// {{{ MIT License
//
// // Copyright 2018 Roland Kaminski, Philipp Wanko, Max Ostrowski
//
// // Permission is hereby granted, free of charge, to any person obtaining a copy
// // of this software and associated documentation files (the "Software"), to
// // deal in the Software without restriction, including without limitation the
// // rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// // sell copies of the Software, and to permit persons to whom the Software is
// // furnished to do so, subject to the following conditions:
//
// // The above copyright notice and this permission notice shall be included in
// // all copies or substantial portions of the Software.
//
// // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// // IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// // FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// // AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// // LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// // FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// // IN THE SOFTWARE.
//
// // }}}

#include <clingo-dl/graph.hh>

namespace ClingoDL {

template <typename T>
DifferenceLogicGraph<T>::DifferenceLogicGraph(DLStats &stats, const std::vector<Edge<T>> &edges, PropagationMode propagate)
    : edges_(edges)
    , propagate_(propagate)
    , stats_(stats) {
    edge_states_.resize(edges_.size(), {1, 1, 0});
    for (int i = 0; i < numeric_cast<int>(edges_.size()); ++i) {
        ensure_index(nodes_, std::max(edges_[i].from, edges_[i].to));
        add_candidate_edge(i);
    }
}

template <typename T>
bool DifferenceLogicGraph<T>::empty() const {
    return nodes_.empty();
}

template <typename T>
bool DifferenceLogicGraph<T>::valid_node(int idx) const { return nodes_.size() > idx; }

template <typename T>
int DifferenceLogicGraph<T>::node_value_defined(int idx) const { return nodes_[idx].defined(); }

template <typename T>
bool DifferenceLogicGraph<T>::has_value(int idx) const { return valid_node(idx) && node_value_defined(idx); }


template <typename T>
T DifferenceLogicGraph<T>::node_value(int idx) const { return -nodes_[idx].potential(); }


template <typename T>
bool DifferenceLogicGraph<T>::edge_is_active(int edge_idx) const { return edge_states_[edge_idx].active; }


template <typename T>
bool DifferenceLogicGraph<T>::can_propagate() const {
    return std::get<4>(changed_trail_.back());
}


template <typename T>
void DifferenceLogicGraph<T>::disable_propagate() {
    std::get<4>(changed_trail_.back()) = false;
}

template <typename T>
void DifferenceLogicGraph<T>::ensure_decision_level(int level, bool enable_propagate) {
    // initialize the trail
    if (changed_trail_.empty() || static_cast<int>(std::get<0>(changed_trail_.back())) < level) {
        bool can_propagate = (changed_trail_.empty() || std::get<4>(changed_trail_.back())) && enable_propagate;
        changed_trail_.emplace_back(level, static_cast<int>(changed_nodes_.size()),
                                           static_cast<int>(changed_edges_.size()),
                                           static_cast<int>(inactive_edges_.size()),
                                           can_propagate);
    }
}

template <typename T>
bool DifferenceLogicGraph<T>::propagate(int xy_idx, Clingo::PropagateControl &ctl) {
    ++stats_.edges_propagated;
    remove_candidate_edge(xy_idx);
    auto &xy = edges_[xy_idx];
    auto &x = nodes_[xy.from];
    auto &y = nodes_[xy.to];
    // BUG: this test is not correct
    // if ((x.incoming.empty() && x.outgoing.size() == 1) || (y.outgoing.empty() && y.incoming.size() == 1)) {
    //    return true;
    //}
    x.relevant_to = true;
    y.relevant_from = true;
    int num_relevant_out_from;
    int num_relevant_in_from;
    int num_relevant_out_to;
    int num_relevant_in_to;
    {
        Timer t{stats_.time_dijkstra};
        std::tie(num_relevant_out_from, num_relevant_in_from) = dijkstra(xy.from, visited_from_, *static_cast<HFM *>(this));
        std::tie(num_relevant_out_to, num_relevant_in_to) = dijkstra(xy.to, visited_to_, *static_cast<HTM *>(this));
    }
#ifdef CROSSCHECK
    int check_relevant_out_from = 0, check_relevant_in_from = 0;
    for (auto &node : visited_from_) {
        if (nodes_[node].relevant_from) {
            for (auto &edge : nodes_[node].candidate_incoming) {
                if (edge_states_[edge].active) {
                    ++check_relevant_in_from;
                }
            }
            for (auto &edge : nodes_[node].candidate_outgoing) {
                if (edge_states_[edge].active) {
                    ++check_relevant_out_from;
                }
            }
        }
    }
    assert(num_relevant_out_from == check_relevant_out_from);
    assert(num_relevant_in_from == check_relevant_in_from);
    int check_relevant_out_to = 0, check_relevant_in_to = 0;
    for (auto &node : visited_to_) {
        if (nodes_[node].relevant_to) {
            for (auto &edge : nodes_[node].candidate_incoming) {
                if (edge_states_[edge].active) {
                    ++check_relevant_in_to;
                }
            }
            for (auto &edge : nodes_[node].candidate_outgoing) {
                if (edge_states_[edge].active) {
                    ++check_relevant_out_to;
                }
            }
        }
    }
    assert(num_relevant_out_to == check_relevant_out_to);
    assert(num_relevant_in_to == check_relevant_in_to);
#endif

    bool forward_from = num_relevant_in_from < num_relevant_out_to;
    bool backward_from = num_relevant_out_from < num_relevant_in_to;

    bool ret = propagate_edges(*static_cast<HFM *>(this), ctl, xy_idx, forward_from, backward_from) && propagate_edges(*static_cast<HTM *>(this), ctl, xy_idx, !forward_from, !backward_from);

    for (auto &x : visited_from_) {
        nodes_[x].visited_from = false;
        nodes_[x].relevant_from = false;
    }
    for (auto &x : visited_to_) {
        nodes_[x].visited_to = false;
        nodes_[x].relevant_to = false;
    }
    visited_from_.clear();
    visited_to_.clear();
    return ret;
}

template <typename T>
bool DifferenceLogicGraph<T>::add_edge(int uv_idx, std::function<bool(std::vector<int>)> f) {
#ifdef CROSSCHECK
    for (auto &node : nodes_) {
        assert(!node.visited_from);
    }
#endif
    assert(visited_from_.empty());
    assert(costs_heap_.empty());
    int level = current_decision_level_();
    auto &uv = edges_[uv_idx];
    // NOTE: would be more efficient if relevant would return statically false here
    //       for the compiler to make comparison cheaper
    auto &m = *static_cast<HFM *>(this);

    // initialize the nodes of the edge to add
    auto &u = nodes_[uv.from];
    auto &v = nodes_[uv.to];
    if (!u.defined()) {
        set_potential(u, level, 0);
    }
    if (!v.defined()) {
        set_potential(v, level, 0);
    }
    v.cost_from = u.potential() + uv.weight - v.potential();
    ++stats_.edges_added;
    if (v.cost_from < 0) {
        costs_heap_.push(m, uv.to);
        visited_from_.emplace_back(uv.to);
        v.visited_from = 1;
        v.path_from = uv_idx;
    }
    else {
        ++stats_.edges_skipped;
    }

    int dfs = 0;
    // detect negative cycles
    while (!costs_heap_.empty() && !u.visited_from) {
        auto s_idx = costs_heap_.pop(m);
        auto &s = nodes_[s_idx];
        assert(s.visited_from);
        s.visited_from = ++dfs;
        set_potential(s, level, s.potential() + s.cost_from);
        for (auto st_idx : s.outgoing) {
            ++stats_.propagate_cost_add;
            assert(st_idx < numeric_cast<int>(edges_.size()));
            auto &st = edges_[st_idx];
            auto &t = nodes_[st.to];
            auto c = s.potential() + st.weight - t.potential();
            if (c < (t.visited_from ? t.cost_from : 0)) {
                assert(c < 0);
                t.path_from = st_idx;
                t.cost_from = c;
                if (!t.visited_from) {
                    t.visited_from = 1;
                    visited_from_.emplace_back(st.to);
                    costs_heap_.push(m, st.to);
                }
                else {
                    costs_heap_.decrease(m, m.offset(st.to));
                }
            }
        }
    }

    bool consistent = true;
    if (!u.visited_from) {
        // add the edge to the graph
        u.outgoing.emplace_back(uv_idx);
        v.incoming.emplace_back(uv_idx);
        changed_edges_.emplace_back(uv_idx);
#ifdef CROSSCHECK
        // NOTE: just a check that will throw if there is a cycle
        bellman_ford(changed_edges_, uv.from);
#endif
    }
    else {
        // gather the edges in the negative cycle
        neg_cycle_.clear();
        neg_cycle_.push_back(v.path_from);
        auto next_idx = edges_[v.path_from].from;
        while (uv.to != next_idx) {
            auto &next = nodes_[next_idx];
            neg_cycle_.push_back(next.path_from);
            next_idx = edges_[next.path_from].from;
        }
#ifdef CROSSCHECK
        T weight = 0;
        for (auto &edge_idx : neg_cycle_) {
            weight += edges_[edge_idx].weight;
        }
        assert(weight < 0);
#endif
        consistent = f(neg_cycle_);
    }

    if (propagate_ >= PropagationMode::Trivial && consistent) {
        if (visited_from_.empty() || propagate_ == PropagationMode::Trivial) {
            consistent = with_incoming(uv.from, f, [&](int t_idx, int ts_idx) {
                auto &ts = edges_[ts_idx];
                if (t_idx == uv.to && uv.weight + ts.weight < 0) {
                    neg_cycle_.emplace_back(uv_idx);
                    neg_cycle_.emplace_back(ts_idx);
                    ++stats_.false_edges_trivial;
                    return true;
                }
                return false;
            });
        }
        else if (propagate_ >= PropagationMode::Weak) {
            consistent = cheap_propagate(uv.from, uv.from, f);
            if (propagate_ >= PropagationMode::WeakPlus && consistent) {
                for (auto &s_idx : visited_from_) {
                    if (!cheap_propagate(uv.from, s_idx, f)) {
                        consistent = false;
                        break;
                    }
                }
            }
        }
    }
    // reset visited flags
    for (auto &x : visited_from_) {
        nodes_[x].visited_from = 0;
    }
    visited_from_.clear();
    costs_heap_.clear();

    return consistent;
}

template <typename T>
template <class P, class F>
bool DifferenceLogicGraph<T>::with_incoming(int s_idx, P p, F f) {
    auto &s = nodes_[s_idx];
    auto &in = s.candidate_incoming;
    auto jt = in.begin();
    for (auto it = jt, ie = in.end(); it != ie; ++it) {
        auto &ts_idx = *it;
        auto &ts = edges_[ts_idx];
        auto t_idx = ts.from;
        auto &t = nodes_[t_idx];
        if (!edge_states_[ts_idx].active) {
            edge_states_[ts_idx].removed_incoming = true;
            continue;
        }
        neg_cycle_.clear();
        if (f(t_idx, ts_idx)) {
            edge_states_[ts_idx].removed_incoming = true;
            remove_candidate_edge(ts_idx);
            if (!p(neg_cycle_)) {
                in.erase(jt, it+1);
                return false;
            }
            continue;
        }
        *jt++ = *it;
    }
    in.erase(jt, in.end());
    return true;
}

template <typename T>
template <class F>
[[nodiscard]] bool DifferenceLogicGraph<T>::cheap_propagate(int u_idx, int s_idx, F f) {
    // we check for the following case:
    // u ->* s -> * t
    //       ^-----/
    //          ts
    return with_incoming(s_idx, f, [&](int t_idx, int ts_idx) {
        auto &s = nodes_[s_idx];
        auto &t = nodes_[t_idx];
        auto &ts = edges_[ts_idx];
        if (s.visited_from < t.visited_from) {
            T weight = t.potential() - s.potential();
            if (weight + ts.weight < 0) {
                T check = 0;
                auto r_idx = t_idx;
                while (u_idx != r_idx && s_idx != r_idx) {
                    auto &r = nodes_[r_idx];
                    auto &rr = edges_[r.path_from];
                    neg_cycle_.emplace_back(r.path_from);
                    r_idx = rr.from;
                    check += rr.weight;
                }
                if (r_idx == s_idx) {
                    if (u_idx == s_idx) {
                        ++stats_.false_edges_weak;
                    }
                    else {
                        ++stats_.false_edges_weak_plus;
                    }
                    assert(weight == check);
                    neg_cycle_.emplace_back(ts_idx);
                    return true;
                }
            }
        }
        return false;
    });
}



template <typename T>
void DifferenceLogicGraph<T>::backtrack() {
    for (auto count = static_cast<int>(changed_nodes_.size()) - std::get<1>(changed_trail_.back()); count > 0; --count) {
        auto &node = nodes_[changed_nodes_.back()];
        node.potential_stack.pop_back();
        changed_nodes_.pop_back();
    }
    for (auto count = static_cast<int>(changed_edges_.size()) - std::get<2>(changed_trail_.back()); count > 0; --count) {
        auto &edge = edges_[changed_edges_.back()];
        nodes_[edge.from].outgoing.pop_back();
        nodes_[edge.to].incoming.pop_back();
        changed_edges_.pop_back();
    }
    int n = std::get<3>(changed_trail_.back());
    for (auto i = inactive_edges_.begin() + n, e = inactive_edges_.end(); i < e; ++i) {
        add_candidate_edge(*i);
    }
    inactive_edges_.resize(n);
    changed_trail_.pop_back();
}

template <typename T>
void DifferenceLogicGraph<T>::remove_candidate_edge(int uv_idx) {
    auto &uv = edges_[uv_idx];
    auto &u = nodes_[uv.from];
    auto &v = nodes_[uv.to];
    --u.degree_out;
    --v.degree_in;
    inactive_edges_.push_back(uv_idx);
    assert(edge_states_[uv_idx].active);
    edge_states_[uv_idx].active = false;
}

template <typename T>
PropagationMode DifferenceLogicGraph<T>::mode() const {
    return propagate_;
}

template <typename T>
void DifferenceLogicGraph<T>::add_candidate_edge(int uv_idx) {
    auto &uv = edges_[uv_idx];
    auto &uv_state = edge_states_[uv_idx];
    auto &u = nodes_[uv.from];
    auto &v = nodes_[uv.to];
    ++u.degree_out;
    ++v.degree_in;
    assert(!uv_state.active);
    uv_state.active = true;
    if (uv_state.removed_outgoing) {
        uv_state.removed_outgoing = false;
        u.candidate_outgoing.emplace_back(uv_idx);
    }
    if (uv_state.removed_incoming) {
        uv_state.removed_incoming = false;
        v.candidate_incoming.emplace_back(uv_idx);
    }
}

template <typename T>
bool DifferenceLogicGraph<T>::propagate_edge_true(int uv_idx, int xy_idx) {
    auto &uv = edges_[uv_idx];
    auto &u = nodes_[uv.from];
    auto &v = nodes_[uv.to];
    assert(u.relevant_to || v.relevant_from);

    if (u.relevant_to && v.relevant_from) {
        auto &xy = edges_[xy_idx];
        auto &x = nodes_[xy.from];
        auto &y = nodes_[xy.to];

        auto a = u.cost_to + y.potential() - u.potential();
        auto b = v.cost_from + v.potential() - x.potential();
        auto d = a + b - xy.weight;
#ifdef CROSSCHECK
        auto bf_costs_from_u = bellman_ford(changed_edges_, uv.from);
        auto bf_costs_from_x = bellman_ford(changed_edges_, xy.from);
        auto aa = bf_costs_from_u.find(xy.to);
        assert(aa != bf_costs_from_u.end());
        assert(aa->second == a);
        auto bb = bf_costs_from_x.find(uv.to);
        assert(bb != bf_costs_from_u.end());
        assert(bb->second == b);
#endif
        if (d <= uv.weight) {
            ++stats_.true_edges;
#ifdef CROSSCHECK
            auto edges = changed_edges_;
            edges.emplace_back(uv_idx);
            // NOTE: throws if there is a cycle
            try {
                bellman_ford(changed_edges_, uv.from);
            }
            catch (...) {
                assert(false && "edge is implied but lead to a conflict :(");
            }
#endif
            remove_candidate_edge(uv_idx);
            return true;
        }
    }
    return false;
}

template <typename T>
bool DifferenceLogicGraph<T>::propagate_edge_false(Clingo::PropagateControl &ctl, int uv_idx, int xy_idx, bool &ret) {
    auto &uv = edges_[uv_idx];
    auto &u = nodes_[uv.from];
    auto &v = nodes_[uv.to];
    assert(v.relevant_to || u.relevant_from);

    if (v.relevant_to && u.relevant_from) {
        auto &xy = edges_[xy_idx];
        auto &x = nodes_[xy.from];
        auto &y = nodes_[xy.to];

        auto a = v.cost_to + y.potential() - v.potential();
        auto b = u.cost_from + u.potential() - x.potential();
        auto d = a + b - xy.weight;
        if (d < -uv.weight) {
            ++stats_.false_edges;
            if (!ctl.assignment().is_false(uv.lit)) {
#ifdef CROSSCHECK
                T sum = uv.weight - xy.weight;
#endif
                std::vector<Clingo::literal_t> clause;
                clause.push_back(-uv.lit);
                // forward
                for (auto next_edge_idx = u.path_from; next_edge_idx >= 0;) {
                    auto &next_edge = edges_[next_edge_idx];
                    auto &next_node = nodes_[next_edge.from];
                    clause.push_back(-next_edge.lit);
#ifdef CROSSCHECK
                    sum += next_edge.weight;
#endif
                    next_edge_idx = next_node.path_from;
                }
                // backward
                for (auto next_edge_idx = v.path_to; next_edge_idx >= 0;) {
                    auto &next_edge = edges_[next_edge_idx];
                    auto &next_node = nodes_[next_edge.to];
                    clause.push_back(-next_edge.lit);
#ifdef CROSSCHECK
                    sum += next_edge.weight;
#endif
                    next_edge_idx = next_node.path_to;
                }
#ifdef CROSSCHECK
                assert(sum < 0);
#endif
                if (!(ret = ctl.add_clause(clause) && ctl.propagate())) {
                    return false;
                }
            }
            remove_candidate_edge(uv_idx);
            return true;
        }
#ifdef CROSSCHECK
        else {
            auto edges = changed_edges_;
            edges.emplace_back(uv_idx);
            // NOTE: throws if there is a cycle
            try {
                bellman_ford(changed_edges_, uv.from);
            }
            catch (...) {
                assert(false && "edge must not cause a conflict");
            }
        }
#endif
    }
    return false;
}

template <typename T>
template <class M>
bool DifferenceLogicGraph<T>::propagate_edges(M &m, Clingo::PropagateControl &ctl, int xy_idx, bool forward, bool backward) {
    if (!forward && !backward) {
        return true;
    }
    for (auto &node : m.visited_set()) {
        if (m.relevant(node)) {
            if (forward) {
                auto &in = m.candidate_incoming(node);
                in.resize(
                    std::remove_if(
                        in.begin(), in.end(),
                        [&](int uv_idx) {
                            if (!edge_states_[uv_idx].active || propagate_edge_true(uv_idx, xy_idx)) {
                                m.remove_incoming(uv_idx);
                                return true;
                            }
                            return false;
                        }) -
                    in.begin());
            }
            if (backward) {
                bool ret = true;
                auto &out = m.candidate_outgoing(node);
                out.resize(
                    std::remove_if(
                        out.begin(), out.end(),
                        [&](int uv_idx) {
                            if (!ret) {
                                return false;
                            }
                            if (!edge_states_[uv_idx].active || propagate_edge_false(ctl, uv_idx, xy_idx, ret)) {
                                m.remove_outgoing(uv_idx);
                                return true;
                            }
                            return false;
                        }) -
                    out.begin());
                if (!ret) {
                    return false;
                }
            }
        }
    }
    return true;
}

template <typename T>
template <class M>
std::pair<int, int> DifferenceLogicGraph<T>::dijkstra(int source_idx, std::vector<int> &visited_set, M &m) {
    int relevant = 0;
    int relevant_degree_out = 0, relevant_degree_in = 0;
    assert(visited_set.empty() && costs_heap_.empty());
    costs_heap_.push(m, source_idx);
    visited_set.push_back(source_idx);
    m.visited(source_idx) = true;
    m.cost(source_idx) = 0;
    m.path(source_idx) = -1;
    while (!costs_heap_.empty()) {
        auto u_idx = costs_heap_.pop(m);
        auto tu = m.path(u_idx);
        if (tu >= 0 && m.relevant(m.from(tu))) {
            m.relevant(u_idx) = true;
            --relevant; // just removed a relevant edge from the queue
        }
        bool relevant_u = m.relevant(u_idx);
        if (relevant_u) {
            relevant_degree_out += nodes_[u_idx].degree_out;
            relevant_degree_in += nodes_[u_idx].degree_in;
        }
        for (auto &uv_idx : m.out(u_idx)) {
            ++m.propagation_cost();
            auto &uv = edges_[uv_idx];
            auto v_idx = m.to(uv_idx);
            // NOTE: explicitely using uv.from and uv.to is intended here
            auto c = m.cost(u_idx) + nodes_[uv.from].potential() + uv.weight - nodes_[uv.to].potential();
            assert(nodes_[uv.from].potential() + uv.weight - nodes_[uv.to].potential() >= 0);
            if (!m.visited(v_idx) || c < m.cost(v_idx)) {
                m.cost(v_idx) = c;
                if (!m.visited(v_idx)) {
                    // node v contributes an edge with a relevant source
                    if (relevant_u) {
                        ++relevant;
                    }
                    visited_set.push_back(m.to(uv_idx));
                    m.visited(v_idx) = true;
                    costs_heap_.push(m, v_idx);
                }
                else {
                    if (m.relevant(m.from(m.path(v_idx)))) {
                        // node v no longer contributes a relevant edge
                        if (!relevant_u) {
                            --relevant;
                        }
                    }
                    // node v contributes a relevant edge now
                    else if (relevant_u) {
                        ++relevant;
                    }
                    costs_heap_.decrease(m, m.offset(v_idx));
                }
                m.path(v_idx) = uv_idx;
            }
        }
        // removed a relevant node from the queue and there are no edges with relevant sources anymore in the queue
        // this condition assumes that initially there is exactly one reachable relevant node in the graph
        if (relevant_u && relevant == 0) {
            costs_heap_.clear();
            break;
        }
    }
    return {relevant_degree_out, relevant_degree_in};
}

#ifdef CROSSCHECK
template <typename T>
std::unordered_map<int, T> bellman_ford(std::vector<int> const &edges, int source) {
    std::unordered_map<int, T> costs;
    costs[source] = 0;
    int nodes = 0;
    for (auto &node : nodes_) {
        if (node.defined()) {
            ++nodes;
        }
    }
    for (int i = 0; i < nodes; ++i) {
        for (auto &uv_idx : edges) {
            auto &uv = edges_[uv_idx];
            auto u_cost = costs.find(uv.from);
            if (u_cost != costs.end()) {
                auto v_cost = costs.find(uv.to);
                auto dist = u_cost->second + uv.weight;
                if (v_cost == costs.end()) {
                    costs[uv.to] = dist;
                }
                else if (dist < v_cost->second) {
                    v_cost->second = dist;
                }
            }
        }
    }
    for (auto &uv_idx : edges) {
        auto &uv = edges_[uv_idx];
        auto u_cost = costs.find(uv.from);
        if (u_cost != costs.end()) {
            auto v_cost = costs.find(uv.to);
            auto dist = u_cost->second + uv.weight;
            if (dist < v_cost->second) {
                throw std::runtime_error("there is a negative cycle!!!");
            }
        }
    }
    return costs;
}
#endif

template <typename T>
void DifferenceLogicGraph<T>::set_potential(DifferenceLogicNode<T> &node, int level, T potential) {
    if (!node.defined() || node.potential_stack.back().first < level) {
        node.potential_stack.emplace_back(level, potential);
        changed_nodes_.emplace_back(numeric_cast<int>(&node - nodes_.data()));
    }
    else {
        node.potential_stack.back().second = potential;
    }
}

template <typename T>
int DifferenceLogicGraph<T>::current_decision_level_() {
    assert(!changed_trail_.empty());
    return std::get<0>(changed_trail_.back());
}

template class DifferenceLogicGraph<int>;
template class DifferenceLogicGraph<double>;

} // namespace ClingoDL
