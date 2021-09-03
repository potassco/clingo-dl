// {{{ MIT License
//
// Copyright Roland Kaminski, Philipp Wanko, and Max Ostrowski
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
//
// }}}

#include <clingo-dl/graph.hh>

namespace ClingoDL {

static constexpr auto invalid_edge_index = std::numeric_limits<vertex_t>::max();

template <typename T>
Graph<T>::Graph(ThreadStatistics &stats, EdgeVec const &edges, PropagationMode propagate)
: edges_{edges}
, stats_{stats}
, propagate_{propagate} {
    edge_states_.resize(edges_.size(), {1, 1, 0});
    for (edge_t i = 0; i < numeric_cast<edge_t>(edges_.size()); ++i) {
        ensure_index(nodes_, std::max(edges_[i].from, edges_[i].to));
        add_candidate_edge_(i);
    }
}

template <typename T>
bool Graph<T>::empty() const {
    return nodes_.empty();
}

template <typename T>
bool Graph<T>::valid_node(vertex_t idx) const { return nodes_.size() > idx; }

template <typename T>
bool Graph<T>::node_value_defined(vertex_t idx) const { return nodes_[idx].defined(); }

template <typename T>
bool Graph<T>::has_value(vertex_t idx) const { return valid_node(idx) && node_value_defined(idx); }


template <typename T>
T Graph<T>::node_value(vertex_t idx) const { return -nodes_[idx].potential(); }


template <typename T>
bool Graph<T>::edge_is_active(edge_t edge_idx) const { return edge_states_[edge_idx].active; }


template <typename T>
bool Graph<T>::can_propagate() const {
    return std::get<4>(changed_trail_.back());
}

template <typename T>
void Graph<T>::disable_propagate() {
    std::get<4>(changed_trail_.back()) = false;
}

template <typename T>
void Graph<T>::ensure_decision_level(level_t level, bool enable_propagate) {
    // initialize the trail
    if (changed_trail_.empty() || numeric_cast<level_t>(std::get<0>(changed_trail_.back())) < level) {
        bool can_propagate = (changed_trail_.empty() || std::get<4>(changed_trail_.back())) && enable_propagate;
        changed_trail_.emplace_back(level, numeric_cast<uint32_t>(changed_nodes_.size()),
                                           numeric_cast<uint32_t>(changed_edges_.size()),
                                           numeric_cast<uint32_t>(inactive_edges_.size()),
                                           can_propagate);
    }
}

template <typename T>
bool Graph<T>::propagate(edge_t xy_idx, Clingo::PropagateControl &ctl) { // NOLINT
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
    vertex_t num_relevant_out_from{0};
    vertex_t num_relevant_in_from{0};
    vertex_t num_relevant_out_to{0};
    vertex_t num_relevant_in_to{0};
    {
        Timer t{stats_.time_dijkstra};
        std::tie(num_relevant_out_from, num_relevant_in_from) = dijkstra_(xy.from, visited_from_, *static_cast<HFM *>(this));
        std::tie(num_relevant_out_to, num_relevant_in_to) = dijkstra_(xy.to, visited_to_, *static_cast<HTM *>(this));
    }
#ifdef CLINGODL_CROSSCHECK
    uint32_t check_relevant_out_from = 0;
    uint32_t check_relevant_in_from = 0;
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
    uint32_t check_relevant_out_to = 0;
    uint32_t check_relevant_in_to = 0;
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

    bool ret = propagate_edges_(*static_cast<HFM *>(this), ctl, xy_idx, forward_from, backward_from) && propagate_edges_(*static_cast<HTM *>(this), ctl, xy_idx, !forward_from, !backward_from);

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
bool Graph<T>::add_edge(edge_t uv_idx, std::function<bool(std::vector<edge_t>)> f) { // NOLINT
#ifdef CLINGODL_CROSSCHECK
    for (auto &node : nodes_) {
        static_cast<void>(node);
        assert(!node.visited_from);
    }
#endif
    assert(visited_from_.empty());
    assert(costs_heap_.empty());
    level_t level = current_decision_level_();
    auto &uv = edges_[uv_idx];
    // NOTE: would be more efficient if relevant would return statically false here
    //       for the compiler to make comparison cheaper
    auto &m = *static_cast<HFM *>(this);

    // initialize the nodes of the edge to add
    auto &u = nodes_[uv.from];
    auto &v = nodes_[uv.to];
    if (!u.defined()) {
        set_potential_(u, level, 0);
    }
    if (!v.defined()) {
        set_potential_(v, level, 0);
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

    uint32_t dfs = 0;
    // detect negative cycles
    while (!costs_heap_.empty() && !u.visited_from) {
        auto s_idx = costs_heap_.pop(m);
        auto &s = nodes_[s_idx];
        assert(s.visited_from);
        s.visited_from = ++dfs;
        set_potential_(s, level, s.potential() + s.cost_from);
        for (auto st_idx : s.outgoing) {
            ++stats_.propagate_cost_add;
            assert(st_idx < numeric_cast<edge_t>(edges_.size()));
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
#ifdef CLINGODL_CROSSCHECK
        // NOTE: just a check that will throw if there is a cycle
        bellman_ford_(changed_edges_, uv.from);
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
#ifdef CLINGODL_CROSSCHECK
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
            consistent = with_incoming_(uv.from, f, [&](vertex_t t_idx, edge_t ts_idx) {
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
            consistent = cheap_propagate_(uv.from, uv.from, f);
            if (propagate_ >= PropagationMode::WeakPlus && consistent) {
                for (auto &s_idx : visited_from_) {
                    if (!cheap_propagate_(uv.from, s_idx, f)) {
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
bool Graph<T>::with_incoming_(vertex_t s_idx, P p, F f) {
    auto &s = nodes_[s_idx];
    auto &in = s.candidate_incoming;
    auto jt = in.begin();
    for (auto it = jt, ie = in.end(); it != ie; ++it) {
        auto &ts_idx = *it;
        auto &ts = edges_[ts_idx];
        auto t_idx = ts.from;
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
[[nodiscard]] bool Graph<T>::cheap_propagate_(vertex_t u_idx, vertex_t s_idx, F f) {
    // we check for the following case:
    // u ->* s -> * t
    //       ^-----/
    //          ts
    return with_incoming_(s_idx, f, [&](vertex_t t_idx, edge_t ts_idx) {
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
void Graph<T>::backtrack() {
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
        add_candidate_edge_(*i);
    }
    inactive_edges_.resize(n);
    changed_trail_.pop_back();
}

template <typename T>
void Graph<T>::remove_candidate_edge(edge_t uv_idx) {
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
PropagationMode Graph<T>::mode() const {
    return propagate_;
}

template <typename T>
void Graph<T>::add_candidate_edge_(edge_t uv_idx) {
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
bool Graph<T>::propagate_edge_true_(edge_t uv_idx, edge_t xy_idx) {
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
#ifdef CLINGODL_CROSSCHECK
        auto bf_costs_from_u = bellman_ford_(changed_edges_, uv.from);
        auto bf_costs_from_x = bellman_ford_(changed_edges_, xy.from);
        auto aa = bf_costs_from_u.find(xy.to);
        static_cast<void>(aa);
        assert(aa != bf_costs_from_u.end());
        assert(aa->second == a);
        auto bb = bf_costs_from_x.find(uv.to);
        static_cast<void>(bb);
        assert(bb != bf_costs_from_u.end());
        assert(bb->second == b);
#endif
        if (d <= uv.weight) {
            ++stats_.true_edges;
#ifdef CLINGODL_CROSSCHECK
            auto edges = changed_edges_;
            edges.emplace_back(uv_idx);
            // NOTE: throws if there is a cycle
            try {
                bellman_ford_(changed_edges_, uv.from);
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
bool Graph<T>::propagate_edge_false_(Clingo::PropagateControl &ctl, edge_t uv_idx, edge_t xy_idx, bool &ret) { // NOLINT
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
#ifdef CLINGODL_CROSSCHECK
                T sum = uv.weight - xy.weight;
#endif
                std::vector<Clingo::literal_t> clause;
                clause.push_back(-uv.lit);
                // forward
                for (auto next_edge_idx = u.path_from; next_edge_idx != invalid_edge_index;) {
                    auto &next_edge = edges_[next_edge_idx];
                    auto &next_node = nodes_[next_edge.from];
                    clause.push_back(-next_edge.lit);
#ifdef CLINGODL_CROSSCHECK
                    sum += next_edge.weight;
#endif
                    next_edge_idx = next_node.path_from;
                }
                // backward
                for (auto next_edge_idx = v.path_to; next_edge_idx != invalid_edge_index;) {
                    auto &next_edge = edges_[next_edge_idx];
                    auto &next_node = nodes_[next_edge.to];
                    clause.push_back(-next_edge.lit);
#ifdef CLINGODL_CROSSCHECK
                    sum += next_edge.weight;
#endif
                    next_edge_idx = next_node.path_to;
                }
#ifdef CLINGODL_CROSSCHECK
                assert(sum < 0);
#endif
                if (!(ret = ctl.add_clause(clause) && ctl.propagate())) {
                    return false;
                }
            }
            remove_candidate_edge(uv_idx);
            return true;
        }
#ifdef CLINGODL_CROSSCHECK
        auto edges = changed_edges_;
        edges.emplace_back(uv_idx);
        // NOTE: throws if there is a cycle
        try {
            bellman_ford_(changed_edges_, uv.from);
        }
        catch (...) {
            assert(false && "edge must not cause a conflict");
        }
#endif
    }
    return false;
}

template <typename T>
template <class M>
bool Graph<T>::propagate_edges_(M &m, Clingo::PropagateControl &ctl, edge_t xy_idx, bool forward, bool backward) { // NOLINT
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
                        [&](edge_t uv_idx) {
                            if (!edge_states_[uv_idx].active || propagate_edge_true_(uv_idx, xy_idx)) {
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
                        [&](edge_t uv_idx) {
                            if (!ret) {
                                return false;
                            }
                            if (!edge_states_[uv_idx].active || propagate_edge_false_(ctl, uv_idx, xy_idx, ret)) {
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
auto Graph<T>::dijkstra_(vertex_t source_idx, std::vector<vertex_t> &visited_set, M &m) -> std::pair<uint32_t, uint32_t> { // NOLINT
    uint32_t relevant = 0;
    uint32_t relevant_degree_out = 0;
    uint32_t relevant_degree_in = 0;
    assert(visited_set.empty() && costs_heap_.empty());
    costs_heap_.push(m, source_idx);
    visited_set.push_back(source_idx);
    m.visited(source_idx) = true;
    m.cost(source_idx) = 0;
    m.path(source_idx) = invalid_edge_index;
    while (!costs_heap_.empty()) {
        auto u_idx = costs_heap_.pop(m);
        auto tu = m.path(u_idx);
        if (tu != invalid_edge_index && m.relevant(m.from(tu))) {
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

#ifdef CLINGODL_CROSSCHECK
template <typename T>
std::unordered_map<int, T> Graph<T>::bellman_ford_(std::vector<vertex_t> const &edges, int source) {
    std::unordered_map<int, T> costs;
    costs[source] = 0;
    int nodes = 0;
    for (auto &node : nodes_) {
        if (node.defined()) {
            ++nodes;
        }
    }
    for (int i = 0; i < nodes; ++i) {
        for (auto const &uv_idx : edges) {
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
    for (auto const &uv_idx : edges) {
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
void Graph<T>::set_potential_(Vertex &node, level_t level, T potential) {
    if (!node.defined() || node.potential_stack.back().first < level) {
        node.potential_stack.emplace_back(level, potential);
        changed_nodes_.emplace_back(numeric_cast<vertex_t>(&node - nodes_.data()));
    }
    else {
        node.potential_stack.back().second = potential;
    }
}

template <typename T>
level_t Graph<T>::current_decision_level_() {
    assert(!changed_trail_.empty());
    return std::get<0>(changed_trail_.back());
}

template class Graph<int>;
template class Graph<double>;

} // namespace ClingoDL
