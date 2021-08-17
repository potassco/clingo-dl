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
