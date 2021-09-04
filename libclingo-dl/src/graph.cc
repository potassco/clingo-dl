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

namespace {

struct From { };
struct To { };

constexpr auto invalid_edge_index = std::numeric_limits<vertex_t>::max();

#ifdef CLINGODL_CROSSCHECK

//! Count relevant incoming/outgoing edges.
template <class M>
std::pair<uint32_t, uint32_t> count_relevant_(M &m) {
    uint32_t relevant_in = 0;
    uint32_t relevant_out = 0;
    for (auto &node : m.visited_set()) {
        if (m.relevant(node)) {
            for (auto &edge : m.candidate_incoming(node)) {
                if (m.active(edge)) {
                    ++relevant_in;
                }
            }
            for (auto &edge : m.candidate_outgoing(node)) {
                if (m.active(edge)) {
                    ++relevant_out;
                }
            }
        }
    }
    return {relevant_in, relevant_out};
}

#endif

} // namespace

// TODO: document + make nested
template <typename T>
struct Graph<T>::Vertex {
    using value_t = T;
    using PotentialStack = std::vector<std::pair<vertex_t, value_t>>;

    //! Return true if the node has a value assigned.
    [[nodiscard]] bool defined() const {
        return !potential_stack.empty();
    }
    //! Return the current value associated with the vertex.
    [[nodiscard]] value_t potential() const {
        return potential_stack.back().second;
    }

    VertexIndexVec outgoing;            //!< Outgoing edges from this vertex that are true.
    VertexIndexVec incoming;            //!< Incoming edges to this vertex that are true.
    VertexIndexVec candidate_incoming;  //!< Edges that might become outgoing edges.
    VertexIndexVec candidate_outgoing;  //!< Edges that might become incoming edges.
    PotentialStack potential_stack;     //!< Vector of pairs of level and potential.
    value_t cost_from{0};
    value_t cost_to{0};
    vertex_t offset{0};
    edge_t path_from{0};
    edge_t path_to{0};
    vertex_t degree_out{0};
    vertex_t degree_in{0};
    vertex_t visited_from{0};
    bool relevant_from{false};
    bool relevant_to{false};
    bool visited_to{false};
};

template <typename T>
struct Graph<T>::EdgeState {
    uint8_t removed_outgoing : 1;
    uint8_t removed_incoming : 1;
    uint8_t active : 1;
};

//! This struct provides functions to access the original and transposed graph
//! based on a template parameter.
//!
//! It furthermore implements functions that can be applied to both kinds of
//! graphs.
//!
//! \note The same effect could be achieved by providing a lot of private
//! template member functions. This implementation has the advantage that the
//! complexity is hidden in the implementation.
template <typename T>
template <typename D>
struct Graph<T>::Impl : Graph {
    using index_t = Heap<4>::index_type;

    //! The index of the vertex in the heap vector.
    index_t &offset(vertex_t idx) {
        return nodes_[idx].offset;
    }

    //! The cost of the vertex.
    value_t &cost(vertex_t idx) {
        if constexpr (std::is_same_v<D, From>) {
            return nodes_[idx].cost_from;
        }
        else {
            return nodes_[idx].cost_to;
        }
    }

    //! The end point of the given edge.
    vertex_t to(edge_t idx) {
        if constexpr (std::is_same_v<D, From>) {
            return edges_[idx].to;
        }
        else {
            return edges_[idx].from;
        }
    }

    //! The starting point of the given edge.
    vertex_t from(edge_t idx) {
        if constexpr (std::is_same_v<D, From>) {
            return edges_[idx].from;
        }
        else {
            return edges_[idx].to;
        }
    }

    //! The outgoing vertices of the given vertex.
    std::vector<vertex_t> &out(vertex_t idx) {
        if constexpr (std::is_same_v<D, From>) {
            return nodes_[idx].outgoing;
        }
        else {
            return nodes_[idx].incoming;
        }
    }

    //! The edge that was used to reach the given vertex.
    edge_t &path(vertex_t idx) {
        if constexpr (std::is_same_v<D, From>) {
            return nodes_[idx].path_from;
        }
        else {
            return nodes_[idx].path_to;
        }
    }

    //! Flag indicating whether the vertex has been visited.
    //!
    //! \note This is a integer here because it is also used as a dfs index
    //! when adding edges.
    std::conditional_t<std::is_same_v<D, From>, vertex_t, bool> &visited(vertex_t idx) {
        if constexpr(std::is_same_v<D, From>) {
            return nodes_[idx].visited_from;
        }
        else {
            return nodes_[idx].visited_to;
        }
    }

    //! Whether the vertex is relevant for propagation.
    bool &relevant(vertex_t idx) {
        if constexpr (std::is_same_v<D, From>) {
            return nodes_[idx].relevant_from;
        }
        else {
            return nodes_[idx].relevant_to;
        }
    }

    //! The set of all visited vertices.
    std::vector<vertex_t> &visited_set() {
        if constexpr (std::is_same_v<D, From>) {
            return visited_from_;
        }
        else {
            return visited_to_;
        }
    }

    //! Outgoing candidate edges that are not false.
    std::vector<vertex_t> &candidate_outgoing(vertex_t idx) {
        if constexpr (std::is_same_v<D, From>) {
            return nodes_[idx].candidate_outgoing;
        }
        else {
            return nodes_[idx].candidate_incoming;
        }
    }

    //! Incoming candidate edges that are not false.
    std::vector<vertex_t> &candidate_incoming(vertex_t idx) {
        if constexpr (std::is_same_v<D, From>) {
            return nodes_[idx].candidate_incoming;
        }
        else {
            return nodes_[idx].candidate_outgoing;
        }
    }

    //! Mark an incoming edge as removed.
    void remove_incoming(edge_t idx) {
        if constexpr (std::is_same_v<D, From>) {
            edge_states_[idx].removed_incoming = true;
        }
        else {
            edge_states_[idx].removed_outgoing = true;
        }
    }

    //! Mark an outgoing edge as removed.
    void remove_outgoing(edge_t idx) {
        if constexpr (std::is_same_v<D, From>) {
            edge_states_[idx].removed_outgoing = true;
        }
        else {
            edge_states_[idx].removed_incoming = true;
        }
    }

    //! Check if the edge is active.
    bool active(edge_t idx) {
        return edge_states_[idx].active;
    }

    //! The cost to propagate the edge.
    uint64_t &propagation_cost() {
        if constexpr (std::is_same_v<D, From>) {
            return stats_.propagate_cost_from;
        }
        else {
            return stats_.propagate_cost_to;
        }
    }

    std::pair<uint32_t, uint32_t> dijkstra(vertex_t source_idx) { // NOLINT
        uint32_t num_relevant = 0;
        uint32_t relevant_degree_out = 0;
        uint32_t relevant_degree_in = 0;
        assert(visited_set().empty() && costs_heap_.empty());
        costs_heap_.push(*this, source_idx);
        visited_set().push_back(source_idx);
        visited(source_idx) = true;
        cost(source_idx) = 0;
        path(source_idx) = invalid_edge_index;
        while (!costs_heap_.empty()) {
            auto u_idx = costs_heap_.pop(*this);
            auto tu = path(u_idx);
            if (tu != invalid_edge_index && relevant(from(tu))) {
                relevant(u_idx) = true;
                --num_relevant; // just removed a relevant edge from the queue
            }
            bool relevant_u = relevant(u_idx);
            if (relevant_u) {
                relevant_degree_out += nodes_[u_idx].degree_out;
                relevant_degree_in += nodes_[u_idx].degree_in;
            }
            for (auto &uv_idx : out(u_idx)) {
                ++propagation_cost();
                auto &uv = edges_[uv_idx];
                auto v_idx = to(uv_idx);
                // NOTE: explicitely using uv.from and uv.to is intended here
                auto c = cost(u_idx) + nodes_[uv.from].potential() + uv.weight - nodes_[uv.to].potential();
                assert(nodes_[uv.from].potential() + uv.weight - nodes_[uv.to].potential() >= 0);
                if (!visited(v_idx) || c < cost(v_idx)) {
                    cost(v_idx) = c;
                    if (!visited(v_idx)) {
                        // node v contributes an edge with a relevant source
                        if (relevant_u) {
                            ++num_relevant;
                        }
                        visited_set().push_back(to(uv_idx));
                        visited(v_idx) = true;
                        costs_heap_.push(*this, v_idx);
                    }
                    else {
                        if (relevant(from(path(v_idx)))) {
                            // node v no longer contributes a relevant edge
                            if (!relevant_u) {
                                --num_relevant;
                            }
                        }
                        // node v contributes a relevant edge now
                        else if (relevant_u) {
                            ++num_relevant;
                        }
                        costs_heap_.decrease(*this, offset(v_idx));
                    }
                    path(v_idx) = uv_idx;
                }
            }
            // removed a relevant node from the queue and there are no edges with relevant sources anymore in the queue
            // this condition assumes that initially there is exactly one reachable relevant node in the graph
            if (relevant_u && num_relevant == 0) {
                costs_heap_.clear();
                break;
            }
        }
        return {relevant_degree_out, relevant_degree_in};
    }

    bool propagate_edges(Clingo::PropagateControl &ctl, edge_t xy_idx, bool forward, bool backward) { // NOLINT
        if (!forward && !backward) {
            return true;
        }
        for (auto &node : visited_set()) {
            if (relevant(node)) {
                if (forward) {
                    auto &in = candidate_incoming(node);
                    in.resize(
                        std::remove_if(
                            in.begin(), in.end(),
                            [&](edge_t uv_idx) {
                                if (!edge_states_[uv_idx].active || propagate_edge_true_(uv_idx, xy_idx)) {
                                    remove_incoming(uv_idx);
                                    return true;
                                }
                                return false;
                            }) -
                        in.begin());
                }
                if (backward) {
                    bool ret = true;
                    auto &out = candidate_outgoing(node);
                    out.resize(
                        std::remove_if(
                            out.begin(), out.end(),
                            [&](edge_t uv_idx) {
                                if (!ret) {
                                    return false;
                                }
                                if (!edge_states_[uv_idx].active || propagate_edge_false_(ctl, uv_idx, xy_idx, ret)) {
                                    remove_outgoing(uv_idx);
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
};

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
Graph<T>::~Graph() = default;

template <typename T>
Graph<T>::Graph(Graph &&other) noexcept = default;

template <typename T>
bool Graph<T>::empty() const {
    return nodes_.empty();
}

template <typename T>
bool Graph<T>::has_value(vertex_t idx) const {
    return idx < nodes_.size() && nodes_[idx].defined();
}

template <typename T>
T Graph<T>::get_value(vertex_t idx) const {
    assert(has_value(idx));
    return -nodes_[idx].potential();
}

template <typename T>
bool Graph<T>::edge_is_enabled(edge_t edge_idx) const {
    return edge_states_[edge_idx].active;
}

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
bool Graph<T>::propagate(edge_t xy_idx, Clingo::PropagateControl &ctl) {
    ++stats_.edges_propagated;
    disable_edge(xy_idx);
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
        std::tie(num_relevant_out_from, num_relevant_in_from) = static_cast<Impl<From>*>(this)->dijkstra(xy.from);
        std::tie(num_relevant_out_to, num_relevant_in_to) = static_cast<Impl<To>*>(this)->dijkstra(xy.to);
    }
#ifdef CLINGODL_CROSSCHECK
    // check if the counts of relevant incoming/outgoing vertices are correct
    assert(std::make_pair(num_relevant_in_from, num_relevant_out_from) == count_relevant_(*static_cast<Impl<From> *>(this)));
    assert(std::make_pair(num_relevant_out_to, num_relevant_in_to) == count_relevant_(*static_cast<Impl<To> *>(this)));
#endif

    bool forward_from = num_relevant_in_from < num_relevant_out_to;
    bool backward_from = num_relevant_out_from < num_relevant_in_to;

    bool ret = static_cast<Impl<From> *>(this)->propagate_edges(ctl, xy_idx, forward_from, backward_from) && static_cast<Impl<To> *>(this)->propagate_edges(ctl, xy_idx, !forward_from, !backward_from);

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
    // This function adds an edge to the graph and returns false if the edge
    // induces a negative cycle.
    //
    // For this to work, the graph must not already have a negative cycle. The
    // function also propagates cheap to propagate edges considering nodes
    // visited during cycle detection. Negative cycles are reported via the
    // given callback function.
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
    auto &m = *static_cast<Impl<From> *>(this);

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
    // traverse over incoming edges
    for (auto it = jt, ie = in.end(); it != ie; ++it) {
        auto &ts_idx = *it;
        auto &ts = edges_[ts_idx];
        auto t_idx = ts.from;
        // remove edges marked inactive
        if (!edge_states_[ts_idx].active) {
            edge_states_[ts_idx].removed_incoming = true;
            continue;
        }
        // visit the incoming vertex and edge
        neg_cycle_.clear();
        if (f(t_idx, ts_idx)) {
            edge_states_[ts_idx].removed_incoming = true;
            disable_edge(ts_idx);
            // add constraint for the negative cycle
            if (!p(neg_cycle_)) {
                // erease edges marked as removed
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
    //   u ->* s ->* t
    //         ^----/
    //           ts
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
void Graph<T>::disable_edge(edge_t uv_idx) {
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
            disable_edge(uv_idx);
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
            disable_edge(uv_idx);
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
