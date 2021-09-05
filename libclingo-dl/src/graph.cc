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

//!< Tag for traversals of the original graph.
struct From { };
//!< Tag for traversals of the transposed graph.
struct To { };

//! An index different from all valid edge indices.
constexpr auto invalid_edge_index = std::numeric_limits<vertex_t>::max();

} // namespace

void ThreadStatistics::reset() {
    *this = ThreadStatistics{};
}

void ThreadStatistics::accu(ThreadStatistics const &x) {
    time_propagate        += x.time_propagate;
    time_undo             += x.time_undo;
    time_dijkstra         += x.time_dijkstra;
    true_edges            += x.true_edges;
    false_edges           += x.false_edges;
    false_edges_trivial   += x.false_edges_trivial;
    false_edges_weak      += x.false_edges_weak;
    false_edges_weak_plus += x.false_edges_weak_plus;
    propagate_cost_add    += x.propagate_cost_add;
    propagate_cost_from   += x.propagate_cost_from;
    propagate_cost_to     += x.propagate_cost_to;
    edges_added           += x.edges_added;
    edges_skipped         += x.edges_skipped;
    edges_propagated      += x.edges_propagated;
}

//!< Thread specific information for vertices.
template <typename T>
struct Graph<T>::Vertex {
    using value_t = T;
    using PotentialStack = std::vector<std::pair<vertex_t, value_t>>;

    //! Return true if the vertex has a value assigned.
    [[nodiscard]] bool defined() const {
        return !potential_stack.empty();
    }
    //! Return the current value associated with the vertex.
    [[nodiscard]] value_t potential() const {
        return potential_stack.back().second;
    }

    VertexIndexVec outgoing;           //!< Outgoing edges from this vertex that are true.
    VertexIndexVec incoming;           //!< Incoming edges to this vertex that are true.
    VertexIndexVec candidate_incoming; //!< Edges that might become incoming edges.
    VertexIndexVec candidate_outgoing; //!< Edges that might become outgoing edges.
    PotentialStack potential_stack;    //!< Vector of pairs of level and potential.
    value_t cost_from{0};              //!< Costs for traversals of the original graph.
    value_t cost_to{0};                //!< Costs for traversals of the transposed graph.
    vertex_t offset{0};                //!< Offset in the cost heap.
    edge_t path_from{0};               //!< Path pointers for traversals of the original graph.
    edge_t path_to{0};                 //!< Path pointers for traversals of the transposed graph.
    vertex_t degree_out{0};            //!< Outgoing degree of candidate edges.
    vertex_t degree_in{0};             //!< Incoming degree of candidate edges.
    vertex_t visited_from{0};          //!< Either a flag to mark the vertex as visited or its depth first index.
    bool visited_to{false};            //!< A flag to mark the vertex as visited for traversals of the transposed graph.
    bool relevant_from{false};         //!< A flag to mark the vertex as visited for traversals of the original graph.
    bool relevant_to{false};           //!< A flag to mark the vertex as visited for traversals of the transposed graph.
};

//!< Thread specific information for edges.
template <typename T>
struct Graph<T>::EdgeState {
    uint8_t removed_outgoing : 1; //!< Flag to mark edges as removed from the candidate_outgoing vector.
    uint8_t removed_incoming : 1; //!< Flag to mark edges as removed from the candidate_incoming vector.
    uint8_t enabled : 1;          //!< Flag to mark the edge as enabled.
};

//!< Struct holding information to backtrack a decision level.
template <typename T>
struct Graph<T>::TrailEntry {
    level_t level;           //!< The corresponding decision level.
    index_t vertex_offset;   //!< Index up to which to backtrack changed vertices.
    index_t edge_offset;     //!< Index up to which to backtrack changed edges.
    index_t disabled_offset; //!< Index up to which to backtrack inactive edges.
    bool can_propagate;      //!< Whether propagation was possible on this level.
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
    using index_t = HeapType::index_type;

    //! The index of the vertex in the heap vector.
    index_t &offset(vertex_t idx) {
        return vertices_[idx].offset;
    }

    //! The cost of the vertex.
    value_t &cost(vertex_t idx) {
        if constexpr (std::is_same_v<D, From>) {
            return vertices_[idx].cost_from;
        }
        else {
            return vertices_[idx].cost_to;
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
            return vertices_[idx].outgoing;
        }
        else {
            return vertices_[idx].incoming;
        }
    }

    //! The edge that was used to reach the given vertex.
    edge_t &path(vertex_t idx) {
        if constexpr (std::is_same_v<D, From>) {
            return vertices_[idx].path_from;
        }
        else {
            return vertices_[idx].path_to;
        }
    }

    //! Flag indicating whether the vertex has been visited.
    //!
    //! \note This is a integer here because it is also used as a dfs index
    //! when adding edges.
    std::conditional_t<std::is_same_v<D, From>, vertex_t, bool> &visited(vertex_t idx) {
        if constexpr(std::is_same_v<D, From>) {
            return vertices_[idx].visited_from;
        }
        else {
            return vertices_[idx].visited_to;
        }
    }

    //! Whether the vertex is relevant for propagation.
    bool &relevant(vertex_t idx) {
        if constexpr (std::is_same_v<D, From>) {
            return vertices_[idx].relevant_from;
        }
        else {
            return vertices_[idx].relevant_to;
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
            return vertices_[idx].candidate_outgoing;
        }
        else {
            return vertices_[idx].candidate_incoming;
        }
    }

    //! Incoming candidate edges that are not false.
    std::vector<vertex_t> &candidate_incoming(vertex_t idx) {
        if constexpr (std::is_same_v<D, From>) {
            return vertices_[idx].candidate_incoming;
        }
        else {
            return vertices_[idx].candidate_outgoing;
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

    //! The cost to propagate the edge.
    uint64_t &propagation_cost() {
        if constexpr (std::is_same_v<D, From>) {
            return stats_.propagate_cost_from;
        }
        else {
            return stats_.propagate_cost_to;
        }
    }

    //! Compute shortests paths starting from the given vertex.
    //!
    //! This function counts relevant nodes and stops as soon as it determines
    //! that there are no more shorted paths through relevant vertices any
    //! more. A vertex is relevant if it was reached via a shortest path
    //! containing a relevant vertex. By construction, the graph always
    //! contains a relevant vertex connected to the starting node via an edge.
    //!
    //! The function also counts the in and out degrees of visited relevant
    //! nodes. The direction with the smaller degree can then be used to try to
    //! find negative cycles.
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
                relevant_degree_out += vertices_[u_idx].degree_out;
                relevant_degree_in += vertices_[u_idx].degree_in;
            }
            for (auto &uv_idx : out(u_idx)) {
                ++propagation_cost();
                auto &uv = edges_[uv_idx];
                auto v_idx = to(uv_idx);
                // NOTE: We have to make sure that weights are positive for the
                // algorithm to run correctly. We can use the potentials to do
                // this. Explicitely using uv.from and uv.to is intended here.
                //
                auto c = cost(u_idx) + vertices_[uv.from].potential() + uv.weight - vertices_[uv.to].potential();
                assert(vertices_[uv.from].potential() + uv.weight - vertices_[uv.to].potential() >= 0);
                if (!visited(v_idx) || c < cost(v_idx)) {
                    cost(v_idx) = c;
                    if (!visited(v_idx)) {
                        // vertex v contributes an edge with a relevant source
                        if (relevant_u) {
                            ++num_relevant;
                        }
                        visited_set().push_back(to(uv_idx));
                        visited(v_idx) = true;
                        costs_heap_.push(*this, v_idx);
                    }
                    else {
                        if (relevant(from(path(v_idx)))) {
                            // vertex v no longer contributes a relevant edge
                            if (!relevant_u) {
                                --num_relevant;
                            }
                        }
                        // vertex v contributes a relevant edge now
                        else if (relevant_u) {
                            ++num_relevant;
                        }
                        costs_heap_.decrease(*this, offset(v_idx));
                    }
                    path(v_idx) = uv_idx;
                }
                else if (v_idx != source_idx && !relevant_u && c == cost(v_idx) && relevant(from(path(v_idx)))) {
                    // similar to the case above where vertex v no longer contributes a relevant edge
                    // we have to specifically handle the case that v is the source vertex which is always unrelevant
                    --num_relevant;
                    path(v_idx) = uv_idx;
                }
            }
            // removed a relevant vertex from the queue and there are no edges with relevant sources anymore in the queue
            // this condition assumes that initially there is exactly one reachable relevant vertex in the graph
            if (relevant_u && num_relevant == 0) {
                costs_heap_.clear();
                break;
            }
        }
        return {relevant_degree_out, relevant_degree_in};
    }

#ifdef CLINGODL_CROSSCHECK
    //! Count relevant incoming/outgoing edges.
    std::pair<uint32_t, uint32_t> count_relevant_() {
        uint32_t relevant_in = 0;
        uint32_t relevant_out = 0;
        for (auto &vertex : visited_set()) {
            if (relevant(vertex)) {
                for (auto &edge : candidate_incoming(vertex)) {
                    if (edge_is_enabled(edge)) {
                        ++relevant_in;
                    }
                }
                for (auto &edge : candidate_outgoing(vertex)) {
                    if (edge_is_enabled(edge)) {
                        ++relevant_out;
                    }
                }
            }
        }
        return {relevant_in, relevant_out};
    }
#endif

    //! Traverse incoming/outgoing edges and disable true edges and propagate
    //! false edges.
    bool propagate_edges(Clingo::PropagateControl &ctl, edge_t xy_idx, bool forward, bool backward) { // NOLINT
        if (!forward && !backward) {
            return true;
        }
        for (auto &vertex : visited_set()) {
            if (relevant(vertex)) {
                if (forward) {
                    auto &in = candidate_incoming(vertex);
                    in.resize(
                        std::remove_if(
                            in.begin(), in.end(),
                            [&](edge_t uv_idx) {
                                if (!edge_states_[uv_idx].enabled || propagate_edge_true_(uv_idx, xy_idx)) {
                                    remove_incoming(uv_idx);
                                    return true;
                                }
                                return false;
                            }) -
                        in.begin());
                }
                if (backward) {
                    bool ret = true;
                    auto &out = candidate_outgoing(vertex);
                    out.resize(
                        std::remove_if(
                            out.begin(), out.end(),
                            [&](edge_t uv_idx) {
                                if (!ret) {
                                    return false;
                                }
                                if (!edge_states_[uv_idx].enabled || propagate_edge_false_(ctl, uv_idx, xy_idx, ret)) {
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
        ensure_index(vertices_, std::max(edges_[i].from, edges_[i].to));
        add_candidate_edge_(i);
    }
}

template <typename T>
Graph<T>::~Graph() = default;

template <typename T>
Graph<T>::Graph(Graph &&other) noexcept = default;

template <typename T>
bool Graph<T>::empty() const {
    return vertices_.empty();
}

template <typename T>
bool Graph<T>::has_value(vertex_t idx) const {
    return idx < vertices_.size() && vertices_[idx].defined();
}

template <typename T>
T Graph<T>::get_value(vertex_t idx) const {
    assert(has_value(idx));
    return -vertices_[idx].potential();
}

template <typename T>
bool Graph<T>::edge_is_enabled(edge_t edge_idx) const {
    return edge_states_[edge_idx].enabled;
}

template <typename T>
bool Graph<T>::can_propagate() const {
    return changed_trail_.empty() || changed_trail_.back().can_propagate;
}

template <typename T>
void Graph<T>::disable_propagate() {
    changed_trail_.back().can_propagate = false;
}

template <typename T>
void Graph<T>::ensure_decision_level(level_t level, bool enable_propagate) {
    // initialize the trail
    if (changed_trail_.empty() || current_decision_level_() < level) {
        changed_trail_.push_back({level,
                                  numeric_cast<index_t>(changed_vertices_.size()),
                                  numeric_cast<index_t>(changed_edges_.size()),
                                  numeric_cast<index_t>(disabled_edges_.size()),
                                  can_propagate() && enable_propagate});
    }
}

template <typename T>
bool Graph<T>::propagate(edge_t xy_idx, Clingo::PropagateControl &ctl) {
    // The function is best understood considering the following example graph:
    //
    //   v ->* x -> y ->* u
    //   ^---------------/
    //
    // We calculate relevant shortest paths from x ->* u.
    // A shorted path is relevant if it contains edge x -> y.
    //
    // Similarly, we calculate relevant shorted paths v *<- y for the transposed graph.
    // Again, a shorted path is relevant if it contains y <- x.
    //
    // There is a negative cycle, if cost(x ->* u) + cost(v *<- y) - cost(x -> y) + cost(u -> v) is negative.
    ++stats_.edges_propagated;
    disable_edge(xy_idx);
    auto &xy = edges_[xy_idx];
    auto &x = vertices_[xy.from];
    auto &y = vertices_[xy.to];
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
    assert(std::make_pair(num_relevant_in_from, num_relevant_out_from) == static_cast<Impl<From> *>(this)->count_relevant_());
    assert(std::make_pair(num_relevant_out_to, num_relevant_in_to) == static_cast<Impl<To> *>(this)->count_relevant_());
#endif

    bool forward_from = num_relevant_in_from < num_relevant_out_to;
    bool backward_from = num_relevant_out_from < num_relevant_in_to;

    bool ret = static_cast<Impl<From> *>(this)->propagate_edges(ctl, xy_idx, forward_from, backward_from) &&
               static_cast<Impl<To> *>(this)->propagate_edges(ctl, xy_idx, !forward_from, !backward_from);

    for (auto &x : visited_from_) {
        vertices_[x].visited_from = false;
        vertices_[x].relevant_from = false;
    }
    for (auto &x : visited_to_) {
        vertices_[x].visited_to = false;
        vertices_[x].relevant_to = false;
    }
    visited_from_.clear();
    visited_to_.clear();
    return ret;
}

template <typename T>
bool Graph<T>::add_edge(Clingo::PropagateControl &ctl, edge_t uv_idx) { // NOLINT
    // This function adds an edge to the graph and returns false if the edge
    // induces a negative cycle.
    //
    // For this to work, the graph must not already have a negative cycle. The
    // function also propagates cheap to propagate edges considering vertices
    // visited during cycle detection. Negative cycles are reported via the
    // given callback function.
#ifdef CLINGODL_CROSSCHECK
    for (auto &vertex : vertices_) {
        static_cast<void>(vertex);
        assert(!vertex.visited_from);
    }
#endif
    assert(visited_from_.empty());
    assert(costs_heap_.empty());
    level_t level = current_decision_level_();
    auto &uv = edges_[uv_idx];
    // NOTE: would be more efficient if relevant would return statically false here
    //       for the compiler to make comparison cheaper
    auto &m = *static_cast<Impl<From> *>(this);

    // initialize the vertices of the edge to add
    auto &u = vertices_[uv.from];
    auto &v = vertices_[uv.to];
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
        auto &s = vertices_[s_idx];
        assert(s.visited_from);
        s.visited_from = ++dfs;
        set_potential_(s, level, s.potential() + s.cost_from);
        for (auto st_idx : s.outgoing) {
            ++stats_.propagate_cost_add;
            assert(st_idx < numeric_cast<edge_t>(edges_.size()));
            auto &st = edges_[st_idx];
            auto &t = vertices_[st.to];
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
        // check that the graph does not have a cycle
        assert(bellman_ford_(changed_edges_, uv.from).has_value());
#endif
    }
    else {
        // gather the edges in the negative cycle
        clause_.clear();
        clause_.emplace_back(-edges_[v.path_from].lit);
        auto next_idx = edges_[v.path_from].from;
#ifdef CLINGODL_CROSSCHECK
        T weight = edges_[v.path_from].weight;
#endif
        while (uv.to != next_idx) {
            auto &vertex = vertices_[next_idx];
            auto &edge = edges_[vertex.path_from];
            clause_.emplace_back(-edge.lit);
            next_idx = edge.from;
#ifdef CLINGODL_CROSSCHECK
            weight += edge.weight;
#endif
        }
#ifdef CLINGODL_CROSSCHECK
        assert(weight < 0);
#endif
        consistent = ctl.add_clause(clause_) && ctl.propagate();
    }

    if (propagate_ >= PropagationMode::Trivial && consistent) {
        if (visited_from_.empty() || propagate_ == PropagationMode::Trivial) {
            consistent = with_incoming_(ctl, uv.from, [&](vertex_t t_idx, edge_t ts_idx) {
                auto &ts = edges_[ts_idx];
                if (t_idx == uv.to && uv.weight + ts.weight < 0) {
                    clause_.emplace_back(-edges_[uv_idx].lit);
                    clause_.emplace_back(-edges_[ts_idx].lit);
                    ++stats_.false_edges_trivial;
                    return true;
                }
                return false;
            });
        }
        else if (propagate_ >= PropagationMode::Weak) {
            consistent = cheap_propagate_(ctl, uv.from, uv.from);
            if (propagate_ >= PropagationMode::WeakPlus && consistent) {
                for (auto &s_idx : visited_from_) {
                    if (!cheap_propagate_(ctl, uv.from, s_idx)) {
                        consistent = false;
                        break;
                    }
                }
            }
        }
    }
    // reset visited flags
    for (auto &x : visited_from_) {
        vertices_[x].visited_from = 0;
    }
    visited_from_.clear();
    costs_heap_.clear();

    return consistent;
}

template <typename T>
template <class F>
bool Graph<T>::with_incoming_(Clingo::PropagateControl &ctl, vertex_t s_idx, F f) {
    auto &s = vertices_[s_idx];
    auto &in = s.candidate_incoming;
    auto jt = in.begin();
    // traverse over incoming edges
    for (auto it = jt, ie = in.end(); it != ie; ++it) {
        auto &ts_idx = *it;
        auto &ts = edges_[ts_idx];
        auto t_idx = ts.from;
        // remove edges marked inactive
        if (!edge_states_[ts_idx].enabled) {
            edge_states_[ts_idx].removed_incoming = true;
            continue;
        }
        // visit the incoming vertex and edge
        clause_.clear();
        if (f(t_idx, ts_idx)) {
            edge_states_[ts_idx].removed_incoming = true;
            disable_edge(ts_idx);
            // add constraint for the negative cycle
            if (!(ctl.add_clause(clause_) && ctl.propagate())) {
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
[[nodiscard]] bool Graph<T>::cheap_propagate_(Clingo::PropagateControl &ctl, vertex_t u_idx, vertex_t s_idx) {
    // we check for the following case:
    //   u ->* s ->* t
    //         ^----/
    //           ts
    return with_incoming_(ctl, s_idx, [&](vertex_t t_idx, edge_t ts_idx) {
        auto &s = vertices_[s_idx];
        auto &t = vertices_[t_idx];
        auto &ts = edges_[ts_idx];
        if (s.visited_from < t.visited_from) {
            T weight = t.potential() - s.potential();
            if (weight + ts.weight < 0) {
                T check = 0;
                auto r_idx = t_idx;
                while (u_idx != r_idx && s_idx != r_idx) {
                    auto &r = vertices_[r_idx];
                    auto &rr = edges_[r.path_from];
                    clause_.emplace_back(-edges_[r.path_from].lit);
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
                    clause_.emplace_back(-edges_[ts_idx].lit);
                    return true;
                }
            }
        }
        return false;
    });
}

template <typename T>
void Graph<T>::backtrack() {
    auto vo = changed_trail_.back().vertex_offset;
    auto eo = changed_trail_.back().edge_offset;
    auto io = changed_trail_.back().disabled_offset;
    for (auto it = changed_vertices_.rbegin(), ie = changed_vertices_.rend() - vo; it != ie; ++it) {
        vertices_[*it].potential_stack.pop_back();
    }
    for (auto it = changed_edges_.rbegin(), ie = changed_edges_.rend() - eo; it != ie; ++it) {
        auto &edge = edges_[*it];
        vertices_[edge.from].outgoing.pop_back();
        vertices_[edge.to].incoming.pop_back();
    }
    for (auto it = disabled_edges_.begin() + io, ie = disabled_edges_.end(); it != ie; ++it) {
        add_candidate_edge_(*it);
    }
    changed_vertices_.resize(vo);
    changed_edges_.resize(eo);
    disabled_edges_.resize(io);
    changed_trail_.pop_back();
}

template <typename T>
void Graph<T>::disable_edge(edge_t uv_idx) {
    auto &uv = edges_[uv_idx];
    auto &u = vertices_[uv.from];
    auto &v = vertices_[uv.to];
    --u.degree_out;
    --v.degree_in;
    disabled_edges_.push_back(uv_idx);
    assert(edge_states_[uv_idx].enabled);
    edge_states_[uv_idx].enabled = false;
}

template <typename T>
PropagationMode Graph<T>::mode() const {
    return propagate_;
}

template <typename T>
void Graph<T>::add_candidate_edge_(edge_t uv_idx) {
    auto &uv = edges_[uv_idx];
    auto &uv_state = edge_states_[uv_idx];
    auto &u = vertices_[uv.from];
    auto &v = vertices_[uv.to];
    ++u.degree_out;
    ++v.degree_in;
    assert(!uv_state.enabled);
    uv_state.enabled = true;
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
    // The function is best understood considering the following example graph:
    //
    //   u ->* x -> y ->* v
    //   \----------------^
    //
    // Using the intermidiate costs calculated in propagate(), we get the
    // length of the shortest path u ->* v. If this path is shorter than
    // u -> v, then we disable the edge.
    auto &uv = edges_[uv_idx];
    auto &u = vertices_[uv.from];
    auto &v = vertices_[uv.to];
    assert(u.relevant_to || v.relevant_from);

    if (u.relevant_to && v.relevant_from) {
        auto &xy = edges_[xy_idx];
        auto &x = vertices_[xy.from];
        auto &y = vertices_[xy.to];

        auto cost_uy = u.cost_to + y.potential() - u.potential();
        auto cost_xv = v.cost_from + v.potential() - x.potential();
        auto cost_uv = cost_uy + cost_xv - xy.weight;
#ifdef CLINGODL_CROSSCHECK
        auto bf_costs_from_u = bellman_ford_(changed_edges_, uv.from);
        auto bf_costs_from_x = bellman_ford_(changed_edges_, xy.from);
        auto bf_cost_uy = bf_costs_from_u->find(xy.to);
        auto bf_cost_xv = bf_costs_from_x->find(uv.to);
        static_cast<void>(bf_cost_uy);
        static_cast<void>(bf_cost_xv);
        assert(bf_cost_uy != bf_costs_from_u->end());
        assert(bf_cost_xv != bf_costs_from_u->end());
        assert(bf_cost_uy->second == cost_uy);
        assert(bf_cost_xv->second == cost_xv);
#endif
        if (cost_uv <= uv.weight) {
#ifdef CLINGODL_CROSSCHECK
            // make sure that the graph does not have a negative cycle even if it contains the edge
            auto edges = changed_edges_;
            edges.emplace_back(uv_idx);
            assert(bellman_ford_(edges, uv.from).has_value());
#endif
            ++stats_.true_edges;
            disable_edge(uv_idx);
            return true;
        }
    }
    return false;
}

template <typename T>
bool Graph<T>::propagate_edge_false_(Clingo::PropagateControl &ctl, edge_t uv_idx, edge_t xy_idx, bool &ret) { // NOLINT
    // The function is best understood considering the following example graph:
    //
    //   v ->* x -> y ->* u
    //   ^----------------/
    //
    // Using the intermidiate costs calculated in propagate(), we get the
    // length of the shortest path v ->* u. If this path extended with u -> v
    // has a negative length, then the edge has to be false.
    auto &uv = edges_[uv_idx];
    auto &u = vertices_[uv.from];
    auto &v = vertices_[uv.to];
    assert(v.relevant_to || u.relevant_from);

    if (v.relevant_to && u.relevant_from) {
        auto &xy = edges_[xy_idx];
        auto &x = vertices_[xy.from];
        auto &y = vertices_[xy.to];

        auto cost_vy = v.cost_to + y.potential() - v.potential();
        auto cost_xu = u.cost_from + u.potential() - x.potential();
        auto cost_vu = cost_vy + cost_xu - xy.weight;
        if (cost_vu + uv.weight < 0) {
            ++stats_.false_edges;
            if (!ctl.assignment().is_false(uv.lit)) {
#ifdef CLINGODL_CROSSCHECK
                value_t sum = uv.weight - xy.weight;
#endif
                clause_.clear();
                clause_.push_back(-uv.lit);
                // forward
                for (auto next_edge_idx = u.path_from; next_edge_idx != invalid_edge_index;) {
                    auto &next_edge = edges_[next_edge_idx];
                    auto &next_vertex = vertices_[next_edge.from];
                    clause_.push_back(-next_edge.lit);
#ifdef CLINGODL_CROSSCHECK
                    sum += next_edge.weight;
#endif
                    next_edge_idx = next_vertex.path_from;
                }
                // backward
                for (auto next_edge_idx = v.path_to; next_edge_idx != invalid_edge_index;) {
                    auto &next_edge = edges_[next_edge_idx];
                    auto &next_vertex = vertices_[next_edge.to];
                    clause_.push_back(-next_edge.lit);
#ifdef CLINGODL_CROSSCHECK
                    sum += next_edge.weight;
#endif
                    next_edge_idx = next_vertex.path_to;
                }
#ifdef CLINGODL_CROSSCHECK
                assert(sum < 0 && sum == cost_vu + uv.weight);
#endif
                if (ret = ctl.add_clause(clause_) && ctl.propagate(); !ret) {
                    return false;
                }
            }
            disable_edge(uv_idx);
            return true;
        }
#ifdef CLINGODL_CROSSCHECK
        // make sure that the graph does not have a negative cycle even if it contains the edge
        auto edges = changed_edges_;
        edges.emplace_back(uv_idx);
        assert(bellman_ford_(edges, uv.from).has_value());
#endif
    }
    return false;
}

#ifdef CLINGODL_CROSSCHECK
template <typename T>
std::optional<std::unordered_map<vertex_t, T>> Graph<T>::bellman_ford_(std::vector<vertex_t> const &edges, vertex_t source) {
    std::unordered_map<vertex_t, T> costs;
    costs[source] = 0;
    vertex_t vertices = 0;
    for (auto &vertex : vertices_) {
        if (vertex.defined()) {
            ++vertices;
        }
    }
    for (vertex_t i = 0; i < vertices; ++i) {
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
                return std::nullopt;
            }
        }
    }
    return costs;
}
#endif

template <typename T>
void Graph<T>::set_potential_(Vertex &vtx, level_t level, T potential) {
    if (!vtx.defined() || vtx.potential_stack.back().first < level) {
        vtx.potential_stack.emplace_back(level, potential);
        changed_vertices_.emplace_back(numeric_cast<vertex_t>(&vtx - vertices_.data()));
    }
    else {
        vtx.potential_stack.back().second = potential;
    }
}

template <typename T>
level_t Graph<T>::current_decision_level_() {
    assert(!changed_trail_.empty());
    return changed_trail_.back().level;
}

template class Graph<int>;
template class Graph<double>;

} // namespace ClingoDL
