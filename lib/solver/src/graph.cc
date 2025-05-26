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
struct From {};
//!< Tag for traversals of the transposed graph.
struct To {};

//! An index different from all valid edge indices.
constexpr auto invalid_edge_index = std::numeric_limits<vertex_t>::max();

} // namespace

void ThreadStatistics::reset() { *this = ThreadStatistics{}; }

void ThreadStatistics::accu(ThreadStatistics const &x) {
    time_propagate += x.time_propagate;
    time_undo += x.time_undo;
    time_dijkstra += x.time_dijkstra;
    true_edges += x.true_edges;
    false_edges += x.false_edges;
    false_edges_trivial += x.false_edges_trivial;
    false_edges_weak += x.false_edges_weak;
    false_edges_weak_plus += x.false_edges_weak_plus;
    propagate_cost_add += x.propagate_cost_add;
    propagate_cost_from += x.propagate_cost_from;
    propagate_cost_to += x.propagate_cost_to;
    edges_added += x.edges_added;
    edges_skipped += x.edges_skipped;
    edges_propagated += x.edges_propagated;
}

//!< Thread specific information for vertices.
template <typename T> struct Graph<T>::Vertex {
    using value_t = T;
    using PotentialStack = std::vector<std::pair<vertex_t, value_t>>;

    //! Return true if the vertex has a value assigned.
    [[nodiscard]] auto defined() const -> bool { return !potential_stack.empty(); }
    //! Return the current value associated with the vertex.
    [[nodiscard]] auto potential() const -> value_t { return potential_stack.back().second; }

    VertexIndexVec outgoing;           //!< Outgoing edges from this vertex that are true.
    VertexIndexVec incoming;           //!< Incoming edges to this vertex that are true.
    VertexIndexVec candidate_incoming; //!< Edges that might become incoming edges.
    VertexIndexVec candidate_outgoing; //!< Edges that might become outgoing edges.
    PotentialStack potential_stack;    //!< Vector of pairs of level and potential.
    value_t cost_from{0};              //!< Costs for traversals of the original graph.
    value_t cost_to{0};                //!< Costs for traversals of the transposed graph.
    value_t bound_lower{0};            //!< The lower bound of a vertex.
    value_t bound_upper{0};            //!< The upper bound of a vertex.
    edge_t path_lower{0};              //!< Path pointers for determining the lower bound.
    edge_t path_upper{0};              //!< Path pointers for determining the upper bound.
    edge_t path_from{0};               //!< Path pointers for traversals of the original graph.
    edge_t path_to{0};                 //!< Path pointers for traversals of the transposed graph.
    vertex_t offset{0};                //!< Offset in the cost heap.
    vertex_t degree_out{0};            //!< Outgoing degree of candidate edges.
    vertex_t degree_in{0};             //!< Incoming degree of candidate edges.
    vertex_t visited_from{0};          //!< Either a flag to mark the vertex as visited or its depth first index.
    bool visited_to{false};            //!< A flag to mark the vertex as visited for traversals of the transposed graph.
    bool visited_lower{false};         //!< A flag to mark the vertex as visited for traversals of the original graph.
    bool visited_upper{false};         //!< A flag to mark the vertex as visited for traversals of the transposed graph.
    bool relevant_from{false};         //!< A flag to mark the vertex as visited for traversals of the original graph.
    bool relevant_to{false};           //!< A flag to mark the vertex as visited for traversals of the transposed graph.
};

//!< Thread specific information for edges.
template <typename T> struct Graph<T>::EdgeState {
    uint8_t removed_outgoing : 1; //!< Flag to mark edges as removed from the candidate_outgoing vector.
    uint8_t removed_incoming : 1; //!< Flag to mark edges as removed from the candidate_incoming vector.
    uint8_t enabled : 1;          //!< Flag to mark the edge as enabled.
};

//!< Struct holding information to backtrack a decision level.
template <typename T> struct Graph<T>::TrailEntry {
    level_t level;                //!< The corresponding decision level.
    index_t vertex_offset;        //!< Index up to which to backtrack changed vertices.
    index_t edge_offset;          //!< Index up to which to backtrack changed edges.
    index_t disabled_offset;      //!< Index up to which to backtrack inactive edges.
    index_t visited_lower_offset; //!< Index up to which to backtrack visited lower bound vertices.
    index_t visited_upper_offset; //!< Index up to which to backtrack visited upper bound vertices.
    index_t lower_value_offset;   //!< Index up to which to backtrack lower bound value changes.
    index_t upper_value_offset;   //!< Index up to which to backtrack upper bound value changes.
    bool can_propagate;           //!< Whether propagation was possible on this level.
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
template <typename T> template <typename D> struct Graph<T>::Impl : Graph {
    using index_t = HeapType::index_type;

    //! The index of the vertex in the heap vector.
    auto offset(vertex_t idx) -> index_t & { return vertices_[idx].offset; }

    //! Compare two vertices by cost/relevant.
    auto less(vertex_t a, vertex_t b) -> bool {
        auto ca = cost(a);
        auto cb = cost(b);
        return ca < cb || (ca == cb && relevant(a) < relevant(b));
    }

    //! Compare a cost with a vertex.
    auto less(value_t ca, bool ra, vertex_t b) -> bool {
        auto cb = cost(b);
        return ca < cb || (ca == cb && ra < relevant(b));
    }

    //! The cost of the vertex.
    auto cost(vertex_t idx) -> value_t & {
        if constexpr (std::is_same_v<D, From>) {
            return vertices_[idx].cost_from;
        } else {
            return vertices_[idx].cost_to;
        }
    }

    //! The bound of a vertex.
    auto bound_value(vertex_t idx) -> value_t & {
        if constexpr (std::is_same_v<D, From>) {
            return vertices_[idx].bound_lower;
        } else {
            return vertices_[idx].bound_upper;
        }
    }

    //! The bound changes during bound propagation.
    auto bound_trail() -> BoundTrailVec & {
        if constexpr (std::is_same_v<D, From>) {
            return lower_trail_;
        } else {
            return upper_trail_;
        }
    }

    //! The end point of the given edge.
    auto to(edge_t idx) -> vertex_t {
        if constexpr (std::is_same_v<D, From>) {
            return edges_[idx].to;
        } else {
            return edges_[idx].from;
        }
    }

    //! The starting point of the given edge.
    auto from(edge_t idx) -> vertex_t {
        if constexpr (std::is_same_v<D, From>) {
            return edges_[idx].from;
        } else {
            return edges_[idx].to;
        }
    }

    //! The outgoing vertices of the given vertex.
    auto out(vertex_t idx) -> std::vector<vertex_t> & {
        if constexpr (std::is_same_v<D, From>) {
            return vertices_[idx].outgoing;
        } else {
            return vertices_[idx].incoming;
        }
    }

    //! The edge that was used to reach the given vertex.
    auto path(vertex_t idx) -> edge_t & {
        if constexpr (std::is_same_v<D, From>) {
            return vertices_[idx].path_from;
        } else {
            return vertices_[idx].path_to;
        }
    }

    //! The edge that was used to reach the given vertex when computing bounds.
    auto bound_path(vertex_t idx) -> edge_t & {
        if constexpr (std::is_same_v<D, From>) {
            return vertices_[idx].path_lower;
        } else {
            return vertices_[idx].path_upper;
        }
    }

    //! Flag indicating whether the vertex has been visited.
    //!
    //! \note This is a integer here because it is also used as a dfs index
    //! when adding edges.
    auto visited(vertex_t idx) -> std::conditional_t<std::is_same_v<D, From>, vertex_t, bool> & {
        if constexpr (std::is_same_v<D, From>) {
            return vertices_[idx].visited_from;
        } else {
            return vertices_[idx].visited_to;
        }
    }

    //! Flag indicating whether the vertex has been visited during bound traversals.
    auto bound_visited(vertex_t idx) -> bool & {
        if constexpr (std::is_same_v<D, From>) {
            return vertices_[idx].visited_lower;
        } else {
            return vertices_[idx].visited_upper;
        }
    }

    //! Whether the vertex is relevant for propagation.
    auto relevant(vertex_t idx) -> bool & {
        if constexpr (std::is_same_v<D, From>) {
            return vertices_[idx].relevant_from;
        } else {
            return vertices_[idx].relevant_to;
        }
    }

    //! The set of all visited vertices.
    auto visited_set() -> std::vector<vertex_t> & {
        if constexpr (std::is_same_v<D, From>) {
            return visited_from_;
        } else {
            return visited_to_;
        }
    }

    //! The set of all visited vertices.
    auto bound_visited_set() -> std::vector<vertex_t> & {
        if constexpr (std::is_same_v<D, From>) {
            return visited_lower_;
        } else {
            return visited_upper_;
        }
    }

    //! Outgoing candidate edges that are not false.
    auto candidate_outgoing(vertex_t idx) -> std::vector<vertex_t> & {
        if constexpr (std::is_same_v<D, From>) {
            return vertices_[idx].candidate_outgoing;
        } else {
            return vertices_[idx].candidate_incoming;
        }
    }

    //! Incoming candidate edges that are not false.
    auto candidate_incoming(vertex_t idx) -> std::vector<vertex_t> & {
        if constexpr (std::is_same_v<D, From>) {
            return vertices_[idx].candidate_incoming;
        } else {
            return vertices_[idx].candidate_outgoing;
        }
    }

    //! Mark an incoming edge as removed.
    void remove_incoming(edge_t idx) {
        if constexpr (std::is_same_v<D, From>) {
            edge_states_[idx].removed_incoming = true;
        } else {
            edge_states_[idx].removed_outgoing = true;
        }
    }

    //! Mark an outgoing edge as removed.
    void remove_outgoing(edge_t idx) {
        if constexpr (std::is_same_v<D, From>) {
            edge_states_[idx].removed_outgoing = true;
        } else {
            edge_states_[idx].removed_incoming = true;
        }
    }

    //! The cost to propagate the edge.
    auto propagation_cost() -> uint64_t & {
        if constexpr (std::is_same_v<D, From>) {
            return stats_.propagate_cost_from;
        } else {
            return stats_.propagate_cost_to;
        }
    }

#ifdef CLINGODL_CROSSCHECK
    //! Compute shortests paths using the Bellman Ford algorithm.
    std::optional<std::unordered_map<vertex_t, value_t>> bellman_ford_(std::vector<edge_t> const &edges,
                                                                       vertex_t source) {
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
                auto u_cost = costs.find(from(uv_idx));
                if (u_cost != costs.end()) {
                    auto u_idx = to(uv_idx);
                    auto v_cost = costs.find(u_idx);
                    auto dist = u_cost->second + uv.weight;
                    if (v_cost == costs.end()) {
                        costs[u_idx] = dist;
                    } else if (dist < v_cost->second) {
                        v_cost->second = dist;
                    }
                }
            }
        }
        for (auto const &uv_idx : edges) {
            auto &uv = edges_[uv_idx];
            auto u_cost = costs.find(from(uv_idx));
            if (u_cost != costs.end()) {
                auto v_cost = costs.find(to(uv_idx));
                auto dist = u_cost->second + uv.weight;
                if (dist < v_cost->second) {
                    return std::nullopt;
                }
            }
        }
        return costs;
    }
#endif

    //! Compute shortests paths starting from the given vertex.
    //!
    //! This function counts relevant nodes and stops as soon as it determines
    //! that there are no more shortest paths through relevant vertices any
    //! more. A vertex is relevant if it was reached via a shortest path
    //! containing a relevant vertex. By construction, the graph always
    //! contains a relevant vertex connected to the starting vertex via an
    //! edge.
    //!
    //! The function also counts the in and out degrees of visited relevant
    //! nodes. The direction with the smaller degree can then be used to try to
    //! find negative cycles.
    std::pair<uint32_t, uint32_t> dijkstra_full(vertex_t source_idx, edge_t relevant_idx) { // NOLINT
        // Note: Initially there is exactly one relevant vertex in the graph,
        // which is guaranteed to be reached by edge relevant_idx. Once this
        // vertex enters the queue, the count corresponds to the number of
        // relevant vertices in the queue. The algorithm can stop, once the
        // count reaches zero.
        uint32_t num_relevant = 1;
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
            bool relevant_u = relevant(u_idx);
            if (relevant_u) {
                --num_relevant;
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
                auto c = cost(u_idx) + vertices_[uv.from].potential() + uv.weight - vertices_[uv.to].potential();
                assert(vertices_[uv.from].potential() + uv.weight - vertices_[uv.to].potential() >= 0);
                if (!visited(v_idx) || less(c, relevant_u, v_idx)) {
                    bool relevant_v = relevant(v_idx);
                    cost(v_idx) = c;
                    if (!visited(v_idx)) {
                        // vertex v became relevant
                        if (relevant_u && !relevant_v) {
                            ++num_relevant;
                        }
                        visited_set().push_back(to(uv_idx));
                        visited(v_idx) = true;
                        relevant(v_idx) = relevant_u || uv_idx == relevant_idx;
                        costs_heap_.push(*this, v_idx);
                    } else {
                        // vertex v became or lost relevance
                        if (relevant_u != relevant_v) {
                            relevant(v_idx) = relevant_u;
                            num_relevant += relevant_u ? 1 : -1;
                        }
                        costs_heap_.decrease(*this, v_idx);
                    }
                    path(v_idx) = uv_idx;
                }
            }
            if (num_relevant == 0) {
                costs_heap_.clear();
                break;
            }
        }
        return {relevant_degree_out, relevant_degree_in};
    }

    //! Compute shortest paths starting from the zero node.
    //!
    //! This function is conceptually similar to dijkstra_full() but computes
    //! paths starting from a zero node. Like this, it can reuse existing paths
    //! from previous calls.
    void dijkstra_bounds(edge_t uv_idx, vertex_t zero_idx) { // NOLINT
        auto u_idx = from(uv_idx);
        // the zero node is reached with zero cost as soon as it is added to the graph
        if (u_idx == zero_idx && !bound_visited(u_idx)) {
            assert(!vertices_[u_idx].visited_lower && !vertices_[u_idx].visited_upper);
            vertices_[u_idx].bound_lower = vertices_[u_idx].bound_upper = 0;
            vertices_[u_idx].path_lower = vertices_[u_idx].path_upper = invalid_edge_index;
            vertices_[u_idx].visited_lower = vertices_[u_idx].visited_upper = true;
            visited_lower_.push_back(u_idx);
            visited_upper_.push_back(u_idx);
        }
        if (!bound_visited(u_idx)) {
            return;
        }
        assert(visited_set().empty() && costs_heap_.empty());
        auto v_idx = to(uv_idx);
        auto bound_uv = bound_value(u_idx) + edges_[uv_idx].weight;
        if (!bound_visited(v_idx) || bound_uv < bound_value(v_idx)) {
            cost(v_idx) = 0;
            costs_heap_.push(*this, v_idx);
            visited(v_idx) = true;
            visited_set().push_back(v_idx);
            if (bound_visited(v_idx)) {
                bound_trail().emplace_back(std::make_tuple(v_idx, bound_path(v_idx), bound_value(v_idx)));
            }
            bound_value(v_idx) = bound_uv;
            bound_path(v_idx) = uv_idx;
        }
        std::vector<value_t> value_set;
        while (!costs_heap_.empty()) {
            auto s_idx = costs_heap_.pop(*this);
            for (auto &st_idx : out(s_idx)) {
                // Note: maybe use a separate counter
                ++propagation_cost();
                auto &st = edges_[st_idx];
                auto t_idx = to(st_idx);
                // NOTE: We have to make sure that weights are positive for the
                // algorithm to run correctly. We can use the potentials to do
                // this. Explicitely using st.from and st.to is intended here.
                auto c = cost(s_idx) + vertices_[st.from].potential() + st.weight - vertices_[st.to].potential();
                auto bound_st = bound_value(s_idx) + st.weight;
                assert(vertices_[st.from].potential() + st.weight - vertices_[st.to].potential() >= 0);
                if ((!visited(t_idx) || less(c, false, t_idx)) &&
                    (!bound_visited(t_idx) || bound_st < bound_value(t_idx))) {
                    cost(t_idx) = c;
                    if (!visited(t_idx)) {
                        visited(t_idx) = true;
                        visited_set().push_back(t_idx);
                        if (bound_visited(t_idx)) {
                            bound_trail().emplace_back(std::make_tuple(t_idx, bound_path(t_idx), bound_value(t_idx)));
                        }
                        costs_heap_.push(*this, t_idx);
                    } else {
                        costs_heap_.decrease(*this, t_idx);
                    }
                    bound_path(t_idx) = st_idx;
                    bound_value(t_idx) = bound_st;
                }
            }
        }
        for (auto &vertex_idx : visited_set()) {
            if (!bound_visited(vertex_idx)) {
                bound_visited(vertex_idx) = true;
                bound_visited_set().push_back(vertex_idx);
            }
            visited(vertex_idx) = false;
        }
#ifdef CLINGODL_CROSSCHECK
        if (!visited_set().empty()) {
            auto costs = bellman_ford_(changed_edges_, zero_idx);
            assert(costs.has_value());
            assert(costs->size() == bound_visited_set().size());
            for (auto [vertex_idx, value] : *costs) {
                static_cast<void>(vertex_idx);
                static_cast<void>(value);
                assert(bound_visited(vertex_idx));
                assert(bound_value(vertex_idx) == value);
            }
        }
#endif
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

    //! Traverse incoming/outgoing edges and disables true edges and propagates
    //! false edges.
    template <bool full>
    bool propagate_edges(Clingo::PropagateControl &ctl, edge_t xy_idx, bool forward, bool backward) { // NOLINT
        if (!forward && !backward) {
            return true;
        }
        for (auto &vertex : visited_set()) {
            if (!full || relevant(vertex)) {
                if (forward) {
                    auto &in = candidate_incoming(vertex);
                    in.resize(std::remove_if(in.begin(), in.end(),
                                             [&](edge_t uv_idx) {
                                                 if (!edge_states_[uv_idx].enabled ||
                                                     propagate_edge_true_<full>(uv_idx, xy_idx)) {
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
                    out.resize(std::remove_if(out.begin(), out.end(),
                                              [&](edge_t uv_idx) {
                                                  if (!ret) {
                                                      return false;
                                                  }
                                                  if (!edge_states_[uv_idx].enabled ||
                                                      propagate_edge_false_<full>(ctl, uv_idx, xy_idx, ret)) {
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
    : edges_{edges}, stats_{stats}, propagate_{propagate} {
    edge_states_.resize(edges_.size(), {1, 1, 0});
    for (edge_t i = 0; i < numeric_cast<edge_t>(edges_.size()); ++i) {
        ensure_index(vertices_, std::max(edges_[i].from, edges_[i].to));
        add_candidate_edge_(i);
    }
}

template <typename T> Graph<T>::~Graph() = default;

template <typename T> Graph<T>::Graph(Graph &&other) noexcept = default;

template <typename T> auto Graph<T>::empty() const -> bool { return vertices_.empty(); }

template <typename T> auto Graph<T>::has_value(vertex_t idx) const -> bool {
    return idx < vertices_.size() && vertices_[idx].defined();
}

template <typename T> auto Graph<T>::get_value(vertex_t idx) const -> T {
    assert(has_value(idx));
    return -vertices_[idx].potential();
}

template <typename T> auto Graph<T>::edge_is_enabled(edge_t edge_idx) const -> bool {
    return edge_states_[edge_idx].enabled;
}

template <typename T> auto Graph<T>::edge_is_negative(edge_t edge_idx) const -> bool {
    auto const &st = edges_[edge_idx];
    auto const &u = vertices_[st.from];
    auto const &v = vertices_[st.to];
    auto up = u.defined() ? u.potential() : 0;
    auto vp = v.defined() ? v.potential() : 0;
    return up + st.weight - vp < 0;
}

template <typename T> auto Graph<T>::can_propagate() const -> bool {
    return changed_trail_.empty() || changed_trail_.back().can_propagate;
}

template <typename T> void Graph<T>::disable_propagate() { changed_trail_.back().can_propagate = false; }

template <typename T> void Graph<T>::ensure_decision_level(level_t level, bool enable_propagate) {
    // initialize the trail
    if (changed_trail_.empty() || current_decision_level_() < level) {
        changed_trail_.push_back(
            {level, numeric_cast<index_t>(changed_vertices_.size()), numeric_cast<index_t>(changed_edges_.size()),
             numeric_cast<index_t>(disabled_edges_.size()), numeric_cast<index_t>(visited_lower_.size()),
             numeric_cast<index_t>(visited_upper_.size()), numeric_cast<index_t>(lower_trail_.size()),
             numeric_cast<index_t>(upper_trail_.size()), can_propagate() && enable_propagate});
    }
}

template <typename T> auto Graph<T>::add_edge(Clingo::PropagateControl &ctl, edge_t uv_idx, vertex_t zero_idx) -> bool {
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

    // cycle check and simple propagation
    auto consistent = check_cycle_(ctl, uv_idx) && propagate_simple_(ctl, uv_idx);

    // reset visited flags
    for (auto &x : visited_from_) {
        vertices_[x].visited_from = 0;
    }
    visited_from_.clear();
    costs_heap_.clear();

    // propagate cycles through zero node
    // (using equality here intentional because full and zero vertex propagation should be exclusive)
    if (mode() == PropagationMode::Zero) {
        consistent = consistent && propagate_zero_(ctl, uv_idx, zero_idx);
    }

    // full propagation
    if (can_propagate()) {
        consistent = consistent && propagate_full_(ctl, uv_idx);
    }

    return consistent;
}

template <typename T>
bool Graph<T>::propagate_zero_(Clingo::PropagateControl &ctl, edge_t uv_idx, vertex_t zero_idx) { // NOLINT
    ++stats_.edges_propagated;
    disable_edge(uv_idx);

    Timer t{stats_.time_dijkstra};
    static_cast<Impl<From> *>(this)->dijkstra_bounds(uv_idx, zero_idx);
    static_cast<Impl<To> *>(this)->dijkstra_bounds(uv_idx, zero_idx);
    t.stop();

    bool ret = static_cast<Impl<From> *>(this)->template propagate_edges<false>(ctl, 0, true, true) &&
               static_cast<Impl<To> *>(this)->template propagate_edges<false>(ctl, 0, true, true);

    visited_from_.clear();
    visited_to_.clear();

#ifdef CLINGODL_CROSSCHECK
    if (ret) {
        auto ass = ctl.assignment();
        auto cost_from_zero = static_cast<Impl<From> *>(this)->bellman_ford_(changed_edges_, zero_idx);
        auto cost_to_zero = static_cast<Impl<To> *>(this)->bellman_ford_(changed_edges_, zero_idx);
        if (cost_from_zero && cost_to_zero) {
            for (auto [x_idx, cost_x] : *cost_from_zero) {
                for (auto &xy : edges_) {
                    if (xy.from != x_idx) {
                        continue;
                    }
                    auto it_y = cost_to_zero->find(xy.to);
                    if (it_y == cost_to_zero->end()) {
                        continue;
                    }
                    auto cost_y = it_y->second;
                    static_cast<void>(cost_y);
                    static_cast<void>(ass);
                    assert(cost_x + cost_y + xy.weight >= 0 || ass.is_false(xy.lit));
                }
            }
        }
    }

#endif

    return ret;
}

template <typename T> bool Graph<T>::check_cycle_(Clingo::PropagateControl &ctl, edge_t uv_idx) { // NOLINT
    // NOTE: would be more efficient if relevant would return statically false here
    //       for the compiler to make comparison cheaper

    auto &m = *static_cast<Impl<From> *>(this);
    level_t level = current_decision_level_();
    auto &uv = edges_[uv_idx];

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
    } else {
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
                } else {
                    costs_heap_.decrease(m, st.to);
                }
            }
        }
    }

    // there no negative cycle
    if (!u.visited_from) {
        // add the edge to the graph
        u.outgoing.emplace_back(uv_idx);
        v.incoming.emplace_back(uv_idx);
        changed_edges_.emplace_back(uv_idx);
#ifdef CLINGODL_CROSSCHECK
        // check that the graph does not have a cycle
        assert(m.bellman_ford_(changed_edges_, uv.from).has_value());
#endif
        return true;
    }

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
    return ctl.add_clause(clause_) && ctl.propagate();
}

template <typename T> bool Graph<T>::propagate_simple_(Clingo::PropagateControl &ctl, edge_t uv_idx) { // NOLINT
    if (propagate_ >= PropagationMode::Trivial) {
        auto &uv = edges_[uv_idx];
        if (visited_from_.empty() || propagate_ == PropagationMode::Trivial) {
            return with_incoming_(ctl, uv.from, [&](vertex_t t_idx, edge_t ts_idx) {
                auto &ts = edges_[ts_idx];
                if (t_idx == uv.to && uv.weight + ts.weight < 0 && !ctl.assignment().is_false(ts.lit)) {
                    clause_.emplace_back(-edges_[uv_idx].lit);
                    clause_.emplace_back(-edges_[ts_idx].lit);
                    ++stats_.false_edges_trivial;
                    return true;
                }
                return false;
            });
        }
        if (propagate_ >= PropagationMode::Weak) {
            if (!cheap_propagate_(ctl, uv.from, uv.from)) {
                return false;
            }
            if (propagate_ >= PropagationMode::WeakPlus) {
                for (auto &s_idx : visited_from_) {
                    if (!cheap_propagate_(ctl, uv.from, s_idx)) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

template <typename T> auto Graph<T>::propagate_full_(Clingo::PropagateControl &ctl, edge_t xy_idx) -> bool {
    // The function is best understood considering the following example graph:
    //
    //   v ->* x -> y ->* u
    //   ^---------------/
    //
    // We calculate relevant shortest paths from x ->* u.
    // A shortesd path is relevant if it contains edge x -> y.
    //
    // Similarly, we calculate relevant shortest paths v *<- y for the transposed graph.
    // Again, a shortest path is relevant if it contains y <- x.
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
        std::tie(num_relevant_out_from, num_relevant_in_from) =
            static_cast<Impl<From> *>(this)->dijkstra_full(xy.from, xy_idx);
        std::tie(num_relevant_out_to, num_relevant_in_to) = static_cast<Impl<To> *>(this)->dijkstra_full(xy.to, xy_idx);
    }
#ifdef CLINGODL_CROSSCHECK
    // check if the counts of relevant incoming/outgoing vertices are correct
    assert(std::make_pair(num_relevant_in_from, num_relevant_out_from) ==
           static_cast<Impl<From> *>(this)->count_relevant_());
    assert(std::make_pair(num_relevant_out_to, num_relevant_in_to) == static_cast<Impl<To> *>(this)->count_relevant_());
#endif

    bool forward_from = num_relevant_in_from < num_relevant_out_to;
    bool backward_from = num_relevant_out_from < num_relevant_in_to;

    bool ret =
        static_cast<Impl<From> *>(this)->template propagate_edges<true>(ctl, xy_idx, forward_from, backward_from) &&
        static_cast<Impl<To> *>(this)->template propagate_edges<true>(ctl, xy_idx, !forward_from, !backward_from);

    for (auto &x : visited_from_) {
        vertices_[x].visited_from = 0;
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
template <class F>
auto Graph<T>::with_incoming_(Clingo::PropagateControl &ctl, vertex_t s_idx, F f) -> bool {
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
                in.erase(jt, it + 1);
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
[[nodiscard]] auto Graph<T>::cheap_propagate_(Clingo::PropagateControl &ctl, vertex_t u_idx, vertex_t s_idx) -> bool {
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
            if (weight + ts.weight < 0 && !ctl.assignment().is_false(ts.lit)) {
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
                    } else {
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

template <typename T> void Graph<T>::backtrack() {
    auto entry = changed_trail_.back();
    for (auto it = changed_vertices_.rbegin(), ie = changed_vertices_.rend() - entry.vertex_offset; it != ie; ++it) {
        vertices_[*it].potential_stack.pop_back();
    }
    for (auto it = changed_edges_.rbegin(), ie = changed_edges_.rend() - entry.edge_offset; it != ie; ++it) {
        auto &edge = edges_[*it];
        vertices_[edge.from].outgoing.pop_back();
        vertices_[edge.to].incoming.pop_back();
    }
    for (auto it = disabled_edges_.begin() + entry.disabled_offset, ie = disabled_edges_.end(); it != ie; ++it) {
        add_candidate_edge_(*it);
    }
    for (auto it = visited_lower_.begin() + entry.visited_lower_offset, ie = visited_lower_.end(); it != ie; ++it) {
        vertices_[*it].visited_lower = false;
    }
    for (auto it = visited_upper_.begin() + entry.visited_upper_offset, ie = visited_upper_.end(); it != ie; ++it) {
        vertices_[*it].visited_upper = false;
    }
    for (auto it = lower_trail_.rbegin(), ie = lower_trail_.rend() - entry.lower_value_offset; it != ie; ++it) {
        auto &[vertex_idx, edge_idx, value] = *it;
        vertices_[vertex_idx].path_lower = edge_idx;
        vertices_[vertex_idx].bound_lower = value;
    }
    for (auto it = upper_trail_.rbegin(), ie = upper_trail_.rend() - entry.upper_value_offset; it != ie; ++it) {
        auto &[vertex_idx, edge_idx, value] = *it;
        vertices_[vertex_idx].path_upper = edge_idx;
        vertices_[vertex_idx].bound_upper = value;
    }
    changed_vertices_.resize(entry.vertex_offset);
    changed_edges_.resize(entry.edge_offset);
    disabled_edges_.resize(entry.disabled_offset);
    visited_lower_.resize(entry.visited_lower_offset);
    visited_upper_.resize(entry.visited_upper_offset);
    lower_trail_.resize(entry.lower_value_offset);
    upper_trail_.resize(entry.upper_value_offset);
    changed_trail_.pop_back();
}

template <typename T> void Graph<T>::disable_edge(edge_t uv_idx) {
    auto &uv = edges_[uv_idx];
    auto &u = vertices_[uv.from];
    auto &v = vertices_[uv.to];
    --u.degree_out;
    --v.degree_in;
    disabled_edges_.push_back(uv_idx);
    assert(edge_states_[uv_idx].enabled);
    edge_states_[uv_idx].enabled = false;
}

template <typename T> auto Graph<T>::mode() const -> PropagationMode { return propagate_; }

template <typename T> void Graph<T>::add_candidate_edge_(edge_t uv_idx) {
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

template <typename T> template <bool full> bool Graph<T>::propagate_edge_true_(edge_t uv_idx, edge_t xy_idx) { // NOLINT
    // The function is best understood considering the following example graph:
    //
    //   u ->* x -> y ->* v
    //   \----------------^
    //
    // Using the intermidiate costs calculated in propagate(), we get the
    // length of the shortest path u ->* v. If this path is shorter than
    // u -> v, then we disable the edge.
    //
    // The case for propagation through the zero node is a bit simpler because
    // we do not have to consider edge x -> y but the zero node instead:
    //
    //   u ->* 0 ->* v
    //   \-----------^
    //
    auto &uv = edges_[uv_idx];
    auto &u = vertices_[uv.from];
    auto &v = vertices_[uv.to];
    bool relevant = false;
    if constexpr (full) {
        assert(u.relevant_to || v.relevant_from);
        relevant = u.relevant_to && v.relevant_from;
    } else {
        assert(u.visited_upper || v.visited_lower);
        relevant = u.visited_upper && v.visited_lower;
    }

    if (relevant) {
        auto &xy = edges_[xy_idx];
        auto &x = vertices_[xy.from];
        auto &y = vertices_[xy.to];

        value_t cost_uv{0};
        if constexpr (full) {
            auto cost_uy = u.cost_to + y.potential() - u.potential();
            auto cost_xv = v.cost_from + v.potential() - x.potential();
            cost_uv = cost_uy + cost_xv - xy.weight;
#ifdef CLINGODL_CROSSCHECK
            auto &m = *static_cast<Impl<From> *>(this);
            auto bf_costs_from_u = m.bellman_ford_(changed_edges_, uv.from);
            auto bf_costs_from_x = m.bellman_ford_(changed_edges_, xy.from);
            auto bf_cost_uy = bf_costs_from_u->find(xy.to);
            auto bf_cost_xv = bf_costs_from_x->find(uv.to);
            static_cast<void>(bf_cost_uy);
            static_cast<void>(bf_cost_xv);
            assert(bf_cost_uy != bf_costs_from_u->end());
            assert(bf_cost_xv != bf_costs_from_u->end());
            assert(bf_cost_uy->second == cost_uy);
            assert(bf_cost_xv->second == cost_xv);
#endif
        } else {
            cost_uv = u.bound_upper + v.bound_lower;
        }
        if (cost_uv <= uv.weight) {
#ifdef CLINGODL_CROSSCHECK
            // make sure that the graph does not have a negative cycle even if it contains the edge
            if constexpr (full) {
                auto edges = changed_edges_;
                edges.emplace_back(uv_idx);
                auto &m = *static_cast<Impl<From> *>(this);
                static_cast<void>(m);
                assert(m.bellman_ford_(edges, uv.from).has_value());
            }
#endif
            ++stats_.true_edges;
            disable_edge(uv_idx);
            return true;
        }
    }
    return false;
}

template <typename T>
template <bool full>
bool Graph<T>::propagate_edge_false_(Clingo::PropagateControl &ctl, edge_t uv_idx, edge_t xy_idx, bool &ret) { // NOLINT
    // The function is best understood considering the following example graph:
    //
    //   v ->* x -> y ->* u
    //   ^----------------/
    //
    // Using the intermediate costs calculated in propagate(), we get the
    // length of the shortest path v ->* u. If this path extended with u -> v
    // has a negative length, then the edge has to be false.
    //
    // The case for propagation through the zero node is a bit simpler because
    // we do not have to consider edge x -> y but the zero node instead:
    //
    //   v ->* 0 ->* u
    //   ^-----------/
    //
    auto &uv = edges_[uv_idx];
    auto &u = vertices_[uv.from];
    auto &v = vertices_[uv.to];
    bool relevant = true;
    if constexpr (full) {
        assert(v.relevant_to || u.relevant_from);
        relevant = v.relevant_to && u.relevant_from;
    } else {
        assert(v.visited_upper || u.visited_lower);
        relevant = v.visited_upper && u.visited_lower;
    }

    if (relevant) {
        auto &xy = edges_[xy_idx];
        auto &x = vertices_[xy.from];
        auto &y = vertices_[xy.to];

        value_t cost_vu{0};
        if constexpr (full) {
            auto cost_vy = v.cost_to + y.potential() - v.potential();
            auto cost_xu = u.cost_from + u.potential() - x.potential();
            cost_vu = cost_vy + cost_xu - xy.weight;
        } else {
            cost_vu = v.bound_upper + u.bound_lower;
        }

        if (cost_vu + uv.weight < 0) {
            ++stats_.false_edges;
            if (!ctl.assignment().is_false(uv.lit)) {
#ifdef CLINGODL_CROSSCHECK
                value_t sum = uv.weight;
                if constexpr (full) {
                    sum -= xy.weight;
                }
#endif
                clause_.clear();
                clause_.push_back(-uv.lit);
                // forward
                for (auto next_edge_idx = full ? u.path_from : u.path_lower; next_edge_idx != invalid_edge_index;) {
                    auto &next_edge = edges_[next_edge_idx];
                    auto &next_vertex = vertices_[next_edge.from];
                    clause_.push_back(-next_edge.lit);
#ifdef CLINGODL_CROSSCHECK
                    sum += next_edge.weight;
#endif
                    next_edge_idx = full ? next_vertex.path_from : next_vertex.path_lower;
                }
                // backward
                for (auto next_edge_idx = full ? v.path_to : v.path_upper; next_edge_idx != invalid_edge_index;) {
                    auto &next_edge = edges_[next_edge_idx];
                    auto &next_vertex = vertices_[next_edge.to];
                    clause_.push_back(-next_edge.lit);
#ifdef CLINGODL_CROSSCHECK
                    sum += next_edge.weight;
#endif
                    next_edge_idx = full ? next_vertex.path_to : next_vertex.path_upper;
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
        if constexpr (full) {
            auto edges = changed_edges_;
            edges.emplace_back(uv_idx);
            auto &m = *static_cast<Impl<From> *>(this);
            static_cast<void>(m);
            assert(m.bellman_ford_(edges, uv.from).has_value());
        }
#endif
    }
    return false;
}

template <typename T> void Graph<T>::set_potential_(Vertex &vtx, level_t level, T potential) {
    if (!vtx.defined() || vtx.potential_stack.back().first < level) {
        vtx.potential_stack.emplace_back(level, potential);
        changed_vertices_.emplace_back(numeric_cast<vertex_t>(&vtx - vertices_.data()));
    } else {
        vtx.potential_stack.back().second = potential;
    }
}

template <typename T> auto Graph<T>::current_decision_level_() -> level_t {
    assert(!changed_trail_.empty());
    return changed_trail_.back().level;
}

template class Graph<int>;
template class Graph<double>;

} // namespace ClingoDL
