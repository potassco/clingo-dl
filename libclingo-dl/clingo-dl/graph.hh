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

#ifndef CLINGODL_GRAPH_HH
#define CLINGODL_GRAPH_HH

#include <clingo.hh>
#include <clingo-dl/theory.hh>
#include <clingo-dl/config.hh>
#include <clingo-dl/util.hh>

namespace ClingoDL {

//! Struct to represent an edge in the difference logic graph.
template <typename T>
struct Edge {
    using value_t = T;
    vertex_t from;         //!< Start vertex index of the edge.
    vertex_t to;           //!< End vertex index of the edge.
    value_t weight;        //!< Weight of the edge.
    Clingo::literal_t lit; //!< Solver literal associated with the edge.
};

//! Struct to capture per thread statistics of the DL propagator/graph.
struct ThreadStatistics {
    //! Reset the object to its initial state.
    void reset();
    //! Add values from given statistics object.
    void accu(ThreadStatistics const &x);

    Duration time_propagate{0};        //!< Total time spend in propagate callback.
    Duration time_undo{0};             //!< Total time spend in the undo function.
    Duration time_dijkstra{0};         //!< Total runtime of the dijkstra algorithm during full propagation.
    uint64_t true_edges{0};            //!< The number of edges that have been made true.
    uint64_t false_edges{0};           //!< The number of edges that have been made false.
    uint64_t false_edges_trivial{0};   //!< Edges made false by inverse propagation.
    uint64_t false_edges_weak{0};      //!< Edges made false by weak propagation.
    uint64_t false_edges_weak_plus{0}; //!< Edges made false by extended weak propagation.
    uint64_t propagate_cost_add{0};    //!< Cost estimate for adding edges.
    uint64_t propagate_cost_from{0};   //!< Cost estimate for propagation on original graph.
    uint64_t propagate_cost_to{0};     //!< Cost estimate for propagation on transposed graph.
    uint64_t edges_added{0};           //!< The number of edges added to the graph.
    uint64_t edges_skipped{0};         //!< Edges added for which no potential updates were necessary.
    uint64_t edges_propagated{0};      //!< The number of edges for which full propagation has been called.
};

//! This class stores the DL graph and provides methods to incrementally add
//! edges, perform propagation, and access the assignment to vertices.
//!
//! The class is meant to be used by the DLPropagator.
template <typename T>
class Graph {
private:
    struct Vertex;
    struct EdgeState;
    struct TrailEntry;
    template <typename D>
    struct Impl;
    using HeapType = Heap<4>;
    using index_t = HeapType::index_type;     //!< Type for indexing (32bit unsigned).
    using value_t = T;                        //!< The value type (integral or floating point).
    using Edge = ClingoDL::Edge<value_t>;     //!< The data structure for edges.
    using EdgeVec = std::vector<Edge>;        //!< A vector of edges.
    using VertexVec = std::vector<Vertex>;    //!< A vector of vertices.
    using TrailVec = std::vector<TrailEntry>; //!< The trail to backtrack per decision level changes.
public:
    Graph(ThreadStatistics &stats, EdgeVec const &edges, PropagationMode propagate);
    Graph(Graph const &other) = delete;
    Graph(Graph &&other) noexcept;
    Graph &operator=(Graph const &other) = delete;
    Graph &operator=(Graph &&other) = delete;
    ~Graph();
    //! Return no edges have been added to the graph yet.
    [[nodiscard]] bool empty() const;
    //! Return true if a potential has been assigned to the vertex.
    [[nodiscard]] bool has_value(vertex_t idx) const;
    //! Return true the (inverted) potential of a vertex.
    [[nodiscard]] value_t get_value(vertex_t idx) const;
    //! Return true if the edge is enabled.
    //!
    //! Disabled edges are ignored during full propagation.
    [[nodiscard]] bool edge_is_enabled(edge_t edge_idx) const;
    //! This disables an edge.
    //!
    //! Any edge whose literal became true or false is disabled. See also
    //! edge_is_enabled(). The state is disabled when backtracking.
    void disable_edge(edge_t uv_idx);
    //! Return true if the propagator can run the full propagate.
    //!
    //! Propagation can be disabled on a decision level. It will stay disabled
    //! on higher decision level and reset during backtracking.
    [[nodiscard]] bool can_propagate() const;
    //! Disable propagation on the current and higher decision levels.
    //!
    //! See also can_propagate().
    void disable_propagate();
    //! Make sure that the trail contains the given decision level.
    void ensure_decision_level(level_t level, bool enable_propagate);
    //! Add an edge to the graph and return false if the edge induces a negative cycle.
    //!
    //! This function assumes that the graph is not conflicting.
    [[nodiscard]] bool add_edge(edge_t uv_idx, std::function<bool(std::vector<edge_t>)> f);
    //! Fully propagates the graph after adding the given edge.
    //!
    //! Afterward any of the remaining edges can be added to the graph without
    //! causing a conflict. The function assumes that the graph was fully
    //! propagated before the edge was added.
    [[nodiscard]] bool propagate(edge_t xy_idx, Clingo::PropagateControl &ctl);
    //! Backtracks the last decision level established with ensure_decision_level().
    void backtrack();
    //! Return the configured propagation mode.
    [[nodiscard]] PropagationMode mode() const;

private:
    //! Traverse the incoming edges of a vertex.
    //!
    //! Edges that became inactive will be removed during the traversal.
    //! Callback f is called for each visited edge and the edge is removed if
    //! the callback returns false. Furthermore, function p is called in this
    //! case, which can be used to add constraints to the graph. If adding a
    //! constraint causes a conflict indicated by the return value of the
    //! callback, the traversal will stop.
    template <class P, class F>
    [[nodiscard]] bool with_incoming_(vertex_t s_idx, P p, F f);
    //! If s has been reached from u, we can use the current potentials to
    //! detect some conflicts involving incoming edges of s.
    template <class F>
    [[nodiscard]] bool cheap_propagate_(vertex_t u_idx, vertex_t s_idx, F f);
    //! Helper to add candidate edges initially and during backtracking.
    void add_candidate_edge_(edge_t uv_idx);
    //! Disable edge u -> v, if there is a shorter path u ->* v.
    [[nodiscard]] bool propagate_edge_true_(edge_t uv_idx, edge_t xy_idx);
    [[nodiscard]] bool propagate_edge_false_(Clingo::PropagateControl &ctl, edge_t uv_idx, edge_t xy_idx, bool &ret);
#ifdef CLINGODL_CROSSCHECK
    std::optional<std::unordered_map<int, value_t>> bellman_ford_(std::vector<edge_t> const &edges, int source);
#endif
    //! Sets the potential of a vertex and ensures that it can be backtracked.
    void set_potential_(Vertex &vtx, level_t level, value_t potential);
    //! Returns the current decision level.
    [[nodiscard]] level_t current_decision_level_();

    HeapType costs_heap_;
    std::vector<vertex_t> visited_from_;
    std::vector<vertex_t> visited_to_;
    EdgeVec const &edges_;
    VertexVec vertices_;
    std::vector<vertex_t> changed_vertices_;
    std::vector<edge_t> changed_edges_;
    TrailVec changed_trail_;
    std::vector<edge_t> inactive_edges_;
    std::vector<EdgeState> edge_states_;
    std::vector<edge_t> neg_cycle_;
    ThreadStatistics &stats_;
    PropagationMode propagate_;
};

} // namespace ClingoDL

#endif // CLINGODL_PROPAGATOR_HH
