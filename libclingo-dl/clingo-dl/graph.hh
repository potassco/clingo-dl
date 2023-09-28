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

#include <clingo-dl/config.hh>
#include <clingo-dl/theory.hh>
#include <clingo-dl/util.hh>
#include <clingo.hh>

namespace ClingoDL {

//! Struct to represent an edge in the difference logic graph.
template <typename T> struct Edge {
    using value_t = T;     //!< The value type (integral or floating point).
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
template <typename T> class Graph {
  private:
    struct Vertex;
    struct EdgeState;
    struct TrailEntry;
    template <typename D> struct Impl;
    using HeapType = Heap<4>;                 //!< The heap type used in (dijkstra-like) algorithms.
    using index_t = HeapType::index_type;     //!< Type for indexing (32bit unsigned).
    using value_t = T;                        //!< The value type (integral or floating point).
    using Edge = ClingoDL::Edge<value_t>;     //!< The data structure for edges.
    using EdgeVec = std::vector<Edge>;        //!< A vector of edges.
    using VertexVec = std::vector<Vertex>;    //!< A vector of vertices.
    using TrailVec = std::vector<TrailEntry>; //!< The trail to backtrack per decision level changes.
    //! Trail of bound changes involving a vertex, path, and value.
    using BoundTrailVec = std::vector<std::tuple<vertex_t, edge_t, value_t>>;

  public:
    Graph(ThreadStatistics &stats, EdgeVec const &edges, PropagationMode propagate);
    Graph(Graph const &other) = delete;
    Graph(Graph &&other) noexcept;
    auto operator=(Graph const &other) -> Graph & = delete;
    auto operator=(Graph &&other) -> Graph & = delete;
    ~Graph();
    //! Return no edges have been added to the graph yet.
    [[nodiscard]] auto empty() const -> bool;
    //! Return true if a potential has been assigned to the vertex.
    [[nodiscard]] auto has_value(vertex_t idx) const -> bool;
    //! Return true the (inverted) potential of a vertex.
    [[nodiscard]] auto get_value(vertex_t idx) const -> value_t;
    //! Return true if the edge is enabled.
    //!
    //! Disabled edges are ignored during full propagation.
    [[nodiscard]] auto edge_is_enabled(edge_t edge_idx) const -> bool;
    //! Check if the edge is negative w.r.t. the current assginment of potentials.
    [[nodiscard]] auto edge_is_negative(edge_t edge_idx) const -> bool;
    //! This disables an edge.
    //!
    //! Any edge whose literal became true or false is disabled. See also
    //! edge_is_enabled(). The state is disabled when backtracking.
    void disable_edge(edge_t uv_idx);
    //! Return true if the propagator can run the full propagate.
    //!
    //! Propagation can be disabled on a decision level. It will stay disabled
    //! on higher decision level and reset during backtracking.
    [[nodiscard]] auto can_propagate() const -> bool;
    //! Disable propagation on the current and higher decision levels.
    //!
    //! See also can_propagate().
    void disable_propagate();
    //! Make sure that the trail contains the given decision level.
    void ensure_decision_level(level_t level, bool enable_propagate);
    //! Add an edge to the graph and return false if the edge induces a negative cycle.
    //!
    //! This function assumes that the graph is not conflicting.
    [[nodiscard]] auto add_edge(Clingo::PropagateControl &ctl, edge_t uv_idx, vertex_t zero_idx) -> bool;
    //! Backtracks the last decision level established with ensure_decision_level().
    void backtrack();
    //! Return the configured propagation mode.
    [[nodiscard]] auto mode() const -> PropagationMode;

  private:
    //! Checks if there is a cycle after adding the given edge.
    //!
    //! The function does not clean up information about found shortest paths.
    //!
    //! \note Has to be called during add_edge().
    auto check_cycle_(Clingo::PropagateControl &ctl, edge_t uv_idx) -> bool;
    //! Perform configured simple propagations after adding the given edge.
    //!
    //! \note Has to be called during add_edge() and uses temporary state
    //! established during check_cycle_().
    auto propagate_simple_(Clingo::PropagateControl &ctl, edge_t uv_idx) -> bool;
    //! Propagates edges to avoid cycles through the zero node.
    //!
    //! \note Has to be called during add_edge().
    auto propagate_zero_(Clingo::PropagateControl &ctl, edge_t uv_idx, vertex_t zero_idx) -> bool;
    //! Fully propagates the graph after adding the given edge.
    //!
    //! Afterward any of the remaining edges can be added to the graph without
    //! causing a conflict. The function assumes that the graph was fully
    //! propagated before the edge was added.
    //!
    //! \note Has to be called during add_edge().
    [[nodiscard]] auto propagate_full_(Clingo::PropagateControl &ctl, edge_t xy_idx) -> bool;
    //! Traverse the incoming edges of a vertex.
    //!
    //! Edges that were disabled will be removed during the traversal. The
    //! clallback is called for each visited edge and the edge is removed if
    //! the callback returns true. Furthermore, the functions assumes that the
    //! callback provides a clause in this case. If adding the clause causes a
    //! conflict, the traversal will stop and the function returns false.
    template <class F> [[nodiscard]] auto with_incoming_(Clingo::PropagateControl &ctl, vertex_t s_idx, F f) -> bool;
    //! If s has been reached from u, we can use the current potentials to
    //! detect some conflicts involving incoming edges of s.
    [[nodiscard]] auto cheap_propagate_(Clingo::PropagateControl &ctl, vertex_t u_idx, vertex_t s_idx) -> bool;
    //! Helper to add candidate edges initially and during backtracking.
    void add_candidate_edge_(edge_t uv_idx);
    //! Disable edge u -> v, if there is a shorter path u ->* v.
    template <bool full> [[nodiscard]] auto propagate_edge_true_(edge_t uv_idx, edge_t xy_idx) -> bool;
    //! Make edge u -> v false, if there is a negative cycle u ->* v -> u.
    template <bool full>
    [[nodiscard]] auto propagate_edge_false_(Clingo::PropagateControl &ctl, edge_t uv_idx, edge_t xy_idx, bool &ret)
        -> bool;
    //! Sets the potential of a vertex and ensures that it can be backtracked.
    void set_potential_(Vertex &vtx, level_t level, value_t potential);
    //! Returns the current decision level.
    [[nodiscard]] auto current_decision_level_() -> level_t;

    HeapType costs_heap_;                    //!< Heap used throught various algorithms.
    std::vector<vertex_t> visited_from_;     //!< Vertices visited during traversal of the original graph.
    std::vector<vertex_t> visited_to_;       //!< Vertices visited during traversal of the transposed graph.
    std::vector<vertex_t> visited_lower_;    //!< Visited vertices with lower bounds.
    std::vector<vertex_t> visited_upper_;    //!< Visited vertices with upper bounds.
    BoundTrailVec lower_trail_;              //!< Trail of changes to lower bounds.
    BoundTrailVec upper_trail_;              //!< Trail of changes to upper bounds.
    EdgeVec const &edges_;                   //!< Reference to the edges.
    VertexVec vertices_;                     //!< Vertex specific information.
    std::vector<vertex_t> changed_vertices_; //!< Vertices changed in chronological order.
    std::vector<edge_t> changed_edges_;      //!< Edges changed in chronological order.
    TrailVec changed_trail_;                 //!< Vector recording backtrackable per decision level information.
    std::vector<edge_t> disabled_edges_;     //!< Disabled edges in chronological order.
    std::vector<EdgeState> edge_states_;     //!< Thread-specific per edge information.
    std::vector<edge_t> neg_cycle_;          //!< Vector for storing negative cycles.
    std::vector<literal_t> clause_;          //!< Vector for storing clauses.
    ThreadStatistics &stats_;                //!< Per thread statistics.
    PropagationMode const propagate_;        //!< The propagation mode.
};

} // namespace ClingoDL

#endif // CLINGODL_PROPAGATOR_HH
