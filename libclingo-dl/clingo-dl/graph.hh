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

struct ThreadStatistics {
    void reset() {
        time_propagate        = std::chrono::steady_clock::duration::zero();
        time_undo             = std::chrono::steady_clock::duration::zero();
        time_dijkstra         = std::chrono::steady_clock::duration::zero();
        true_edges            = 0;
        false_edges           = 0;
        false_edges_trivial   = 0;
        false_edges_weak      = 0;
        false_edges_weak_plus = 0;
        propagate_cost_add    = 0;
        propagate_cost_from   = 0;
        propagate_cost_to     = 0;
        edges_added           = 0;
        edges_skipped         = 0;
        edges_propagated      = 0;
    }
    void accu(ThreadStatistics const &x) {
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
    Duration time_propagate = Duration{0};
    Duration time_undo = Duration{0};
    Duration time_dijkstra = Duration{0};
    uint64_t true_edges{0};
    uint64_t false_edges{0};
    uint64_t false_edges_trivial{0};
    uint64_t false_edges_weak{0};
    uint64_t false_edges_weak_plus{0};
    uint64_t propagate_cost_add{0};
    uint64_t propagate_cost_from{0};
    uint64_t propagate_cost_to{0};
    uint64_t edges_added{0};
    uint64_t edges_skipped{0};
    uint64_t edges_propagated{0};
};

template <typename T>
class Graph {
private:
    struct Vertex;
    struct EdgeState;
    template <class D>
    struct Impl;
    using value_t = T;                     //!< The value type (integral or floating point).
    using Edge = ClingoDL::Edge<value_t>;  //!< The data structure for edges.
    using EdgeVec = std::vector<Edge>;     //!< A vector of edges.
    using VertexVec = std::vector<Vertex>; //!< A vector of vertices.

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
    //! Traverse the incoming edges of a node.
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
    void add_candidate_edge_(edge_t uv_idx);
    [[nodiscard]] bool propagate_edge_true_(edge_t uv_idx, edge_t xy_idx);
    [[nodiscard]] bool propagate_edge_false_(Clingo::PropagateControl &ctl, edge_t uv_idx, edge_t xy_idx, bool &ret);
#ifdef CLINGODL_CROSSCHECK
    std::unordered_map<int, value_t> bellman_ford_(std::vector<edge_t> const &edges, int source);
#endif
    void set_potential_(Vertex &node, level_t level, value_t potential);
    [[nodiscard]] level_t current_decision_level_();

    Heap<4> costs_heap_;
    std::vector<vertex_t> visited_from_;
    std::vector<vertex_t> visited_to_;
    EdgeVec const &edges_;
    VertexVec nodes_;
    std::vector<vertex_t> changed_nodes_;
    std::vector<edge_t> changed_edges_;
    std::vector<std::tuple<level_t, uint32_t, uint32_t, uint32_t, bool>> changed_trail_;
    std::vector<edge_t> inactive_edges_;
    std::vector<EdgeState> edge_states_;
    std::vector<edge_t> neg_cycle_;
    ThreadStatistics &stats_;
    PropagationMode propagate_;
};

} // namespace ClingoDL

#endif // CLINGODL_PROPAGATOR_HH
