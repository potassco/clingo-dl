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

template <class T, class P>
struct HeapFromM {
    using index_t = Heap<0>::index_type;
    using value_t = T;
    index_t &offset(vertex_t idx) { return static_cast<P *>(this)->nodes_[idx].offset; }
    value_t &cost(vertex_t idx) { return static_cast<P *>(this)->nodes_[idx].cost_from; }
    vertex_t to(edge_t idx) { return static_cast<P *>(this)->edges_[idx].to; }
    vertex_t from(edge_t idx) { return static_cast<P *>(this)->edges_[idx].from; }
    std::vector<vertex_t> &out(vertex_t idx) { return static_cast<P *>(this)->nodes_[idx].outgoing; }
    edge_t &path(vertex_t idx) { return static_cast<P *>(this)->nodes_[idx].path_from; }
    vertex_t &visited(vertex_t idx) { return static_cast<P *>(this)->nodes_[idx].visited_from; }
    bool &relevant(vertex_t idx) { return static_cast<P *>(this)->nodes_[idx].relevant_from; }
    std::vector<vertex_t> &visited_set() { return static_cast<P *>(this)->visited_from_; }
    std::vector<vertex_t> &candidate_outgoing(vertex_t idx) { return static_cast<P *>(this)->nodes_[idx].candidate_outgoing; }
    std::vector<vertex_t> &candidate_incoming(vertex_t idx) { return static_cast<P *>(this)->nodes_[idx].candidate_incoming; }
    void remove_incoming(edge_t idx) { static_cast<P *>(this)->edge_states_[idx].removed_incoming = true; }
    void remove_outgoing(edge_t idx) { static_cast<P *>(this)->edge_states_[idx].removed_outgoing = true; }
    uint64_t &propagation_cost() { return static_cast<P *>(this)->stats_.propagate_cost_from; }
};

template <class T, class P>
struct HeapToM {
    using index_t = Heap<0>::index_type;
    using value_t = T;
    index_t &offset(vertex_t idx) { return static_cast<P *>(this)->nodes_[idx].offset; }
    value_t &cost(vertex_t idx) { return static_cast<P *>(this)->nodes_[idx].cost_to; }
    vertex_t to(edge_t idx) { return static_cast<P *>(this)->edges_[idx].from; }
    vertex_t from(edge_t idx) { return static_cast<P *>(this)->edges_[idx].to; }
    std::vector<vertex_t> &out(vertex_t idx) { return static_cast<P *>(this)->nodes_[idx].incoming; }
    edge_t &path(vertex_t idx) { return static_cast<P *>(this)->nodes_[idx].path_to; }
    bool &visited(vertex_t idx) { return static_cast<P *>(this)->nodes_[idx].visited_to; }
    bool &relevant(vertex_t idx) { return static_cast<P *>(this)->nodes_[idx].relevant_to; }
    std::vector<vertex_t> &visited_set() { return static_cast<P *>(this)->visited_to_; }
    std::vector<vertex_t> &candidate_outgoing(vertex_t idx) { return static_cast<P *>(this)->nodes_[idx].candidate_incoming; }
    std::vector<vertex_t> &candidate_incoming(vertex_t idx) { return static_cast<P *>(this)->nodes_[idx].candidate_outgoing; }
    void remove_incoming(edge_t idx) { static_cast<P *>(this)->edge_states_[idx].removed_outgoing = true; }
    void remove_outgoing(edge_t idx) { static_cast<P *>(this)->edge_states_[idx].removed_incoming = true; }
    uint64_t &propagation_cost() { return static_cast<P *>(this)->stats_.propagate_cost_to; }
};

//! Struct to represent an edge in the difference logic graph.
template <typename T>
struct Edge {
    using value_t = T;
    vertex_t from;         //!< Start vertex index of the edge.
    vertex_t to;           //!< End vertex index of the edge.
    value_t weight;        //!< Weight of the edge.
    Clingo::literal_t lit; //!< Solver literal associated with the edge.
};

template <typename T>
struct Vertex {
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
    PotentialStack potential_stack; //!< Vector of pairs of level and potential.
    value_t cost_from{0};
    value_t cost_to{0};
    uint32_t offset{0};
    edge_t path_from{0};
    edge_t path_to{0};
    vertex_t degree_out{0};
    vertex_t degree_in{0};
    vertex_t visited_from{0};
    bool relevant_from{false};
    bool relevant_to{false};
    bool visited_to{false};
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

struct EdgeState {
    uint8_t removed_outgoing : 1;
    uint8_t removed_incoming : 1;
    uint8_t active : 1;
};

template <typename T>
class Graph : private HeapToM<T, Graph<T>>, private HeapFromM<T, Graph<T>> {
    using value_t = T;
    using Edge = ClingoDL::Edge<value_t>;
    using EdgeVec = std::vector<Edge>;
    using Vertex = ClingoDL::Vertex<value_t>;
    using VertexVec = std::vector<Vertex>;
    using HTM = HeapToM<value_t, Graph<value_t>>;
    using HFM = HeapFromM<value_t, Graph<value_t>>;
    friend HTM;
    friend HFM;

public:
    Graph(ThreadStatistics &stats, EdgeVec const &edges, PropagationMode propagate);
    [[nodiscard]] bool empty() const;
    [[nodiscard]] bool valid_node(vertex_t idx) const;
    [[nodiscard]] bool node_value_defined(vertex_t idx) const;
    [[nodiscard]] bool has_value(vertex_t idx) const;
    [[nodiscard]] value_t node_value(vertex_t idx) const;
    [[nodiscard]] bool edge_is_active(edge_t edge_idx) const;
    [[nodiscard]] bool can_propagate() const;
    void disable_propagate();
    void ensure_decision_level(level_t level, bool enable_propagate);
    [[nodiscard]] bool add_edge(edge_t uv_idx, std::function<bool(std::vector<edge_t>)> f);
    [[nodiscard]] bool propagate(edge_t xy_idx, Clingo::PropagateControl &ctl);
    void backtrack();
    void remove_candidate_edge(edge_t uv_idx);
    [[nodiscard]] PropagationMode mode() const;

private:
    template <class P, class F>
    [[nodiscard]] bool with_incoming_(vertex_t s_idx, P p, F f);
    template <class F>
    [[nodiscard]] bool cheap_propagate_(vertex_t u_idx, vertex_t s_idx, F f);
    void add_candidate_edge_(edge_t uv_idx);
    [[nodiscard]] bool propagate_edge_true_(edge_t uv_idx, edge_t xy_idx);
    [[nodiscard]] bool propagate_edge_false_(Clingo::PropagateControl &ctl, edge_t uv_idx, edge_t xy_idx, bool &ret);
    template <class M>
    [[nodiscard]] bool propagate_edges_(M &m, Clingo::PropagateControl &ctl, edge_t xy_idx, bool forward, bool backward);
    template <class M>
    [[nodiscard]] std::pair<vertex_t, vertex_t> dijkstra_(vertex_t source_idx, std::vector<vertex_t> &visited_set, M &m);
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
