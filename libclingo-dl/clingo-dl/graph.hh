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

#ifndef CLINGODL_GRAPH_HH
#define CLINGODL_GRAPH_HH

#include <clingo.hh>
#include <clingo-dl/theory.hh>
#include <clingo-dl/config.hh>
#include <clingo-dl/util.hh>

namespace ClingoDL {

template <class T, class P>
struct HeapFromM {
    int &offset(int idx) { return static_cast<P *>(this)->nodes_[idx].offset; }
    T &cost(int idx) { return static_cast<P *>(this)->nodes_[idx].cost_from; }
    int to(int idx) { return static_cast<P *>(this)->edges_[idx].to; }
    int from(int idx) { return static_cast<P *>(this)->edges_[idx].from; }
    std::vector<int> &out(int idx) { return static_cast<P *>(this)->nodes_[idx].outgoing; }
    int &path(int idx) { return static_cast<P *>(this)->nodes_[idx].path_from; }
    int &visited(int idx) { return static_cast<P *>(this)->nodes_[idx].visited_from; }
    bool &relevant(int idx) { return static_cast<P *>(this)->nodes_[idx].relevant_from; }
    std::vector<int> &visited_set() { return static_cast<P *>(this)->visited_from_; }
    std::vector<int> &candidate_outgoing(int idx) { return static_cast<P *>(this)->nodes_[idx].candidate_outgoing; }
    std::vector<int> &candidate_incoming(int idx) { return static_cast<P *>(this)->nodes_[idx].candidate_incoming; }
    void remove_incoming(int idx) { static_cast<P *>(this)->edge_states_[idx].removed_incoming = true; }
    void remove_outgoing(int idx) { static_cast<P *>(this)->edge_states_[idx].removed_outgoing = true; }
    uint64_t &propagation_cost() {return static_cast<P *>(this)->stats_.propagate_cost_from; }
};

template <class T, class P>
struct HeapToM {
    int &offset(int idx) { return static_cast<P *>(this)->nodes_[idx].offset; }
    T &cost(int idx) { return static_cast<P *>(this)->nodes_[idx].cost_to; }
    int to(int idx) { return static_cast<P *>(this)->edges_[idx].from; }
    int from(int idx) { return static_cast<P *>(this)->edges_[idx].to; }
    std::vector<int> &out(int idx) { return static_cast<P *>(this)->nodes_[idx].incoming; }
    int &path(int idx) { return static_cast<P *>(this)->nodes_[idx].path_to; }
    bool &visited(int idx) { return static_cast<P *>(this)->nodes_[idx].visited_to; }
    bool &relevant(int idx) { return static_cast<P *>(this)->nodes_[idx].relevant_to; }
    std::vector<int> &visited_set() { return static_cast<P *>(this)->visited_to_; }
    std::vector<int> &candidate_outgoing(int idx) { return static_cast<P *>(this)->nodes_[idx].candidate_incoming; }
    std::vector<int> &candidate_incoming(int idx) { return static_cast<P *>(this)->nodes_[idx].candidate_outgoing; }
    void remove_incoming(int idx) { static_cast<P *>(this)->edge_states_[idx].removed_outgoing = true; }
    void remove_outgoing(int idx) { static_cast<P *>(this)->edge_states_[idx].removed_incoming = true; }
    uint64_t &propagation_cost() {return static_cast<P *>(this)->stats_.propagate_cost_to; }
};

template <typename T>
struct Edge {
    vertex_t from;
    vertex_t to;
    T weight;
    Clingo::literal_t lit;
};

template <typename T>
struct DifferenceLogicNode {
    [[nodiscard]] bool defined() const { return !potential_stack.empty(); }
    [[nodiscard]] T potential() const { return potential_stack.back().second; }

    std::vector<int> outgoing;
    std::vector<int> incoming;
    std::vector<int> candidate_incoming;
    std::vector<int> candidate_outgoing;
    std::vector<std::pair<int, T>> potential_stack; // [(level,potential)]
    T cost_from = 0;
    T cost_to = 0;
    int offset = 0;
    int path_from = 0;
    int path_to = 0;
    int degree_out = 0;
    int degree_in = 0;
    int visited_from = 0;
    bool relevant_from = false;
    bool relevant_to = false;
    bool visited_to = false;
};

struct DLStats {
    void reset() {
        time_propagate = std::chrono::steady_clock::duration::zero();
        time_undo      = std::chrono::steady_clock::duration::zero();
        time_dijkstra  = std::chrono::steady_clock::duration::zero();
        true_edges            = 0;
        false_edges           = 0;
        false_edges_trivial   = 0;
        false_edges_weak      = 0;
        false_edges_weak_plus = 0;
        propagate_cost_add  = 0;
        propagate_cost_from = 0;
        propagate_cost_to   = 0;
        edges_added      = 0;
        edges_skipped    = 0;
        edges_propagated = 0;
    }
    void accu(DLStats const &x) {
        time_propagate+= x.time_propagate;
        time_undo     += x.time_undo;
        time_dijkstra += x.time_dijkstra;
        true_edges    += x.true_edges;
        false_edges   += x.false_edges;
        false_edges_trivial  += x.false_edges_trivial;
        false_edges_weak     += x.false_edges_weak;
        false_edges_weak_plus+= x.false_edges_weak_plus;
        propagate_cost_add += x.propagate_cost_add;
        propagate_cost_from+= x.propagate_cost_from;
        propagate_cost_to  += x.propagate_cost_to;
        edges_added      += x.edges_added;
        edges_skipped    += x.edges_skipped;
        edges_propagated += x.edges_propagated;
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
class DifferenceLogicGraph : private HeapToM<T, DifferenceLogicGraph<T>>, private HeapFromM<T, DifferenceLogicGraph<T>> {
    using HTM = HeapToM<T, DifferenceLogicGraph<T>>;
    using HFM = HeapFromM<T, DifferenceLogicGraph<T>>;
    friend HTM;
    friend HFM;

public:
    DifferenceLogicGraph(DLStats &stats, const std::vector<Edge<T>> &edges, PropagationMode propagate);
    [[nodiscard]] bool empty() const;
    [[nodiscard]] bool valid_node(int idx) const;
    [[nodiscard]] int node_value_defined(int idx) const;
    [[nodiscard]] bool has_value(int idx) const;
    [[nodiscard]] T node_value(int idx) const;
    [[nodiscard]] bool edge_is_active(int edge_idx) const;
    [[nodiscard]] bool can_propagate() const;
    void disable_propagate();
    void ensure_decision_level(int level, bool enable_propagate);
    [[nodiscard]] bool add_edge(int uv_idx, std::function<bool(std::vector<int>)> f);
    [[nodiscard]] bool propagate(int xy_idx, Clingo::PropagateControl &ctl);
    void backtrack();
    void remove_candidate_edge(int uv_idx);
    [[nodiscard]] PropagationMode mode() const;

private:
    template <class P, class F>
    [[nodiscard]] bool with_incoming(int s_idx, P p, F f);

    template <class F>
    [[nodiscard]] bool cheap_propagate(int u_idx, int s_idx, F f);

    void add_candidate_edge(int uv_idx);
    [[nodiscard]] bool propagate_edge_true(int uv_idx, int xy_idx);
    [[nodiscard]] bool propagate_edge_false(Clingo::PropagateControl &ctl, int uv_idx, int xy_idx, bool &ret);
    template <class M>
    [[nodiscard]] bool propagate_edges(M &m, Clingo::PropagateControl &ctl, int xy_idx, bool forward, bool backward);
    template <class M>
    [[nodiscard]] std::pair<int, int> dijkstra(int source_idx, std::vector<int> &visited_set, M &m);
#ifdef CROSSCHECK
    [[nodiscard]] std::unordered_map<int, T> bellman_ford(std::vector<int> const &edges, int source);
#endif
    void set_potential(DifferenceLogicNode<T> &node, int level, T potential);
    [[nodiscard]] int current_decision_level_();

    Heap<4> costs_heap_;
    std::vector<int> visited_from_;
    std::vector<int> visited_to_;
    std::vector<Edge<T>> const &edges_;
    std::vector<DifferenceLogicNode<T>> nodes_;
    std::vector<int> changed_nodes_;
    std::vector<int> changed_edges_;
    std::vector<std::tuple<int, int, int, int, bool>> changed_trail_;
    std::vector<int> inactive_edges_;
    std::vector<EdgeState> edge_states_;
    std::vector<int> neg_cycle_;
    DLStats &stats_;
    PropagationMode propagate_;
};

} // namespace ClingoDL

#endif // CLINGODL_PROPAGATOR_HH
