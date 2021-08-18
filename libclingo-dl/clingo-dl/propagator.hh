// {{{ MIT License
//
// Copyright 2018 Roland Kaminski, Philipp Wanko, Max Ostrowski
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

#ifndef CLINGODL_PROPAGATOR_HH
#define CLINGODL_PROPAGATOR_HH

#include <clingo.hh>
#include <clingo-dl/util.hh>
#include <clingo-dl/parsing.hh>
#include <clingo-dl/graph.hh>

#include <unordered_map>

#define CLINGODL_CHECKSOLUTION

namespace ClingoDL {

struct Statistics {
    void reset() {
        time_init  = std::chrono::steady_clock::duration::zero();
        ccs = 0;
        mutexes = 0;
        edges = 0;
        variables = 0;
        for (auto& i : dl_stats) {
            i.reset();
        }
    }
    void accu(Statistics const &x) {
        time_init += x.time_init;
        ccs = x.ccs;
        mutexes += x.mutexes;
        edges = x.edges;
        variables = x.variables;
        if (dl_stats.size() < x.dl_stats.size()) {
            dl_stats.resize(x.dl_stats.size());
        }
        auto it = x.dl_stats.begin();
        for (auto &y : dl_stats) {
            y.accu(*it++);
        }
    }
    Duration time_init = Duration{0};
    uint64_t ccs{0};
    uint64_t mutexes{0};
    uint64_t edges{0};
    uint64_t variables{0};
    std::vector<ThreadStatistics> dl_stats;
};

//! A propagator for difference constraints.
template <typename T>
class DLPropagator : public Clingo::Propagator {
private:
    struct NodeInfo;
    struct ThreadState;
    struct FactState;
    using value_t = T;
    using Edge = ClingoDL::Edge<value_t>;
    using Graph = ClingoDL::Graph<value_t>;
    using CoVarVec = ClingoDL::CoVarVec<value_t>; //!< Vector of coefficients and variables.

public:
    DLPropagator(Statistics &stats, PropagatorConfig conf);
    DLPropagator(DLPropagator const &other) = delete;
    DLPropagator(DLPropagator &&other) noexcept;
    DLPropagator &operator=(DLPropagator const &other) = delete;
    DLPropagator &operator=(DLPropagator &&other) = delete;
    ~DLPropagator() override;

    //! Get the number of vertices in the graph.
    size_t num_vertices() const;
    //! Get the symbol associated with a vertex index.
    Clingo::Symbol symbol(size_t index) const;
    //! Lookup the index of a vertex.
    uint32_t lookup(clingo_symbol_t symbol);
    //! Check if teh given vertex has a lower bound in the given thread.
    bool has_lower_bound(uint32_t thread_id, size_t index) const;
    //! Get the lower bound of a vertex in the given thread.
    T lower_bound(uint32_t thread_id, size_t index) const;
    //! Extend the model with vertex assignments.
    void extend_model(Clingo::Model &model);

    // propagator interface
    //! Initialize the propagator.
    void init(Clingo::PropagateInit &init) override;
    //! Propagate edges.
    void propagate(Clingo::PropagateControl &ctl, Clingo::LiteralSpan changes) override;
    //! Undo propgated edges.
    void undo(Clingo::PropagateControl const &ctl, Clingo::LiteralSpan changes) noexcept override;
    //! Check a propagation fixed point.
    void check(Clingo::PropagateControl &ctl) override;

private:
    // initialization functions
    //! Map a symbol to an integer.
    int map_vertex_(Clingo::Symbol v);
    //! Add constraints in the theory data.
    void add_constraints_(Clingo::PropagateInit &init);
    //! Normalize constraints to individual edges over `<=`.
    bool normalize_constraint_(Clingo::PropagateInit &init, int literal, CoVarVec const &elements, char const *op, T rhs, bool strict);
    //! Add up to two edges for the given constraint if the has at most 2 variables and suitable coefficients.
    bool add_edges_(Clingo::PropagateInit& init, int literal, CoVarVec const &covec, T rhs, bool strict);
    //! Add up to two edges for a constraint.
    void add_edges_(Clingo::PropagateInit& init, int u_id, int v_id, T weight, int lit, bool strict);
    //! Add (up to one) edge for a constraint.
    void add_edge_(Clingo::PropagateInit &init, int u_id, int v_id, T weight, int lit);
    //! Reset connected components (to be recalculated with the next call to cc_calculate_).
    void cc_reset_();
    //! Mark a vertex in a component as visited.
    bool cc_visited_(vertex_t vertex) const;
    //! Check if the given vertex is a zero node.
    bool is_zero_(vertex_t vertex) const;
    //! Calculate the connected components.
    void cc_calculate_(std::unordered_multimap<int, int> &outgoing, std::unordered_multimap<int, int> &incoming);
    //! Calculate mutually exclusive edges.
    void calculate_mutexes_(Clingo::PropagateInit &init, int edge_start, std::unordered_multimap<int, int> &outgoing);
    //! Finalize initialization by initializing thread states for propagation.
    void initialize_states_(Clingo::PropagateInit &init);

    // propagation functions
    //! Disable all edges associated with the given literal in the thread state.
    void disable_edge_by_lit(ThreadState &state, Clingo::literal_t lit);
    //! Compute the current potential of a vertex in the given graph.
    int get_potential_(Graph const &graph, vertex_t vertex);
    //! Compute the current cost of a vertex in the given graph.
    int cost_(Graph const &graph, Edge const &edge);
    //! Compute a weight to sort vertices before propagation.
    int cost_(SortMode mode, Graph const &graph, int i);
    //! Sort vertices in the propagation queue according to the given mode.
    void sort_edges(SortMode mode, ThreadState &state);
    //! Propagate the given literals.
    void do_propagate(Clingo::PropagateControl &ctl, Clingo::LiteralSpan changes);

    using LitEdgeMap = std::unordered_multimap<Clingo::literal_t, edge_t>;
    using SymVertexMap = std::unordered_map<Clingo::Symbol, vertex_t>;

    std::vector<ThreadState> states_;      //!< Thread specific state.
    std::vector<FactState> facts_;         //!< Thread specific state for fact propagation.
    LitEdgeMap lit_to_edges_;              //!< Map from literals to edges.
    LitEdgeMap false_lit_to_edges_;        //!< Map from (inverted) literals to edges. (looks unnecessary)
    std::vector<Edge> edges_;              //!< Vector holding all edges.
    std::vector<Clingo::Symbol> vert_map_; //!< Mapping from vertex indices to symbols. (can go to node info)
    SymVertexMap vert_map_inv_;            //!< Mapping from symbols to vertex indices.
    std::vector<NodeInfo> node_info_;      //!< Node specific information.
    std::vector<vertex_t> zero_nodes_;     //!< Vector holding all zero nodes.
    Statistics &stats_;                    //!< Statistics for the propagator.
    PropagatorConfig conf_;                //!< Configuration of the propagator.
};

} // namespace ClingoDL

#endif // CLINGODL_PROPAGATOR_HH
