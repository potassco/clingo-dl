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

#ifndef CLINGODL_PROPAGATOR_HH
#define CLINGODL_PROPAGATOR_HH

#include <clingo.hh>
#include <clingo-dl/util.hh>
#include <clingo-dl/parsing.hh>
#include <clingo-dl/graph.hh>

#include <unordered_map>

namespace ClingoDL {

//! Struct with statistics for the DLPropagator.
struct Statistics {
    //! Reset the statistics object.
    void reset();
    //! Accumulate values from another statistics object.
    void accu(Statistics const &x);

    //! The duration of the initialization step.
    Duration time_init = Duration{0};
    //! The number of connected componenets in the difference logic graph.
    uint64_t ccs{0};
    //! The number of mutexes found during initialisation.
    uint64_t mutexes{0};
    //! The number of edges in the difference logic graph.
    uint64_t edges{0};
    //! The number of distinct variables appearing in the difference constraints.
    uint64_t variables{0};
    //! Per thread statistics.
    std::vector<ThreadStatistics> thread_statistics;
};

//! A propagator for difference constraints.
template <typename T>
class DLPropagator : public Clingo::Propagator {
public:
    using value_t = T;

    // construction
    DLPropagator(Statistics &stats, PropagatorConfig conf);
    DLPropagator(DLPropagator const &other) = delete;
    DLPropagator(DLPropagator &&other) = delete;
    DLPropagator &operator=(DLPropagator const &other) = delete;
    DLPropagator &operator=(DLPropagator &&other) = delete;
    ~DLPropagator() override;

    //! Get the number of vertices in the graph.
    [[nodiscard]] vertex_t num_vertices() const;
    //! Get the symbol associated with a vertex index.
    [[nodiscard]] Clingo::Symbol symbol(vertex_t index) const;
    //! Lookup the index of a vertex.
    vertex_t lookup(Clingo::Symbol symbol);
    //! Check if the given vertex has a lower bound in the given thread.
    [[nodiscard]] bool has_lower_bound(Clingo::id_t thread_id, vertex_t index) const;
    //! Get the lower bound of a vertex in the given thread.
    [[nodiscard]] value_t lower_bound(Clingo::id_t thread_id, vertex_t index) const;
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
    struct VertexInfo;
    struct ThreadState;
    struct FactState;
    using Edge = ClingoDL::Edge<value_t>;
    using Graph = ClingoDL::Graph<value_t>;
    //! Vector of coefficients and variables.
    using CoVarVec = ClingoDL::CoVarVec<value_t>;
    //! Mapping from literals to edge indices.
    using LitEdgeMap = std::unordered_multimap<literal_t, edge_t>;
    using SymVertexMap = std::unordered_map<Clingo::Symbol, vertex_t>;
    using VertexInfoVec = std::vector<VertexInfo>;
    using ThreadStateVec = std::vector<ThreadState>;
    using FactStateVec = std::vector<FactState>;
    using AdjacencyMap = std::unordered_multimap<vertex_t, vertex_t>;

    // initialization functions
    //! Map a symbol to an integer.
    vertex_t map_vertex_(Clingo::Symbol symbol);
    //! Add constraints in the theory data.
    [[nodiscard]] bool add_constraints_(Clingo::PropagateInit &init);
    //! Normalize constraints to individual edges over `<=`.
    [[nodiscard]] bool normalize_constraint_(Clingo::PropagateInit &init, literal_t literal, CoVarVec const &elements, char const *op, value_t rhs, bool strict);
    //! Add up to two edges for the given constraint if the has at most 2 variables and suitable coefficients.
    [[nodiscard]] bool add_edges_(Clingo::PropagateInit& init, literal_t literal, CoVarVec const &covec, value_t rhs, bool strict);
    //! Add up to two edges for a constraint.
    void add_edges_(Clingo::PropagateInit& init, vertex_t u_id, vertex_t v_id, value_t weight, literal_t lit, bool strict);
    //! Add (up to one) edge for a constraint.
    void add_edge_(Clingo::PropagateInit &init, vertex_t u_id, vertex_t v_id, value_t weight, literal_t lit);
    //! Reset connected components (to be recalculated with the next call to cc_calculate_).
    void cc_reset_();
    //! Mark a vertex in a component as visited.
    [[nodiscard]] bool cc_visited_(vertex_t index) const;
    //! Check if the given vertex is a zero vertex.
    [[nodiscard]] bool is_zero_(vertex_t index) const;
    //! Calculate the connected components.
    void cc_calculate_(AdjacencyMap &outgoing, AdjacencyMap &incoming);
    //! Calculate mutually exclusive edges.
    void calculate_mutexes_(Clingo::PropagateInit &init, edge_t edge_start, AdjacencyMap &outgoing);
    //! Finalize initialization by initializing thread states for propagation.
    void initialize_states_(Clingo::PropagateInit &init);

    // propagation functions
    //! Disable all edges associated with the given literal in the thread state.
    void disable_edge_by_lit(ThreadState &state, literal_t lit);
    //! Compute the current potential of a vertex in the given graph.
    [[nodiscard]] value_t get_potential_(Graph const &graph, vertex_t index) const;
    //! Compute the current cost of a vertex in the given graph.
    [[nodiscard]] value_t cost_(Graph const &graph, Edge const &edge) const;
    //! Compute a weight to sort vertices before propagation.
    [[nodiscard]] value_t cost_(SortMode mode, Graph const &graph, edge_t index) const;
    //! Sort vertices in the propagation queue according to the given mode.
    void sort_edges(SortMode mode, ThreadState &state);
    //! Propagate the given literals.
    void do_propagate(Clingo::PropagateControl &ctl, Clingo::LiteralSpan changes);

    ThreadStateVec states_;        //!< Thread specific state.
    FactStateVec facts_;           //!< Thread specific state for fact propagation.
    LitEdgeMap lit_to_edges_;      //!< Map from literals to edges.
    std::vector<Edge> edges_;      //!< Vector holding all edges.
    SymVertexMap vert_map_inv_;    //!< Mapping from symbols to vertex indices.
    VertexInfoVec vertex_info_;    //!< Node specific information.
    VertexIndexVec zero_vertices_; //!< Vector holding all zero vertices.
    Statistics &stats_;            //!< Statistics for the propagator.
    PropagatorConfig conf_;        //!< Configuration of the propagator.
    bool disable_edges_{false};    //!< Whether to disable edges.
};

} // namespace ClingoDL

#endif // CLINGODL_PROPAGATOR_HH
