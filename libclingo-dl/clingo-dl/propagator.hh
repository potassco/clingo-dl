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

template <typename T>
class DLPropagator : public Clingo::Propagator {
private:
    struct NodeInfo;
    struct ThreadState;
    struct FactState;
    using CoVarVec = ClingoDL::CoVarVec<T>; // vector of coefficients and variables

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

    void cc_reset_();

    bool cc_visited_(int node) const;

    bool is_zero_(int node) const;

    void cc_calculate_(std::unordered_multimap<int, int> &outgoing, std::unordered_multimap<int, int> &incoming);

    void calculate_mutexes_(Clingo::PropagateInit &init, int edge_start, std::unordered_multimap<int, int> &outgoing);

    void initialize_states_(Clingo::PropagateInit &init);

    // propagation

    void disable_edge_by_lit(ThreadState &state, Clingo::literal_t lit);

    int get_potential_(Graph<T> const &graph, int idx);

    int cost_(Graph<T> const &graph, Edge<T> const &edge);

    int cost_(SortMode mode, Graph<T> const &graph, int i);

    void sort_edges(SortMode mode, ThreadState &state);

    void do_propagate(Clingo::PropagateControl &ctl, Clingo::LiteralSpan changes);

    std::vector<ThreadState> states_;
    std::vector<FactState> facts_;
    std::unordered_multimap<Clingo::literal_t, int> lit_to_edges_;
    std::unordered_multimap<Clingo::literal_t, int> false_lit_to_edges_;
    std::vector<Edge<T>> edges_;
    std::vector<Clingo::Symbol> vert_map_;
    std::unordered_map<Clingo::Symbol, int> vert_map_inv_;
    std::vector<NodeInfo> node_info_;
    std::vector<int> zero_nodes_;
    Statistics &stats_;
    PropagatorConfig conf_;
};

} // namespace ClingoDL

#endif // CLINGODL_PROPAGATOR_HH
