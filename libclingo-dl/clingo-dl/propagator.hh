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

struct FactState {
    std::vector<Clingo::literal_t> lits;
    size_t limit{0};
};

template <typename T>
struct ThreadState {
    ThreadState(ThreadStatistics &stats, const std::vector<Edge<T>> &edges, PropagationMode propagate, uint64_t propagate_root, uint64_t propagate_budget)
        : stats(stats)
        , dl_graph(stats, edges, propagate)
        , propagate_root{propagate_root}
        , propagate_budget{propagate_budget} { }
    ThreadStatistics &stats;
    Graph<T> dl_graph;
    std::vector<Clingo::literal_t> false_lits;
    std::vector<int> todo_edges;
    uint64_t propagate_root;
    uint64_t propagate_budget;
};

template <typename T>
Clingo::Symbol to_symbol(T value);

template <>
inline Clingo::Symbol to_symbol(int value) {
    return Clingo::Number(value);
}

template <>
inline Clingo::Symbol to_symbol(double value) {
    return Clingo::String(std::to_string(value).c_str());
}

struct NodeInfo {
    NodeInfo(uint32_t cc = 0, bool visited = false) : cc(cc), visited(visited) {}
    uint32_t cc : 31;
    uint32_t visited : 1;
};

template <typename T>
class DifferenceLogicPropagator : public Clingo::Propagator {
private:
    using CoVarVec = ClingoDL::CoVarVec<T>; // vector of coefficients and variables

public:
    DifferenceLogicPropagator(Statistics &stats, PropagatorConfig conf)
    : stats_{stats}
    , conf_{std::move(conf)} {
        zero_nodes_.emplace_back(map_vertex_(Clingo::Number(0)));
        cc_reset_();
    }

    //! Get the number of vertices in the graph.
    size_t num_vertices() const {
        return vert_map_.size();
    }

    //! Get the symbol associated with a vertex index.
    Clingo::Symbol symbol(size_t index) const {
        return vert_map_[index];
    }

    //! Lookup the index of a vertex.
    uint32_t lookup(clingo_symbol_t symbol) {
        auto it = vert_map_inv_.find(Clingo::Symbol(symbol));
        return it != vert_map_inv_.end()
            ? it->second
            : num_vertices();
    }

    //! Check if teh given vertex has a lower bound in the given thread.
    bool has_lower_bound(uint32_t thread_id, size_t index) const {
        return index < vert_map_.size() && !is_zero_(index) && states_[thread_id].dl_graph.has_value(index);
    }

    //! Get the lower bound of a vertex in the given thread.
    T lower_bound(uint32_t thread_id, size_t index) const {
        assert(has_lower_bound(thread_id, index));
        auto &state = states_[thread_id];
        T adjust = 0;
        auto cc = node_info_[index].cc;
        auto zero_node = zero_nodes_[cc];

        if (state.dl_graph.has_value(zero_node)) {
            adjust = state.dl_graph.node_value(zero_node);
        }
        return state.dl_graph.node_value(index) - adjust;
    }

    //! Extend the model with vertex assignments.
    void extend_model(Clingo::Model &model) {
        auto &state = states_[model.thread_id()];
        std::vector<T> adjust;
        adjust.reserve(zero_nodes_.size());
        for (auto node : zero_nodes_) {
            adjust.emplace_back(state.dl_graph.has_value(node) ? state.dl_graph.node_value(node) : 0);
        }

        Clingo::SymbolVector vec;
        for (auto idx = 0; idx < vert_map_.size(); ++idx) {
            if (!is_zero_(idx) && state.dl_graph.has_value(idx)) {
                Clingo::SymbolVector params;
                params.emplace_back(vert_map_[idx]);
                auto cc = node_info_[idx].cc;
                params.emplace_back(to_symbol<T>(state.dl_graph.node_value(idx) - adjust[cc]));
                vec.emplace_back(Function("dl", params));
            }
        }
        model.extend(vec);
    }

    // propagator interface

    //! Initialize the propagator.
    void init(Clingo::PropagateInit &init) override {
        if (!edges_.empty()) {
            init.set_check_mode(Clingo::PropagatorCheckMode::Partial);
        }

        int edge_start = edges_.size();

        Timer t{stats_.time_init};
        add_constraints_(init);

        // build adjacency list
        std::unordered_multimap<int, int> outgoing;
        std::unordered_multimap<int, int> incoming;
        for (int edge_id = 0, size = edges_.size(); edge_id < size; ++edge_id) {
            outgoing.emplace(std::make_pair(edges_[edge_id].from, edge_id));
            incoming.emplace(std::make_pair(edges_[edge_id].to, edge_id));
        }

        cc_calculate_(outgoing, incoming);

        stats_.edges = edges_.size();
        stats_.variables = num_vertices();

        if (conf_.mutex_size > 0 && conf_.mutex_cutoff > 0) {
            calculate_mutexes_(init, edge_start, outgoing);
        }
        initialize_states_(init);
    }

    //! Propagate edges.
    void propagate(Clingo::PropagateControl &ctl, Clingo::LiteralSpan changes) override {
        if (ctl.assignment().decision_level() == 0) {
            auto &facts = facts_[ctl.thread_id()];
            facts.lits.insert(facts.lits.end(), changes.begin(), changes.end());
        }
        do_propagate(ctl, changes);
    }

    //! Undo propgated edges.
    void undo(Clingo::PropagateControl const &ctl, Clingo::LiteralSpan changes) noexcept override {
        static_cast<void>(changes);
        auto &state = states_[ctl.thread_id()];
        Timer t{state.stats.time_undo};
        state.dl_graph.backtrack();
    }

    void check(Clingo::PropagateControl &ctl) override {
        ThreadState<T> &state = states_[ctl.thread_id()];
        auto &facts = facts_[ctl.thread_id()];
        auto assignment = ctl.assignment();
        if (assignment.decision_level() == 0 && facts.limit > 0) {
            do_propagate(ctl, {facts.lits.data(), facts.lits.data() + facts.limit}); // NOLINT
            facts.limit = 0;
        }
#if defined(CLINGODL_CHECKSOLUTION) || defined(CLINGODL_CROSSCHECK)
        if (ctl.assignment().is_total()) {
            for (auto &x : edges_) {
                if (ctl.assignment().is_true(x.lit)) {
                    if (!state.dl_graph.node_value_defined(x.from) || !state.dl_graph.node_value_defined(x.to) || !(state.dl_graph.node_value(x.from) - state.dl_graph.node_value(x.to) <= x.weight)) {
                        throw std::logic_error("not a valid solution");
                    }
                }
            }
        }
#endif
    }

private:
    // initialization functions

    //! Map a symbol to an integer.
    int map_vertex_(Clingo::Symbol v) {
        auto ret = vert_map_inv_.emplace(v, static_cast<int>(vert_map_.size()));
        if (ret.second) {
            vert_map_.emplace_back(ret.first->first);
        }
        return ret.first->second;
    }

    //! Add constraints in the theory data.
    void add_constraints_(Clingo::PropagateInit &init) {
        for (auto atom : init.theory_atoms()) {
            auto term = atom.term();
            if (match(term, "__diff_h", 0) || match(term, "__diff_b", 0)) {
                auto edge = parse<T>(atom, [this](Clingo::Symbol const &sym) { return map_vertex_(sym); });
                int lit = init.solver_literal(atom.literal());
                normalize_constraint_(init, lit, edge.lhs, edge.rel, edge.rhs, edge.strict);
            }
        }
    }

    //! Normalize constraints to individual edges over `<=`.
    bool normalize_constraint_(Clingo::PropagateInit &init, int literal, CoVarVec const &elements, char const *op, T rhs, bool strict) {
        // rewrite '>', '<', and '>=' into '<='
        if (std::strcmp(op, ">") == 0) {
            op = ">=";
            rhs = safe_add<T>(rhs, epsilon<T>());
        }
        else if (std::strcmp(op, "<") == 0) {
            op = "<=";
            rhs = safe_sub<T>(rhs, epsilon<T>());
        }

        if (std::strcmp(op, ">=") == 0) {
            CoVarVec copy;
            copy.reserve(elements.size());
            for (auto &[co, var] : elements) {
                copy.emplace_back(safe_inv<T>(co), var);
            }
            return normalize_constraint_(init, literal, copy, "<=", safe_inv<T>(rhs), strict);
        }

        // hanle remainig '<=', '=', and '!='
        if (std::strcmp(op, "<=") == 0) {
            if (!init.assignment().is_true(-literal) && !add_edges_(init, literal, elements, rhs, false)) {
                return false;
            }
        }
        else if (std::strcmp(op, "=") == 0) {
            int a{0};
            int b{0};
            if (strict) {
                if (init.assignment().is_true(literal)) {
                    a = b = 1;
                }
                else {
                    a = init.add_literal();
                    b = init.add_literal();
                }

                // Note: this cannot fail because constraint normalization does not propagate
                if (!init.add_clause({-literal, a})) {
                    return false;
                }
                if (!init.add_clause({-literal, b})) {
                    return false;
                }
                if (!init.add_clause({-a, -b, literal})) {
                    return false;
                }
            }
            else {
                a = b = literal;
            }

            if (!normalize_constraint_(init, a, elements, "<=", rhs, strict)) {
                return false;
            }
            if (!normalize_constraint_(init, b, elements, ">=", rhs, strict)) {
                return false;
            }

            if (strict) {
                return true;
            }
        }
        else if (std::strcmp(op, "!=") == 0) {
            if (strict) {
                return normalize_constraint_(init, -literal, elements, "=", rhs, true);
            }

            auto a = init.add_literal();
            auto b = init.add_literal();

            if (!init.add_clause({a, b, -literal})) {
                return false;
            }
            if (!init.add_clause({-a, -b})) {
                return false;
            }
            if (!init.add_clause({literal, -a})) {
                return false;
            }
            if (!init.add_clause({literal, -b})) {
                return false;
            }

            if (!normalize_constraint_(init, a, elements, "<", rhs, false)) {
                return false;
            }
            if (!normalize_constraint_(init, b, elements, ">", rhs, false)) {
                return false;
            }
        }

        if (strict) {
            assert(std::strcmp(op, "=") != 0);

            if (std::strcmp(op, "<=") == 0) {
                op = ">";
            }
            else if (std::strcmp(op, "!=") == 0) {
                op = "=";
            }

            if (!normalize_constraint_(init, -literal, elements, op, rhs, false)) {
                return false;
            }
        }

        return true;
    }

    //! Add up to two edges for the given constraint if the has at most 2 variables and suitable coefficients.
    bool add_edges_(Clingo::PropagateInit& init, int literal, CoVarVec const &covec, T rhs, bool strict) {
        char const *msg = "normalizing difference constraint failed: only constraints of form &diff {u - v} <= b are accepted";
        if (strict && init.assignment().is_false(literal)) {
            return true;
        }
        if (covec.size() > 2) {
            throw std::runtime_error(msg);
        }
        auto u_id = map_vertex_(Clingo::Number(0));
        auto v_id = map_vertex_(Clingo::Number(0));
        if (covec.size() == 0) {
            if (rhs < 0) {
                return init.add_clause({-literal});
            }
            return !strict || init.add_clause({literal});
        }

        if (covec.size() == 1) {
            if (covec[0].first == 1) {
                u_id = covec[0].second;
            }
            else if (covec[0].first == -1) {
                v_id = covec[0].second;
            }
            else {
                throw std::runtime_error(msg);
            }
        }
        else if (covec.size() == 2) {
            if (covec[0].first == 1) {
                u_id = covec[0].second;
                if (covec[1].first == -1) {
                    v_id = covec[1].second;
                }
                else {
                    throw std::runtime_error(msg);
                }
            }
            else if (covec[0].first == -1) {
                v_id = covec[0].second;
                if (covec[1].first == 1) {
                    u_id = covec[1].second;
                }
                else {
                    throw std::runtime_error(msg);
                }
            }
            else {
                throw std::runtime_error(msg);
            }
        }
        add_edges_(init, u_id, v_id, rhs, literal, strict);
        return true;
    }

    //! Add up to two edges for a constraint.
    void add_edges_(Clingo::PropagateInit& init, int u_id, int v_id, T weight, int lit, bool strict) {
        add_edge_(init, u_id, v_id, weight, lit);
        if (strict) {
            add_edge_(init, v_id, u_id, -weight-1, -lit);
        }
    }

    //! Add (up to one) edge for a constraint.
    void add_edge_(Clingo::PropagateInit &init, int u_id, int v_id, T weight, int lit) {
        auto id = numeric_cast<int>(edges_.size());
        edges_.push_back({u_id, v_id, weight, lit});
        lit_to_edges_.emplace(lit, id);
        bool add = false;
        for (int i = 0; i < init.number_of_threads(); ++i) {
            init.add_watch(lit, i);
            if (conf_.get_propagate_mode(i) >= PropagationMode::Strong || conf_.get_propagate_root(i) > 0 || conf_.get_propagate_budget(i) > 0) {
                add = true;
                init.add_watch(-lit, i);
            }
        }
        if (add) {
            false_lit_to_edges_.emplace(-lit, id);
        }
    }

    void cc_reset_() {
        node_info_.clear();
        node_info_.resize(vert_map_.size());
        for (unsigned int i = 0; i < zero_nodes_.size(); ++i) {
            node_info_[zero_nodes_[i]] = NodeInfo(i, true);
        }
    }

    bool cc_visited_(int node) const {
        return node_info_[node].visited;
    }

    bool is_zero_(int node) const {
        assert(node_info_[node].cc < zero_nodes_.size());
        return zero_nodes_[node_info_[node].cc] == node;
    }

    void cc_calculate_(std::unordered_multimap<int, int> &outgoing, std::unordered_multimap<int, int> &incoming) {
        uint32_t cc = 0;
        // Note that this marks zero nodes as visited.
        cc_reset_();

        std::vector<int> node_stack;
        for (int node = 0; node < vert_map_.size(); ++node) {
            if (cc_visited_(node)) {
                continue;
            }
            node_info_[node] = NodeInfo(cc, true);
            node_stack.emplace_back(node);
            while (!node_stack.empty()) {
                auto node = node_stack.back();
                node_stack.pop_back();
                auto edges = outgoing.equal_range(node);
                for (auto edge = edges.first; edge != edges.second; ++edge) {
                    auto add_node = edges_[edge->second].to;
                    if (!cc_visited_(add_node)) {
                        node_info_[add_node] = NodeInfo(cc, true);
                        node_stack.emplace_back(add_node);
                    }
                }
                edges = incoming.equal_range(node);
                for (auto edge = edges.first; edge != edges.second; ++edge) {
                    auto add_node = edges_[edge->second].from;
                    if (!cc_visited_(add_node)) {
                        node_info_[add_node] = NodeInfo(cc, true);
                        node_stack.emplace_back(add_node);
                    }
                }
            }
            ++cc;
        }
        stats_.ccs = cc;

        zero_nodes_.reserve(cc);
        for (int i = zero_nodes_.size(); i < cc; ++i) {
            auto node = map_vertex_(Clingo::Function("__null", {Clingo::Number(i)}));
            zero_nodes_.emplace_back(node);
            node_info_.resize(std::max(node_info_.size(), static_cast<size_t>(node + 1)));
            node_info_[node] = NodeInfo(i, true);
        }

        std::vector<std::pair<int, int>> outgoing_change;
        std::vector<std::pair<int, int>> incoming_change;
        for (auto zero_node : zero_nodes_) {
            auto range = outgoing.equal_range(zero_node);
            for (auto edge = range.first; edge != range.second; ++edge) {
                auto &e = edges_[edge->second];
                auto cc = node_info_[e.to].cc;
                e.from = zero_nodes_[cc];
                outgoing_change.emplace_back(zero_nodes_[cc], edge->second);
            }
            outgoing.erase(range.first, range.second);
            range = incoming.equal_range(zero_node);
            for (auto edge = range.first; edge != range.second; ++edge) {
                auto &e = edges_[edge->second];
                auto cc = node_info_[e.from].cc;
                e.to = zero_nodes_[cc];
                incoming_change.emplace_back(zero_nodes_[cc], edge->second);
            }
            incoming.erase(range.first, range.second);
        }
        outgoing.insert(outgoing_change.begin(), outgoing_change.end());
        incoming.insert(incoming_change.begin(), incoming_change.end());
    }

    void calculate_mutexes_(Clingo::PropagateInit &init, int edge_start, std::unordered_multimap<int, int> &outgoing) {
        // let r and s be edge literals and T be the true literal:
        //
        //            r     T
        //         o --> o --> o
        //         ^           |
        //       T |           | r
        //         |           v
        //         o <-- o <-- o
        //            s     s
        //
        // if the above graph gives rise to a negative cycle, then r and s are mutually exclusive
        // the algorithm below adds clauses excluding such mutexes
        struct State {
            T weight;
            int id;
            int n;
            int previous;
        };
        std::vector<State> queue;
        std::vector<Clingo::literal_t> clause;

        auto ass = init.assignment();

        // traverse graph starting from each edge
        for (int start_id = edge_start, size = edges_.size(); start_id < size; ++start_id) {
            auto &start = edges_[start_id];
            // skipping over true literals forgoes some mutexes in the incremental case
            // but makes the algorithm much faster when there are many static edges
            // which are checked on level zero anyway
            if (ass.truth_value(start.lit) != Clingo::TruthValue::Free) {
                continue;
            }

            queue.emplace_back(State{start.weight, start_id, 1, -1});
            for (int queue_offset = 0; queue_offset < queue.size(); ++queue_offset) {
                auto rs_state = queue[queue_offset];
                auto rs = edges_[rs_state.id];
                auto out = outgoing.equal_range(rs.to);
                for (auto it = out.first; it != out.second; ++it) {
                    auto st_id = it->second;
                    auto &st = edges_[st_id];
                    auto st_truth = ass.truth_value(st.lit);
                    if ((st_id > start_id && st_truth == Clingo::TruthValue::Free) || st_truth == Clingo::TruthValue::False) {
                        continue;
                    }
                    auto w = rs_state.weight + st.weight;
                    auto n = rs_state.n;
                    auto c = queue_offset;
                    int found = 0;
                    while (c != -1) {
                        auto &cc = edges_[queue[c].id];
                        if (cc.lit == -st.lit || queue[c].id == st_id) {
                            found = 2;
                            break;
                        }
                        if (cc.lit == st.lit) {
                            found = 1;
                        }
                        c = queue[c].previous;
                    }
                    if (found == 2) {
                        continue;
                    }
                    if (found == 0 && st_truth == Clingo::TruthValue::Free) {
                        n += 1;
                    }
                    if (st.to == start.from && w < 0) {
                        ++stats_.mutexes;
                        clause.emplace_back(-st.lit);
                        auto c = queue_offset;
                        while (c != -1) {
                            auto &cc = edges_[queue[c].id];
                            clause.emplace_back(-cc.lit);
                            c = queue[c].previous;
                        }
                        if (!init.add_clause(clause)) {
                            return;
                        }
                        clause.clear();
                    }
                    else if (n < conf_.mutex_size && queue.size() < conf_.mutex_cutoff) {
                        queue.emplace_back(State{w, st_id, n, queue_offset});
                    }
                }
            }
            queue.clear();
        }
    }

    void initialize_states_(Clingo::PropagateInit &init) {
        stats_.dl_stats.resize(init.number_of_threads());
        states_.clear();
        if (facts_.size() < init.number_of_threads()) {
            facts_.resize(init.number_of_threads());
        }
        for (int i = 0; i < init.number_of_threads(); ++i) {
            states_.emplace_back(stats_.dl_stats[i], edges_, conf_.get_propagate_mode(i), conf_.get_propagate_root(i), conf_.get_propagate_budget(i));
            facts_[i].limit = facts_[i].lits.size();
        }
    }

    // propagation

    void disable_edge_by_lit(ThreadState<T> &state, Clingo::literal_t lit) {
        for (auto it = false_lit_to_edges_.find(lit), ie = false_lit_to_edges_.end(); it != ie && it->first == lit; ++it) {
            if (state.dl_graph.edge_is_active(it->second)) {
                state.dl_graph.remove_candidate_edge(it->second);
            }
        }
    }

    int get_potential_(Graph<T> const &graph, int idx) {
        return graph.node_value_defined(idx) ? -graph.node_value(idx) : 0;
    };

    int cost_(Graph<T> const &graph, Edge<T> const &edge) {
        return get_potential_(graph, edge.from) + edge.weight - get_potential_(graph, edge.to);
    };

    int cost_(SortMode mode, Graph<T> const &graph, int i) {
        switch(mode) {
            case SortMode::Weight: {
                return edges_[i].weight;
            }
            case SortMode::WeightRev: {
                return -edges_[i].weight;
            }
            case SortMode::Potential: {
                return cost_(graph, edges_[i]);
            }
            case SortMode::PotentialRev: {
                return -cost_(graph, edges_[i]);
            }
            case SortMode::No: {
                break;
            }
        }
        return 0;
    }

    void sort_edges(SortMode mode, ThreadState<T> &state) {
        std::sort(state.todo_edges.begin(), state.todo_edges.end(), [&](int l, int r) {
            return cost_(mode, state.dl_graph, l) < cost_(mode, state.dl_graph, r);
        });
    }

    void do_propagate(Clingo::PropagateControl &ctl, Clingo::LiteralSpan changes) {
        auto thread_id = ctl.thread_id();
        ThreadState<T> &state = states_[thread_id];
        Timer t{state.stats.time_propagate};
        auto level = ctl.assignment().decision_level();
        bool enable_propagate = state.dl_graph.mode() >= PropagationMode::Strong || level < state.propagate_root || state.propagate_budget > 0;
        state.dl_graph.ensure_decision_level(level, enable_propagate);
        if (state.dl_graph.can_propagate()) {
            for (auto &lit : state.false_lits) {
                if (ctl.assignment().is_true(lit)) { disable_edge_by_lit(state, lit); }
                ctl.add_watch(lit);
            }
            state.false_lits.clear();
        }
        state.todo_edges.clear();
        for (auto lit : changes) {
            auto it = lit_to_edges_.find(lit);
            auto ie = lit_to_edges_.end();
            if (state.dl_graph.can_propagate()) { disable_edge_by_lit(state, lit); }
            else if (it == ie) {
                state.false_lits.emplace_back(lit);
                ctl.remove_watch(lit);
            }
            for (; it != ie && it->first == lit; ++it) {
                if (state.dl_graph.edge_is_active(it->second)) {
                    state.todo_edges.push_back(it->second);
                }
            }
        }
        sort_edges(conf_.get_sort_mode(thread_id), state);
        for (auto edge : state.todo_edges) {
            if (state.dl_graph.edge_is_active(edge)) {
                auto ret = state.dl_graph.add_edge(edge, [&](std::vector<int> const &neg_cycle) {
                    std::vector<Clingo::literal_t> clause;
                    for (auto eid : neg_cycle) {
                        auto lit = -edges_[eid].lit;
                        if (ctl.assignment().is_true(lit)) { return true; }
                        clause.emplace_back(lit);
                    }
                    return ctl.add_clause(clause) && ctl.propagate();
                });
                if (!ret) { return; }
                bool propagate = (state.dl_graph.mode() >= PropagationMode::Strong) ||
                    (level < state.propagate_root) || (
                        state.propagate_budget > 0 &&
                        state.dl_graph.can_propagate() &&
                        state.stats.propagate_cost_add + state.propagate_budget > state.stats.propagate_cost_from + state.stats.propagate_cost_to);
                if (!propagate) { state.dl_graph.disable_propagate(); }
                // if !propgate -> can no longer propagate!
                if (propagate && !state.dl_graph.propagate(edge, ctl)) { return; }
            }
        }
    }

    std::vector<ThreadState<T>> states_;
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
