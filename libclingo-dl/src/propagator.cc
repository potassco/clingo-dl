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

#include <clingo-dl/propagator.hh>

#include <unordered_set>

#define CLINGODL_CHECKSOLUTION

namespace ClingoDL {

namespace {

template <typename T, typename std::enable_if<std::is_integral_v<T>, bool>::type = true>
inline auto to_symbol(T value) -> Clingo::Symbol {
    return Clingo::Number(value);
}

template <typename T, typename std::enable_if<std::is_floating_point_v<T>, bool>::type = true>
inline auto to_symbol(T value) -> Clingo::Symbol {
    return Clingo::String(std::to_string(value).c_str());
}

} // namespace

void Statistics::reset() {
    time_init = std::chrono::steady_clock::duration::zero();
    ccs = 0;
    mutexes = 0;
    edges = 0;
    variables = 0;
    for (auto &i : thread_statistics) {
        i.reset();
    }
}

void Statistics::accu(Statistics const &x) {
    time_init += x.time_init;
    ccs = x.ccs;
    mutexes += x.mutexes;
    edges = x.edges;
    variables = x.variables;
    if (thread_statistics.size() < x.thread_statistics.size()) {
        thread_statistics.resize(x.thread_statistics.size());
    }
    auto it = x.thread_statistics.begin();
    for (auto &y : thread_statistics) {
        y.accu(*it++);
    }
}

//! Struct to store vertex specific information.
template <typename T> struct DLPropagator<T>::VertexInfo {
    VertexInfo(Clingo::Symbol symbol) : symbol{symbol}, cc{0}, visited{0} {}
    //! Mark the vertex visited in the given connected component.
    void set_visited(index_t cc, bool visited) {
        this->cc = cc;
        this->visited = visited;
    }

    Clingo::Symbol symbol;                //!< The symbol associated with the vertex.
    index_t cc : sizeof(index_t) * 8 - 1; //!< The connected component the vertex is in.
    index_t visited : 1;                  //!< Whether the vertex has been visited.
};

//! Struct to store thread specific state.
template <typename T> struct DLPropagator<T>::ThreadState {
    ThreadState(ThreadStatistics &stats, const std::vector<Edge> &edges, PropagationMode propagate,
                level_t propagate_root, uint64_t propagate_budget)
        : stats{stats}, graph{stats, edges, propagate}, propagate_root{propagate_root},
          propagate_budget{propagate_budget} {}

    ThreadStatistics &stats;               //!< Thread specific statistics.
    Graph graph;                           //!< The incremental graph associated with the thread.
    std::vector<literal_t> removed_watchs; //!< Literals from which watches have been removed.
    std::vector<edge_t> todo_edges;        //!< Edges that have to be propagated.
    uint64_t propagate_root;               //!< Propagation is disabled above this level.
    uint64_t propagate_budget;             //!< The maximum budget invested into propagation.
};

//! Struct to store facts to repropagate.
//!
//! The solver never backtracks the top-level level. When solving a new step,
//! we simply rebuild the per thread states and repropagate the facts stored
//! here.
template <typename T> struct DLPropagator<T>::FactState {
    std::vector<literal_t> lits; //!< Literals that are true on level 0.
    size_t limit{0};             //!< Offset until which literals have to be repropagated.
};

template <typename T>
DLPropagator<T>::DLPropagator(Statistics &stats, PropagatorConfig conf) : stats_{stats}, conf_{std::move(conf)} {
    zero_vertices_.emplace_back(map_vertex_(Clingo::Number(0)));
    cc_reset_();
}

template <typename T> DLPropagator<T>::~DLPropagator() = default;

template <typename T> auto DLPropagator<T>::num_vertices() const -> vertex_t { return vertex_info_.size(); }

template <typename T> auto DLPropagator<T>::symbol(vertex_t index) const -> Clingo::Symbol {
    return vertex_info_[index].symbol;
}

template <typename T> auto DLPropagator<T>::lookup(Clingo::Symbol symbol) -> vertex_t {
    auto it = vert_map_inv_.find(symbol);
    return it != vert_map_inv_.end() ? it->second : num_vertices();
}

template <typename T> auto DLPropagator<T>::has_lower_bound(Clingo::id_t thread_id, vertex_t index) const -> bool {
    return index < vertex_info_.size() && !is_zero_(index) && states_[thread_id].graph.has_value(index);
}

template <typename T> auto DLPropagator<T>::lower_bound(Clingo::id_t thread_id, vertex_t index) const -> value_t {
    assert(has_lower_bound(thread_id, index));
    auto &state = states_[thread_id];
    auto zero_vertex = zero_vertices_[vertex_info_[index].cc];
    T adjust = state.graph.has_value(zero_vertex) ? state.graph.get_value(zero_vertex) : 0;
    return state.graph.get_value(index) - adjust;
}

template <typename T> void DLPropagator<T>::extend_model(Clingo::Model &model) {
    auto &state = states_[model.thread_id()];
    Clingo::SymbolVector vec;
    for (vertex_t idx = 0; idx < numeric_cast<vertex_t>(vertex_info_.size()); ++idx) {
        if (!is_zero_(idx) && state.graph.has_value(idx)) {
            Clingo::SymbolVector params;
            auto zero_vertex = zero_vertices_[vertex_info_[idx].cc];
            T adjust = state.graph.has_value(zero_vertex) ? state.graph.get_value(zero_vertex) : 0;
            params.emplace_back(vertex_info_[idx].symbol);
            params.emplace_back(to_symbol<T>(state.graph.get_value(idx) - adjust));
            vec.emplace_back(Function("dl", params));
        }
    }
    model.extend(vec);
}

template <typename T> void DLPropagator<T>::init(Clingo::PropagateInit &init) {
    if (!edges_.empty()) {
        init.set_check_mode(Clingo::PropagatorCheckMode::Partial);
    }

    edge_t edge_start = edges_.size();

    Timer t{stats_.time_init};
    if (!add_constraints_(init)) {
        return;
    }

    // build adjacency list
    AdjacencyMap outgoing;
    AdjacencyMap incoming;
    for (edge_t edge_id = 0, size = numeric_cast<edge_t>(edges_.size()); edge_id < size; ++edge_id) {
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

template <typename T> void DLPropagator<T>::propagate(Clingo::PropagateControl &ctl, Clingo::LiteralSpan changes) {
    // add facts for propagation at the next step
    if (ctl.assignment().decision_level() == 0) {
        auto &facts = facts_[ctl.thread_id()];
        facts.lits.insert(facts.lits.end(), changes.begin(), changes.end());
    }
    do_propagate(ctl, changes);
}

template <typename T>
void DLPropagator<T>::undo(Clingo::PropagateControl const &ctl, Clingo::LiteralSpan changes) noexcept {
    static_cast<void>(changes);
    auto &state = states_[ctl.thread_id()];
    Timer t{state.stats.time_undo};
    state.graph.backtrack();
}

template <typename T> void DLPropagator<T>::check(Clingo::PropagateControl &ctl) {
    ThreadState &state = states_[ctl.thread_id()];
    auto &facts = facts_[ctl.thread_id()];
    auto assignment = ctl.assignment();
    // propagate facts from previous step
    if (assignment.decision_level() == 0 && facts.limit > 0) {
        do_propagate(ctl, {facts.lits.data(), facts.lits.data() + facts.limit}); // NOLINT
        facts.limit = 0;
    }
#if defined(CLINGODL_CHECKSOLUTION) || defined(CLINGODL_CROSSCHECK)
    if (ctl.assignment().is_total()) {
        for (auto &x : edges_) {
            if (ctl.assignment().is_true(x.lit)) {
                if (!state.graph.has_value(x.from) || !state.graph.has_value(x.to) ||
                    !(state.graph.get_value(x.from) - state.graph.get_value(x.to) <= x.weight)) {
                    throw std::logic_error("not a valid solution");
                }
            }
        }
    }
#endif
}

template <typename T> auto DLPropagator<T>::map_vertex_(Clingo::Symbol symbol) -> vertex_t {
    auto [it, ins] = vert_map_inv_.emplace(symbol, numeric_cast<vertex_t>(vertex_info_.size()));
    if (ins) {
        vertex_info_.emplace_back(it->first);
    }
    return it->second;
}

template <typename T> auto DLPropagator<T>::add_constraints_(Clingo::PropagateInit &init) -> bool {
    for (auto atom : init.theory_atoms()) {
        auto term = atom.term();
        if (match(term, "__diff_h", 0) || match(term, "__diff_b", 0)) {
            auto edge = parse<T>(atom, [this](Clingo::Symbol const &sym) { return map_vertex_(sym); });
            literal_t lit = init.solver_literal(atom.literal());
            if (!normalize_constraint_(init, lit, edge.lhs, edge.rel, edge.rhs, edge.strict)) {
                return false;
            }
        }
    }
    return true;
}

template <typename T>
auto DLPropagator<T>::normalize_constraint_(Clingo::PropagateInit &init, literal_t literal, CoVarVec const &elements,
                                            char const *op, T rhs, bool strict) -> bool { // NOLINT
    // rewrite '>', '<', and '>=' into '<='
    if (std::strcmp(op, ">") == 0) {
        op = ">=";
        rhs = safe_add<T>(rhs, epsilon<T>());
    } else if (std::strcmp(op, "<") == 0) {
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

    // hanle remaining '<=', '=', and '!='
    if (std::strcmp(op, "<=") == 0) {
        if (!init.assignment().is_true(-literal) && !add_edges_(init, literal, elements, rhs, false)) {
            return false;
        }
    } else if (std::strcmp(op, "=") == 0) {
        literal_t a = 0;
        literal_t b = 0;
        if (strict) {
            if (init.assignment().is_true(literal)) {
                a = b = 1;
            } else {
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
        } else {
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
    } else if (std::strcmp(op, "!=") == 0) {
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
        } else if (std::strcmp(op, "!=") == 0) {
            op = "=";
        }

        if (!normalize_constraint_(init, -literal, elements, op, rhs, false)) {
            return false;
        }
    }

    return true;
}

template <typename T>
auto DLPropagator<T>::add_edges_(Clingo::PropagateInit &init, literal_t literal, CoVarVec const &covec, T rhs,
                                 bool strict) -> bool {
    char const *msg =
        "normalizing difference constraint failed: only constraints of form &diff {u - v} <= b are accepted";
    if (strict && init.assignment().is_false(literal)) {
        return true;
    }
    if (covec.size() > 2) {
        throw std::runtime_error(msg);
    }
    auto u_id = map_vertex_(Clingo::Number(0));
    auto v_id = map_vertex_(Clingo::Number(0));
    if (covec.empty()) {
        if (rhs < 0) {
            return init.add_clause({-literal});
        }
        return !strict || init.add_clause({literal});
    }
    if (covec.size() == 1) {
        if (covec[0].first == 1) {
            u_id = covec[0].second;
        } else if (covec[0].first == -1) {
            v_id = covec[0].second;
        } else {
            throw std::runtime_error(msg);
        }
    } else if (covec.size() == 2) {
        if (covec[0].first == 1) {
            u_id = covec[0].second;
            if (covec[1].first == -1) {
                v_id = covec[1].second;
            } else {
                throw std::runtime_error(msg);
            }
        } else if (covec[0].first == -1) {
            v_id = covec[0].second;
            if (covec[1].first == 1) {
                u_id = covec[1].second;
            } else {
                throw std::runtime_error(msg);
            }
        } else {
            throw std::runtime_error(msg);
        }
    }
    add_edges_(init, u_id, v_id, rhs, literal, strict);
    return true;
}

template <typename T>
void DLPropagator<T>::add_edges_(Clingo::PropagateInit &init, vertex_t u_id, vertex_t v_id, value_t weight,
                                 literal_t lit, bool strict) {
    add_edge_(init, u_id, v_id, weight, lit);
    if (strict) {
        add_edge_(init, v_id, u_id, -weight - 1, -lit);
    }
}

template <typename T>
void DLPropagator<T>::add_edge_(Clingo::PropagateInit &init, vertex_t u_id, vertex_t v_id, value_t weight,
                                literal_t lit) {
    auto id = numeric_cast<edge_t>(edges_.size());
    edges_.push_back({u_id, v_id, weight, lit});
    lit_to_edges_.emplace(lit, id);
    for (int i = 0; i < init.number_of_threads(); ++i) {
        init.add_watch(lit, i);
        if (conf_.get_propagate_mode(i) >= PropagationMode::Zero || conf_.get_propagate_root(i) > 0 ||
            conf_.get_propagate_budget(i) > 0) {
            disable_edges_ = true;
            init.add_watch(-lit, i);
        }
    }
}

template <typename T> void DLPropagator<T>::cc_reset_() {
    for (auto &info : vertex_info_) {
        info.set_visited(0, false);
    }
    for (size_t i = 0; i < zero_vertices_.size(); ++i) {
        vertex_info_[zero_vertices_[i]].set_visited(i, true);
    }
}

template <typename T> auto DLPropagator<T>::cc_visited_(vertex_t index) const -> bool {
    return vertex_info_[index].visited;
}

template <typename T> auto DLPropagator<T>::is_zero_(vertex_t index) const -> bool {
    assert(vertex_info_[index].cc < zero_vertices_.size());
    return zero_vertices_[vertex_info_[index].cc] == index;
}

template <typename T> void DLPropagator<T>::cc_calculate_(AdjacencyMap &outgoing, AdjacencyMap &incoming) {
    index_t cc = 0;
    // Note that this marks zero vertices as visited.
    cc_reset_();

    std::vector<vertex_t> vertex_stack;
    for (vertex_t start_vertex = 0; start_vertex < numeric_cast<vertex_t>(vertex_info_.size()); ++start_vertex) {
        if (cc_visited_(start_vertex)) {
            continue;
        }
        if (!conf_.calculate_cc) {
            vertex_info_[start_vertex].set_visited(cc, true);
            continue;
        }

        vertex_info_[start_vertex].set_visited(cc, true);
        vertex_stack.emplace_back(start_vertex);
        while (!vertex_stack.empty()) {
            auto vertex = vertex_stack.back();
            vertex_stack.pop_back();
            auto edges = outgoing.equal_range(vertex);
            for (auto edge = edges.first; edge != edges.second; ++edge) {
                auto add_vertex = edges_[edge->second].to;
                if (!cc_visited_(add_vertex)) {
                    vertex_info_[add_vertex].set_visited(cc, true);
                    vertex_stack.emplace_back(add_vertex);
                }
            }
            edges = incoming.equal_range(vertex);
            for (auto edge = edges.first; edge != edges.second; ++edge) {
                auto add_vertex = edges_[edge->second].from;
                if (!cc_visited_(add_vertex)) {
                    vertex_info_[add_vertex].set_visited(cc, true);
                    vertex_stack.emplace_back(add_vertex);
                }
            }
        }
        ++cc;
    }
    stats_.ccs = cc;

    zero_vertices_.reserve(cc);
    for (auto i = numeric_cast<index_t>(zero_vertices_.size()); i < cc; ++i) {
        auto vertex = map_vertex_(Clingo::Function("__null", {Clingo::Number(numeric_cast<int>(i))}));
        zero_vertices_.emplace_back(vertex);
        vertex_info_[vertex].set_visited(i, true);
    }

    std::vector<std::pair<vertex_t, vertex_t>> outgoing_change;
    std::vector<std::pair<vertex_t, vertex_t>> incoming_change;
    for (auto zero_vertex : zero_vertices_) {
        auto range = outgoing.equal_range(zero_vertex);
        for (auto edge = range.first; edge != range.second; ++edge) {
            auto &e = edges_[edge->second];
            auto cc = vertex_info_[e.to].cc;
            e.from = zero_vertices_[cc];
            outgoing_change.emplace_back(zero_vertices_[cc], edge->second);
        }
        outgoing.erase(range.first, range.second);
        range = incoming.equal_range(zero_vertex);
        for (auto edge = range.first; edge != range.second; ++edge) {
            auto &e = edges_[edge->second];
            auto cc = vertex_info_[e.from].cc;
            e.to = zero_vertices_[cc];
            incoming_change.emplace_back(zero_vertices_[cc], edge->second);
        }
        incoming.erase(range.first, range.second);
    }
    outgoing.insert(outgoing_change.begin(), outgoing_change.end());
    incoming.insert(incoming_change.begin(), incoming_change.end());
}

template <typename T>
void DLPropagator<T>::calculate_mutexes_(Clingo::PropagateInit &init, edge_t edge_start,
                                         AdjacencyMap &outgoing) { // NOLINT
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
        value_t weight;
        edge_t id;
        index_t n;
        index_t previous;
    };
    static constexpr auto invalid = std::numeric_limits<index_t>::max();
    std::vector<State> queue;
    std::vector<literal_t> clause;

    auto ass = init.assignment();

    // traverse graph starting from each edge
    for (edge_t start_id = edge_start, size = numeric_cast<edge_t>(edges_.size()); start_id < size; ++start_id) {
        auto &start = edges_[start_id];
        // skipping over true literals forgoes some mutexes in the incremental case
        // but makes the algorithm much faster when there are many static edges
        // which are checked on level zero anyway
        if (ass.truth_value(start.lit) != Clingo::TruthValue::Free) {
            continue;
        }

        queue.emplace_back(State{start.weight, start_id, 1, invalid});
        for (index_t queue_offset = 0; queue_offset < numeric_cast<index_t>(queue.size()); ++queue_offset) {
            auto rs_state = queue[queue_offset];
            auto rs = edges_[rs_state.id];
            auto out = outgoing.equal_range(rs.to);
            for (auto it = out.first; it != out.second; ++it) {
                auto st_id = it->second;
                auto &st = edges_[st_id];
                auto st_truth = ass.truth_value(st.lit);
                if ((st_id > start_id && st_truth == Clingo::TruthValue::Free) ||
                    st_truth == Clingo::TruthValue::False) {
                    continue;
                }
                auto w = rs_state.weight + st.weight;
                auto n = rs_state.n;
                auto c = queue_offset;
                int found = 0;
                while (c != invalid) {
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
                    while (c != invalid) {
                        auto &cc = edges_[queue[c].id];
                        clause.emplace_back(-cc.lit);
                        c = queue[c].previous;
                    }
                    if (!init.add_clause(clause)) {
                        return;
                    }
                    clause.clear();
                } else if (n < conf_.mutex_size && queue.size() < conf_.mutex_cutoff) {
                    queue.emplace_back(State{w, st_id, n, queue_offset});
                }
            }
        }
        queue.clear();
    }
}

template <typename T> void DLPropagator<T>::initialize_states_(Clingo::PropagateInit &init) {
    states_.clear();
    stats_.thread_statistics.resize(init.number_of_threads());
    if (facts_.size() < numeric_cast<size_t>(init.number_of_threads())) {
        facts_.resize(init.number_of_threads());
    }
    for (Clingo::id_t i = 0; i < numeric_cast<Clingo::id_t>(init.number_of_threads()); ++i) {
        states_.emplace_back(stats_.thread_statistics[i], edges_, conf_.get_propagate_mode(i),
                             conf_.get_propagate_root(i), conf_.get_propagate_budget(i));
        facts_[i].limit = facts_[i].lits.size();
    }
}

template <typename T> void DLPropagator<T>::disable_edge_by_lit(ThreadState &state, literal_t lit) {
    if (disable_edges_) {
        for (auto it = lit_to_edges_.find(-lit), ie = lit_to_edges_.end(); it != ie && it->first == -lit; ++it) {
            if (state.graph.edge_is_enabled(it->second)) {
                state.graph.disable_edge(it->second);
            }
        }
    }
}

template <typename T> auto DLPropagator<T>::get_potential_(Graph const &graph, vertex_t index) const -> value_t {
    return graph.has_value(index) ? -graph.get_value(index) : 0;
}

template <typename T> auto DLPropagator<T>::cost_(Graph const &graph, Edge const &edge) const -> value_t {
    return get_potential_(graph, edge.from) + edge.weight - get_potential_(graph, edge.to);
}

template <typename T> auto DLPropagator<T>::cost_(SortMode mode, Graph const &graph, edge_t index) const -> value_t {
    switch (mode) {
        case SortMode::Weight: {
            return edges_[index].weight;
        }
        case SortMode::WeightRev: {
            return -edges_[index].weight;
        }
        case SortMode::Potential: {
            return cost_(graph, edges_[index]);
        }
        case SortMode::PotentialRev: {
            return -cost_(graph, edges_[index]);
        }
        case SortMode::No: {
            break;
        }
    }
    return 0;
}

template <typename T> void DLPropagator<T>::sort_edges(SortMode mode, ThreadState &state) {
    std::sort(state.todo_edges.begin(), state.todo_edges.end(),
              [&](edge_t l, edge_t r) { return cost_(mode, state.graph, l) < cost_(mode, state.graph, r); });
}

template <typename T>
void DLPropagator<T>::do_propagate(Clingo::PropagateControl &ctl, Clingo::LiteralSpan changes) { // NOLINT
    // This function checks for conflicts and propagates edges if enabled. If
    // propagation is enabled, the graph has to be propagated after each edge
    // added. If limited propagation is enabled and the limit is reached,
    // propagation will be disabled for all further calls below the current
    // decision level.
    auto thread_id = ctl.thread_id();
    auto ass = ctl.assignment();
    ThreadState &state = states_[thread_id];
    Timer timer{state.stats.time_propagate};
    auto level = ass.decision_level();
    bool propagate = state.graph.mode() >= PropagationMode::Strong || level < state.propagate_root;
    state.graph.ensure_decision_level(level, propagate || state.propagate_budget > 0);

    // re-enable removed watches
    if (state.graph.can_propagate()) {
        // Note: If propagation is re-enabled, we re-add watches for literals
        // that have been removed. Since we removed the watches, some of them
        // might have become true unnoticed on earlier decision levels. We keep
        // such edges in the vector and simply disable them again. It would be
        // more efficient to re-add watches in the undo function.
        // Unfortunately, this is not supported by the current API.
        auto it = state.removed_watchs.begin();
        auto ie = state.removed_watchs.end();
        for (auto jt = it; jt != ie; ++jt) {
            auto truth = ass.truth_value(*jt);
            if (truth == Clingo::TruthValue::True) {
                disable_edge_by_lit(state, *jt);
            }
            if (truth == Clingo::TruthValue::Free || ass.level(*jt) == level) {
                ctl.add_watch(*jt);
                if (it != jt) {
                    *it = *jt;
                }
                ++it;
            }
        }
        state.removed_watchs.erase(it, ie);
    }

    // fill the todo queue
    state.todo_edges.clear();
    for (auto lit : changes) {
        auto it = lit_to_edges_.find(lit);
        auto ie = lit_to_edges_.end();
        if (state.graph.can_propagate()) {
            disable_edge_by_lit(state, lit);
        } else if (it == ie) {
            state.removed_watchs.emplace_back(lit);
            ctl.remove_watch(lit);
        }
        for (; it != ie && it->first == lit; ++it) {
            if (state.graph.edge_is_enabled(it->second)) {
                state.todo_edges.push_back(it->second);
            }
        }
    }

    // process edges in the todo queue
    sort_edges(conf_.get_sort_mode(thread_id), state);
    for (auto edge_idx : state.todo_edges) {
        if (state.graph.edge_is_enabled(edge_idx)) {
            // disable propagation once budget exceeded
            bool has_budget =
                state.propagate_budget > 0 && state.stats.propagate_cost_add + state.propagate_budget >
                                                  state.stats.propagate_cost_from + state.stats.propagate_cost_to;
            if (!propagate && !has_budget) {
                state.graph.disable_propagate();
            }
            auto &edge = edges_[edge_idx];
            auto &info = vertex_info_[edge.from];
            // check for conflicts
            if (!state.graph.add_edge(ctl, edge_idx, info.cc)) {
                return;
            }
        }
    }
}

template <typename T>
auto DLPropagator<T>::decide(id_t thread_id, Clingo::Assignment const &assign, literal_t fallback) -> literal_t {
    static_cast<void>(assign);
    if (conf_.decision_mode == DecisionMode::Disabled) {
        return fallback;
    }
    bool phase = conf_.decision_mode == DecisionMode::MinConflict;
    ThreadState &state = states_[thread_id];
    auto it = lit_to_edges_.find(fallback);
    if (it != lit_to_edges_.end() && state.graph.edge_is_negative(it->second) == phase) {
        return -fallback;
    }
    it = lit_to_edges_.find(-fallback);
    if (it != lit_to_edges_.end() && state.graph.edge_is_negative(it->second) != phase) {
        return -fallback;
    }
    return fallback;
}

template class DLPropagator<int>;
template class DLPropagator<double>;

} // namespace ClingoDL
