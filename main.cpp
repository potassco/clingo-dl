// {{{ GPL License

// This file is part of gringo - a grounder for logic programs.
// Copyright Roland Kaminski

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// }}}

#include <clingo.hh>
#include <iostream>
#include <string>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <queue>
#include <limits>
#include <chrono>
#include <iomanip>
#include <stdlib.h>

using namespace Clingo;

namespace Detail {

template <int X>
using int_type = std::integral_constant<int, X>;
template <class T, class S>
inline void nc_check(S s, int_type<0>) { // same sign
    (void)s;
    assert((std::is_same<T, S>::value) || (s >= std::numeric_limits<T>::min() && s <= std::numeric_limits<T>::max()));
}
template <class T, class S>
inline void nc_check(S s, int_type<-1>) { // Signed -> Unsigned
    (void)s;
    assert(s >= 0 && static_cast<S>(static_cast<T>(s)) == s);
}
template <class T, class S>
inline void nc_check(S s, int_type<1>) { // Unsigned -> Signed
    (void)s;
    assert(!(s > std::numeric_limits<T>::max()));
}

} // namespace Detail

template <class T, class S>
inline T numeric_cast(S s) {
    constexpr int sv = int(std::numeric_limits<T>::is_signed) - int(std::numeric_limits<S>::is_signed);
    ::Detail::nc_check<T>(s, ::Detail::int_type<sv>());
    return static_cast<T>(s);
}

template <typename T>
struct Edge {
    int from;
    int to;
    T weight;
    literal_t lit;
};

template <class K, class V>
std::ostream &operator<<(std::ostream &out, std::unordered_map<K, V> const &map);
template <class T>
std::ostream &operator<<(std::ostream &out, std::vector<T> const &vec);
template <class K, class V>
std::ostream &operator<<(std::ostream &out, std::pair<K, V> const &pair);

template <class T>
std::ostream &operator<<(std::ostream &out, std::vector<T> const &vec) {
    out << "{";
    for (auto &x : vec) {
        out << " " << x;
    }
    out << " }";
    return out;
}

template <class K, class V>
std::ostream &operator<<(std::ostream &out, std::unordered_map<K, V> const &map) {
    using T = std::pair<K, V>;
    std::vector<T> vec;
    vec.assign(map.begin(), map.end());
    std::sort(vec.begin(), vec.end(), [](T const &a, T const &b) { return a.first < b.first; });
    out << vec;
    return out;
}

template <class K, class V>
std::ostream &operator<<(std::ostream &out, std::pair<K, V> const &pair) {
    out << "( " << pair.first << " " << pair.second << " )";
    return out;
}

template <class C>
void ensure_index(C &c, size_t index) {
    if (index >= c.size()) {
        c.resize(index + 1);
    }
}

using Duration = std::chrono::duration<double>;

class Timer {
public:
    Timer(Duration &elapsed)
        : elapsed_(elapsed)
        , start_(std::chrono::steady_clock::now()) {}
    ~Timer() { elapsed_ += std::chrono::steady_clock::now() - start_; }

private:
    Duration &elapsed_;
    std::chrono::time_point<std::chrono::steady_clock> start_;
};

template <typename T>
struct DifferenceLogicNode {
    bool defined() const { return !potential_stack.empty(); }
    T potential() const { return potential_stack.back().second; }
    std::vector<int> outgoing;
    std::vector<int> incoming;
    std::vector<std::pair<int, T>> potential_stack; // [(level,potential)]
    int last_edge = 0;
    T gamma = 0;
    bool changed = false;
};

template <typename T>
struct DifferenceLogicNodeUpdate {
    int node_idx;
    T gamma;
};

template <class T>
bool operator<(DifferenceLogicNodeUpdate<T> const &a, DifferenceLogicNodeUpdate<T> const &b) {
    return a.gamma > b.gamma;
}

template <typename T>
class DifferenceLogicGraph {
public:
    DifferenceLogicGraph(const std::vector<Edge<T>> &edges)
        : edges_(edges)
        , offset_active_edges_(0) {
        for (int i = 0; i < numeric_cast<int>(edges_.size()); ++i) {
            edge_partition_.emplace_back(i);
        }
        active_edges_.resize(edge_partition_.size(), true);
    }

    bool empty() const { return nodes_.empty(); }

    int node_value_defined(int idx) const { return nodes_[idx].defined(); }
    T node_value(int idx) const { return -nodes_[idx].potential(); }

    bool edge_is_active(int edge_idx) const {
        return active_edges_[edge_idx];
    }

    std::vector<int> add_edge(int level, int uv_idx) {
        auto &uv = edges_[uv_idx];

        // initialize the trail
        if (changed_trail_.empty() || std::get<0>(changed_trail_.back()) < level) {
            assert(changed_trail_.empty() || std::get<3>(changed_trail_.back()) <= offset_active_edges_);
            changed_trail_.emplace_back(level, changed_nodes_.size(), changed_edges_.size(), offset_active_edges_);
        }

        // initialize the nodes of the edge to add
        ensure_index(nodes_, std::max(uv.from, uv.to));
        auto &u = nodes_[uv.from];
        auto &v = nodes_[uv.to];
        if (!u.defined()) {
            set_potential(u, level, 0);
        }
        if (!v.defined()) {
            set_potential(v, level, 0);
        }
        v.gamma = u.potential() + uv.weight - v.potential();
        if (v.gamma < 0) {
            gamma_.push({uv.to, v.gamma});
            v.last_edge = uv_idx;
        }

        // detect negative cycles
        while (!gamma_.empty() && u.gamma == 0) {
            auto s_idx = gamma_.top().node_idx;
            auto &s = nodes_[s_idx];
            if (!s.changed) {
                assert(s.gamma == gamma_.top().gamma);
                set_potential(s, level, s.potential() + s.gamma);
                s.gamma = 0;
                s.changed = true;
                changed_.emplace_back(s_idx);
                for (auto st_idx : s.outgoing) {
                    assert(st_idx < numeric_cast<int>(edges_.size()));
                    auto &st = edges_[st_idx];
                    auto &t = nodes_[st.to];
                    if (!t.changed) {
                        auto gamma = s.potential() + st.weight - t.potential();
                        if (gamma < t.gamma) {
                            t.gamma = gamma;
                            gamma_.push({st.to, t.gamma});
                            t.last_edge = st_idx;
                        }
                    }
                }
            }
            gamma_.pop();
        }

        std::vector<int> neg_cycle;
        if (u.gamma < 0) {
            // gather the edges in the negative cycle
            neg_cycle.push_back(v.last_edge);
            auto next_idx = edges_[v.last_edge].from;
            while (uv.to != next_idx) {
                auto &next = nodes_[next_idx];
                neg_cycle.push_back(next.last_edge);
                next_idx = edges_[next.last_edge].from;
            }
#ifdef CROSSCHECK
            T weight = 0;
            for (auto &edge_idx : neg_cycle) {
                weight += edges_[edge_idx].weight;
            }
            assert(weight < 0);
#endif
        }
        else {
            // add the edge to the graph
            u.outgoing.emplace_back(uv_idx);
            v.incoming.emplace_back(uv_idx);
            changed_edges_.emplace_back(uv_idx);
#ifdef CROSSCHECK
            // TODO: just a check that will throw if there is a cycle
            std::unordered_map<int, T> costs;
            bellman_ford(costs, uv.from);
#endif
        }

        // reset gamma and changed flags
        v.gamma = 0;
        while (!gamma_.empty()) {
            nodes_[gamma_.top().node_idx].gamma = 0;
            gamma_.pop();
        }
        for (auto x : changed_) {
            nodes_[x].changed = false;
        }
        changed_.clear();

        return neg_cycle;
    }

    bool propagate(int xy_idx, Clingo::PropagateControl &ctl) {
        bool ret = true;
        // these maps are of course completely unnecessary
        // the costs can be stored in the nodes as well
        // also the changed flags from the algorithm above can be reused!
        // then there is also the relevancy criterion described in the paper
        auto &xy = edges_[xy_idx];
        auto &x = nodes_[xy.from];
        auto &y = nodes_[xy.to];
        std::unordered_map<int, T> costs_from_y, costs_to_x;
        dijkstra(costs_from_y, xy.to, [](Edge<T> const &edge) { return edge.to; }, [](DifferenceLogicNode<T> const &node) { return node.outgoing; });
        dijkstra(costs_to_x, xy.from, [](Edge<T> const &edge) { return edge.from; }, [](DifferenceLogicNode<T> const &node) { return node.incoming; });
        offset_active_edges_ = std::partition(edge_partition_.begin() + offset_active_edges_, edge_partition_.end(), [&](int uv_idx) {
            if (uv_idx == xy_idx) {
                active_edges_[xy_idx] = false;
                return true;
            }
            assert(edge_is_active(uv_idx));
            auto &uv = edges_[uv_idx];
            auto yv_it = costs_from_y.find(uv.to), ux_it = costs_to_x.find(uv.from);
            auto &u = nodes_[uv.from];
            auto &v = nodes_[uv.to];
            if (yv_it != costs_from_y.end() && ux_it != costs_to_x.end()) {
                auto a = ux_it->second + x.potential() - u.potential();
                auto b = yv_it->second + v.potential() - y.potential();
                auto d = a + b + xy.weight;
#ifdef CROSSCHECK
                std::unordered_map<int, T> bf_costs_from_u, bf_costs_from_y;
                bellman_ford(bf_costs_from_u, uv.from);
                bellman_ford(bf_costs_from_y, xy.to);
                auto aa = bf_costs_from_u.find(xy.from);
                assert (aa != bf_costs_from_u.end());
                assert(aa->second == a);
                auto bb = bf_costs_from_y.find(uv.to);
                assert (bb != bf_costs_from_u.end());
                assert(bb->second == b);
#endif

                if (d <= uv.weight) {
                    active_edges_[uv_idx] = false;
#ifdef CROSSCHECK
                    auto cycle = add_edge(std::numeric_limits<int>::max(), uv_idx);
                    backtrack();
                    if (!cycle.empty()) { throw std::runtime_error("edge is implied but lead to a conflict :("); }
#endif
                    return true;
                }
            }
            yv_it = costs_from_y.find(uv.from), ux_it = costs_to_x.find(uv.to);
            if (yv_it != costs_from_y.end() && ux_it != costs_to_x.end()) {
                auto a = ux_it->second + x.potential() - v.potential();
                auto b = yv_it->second + u.potential() - y.potential();
                auto d = a + b + xy.weight;
                if (d <= -uv.weight - 1) {
                    active_edges_[uv_idx] = false;
                    // TODO: the cycle can be obtained from the dijkstra too
                    auto cycle = add_edge(std::numeric_limits<int>::max(), uv_idx);
                    backtrack();
                    assert (!cycle.empty());
                    std::vector<literal_t> clause;
                    for (auto eid : cycle) {
                        clause.emplace_back(-edges_[eid].lit);
                    }
                    ret = ret && ctl.add_clause(clause) && ctl.propagate();
                    return true;
                }
#ifdef CROSSCHECK
                else {
                    auto cycle = add_edge(std::numeric_limits<int>::max(), uv_idx);
                    backtrack();
                    if (!cycle.empty()) { throw std::runtime_error("edge should not cause a conflict!"); }
                }
#endif
            }
            return false;
        }) - edge_partition_.begin();
        return ret;
    }
    template <class To, class Out>
    void dijkstra(std::unordered_map<int, T> &costs, int source, To to, Out getOut) {
        std::priority_queue<DifferenceLogicNodeUpdate<T>> queue;
        queue.push({source, 0});
        costs[source] = 0;
        while (!queue.empty()) {
            auto top = queue.top();
            if (top.gamma == costs[top.node_idx]) {
                auto &u = nodes_[top.node_idx];
                for (auto &uv_idx : getOut(u)) {
                    auto &uv = edges_[uv_idx];
                    auto v_it = costs.find(to(uv));
                    // NOTE: explicitely using uv.from and uv.to is intended here
                    auto gamma = top.gamma + nodes_[uv.from].potential() + uv.weight - nodes_[uv.to].potential();
                    assert(nodes_[uv.from].potential() + uv.weight - nodes_[uv.to].potential() >= 0);
                    if (v_it == costs.end() || gamma < v_it->second) {
                        queue.push({to(uv), gamma});
                        costs[to(uv)] = gamma;
                    }
                }
            }
            queue.pop();
        }
    }

#ifdef CROSSCHECK
    void bellman_ford(std::unordered_map<int, T> &costs, int source) {
        costs[source] = 0;
        int nodes = 0;
        for (auto &node : nodes_) {
            if (node.defined()) {
                ++nodes;
            }
        }
        for (int i = 0; i < nodes; ++i) {
            for (auto &uv_idx : changed_edges_) {
                auto &uv = edges_[uv_idx];
                auto u_cost = costs.find(uv.from);
                if (u_cost != costs.end()) {
                    auto v_cost = costs.find(uv.to);
                    auto dist = u_cost->second + uv.weight;
                    if (v_cost == costs.end()) {
                        costs[uv.to] = dist;
                    }
                    else if (dist < v_cost->second) {
                        v_cost->second = dist;
                    }
                }
            }
        }
        for (auto &uv_idx : changed_edges_) {
            auto &uv = edges_[uv_idx];
            auto u_cost = costs.find(uv.from);
            if (u_cost != costs.end()) {
                auto v_cost = costs.find(uv.to);
                auto dist = u_cost->second + uv.weight;
                if (dist < v_cost->second) {
                    throw std::runtime_error("there is a negative cycle!!!");
                }
            }
        }
    }
#endif

    void backtrack() {
        for (int count = changed_nodes_.size() - std::get<1>(changed_trail_.back()); count > 0; --count) {
            auto &node = nodes_[changed_nodes_.back()];
            node.potential_stack.pop_back();
            changed_nodes_.pop_back();
        }
        for (int count = changed_edges_.size() - std::get<2>(changed_trail_.back()); count > 0; --count) {
            auto &edge = edges_[changed_edges_.back()];
            nodes_[edge.from].outgoing.pop_back();
            nodes_[edge.to].incoming.pop_back();
            changed_edges_.pop_back();
        }
        for (; offset_active_edges_ > std::get<3>(changed_trail_.back()); --offset_active_edges_) {
            active_edges_[edge_partition_[offset_active_edges_ - 1]] = true;
        }
        assert(offset_active_edges_ == std::get<3>(changed_trail_.back()));
        changed_trail_.pop_back();
    }

private:
    void set_potential(DifferenceLogicNode<T> &node, int level, T potential) {
        if (!node.defined() || node.potential_stack.back().first < level) {
            node.potential_stack.emplace_back(level, potential);
            changed_nodes_.emplace_back(numeric_cast<int>(&node - nodes_.data()));
        }
        else {
            node.potential_stack.back().second = potential;
        }
    }

private:
    std::priority_queue<DifferenceLogicNodeUpdate<T>> gamma_;
    std::vector<int> changed_;
    const std::vector<Edge<T>> &edges_;
    std::vector<int> edge_partition_;
    std::vector<bool> active_edges_;
    std::vector<DifferenceLogicNode<T>> nodes_;
    std::vector<int> changed_nodes_;
    std::vector<int> changed_edges_;
    std::vector<std::tuple<int, int, int, int>> changed_trail_;
    int offset_active_edges_;
};

struct DLStats {
    Duration time_propagate = Duration{0};
    Duration time_undo = Duration{0};
};

struct Stats {
    Duration time_total = Duration{0};
    Duration time_init = Duration{0};
    std::vector<DLStats> dl_stats;
    int64_t conflicts{0};
    int64_t choices{0};
    int64_t restarts{0};
};

template <typename T>
struct DLState {
    DLState(DLStats &stats, const std::vector<Edge<T>> &edges)
        : stats(stats)
        , dl_graph(edges) {}
    DLStats &stats;
    std::vector<int> edge_trail;
    DifferenceLogicGraph<T> dl_graph;
    int propagated = 0;
};

template <typename T>
T get_weight(TheoryAtom const &atom);
template <>
int get_weight(TheoryAtom const &atom) {
    return atom.guard().second.arguments().empty() ? atom.guard().second.number() : -atom.guard().second.arguments()[0].number();
}
template <>
double get_weight(TheoryAtom const &atom) {
    static const std::string chars = "()\"";
    // TODO: why not simply atom.guard().second.name();
    std::string guard = atom.guard().second.to_string();
    guard.erase(std::remove_if(guard.begin(), guard.end(), [](char c) { return chars.find(c) != std::string::npos; }), guard.end());
    return std::stod(guard);
}

template <typename T>
class DifferenceLogicPropagator : public Propagator {
public:
    DifferenceLogicPropagator(Stats &stats, bool strict)
        : stats_(stats)
        , strict_(strict) {}

    void print_assignment(int thread) const {
        auto &state = states_[thread];
        T adjust = 0;
        int idx = 0;
        for (std::string const &name : vert_map_) {
            if (name == "0") {
                adjust = state.dl_graph.node_value(idx);
                break;
            }
            ++idx;
        }

        std::cout << "with assignment:\n";
        idx = 0;
        for (std::string const &name : vert_map_) {
            if (state.dl_graph.node_value_defined(idx) && name != "0") {
                std::cout << name << ":" << adjust + state.dl_graph.node_value(idx) << " ";
            }
            ++idx;
        }
        std::cout << "\n";
    }

private:
    // initialization

    void init(PropagateInit &init) override {
        Timer t{stats_.time_init};
        for (auto atom : init.theory_atoms()) {
            auto term = atom.term();
            if (term.to_string() == "diff") {
                add_edge_atom(init, atom);
            }
        }
        initialize_states(init);
    }

    void add_edge_atom(PropagateInit &init, TheoryAtom const &atom) {
        int lit = init.solver_literal(atom.literal());
        T weight = get_weight<T>(atom);
        auto u_id = map_vert(atom.elements()[0].tuple()[0].arguments()[0].to_string());
        auto v_id = map_vert(atom.elements()[0].tuple()[0].arguments()[1].to_string());
        auto id = numeric_cast<int>(edges_.size());
        edges_.push_back({u_id, v_id, weight, lit});
        lit_to_edges_.emplace(lit, id);
        init.add_watch(lit);
        if (strict_) {
            auto id = numeric_cast<int>(edges_.size());
            edges_.push_back({v_id, u_id, -weight - 1, -lit});
            lit_to_edges_.emplace(-lit, id);
            init.add_watch(-lit);
        }
    }

    int map_vert(std::string v) {
        auto ret = vert_map_inv_.emplace(std::move(v), vert_map_.size());
        if (ret.second) {
            vert_map_.emplace_back(ret.first->first);
        }
        return ret.first->second;
    }

    void initialize_states(PropagateInit &init) {
        stats_.dl_stats.resize(init.number_of_threads());
        for (int i = 0; i < init.number_of_threads(); ++i) {
            states_.emplace_back(stats_.dl_stats[i], edges_);
        }
    }

    // propagation

    void propagate(PropagateControl &ctl, LiteralSpan changes) override {
        auto &state = states_[ctl.thread_id()];
        Timer t{state.stats.time_propagate};
        auto level = ctl.assignment().decision_level();
        // NOTE: vector copy only because clasp bug
        //       can only be triggered with propagation
        for (auto lit : std::vector<Clingo::literal_t>(changes.begin(), changes.end())) {
            for (auto it = lit_to_edges_.find(lit), ie = lit_to_edges_.end(); it != ie && it->first == lit; ++it) {
                if (state.dl_graph.edge_is_active(it->second)) {
                    auto neg_cycle = state.dl_graph.add_edge(level, it->second);
                    if (!neg_cycle.empty()) {
                        std::vector<literal_t> clause;
                        for (auto eid : neg_cycle) {
                            clause.emplace_back(-edges_[eid].lit);
                        }
                        if (!ctl.add_clause(clause) || !ctl.propagate()) {
                            return;
                        }
                        assert(false && "must not happen");
                    }
                    else {
                        if (!state.dl_graph.propagate(it->second, ctl)) {
                            return;
                        }
                    }
                }
            }
        }
    }

    // undo

void undo(PropagateControl const &ctl, LiteralSpan changes) override {
    static_cast<void>(changes);
        auto &state = states_[ctl.thread_id()];
        Timer t{state.stats.time_undo};
        state.dl_graph.backtrack();
    }

private:
    std::vector<DLState<T>> states_;
    std::unordered_multimap<literal_t, int> lit_to_edges_;
    std::vector<Edge<T>> edges_;
    std::vector<std::reference_wrapper<const std::string>> vert_map_;
    std::unordered_map<std::string, int> vert_map_inv_;
    Stats &stats_;
    bool strict_;
};

int get_int(std::string constname, Control &ctl, int def) {
    Symbol val = ctl.get_const(constname.c_str());
    if (val.to_string() == constname.c_str()) {
        return def;
    }
    return val.number();
}

template <typename T>
void solve(Stats &stats, Control &ctl, bool strict) {
    DifferenceLogicPropagator<T> p{stats, strict};
    ctl.register_propagator(p);
    ctl.ground({{"base", {}}});
    int i = 0;
    for (auto m : ctl.solve()) {
        i++;
        std::cout << "Answer " << i << "\n";
        std::cout << m << "\n";
        p.print_assignment(m.thread_id());
    }
    if (i == 0) {
        std::cout << "UNSATISFIABLE\n";
    }
    else {
        std::cout << "SATISFIABLE\n";
    }
}

int main(int argc, char *argv[]) {
    Stats stats;
    {
        Timer t{stats.time_total};
        auto argb = argv + 1, arge = argb;
        for (; *argb; ++argb, ++arge) {
            if (std::strcmp(*argb, "--") == 0) {
                ++argb;
                break;
            }
        }
        Control ctl{{argb, numeric_cast<size_t>(argv + argc - argb)}};
        ctl.add("base", {}, R"(#theory dl {
    term{};
    constant {- : 1, unary};
    diff_term {- : 1, binary, left};
    &diff/0 : diff_term, {<=}, constant, any;
    &show_assignment/0 : term, directive
}.)");
        for (auto arg = argv + 1; arg != arge; ++arg) {
            ctl.load(*arg);
        }
        // configure strict/non-strict mode
        auto strict = get_int("strict", ctl, 0) != 0;
        // configure IDL/RDL mode
        auto rdl = get_int("rdl", ctl, 0) != 0;

        if (rdl) {
            if (strict) {
                std::cout << "Real difference logic not available with strict semantics!" << std::endl;
                exit(EXIT_FAILURE);
            }
            solve<double>(stats, ctl, strict);
        }
        else {
            solve<int>(stats, ctl, strict);
        }
        auto solvers = ctl.statistics()["solving"]["solvers"];
        stats.choices = solvers["choices"];
        stats.conflicts = solvers["conflicts"];
        stats.restarts = solvers["restarts"];
    }

    std::cout << "total: " << stats.time_total.count() << "s\n";
    std::cout << "  init: " << stats.time_init.count() << "s\n";
    int thread = 0;
    for (auto &stat : stats.dl_stats) {
        std::cout << "  total[" << thread << "]: ";
        std::cout << (stat.time_undo + stat.time_propagate).count() << "s\n";
        std::cout << "    propagate: " << stat.time_propagate.count() << "s\n";
        std::cout << "    undo     : " << stat.time_undo.count() << "s\n";
        ++thread;
    }
    std::cout << "conflicts: " << stats.conflicts << "\n";
    std::cout << "choices  : " << stats.choices << "\n";
    std::cout << "restarts : " << stats.restarts << "\n";
}
