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

using namespace Clingo;

namespace Detail {

template <int X> using int_type = std::integral_constant<int, X>;
template <class T, class S> inline void nc_check(S s, int_type<0>) { // same sign
    (void)s;
    assert((std::is_same<T, S>::value) || (s >= std::numeric_limits<T>::min() && s <= std::numeric_limits<T>::max()));
}
template <class T, class S> inline void nc_check(S s, int_type<-1>) { // Signed -> Unsigned
    (void)s;
    assert(s >= 0 && static_cast<S>(static_cast<T>(s)) == s);
}
template <class T, class S> inline void nc_check(S s, int_type<1>) { // Unsigned -> Signed
    (void)s;
    assert(!(s > std::numeric_limits<T>::max()));
}

} // namespace Detail

template <class T, class S> inline T numeric_cast(S s) {
    constexpr int sv = int(std::numeric_limits<T>::is_signed) - int(std::numeric_limits<S>::is_signed);
    ::Detail::nc_check<T>(s, ::Detail::int_type<sv>());
    return static_cast<T>(s);
}

struct Edge {
    Edge(int i, int x, int y, int k) : id(i), from(x), to(y), weight(k) {}
    int id;
    int from;
    int to;
    int weight;
};

template <class K, class V> std::ostream &operator<<(std::ostream &out, std::unordered_map<K, V> const &map);
template <class T> std::ostream &operator<<(std::ostream &out, std::vector<T> const &vec);
template <class K, class V> std::ostream &operator<<(std::ostream &out, std::pair<K, V> const &pair);

template <class T> std::ostream &operator<<(std::ostream &out, std::vector<T> const &vec) {
    out << "{";
    for (auto &x : vec) {
        out << " " << x;
    }
    out << " }";
    return out;
}

template <class K, class V> std::ostream &operator<<(std::ostream &out, std::unordered_map<K, V> const &map) {
    using T = std::pair<K, V>;
    std::vector<T> vec;
    vec.assign(map.begin(), map.end());
    std::sort(vec.begin(), vec.end(), [](T const &a, T const &b) { return a.first < b.first; });
    out << vec;
    return out;
}

template <class K, class V> std::ostream &operator<<(std::ostream &out, std::pair<K, V> const &pair) {
    out << "( " << pair.first << " " << pair.second << " )";
    return out;
}

class DifferenceLogicGraph {
  private:
    class Graph {
      public:
        std::vector<std::vector<int>> outgoing;
        const std::vector<Edge> &edges;

        Graph(const std::vector<Edge> &edges) : edges(edges) {}

        void add_edge(int edge) {
            assert(edge < numeric_cast<int>(edges.size()));
            int maxid = std::max(edges[edge].from, edges[edge].to);
            if (maxid >= numeric_cast<int>(outgoing.size())) {
                outgoing.resize(maxid + 1, {});
            }
            outgoing[edges[edge].from].emplace_back(edge);
        }

        void reset() { outgoing.clear(); }
    };

    class PairComp {
        bool reverse;

      public:
        PairComp(const bool &revparam = false) { reverse = revparam; }
        bool operator()(const std::pair<int, int> &lhs, const std::pair<int, int> &rhs) const {
            if (reverse)
                return (lhs.second < rhs.second);
            else
                return (lhs.second > rhs.second);
        }
    };

    bool valid;
    Graph graph;
    std::vector<int> potential;
    std::vector<bool> changed;

  public:
    DifferenceLogicGraph(const std::vector<Edge> &edges) : valid(false), graph(edges) {}

    bool is_valid() { return this->valid; }

    std::unordered_map<int, int> get_assignment() {
        std::unordered_map<int, int> ass;
        int idx = 0;
        for (int pot : potential) {
            if (potential[idx] != std::numeric_limits<int>::max()) {
                ass[idx] = -pot;
            }
            idx++;
        }
        return ass;
    }

    std::vector<int> add_edge(int edge) {
        int u = graph.edges[edge].from;
        int v = graph.edges[edge].to;
        int d = graph.edges[edge].weight;
        graph.add_edge(edge);

        typedef std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, PairComp> PairPrioQueue;

        if (numeric_cast<int>(potential.size()) <= u) {
            potential.resize(u + 1, std::numeric_limits<int>::max());
            potential[u] = 0;
        } else if (potential[u] == std::numeric_limits<int>::max()) {
            potential[u] = 0;
        }

        if (numeric_cast<int>(potential.size()) <= v) {
            potential.resize(v + 1, std::numeric_limits<int>::max());
            potential[v] = potential[u] + d;
        } else if (potential[v] == std::numeric_limits<int>::max()) {
            potential[v] = potential[u] + d;
        }

        if (!valid) {
            valid = true;
            return {};
        }

        changed.clear();
        changed.resize(potential.size(), false);

        PairPrioQueue gamma;
        std::vector<int> last_edges(potential.size());

        int v_gamma = potential[u] + d - potential[v];
        if (v_gamma >= 0) {
            return {};
        } else {
            last_edges[v] = edge;
        }

        std::vector<int> gamma_vec(potential.size(), 0);
        gamma.push(std::make_pair(v, v_gamma));
        gamma_vec[v] = v_gamma;

        bool inconsistent = false;

        while (!gamma.empty() && gamma.top().second < 0) {
            int s = gamma.top().first;
            if (s == u) {
                inconsistent = true;
                break;
            }
            potential[s] += gamma.top().second;
            changed[s] = true;
            gamma.pop();
            gamma_vec[s] = 0;
            if (!graph.outgoing[s].empty()) {
                for (auto eid : graph.outgoing[s]) {
                    assert(eid < numeric_cast<int>(graph.edges.size()));
                    int t = graph.edges[eid].to;
                    if (!changed[t]) {
                        int weight = graph.edges[eid].weight;
                        int t_gamma = std::min(gamma_vec[t], potential[s] + weight - potential[t]);
                        gamma_vec[t] = t_gamma;
                        if (t_gamma < 0) {
                            gamma.push(std::make_pair(t, t_gamma));
                            last_edges[t] = eid;
                        }
                    }
                }
            }
        }

        if (inconsistent) {
            std::vector<int> neg_cycle;
            int begin = v;
            neg_cycle.push_back(graph.edges[last_edges[v]].id);
            int next = graph.edges[last_edges[v]].from;
            while (begin != next) {
                neg_cycle.push_back(graph.edges[last_edges[next]].id);
                next = graph.edges[last_edges[next]].from;
            }
            return neg_cycle;
        }

        return {};
    }

    void reset() {
        this->valid = false;
        this->graph.reset();
        this->potential.clear();
    }
};

class DifferenceLogicPropagator : public Propagator {

  private:
    struct StackItem {
        StackItem(uint32_t dl, int si) : decision_level(dl), trail_index(si) {}
        uint32_t decision_level;
        int trail_index;
    };

    struct State {
        State(const std::vector<Edge> &edges) : dl_graph(edges) {}
        std::vector<StackItem> stack;
        std::vector<int> edge_trail;
        DifferenceLogicGraph dl_graph;
    };

    std::vector<State> states;
    std::vector<literal_t> edges_to_lit;
    std::unordered_map<literal_t, std::vector<int>> lit_to_edges;
    std::vector<Edge> edges;
    std::vector<std::string> vert_map;
    std::unordered_map<std::string, int> vert_map_inv;

    int map_vert(std::string v) {
        auto it = std::find(vert_map.begin(), vert_map.end(), v);
        if (it != vert_map.end()) {
            return it - vert_map.begin();
        }
        vert_map.emplace_back(v);
        vert_map_inv[v] = vert_map.size() - 1;
        return vert_map.size() - 1;
    }

    void add_edge_atom(PropagateInit &init, TheoryAtom const &atom) {
        int lit = init.solver_literal(atom.literal());
        int weight = 0;
        std::string u, v;
        if (atom.guard().second.arguments().empty()) { // Check if constant is  negated
            weight = atom.guard().second.number();
        } else {
            weight = -atom.guard().second.arguments()[0].number();
        }
        u = atom.elements()[0].tuple()[0].arguments()[0].to_string();
        v = atom.elements()[0].tuple()[0].arguments()[1].to_string();
        int u_id = map_vert(u);
        int v_id = map_vert(v);
        Edge edge = Edge(edges.size(), u_id, v_id, weight);
        edges.emplace_back(edge);
        edges_to_lit.emplace_back(lit);
        lit_to_edges[lit].emplace_back(edge.id);
        init.add_watch(lit);
    }

    void initialize_states(PropagateInit &init) {
        for (int i = 0; i < init.number_of_threads(); ++i) {
            states.emplace_back(edges);
        }
    }

    std::vector<int> check_neg_cycle(State &state, int offset) {
        int eid;
        if (state.dl_graph.is_valid()) {
            eid = state.stack.back().trail_index + offset;
        } else {
            eid = 0;
        }

        std::vector<int> neg_cycle;
        for (int n = eid; n < numeric_cast<int>(state.edge_trail.size()); n++) {
            for (auto &id : lit_to_edges[state.edge_trail[n]]) {
                assert(id < numeric_cast<int>(edges.size()));
                // assert(id < state.dl_graph.graph.edges.size());
                neg_cycle = state.dl_graph.add_edge(id);
                if (!neg_cycle.empty()) {
                    return neg_cycle;
                }
            }
        }

        return {};
    }

    bool check_consistency(PropagateControl &ctl, State &state, int offset) {
        std::vector<int> neg_cycle = check_neg_cycle(state, offset);
        if (!neg_cycle.empty()) {
            std::vector<literal_t> clause;
            for (auto eid : neg_cycle) {
                clause.emplace_back(-edges_to_lit[eid]);
            }
            return ctl.add_clause(clause) && ctl.propagate();
        }
        return true;
    }

  public:
    DifferenceLogicPropagator() {}

    void init(PropagateInit &init) override {
        for (auto atom : init.theory_atoms()) {
            auto term = atom.term();
            if (term.to_string() == "diff") {
                add_edge_atom(init, atom);
            }
        }

        initialize_states(init);
    }

    void propagate(PropagateControl &ctl, LiteralSpan changes) override {
        auto &state = states[ctl.thread_id()];
        uint32_t dl = ctl.assignment().decision_level();
        uint32_t old_dl = 0;
        if (!state.stack.empty()) {
            old_dl = state.stack.back().decision_level;
        }
        if (state.stack.empty() || old_dl < dl) {
            state.stack.emplace_back(dl, static_cast<int>(state.edge_trail.size()));
        }
        for (auto &lit : changes) {
            state.edge_trail.emplace_back(lit);
        }
        int offset = 0;
        if (!state.stack.empty() && old_dl == dl) {
            offset = state.edge_trail.size() - state.stack.back().trail_index;
        }
        check_consistency(ctl, state, offset);

        return;
    }

    void undo(PropagateControl const &ctl, LiteralSpan changes) override {
        static_cast<void>(changes);
        auto &state = states[ctl.thread_id()];
        int sid = state.stack.back().trail_index;
        auto ib = state.edge_trail.begin() + sid, ie = state.edge_trail.end();
        state.edge_trail.erase(ib, ie);
        state.stack.pop_back();
        state.dl_graph.reset();
    }

    void check(PropagateControl &ctl) override {
        auto &state = states[ctl.thread_id()];
        std::cout << "Valid assignment found:" << std::endl;
        std::unordered_map<int, int> assignment = state.dl_graph.get_assignment();
        for (auto &it : assignment) {
            if (vert_map[it.first] != "0") {
                std::cout << vert_map[it.first] << ":" << it.second << " ";
            }
        }
        std::cout << std::endl;
    }
};

int get_int(std::string constname, Control &ctl, int def) {
    Symbol val = ctl.get_const(constname.c_str());
    if (val.to_string() == constname.c_str()) {
        return def;
    }
    return val.number();
}

int main(int argc, char *argv[]) {
    Control ctl{{argv + 1, static_cast<size_t>(argc - 2)}, nullptr, 20};
    // TODO: configure strict/non-strict mode
    int c = 0;
    if (argc > 1) {
        std::string s = argv[argc - 1];
        size_t pos = 0;
        std::string token;
        while ((pos = s.find(' ')) != std::string::npos) {
            token = s.substr(0, pos);
            ctl.load(token.c_str());
            s.erase(0, pos + 1);
        }
        ctl.load(s.c_str());
        c = get_int("strict", ctl, 0);
    }

    DifferenceLogicPropagator p;
    ctl.register_propagator(p, false);
    ctl.ground({{"base", {}}}, nullptr);
    int i = 0;
    for (auto m : ctl.solve()) {
        i++;
        std::cout << "Answer " << i << std::endl;
        std::cout << m << std::endl << std::endl;
    }
    if (i == 0) {
        std::cout << "UNSATISFIABLE " << std::endl;
    } else {
        std::cout << "SATISFIABLE " << std::endl;
    }
}
