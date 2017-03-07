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
    Edge(int x, int y, int k) :from(x), to(y), weight(k) {}
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

template <class C>
void ensure_index(C &c, size_t index) {
    if (index >= c.size()) {
        c.resize(index + 1);
    }
}

constexpr int undefined_potential = std::numeric_limits<int>::max();
struct DifferenceLogicNode {
    std::vector<int> outgoing;
    int potential = undefined_potential;
    int last_edge = 0;
    int gamma = 0;
    bool changed = false;
};

class DLPairComp {
public:
    // Note: reverse == true is never used
    DLPairComp(bool reverse = false) : reverse(reverse) { }
    bool operator()(const std::pair<int, int> &lhs, const std::pair<int, int> &rhs) const {
        return reverse ? lhs.second < rhs.second : lhs.second > rhs.second;
    }
private:
    bool reverse;
};

class DifferenceLogicGraph {
public:
    DifferenceLogicGraph(const std::vector<Edge> &edges) : edges_(edges) {}

    bool is_valid() { return valid_; }

    std::unordered_map<int, int> get_assignment() {
        std::unordered_map<int, int> ass;
        int idx = 0;
        for (auto &&node : nodes_) {
            if (node.potential != undefined_potential) {
                ass[idx] = -node.potential;
            }
            idx++;
        }
        return ass;
    }

    std::vector<int> add_edge(int uv_idx) {
        auto &&uv = edges_[uv_idx];

        ensure_index(nodes_, std::max(uv.from, uv.to));
        auto &&u = nodes_[uv.from], &&v = nodes_[uv.to];
        if (u.potential == undefined_potential) {
            u.potential = 0;
        }
        u.outgoing.emplace_back(uv_idx);
        if (v.potential == undefined_potential) {
            v.potential = u.potential + uv.weight;
        }
        if (!valid_) {
            valid_ = true;
            return {};
        }

        v.gamma = u.potential + uv.weight - v.potential;
        if (v.gamma >= 0) {
            return {};
        } else {
            v.last_edge = uv_idx;
        }

        while (!gamma_.empty()) { gamma_.pop(); }
        if (v.gamma < 0) {
            gamma_.push(std::make_pair(uv.to, v.gamma));
        }

        bool inconsistent = false;

        while (!gamma_.empty()) {
            int s_idx = gamma_.top().first;
            if (s_idx == uv.from) {
                inconsistent = true;
                break;
            }
            auto &&s = nodes_[s_idx];
            s.potential += gamma_.top().second;
            s.changed = true;
            changed_.emplace_back(s_idx);
            gamma_.pop();
            s.gamma = 0;
            for (auto st_idx : s.outgoing) {
                assert(st_idx < numeric_cast<int>(edges_.size()));
                auto &&st = edges_[st_idx];
                auto &&t = nodes_[st.to];
                if (!t.changed) {
                    t.gamma = std::min(t.gamma, s.potential + st.weight - t.potential);
                    if (t.gamma < 0) {
                        gamma_.push(std::make_pair(st.to, t.gamma));
                        t.last_edge = st_idx;
                    }
                }
            }
        }
        for (auto x : changed_) {
            nodes_[x].changed = false;
        }
        changed_.clear();

        std::vector<int> neg_cycle;
        if (inconsistent) {
            neg_cycle.push_back(v.last_edge);
            auto next_idx = edges_[v.last_edge].from;
            while (uv.to != next_idx) {
                auto next = nodes_[next_idx];
                neg_cycle.push_back(next.last_edge);
                next_idx = edges_[next.last_edge].from;
            }
        }
        return neg_cycle;
    }

    void reset() {
        valid_ = false;
        nodes_.clear();
    }

private:
    using PairPrioQueue = std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, DLPairComp>;

    PairPrioQueue gamma_;
    std::vector<int> changed_;
    const std::vector<Edge> &edges_;
    std::vector<DifferenceLogicNode> nodes_;
    bool valid_ = false;
};

struct DLStackItem {
    DLStackItem(uint32_t dl, int si) : decision_level(dl), trail_index(si) {}
    uint32_t decision_level;
    int trail_index;
};

struct DLState {
    DLState(const std::vector<Edge> &edges) : dl_graph(edges) {}
    std::vector<DLStackItem> stack;
    std::vector<int> edge_trail;
    DifferenceLogicGraph dl_graph;
};

class DifferenceLogicPropagator : public Propagator {
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
        int id = edges.size();
        Edge edge = Edge(u_id, v_id, weight);
        edges.emplace_back(edge);
        edges_to_lit.emplace_back(lit);
        lit_to_edges[lit].emplace_back(id);
        init.add_watch(lit);
    }

    void initialize_states(PropagateInit &init) {
        for (int i = 0; i < init.number_of_threads(); ++i) {
            states.emplace_back(edges);
        }
    }

    std::vector<int> check_neg_cycle(DLState &state, int offset) {
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

    bool check_consistency(PropagateControl &ctl, DLState &state, int offset) {
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

private:
    std::vector<DLState> states;
    std::vector<literal_t> edges_to_lit;
    std::unordered_map<literal_t, std::vector<int>> lit_to_edges;
    std::vector<Edge> edges;
    std::vector<std::string> vert_map;
    std::unordered_map<std::string, int> vert_map_inv;
};

int get_int(std::string constname, Control &ctl, int def) {
    Symbol val = ctl.get_const(constname.c_str());
    if (val.to_string() == constname.c_str()) {
        return def;
    }
    return val.number();
}

int main(int argc, char *argv[]) {
    char **argv2 = argv + 1;
    int argc2 = argc - 1;
    for (; *argv2; ++argv2, --argc2) {
        if (std::strcmp(*argv2, "--") == 0) {
            ++argv2;
            --argc2;
            break;
        }
    }
    Control ctl{{argv2, numeric_cast<size_t>(argc2)}, nullptr, 20};
    for (char **arg = argv + 1; arg != argv + (argc - argc2 - 1); ++arg) {
        ctl.load(*arg);
    }
    // TODO: configure strict/non-strict mode
    // int c = get_int("strict", ctl, 0);

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
