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
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <limits>
#include <chrono>
#include <iomanip>
#include <stdlib.h>

//#define CROSSCHECK
#define CHECKSOLUTION

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

// NOTE: according to papers 4-ary heaps perform best for dijkstra and are quite easy to implement
class BinaryHeap {
public:
    template <class M>
    void push(M &m, int item) {
        auto i = m.offset(item) = heap_.size();
        heap_.push_back(item);
        decrease(m, i);
    }
    template <class M>
    int pop(M &m) {
        assert(!heap_.empty());
        auto ret = heap_[0];
        if (heap_.size() > 1) {
            heap_[0] = heap_.back();
            m.offset(heap_[0]) = 0;
            heap_.pop_back();
            increase(m, 0);
        }
        else { heap_.pop_back(); }
        return ret;
    }

    template <class M>
    void decrease(M &m, int i) {
        while (i > 0) {
            int p = parent_(i);
            if (m.cost(heap_[p]) > m.cost(heap_[i])) {
                swap_(m, i, p);
                i = p;
            }
            else { break; }
        }
    }
    template <class M>
    void increase(M &m, int i) {
        for (int p = i, l = left_(p), s = numeric_cast<int>(heap_.size()); l < s; l = left_(p)) {
            int r = right_(p);
            if (r < s) {
                if (less_(m, r, l)) {
                    if (less_(m, r, p)) {
                        swap_(m, p, r);
                        p = r;
                        continue;
                    }
                }
            }
            if (less_(m, l, p)) {
                swap_(m, p, l);
                p = l;
                continue;

            }
            break;
        }
    }
    int size() {
        return heap_.size();
    }
    bool empty() {
        return heap_.empty();
    }
    void clear() {
        heap_.clear();
    }
private:
    template <class M>
    void swap_(M &m, int i, int j) {
        m.offset(heap_[j]) = i;
        m.offset(heap_[i]) = j;
        std::swap(heap_[i], heap_[j]);
    }
    int left_(int offset) {
        return 2 * offset + 1;
    }
    int right_(int offset) {
        return left_(offset) + 1;
    }
    int parent_(int offset) {
        return (offset - 1) / 2;
    }
    template <class M>
    bool less_(M &m, int a, int b) {
        a = heap_[a], b = heap_[b];
        auto ca = m.cost(a), cb = m.cost(b);
        return ca < cb || (ca == cb && m.relevant(a) < m.relevant(b));
    }

private:
    std::vector<int> heap_;
};

template <class T, class P>
struct HeapFromM {
    int &offset(int idx)                          { return static_cast<P*>(this)->nodes_[idx].offset; }
    T &cost(int idx)                              { return static_cast<P*>(this)->nodes_[idx].cost_from; }
    int to(int idx)                               { return static_cast<P*>(this)->edges_[idx].to; }
    int from(int idx)                             { return static_cast<P*>(this)->edges_[idx].from; }
    std::vector<int> &out(int idx)                { return static_cast<P*>(this)->nodes_[idx].outgoing; }
    int &path(int idx)                            { return static_cast<P*>(this)->nodes_[idx].path_from; }
    bool &visited(int idx)                        { return static_cast<P*>(this)->nodes_[idx].visited_from; }
    bool &relevant(int idx)                       { return static_cast<P*>(this)->nodes_[idx].relevant_from; }
    std::vector<int> &visited_set()               { return static_cast<P*>(this)->visited_from_; }
    std::vector<int> &candidate_outgoing(int idx) { return static_cast<P*>(this)->nodes_[idx].candidate_outgoing; }
    std::vector<int> &candidate_incoming(int idx) { return static_cast<P*>(this)->nodes_[idx].candidate_incoming; }
    void remove_incoming(int idx)                 { static_cast<P*>(this)->edge_states_[idx].removed_incoming = true; }
    void remove_outgoing(int idx)                 { static_cast<P*>(this)->edge_states_[idx].removed_outgoing = true; }
};

template <class T, class P>
struct HeapToM {
    int &offset(int idx)                          { return static_cast<P*>(this)->nodes_[idx].offset; }
    T &cost(int idx)                              { return static_cast<P*>(this)->nodes_[idx].cost_to; }
    int to(int idx)                               { return static_cast<P*>(this)->edges_[idx].from; }
    int from(int idx)                             { return static_cast<P*>(this)->edges_[idx].to; }
    std::vector<int> &out(int idx)                { return static_cast<P*>(this)->nodes_[idx].incoming; }
    int &path(int idx)                            { return static_cast<P*>(this)->nodes_[idx].path_to; }
    bool &visited(int idx)                        { return static_cast<P*>(this)->nodes_[idx].visited_to; }
    bool &relevant(int idx)                       { return static_cast<P*>(this)->nodes_[idx].relevant_to; }
    std::vector<int> &visited_set()               { return static_cast<P*>(this)->visited_to_; }
    std::vector<int> &candidate_outgoing(int idx) { return static_cast<P*>(this)->nodes_[idx].candidate_incoming; }
    std::vector<int> &candidate_incoming(int idx) { return static_cast<P*>(this)->nodes_[idx].candidate_outgoing; }
    void remove_incoming(int idx)                 { static_cast<P*>(this)->edge_states_[idx].removed_outgoing = true; }
    void remove_outgoing(int idx)                 { static_cast<P*>(this)->edge_states_[idx].removed_incoming = true; }
};

template <typename T>
struct DifferenceLogicNode {
    bool defined() const { return !potential_stack.empty(); }
    T potential() const { return potential_stack.back().second; }
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
    int degree = 0;
    bool relevant_from = false;
    bool relevant_to = false;
    bool visited_from = false;
    bool visited_to = false;
};

struct DLStats {
    Duration time_propagate = Duration{0};
    Duration time_undo = Duration{0};
    Duration time_dijkstra = Duration{0};
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
    friend struct HeapToM<T, DifferenceLogicGraph<T>>;
    friend struct HeapFromM<T, DifferenceLogicGraph<T>>;
public:
    DifferenceLogicGraph(DLStats &stats, const std::vector<Edge<T>> &edges)
        : edges_(edges)
        , stats_(stats) {
        edge_states_.resize(edges_.size(), {1, 1, 1});
        for (int i = 0; i < numeric_cast<int>(edges_.size()); ++i) {
            ensure_index(nodes_, std::max(edges_[i].from, edges_[i].to));
            add_candidate_edge(i);
        }
    }

    bool empty() const { return nodes_.empty(); }

    int node_value_defined(int idx) const { return nodes_[idx].defined(); }
    T node_value(int idx) const { return -nodes_[idx].potential(); }

    bool edge_is_active(int edge_idx) const {
        return edge_states_[edge_idx].active;
    }

    void ensure_decision_level(int level) {
        // initialize the trail
        if (changed_trail_.empty() || std::get<0>(changed_trail_.back()) < level) {
            changed_trail_.emplace_back(level, changed_nodes_.size(), changed_edges_.size(), inactive_edges_.size());
        }
    }

    std::vector<int> add_edge(int uv_idx) {
#ifdef CROSSCHECK
        for (auto &node : nodes_) {
            assert(!node.visited_from);
        }
#endif
        assert(visited_from_.empty());
        assert(costs_heap_.empty());
        int level = current_decision_level_();
        auto &uv = edges_[uv_idx];
        // NOTE: would be more efficient if relevant would return statically false here
        //       for the compiler to make comparison cheaper
        auto &m = *static_cast<HFM*>(this);

        // initialize the nodes of the edge to add
        auto &u = nodes_[uv.from];
        auto &v = nodes_[uv.to];
        if (!u.defined()) { set_potential(u, level, 0); }
        if (!v.defined()) { set_potential(v, level, 0); }
        v.cost_from = u.potential() + uv.weight - v.potential();
        if (v.cost_from < 0) {
            costs_heap_.push(m, uv.to);
            visited_from_.emplace_back(uv.to);
            v.visited_from = true;
            v.path_from = uv_idx;
        }

        // detect negative cycles
        while (!costs_heap_.empty() && !u.visited_from) {
            auto s_idx = costs_heap_.pop(m);
            auto &s = nodes_[s_idx];
            assert(s.visited_from);
            set_potential(s, level, s.potential() + s.cost_from);
            for (auto st_idx : s.outgoing) {
                assert(st_idx < numeric_cast<int>(edges_.size()));
                auto &st = edges_[st_idx];
                auto &t = nodes_[st.to];
                auto c = s.potential() + st.weight - t.potential();
                if (c < (t.visited_from ? t.cost_from : 0)) {
                    assert(c < 0);
                    t.path_from = st_idx;
                    t.cost_from = c;
                    if (!t.visited_from) {
                        t.visited_from = true;
                        visited_from_.emplace_back(st.to);
                        costs_heap_.push(m, st.to);
                    }
                    else { costs_heap_.decrease(m, m.offset(st.to)); }
                }
            }
        }

        std::vector<int> neg_cycle;
        if (u.visited_from) {
            // gather the edges in the negative cycle
            neg_cycle.push_back(v.path_from);
            auto next_idx = edges_[v.path_from].from;
            while (uv.to != next_idx) {
                auto &next = nodes_[next_idx];
                neg_cycle.push_back(next.path_from);
                next_idx = edges_[next.path_from].from;
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
            // NOTE: just a check that will throw if there is a cycle
            bellman_ford(changed_edges_, uv.from);
#endif
        }

        // reset visited flags
        for (auto &x : visited_from_) { nodes_[x].visited_from = false; }
        visited_from_.clear();
        costs_heap_.clear();

        return neg_cycle;
    }

    bool propagate(int xy_idx, Clingo::PropagateControl &ctl) {
        static int props = 0;
        ++props;
        remove_candidate_edge(xy_idx);
        auto &xy = edges_[xy_idx];
        auto &x = nodes_[xy.from];
        auto &y = nodes_[xy.to];
        x.relevant_to = true;
        y.relevant_from = true;
        int relevant_from;
        int relevant_to;
        {
            Timer t{stats_.time_dijkstra};
            relevant_from = dijkstra(xy.from, visited_from_, *static_cast<HFM*>(this));
            relevant_to = dijkstra(xy.to, visited_to_, *static_cast<HTM*>(this));
        }

        bool ret = relevant_from < relevant_to
            ? propagate_edges(*static_cast<HFM*>(this), ctl, xy_idx)
            : propagate_edges(*static_cast<HTM*>(this), ctl, xy_idx);

        for (auto &x : visited_from_) {
            nodes_[x].visited_from = false;
            nodes_[x].relevant_from = false;
        }
        for (auto &x : visited_to_) {
            nodes_[x].visited_to = false;
            nodes_[x].relevant_to = false;
        }
        visited_from_.clear();
        visited_to_.clear();
        return ret;
    }

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
        int n = std::get<3>(changed_trail_.back());
        for (auto i = inactive_edges_.begin() + n, e = inactive_edges_.end(); i < e; ++i) {
            add_candidate_edge(*i);
        }
        inactive_edges_.resize(n);
        changed_trail_.pop_back();
    }

private:
    void add_candidate_edge(int uv_idx) {
        auto &uv = edges_[uv_idx];
        auto &u = nodes_[uv.from];
        auto &v = nodes_[uv.to];
        ++u.degree; ++v.degree;
        edge_states_[uv_idx].active = true;
        if (edge_states_[uv_idx].removed_outgoing) {
            edge_states_[uv_idx].removed_outgoing = false;
            u.candidate_outgoing.emplace_back(uv_idx);
        }
        if (edge_states_[uv_idx].removed_incoming) {
            edge_states_[uv_idx].removed_incoming = false;
            v.candidate_incoming.emplace_back(uv_idx);
        }
    }

    void remove_candidate_edge(int uv_idx) {
        auto &uv = edges_[uv_idx];
        auto &u = nodes_[uv.from];
        auto &v = nodes_[uv.to];
        --u.degree;
        --v.degree;
        inactive_edges_.push_back(uv_idx);
        edge_states_[uv_idx].active = false;
    }

    bool propagate_edge_true(int uv_idx, int xy_idx) {
        auto &uv = edges_[uv_idx];
        auto &u = nodes_[uv.from];
        auto &v = nodes_[uv.to];
        assert(u.relevant_to || v.relevant_from);

        if (u.relevant_to && v.relevant_from) {
            auto &xy = edges_[xy_idx];
            auto &x = nodes_[xy.from];
            auto &y = nodes_[xy.to];

            auto a = u.cost_to + y.potential() - u.potential();
            auto b = v.cost_from + v.potential() - x.potential();
            auto d = a + b - xy.weight;
#ifdef CROSSCHECK
            auto bf_costs_from_u = bellman_ford(changed_edges_, uv.from);
            auto bf_costs_from_x = bellman_ford(changed_edges_, xy.from);
            auto aa = bf_costs_from_u.find(xy.to);
            assert (aa != bf_costs_from_u.end());
            assert(aa->second == a);
            auto bb = bf_costs_from_x.find(uv.to);
            assert (bb != bf_costs_from_u.end());
            assert(bb->second == b);
#endif
            if (d <= uv.weight) {
#ifdef CROSSCHECK
                auto edges = changed_edges_;
                edges.emplace_back(uv_idx);
                // NOTE: throws if there is a cycle
                try         { bellman_ford(changed_edges_, uv.from); }
                catch (...) { assert(false && "edge is implied but lead to a conflict :("); }
#endif
                remove_candidate_edge(uv_idx);
                return true;
            }
        }
        return false;
    }

    bool propagate_edge_false(Clingo::PropagateControl &ctl, int uv_idx, int xy_idx, bool &ret) {
        auto &uv = edges_[uv_idx];
        auto &u = nodes_[uv.from];
        auto &v = nodes_[uv.to];
        assert (v.relevant_to || u.relevant_from);

        if (v.relevant_to && u.relevant_from) {
            auto &xy = edges_[xy_idx];
            auto &x = nodes_[xy.from];
            auto &y = nodes_[xy.to];

            auto a = v.cost_to + y.potential() - v.potential();
            auto b = u.cost_from + u.potential() - x.potential();
            auto d = a + b - xy.weight;
            if (d < -uv.weight) {
                if (!ctl.assignment().is_false(uv.lit)) {
#ifdef CROSSCHECK
                    T sum = uv.weight - xy.weight;
#endif
                    std::vector<literal_t> clause;
                    clause.push_back(-uv.lit);
                    // forward
                    for (auto next_edge_idx = u.path_from; next_edge_idx >= 0; ) {
                        auto &next_edge = edges_[next_edge_idx];
                        auto &next_node = nodes_[next_edge.from];
                        clause.push_back(-next_edge.lit);
#ifdef CROSSCHECK
                        sum+= next_edge.weight;
#endif
                        next_edge_idx = next_node.path_from;
                    }
                    // backward
                    for (auto next_edge_idx = v.path_to; next_edge_idx >= 0; ) {
                        auto &next_edge = edges_[next_edge_idx];
                        auto &next_node = nodes_[next_edge.to];
                        clause.push_back(-next_edge.lit);
#ifdef CROSSCHECK
                        sum+= next_edge.weight;
#endif
                        next_edge_idx = next_node.path_to;
                    }
#ifdef CROSSCHECK
                    assert(sum < 0);
#endif
                    if (!(ret = ctl.add_clause(clause) && ctl.propagate())) {
                        return false;
                    }
                }
                remove_candidate_edge(uv_idx);
                return true;
            }
#ifdef CROSSCHECK
            else {
                auto edges = changed_edges_;
                edges.emplace_back(uv_idx);
                // NOTE: throws if there is a cycle
                try         { bellman_ford(changed_edges_, uv.from); }
                catch (...) { assert(false && "edge must not cause a conflict"); }
            }
#endif
        }
        return false;
    }

    template <class M>
    bool propagate_edges(M &m, Clingo::PropagateControl &ctl, int xy_idx) {
        for (auto &node : m.visited_set()) {
            if (m.relevant(node)) {
                auto &in = m.candidate_incoming(node);
                in.resize(std::remove_if(in.begin(), in.end(), [&](int uv_idx) {
                    if (!edge_states_[uv_idx].active || propagate_edge_true(uv_idx, xy_idx)) {
                        m.remove_incoming(uv_idx);
                        return true;
                    }
                    return false;
                }) - in.begin());
                bool ret = true;
                auto &out = m.candidate_outgoing(node);
                out.resize(std::remove_if(out.begin(), out.end(), [&](int uv_idx) {
                    if (!ret) { return false; }
                    if (!edge_states_[uv_idx].active || propagate_edge_false(ctl, uv_idx, xy_idx, ret)) {
                        m.remove_outgoing(uv_idx);
                        return true;
                    }
                    return false;
                }) - out.begin());
                if (!ret) { return false; }
            }
        }
        return true;
    }

    template <class M>
    int dijkstra(int source_idx, std::vector<int> &visited_set, M &m) {
        // TODO: the paper argues that the SSSP algorithm in "Shortest path algorithms: Engineering aspects" is faster
        //       than a simple dijkstra. Maybe there are even better ones nowadays.
        int relevant = 0;
        int relevant_degree = 0;
        assert(visited_set.empty() && costs_heap_.empty());
        costs_heap_.push(m, source_idx);
        visited_set.push_back(source_idx);
        m.visited(source_idx) = true;
        m.cost(source_idx) = 0;
        m.path(source_idx) = -1;
        while (!costs_heap_.empty()) {
            auto u_idx = costs_heap_.pop(m);
            auto tu = m.path(u_idx);
            if (tu >= 0 && m.relevant(m.from(tu))) {
                m.relevant(u_idx) = true;
                --relevant; // just removed a relevant edge from the queue
            }
            bool relevant_u = m.relevant(u_idx);
            if (relevant_u) { relevant_degree += nodes_[u_idx].degree; }
            for (auto &uv_idx : m.out(u_idx)) {
                auto &uv = edges_[uv_idx];
                auto v_idx = m.to(uv_idx);
                // NOTE: explicitely using uv.from and uv.to is intended here
                auto c = m.cost(u_idx) + nodes_[uv.from].potential() + uv.weight - nodes_[uv.to].potential();
                assert(nodes_[uv.from].potential() + uv.weight - nodes_[uv.to].potential() >= 0);
                if (!m.visited(v_idx) || c < m.cost(v_idx)) {
                    m.cost(v_idx) = c;
                    if (!m.visited(v_idx)) {
                        // node v contributes an edge with a relevant source
                        if (relevant_u) { ++relevant; }
                        visited_set.push_back(m.to(uv_idx));
                        m.visited(v_idx) = true;
                        costs_heap_.push(m, v_idx);
                    }
                    else {
                        if (m.relevant(m.from(m.path(v_idx)))) {
                            // node v no longer contributes a relevant edge
                            if (!relevant_u) { --relevant; }
                        }
                        // node v contributes a relevant edge now
                        else if (relevant_u) { ++relevant; }
                        costs_heap_.decrease(m, m.offset(v_idx));
                    }
                    m.path(v_idx) = uv_idx;
                }
            }
            // removed a relevant node from the queue and there are no edges with relevant sources anymore in the queue
            // this condition assumes that initially there is exactly one reachable relevant node in the graph
            if (relevant_u && relevant == 0) {
                costs_heap_.clear();
                break;
            }
        }
        return relevant_degree;
    }

#ifdef CROSSCHECK
    std::unordered_map<int, T> bellman_ford(std::vector<int> const &edges, int source) {
        std::unordered_map<int, T> costs;
        costs[source] = 0;
        int nodes = 0;
        for (auto &node : nodes_) {
            if (node.defined()) {
                ++nodes;
            }
        }
        for (int i = 0; i < nodes; ++i) {
            for (auto &uv_idx : edges) {
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
        for (auto &uv_idx : edges) {
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
        return costs;
    }
#endif

    void set_potential(DifferenceLogicNode<T> &node, int level, T potential) {
        if (!node.defined() || node.potential_stack.back().first < level) {
            node.potential_stack.emplace_back(level, potential);
            changed_nodes_.emplace_back(numeric_cast<int>(&node - nodes_.data()));
        }
        else {
            node.potential_stack.back().second = potential;
        }
    }

    int current_decision_level_() {
        assert(!changed_trail_.empty());
        return std::get<0>(changed_trail_.back());
    }

private:
    BinaryHeap costs_heap_;
    std::vector<int> visited_from_;
    std::vector<int> visited_to_;
    std::vector<Edge<T>> const &edges_;
    std::vector<DifferenceLogicNode<T>> nodes_;
    std::vector<int> changed_nodes_;
    std::vector<int> changed_edges_;
    std::vector<std::tuple<int, int, int, int>> changed_trail_;
    std::vector<int> inactive_edges_;
    std::vector<EdgeState> edge_states_;
    DLStats &stats_;
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
        , dl_graph(stats, edges) {}
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
    DifferenceLogicPropagator(Stats &stats, bool strict, bool propagate)
        : stats_(stats)
        , strict_(strict)
        , propagate_(propagate) { }

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
        state.dl_graph.ensure_decision_level(level);
        // NOTE: vector copy only because clasp bug
        //       can only be triggered with propagation
        for (auto lit : std::vector<Clingo::literal_t>(changes.begin(), changes.end())) {
            for (auto it = lit_to_edges_.find(lit), ie = lit_to_edges_.end(); it != ie && it->first == lit; ++it) {
                if (state.dl_graph.edge_is_active(it->second)) {
                    auto neg_cycle = state.dl_graph.add_edge(it->second);
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
                    else if (propagate_) {
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

#if defined(CHECKSOLUTION) || defined(CROSSCHECK)
    void check(PropagateControl &ctl) override {
        auto &state = states_[ctl.thread_id()];
        for (auto &x : edges_) {
            if (ctl.assignment().is_true(x.lit)) {
                if (!state.dl_graph.node_value_defined(x.from) ||
                    !state.dl_graph.node_value_defined(x.to) ||
                    !(state.dl_graph.node_value(x.from) - state.dl_graph.node_value(x.to) <= x.weight)) {
                    throw std::logic_error("not a valid solution");
                }
            }
        }
    }
#endif

private:
    std::vector<DLState<T>> states_;
    std::unordered_multimap<literal_t, int> lit_to_edges_;
    std::vector<Edge<T>> edges_;
    std::vector<std::reference_wrapper<const std::string>> vert_map_;
    std::unordered_map<std::string, int> vert_map_inv_;
    Stats &stats_;
    bool strict_;
    bool propagate_;
};

int get_int(std::string constname, Control &ctl, int def) {
    Symbol val = ctl.get_const(constname.c_str());
    if (val.to_string() == constname.c_str()) {
        return def;
    }
    return val.number();
}

template <typename T>
void solve(Stats &stats, Control &ctl, bool strict, bool propagate) {
    DifferenceLogicPropagator<T> p{stats, strict, propagate};
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
        bool propagate = false;
        for (auto arg = argv + 1; arg != arge; ++arg) {
            if (std::strcmp(*arg, "-p") == 0) {
                propagate = true;
            }
            else {
                ctl.load(*arg);
            }
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
            solve<double>(stats, ctl, strict, propagate);
        }
        else {
            solve<int>(stats, ctl, strict, propagate);
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
        std::cout << "      dijkstra: " << stat.time_dijkstra.count() << "s\n";
        std::cout << "    undo     : " << stat.time_undo.count() << "s\n";
        ++thread;
    }
    std::cout << "conflicts: " << stats.conflicts << "\n";
    std::cout << "choices  : " << stats.choices << "\n";
    std::cout << "restarts : " << stats.restarts << "\n";
}
