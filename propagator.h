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

#ifndef PROPAGATOR_H
#define PROPAGATOR_H


#include <clingo.hh>
#include <unordered_map>
#include <chrono>

//#define CROSSCHECK
#define CHECKSOLUTION
using namespace Clingo;
using std::chrono::steady_clock;

namespace std {

template<>
struct hash<std::pair<int, int>> {
    size_t operator()(std::pair<int, int> const &p) const {
        return static_cast<size_t>(p.first) + (static_cast<size_t>(p.second) >> 32);
    }
};

}

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

template <int N>
class Heap {
public:
    template <class M>
    void push(M &m, int item) {
        auto i = m.offset(item) = static_cast<int>(heap_.size());
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
        else {
            heap_.pop_back();
        }
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
            else {
                break;
            }
        }
    }
    template <class M>
    void increase(M &m, int i) {
        for (int p = i, j = children_(p), s = numeric_cast<int>(heap_.size()); j < s; j = children_(p)) {
            int min = j;
            for (int k = j + 1; k < j + N; ++k) {
                if (k < s && less_(m, k, min)) {
                    min = k;
                }
            }
            if (less_(m, min, p)) {
                swap_(m, min, p);
                p = min;
            }
            else {
                return;
            }
        }
    }
    int size() { return heap_.size(); }
    bool empty() { return heap_.empty(); }
    void clear() { heap_.clear(); }

private:
    template <class M>
    void swap_(M &m, int i, int j) {
        m.offset(heap_[j]) = i;
        m.offset(heap_[i]) = j;
        std::swap(heap_[i], heap_[j]);
    }
    int parent_(int offset) { return (offset - 1) / N; }
    int children_(int offset) { return N * offset + 1; }
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
    int &offset(int idx) { return static_cast<P *>(this)->nodes_[idx].offset; }
    T &cost(int idx) { return static_cast<P *>(this)->nodes_[idx].cost_from; }
    int to(int idx) { return static_cast<P *>(this)->edges_[idx].to; }
    int from(int idx) { return static_cast<P *>(this)->edges_[idx].from; }
    std::vector<int> &out(int idx) { return static_cast<P *>(this)->nodes_[idx].outgoing; }
    int &path(int idx) { return static_cast<P *>(this)->nodes_[idx].path_from; }
    bool &visited(int idx) { return static_cast<P *>(this)->nodes_[idx].visited_from; }
    bool &relevant(int idx) { return static_cast<P *>(this)->nodes_[idx].relevant_from; }
    std::vector<int> &visited_set() { return static_cast<P *>(this)->visited_from_; }
    std::vector<int> &candidate_outgoing(int idx) { return static_cast<P *>(this)->nodes_[idx].candidate_outgoing; }
    std::vector<int> &candidate_incoming(int idx) { return static_cast<P *>(this)->nodes_[idx].candidate_incoming; }
    void remove_incoming(int idx) { static_cast<P *>(this)->edge_states_[idx].removed_incoming = true; }
    void remove_outgoing(int idx) { static_cast<P *>(this)->edge_states_[idx].removed_outgoing = true; }
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
    int degree_out = 0;
    int degree_in = 0;
    bool relevant_from = false;
    bool relevant_to = false;
    bool visited_from = false;
    bool visited_to = false;
};

struct DLStats {
    void reset() {
        time_propagate = steady_clock::duration::zero();
        time_undo      = steady_clock::duration::zero();
        time_dijkstra  = steady_clock::duration::zero();
        true_edges     = 0;
        false_edges    = 0;
    }
    Duration time_propagate = Duration{0};
    Duration time_undo = Duration{0};
    Duration time_dijkstra = Duration{0};
    uint64_t true_edges{0};
    uint64_t false_edges{0};
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
        edge_states_.resize(edges_.size(), {1, 1, 0});
        for (int i = 0; i < numeric_cast<int>(edges_.size()); ++i) {
            ensure_index(nodes_, std::max(edges_[i].from, edges_[i].to));
            add_candidate_edge(i);
        }
    }

    bool empty() const { return nodes_.empty(); }

    int node_value_defined(int idx) const { return nodes_[idx].defined(); }
    T node_value(int idx) const { return -nodes_[idx].potential(); }

    bool edge_is_active(int edge_idx) const { return edge_states_[edge_idx].active; }

    void ensure_decision_level(int level) {
        // initialize the trail
        if (changed_trail_.empty() || static_cast<int>(std::get<0>(changed_trail_.back())) < level) {
            changed_trail_.emplace_back(level, static_cast<int>(changed_nodes_.size()),
                                               static_cast<int>(changed_edges_.size()),
                                               static_cast<int>(inactive_edges_.size()));
        }
    }

    template <class F>
    bool cheap_propagate(int u_idx, F f) {
        // this function can be called during the dijkstra algorithm
        // a path together with its cost can be extracted from the state
        auto &u = nodes_[u_idx];
        auto &in = u.candidate_incoming;
        auto jt = in.begin();
        for (auto it = jt, ie = in.end(); it != ie; ++it) {
            auto vu_idx = *it;
            auto &vu = edges_[vu_idx];
            assert(u_idx == vu.to);
            auto v_idx = vu.from;
            auto &v = nodes_[v_idx];
            T weight = v.visited_from ? v.potential() - u.potential() : -vu.weight;
            if (!edge_states_[vu_idx].active) {
                edge_states_[vu_idx].removed_incoming = true;
            }
            else if (vu.weight + weight < 0) {
                edge_states_[vu_idx].removed_incoming = true;
                remove_candidate_edge(vu_idx);
                ++stats_.false_edges;
                T check = 0;
                auto t_idx = v_idx;
                std::vector<int> neg_cycle;
                while (u_idx != t_idx) {
                    auto &t = nodes_[t_idx];
                    auto &st = edges_[t.path_from];
                    neg_cycle.emplace_back(t.path_from);
                    t_idx = st.from;
                    check += st.weight;
                }
                assert(weight == check);
                neg_cycle.emplace_back(vu_idx);
                if (!f(neg_cycle)) {
                    in.erase(jt, it+1);
                    return false;
                }
            }
            else {
                *jt++ = *it;
            }
        }
        in.erase(jt, in.end());
        return true;
    }

    template <class F>
    bool add_edge(int uv_idx, F f) {
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
        auto &m = *static_cast<HFM *>(this);

        // initialize the nodes of the edge to add
        auto &u = nodes_[uv.from];
        auto &v = nodes_[uv.to];
        if (!u.defined()) {
            set_potential(u, level, 0);
        }
        if (!v.defined()) {
            set_potential(v, level, 0);
        }
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
            // NOTE: here we should have a cost from u to t
            //       and can check for an edge t to v
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
                    else {
                        costs_heap_.decrease(m, m.offset(st.to));
                    }
                }
            }
        }

        bool consistent = true;
        if (!u.visited_from) {
            // add the edge to the graph
            u.outgoing.emplace_back(uv_idx);
            v.incoming.emplace_back(uv_idx);
            changed_edges_.emplace_back(uv_idx);
#ifdef CROSSCHECK
            // NOTE: just a check that will throw if there is a cycle
            bellman_ford(changed_edges_, uv.from);
#endif
        }
        else {
            // gather the edges in the negative cycle
            std::vector<int> neg_cycle;
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
            consistent = f(neg_cycle);
        }

        consistent = consistent && cheap_propagate(uv.from, f);
        // reset visited flags
        for (auto &x : visited_from_) {
            nodes_[x].visited_from = false;
        }
        visited_from_.clear();
        costs_heap_.clear();

        return consistent;
    }

    bool propagate(int xy_idx, Clingo::PropagateControl &ctl) {
        remove_candidate_edge(xy_idx);
        auto &xy = edges_[xy_idx];
        auto &x = nodes_[xy.from];
        auto &y = nodes_[xy.to];
        // BUG: this test is not correct
        // if ((x.incoming.empty() && x.outgoing.size() == 1) || (y.outgoing.empty() && y.incoming.size() == 1)) {
        //    return true;
        //}
        x.relevant_to = true;
        y.relevant_from = true;
        int num_relevant_out_from;
        int num_relevant_in_from;
        int num_relevant_out_to;
        int num_relevant_in_to;
        {
            Timer t{stats_.time_dijkstra};
            std::tie(num_relevant_out_from, num_relevant_in_from) = dijkstra(xy.from, visited_from_, *static_cast<HFM *>(this));
            std::tie(num_relevant_out_to, num_relevant_in_to) = dijkstra(xy.to, visited_to_, *static_cast<HTM *>(this));
        }
#ifdef CROSSCHECK
        int check_relevant_out_from = 0, check_relevant_in_from = 0;
        for (auto &node : visited_from_) {
            if (nodes_[node].relevant_from) {
                for (auto &edge : nodes_[node].candidate_incoming) {
                    if (edge_states_[edge].active) {
                        ++check_relevant_in_from;
                    }
                }
                for (auto &edge : nodes_[node].candidate_outgoing) {
                    if (edge_states_[edge].active) {
                        ++check_relevant_out_from;
                    }
                }
            }
        }
        assert(num_relevant_out_from == check_relevant_out_from);
        assert(num_relevant_in_from == check_relevant_in_from);
        int check_relevant_out_to = 0, check_relevant_in_to = 0;
        for (auto &node : visited_to_) {
            if (nodes_[node].relevant_to) {
                for (auto &edge : nodes_[node].candidate_incoming) {
                    if (edge_states_[edge].active) {
                        ++check_relevant_in_to;
                    }
                }
                for (auto &edge : nodes_[node].candidate_outgoing) {
                    if (edge_states_[edge].active) {
                        ++check_relevant_out_to;
                    }
                }
            }
        }
        assert(num_relevant_out_to == check_relevant_out_to);
        assert(num_relevant_in_to == check_relevant_in_to);
#endif

        bool forward_from = num_relevant_in_from < num_relevant_out_to;
        bool backward_from = num_relevant_out_from < num_relevant_in_to;

        bool ret = propagate_edges(*static_cast<HFM *>(this), ctl, xy_idx, forward_from, backward_from) && propagate_edges(*static_cast<HTM *>(this), ctl, xy_idx, !forward_from, !backward_from);

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
        for (auto count = static_cast<int>(changed_nodes_.size()) - std::get<1>(changed_trail_.back()); count > 0; --count) {
            auto &node = nodes_[changed_nodes_.back()];
            node.potential_stack.pop_back();
            changed_nodes_.pop_back();
        }
        for (auto count = static_cast<int>(changed_edges_.size()) - std::get<2>(changed_trail_.back()); count > 0; --count) {
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

    void remove_candidate_edge(int uv_idx) {
        auto &uv = edges_[uv_idx];
        auto &u = nodes_[uv.from];
        auto &v = nodes_[uv.to];
        --u.degree_out;
        --v.degree_in;
        inactive_edges_.push_back(uv_idx);
        assert(edge_states_[uv_idx].active);
        edge_states_[uv_idx].active = false;
    }

private:
    void add_candidate_edge(int uv_idx) {
        auto &uv = edges_[uv_idx];
        auto &uv_state = edge_states_[uv_idx];
        auto &u = nodes_[uv.from];
        auto &v = nodes_[uv.to];
        ++u.degree_out;
        ++v.degree_in;
        assert(!uv_state.active);
        uv_state.active = true;
        if (uv_state.removed_outgoing) {
            uv_state.removed_outgoing = false;
            u.candidate_outgoing.emplace_back(uv_idx);
        }
        if (uv_state.removed_incoming) {
            uv_state.removed_incoming = false;
            v.candidate_incoming.emplace_back(uv_idx);
        }
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
            assert(aa != bf_costs_from_u.end());
            assert(aa->second == a);
            auto bb = bf_costs_from_x.find(uv.to);
            assert(bb != bf_costs_from_u.end());
            assert(bb->second == b);
#endif
            if (d <= uv.weight) {
                ++stats_.true_edges;
#ifdef CROSSCHECK
                auto edges = changed_edges_;
                edges.emplace_back(uv_idx);
                // NOTE: throws if there is a cycle
                try {
                    bellman_ford(changed_edges_, uv.from);
                }
                catch (...) {
                    assert(false && "edge is implied but lead to a conflict :(");
                }
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
        assert(v.relevant_to || u.relevant_from);

        if (v.relevant_to && u.relevant_from) {
            auto &xy = edges_[xy_idx];
            auto &x = nodes_[xy.from];
            auto &y = nodes_[xy.to];

            auto a = v.cost_to + y.potential() - v.potential();
            auto b = u.cost_from + u.potential() - x.potential();
            auto d = a + b - xy.weight;
            if (d < -uv.weight) {
                ++stats_.false_edges;
                if (!ctl.assignment().is_false(uv.lit)) {
#ifdef CROSSCHECK
                    T sum = uv.weight - xy.weight;
#endif
                    std::vector<literal_t> clause;
                    clause.push_back(-uv.lit);
                    // forward
                    for (auto next_edge_idx = u.path_from; next_edge_idx >= 0;) {
                        auto &next_edge = edges_[next_edge_idx];
                        auto &next_node = nodes_[next_edge.from];
                        clause.push_back(-next_edge.lit);
#ifdef CROSSCHECK
                        sum += next_edge.weight;
#endif
                        next_edge_idx = next_node.path_from;
                    }
                    // backward
                    for (auto next_edge_idx = v.path_to; next_edge_idx >= 0;) {
                        auto &next_edge = edges_[next_edge_idx];
                        auto &next_node = nodes_[next_edge.to];
                        clause.push_back(-next_edge.lit);
#ifdef CROSSCHECK
                        sum += next_edge.weight;
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
                try {
                    bellman_ford(changed_edges_, uv.from);
                }
                catch (...) {
                    assert(false && "edge must not cause a conflict");
                }
            }
#endif
        }
        return false;
    }

    template <class M>
    bool propagate_edges(M &m, Clingo::PropagateControl &ctl, int xy_idx, bool forward, bool backward) {
        if (!forward && !backward) {
            return true;
        }
        for (auto &node : m.visited_set()) {
            if (m.relevant(node)) {
                if (forward) {
                    auto &in = m.candidate_incoming(node);
                    in.resize(
                        std::remove_if(
                            in.begin(), in.end(),
                            [&](int uv_idx) {
                                if (!edge_states_[uv_idx].active || propagate_edge_true(uv_idx, xy_idx)) {
                                    m.remove_incoming(uv_idx);
                                    return true;
                                }
                                return false;
                            }) -
                        in.begin());
                }
                if (backward) {
                    bool ret = true;
                    auto &out = m.candidate_outgoing(node);
                    out.resize(
                        std::remove_if(
                            out.begin(), out.end(),
                            [&](int uv_idx) {
                                if (!ret) {
                                    return false;
                                }
                                if (!edge_states_[uv_idx].active || propagate_edge_false(ctl, uv_idx, xy_idx, ret)) {
                                    m.remove_outgoing(uv_idx);
                                    return true;
                                }
                                return false;
                            }) -
                        out.begin());
                    if (!ret) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    template <class M>
    std::pair<int, int> dijkstra(int source_idx, std::vector<int> &visited_set, M &m) {
        int relevant = 0;
        int relevant_degree_out = 0, relevant_degree_in = 0;
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
            if (relevant_u) {
                relevant_degree_out += nodes_[u_idx].degree_out;
                relevant_degree_in += nodes_[u_idx].degree_in;
            }
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
                        if (relevant_u) {
                            ++relevant;
                        }
                        visited_set.push_back(m.to(uv_idx));
                        m.visited(v_idx) = true;
                        costs_heap_.push(m, v_idx);
                    }
                    else {
                        if (m.relevant(m.from(m.path(v_idx)))) {
                            // node v no longer contributes a relevant edge
                            if (!relevant_u) {
                                --relevant;
                            }
                        }
                        // node v contributes a relevant edge now
                        else if (relevant_u) {
                            ++relevant;
                        }
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
        return {relevant_degree_out, relevant_degree_in};
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
    Heap<4> costs_heap_;
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
    void reset() {
        time_total = steady_clock::duration::zero();
        time_init  = steady_clock::duration::zero();
        conflicts  = 0;
        choices    = 0;
        restarts   = 0;
        for (auto& i : dl_stats) {
            i.reset();
        }
    }
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
    DifferenceLogicGraph<T> dl_graph;
};

template <typename T>
T evaluate_binary(char const *op, T left, T right) {
    if (std::strcmp(op, "+") == 0) {
        return left + right;
    }
    if (std::strcmp(op, "-") == 0) {
        return left - right;
    }
    if (std::strcmp(op, "*") == 0) {
        return left * right;
    }
    if (std::strcmp(op, "/") == 0) {
        if (std::is_integral<T>::value && right == 0) {
            throw std::runtime_error("could not evaluate term: division by zero");
        }
        return left / right;
    }
    throw std::runtime_error("could not evaluate term: unknown binary operator");
}

template <typename T>
T evaluate(Clingo::TheoryTerm term);

template <typename T>
T get_weight(TheoryAtom const &atom);

Clingo::Symbol evaluate_term(Clingo::TheoryTerm term);

class ExtendedPropagator {
public:
    virtual bool extend_model(Model &model) = 0;
    void addName(const std::string& s) {
        names_.emplace_back(s);
    }
    char const* name(size_t index) { return names_[index].c_str(); }
    virtual size_t numVars() const = 0;
    virtual bool hasLowerBound(uint32_t threadId, size_t index) = 0;
    virtual double lowerBound(uint32_t threadId, size_t index) = 0;
private:
    std::vector<std::string> names_;
};
template <typename T>
class DifferenceLogicPropagator : public Propagator, public ExtendedPropagator {
public:
    DifferenceLogicPropagator(Stats &stats, bool strict, bool propagate)
        : stats_(stats)
        , strict_(strict)
        , propagate_(propagate)
        {
            map_vert(Clingo::Number(0));
        }

public:
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
        char const *msg = "parsing difference constraint failed: only constraints of form &diff {u - v} <= b are accepted";
        int lit = init.solver_literal(atom.literal());
        if (!atom.has_guard()) {
            throw std::runtime_error(msg);
        }
        T weight = get_weight<T>(atom);
        auto elems = atom.elements();
        if (elems.size() != 1) {
            throw std::runtime_error(msg);
        }
        auto tuple = elems[0].tuple();
        if (tuple.size() != 1) {
            throw std::runtime_error(msg);
        }
        auto term = tuple[0];
        if (term.type() != Clingo::TheoryTermType::Function || std::strcmp(term.name(), "-") != 0) {
            throw std::runtime_error(msg);
        }
        auto args = term.arguments();
        if (args.size() != 2) {
            throw std::runtime_error(msg);
        }
        auto x = evaluate_term(args[0]);
        addName(x.to_string());
        auto u_id = map_vert(x);
        x = evaluate_term(args[1]);
        addName(x.to_string());
        auto v_id = map_vert(x);
        auto id = numeric_cast<int>(edges_.size());
        edges_.push_back({u_id, v_id, weight, lit});
        lit_to_edges_.emplace(lit, id);
        init.add_watch(lit);
        if (propagate_) {
            false_lit_to_edges_.emplace(-lit, id);
            init.add_watch(-lit);
        }
        if (strict_) {
            auto id = numeric_cast<int>(edges_.size());
            edges_.push_back({v_id, u_id, -weight - 1, -lit});
            lit_to_edges_.emplace(-lit, id);
            if (propagate_) {
                false_lit_to_edges_.emplace(lit, id);
            }
            else {
                init.add_watch(-lit);
            }
        }
    }

    int map_vert(Clingo::Symbol v) {
        auto ret = vert_map_inv_.emplace(v, static_cast<int>(vert_map_.size()));
        if (ret.second) {
            vert_map_.emplace_back(ret.first->first);
        }
        return ret.first->second;
    }

    void initialize_states(PropagateInit &init) {
        stats_.dl_stats.resize(init.number_of_threads());
        states_.clear();
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
        for (auto lit : changes) {
            for (auto it = false_lit_to_edges_.find(lit), ie = false_lit_to_edges_.end(); it != ie && it->first == lit; ++it) {
                if (state.dl_graph.edge_is_active(it->second)) {
                    state.dl_graph.remove_candidate_edge(it->second);
                }
            }
            for (auto it = lit_to_edges_.find(lit), ie = lit_to_edges_.end(); it != ie && it->first == lit; ++it) {
                if (state.dl_graph.edge_is_active(it->second)) {
                    if (!state.dl_graph.add_edge(it->second, [&](std::vector<int> const &neg_cycle) {
                        std::vector<literal_t> clause;
                        for (auto eid : neg_cycle) {
                            clause.emplace_back(-edges_[eid].lit);
                        }
                        return ctl.add_clause(clause) && ctl.propagate();
                    })) {
                        return;
                    }
                    else if (propagate_ && !state.dl_graph.propagate(it->second, ctl)) {
                        return;
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
                if (!state.dl_graph.node_value_defined(x.from) || !state.dl_graph.node_value_defined(x.to) || !(state.dl_graph.node_value(x.from) - state.dl_graph.node_value(x.to) <= x.weight)) {
                    throw std::logic_error("not a valid solution");
                }
            }
        }
    }
#endif
    bool extend_model(Model &model) override {
        auto &state = states_[model.thread_id()];
        T adjust = 0;
        assert(vert_map_[0] == Clingo::Number(0));
        if (!state.dl_graph.empty() && state.dl_graph.node_value_defined(0)) {
            adjust = state.dl_graph.node_value(0);
        }

        SymbolVector vec;
        for (auto idx = 1; idx < vert_map_.size(); ++idx) {
            if (state.dl_graph.node_value_defined(idx)) {
                SymbolVector params;
                params.emplace_back(vert_map_[idx]);
                params.emplace_back(String(std::to_string(adjust + state.dl_graph.node_value(idx)).c_str()));
                vec.emplace_back(Function("dl",params));
            }
        }
        model.extend(vec);
        return true;
    }

    size_t numVars() const {
        return vert_map_.size();
    }

    bool hasLowerBound(uint32_t threadId, size_t index) {
        assert(index < vert_map_.size());
        auto &state = states_[threadId];
        return index > 0 && states_[threadId].dl_graph.node_value_defined(index);

    }
    virtual double lowerBound(uint32_t threadId, size_t index) {
        assert(hasLowerBound(threadId, index));
        auto &state = states_[threadId];
        T adjust = 0;
        assert(vert_map_[0] == Clingo::Number(0));
        if (state.dl_graph.node_value_defined(0)) {
            adjust = state.dl_graph.node_value(0);
        }
        return state.dl_graph.node_value(index) + adjust;
    }


private:
    std::vector<DLState<T>> states_;
    std::unordered_multimap<literal_t, int> lit_to_edges_;
    std::unordered_multimap<literal_t, int> false_lit_to_edges_;
    std::vector<Edge<T>> edges_;
    std::vector<Clingo::Symbol> vert_map_;
    std::unordered_map<Clingo::Symbol, int> vert_map_inv_;
    Stats &stats_;
    bool strict_;
    bool propagate_;
};


#endif
