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

#ifndef CLINGODL_GRAPH_HH
#define CLINGODL_GRAPH_HH

#include <clingo.hh>
#include <clingo-dl/theory.hh>
#include <clingo-dl/util.hh>

namespace ClingoDL {

template <class T, class P>
struct HeapFromM {
    int &offset(int idx) { return static_cast<P *>(this)->nodes_[idx].offset; }
    T &cost(int idx) { return static_cast<P *>(this)->nodes_[idx].cost_from; }
    int to(int idx) { return static_cast<P *>(this)->edges_[idx].to; }
    int from(int idx) { return static_cast<P *>(this)->edges_[idx].from; }
    std::vector<int> &out(int idx) { return static_cast<P *>(this)->nodes_[idx].outgoing; }
    int &path(int idx) { return static_cast<P *>(this)->nodes_[idx].path_from; }
    int &visited(int idx) { return static_cast<P *>(this)->nodes_[idx].visited_from; }
    bool &relevant(int idx) { return static_cast<P *>(this)->nodes_[idx].relevant_from; }
    std::vector<int> &visited_set() { return static_cast<P *>(this)->visited_from_; }
    std::vector<int> &candidate_outgoing(int idx) { return static_cast<P *>(this)->nodes_[idx].candidate_outgoing; }
    std::vector<int> &candidate_incoming(int idx) { return static_cast<P *>(this)->nodes_[idx].candidate_incoming; }
    void remove_incoming(int idx) { static_cast<P *>(this)->edge_states_[idx].removed_incoming = true; }
    void remove_outgoing(int idx) { static_cast<P *>(this)->edge_states_[idx].removed_outgoing = true; }
    uint64_t &propagation_cost() {return static_cast<P *>(this)->stats_.propagate_cost_from; }
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
    uint64_t &propagation_cost() {return static_cast<P *>(this)->stats_.propagate_cost_to; }
};

template <typename T>
struct Edge {
    vertex_t from;
    vertex_t to;
    T weight;
    Clingo::literal_t lit;
};

template <typename T>
struct DifferenceLogicNode {
    [[nodiscard]] bool defined() const { return !potential_stack.empty(); }
    [[nodiscard]] T potential() const { return potential_stack.back().second; }

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
    int visited_from = 0;
    bool relevant_from = false;
    bool relevant_to = false;
    bool visited_to = false;
};

struct DLStats {
    void reset() {
        time_propagate = std::chrono::steady_clock::duration::zero();
        time_undo      = std::chrono::steady_clock::duration::zero();
        time_dijkstra  = std::chrono::steady_clock::duration::zero();
        true_edges            = 0;
        false_edges           = 0;
        false_edges_trivial   = 0;
        false_edges_weak      = 0;
        false_edges_weak_plus = 0;
        propagate_cost_add  = 0;
        propagate_cost_from = 0;
        propagate_cost_to   = 0;
        edges_added      = 0;
        edges_skipped    = 0;
        edges_propagated = 0;
    }
    void accu(DLStats const &x) {
        time_propagate+= x.time_propagate;
        time_undo     += x.time_undo;
        time_dijkstra += x.time_dijkstra;
        true_edges    += x.true_edges;
        false_edges   += x.false_edges;
        false_edges_trivial  += x.false_edges_trivial;
        false_edges_weak     += x.false_edges_weak;
        false_edges_weak_plus+= x.false_edges_weak_plus;
        propagate_cost_add += x.propagate_cost_add;
        propagate_cost_from+= x.propagate_cost_from;
        propagate_cost_to  += x.propagate_cost_to;
        edges_added      += x.edges_added;
        edges_skipped    += x.edges_skipped;
        edges_propagated += x.edges_propagated;
    }
    Duration time_propagate = Duration{0};
    Duration time_undo = Duration{0};
    Duration time_dijkstra = Duration{0};
    uint64_t true_edges{0};
    uint64_t false_edges{0};
    uint64_t false_edges_trivial{0};
    uint64_t false_edges_weak{0};
    uint64_t false_edges_weak_plus{0};
    uint64_t propagate_cost_add{0};
    uint64_t propagate_cost_from{0};
    uint64_t propagate_cost_to{0};
    uint64_t edges_added{0};
    uint64_t edges_skipped{0};
    uint64_t edges_propagated{0};
};

struct EdgeState {
    uint8_t removed_outgoing : 1;
    uint8_t removed_incoming : 1;
    uint8_t active : 1;
};

enum class PropagationMode { Check = 0, Trivial = 1, Weak = 2, WeakPlus = 3, Strong = 4 };
enum class SortMode { No = 0, Weight = 1, WeightRev = 2, Potential = 3 , PotentialRev = 4};

struct ThreadConfig {
    std::pair<bool,uint64_t> propagate_root{false,0};
    std::pair<bool,uint64_t> propagate_budget{false,0};
    std::pair<bool,PropagationMode> mode{false,PropagationMode::Check};
    std::pair<bool,SortMode> sort_edges{false,SortMode::Weight};
};

template <typename T>
class DifferenceLogicGraph : private HeapToM<T, DifferenceLogicGraph<T>>, private HeapFromM<T, DifferenceLogicGraph<T>> {
    using HTM = HeapToM<T, DifferenceLogicGraph<T>>;
    using HFM = HeapFromM<T, DifferenceLogicGraph<T>>;
    friend struct HeapToM<T, DifferenceLogicGraph<T>>;
    friend struct HeapFromM<T, DifferenceLogicGraph<T>>;

public:
    DifferenceLogicGraph(DLStats &stats, const std::vector<Edge<T>> &edges, PropagationMode propagate);
    [[nodiscard]] bool empty() const;
    [[nodiscard]] bool valid_node(int idx) const;
    [[nodiscard]] int node_value_defined(int idx) const;
    [[nodiscard]] bool has_value(int idx) const;

    [[nodiscard]] T node_value(int idx) const;

    [[nodiscard]] bool edge_is_active(int edge_idx) const;

    [[nodiscard]] bool can_propagate() const;
    void disable_propagate();
    void ensure_decision_level(int level, bool enable_propagate);

    template <class F>
    [[nodiscard]] bool add_edge(int uv_idx, F f) {
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
        ++stats_.edges_added;
        if (v.cost_from < 0) {
            costs_heap_.push(m, uv.to);
            visited_from_.emplace_back(uv.to);
            v.visited_from = 1;
            v.path_from = uv_idx;
        }
        else {
            ++stats_.edges_skipped;
        }

        int dfs = 0;
        // detect negative cycles
        while (!costs_heap_.empty() && !u.visited_from) {
            auto s_idx = costs_heap_.pop(m);
            auto &s = nodes_[s_idx];
            assert(s.visited_from);
            s.visited_from = ++dfs;
            set_potential(s, level, s.potential() + s.cost_from);
            for (auto st_idx : s.outgoing) {
                ++stats_.propagate_cost_add;
                assert(st_idx < numeric_cast<int>(edges_.size()));
                auto &st = edges_[st_idx];
                auto &t = nodes_[st.to];
                auto c = s.potential() + st.weight - t.potential();
                if (c < (t.visited_from ? t.cost_from : 0)) {
                    assert(c < 0);
                    t.path_from = st_idx;
                    t.cost_from = c;
                    if (!t.visited_from) {
                        t.visited_from = 1;
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
            neg_cycle_.clear();
            neg_cycle_.push_back(v.path_from);
            auto next_idx = edges_[v.path_from].from;
            while (uv.to != next_idx) {
                auto &next = nodes_[next_idx];
                neg_cycle_.push_back(next.path_from);
                next_idx = edges_[next.path_from].from;
            }
#ifdef CROSSCHECK
            T weight = 0;
            for (auto &edge_idx : neg_cycle_) {
                weight += edges_[edge_idx].weight;
            }
            assert(weight < 0);
#endif
            consistent = f(neg_cycle_);
        }

        if (propagate_ >= PropagationMode::Trivial && consistent) {
            if (visited_from_.empty() || propagate_ == PropagationMode::Trivial) {
                consistent = with_incoming(uv.from, f, [&](int t_idx, int ts_idx) {
                    auto &ts = edges_[ts_idx];
                    if (t_idx == uv.to && uv.weight + ts.weight < 0) {
                        neg_cycle_.emplace_back(uv_idx);
                        neg_cycle_.emplace_back(ts_idx);
                        ++stats_.false_edges_trivial;
                        return true;
                    }
                    return false;
                });
            }
            else if (propagate_ >= PropagationMode::Weak) {
                consistent = cheap_propagate(uv.from, uv.from, f);
                if (propagate_ >= PropagationMode::WeakPlus && consistent) {
                    for (auto &s_idx : visited_from_) {
                        if (!cheap_propagate(uv.from, s_idx, f)) {
                            consistent = false;
                            break;
                        }
                    }
                }
            }
        }
        // reset visited flags
        for (auto &x : visited_from_) {
            nodes_[x].visited_from = 0;
        }
        visited_from_.clear();
        costs_heap_.clear();

        return consistent;
    }

    [[nodiscard]] bool propagate(int xy_idx, Clingo::PropagateControl &ctl);
    void backtrack();
    void remove_candidate_edge(int uv_idx);
    [[nodiscard]] PropagationMode mode() const;

private:
    template <class P, class F>
    [[nodiscard]] bool with_incoming(int s_idx, P p, F f) {
        auto &s = nodes_[s_idx];
        auto &in = s.candidate_incoming;
        auto jt = in.begin();
        for (auto it = jt, ie = in.end(); it != ie; ++it) {
            auto &ts_idx = *it;
            auto &ts = edges_[ts_idx];
            auto t_idx = ts.from;
            auto &t = nodes_[t_idx];
            if (!edge_states_[ts_idx].active) {
                edge_states_[ts_idx].removed_incoming = true;
                continue;
            }
            neg_cycle_.clear();
            if (f(t_idx, ts_idx)) {
                edge_states_[ts_idx].removed_incoming = true;
                remove_candidate_edge(ts_idx);
                if (!p(neg_cycle_)) {
                    in.erase(jt, it+1);
                    return false;
                }
                continue;
            }
            *jt++ = *it;
        }
        in.erase(jt, in.end());
        return true;
    }

    template <class F>
    [[nodiscard]] bool cheap_propagate(int u_idx, int s_idx, F f) {
        // we check for the following case:
        // u ->* s -> * t
        //       ^-----/
        //          ts
        return with_incoming(s_idx, f, [&](int t_idx, int ts_idx) {
            auto &s = nodes_[s_idx];
            auto &t = nodes_[t_idx];
            auto &ts = edges_[ts_idx];
            if (s.visited_from < t.visited_from) {
                T weight = t.potential() - s.potential();
                if (weight + ts.weight < 0) {
                    T check = 0;
                    auto r_idx = t_idx;
                    while (u_idx != r_idx && s_idx != r_idx) {
                        auto &r = nodes_[r_idx];
                        auto &rr = edges_[r.path_from];
                        neg_cycle_.emplace_back(r.path_from);
                        r_idx = rr.from;
                        check += rr.weight;
                    }
                    if (r_idx == s_idx) {
                        if (u_idx == s_idx) {
                            ++stats_.false_edges_weak;
                        }
                        else {
                            ++stats_.false_edges_weak_plus;
                        }
                        assert(weight == check);
                        neg_cycle_.emplace_back(ts_idx);
                        return true;
                    }
                }
            }
            return false;
        });
    }

    void add_candidate_edge(int uv_idx);
    [[nodiscard]] bool propagate_edge_true(int uv_idx, int xy_idx);
    [[nodiscard]] bool propagate_edge_false(Clingo::PropagateControl &ctl, int uv_idx, int xy_idx, bool &ret);
    template <class M>
    [[nodiscard]] bool propagate_edges(M &m, Clingo::PropagateControl &ctl, int xy_idx, bool forward, bool backward);
    template <class M>
    [[nodiscard]] std::pair<int, int> dijkstra(int source_idx, std::vector<int> &visited_set, M &m);
#ifdef CROSSCHECK
    [[nodiscard]] std::unordered_map<int, T> bellman_ford(std::vector<int> const &edges, int source);
#endif
    void set_potential(DifferenceLogicNode<T> &node, int level, T potential);
    [[nodiscard]] int current_decision_level_();

    Heap<4> costs_heap_;
    std::vector<int> visited_from_;
    std::vector<int> visited_to_;
    std::vector<Edge<T>> const &edges_;
    std::vector<DifferenceLogicNode<T>> nodes_;
    std::vector<int> changed_nodes_;
    std::vector<int> changed_edges_;
    std::vector<std::tuple<int, int, int, int, bool>> changed_trail_;
    std::vector<int> inactive_edges_;
    std::vector<EdgeState> edge_states_;
    std::vector<int> neg_cycle_;
    DLStats &stats_;
    PropagationMode propagate_;
};

} // namespace ClingoDL

#endif // CLINGODL_PROPAGATOR_HH
