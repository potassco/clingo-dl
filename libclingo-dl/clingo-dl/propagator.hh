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

#ifndef CLINGODL_PROPAGATOR_HH
#define CLINGODL_PROPAGATOR_HH

#include <clingo.hh>
#include <clingo-dl/util.hh>

#include <unordered_map>

#if CLINGO_VERSION_MAJOR*1000 + CLINGO_VERSION_MINOR >= 5005
#   define CLINGODL_UNDO_NOEXCEPT noexcept
#else
#   define CLINGODL_UNDO_NOEXCEPT
#endif

//#define CROSSCHECK
#define CHECKSOLUTION

namespace {
inline std::string quote(Clingo::StringSpan str) {
    std::string res;
    for (auto c : str) {
        switch (*c) {
            case '\n': {
                res.push_back('\\');
                res.push_back('n');
                break;
            }
            case '\\': {
                res.push_back('\\');
                res.push_back('\\');
                break;
            }
            case '"': {
                res.push_back('\\');
                res.push_back('"');
                break;
            }
            default: {
                res.push_back(*c);
                break;
            }
        }
    }
    return res;
}

inline std::string unquote(Clingo::StringSpan str) {
    std::string res;
    bool slash = false;
    for (auto c : str) {
        if (slash) {
            switch (*c) {
                case 'n': {
                    res.push_back('\n');
                    break;
                }
                case '\\': {
                    res.push_back('\\');
                    break;
                }
                case '"': {
                    res.push_back('"');
                    break;
                }
                default: {
                    assert(false);
                    break;
                }
            }
            slash = false;
        }
        else if (*c == '\\') { slash = true; }
        else { res.push_back(*c); }
    }
    return res;
}
}
namespace ClingoDL {
template <typename T=void>
inline T throw_syntax_error(char const *message="Invalid Syntax") {
    throw std::runtime_error(message);
}

inline void check_syntax(bool condition, char const *message="Invalid Syntax") {
    if (!condition) {
        throw_syntax_error(message);
    }
}

char const *negate_relation(char const *op);

using NodeCallback = std::function<void(Clingo::AST::Node &&ast)>;

//! Transform the given statement with dl constraints and pass it on to the
//! given callback.
//!
//! Optionally shifts constraints from rule bodies into heads of integrity
//! constraints if possible.
void transform(Clingo::AST::Node const &ast, NodeCallback const &cb, bool shift);

template <typename T>
struct Edge {
    int from;
    int to;
    T weight;
    Clingo::literal_t lit;
};

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

struct PropagatorConfig {
    SortMode sort_edges{SortMode::Weight};
    uint64_t mutex_size{0};
    uint64_t mutex_cutoff{10};
    uint64_t propagate_root{0};
    uint64_t propagate_budget{0};
    PropagationMode mode{PropagationMode::Check};
    std::vector<ThreadConfig> thread_config;

    uint64_t get_propagate_root(Clingo::id_t thread_id) {
        if (thread_id < thread_config.size() && thread_config[thread_id].propagate_root.first) {
            return thread_config[thread_id].propagate_root.second;
        }
        return propagate_root;
    }
    uint64_t get_propagate_budget(Clingo::id_t thread_id) {
        if (thread_id < thread_config.size() && thread_config[thread_id].propagate_budget.first) {
            return thread_config[thread_id].propagate_budget.second;
        }
        return propagate_budget;
    }
    PropagationMode get_propagate_mode(Clingo::id_t thread_id) {
        if (thread_id < thread_config.size() && thread_config[thread_id].mode.first) {
            return thread_config[thread_id].mode.second;
        }
        return mode;
    }
    SortMode get_sort_mode(Clingo::id_t thread_id) {
        if (thread_id < thread_config.size() && thread_config[thread_id].sort_edges.first) {
            return thread_config[thread_id].sort_edges.second;
        }
        return sort_edges;
    }


    ThreadConfig &ensure(Clingo::id_t thread_id) {
        if (thread_config.size() < thread_id + 1) {
            thread_config.resize(thread_id + 1);
        }
        return thread_config[thread_id];
    }
};


template <typename T>
class DifferenceLogicGraph : private HeapToM<T, DifferenceLogicGraph<T>>, private HeapFromM<T, DifferenceLogicGraph<T>> {
    using HTM = HeapToM<T, DifferenceLogicGraph<T>>;
    using HFM = HeapFromM<T, DifferenceLogicGraph<T>>;
    friend struct HeapToM<T, DifferenceLogicGraph<T>>;
    friend struct HeapFromM<T, DifferenceLogicGraph<T>>;

public:
    DifferenceLogicGraph(DLStats &stats, const std::vector<Edge<T>> &edges, PropagationMode propagate)
        : edges_(edges)
        , propagate_(propagate)
        , stats_(stats) {
        edge_states_.resize(edges_.size(), {1, 1, 0});
        for (int i = 0; i < numeric_cast<int>(edges_.size()); ++i) {
            ensure_index(nodes_, std::max(edges_[i].from, edges_[i].to));
            add_candidate_edge(i);
        }
    }

    bool empty() const { return nodes_.empty(); }
    bool valid_node(int idx) const { return nodes_.size() > idx; }
    int node_value_defined(int idx) const { return nodes_[idx].defined(); }
    bool has_value(int idx) const { return valid_node(idx) && node_value_defined(idx); }

    T node_value(int idx) const { return -nodes_[idx].potential(); }

    bool edge_is_active(int edge_idx) const { return edge_states_[edge_idx].active; }

    bool can_propagate() const {
        return std::get<4>(changed_trail_.back());
    }

    void disable_propagate() {
        std::get<4>(changed_trail_.back()) = false;
    }

    void ensure_decision_level(int level, bool enable_propagate) {
        // initialize the trail
        if (changed_trail_.empty() || static_cast<int>(std::get<0>(changed_trail_.back())) < level) {
            bool can_propagate = (changed_trail_.empty() || std::get<4>(changed_trail_.back())) && enable_propagate;
            changed_trail_.emplace_back(level, static_cast<int>(changed_nodes_.size()),
                                               static_cast<int>(changed_edges_.size()),
                                               static_cast<int>(inactive_edges_.size()),
                                               can_propagate);
        }
    }

    std::vector<int> neg_cycle;
    template <class P, class F>
    bool with_incoming(int s_idx, P p, F f) {
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
            neg_cycle.clear();
            if (f(t_idx, ts_idx)) {
                edge_states_[ts_idx].removed_incoming = true;
                remove_candidate_edge(ts_idx);
                if (!p(neg_cycle)) {
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
    bool cheap_propagate(int u_idx, int s_idx, F f) {
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
                        neg_cycle.emplace_back(r.path_from);
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
                        neg_cycle.emplace_back(ts_idx);
                        return true;
                    }
                }
            }
            return false;
        });
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
            neg_cycle.clear();
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

        if (propagate_ >= PropagationMode::Trivial && consistent) {
            if (visited_from_.empty() || propagate_ == PropagationMode::Trivial) {
                consistent = with_incoming(uv.from, f, [&](int t_idx, int ts_idx) {
                    auto &ts = edges_[ts_idx];
                    if (t_idx == uv.to && uv.weight + ts.weight < 0) {
                        neg_cycle.emplace_back(uv_idx);
                        neg_cycle.emplace_back(ts_idx);
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

    bool propagate(int xy_idx, Clingo::PropagateControl &ctl) {
        ++stats_.edges_propagated;
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

    PropagationMode mode() const {
        return propagate_;
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
                    std::vector<Clingo::literal_t> clause;
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
                ++m.propagation_cost();
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
    std::vector<std::tuple<int, int, int, int, bool>> changed_trail_;
    std::vector<int> inactive_edges_;
    std::vector<EdgeState> edge_states_;
    DLStats &stats_;
    PropagationMode propagate_;
};

struct Stats {
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
    void accu(Stats const &x) {
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
    std::vector<DLStats> dl_stats;
};

struct FactState {
    std::vector<Clingo::literal_t> lits;
    size_t limit{0};
};

template <typename T>
struct DLState {
    DLState(DLStats &stats, const std::vector<Edge<T>> &edges, PropagationMode propagate, uint64_t propagate_root, uint64_t propagate_budget)
        : stats(stats)
        , dl_graph(stats, edges, propagate)
        , propagate_root{propagate_root}
        , propagate_budget{propagate_budget} { }
    DLStats &stats;
    DifferenceLogicGraph<T> dl_graph;
    std::vector<Clingo::literal_t> false_lits;
    std::vector<int> todo_edges;
    uint64_t propagate_root;
    uint64_t propagate_budget;
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

constexpr int INVALID_VAR{std::numeric_limits<int>::max()};
bool match(Clingo::TheoryTerm const &term, char const *name, size_t arity);

//! Test whether a variable is valid.
inline bool is_valid_var(int var) {
    return var < INVALID_VAR;
}



template <typename T>
class DifferenceLogicPropagator : public Clingo::Propagator {
private:

    using CoVarVec = std::vector<std::pair<T,int>>; // vector of coefficients and variables
public:
    DifferenceLogicPropagator(Stats &stats, PropagatorConfig const &conf)
    : stats_(stats)
    , conf_{conf} {
        zero_nodes_.emplace_back(map_vert(Clingo::Number(0)));
        cc_reset();
    }

public:
    // initialization

    void init(Clingo::PropagateInit &init) override {
        if (!edges_.empty()) {
            init.set_check_mode(Clingo::PropagatorCheckMode::Partial);
        }

        int edge_start = edges_.size();

        Timer t{stats_.time_init};
        for (auto atom : init.theory_atoms()) {
            auto term = atom.term();
            if (term.to_string() == "__diff_h" || term.to_string() == "__diff_b") {
                add_edge_atom(init, atom);
            }
        }

        std::unordered_multimap<int, int> outgoing;
        std::unordered_multimap<int, int> incoming;
        // build adjacency list
        for (int edge_id = 0, size = edges_.size(); edge_id < size; ++edge_id) {
            outgoing.emplace(std::make_pair(edges_[edge_id].from, edge_id));
            incoming.emplace(std::make_pair(edges_[edge_id].to, edge_id));
        }

        cc(outgoing, incoming);

        stats_.edges = edges_.size();
        stats_.variables = num_vertices();

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
        if (conf_.mutex_size > 0 && conf_.mutex_cutoff > 0) {
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
                            else if (cc.lit == st.lit) {
                                found = 1;
                            }
                            c = queue[c].previous;
                        }
                        if (found == 2) { continue; }
                        else if (found == 0 && st_truth == Clingo::TruthValue::Free) { n += 1; }
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

        initialize_states(init);
    }

    template <class N, typename std::enable_if<std::is_integral<N>::value, int>::type = 0>
    T round(double val) {
        if (ceilf(val) == val) {
            return static_cast<T>(val);
        }
        throw std::runtime_error("could not evaluate term: for real numbers use option rdl");
    }

    template <class N, typename std::enable_if<std::is_floating_point<N>::value, int>::type = 0>
    T round(T val) {
        return static_cast<T>(val);
    }


    bool evaluate_real(char const *name, T &val) {
        static const std::string chars = "\"";
        auto len = std::strlen(name);
        if (len < 2 || name[0] != '"' || name[len - 1] != '"') {
            return false;
        }
        char *parsed = nullptr;
        auto ret = std::strtod(name + 1, &parsed);
        if (parsed != name + len - 1) {
            return false;
        }
        val = round<T>(ret);
        return true;
    }

    void parse_constraint_elem(Clingo::TheoryTerm const &term, CoVarVec &res) {
        if (term.type() == Clingo::TheoryTermType::Number) {
            res.emplace_back(term.number(), INVALID_VAR);
        }
        else if (match(term, "+", 2)) {
            auto args = term.arguments();
            parse_constraint_elem(args.front(), res);
            parse_constraint_elem(args.back(), res);
        }
        else if (match(term, "-", 2)) {
            auto args = term.arguments();
            parse_constraint_elem(args.front(), res);
            auto pos = res.size();
            parse_constraint_elem(args.back(), res);
            for (auto it = res.begin() + pos, ie = res.end(); it != ie; ++it) {
                it->first = safe_inv(it->first);
            }
        }
        else if (match(term, "-", 1)) {
            auto pos = res.size();
            parse_constraint_elem(term.arguments().front(), res);
            for (auto it = res.begin() + pos, ie = res.end(); it != ie; ++it) {
                it->first = safe_inv(it->first);
            }
        }
        else if (match(term, "+", 1)) {
            parse_constraint_elem(term.arguments().front(), res);
        }
        else if (match(term, "*", 2)) {
            auto args = term.arguments();
            CoVarVec lhs, rhs; // NOLINT
            parse_constraint_elem(args.front(), lhs);
            parse_constraint_elem(args.back(), rhs);
            for (auto &l : lhs) {
                for (auto &r : rhs) {
                    if (!is_valid_var(l.second)) {
                        res.emplace_back(safe_mul(l.first, r.first), r.second);
                    }
                    else if (!is_valid_var(r.second)) {
                        res.emplace_back(safe_mul(l.first, r.first), l.second);
                    }
                    else {
                        throw_syntax_error("Invalid Syntax: only linear difference constraints are supported");
                    }
                }
            }
        }
        else if (term.type() == Clingo::TheoryTermType::Symbol) {
            T val;
            if (evaluate_real(term.name(), val)) {
                res.emplace_back(val, INVALID_VAR);
            }
            else {
                res.emplace_back(1, map_vert(evaluate(term)));
            }
        }
        else if (term.type() == Clingo::TheoryTermType::Function || term.type() == Clingo::TheoryTermType::Tuple) {
            res.emplace_back(1, map_vert(evaluate(term)));
        }
        else {
            throw_syntax_error("Invalid Syntax: invalid diff constraint");
        }
    }

    T simplify(CoVarVec &vec) const {
        static thread_local std::unordered_map<int, typename CoVarVec::iterator> seen;
        T rhs = 0;
        seen.clear();

        auto jt = vec.begin();
        for (auto it = jt, ie = vec.end(); it != ie; ++it) {
            auto &co = it->first;
            auto &var = it->second;
            if (co == 0) {
                continue;
            }
            if (!is_valid_var(var)) {
                rhs = safe_sub<T>(rhs, co);
            }
            else {
                auto r = seen.emplace(var, jt);
                auto kt = r.first;
                auto ins = r.second;
                if (!ins) {
                    kt->second->first = safe_add<T>(kt->second->first, co);
                }
                else {
                    if (it != jt) {
                        *jt = *it;
                    }
                    ++jt;
                }
            }
        }

        jt = std::remove_if(vec.begin(), jt, [](auto &co_var) { return co_var.first == 0; } );
        vec.erase(jt, vec.end());
        return rhs;
    }

    void add_edge_atom(Clingo::PropagateInit &init, Clingo::TheoryAtom const &atom) {
        char const *msg = "parsing difference constraint failed: only constraints of form &diff {u - v} <= b are accepted";
        int lit = init.solver_literal(atom.literal());
        if (!atom.has_guard()) {
            throw std::runtime_error(msg);
        }

        auto term = atom.term();
        bool strict = (term.to_string() == "__diff_b");
        if (strict && std::is_floating_point<T>::value) {
            std::runtime_error("real difference logic not available with strict semantics in the body of a rule");
        }
        CoVarVec covec;
        parse_constraint_elem(atom.guard().second, covec);
        for (auto it = covec.begin(), ie = covec.end(); it != ie; ++it) {
            it->first = safe_inv<T>(it->first);
        }
        auto rel = atom.guard().first;

        auto elems = atom.elements();
        if (elems.size() > 1) {
            throw std::runtime_error(msg);
        }
        for (auto const &element : elems) {
            auto tuple = element.tuple();
            check_syntax(!tuple.empty() && element.condition().empty(), "Invalid Syntax: invalid sum constraint");
            parse_constraint_elem(element.tuple().front(), covec);
        }
        auto rhs = simplify(covec);
        normalize_constraint(init, lit, covec, rel, rhs, strict);
    }

    bool add_edge(Clingo::PropagateInit& init, int literal, CoVarVec const &covec, T rhs, bool strict) {
        char const *msg = "parsing difference constraint failed: only constraints of form &diff {u - v} <= b are accepted";
        if (strict && init.assignment().is_false(literal)) {
            return true;
        }
        if (covec.size() > 2) {
            throw std::runtime_error(msg);
        }
        auto u_id = map_vert(Clingo::Number(0));
        auto v_id = map_vert(Clingo::Number(0));
        if (covec.size() == 0) {
            if (rhs < 0) {
                return init.add_clause({-literal});
            }
            return !strict || init.add_clause({literal});
        }
        else if (covec.size() == 1) {
            if (covec[0].first == 1) {
                u_id = covec[0].second;
            }
            else if (covec[0].first == -1) {
                v_id = covec[0].second;
            }
            else throw std::runtime_error(msg);
        }
        else if (covec.size() == 2) {
            if (covec[0].first == 1) {
                u_id = covec[0].second;
                if (covec[1].first == -1) {
                    v_id = covec[1].second;
                }
                else throw std::runtime_error(msg);
            }
            else if (covec[0].first == -1) {
                v_id = covec[0].second;
                if (covec[1].first == 1) {
                    u_id = covec[1].second;
                }
                else throw std::runtime_error(msg);
            }
            else throw std::runtime_error(msg);
        }
        add_edge(init, u_id, v_id, rhs, literal, strict);
        return true;
    }

    void add_edge(Clingo::PropagateInit& init, int u_id, int v_id, T weight, int lit, bool strict) {
        add_edge(init, u_id, v_id, weight, lit);
        if (strict) {
            add_edge(init, v_id, u_id, -weight-1, -lit);
        }
    }

    void add_edge(Clingo::PropagateInit &init, int u_id, int v_id, T weight, int lit) {
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

    bool normalize_constraint(Clingo::PropagateInit &init, int literal, CoVarVec const &elements, char const *op, T rhs, bool strict) {
        CoVarVec copy;
        CoVarVec const *elems = &elements;

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
            op = "<=";
            rhs = safe_inv<T>(rhs);
            copy.reserve(elements.size());
            for (auto const &covar : elements) {
                copy.emplace_back(safe_inv<T>(covar.first), covar.second);
            }
            elems = &copy;
        }

        // hanle remainig '<=', '=', and '!='
        if (std::strcmp(op, "<=") == 0) {
            if (!init.assignment().is_true(-literal) && !add_edge(init, literal, *elems, rhs, false)) {
                return false;
            }
        }
        else if (std::strcmp(op, "=") == 0) {
            int a, b;
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

            if (!normalize_constraint(init, a, *elems, "<=", rhs, strict)) {
                return false;
            }
            if (!normalize_constraint(init, b, *elems, ">=", rhs, strict)) {
                return false;
            }

            if (strict) {
                return true;
            }
        }
        else if (std::strcmp(op, "!=") == 0) {
            if (strict) {
                return normalize_constraint(init, -literal, *elems, "=", rhs, true);
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

            if (!normalize_constraint(init, a, *elems, "<", rhs, false)) {
                return false;
            }
            if (!normalize_constraint(init, b, *elems, ">", rhs, false)) {
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

            if (!normalize_constraint(init, -literal, *elems, op, rhs, false)) {
                return false;
            }
        }

        return true;
    }

    void cc_reset() {
        node_info_.clear();
        node_info_.resize(vert_map_.size());
        for (unsigned int i = 0; i < zero_nodes_.size(); ++i) {
            node_info_[zero_nodes_[i]] = NodeInfo(i, true);
        }
    }

    bool cc_visited(int node) const {
        return node_info_[node].visited;
    }

    bool is_zero(int node) const {
        assert(node_info_[node].cc < zero_nodes_.size());
        return zero_nodes_[node_info_[node].cc] == node;
    }

    void cc(std::unordered_multimap<int, int> &outgoing, std::unordered_multimap<int, int> &incoming) {
        uint32_t cc = 0;
        // Note that this marks zero nodes as visited.
        cc_reset();

        std::vector<int> node_stack;
        for (int node = 0; node < vert_map_.size(); ++node) {
            if (cc_visited(node)) {
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
                    if (!cc_visited(add_node)) {
                        node_info_[add_node] = NodeInfo(cc, true);
                        node_stack.emplace_back(add_node);
                    }
                }
                edges = incoming.equal_range(node);
                for (auto edge = edges.first; edge != edges.second; ++edge) {
                    auto add_node = edges_[edge->second].from;
                    if (!cc_visited(add_node)) {
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
            auto node = map_vert(Clingo::Function("__null", {Clingo::Number(i)}));
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

    int map_vert(Clingo::Symbol v) {
        auto ret = vert_map_inv_.emplace(v, static_cast<int>(vert_map_.size()));
        if (ret.second) {
            vert_map_.emplace_back(ret.first->first);
        }
        return ret.first->second;
    }

    void initialize_states(Clingo::PropagateInit &init) {
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

    void check(Clingo::PropagateControl &ctl) override {
        DLState<T> &state = states_[ctl.thread_id()];
        auto &facts = facts_[ctl.thread_id()];
        auto assignment = ctl.assignment();
        if (assignment.decision_level() == 0 && facts.limit > 0) {
            do_propagate(ctl, {facts.lits.data(), facts.lits.data() + facts.limit});
            facts.limit = 0;
        }
#if defined(CHECKSOLUTION) || defined(CROSSCHECK)
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

    void disable_edge_by_lit(DLState<T> &state, Clingo::literal_t lit) {
        for (auto it = false_lit_to_edges_.find(lit), ie = false_lit_to_edges_.end(); it != ie && it->first == lit; ++it) {
            if (state.dl_graph.edge_is_active(it->second)) {
                state.dl_graph.remove_candidate_edge(it->second);
            }
        }
    }

    int get_potential_(DifferenceLogicGraph<T> const &graph, int idx) {
        return graph.node_value_defined(idx) ? -graph.node_value(idx) : 0;
    };

    int cost_(DifferenceLogicGraph<T> const &graph, Edge<T> const &edge) {
        return get_potential_(graph, edge.from) + edge.weight - get_potential_(graph, edge.to);
    };

    int cost_(SortMode mode, DifferenceLogicGraph<T> const &graph, int i) {
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

    void sort_edges(SortMode mode, DLState<T> &state) {
        std::sort(state.todo_edges.begin(), state.todo_edges.end(), [&](int l, int r) {
            return cost_(mode, state.dl_graph, l) < cost_(mode, state.dl_graph, r);
        });
    }

    void do_propagate(Clingo::PropagateControl &ctl, Clingo::LiteralSpan changes) {
        auto thread_id = ctl.thread_id();
        DLState<T> &state = states_[thread_id];
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
            auto it = lit_to_edges_.find(lit), ie = lit_to_edges_.end();
            if (state.dl_graph.can_propagate()) { disable_edge_by_lit(state, lit); }
            else if (it == ie) {
                state.false_lits.emplace_back(lit);
                ctl.remove_watch(lit);
            }
            for (; it != ie && it->first == lit; ++it) {
                if (state.dl_graph.edge_is_active(it->second)) state.todo_edges.push_back(it->second);
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

    void propagate(Clingo::PropagateControl &ctl, Clingo::LiteralSpan changes) override {
        if (ctl.assignment().decision_level() == 0) {
            auto &facts = facts_[ctl.thread_id()];
            facts.lits.insert(facts.lits.end(), changes.begin(), changes.end());
        }
        do_propagate(ctl, changes);
    }

    // undo

    void undo(Clingo::PropagateControl const &ctl, Clingo::LiteralSpan changes) CLINGODL_UNDO_NOEXCEPT override {
        static_cast<void>(changes);
        auto &state = states_[ctl.thread_id()];
        Timer t{state.stats.time_undo};
        state.dl_graph.backtrack();
    }

    void extend_model(Clingo::Model &model) {
        auto &state = states_[model.thread_id()];
        std::vector<T> adjust;
        adjust.reserve(zero_nodes_.size());
        for (auto node : zero_nodes_) {
            adjust.emplace_back(state.dl_graph.has_value(node) ? state.dl_graph.node_value(node) : 0);
        }

        Clingo::SymbolVector vec;
        for (auto idx = 0; idx < vert_map_.size(); ++idx) {
            if (!is_zero(idx) && state.dl_graph.has_value(idx)) {
                Clingo::SymbolVector params;
                params.emplace_back(vert_map_[idx]);
                auto cc = node_info_[idx].cc;
                params.emplace_back(to_symbol<T>(state.dl_graph.node_value(idx) - adjust[cc]));
                vec.emplace_back(Function("dl", params));
            }
        }
        model.extend(vec);
    }

    size_t num_vertices() const {
        return vert_map_.size();
    }

    Clingo::Symbol symbol(size_t index) const {
        return vert_map_[index];
    }

    uint32_t lookup(clingo_symbol_t symbol) {
        auto it = vert_map_inv_.find(Clingo::Symbol(symbol));
        return it != vert_map_inv_.end()
            ? it->second
            : num_vertices();
    }

    bool has_lower_bound(uint32_t thread_id, size_t index) const {
        return index < vert_map_.size() && !is_zero(index) && states_[thread_id].dl_graph.has_value(index);
    }

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
private:

    T to_T(Clingo::Symbol const &a) const {
        if (a.type() == Clingo::SymbolType::Number) {
            return static_cast<T>(a.number());
        }
        if (a.type() == Clingo::SymbolType::String) {
            return std::stod(a.string());
        }
        return throw_syntax_error<T>();
    }

    template <class F, class N, typename std::enable_if<std::is_integral<N>::value, int>::type = 0>
    Clingo::Symbol evaluate(Clingo::TheoryTerm const &a, Clingo::TheoryTerm const &b, F f) const {
        auto ea = evaluate(a);
        check_syntax(ea.type() == Clingo::SymbolType::Number);
        auto eb = evaluate(b);
        check_syntax(eb.type() == Clingo::SymbolType::Number);
        return Clingo::Number(f(to_T(ea), to_T(eb)));
    }

    template <class F, class N, typename std::enable_if<std::is_floating_point<N>::value, int>::type = 0>
    Clingo::Symbol evaluate(Clingo::TheoryTerm const &a, Clingo::TheoryTerm const &b, F f) const {
        auto ea = evaluate(a);
        auto eb = evaluate(b);
        return Clingo::String(std::to_string(f(to_T(ea), to_T(eb))).c_str());
    }

    template <class F>
    Clingo::Symbol evaluate(Clingo::TheoryTerm const &a, Clingo::TheoryTerm const &b, F &&f) const {
        return evaluate<F, T>(a, b, std::forward<F>(f));
    }

    template <class N, typename std::enable_if<std::is_integral<N>::value, int>::type = 0>
    T epsilon() const {
        return 1;
    }

    template <class N, typename std::enable_if<std::is_floating_point<N>::value, int>::type = 0>
    T epsilon() const {
        return 0.00001;
    }

    Clingo::Symbol evaluate(Clingo::TheoryTerm const &term) const {
        if (term.type() == Clingo::TheoryTermType::Symbol) {
            const char* const name = term.name();
            if (name[0]=='\"' && name[strlen(name)-1]=='\"') {
                return Clingo::String(quote({name,strlen(name)}).c_str());
            }
            else
                return Clingo::Function(unquote({name,strlen(name)}).c_str(), {});
        }

        if (term.type() == Clingo::TheoryTermType::Number) {
            return Clingo::Number(term.number());
        }

        if (match(term, "+", 2)) {
            return evaluate(term.arguments().front(), term.arguments().back(), safe_add<T>);
        }
        if (match(term, "-", 2)) {
            return evaluate(term.arguments().front(), term.arguments().back(), safe_sub<T>);
        }
        if (match(term, "*", 2)) {
            return evaluate(term.arguments().front(), term.arguments().back(), safe_mul<T>);
        }
        if (match(term, "/", 2)) {
            return evaluate(term.arguments().front(), term.arguments().back(), safe_div<T>);
        }
        if (match(term, "\\", 2)) {
            return evaluate(term.arguments().front(), term.arguments().back(), safe_mod<T>);
        }
        if (match(term, "**", 2)) {
            return evaluate(term.arguments().front(), term.arguments().back(), safe_pow<T>);
        }

        if (match(term, "-", 1)) {
            auto ea = evaluate(term.arguments().front());
            if (ea.type() == Clingo::SymbolType::Number) {
                return Clingo::Number(safe_inv(ea.number()));
            }
            if (ea.type() == Clingo::SymbolType::Function && std::strlen(ea.name()) > 0) {
                return Clingo::Function(ea.name(), ea.arguments(), !ea.is_positive());
            }
            return throw_syntax_error<Clingo::Symbol>();
        }

        check_syntax(!match(term, "..", 2));

        if (term.type() == Clingo::TheoryTermType::Tuple || term.type() == Clingo::TheoryTermType::Function) {
            std::vector<Clingo::Symbol> args;
            args.reserve(term.arguments().size());
            for (auto const &arg : term.arguments()) {
                args.emplace_back(evaluate(arg));
            }
            return Clingo::Function(term.type() == Clingo::TheoryTermType::Function ? term.name() : "", args);
        }
        return throw_syntax_error<Clingo::Symbol>();
    }

private:

    std::vector<DLState<T>> states_;
    std::vector<FactState> facts_;
    std::unordered_multimap<Clingo::literal_t, int> lit_to_edges_;
    std::unordered_multimap<Clingo::literal_t, int> false_lit_to_edges_;
    std::vector<Edge<T>> edges_;
    std::vector<Clingo::Symbol> vert_map_;
    std::unordered_map<Clingo::Symbol, int> vert_map_inv_;
    std::vector<NodeInfo> node_info_;
    std::vector<int> zero_nodes_;
    Stats &stats_;
    PropagatorConfig conf_;
};

} // namespace ClingoDL

#endif // CLINGODL_PROPAGATOR_HH
