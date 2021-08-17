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

#ifndef CLINGODL_CONFIG_HH
#define CLINGODL_CONFIG_HH

#include <clingo.hh>
#include <utility>
#include <vector>

namespace ClingoDL {

enum class PropagationMode {
    Check = 0,
    Trivial = 1,
    Weak = 2,
    WeakPlus = 3,
    Strong = 4
};

enum class SortMode {
    No = 0,
    Weight = 1,
    WeightRev = 2,
    Potential = 3,
    PotentialRev = 4
};

static constexpr SortMode SORT_EDGES{SortMode::Weight};
static constexpr uint64_t MUTEX_SIZE{0};
static constexpr uint64_t MUTEX_CUTOFF{10};
static constexpr uint64_t PROPAGATE_ROOT{0};
static constexpr uint64_t PROPAGATE_BUDGET{0};
static constexpr PropagationMode PROPAGATE_MODE{PropagationMode::Check};

struct ThreadConfig {
    std::optional<uint64_t> propagate_root;
    std::optional<uint64_t> propagate_budget;
    std::optional<PropagationMode> mode;
    std::optional<SortMode> sort_edges;
};

struct PropagatorConfig {
    SortMode sort_edges{SORT_EDGES};
    uint64_t mutex_size{MUTEX_SIZE};
    uint64_t mutex_cutoff{MUTEX_CUTOFF};
    uint64_t propagate_root{PROPAGATE_ROOT};
    uint64_t propagate_budget{PROPAGATE_BUDGET};
    PropagationMode mode{PROPAGATE_MODE};
    std::vector<ThreadConfig> thread_config;

    uint64_t get_propagate_root(Clingo::id_t thread_id) {
        if (thread_id < thread_config.size()) {
            return thread_config[thread_id].propagate_root.value_or(propagate_root);
        }
        return propagate_root;
    }
    uint64_t get_propagate_budget(Clingo::id_t thread_id) {
        if (thread_id < thread_config.size()) {
            return thread_config[thread_id].propagate_budget.value_or(propagate_budget);
        }
        return propagate_budget;
    }
    PropagationMode get_propagate_mode(Clingo::id_t thread_id) {
        if (thread_id < thread_config.size()) {
            return thread_config[thread_id].mode.value_or(mode);
        }
        return mode;
    }
    SortMode get_sort_mode(Clingo::id_t thread_id) {
        if (thread_id < thread_config.size()) {
            return thread_config[thread_id].sort_edges.value_or(sort_edges);
        }
        return sort_edges;
    }

    ThreadConfig &ensure(Clingo::id_t thread_id) {
        if (thread_config.size() < thread_id + 1) {
            thread_config.resize(thread_id + 1);
        }
        return thread_config[thread_id];
    }

private:
    /*
    template <typename T>
    std::optional<T> get(Clingo::id_t thread_id) {
        if (thread_id < thread_config.size()) {
            return thread_config[thread_id];
        }
        return std::nullopt;
    }
    */
};

} // namespace ClingoDL

#endif // CLINGODL_CONFIG_HH
