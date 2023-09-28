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

#ifndef CLINGODL_CONFIG_HH
#define CLINGODL_CONFIG_HH

#include <clingo.hh>
#include <optional>
#include <utility>
#include <vector>

namespace ClingoDL {

//! Enumeration to configure propagation strength.
enum class PropagationMode {
    Check = 0,    //!< Only check for conflicting assignments.
    Trivial = 1,  //!< Check inverse constraits.
    Weak = 2,     //!< Perform weak propagation.
    WeakPlus = 3, //!< Perform weak propagation with some extra effort.
    Zero = 4,     //!< Perform propagation through zero node.
    Strong = 5    //!< Perform full propagation.
};

//! Enumeration to configure sorting of edges before propagation.
enum class SortMode {
    No = 0,          //! Do not sort edges.
    Weight = 1,      //! Sort edges by weight in ascending order.
    WeightRev = 2,   //! Sort edges by weight in descending order.
    Potential = 3,   //! Sort edges by potential in ascending order.
    PotentialRev = 4 //! Sort edges by potential in descending order.
};

//! Enumeration to configure sorting of edges before propagation.
enum class DecisionMode {
    Disabled = 0,    //! Do not modify decision heuristic.
    MinConflict = 1, //! Try to minimize conflicts.
    MaxConflict = 2, //! Try to maximize conflicts.
};

//! Default value for PropagatorConfig::sort_mode.
static constexpr SortMode SORT_EDGES{SortMode::Weight};
//! Default value for PropagatorConfig::decision_mode.
static constexpr DecisionMode DECISON_MODE{DecisionMode::Disabled};
//! Default value for PropagatorConfig::mutex_size.
static constexpr uint64_t MUTEX_SIZE{0};
//! Default value for PropagatorConfig::mutex_cutoff.
static constexpr uint64_t MUTEX_CUTOFF{10};
//! Default value for PropagatorConfig::sort_mode.
static constexpr uint64_t PROPAGATE_ROOT{0};
//! Default value for PropagatorConfig::propagate_budget.
static constexpr uint64_t PROPAGATE_BUDGET{0};
//! Default value for PropagatorConfig::propagate_mode.
static constexpr PropagationMode PROPAGATE_MODE{PropagationMode::Check};
//! Default value for PropagatorConfig::calculate_cc.
static constexpr bool CALCULATE_CC{true};

//! Struct to configure per thread options.
struct ThreadConfig {
    //! See PropagatorConfig::propagate_root.
    std::optional<uint64_t> propagate_root;
    //! See PropagatorConfig::propagate_budget.
    std::optional<uint64_t> propagate_budget;
    //! See PropagatorConfig::propagate_mode.
    std::optional<PropagationMode> propagate_mode;
    //! See PropagatorConfig::sort_mode.
    std::optional<SortMode> sort_mode;
};

//! Struct to configure a propagator.
struct PropagatorConfig {
    //! Configure sorting of edges before propagation.
    SortMode sort_mode{SORT_EDGES};
    //! Configure sorting of edges before propagation.
    DecisionMode decision_mode{DECISON_MODE};
    //! Maximum size of mutexes to add during preprocessing.
    uint64_t mutex_size{MUTEX_SIZE};
    //! Maximum budget to calculate a mutex.
    uint64_t mutex_cutoff{MUTEX_CUTOFF};
    //! Enable full propagation below this level.
    uint64_t propagate_root{PROPAGATE_ROOT};
    //! Maximum budget to spend before disabling propagation.
    uint64_t propagate_budget{PROPAGATE_BUDGET};
    //! Configure propagation strength.
    PropagationMode propagate_mode{PROPAGATE_MODE};
    //! Per thread configuration.
    std::vector<ThreadConfig> thread_config;
    //! Enable component optimization.
    bool calculate_cc{CALCULATE_CC};

    //! Get per thread propagate_root if present or global value.
    [[nodiscard]] auto get_propagate_root(Clingo::id_t thread_id) const -> uint64_t {
        return get_prop(thread_id, propagate_root, &ThreadConfig::propagate_root);
    }
    //! Get per thread propagate_budget if present or global value.
    [[nodiscard]] auto get_propagate_budget(Clingo::id_t thread_id) const -> uint64_t {
        return get_prop(thread_id, propagate_budget, &ThreadConfig::propagate_budget);
    }
    //! Get per thread propagate_mode if present or global value.
    [[nodiscard]] auto get_propagate_mode(Clingo::id_t thread_id) const -> PropagationMode {
        return get_prop(thread_id, propagate_mode, &ThreadConfig::propagate_mode);
    }
    //! Get per thread sort_mode if present or global value.
    [[nodiscard]] auto get_sort_mode(Clingo::id_t thread_id) const -> SortMode {
        return get_prop(thread_id, sort_mode, &ThreadConfig::sort_mode);
    }
    //! Return the thread config for the given thread or return a default
    //! constructed object if absent.
    auto ensure(Clingo::id_t thread_id) -> ThreadConfig & {
        if (thread_config.size() < thread_id + 1) {
            thread_config.resize(thread_id + 1);
        }
        return thread_config[thread_id];
    }

  private:
    //! Helper to access per thread or global properties.
    template <class T, class P> [[nodiscard]] auto get_prop(Clingo::id_t thread_id, T &&def, P &&prop) const -> T {
        if (thread_id < thread_config.size() && thread_config[thread_id].*prop) {
            return *(thread_config[thread_id].*prop); // NOLINT(bugprone-unchecked-optional-access)
        }
        return def;
    }
};

} // namespace ClingoDL

#endif // CLINGODL_CONFIG_HH
