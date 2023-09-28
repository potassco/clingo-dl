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

#include <clingo-dl.h>
#include <clingo-dl/propagator.hh>

#include <clingo.h>
#include <clingo.hh>
#include <sstream>

#define CLINGODL_TRY try // NOLINT
#define CLINGODL_CATCH                                                                                                 \
    catch (...) {                                                                                                      \
        Clingo::Detail::handle_cxx_error();                                                                            \
        return false;                                                                                                  \
    }                                                                                                                  \
    return true // NOLINT

using namespace ClingoDL;

namespace {

using Clingo::Detail::handle_error;

//! C initialization callback for the DL propagator.
template <typename T> auto init(clingo_propagate_init_t *i, void *data) -> bool {
    CLINGODL_TRY {
        Clingo::PropagateInit in(i);
        static_cast<DLPropagator<T> *>(data)->init(in);
    }
    CLINGODL_CATCH;
}

//! C propagation callback for the DL propagator.
template <typename T>
auto propagate(clingo_propagate_control_t *i, const clingo_literal_t *changes, size_t size, void *data) -> bool {
    CLINGODL_TRY {
        Clingo::PropagateControl in(i);
        static_cast<DLPropagator<T> *>(data)->propagate(in, {changes, size});
    }
    CLINGODL_CATCH;
}

//! C undo callback for the DL propagator.
template <typename T>
void undo(clingo_propagate_control_t const *i, const clingo_literal_t *changes, size_t size, void *data) {
    Clingo::PropagateControl in(const_cast<clingo_propagate_control_t *>(i)); // NOLINT
    static_cast<DLPropagator<T> *>(data)->undo(in, {changes, size});
}

//! C check callback for the DL propagator.
template <typename T> auto check(clingo_propagate_control_t *i, void *data) -> bool {
    CLINGODL_TRY {
        Clingo::PropagateControl in(i);
        static_cast<DLPropagator<T> *>(data)->check(in);
    }
    CLINGODL_CATCH;
}

//! C decide callback for the DL heuristic.
template <typename T>
auto decide(clingo_id_t thread_id, clingo_assignment_t const *assignment, clingo_literal_t fallback, void *data,
            clingo_literal_t *decision) -> bool {
    CLINGODL_TRY {
        Clingo::Assignment ass(assignment);
        *decision = static_cast<DLPropagator<T> *>(data)->decide(thread_id, ass, fallback);
    }
    CLINGODL_CATCH;
}

//! High level interface to use the DL propagator hiding the value type.
class PropagatorFacade {
  public:
    PropagatorFacade() = default;
    PropagatorFacade(PropagatorFacade const &other) = default;
    PropagatorFacade(PropagatorFacade &&other) = default;
    auto operator=(PropagatorFacade const &other) -> PropagatorFacade & = default;
    auto operator=(PropagatorFacade &&other) noexcept -> PropagatorFacade & = default;
    virtual ~PropagatorFacade() = default;

    //! Look up the index of a symbol.
    //!
    //! The function returns false if the symbol could not be found.
    virtual auto lookup_symbol(clingo_symbol_t name, size_t *index) -> bool = 0;
    //! Get the symbol associated with an index.
    virtual auto get_symbol(size_t index) -> clingo_symbol_t = 0;
    //! Check if a symbol has a value in a thread.
    virtual auto has_value(uint32_t thread_id, size_t index) -> bool = 0;
    //! Get the value of a symbol in a thread.
    virtual void get_value(uint32_t thread_id, size_t index, clingodl_value_t *value) = 0;
    //! Function to iterato over the thread specific assignment of symbols and values.
    //!
    //! Argument current should initially be set to 0. The function returns
    //! false if no more values are available.
    virtual auto next(uint32_t thread_id, size_t *current) -> bool = 0;
    //! Extend the given model with the assignment stored in the propagator.
    virtual void extend_model(Clingo::Model &m) = 0;
    //! Add the propagator statistics to clingo's statistics.
    virtual void on_statistics(Clingo::UserStatistics &step, Clingo::UserStatistics &accu) = 0;
};

//! Set variant to an integer value.
void set_value(clingodl_value_t *variant, int value) {
    variant->type = clingodl_value_type_int;
    variant->int_number = value; // NOLINT
}

//! Set variant to a double value.
void set_value(clingodl_value_t *variant, double value) {
    variant->type = clingodl_value_type_double;
    variant->double_number = value; // NOLINT
}

//! High level interface to use the DL propagator.
template <typename T> class DLPropagatorFacade : public PropagatorFacade {
  public:
    DLPropagatorFacade(clingo_control_t *control, PropagatorConfig const &conf) : prop_{step_, conf} {
        handle_error(clingo_control_add(control, "base", nullptr, 0, THEORY));
        static clingo_propagator_t prop = {init<T>, propagate<T>, undo<T>, check<T>,
                                           conf.decision_mode != DecisionMode::Disabled ? decide<T> : nullptr};
        handle_error(clingo_control_register_propagator(control, &prop, &prop_, false));
    }

    auto lookup_symbol(clingo_symbol_t name, size_t *index) -> bool override {
        *index = prop_.lookup(Clingo::Symbol{name}) + 1;
        return *index <= prop_.num_vertices();
    }

    auto get_symbol(size_t index) -> clingo_symbol_t override {
        return prop_.symbol(numeric_cast<vertex_t>(index - 1)).to_c();
    }

    auto has_value(uint32_t thread_id, size_t index) -> bool override {
        return prop_.has_lower_bound(thread_id, numeric_cast<vertex_t>(index - 1));
    }

    void get_value(uint32_t thread_id, size_t index, clingodl_value_t *value) override {
        assert(index > 0 && index <= prop_.num_vertices());
        set_value(value, prop_.lower_bound(thread_id, numeric_cast<vertex_t>(index - 1)));
    }

    auto next(uint32_t thread_id, size_t *current) -> bool override {
        for (++*current; *current <= prop_.num_vertices(); ++*current) {
            if (prop_.has_lower_bound(thread_id, numeric_cast<vertex_t>(*current - 1))) {
                return true;
            }
        }
        return false;
    }

    void extend_model(Clingo::Model &m) override { prop_.extend_model(m); }

    void on_statistics(Clingo::UserStatistics &step, Clingo::UserStatistics &accu) override {
        accu_.accu(step_);
        add_statistics_(step, step_);
        add_statistics_(accu, accu_);
        step_.reset();
    }

  private:
    //! Add an integral value to the statistics.
    template <class V, std::enable_if_t<std::is_integral_v<V>, bool> = true>
    static void add_subkey_(Clingo::UserStatistics &root, char const *name, V value) {
        root.add_subkey(name, Clingo::StatisticsType::Value).set_value(static_cast<double>(value));
    }
    //! Add an floating point value to the statistics.
    template <class V, std::enable_if_t<std::is_floating_point_v<V>, bool> = true>
    static void add_subkey_(Clingo::UserStatistics &root, char const *name, V value) {
        root.add_subkey(name, Clingo::StatisticsType::Value).set_value(value);
    }

    //!< Helper function to add the DL statistics to clingo's statistics.
    void add_statistics_(Clingo::UserStatistics &root, Statistics const &stats) {
        Clingo::UserStatistics diff = root.add_subkey("DifferenceLogic", Clingo::StatisticsType::Map);
        add_subkey_(diff, "Time init(s)", stats.time_init.count());
        add_subkey_(diff, "CCs", stats.ccs);
        add_subkey_(diff, "Mutexes", stats.mutexes);
        add_subkey_(diff, "Edges", stats.edges);
        add_subkey_(diff, "Variables", stats.variables);
        Clingo::UserStatistics threads = diff.add_subkey("Thread", Clingo::StatisticsType::Array);
        threads.ensure_size(stats.thread_statistics.size(), Clingo::StatisticsType::Map);
        auto it = threads.begin();
        for (auto const &stat : stats.thread_statistics) {
            auto thread = *it++;
            add_subkey_(thread, "Propagation(s)", stat.time_propagate.count());
            add_subkey_(thread, "Dijkstra(s)", stat.time_dijkstra.count());
            add_subkey_(thread, "Undo(s)", stat.time_undo.count());
            add_subkey_(thread, "True edges", stat.true_edges);
            add_subkey_(thread, "False edges", stat.false_edges);
            add_subkey_(thread, "False edges (inverse)", stat.false_edges_trivial);
            add_subkey_(thread, "False edges (partial)", stat.false_edges_weak);
            add_subkey_(thread, "False edges (partial+)", stat.false_edges_weak_plus);
            add_subkey_(thread, "Edges added", stat.edges_added);
            add_subkey_(thread, "Edges skipped", stat.edges_skipped);
            add_subkey_(thread, "Edges propagated", stat.edges_propagated);
            add_subkey_(thread, "Cost consistency", stat.propagate_cost_add);
            add_subkey_(thread, "Cost forward", stat.propagate_cost_from);
            add_subkey_(thread, "Cost backward", stat.propagate_cost_to);
        }
    }

    Statistics step_;      //!< Per step statistics.
    Statistics accu_;      //!< Accumulated statistics over all steps.
    DLPropagator<T> prop_; //!< The underlying difference logic propagator.
};

//! Check if b is a lower case prefix of a returning a pointer to the remainder of a.
auto iequals_pre(char const *a, char const *b) -> char const * {
    for (; *a && *b; ++a, ++b) { // NOLINT
        if (tolower(*a) != tolower(*b)) {
            return nullptr;
        }
    }
    return *b != '\0' ? nullptr : a;
}

//! Check if two strings are lower case equal.
auto iequals(char const *a, char const *b) -> bool {
    a = iequals_pre(a, b);
    return a != nullptr && *a == '\0';
}

//! Turn the largest prefix of value into an unsigned integer and return the remainder.
//!
//! The function returns a nullpointer if there are no leading digits. The
//! result is stored in data which is assumed to be a pointer to an uint64_t.
auto parse_uint64_pre(const char *value, void *data) -> char const * {
    auto &res = *static_cast<uint64_t *>(data);
    char const *it = value;
    res = 0;

    for (; *it != '\0'; ++it) { // NOLINT
        if ('0' <= *it && *it <= '9') {
            auto tmp = res;
            res *= 10; // NOLINT
            res += *it - '0';
            if (res < tmp) {
                return nullptr;
            }
        } else {
            break;
        }
    }

    return value != it ? it : nullptr;
}

//! Turn the value into an uint64_t assuming that data is a pointer to an
//! uint64_t.
auto parse_uint64(const char *value, void *data) -> bool {
    value = parse_uint64_pre(value, data);
    return value != nullptr && *value == '\0';
}

//! Parse thread-specific option via a callback.
//!
//! The thread number is optional and can follow separated with a comma.
template <typename F, typename G> auto set_config(char const *value, void *data, F f, G g) -> bool {
    try {
        auto &config = *static_cast<PropagatorConfig *>(data);
        uint64_t id = 0;
        if (*value == '\0') {
            f(config);
            return true;
        }
        if (*value == ',' && parse_uint64(value + 1, &id) && id < 64) { // NOLINT
            g(config.ensure(id));
            return true;
        }
    } catch (...) {
    }
    return false;
}

//! Parse a level to limit full propagation.
//!
//! Return false if there is a parse error.
auto parse_root(const char *value, void *data) -> bool {
    uint64_t x = 0;
    return (value = parse_uint64_pre(value, &x)) != nullptr &&
           set_config(
               value, data, [x](PropagatorConfig &config) { config.propagate_root = x; },
               [x](ThreadConfig &config) { config.propagate_root = x; });
}

//! Parse the propagation budget and store it data.
//!
//! Return false if there is a parse error.
auto parse_budget(const char *value, void *data) -> bool {
    uint64_t x = 0;
    return (value = parse_uint64_pre(value, &x)) != nullptr &&
           set_config(
               value, data, [x](PropagatorConfig &config) { config.propagate_budget = x; },
               [x](ThreadConfig &config) { config.propagate_budget = x; });
}

//! Parse the mutex detection mode and store it data.
//!
//! Return false if there is a parse error.
auto parse_mutex(const char *value, void *data) -> bool {
    auto &pc = *static_cast<PropagatorConfig *>(data);
    uint64_t x = 0;
    if ((value = parse_uint64_pre(value, &x)) == nullptr) {
        return false;
    }
    pc.mutex_size = x;
    if (*value == '\0') {
        pc.mutex_cutoff = 10 * x; // NOLINT
        return true;
    }
    if (*value == ',') {
        if (!parse_uint64(value + 1, &x)) {
            return false;
        } // NOLINT
        pc.mutex_cutoff = x;
    }
    return true;
}

//! Parse the propagation mode and store it data.
//!
//! Return false if there is a parse error.
auto parse_mode(const char *value, void *data) -> bool {
    PropagationMode mode = PropagationMode::Check;
    char const *rem = nullptr;
    if (rem = iequals_pre(value, "no"); rem != nullptr) {
        mode = PropagationMode::Check;
    } else if (rem = iequals_pre(value, "inverse"); rem != nullptr) {
        mode = PropagationMode::Trivial;
    } else if (rem = iequals_pre(value, "partial+"); rem != nullptr) {
        mode = PropagationMode::WeakPlus;
    } else if (rem = iequals_pre(value, "partial"); rem != nullptr) {
        mode = PropagationMode::Weak;
    } else if (rem = iequals_pre(value, "zero"); rem != nullptr) {
        mode = PropagationMode::Zero;
    } else if (rem = iequals_pre(value, "full"); rem != nullptr) {
        mode = PropagationMode::Strong;
    }
    return rem != nullptr && set_config(
                                 rem, data, [mode](PropagatorConfig &config) { config.propagate_mode = mode; },
                                 [mode](ThreadConfig &config) { config.propagate_mode = mode; });
}

//! Parse the sort mode and store it data.
//!
//! Return false if there is a parse error.
auto parse_sort(const char *value, void *data) -> bool {
    SortMode sort = SortMode::Weight;
    char const *rem = nullptr;
    if (rem = iequals_pre(value, "no"); rem != nullptr) {
        sort = SortMode::No;
    } else if (rem = iequals_pre(value, "weight-reversed"); rem != nullptr) {
        sort = SortMode::WeightRev;
    } else if (rem = iequals_pre(value, "weight"); rem != nullptr) {
        sort = SortMode::Weight;
    } else if (rem = iequals_pre(value, "potential-reversed"); rem != nullptr) {
        sort = SortMode::PotentialRev;
    } else if (rem = iequals_pre(value, "potential"); rem != nullptr) {
        sort = SortMode::Potential;
    }
    return rem != nullptr && set_config(
                                 rem, data, [sort](PropagatorConfig &config) { config.sort_mode = sort; },
                                 [sort](ThreadConfig &config) { config.sort_mode = sort; });
}

//! Parse the decision mode.
//!
//! Return false if there is a parse error.
auto parse_decide(const char *value, void *data) -> bool {
    DecisionMode mode = DecisionMode::Disabled;
    if (iequals(value, "no")) {
        mode = DecisionMode::Disabled;
    } else if (iequals(value, "min")) {
        mode = DecisionMode::MinConflict;
    } else if (iequals(value, "max")) {
        mode = DecisionMode::MaxConflict;
    }
    static_cast<PropagatorConfig *>(data)->decision_mode = mode;
    return true;
}

//! Parse a Boolean and store it in data.
//!
//! Return false if there is a parse error.
auto parse_bool(const char *value, void *data) -> bool {
    auto &result = *static_cast<bool *>(data);
    if (iequals(value, "no") || iequals(value, "off") || iequals(value, "0")) {
        result = false;
        return true;
    }
    if (iequals(value, "yes") || iequals(value, "on") || iequals(value, "1")) {
        result = true;
        return true;
    }
    return false;
}

//! Set the given error message if the Boolean is false.
//!
//! Return false if there is a parse error.
auto check_parse(char const *key, bool ret) -> bool {
    if (!ret) {
        std::ostringstream msg;
        msg << "invalid value for '" << key << "'";
        clingo_set_error(clingo_error_runtime, msg.str().c_str());
    }
    return ret;
}

} // namespace

struct clingodl_theory {
    std::unique_ptr<PropagatorFacade> clingodl{nullptr};
    PropagatorConfig config;
    bool rdl;
    bool shift_constraints{false};
};

extern "C" void clingodl_version(int *major, int *minor, int *patch) {
    if (major != nullptr) {
        *major = CLINGODL_VERSION_MAJOR;
    }
    if (minor != nullptr) {
        *minor = CLINGODL_VERSION_MINOR;
    }
    if (patch != nullptr) {
        *patch = CLINGODL_VERSION_REVISION;
    }
}

extern "C" auto clingodl_create(clingodl_theory_t **theory) -> bool {
    CLINGODL_TRY { *theory = new clingodl_theory{}; } // NOLINT
    CLINGODL_CATCH;
}

extern "C" auto clingodl_register(clingodl_theory_t *theory, clingo_control_t *control) -> bool {
    CLINGODL_TRY {
        if (!theory->rdl) {
            theory->clingodl = std::make_unique<DLPropagatorFacade<int>>(control, theory->config);
        } else {
            theory->clingodl = std::make_unique<DLPropagatorFacade<double>>(control, theory->config);
        }
    }
    CLINGODL_CATCH;
}

extern "C" auto clingodl_rewrite_ast(clingodl_theory_t *theory, clingo_ast_t *ast, clingodl_ast_callback_t add,
                                     void *data) -> bool {
    CLINGODL_TRY {
        clingo_ast_acquire(ast);
        Clingo::AST::Node ast_cpp{ast};
        transform(
            ast_cpp, [add, data](Clingo::AST::Node &&ast_trans) { handle_error(add(ast_trans.to_c(), data)); },
            theory->shift_constraints);
    }
    CLINGODL_CATCH;
}

extern "C" auto clingodl_prepare(clingodl_theory_t *theory, clingo_control_t *control) -> bool {
    static_cast<void>(theory);
    static_cast<void>(control);
    return true;
}

extern "C" auto clingodl_destroy(clingodl_theory_t *theory) -> bool {
    CLINGODL_TRY { delete theory; } // NOLINT
    CLINGODL_CATCH;
}

extern "C" auto clingodl_configure(clingodl_theory_t *theory, char const *key, char const *value) -> bool {
    CLINGODL_TRY {
        if (strcmp(key, "propagate") == 0) {
            return check_parse("propagate", parse_mode(value, &theory->config));
        }
        if (strcmp(key, "propagate-root") == 0) {
            return check_parse("propagate-root", parse_root(value, &theory->config));
        }
        if (strcmp(key, "propagate-budget") == 0) {
            return check_parse("propgate-budget", parse_budget(value, &theory->config));
        }
        if (strcmp(key, "add-mutexes") == 0) {
            return check_parse("add-mutexes", parse_mutex(value, &theory->config));
        }
        if (strcmp(key, "sort-edges") == 0) {
            return check_parse("sort-edges", parse_sort(value, &theory->config));
        }
        if (strcmp(key, "rdl") == 0) {
            return check_parse("rdl", parse_bool(value, &theory->rdl));
        }
        if (strcmp(key, "dl-heuristic") == 0) {
            return check_parse("dl-heuristic", parse_decide(value, &theory->config));
        }
        if (strcmp(key, "shift-constraints") == 0) {
            return check_parse("shift-constraints", parse_bool(value, &theory->shift_constraints));
        }
        if (strcmp(key, "compute-components") == 0) {
            return check_parse("compute-components", parse_bool(value, &theory->config.calculate_cc));
        }
        std::ostringstream msg;
        msg << "invalid configuration key '" << key << "'";
        clingo_set_error(clingo_error_runtime, msg.str().c_str());
        return false;
    }
    CLINGODL_CATCH;
}

extern "C" auto clingodl_register_options(clingodl_theory_t *theory, clingo_options_t *options) -> bool {
    CLINGODL_TRY {
        char const *group = "Clingo.DL Options";
        handle_error(clingo_options_add(options, group, "propagate",
                                        "Set propagation mode [no]\n"
                                        "      <mode>  : {no,inverse,partial,partial+,zero,full}[,<thread>]\n"
                                        "        no      : No propagation; only detect conflicts\n"
                                        "        inverse : Check inverse constraints\n"
                                        "        partial : Detect some conflicts\n"
                                        "        partial+: Detect some more conflicts\n"
                                        "        zero    : Detect all immediate conflicts through zero nodes\n"
                                        "        full    : Detect all immediate conflicts\n"
                                        "      <thread>: Restrict to thread",
                                        &parse_mode, &theory->config, true, "<mode>"));
        handle_error(clingo_options_add(options, group, "propagate-root",
                                        "Enable full propagation below decision level [0]\n"
                                        "      <arg>   : <n>[,<thread>]\n"
                                        "      <n>     : Upper bound for decision level\n"
                                        "      <thread>: Restrict to thread",
                                        &parse_root, &theory->config, true, "<arg>"));
        handle_error(clingo_options_add(options, group, "propagate-budget",
                                        "Enable full propagation limiting to budget [0]\n"
                                        "      <arg>   : <n>[,<thread>]\n"
                                        "      <n>     : Budget roughly corresponding to cost of consistency checks\n"
                                        "                (if possible use with --propagate-root greater 0)\n"
                                        "      <thread>: Restrict to thread",
                                        &parse_budget, &theory->config, true, "<arg>"));
        handle_error(clingo_options_add(options, group, "add-mutexes",
                                        "Add mutexes in a preprocessing step [0]\n"
                                        "      <arg>: <max>[,<cut>]\n"
                                        "      <max>: Maximum size of mutexes to add\n"
                                        "      <cut>: Limit costs to calculate mutexes",
                                        &parse_mutex, &theory->config, true, "<arg>"));
        handle_error(clingo_options_add(options, group, "sort-edges",
                                        "Sort edges for propagation [weight]\n"
                                        "      <arg>: {no, weight, weight-reversed, potential, potential-reversed}\n"
                                        "        no                : No sorting\n"
                                        "        weight            : Sort by edge weight\n"
                                        "        weight-reversed   : Sort by negative edge weight\n"
                                        "        potential         : Sort by relative potential\n"
                                        "        potential-reversed: Sort by relative negative potential",
                                        &parse_sort, &theory->config, true, "<arg>"));
        handle_error(clingo_options_add(options, group, "dl-heuristic",
                                        "Decision heuristic for difference constraints\n"
                                        "      <arg>: {none, min, max}\n"
                                        "        no : Use default decision heuristic\n"
                                        "        min: Try to minimize conflicts\n"
                                        "        max: Try to maximize conflicts",
                                        &parse_decide, &theory->config, false, "<arg>"));
        handle_error(
            clingo_options_add_flag(options, group, "rdl", "Enable support for real numbers [no]", &theory->rdl));
        handle_error(clingo_options_add_flag(options, group, "shift-constraints",
                                             "Shift constraints into head of integrity constraints [no]",
                                             &theory->shift_constraints));
        handle_error(clingo_options_add_flag(options, group, "compute-components", "Compute connected components [yes]",
                                             &theory->config.calculate_cc));
    }
    CLINGODL_CATCH;
}

extern "C" auto clingodl_validate_options(clingodl_theory_t *theory) -> bool {
    static_cast<void>(theory);
    return true;
}

extern "C" auto clingodl_on_model(clingodl_theory_t *theory, clingo_model_t *model) -> bool {
    CLINGODL_TRY {
        Clingo::Model m(model);
        theory->clingodl->extend_model(m);
    }
    CLINGODL_CATCH;
}

extern "C" auto clingodl_lookup_symbol(clingodl_theory_t *theory, clingo_symbol_t symbol, size_t *index) -> bool {
    return theory->clingodl->lookup_symbol(symbol, index);
}

extern "C" auto clingodl_get_symbol(clingodl_theory_t *theory, size_t index) -> clingo_symbol_t {
    return theory->clingodl->get_symbol(index);
}

extern "C" void clingodl_assignment_begin(clingodl_theory_t *theory, uint32_t thread_id, size_t *index) {
    // Note: the first vertex is always 0 and can be skipped because its value is 0
    static_cast<void>(theory);
    static_cast<void>(thread_id);
    *index = 1;
}

extern "C" auto clingodl_assignment_next(clingodl_theory_t *theory, uint32_t thread_id, size_t *index) -> bool {
    return theory->clingodl->next(thread_id, index);
}

extern "C" auto clingodl_assignment_has_value(clingodl_theory_t *theory, uint32_t thread_id, size_t index) -> bool {
    return theory->clingodl->has_value(thread_id, index);
}

extern "C" void clingodl_assignment_get_value(clingodl_theory_t *theory, uint32_t thread_id, size_t index,
                                              clingodl_value_t *value) {
    theory->clingodl->get_value(thread_id, index, value);
}

extern "C" auto clingodl_on_statistics(clingodl_theory_t *theory, clingo_statistics_t *step, clingo_statistics_t *accu)
    -> bool {
    CLINGODL_TRY {
        uint64_t root_s{0};
        uint64_t root_a{0};
        handle_error(clingo_statistics_root(step, &root_s));
        handle_error(clingo_statistics_root(accu, &root_a));
        Clingo::UserStatistics s(step, root_s);
        Clingo::UserStatistics a(accu, root_a);
        theory->clingodl->on_statistics(s, a);
    }
    CLINGODL_CATCH;
}

#undef CLINGODL_TRY
#undef CLINGODL_CATCH
