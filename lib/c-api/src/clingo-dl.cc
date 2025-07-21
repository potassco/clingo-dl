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
#include <clingo/propagate.hh>

#include <sstream>

using namespace ClingoDL;

namespace {

using Clingo::Detail::handle_error;

//! C initialization callback for the DL propagator.
template <typename T> auto init(clingo_assignment_t const *ass, clingo_propagate_init_t *init, void *data) -> bool {
    CLINGO_TRY { static_cast<DLPropagator<T> *>(data)->init(Clingo::Assignment{ass}, Clingo::PropagateInit{init}); }
    CLINGO_CATCH;
}

//! C propagation callback for the DL propagator.
template <typename T>
auto propagate(clingo_assignment_t const *ass, clingo_propagate_control_t *ctl, const clingo_literal_t *changes,
               size_t size, void *data) -> bool {
    CLINGO_TRY {

        static_cast<DLPropagator<T> *>(data)->propagate(Clingo::Assignment{ass}, Clingo::PropagateControl{ctl},
                                                        {changes, size});
    }
    CLINGO_CATCH;
}

//! C undo callback for the DL propagator.
template <typename T>
void undo(clingo_assignment_t const *ass, clingo_literal_t const *changes, size_t size, void *data) {
    try {
        static_cast<DLPropagator<T> *>(data)->undo(Clingo::Assignment{ass}, {changes, size});
    } catch (std::exception const &e) {
        printf("panic: %s\n", e.what());
    }
}

//! C check callback for the DL propagator.
template <typename T> auto check(clingo_assignment_t const *ass, clingo_propagate_control_t *ctl, void *data) -> bool {
    CLINGO_TRY { static_cast<DLPropagator<T> *>(data)->check(Clingo::Assignment{ass}, Clingo::PropagateControl{ctl}); }
    CLINGO_CATCH;
}

//! C decide callback for the DL heuristic.
template <typename T>
auto decide(clingo_assignment_t const *assignment, clingo_literal_t fallback, void *data, clingo_literal_t *decision)
    -> bool {
    CLINGO_TRY {
        Clingo::Assignment ass(assignment);
        *decision = static_cast<DLPropagator<T> *>(data)->decide(ass, fallback);
    }
    CLINGO_CATCH;
}

//! High level interface to use the DL propagator hiding the value type.
class PropagatorFacade {
  public:
    PropagatorFacade() = default;
    PropagatorFacade(PropagatorFacade &&other) = default;
    virtual ~PropagatorFacade() = default;

    //! Look up the index of a symbol.
    //!
    //! The function returns false if the symbol could not be found.
    virtual auto lookup_symbol(Clingo::Symbol const &name, size_t &index) -> bool = 0;
    //! Get the symbol associated with an index.
    virtual auto get_symbol(size_t index) -> Clingo::Symbol const & = 0;
    //! Check if a symbol has a value in a thread.
    virtual auto has_value(uint32_t thread_id, size_t index) -> bool = 0;
    //! Get the value of a symbol in a thread.
    virtual void get_value(uint32_t thread_id, size_t index, clingo_theory_value_t &value) = 0;
    //! Function to iterato over the thread specific assignment of symbols and values.
    //!
    //! Argument current should initially be set to 0. The function returns
    //! false if no more values are available.
    virtual auto next(uint32_t thread_id, size_t &current) -> bool = 0;
    //! Extend the given model with the assignment stored in the propagator.
    virtual void extend_model(Clingo::Model m) = 0;
    //! Add the propagator statistics to clingo's statistics.
    virtual void on_statistics(Clingo::StatsMap step, Clingo::StatsMap accu) = 0;
};

//! Set variant to an integer value.
void set_value(clingo_theory_value_t &variant, int value) {
    variant.type = clingo_theory_value_type_int;
    variant.int_number = value; // NOLINT
}

//! Set variant to a double value.
void set_value(clingo_theory_value_t &variant, double value) {
    variant.type = clingo_theory_value_type_double;
    variant.double_number = value; // NOLINT
}

//! High level interface to use the DL propagator.
template <typename T> class DLPropagatorFacade : public PropagatorFacade {
  public:
    DLPropagatorFacade(clingo_lib_t *lib, clingo_control_t *control, PropagatorConfig const &conf)
        : prop_{Clingo::Library{lib, true}, step_, conf} {
        handle_error(clingo_control_parse_string(control, THEORY, std::strlen(THEORY)));
        static clingo_propagator_t prop = {
            init<T>, nullptr,  propagate<T>,
            undo<T>, check<T>, conf.decision_mode != DecisionMode::Disabled ? decide<T> : nullptr,
            nullptr};
        handle_error(clingo_control_register_propagator(control, &prop, &prop_));
    }

    auto lookup_symbol(Clingo::Symbol const &name, size_t &index) -> bool override {
        index = prop_.lookup(name) + 1;
        return index <= prop_.num_vertices();
    }

    auto get_symbol(size_t index) -> Clingo::Symbol const & override {
        return prop_.symbol(numeric_cast<vertex_t>(index - 1));
    }

    auto has_value(uint32_t thread_id, size_t index) -> bool override {
        return prop_.has_lower_bound(thread_id, numeric_cast<vertex_t>(index - 1));
    }

    void get_value(uint32_t thread_id, size_t index, clingo_theory_value_t &value) override {
        assert(index > 0 && index <= prop_.num_vertices());
        set_value(value, prop_.lower_bound(thread_id, numeric_cast<vertex_t>(index - 1)));
    }

    auto next(uint32_t thread_id, size_t &current) -> bool override {
        for (++current; current <= prop_.num_vertices(); ++current) {
            if (prop_.has_lower_bound(thread_id, numeric_cast<vertex_t>(current - 1))) {
                return true;
            }
        }
        return false;
    }

    void extend_model(Clingo::Model m) override { prop_.extend_model(m); }

    void on_statistics(Clingo::StatsMap step, Clingo::StatsMap accu) override {
        accu_.accu(step_);
        add_statistics_(step, step_);
        add_statistics_(accu, accu_);
        step_.reset();
    }

  private:
    //! Add an integral value to the statistics.
    template <class V, std::enable_if_t<std::is_integral_v<V>, bool> = true>
    static void add_subkey_(Clingo::StatsMap root, char const *name, V value) {
        root.insert(name, Clingo::StatsType::value).value(static_cast<double>(value));
    }
    //! Add an floating point value to the statistics.
    template <class V, std::enable_if_t<std::is_floating_point_v<V>, bool> = true>
    static void add_subkey_(Clingo::StatsMap root, char const *name, V value) {
        root.insert(name, Clingo::StatsType::value).value(value);
    }

    //!< Helper function to add the DL statistics to clingo's statistics.
    void add_statistics_(Clingo::StatsMap root, Statistics const &stats) {
        auto diff = root.insert("DifferenceLogic", Clingo::StatsType::map).map();
        add_subkey_(diff, "Time init(s)", stats.time_init.count());
        add_subkey_(diff, "CCs", stats.ccs);
        add_subkey_(diff, "Mutexes", stats.mutexes);
        add_subkey_(diff, "Edges", stats.edges);
        add_subkey_(diff, "Variables", stats.variables);
        auto threads = diff.insert("Thread", Clingo::StatsType::array).array();
        std::ignore = threads.ensure(stats.thread_statistics.size() - 1, Clingo::StatsType::map);
        auto it = threads.begin();
        for (auto const &stat : stats.thread_statistics) {
            auto thread = (*it++).map();
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

//! Ascii tolower conversion.
constexpr auto tolower(char c) -> char { return (c >= 'A' && c <= 'Z') ? c + ('a' - 'A') : c; }

//! Check if b is a lower case prefix of a returning a string_view to the remainder of a.
auto iequals_pre(std::string_view a, std::string_view b) -> std::optional<std::string_view> {
    if (a.size() < b.size()) {
        return std::nullopt;
    }
    auto cmp = [](char ac, char bc) { return tolower(ac) == tolower(bc); };
    auto res = std::ranges::mismatch(b, a, cmp);
    return res.in1 == b.end() ? std::optional{a.substr(b.size())} : std::nullopt;
}

//! Check if two strings are lower case equal.
auto iequals(std::string_view a, std::string_view b) -> bool {
    auto res = iequals_pre(a, b);
    return res && res->empty();
}

//! Turn the largest prefix of value into an unsigned integer and return the value and remainder.
//!
//! The function returns a nullopt if there are no leading digits.
auto parse_uint64_pre(std::string_view value) -> std::optional<std::pair<uint64_t, std::string_view>> {
    uint64_t res = 0;
    auto const *first = value.data();
    auto const *last = value.data() + value.size();
    auto [ptr, err] = std::from_chars(first, last, res);
    return err == std::errc{} ? std::make_optional(std::make_pair(res, std::string_view(ptr, last))) : std::nullopt;
}

//! Turn the value into an uint64_t and return it as optional.
auto parse_uint64(std::string_view value) -> std::optional<uint64_t> {
    auto opt = parse_uint64_pre(value);
    return (opt && opt->second.empty()) ? std::optional{opt->first} : std::nullopt;
}

//! Parse thread-specific option via a callback.
//!
//! The thread number is optional and can follow separated with a comma.
template <typename F, typename G> auto set_config(std::string_view value, void *data, F f, G g) -> bool {
    auto &config = *static_cast<PropagatorConfig *>(data);
    if (value.empty()) {
        f(config);
        return true;
    }
    if (auto opt = value.starts_with(',') ? parse_uint64(value.substr(1)) : std::nullopt; opt && *opt < 64) {
        g(config.ensure(*opt));
        return true;
    }
    return false;
}

//! Parse a level to limit full propagation.
auto parse_root(char const *value, size_t size, void *data, bool *result) -> bool {
    CLINGO_TRY {
        auto res = parse_uint64_pre({value, size});
        *result = res && set_config(
                             res->second, data, [&](PropagatorConfig &config) { config.propagate_root = res->first; },
                             [&](ThreadConfig &config) { config.propagate_root = res->first; });
    }
    CLINGO_CATCH;
}

//! Parse the propagation budget and store it in data.
auto parse_budget(const char *value, size_t size, void *data, bool *result) -> bool {
    CLINGO_TRY {
        auto res = parse_uint64_pre({value, size});
        *result = res && set_config(
                             res->second, data, [&](PropagatorConfig &config) { config.propagate_budget = res->first; },
                             [&](ThreadConfig &config) { config.propagate_budget = res->first; });
    }
    CLINGO_CATCH;
}

//! Parse the mutex detection mode and store it in data.
auto parse_mutex(const char *value, size_t size, void *data, bool *result) -> bool {
    CLINGO_TRY {
        auto &pc = *static_cast<PropagatorConfig *>(data);
        *result = false;
        if (auto fst = parse_uint64_pre({value, size})) {
            pc.mutex_size = fst->first;
            if (fst->second.empty()) {
                pc.mutex_cutoff = 10 * fst->first;
                *result = true;
            } else if (auto snd = fst->second.starts_with(',') ? parse_uint64(fst->second.substr(1)) : std::nullopt) {
                pc.mutex_cutoff = *snd;
                *result = true;
            }
        }
    }
    CLINGO_CATCH;
}

//! Parse the propagation mode and store it in data.
auto parse_mode(char const *value, size_t size, void *data, bool *result) -> bool {
    CLINGO_TRY {
        auto mode = PropagationMode::Check;
        auto str = std::string_view{value, size};
        auto res = std::optional<std::string_view>{};
        if (res = iequals_pre(str, "no"), res) {
            mode = PropagationMode::Check;
        } else if (res = iequals_pre(str, "inverse"); res) {
            mode = PropagationMode::Trivial;
        } else if (res = iequals_pre(str, "partial+"); res) {
            mode = PropagationMode::WeakPlus;
        } else if (res = iequals_pre(str, "partial"); res) {
            mode = PropagationMode::Weak;
        } else if (res = iequals_pre(str, "zero"); res) {
            mode = PropagationMode::Zero;
        } else if (res = iequals_pre(str, "full"); res) {
            mode = PropagationMode::Strong;
        }
        *result = res && set_config(
                             *res, data, [mode](PropagatorConfig &config) { config.propagate_mode = mode; },
                             [mode](ThreadConfig &config) { config.propagate_mode = mode; });
    }
    CLINGO_CATCH;
}

//! Parse the sort mode and store it data.
//!
//! Return false if there is a parse error.
auto parse_sort(char const *value, size_t size, void *data, bool *result) -> bool {
    CLINGO_TRY {
        auto sort = SortMode::Weight;
        auto str = std::string_view{value, size};
        auto res = std::optional<std::string_view>{};
        if (res = iequals_pre(str, "no"); res) {
            sort = SortMode::No;
        } else if (res = iequals_pre(str, "weight-reversed"); res) {
            sort = SortMode::WeightRev;
        } else if (res = iequals_pre(str, "weight"); res) {
            sort = SortMode::Weight;
        } else if (res = iequals_pre(str, "potential-reversed"); res) {
            sort = SortMode::PotentialRev;
        } else if (res = iequals_pre(str, "potential"); res) {
            sort = SortMode::Potential;
        }
        *result = res && set_config(
                             *res, data, [sort](PropagatorConfig &config) { config.sort_mode = sort; },
                             [sort](ThreadConfig &config) { config.sort_mode = sort; });
    }
    CLINGO_CATCH;
}

//! Parse the decision mode.
auto parse_decide(char const *value, size_t size, void *data, bool *result) -> bool {
    CLINGO_TRY {
        auto &mode = static_cast<PropagatorConfig *>(data)->decision_mode;
        auto str = std::string_view{value, size};
        auto res = false;
        if (res = iequals(str, "no"); res) {
            mode = DecisionMode::Disabled;
        } else if (res = iequals(str, "min"); res) {
            mode = DecisionMode::MinConflict;
        } else if (res = iequals(str, "max"); res) {
            mode = DecisionMode::MaxConflict;
        }
        *result = res;
    }
    CLINGO_CATCH;
}

//! Parse a Boolean and store it in data.
auto parse_bool(const char *value, size_t size, void *data, bool *result) -> bool {
    CLINGO_TRY {
        auto &flag = *static_cast<bool *>(data);
        auto str = std::string_view{value, size};
        auto res = false;
        if (res = iequals(str, "no") || iequals(str, "off") || iequals(str, "0"); res) {
            flag = false;
        } else if (res = iequals(str, "yes") || iequals(str, "on") || iequals(str, "1"); res) {
            flag = true;
        }
        *result = res;
    }
    CLINGO_CATCH;
}

//! Set the given error message if the Boolean is false.
//!
//! Return false if there is a parse error.
template <class F, class... As> auto check_parse(char const *key, F fun, As &&...as) -> bool {
    CLINGO_TRY {
        bool res = false;
        if (!fun(std::forward<As>(as)..., &res)) {
            return false;
        }
        if (!res) {
            auto msg = std::ostringstream{};
            msg << "invalid value for '" << key << "'";
            clingo_set_error(clingo_result_invalid, msg.view().data(), msg.view().size());
        }
    }
    CLINGO_CATCH;
}

struct clingodl_theory {
    clingodl_theory(clingo_lib_t *lib) : lib{lib, true} {}
    Clingo::Library lib;
    std::unique_ptr<PropagatorFacade> clingodl;
    PropagatorConfig config;
    bool rdl{false};
    bool shift_constraints{false};

    static auto info([[maybe_unused]] void *self, clingo_string_t *name, int *major, int *minor, int *patch) -> bool {
        using namespace std::string_view_literals;
        CLINGO_TRY {
            if (name != nullptr) {
                constexpr auto str = "clingo-dl"sv;
                name->data = str.data();
                name->size = str.size();
            }
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
        CLINGO_CATCH;
    }

    static auto register_(void *self, clingo_control_t *control) -> bool {
        auto theory = static_cast<clingodl_theory *>(self);
        CLINGO_TRY {
            if (!theory->rdl) {
                theory->clingodl =
                    std::make_unique<DLPropagatorFacade<int>>(c_cast(theory->lib), control, theory->config);
            } else {
                theory->clingodl =
                    std::make_unique<DLPropagatorFacade<double>>(c_cast(theory->lib), control, theory->config);
            }
        }
        CLINGO_CATCH;
    }

    static auto rewrite_ast(void *self, clingo_ast_t *ast, clingo_theory_ast_callback_t add, void *data) -> bool {
        auto theory = static_cast<clingodl_theory *>(self);
        CLINGO_TRY {
            rewrite(
                theory->lib, Clingo::AST::Node{ast, true},
                [add, data](Clingo::AST::Node ast) { handle_error(add(c_cast(ast), data)); },
                theory->shift_constraints);
        }
        CLINGO_CATCH;
    }

    static auto prepare([[maybe_unused]] void *self, [[maybe_unused]] clingo_control_t *control) -> bool {
        return true;
    }

    static void destroy(void *self) {
        auto theory = static_cast<clingodl_theory *>(self);
        std::unique_ptr<clingodl_theory>{theory};
    }

    static auto configure(void *self, char const *key, size_t key_size, char const *value, size_t value_size) -> bool {
        CLINGO_TRY {
            auto theory = static_cast<clingodl_theory *>(self);
            auto sv_key = std::string_view{key, key_size};
            if (sv_key == "propagate") {
                return check_parse("propagate", parse_mode, value, value_size, &theory->config);
            }
            if (sv_key == "propagate-root") {
                return check_parse("propagate-root", parse_root, value, value_size, &theory->config);
            }
            if (sv_key == "propagate-budget") {
                return check_parse("propgate-budget", parse_budget, value, value_size, &theory->config);
            }
            if (sv_key == "add-mutexes") {
                return check_parse("add-mutexes", parse_mutex, value, value_size, &theory->config);
            }
            if (sv_key == "sort-edges") {
                return check_parse("sort-edges", parse_sort, value, value_size, &theory->config);
            }
            if (sv_key == "rdl") {
                return check_parse("rdl", parse_bool, value, value_size, &theory->rdl);
            }
            if (sv_key == "dl-heuristic") {
                return check_parse("dl-heuristic", parse_decide, value, value_size, &theory->config);
            }
            if (sv_key == "shift-constraints") {
                return check_parse("shift-constraints", parse_bool, value, value_size, &theory->shift_constraints);
            }
            if (sv_key == "compute-components") {
                return check_parse("compute-components", parse_bool, value, value_size, &theory->config.calculate_cc);
            }
            std::ostringstream msg;
            msg << "invalid configuration key '" << key << "'";
            clingo_set_error(clingo_result_invalid, msg.view().data(), msg.view().size());
            return false;
        }
        CLINGO_CATCH;
    }

    static auto register_options(void *self, clingo_options_t *options) -> bool {
        auto theory = static_cast<clingodl_theory *>(self);
        CLINGO_TRY {
            using namespace std::string_view_literals;
            auto group = "Clingo.DL Options"sv;
            auto opt = [&](std::string_view name, std::string_view desc, clingo_option_parser_t parser,
                           bool multi = false, std::string_view arg = {}) {
                handle_error(clingo_options_add(options, group.data(), group.size(), name.data(), name.size(),
                                                desc.data(), desc.size(), parser, &theory->config, multi,
                                                arg.empty() ? nullptr : arg.data(), arg.size()));
            };
            auto flag = [&](std::string_view name, std::string_view desc, bool &target) {
                handle_error(clingo_options_add_flag(options, group.data(), group.size(), name.data(), name.size(),
                                                     desc.data(), desc.size(), &target));
            };
            opt("propagate",
                "Set propagation mode [no]\n"
                "      <mode>  : {no,inverse,partial,partial+,zero,full}[,<thread>]\n"
                "        no      : No propagation; only detect conflicts\n"
                "        inverse : Check inverse constraints\n"
                "        partial : Detect some conflicts\n"
                "        partial+: Detect some more conflicts\n"
                "        zero    : Detect all immediate conflicts through zero nodes\n"
                "        full    : Detect all immediate conflicts\n"
                "      <thread>: Restrict to thread",
                &parse_mode, true, "<mode>");
            opt("propagate-root",
                "Enable full propagation below decision level [0]\n"
                "      <arg>   : <n>[,<thread>]\n"
                "      <n>     : Upper bound for decision level\n"
                "      <thread>: Restrict to thread",
                &parse_root, true, "<arg>");
            opt("propagate-budget",
                "Enable full propagation limiting to budget [0]\n"
                "      <arg>   : <n>[,<thread>]\n"
                "      <n>     : Budget roughly corresponding to cost of consistency checks\n"
                "                (if possible use with --propagate-root greater 0)\n"
                "      <thread>: Restrict to thread",
                &parse_budget, true, "<arg>");
            opt("add-mutexes",
                "Add mutexes in a preprocessing step [0]\n"
                "      <arg>: <max>[,<cut>]\n"
                "      <max>: Maximum size of mutexes to add\n"
                "      <cut>: Limit costs to calculate mutexes",
                &parse_mutex, true, "<arg>");
            opt("sort-edges",
                "Sort edges for propagation [weight]\n"
                "      <arg>: {no, weight, weight-reversed, potential, potential-reversed}\n"
                "        no                : No sorting\n"
                "        weight            : Sort by edge weight\n"
                "        weight-reversed   : Sort by negative edge weight\n"
                "        potential         : Sort by relative potential\n"
                "        potential-reversed: Sort by relative negative potential",
                &parse_sort, true, "<arg>");
            opt("dl-heuristic",
                "Decision heuristic for difference constraints\n"
                "      <arg>: {none, min, max}\n"
                "        no : Use default decision heuristic\n"
                "        min: Try to minimize conflicts\n"
                "        max: Try to maximize conflicts",
                &parse_decide, false, "<arg>");
            flag("rdl", "Enable support for real numbers [no]", theory->rdl);
            flag("shift-constraints", "Shift constraints into head of integrity constraints [no]",
                 theory->shift_constraints);
            flag("compute-components", "Compute connected components [yes]", theory->config.calculate_cc);
        }
        CLINGO_CATCH;
    }

    static auto validate_options([[maybe_unused]] void *self) -> bool { return true; }

    static auto on_model(void *self, clingo_model_t *model) -> bool {
        CLINGO_TRY {
            auto theory = static_cast<clingodl_theory *>(self);
            Clingo::Model m(model);
            theory->clingodl->extend_model(m);
        }
        CLINGO_CATCH;
    }

    static auto lookup_symbol(void *self, clingo_symbol_t symbol, size_t *index, bool *found) -> bool {
        CLINGO_TRY {
            auto theory = static_cast<clingodl_theory *>(self);
            *found = theory->clingodl->lookup_symbol(Clingo::Symbol{symbol, true}, *index);
        }
        CLINGO_CATCH;
    }

    static auto assignment_next(void *self, uint32_t thread_id, bool *init, size_t *index, bool *has_value) -> bool {
        CLINGO_TRY {
            auto theory = static_cast<clingodl_theory *>(self);
            if (std::exchange(*init, false)) {
                *index = 0;
            }
            *has_value = theory->clingodl->next(thread_id, *index);
        }
        CLINGO_CATCH;
    }

    static auto assignment_get_value(void *self, uint32_t thread_id, size_t index, clingo_symbol_t *symbol,
                                     clingo_theory_value_t *value, bool *has_value) -> bool {
        CLINGO_TRY {
            auto theory = static_cast<clingodl_theory *>(self);
            bool hv = theory->clingodl->has_value(thread_id, index);
            if (has_value != nullptr) {
                *has_value = hv;
            }
            if (hv) {
                if (symbol != nullptr) {
                    *symbol = c_cast(theory->clingodl->get_symbol(index));
                    clingo_symbol_acquire(*symbol);
                }
                if (value != nullptr) {
                    theory->clingodl->get_value(thread_id, index, *value);
                }
            }
        }
        CLINGO_CATCH;
    }

    static auto on_statistics(void *self, clingo_stats_t *stats) -> bool {
        CLINGO_TRY {
            auto theory = static_cast<clingodl_theory *>(self);
            uint64_t root = 0;
            handle_error(clingo_stats_root(stats, &root));
            auto cpp_stats = Clingo::Stats{stats, root};
            theory->clingodl->on_statistics(cpp_stats["user_step"].map(), cpp_stats["user_accu"].map());
        }
        CLINGO_CATCH;
    }
};

} // namespace

extern "C" bool clingodl_create(clingo_lib_t *lib, clingo_theory_t *theory) {
    CLINGO_TRY {
        *theory = clingo_theory_t{
            clingodl_theory::info,
            clingodl_theory::destroy,
            clingodl_theory::register_,
            clingodl_theory::rewrite_ast,
            clingodl_theory::prepare,
            clingodl_theory::register_options,
            clingodl_theory::validate_options,
            clingodl_theory::configure,
            clingodl_theory::on_model,
            clingodl_theory::on_statistics,
            clingodl_theory::lookup_symbol,
            clingodl_theory::assignment_next,
            clingodl_theory::assignment_get_value,
            nullptr,
        };
        theory->self = std::make_unique<clingodl_theory>(lib).release();
    }
    CLINGO_CATCH;
}
