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

#include <clingo-dl.h>
#include <clingo-dl/propagator.hh>

#include <sstream>

#define CLINGODL_TRY try
#define CLINGODL_CATCH catch (...){ Clingo::Detail::handle_cxx_error(); return false; } return true

#define CLINGO_CALL(x) Clingo::Detail::handle_error(x)

template <typename T>
bool init(clingo_propagate_init_t* i, void* data)
{
    CLINGODL_TRY {
        PropagateInit in(i);
        static_cast<DifferenceLogicPropagator<T>*>(data)->init(in);
    }
    CLINGODL_CATCH;
}

template <typename T>
bool propagate(clingo_propagate_control_t* i, const clingo_literal_t *changes, size_t size, void* data)
{
    CLINGODL_TRY {
        PropagateControl in(i);
        static_cast<DifferenceLogicPropagator<T>*>(data)->propagate(in, {changes, size});
    }
    CLINGODL_CATCH;
}

#if CLINGO_VERSION_MAJOR*1000 + CLINGO_VERSION_MINOR >= 5005
template <typename T>
void undo(clingo_propagate_control_t const* i, const clingo_literal_t *changes, size_t size, void* data)
{
    PropagateControl in(const_cast<clingo_propagate_control_t *>(i));
    static_cast<DifferenceLogicPropagator<T>*>(data)->undo(in, {changes, size});
}
#else
template <typename T>
bool undo(clingo_propagate_control_t const* i, const clingo_literal_t *changes, size_t size, void* data)
{
    CLINGODL_TRY {
        PropagateControl in(const_cast<clingo_propagate_control_t *>(i));
        static_cast<DifferenceLogicPropagator<T>*>(data)->undo(in, {changes, size});
    }
    CLINGODL_CATCH;
}
#endif

template <typename T>
bool check(clingo_propagate_control_t* i, void* data)
{
    CLINGODL_TRY {
        PropagateControl in(i);
        static_cast<DifferenceLogicPropagator<T>*>(data)->check(in);
    }
    CLINGODL_CATCH;
}

struct PropagatorFacade {
public:
    virtual ~PropagatorFacade() {};
    virtual bool lookup_symbol(clingo_symbol_t name, size_t *index) = 0;
    virtual clingo_symbol_t get_symbol(size_t index) = 0;
    virtual bool has_value(uint32_t thread_id, size_t index) = 0;
    virtual void get_value(uint32_t thread_id, size_t index, clingodl_value_t *value) = 0;
    virtual bool next(uint32_t thread_id, size_t *current) = 0;
    virtual void extend_model(Model &m) = 0;
    virtual void on_statistics(UserStatistics& step, UserStatistics &accu) = 0;
};

template<typename T>
void set_value(clingodl_value_t *variant, T value);

template<>
void set_value<int>(clingodl_value_t *variant, int value) {
    variant->type = clingodl_value_type_int;
    variant->int_number = value;
}

template<>
void set_value<double>(clingodl_value_t *variant, double value) {
    variant->type = clingodl_value_type_double;
    variant->double_number = value;
}

template<typename T>
class DLPropagatorFacade : public PropagatorFacade {
public:
    DLPropagatorFacade(clingo_control_t *ctl, PropagatorConfig const &conf)
    : prop_{step_, conf} {
        CLINGO_CALL(clingo_control_add(ctl,"base", nullptr, 0, R"(#theory dl {
term {
  + : 1, binary, left;
  - : 1, binary, left;
  * : 2, binary, left;
  / : 2, binary, left;
  - : 3, unary
};
&__diff_h/0 : term, {<=,>=,<,>,=,!=}, term, head;
&__diff_b/0 : term, {<=,>=,<,>,=,!=}, term, body;
&show_assignment/0 : term, directive
}.)"));
        static clingo_propagator_t prop = {
            init<T>,
            propagate<T>,
            undo<T>,
            check<T>,
            nullptr
        };
        CLINGO_CALL(clingo_control_register_propagator(ctl, &prop, &prop_, false));
    }

    bool lookup_symbol(clingo_symbol_t name, size_t *index) override {
        *index = prop_.lookup(name) + 1;
        return *index <= prop_.num_vertices();
    }

    clingo_symbol_t get_symbol(size_t index) override {
        return prop_.symbol(index - 1).to_c();
    }

    bool has_value(uint32_t thread_id, size_t index) override {
        return prop_.has_lower_bound(thread_id, index - 1);
    }
    void get_value(uint32_t thread_id, size_t index, clingodl_value_t *value) override {
        assert(index > 0 && index <= prop_.num_vertices());
        set_value(value, prop_.lower_bound(thread_id, index - 1));
    }

    bool next(uint32_t thread_id, size_t *current) override {
        for (++*current; *current <= prop_.num_vertices(); ++*current) {
            if (prop_.has_lower_bound(thread_id, *current - 1)) {
                return true;
            }
        }
        return false;
    }
    void extend_model(Model &m) override {
        prop_.extend_model(m);
    }
    void on_statistics(UserStatistics& step, UserStatistics &accu) override {
        accu_.accu(step_);
        add_statistics(step, step_);
        add_statistics(accu, accu_);
        step_.reset();
    }

    void add_statistics(UserStatistics& root, Stats const &stats) {
        UserStatistics diff = root.add_subkey("DifferenceLogic", StatisticsType::Map);
        diff.add_subkey("Time init(s)", StatisticsType::Value).set_value(stats.time_init.count());
        diff.add_subkey("CCs", StatisticsType::Value).set_value(stats.ccs);
        diff.add_subkey("Mutexes", StatisticsType::Value).set_value(stats.mutexes);
        diff.add_subkey("Edges", StatisticsType::Value).set_value(stats.edges);
        diff.add_subkey("Variables", StatisticsType::Value).set_value(stats.variables);
        UserStatistics threads = diff.add_subkey("Thread", StatisticsType::Array);
        threads.ensure_size(stats.dl_stats.size(), StatisticsType::Map);
        auto it = threads.begin();
        for (DLStats const& stat : stats.dl_stats) {
            auto thread = *it++;
            thread.add_subkey("Propagation(s)", StatisticsType::Value).set_value(stat.time_propagate.count());
            thread.add_subkey("Dijkstra(s)", StatisticsType::Value).set_value(stat.time_dijkstra.count());
            thread.add_subkey("Undo(s)", StatisticsType::Value).set_value(stat.time_undo.count());
            thread.add_subkey("True edges", StatisticsType::Value).set_value(stat.true_edges);
            thread.add_subkey("False edges", StatisticsType::Value).set_value(stat.false_edges);
            thread.add_subkey("False edges (inverse)", StatisticsType::Value).set_value(stat.false_edges_trivial);
            thread.add_subkey("False edges (partial)", StatisticsType::Value).set_value(stat.false_edges_weak);
            thread.add_subkey("False edges (partial+)", StatisticsType::Value).set_value(stat.false_edges_weak_plus);
            thread.add_subkey("Edges added", StatisticsType::Value).set_value(stat.edges_added);
            thread.add_subkey("Edges skipped", StatisticsType::Value).set_value(stat.edges_skipped);
            thread.add_subkey("Edges propagated", StatisticsType::Value).set_value(stat.edges_propagated);
            thread.add_subkey("Cost consistency", StatisticsType::Value).set_value(stat.propagate_cost_add);
            thread.add_subkey("Cost forward", StatisticsType::Value).set_value(stat.propagate_cost_from);
            thread.add_subkey("Cost backward", StatisticsType::Value).set_value(stat.propagate_cost_to);
        }
    }

private:
    Stats step_;
    Stats accu_;
    DifferenceLogicPropagator<T> prop_;
};

struct clingodl_theory {
    std::unique_ptr<PropagatorFacade> clingodl{nullptr};
    PropagatorConfig config;
    bool rdl;
    bool shift_constraints{true};
};

extern "C" bool clingodl_create(clingodl_theory_t **theory) {
    CLINGODL_TRY { *theory = new clingodl_theory{}; }
    CLINGODL_CATCH;
}

extern "C" bool clingodl_register(clingodl_theory_t *theory, clingo_control_t* ctl) {
    CLINGODL_TRY {
        if (!theory->rdl) {
            theory->clingodl = std::make_unique<DLPropagatorFacade<int>>(ctl, theory->config);
        }
        else {
            theory->clingodl = std::make_unique<DLPropagatorFacade<double>>(ctl, theory->config);
        }
    }
    CLINGODL_CATCH;
}

extern "C" bool clingodl_rewrite_statement(clingodl_theory_t *theory, clingo_ast_statement_t const *stm, clingodl_rewrite_callback_t add, void *data) {
    CLINGODL_TRY {
        Clingo::StatementCallback cb = [&](Clingo::AST::Statement &&stm) {
            transform(std::move(stm), [add, data](Clingo::AST::Statement &&stm){
                Clingo::AST::Detail::ASTToC visitor;
                auto x = stm.data.accept(visitor);
                x.location = stm.location;
                CLINGO_CALL(add(&x, data));
            }, theory->shift_constraints);
        };
        Clingo::AST::Detail::convStatement(stm, cb);
    }
    CLINGODL_CATCH;
}

extern "C" bool clingodl_prepare(clingodl_theory_t *, clingo_control_t *) {
    return true;
}

extern "C" bool clingodl_destroy(clingodl_theory_t *theory) {
    CLINGODL_TRY { delete theory; }
    CLINGODL_CATCH;
}

static char const *iequals_pre(char const *a, char const *b) {
    for (; *a && *b; ++a, ++b) {
        if (tolower(*a) != tolower(*b)) { return nullptr; }
    }
    return *b ? nullptr : a;
}
static bool iequals(char const *a, char const *b) {
    a = iequals_pre(a, b);
    return a && !*a;
}
static char const *parse_uint64_pre(const char *value, void *data) {
    auto &res = *static_cast<uint64_t*>(data);
    char const *it = value;
    res = 0;

    for (; *it; ++it) {
        if ('0' <= *it && *it <= '9') {
            auto tmp = res;
            res *= 10;
            res += *it - '0';
            if (res < tmp) { return nullptr; }
        }
        else { break; }
    }

    return value != it ? it : nullptr;
}
static bool parse_uint64(const char *value, void *data) {
    value = parse_uint64_pre(value, data);
    return value && !*value;
}

template <typename F, typename G>
bool set_config(char const *value, void *data, F f, G g) {
    try {
        auto &config = *static_cast<PropagatorConfig*>(data);
        uint64_t id = 0;
        if (*value == '\0') {
            f(config);
            return true;
        }
        else if (*value == ',' && parse_uint64(value + 1, &id) && id < 64) {
            g(config.ensure(id));
            return true;
        }
    }
    catch (...) { }
    return false;
}

static bool parse_root(const char *value, void *data) {
    uint64_t x = 0;
    return (value = parse_uint64_pre(value, &x)) && set_config(value, data,
        [x](PropagatorConfig &config) { config.propagate_root = x; },
        [x](ThreadConfig &config) { config.propagate_root = {true, x}; });
}
static bool parse_budget(const char *value, void *data) {
    uint64_t x = 0;
    return (value = parse_uint64_pre(value, &x)) && set_config(value, data,
        [x](PropagatorConfig &config) { config.propagate_budget = x; },
        [x](ThreadConfig &config) { config.propagate_budget = {true, x}; });
}
static bool parse_mutex(const char *value, void *data) {
    auto &pc = *static_cast<PropagatorConfig*>(data);
    uint64_t x = 0;
    if (!(value = parse_uint64_pre(value, &x))) { return false; }
    pc.mutex_size = x;
    if (*value == '\0') {
        pc.mutex_cutoff = 10 * x;
        return true;
    }
    if (*value == ',') {
        if (!parse_uint64(value+1, &x)) { return false; }
        pc.mutex_cutoff = x;
    }
    return true;
}
static bool parse_mode(const char *value, void *data) {
    PropagationMode mode = PropagationMode::Check;
    char const *rem = nullptr;
    if ((rem = iequals_pre(value, "no"))) {
        mode = PropagationMode::Check;
    }
    else if ((rem = iequals_pre(value, "inverse"))) {
        mode = PropagationMode::Trivial;
    }
    else if ((rem = iequals_pre(value, "partial+"))) {
        mode = PropagationMode::WeakPlus;
    }
    else if ((rem = iequals_pre(value, "partial"))) {
        mode = PropagationMode::Weak;
    }
    else if ((rem = iequals_pre(value, "full"))) {
        mode = PropagationMode::Strong;
    }
    return rem && set_config(rem, data,
        [mode](PropagatorConfig &config) { config.mode = mode; },
        [mode](ThreadConfig &config) { config.mode = {true, mode}; });
}
static bool parse_sort(const char *value, void *data) {
    SortMode sort = SortMode::Weight;
    char const *rem = nullptr;
    if ((rem = iequals_pre(value, "no"))) {
        sort = SortMode::No;
    }
    else if ((rem = iequals_pre(value, "weight-reversed"))) {
        sort = SortMode::WeightRev;
    }
    else if ((rem = iequals_pre(value, "weight"))) {
        sort = SortMode::Weight;
    }
    else if ((rem = iequals_pre(value, "potential-reversed"))) {
        sort = SortMode::PotentialRev;
    }
    else if ((rem = iequals_pre(value, "potential"))) {
        sort = SortMode::Potential;
    }
    return rem && set_config(rem, data,
        [sort](PropagatorConfig &config) { config.sort_edges = sort; },
        [sort](ThreadConfig &config) { config.sort_edges = {true, sort}; });
}

static bool parse_bool(const char *value, void *data) {
    auto &result = *static_cast<bool*>(data);
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

static bool check_parse(char const *key, bool ret) {
    if (!ret) {
        std::ostringstream msg;
        msg << "invalid value for '" << key << "'";
        clingo_set_error(clingo_error_runtime, msg.str().c_str());
    }
    return ret;
}

extern "C" bool clingodl_configure(clingodl_theory_t *theory, char const *key, char const *value) {
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
        std::ostringstream msg;
        msg << "invalid configuration key '" << key << "'";
        clingo_set_error(clingo_error_runtime, msg.str().c_str());
        return false;
    }
    CLINGODL_CATCH;
}

extern "C" bool clingodl_register_options(clingodl_theory_t *theory, clingo_options_t* options) {
    CLINGODL_TRY {
        char const * group = "Clingo.DL Options";
        CLINGO_CALL(clingo_options_add(options, group, "propagate",
            "Set propagation mode [no]\n"
            "      <mode>  : {no,inverse,partial,partial+,full}[,<thread>]\n"
            "        no      : No propagation; only detect conflicts\n"
            "        inverse : Check inverse constraints\n"
            "        partial : Detect some conflicting constraints\n"
            "        partial+: Detect some more conflicting constraints\n"
            "        full    : Detect all conflicting constraints\n"
            "      <thread>: Restrict to thread",
            &parse_mode, &theory->config, true, "<mode>"));
        CLINGO_CALL(clingo_options_add(options, group, "propagate-root",
            "Enable full propagation below decision level [0]\n"
            "      <arg>   : <n>[,<thread>]\n"
            "      <n>     : Upper bound for decision level\n"
            "      <thread>: Restrict to thread",
            &parse_root, &theory->config, true, "<arg>"));
        CLINGO_CALL(clingo_options_add(options, group, "propagate-budget",
            "Enable full propagation limiting to budget [0]\n"
            "      <arg>   : <n>[,<thread>]\n"
            "      <n>     : Budget roughly corresponding to cost of consistency checks\n"
            "                (if possible use with --propagate-root greater 0)\n"
            "      <thread>: Restrict to thread",
            &parse_budget, &theory->config, true, "<arg>"));
        CLINGO_CALL(clingo_options_add(options, group, "add-mutexes",
            "Add mutexes in a preprocessing step [0]\n"
            "      <arg>   : <max>[,<cut>]\n"
            "      <max>   : Maximum size of mutexes to add\n"
            "      <cut>   : Limit costs to calculate mutexes",
            &parse_mutex, &theory->config, true, "<arg>"));
        CLINGO_CALL(clingo_options_add(options, group, "sort-edges",
            "Sort edges for propagation [weight]\n"
            "      <arg>   : {no, weight, weight-reversed, potential, potential-reversed}\n"
            "        no                 : No sorting\n"
            "        weight             : Sort by edge weight\n"
            "        weight-reversed    : Sort by negative edge weight\n"
            "        potential          : Sort by relative potential\n"
            "        potential-reversed : Sort by relative negative potential",
            &parse_sort, &theory->config, true, "<arg>"));
        CLINGO_CALL(clingo_options_add_flag(options, group, "rdl", "Enable support for real numbers", &theory->rdl));
    }
    CLINGODL_CATCH;
}

extern "C" bool clingodl_validate_options(clingodl_theory_t *theory) {
    return true;
}

extern "C" bool clingodl_on_model(clingodl_theory_t *theory, clingo_model_t* model) {
    CLINGODL_TRY {
        Model m(model);
        theory->clingodl->extend_model(m);
    }
    CLINGODL_CATCH;
}

extern "C" bool clingodl_lookup_symbol(clingodl_theory_t *theory, clingo_symbol_t symbol, size_t *index) {
    return theory->clingodl->lookup_symbol(symbol, index);
}

extern "C" clingo_symbol_t clingodl_get_symbol(clingodl_theory_t *theory, size_t index) {
    return theory->clingodl->get_symbol(index);
}

extern "C" void clingodl_assignment_begin(clingodl_theory_t *, uint32_t, size_t *current) {
    // Note: the first vertex is always 0 and can be skipped because its value is 0
    *current = 1;
}

extern "C" bool clingodl_assignment_next(clingodl_theory_t *theory, uint32_t thread_id, size_t *index) {
    return theory->clingodl->next(thread_id, index);
}

extern "C" bool clingodl_assignment_has_value(clingodl_theory_t *theory, uint32_t thread_id, size_t index) {
    return theory->clingodl->has_value(thread_id, index);
}

extern "C" void clingodl_assignment_get_value(clingodl_theory_t *theory, uint32_t thread_id, size_t index, clingodl_value_t *value) {
    theory->clingodl->get_value(thread_id, index, value);
}

extern "C" bool clingodl_on_statistics(clingodl_theory_t *theory, clingo_statistics_t* step, clingo_statistics_t* accu) {
    CLINGODL_TRY {
        uint64_t root_s, root_a;
        CLINGO_CALL(clingo_statistics_root(step, &root_s));
        CLINGO_CALL(clingo_statistics_root(accu, &root_a));
        UserStatistics s(step, root_s);
        UserStatistics a(accu, root_a);
        theory->clingodl->on_statistics(s, a);
    }
    CLINGODL_CATCH;
}

#undef CLINGODL_TRY
#undef CLINGODL_CATCH
#undef CLINGO_CALL
