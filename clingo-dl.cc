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
#include <propagator.hh>

#include <cerrno>

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

template <typename T>
bool undo(clingo_propagate_control_t const* i, const clingo_literal_t *changes, size_t size, void* data)
{
    CLINGODL_TRY {
        PropagateControl in(const_cast<clingo_propagate_control_t *>(i));
        static_cast<DifferenceLogicPropagator<T>*>(data)->undo(in, {changes, size});
    }
    CLINGODL_CATCH;
}

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
T *get_value(clingodl_value_t *value);

template<>
int *get_value<int>(clingodl_value_t *value) {
    return &value->int_number;
}

template<>
double *get_value<double>(clingodl_value_t *value) {
    return &value->double_number;
}

template<typename T>
class DLPropagatorFacade : public PropagatorFacade {
public:
    DLPropagatorFacade(clingo_control_t *ctl, bool strict, uint64_t propagate_root, uint64_t propagate_budget, PropagationMode mode)
    : prop_{step_, strict, propagate_root, propagate_budget, mode} {
        CLINGO_CALL(clingo_control_add(ctl,"base", nullptr, 0, R"(#theory dl {
term{};
constant {
  + : 1, binary, left;
  - : 1, binary, left;
  * : 2, binary, left;
  / : 2, binary, left;
  - : 3, unary
};
diff_term {- : 1, binary, left};
&diff/0 : diff_term, {<=}, constant, any;
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
        assert(index > 0 && index <= prop_.num_vertices());
        return prop_.has_lower_bound(thread_id, index);
    }
    void get_value(uint32_t thread_id, size_t index, clingodl_value_t *value) override {
        assert(index > 0 && index <= prop_.num_vertices());
        *::get_value<T>(value) = prop_.lower_bound(thread_id, index - 1);
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

struct clingodl_propagator {
    std::unique_ptr<PropagatorFacade> clingodl{nullptr};
    bool strict{false};
    bool rdl{false};
    uint64_t propagate_root{0};
    uint64_t propagate_budget{0};
    PropagationMode mode{PropagationMode::Check};
};

extern "C" bool clingodl_create_propagator(clingodl_propagator_t **prop) {
    CLINGODL_TRY { *prop = new clingodl_propagator{}; }
    CLINGODL_CATCH;
}

extern "C" bool clingodl_register_propagator(clingodl_propagator_t *prop, clingo_control_t* ctl) {
    CLINGODL_TRY {
        if (!prop->rdl) {
            prop->clingodl = std::make_unique<DLPropagatorFacade<int>>(ctl, prop->strict, prop->propagate_root, prop->propagate_budget, prop->mode);
        }
        else {
            prop->clingodl = std::make_unique<DLPropagatorFacade<double>>(ctl, prop->strict, prop->propagate_root, prop->propagate_budget, prop->mode);
        }
    }
    CLINGODL_CATCH;
}

extern "C" bool clingodl_destroy_propagator(clingodl_propagator_t *prop) {
    CLINGODL_TRY { delete prop; }
    CLINGODL_CATCH;
}

static bool iequals(char const *a, char const *b) {
    for (; *a && *b; ++a, ++b) {
        if (tolower(*a) != tolower(*b)) { return false; }
    }
    return !*a && !*b;
}
static bool parse_uint64(const char *value, void *data) {
    auto &root = *static_cast<uint64_t*>(data);
    char *end = nullptr;
    root = std::strtoull(value, &end, 10);
    return end != value && static_cast<size_t>(end - value) == std::strlen(value) && errno == 0;
}
static bool parse_mode(const char *value, void *data) {
    auto &mode = *static_cast<PropagationMode*>(data);
    if (iequals(value, "no")) {
        mode = PropagationMode::Check;
        return true;
    }
    else if (iequals(value, "inverse")) {
        mode = PropagationMode::Trivial;
        return true;
    }
    else if (iequals(value, "partial")) {
        mode = PropagationMode::Weak;
        return true;
    }
    else if (iequals(value, "partial+")) {
        mode = PropagationMode::WeakPlus;
        return true;
    }
    else if (iequals(value, "full")) {
        mode = PropagationMode::Strong;
        return true;
    }
    return false;
}

extern "C" bool clingodl_register_options(clingodl_propagator_t *prop, clingo_options_t* options) {
    CLINGODL_TRY {
        char const * group = "Clingo.DL Options";
        CLINGO_CALL(clingo_options_add(options, group, "propagate",
            "Set propagation mode [no]\n"
            "      <mode>: {no,inverse,partial,partial+,full}\n"
            "        no      : No propagation; only detect conflicts\n"
            "        inverse : Check inverse constraints\n"
            "        partial : Detect some conflicting constraints\n"
            "        partial+: Detect some more conflicting constraints\n"
            "        full    : Detect all conflicting constraints",
            &parse_mode, &prop->mode, false, "<mode>"));
        CLINGO_CALL(clingo_options_add(options, group, "propagate-root",
            "Enable full propagation below decision level <n> [0]",
            &parse_uint64, &prop->propagate_root, false, "<n>"));
        CLINGO_CALL(clingo_options_add(options, group, "propagate-budget",
            "Enable full propagation limiting to budget <n> [0]\n"
            "                            (if possible use with --propagate-root=1)\n",
            &parse_uint64, &prop->propagate_budget, false, "<n>"));
        CLINGO_CALL(clingo_options_add_flag(options, group, "rdl", "Enable support for real numbers.", &prop->rdl));
        CLINGO_CALL(clingo_options_add_flag(options, group, "strict", "Enable strict mode.", &prop->strict));
    }
    CLINGODL_CATCH;
}

extern "C" bool clingodl_validate_options(clingodl_propagator_t *prop) {
    CLINGODL_TRY {
        if (prop->strict && prop->rdl) {
            throw std::runtime_error("real difference logic not available with strict semantics");
        }
    }
    CLINGODL_CATCH;
}

extern "C" bool clingodl_on_model(clingodl_propagator_t *prop, clingo_model_t* model) {
    CLINGODL_TRY {
        Model m(model);
        prop->clingodl->extend_model(m);
    }
    CLINGODL_CATCH;
}

extern "C" bool clingodl_lookup_symbol(clingodl_propagator_t *prop_, clingo_symbol_t symbol, size_t *index) {
    return prop_->clingodl->lookup_symbol(symbol, index);
}

extern "C" clingo_symbol_t clingodl_get_symbol(clingodl_propagator_t *prop, size_t index) {
    return prop->clingodl->get_symbol(index);
}

extern "C" void clingodl_assignment_begin(clingodl_propagator_t *, uint32_t, size_t *current) {
    // Note: the first vertex is always 0 and can be skipped because its value is 0
    *current = 1;
}

extern "C" bool clingodl_assignment_next(clingodl_propagator_t *prop, uint32_t thread_id, size_t *index) {
    return prop->clingodl->next(thread_id, index);
}

extern "C" bool clingodl_assignment_has_value(clingodl_propagator_t *prop, uint32_t thread_id, size_t index) {
    return prop->clingodl->has_value(thread_id, index);
}

extern "C" void clingodl_assignment_get_value(clingodl_propagator_t *prop, uint32_t thread_id, size_t index, clingodl_value_t *value) {
    prop->clingodl->get_value(thread_id, index, value);
}

extern "C" bool clingodl_on_statistics(clingodl_propagator_t *prop, clingo_statistics_t* step, clingo_statistics_t* accu) {
    CLINGODL_TRY {
        uint64_t root_s, root_a;
        CLINGO_CALL(clingo_statistics_root(step, &root_s));
        CLINGO_CALL(clingo_statistics_root(accu, &root_a));
        UserStatistics s(step, root_s);
        UserStatistics a(accu, root_a);
        prop->clingodl->on_statistics(s, a);
    }
    CLINGODL_CATCH;
}

#undef CLINGODL_TRY
#undef CLINGODL_CATCH
#undef CLINGO_CALL
