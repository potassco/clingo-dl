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

struct Storage {
public:
    virtual ~Storage() {};
    virtual void extend_model(Model &m) = 0;
    virtual size_t num_vertices() const = 0;
    virtual Symbol symbol(size_t idx) const = 0;
    virtual bool has_lower_bound(uint32_t thread_id, size_t index) const = 0;
    virtual double lower_bound(uint32_t thread_id, size_t index) const = 0;
    virtual void on_statistics(UserStatistics& step, UserStatistics &accu) = 0;
};

template<typename T>
class PropagatorStorage : public Storage {
public:
    PropagatorStorage(clingo_control_t *ctl, bool strict, PropagationMode mode)
    : prop_{step_, strict, mode} {
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

    size_t num_vertices() const override {
        return prop_.num_vertices();
    }
    Symbol symbol(size_t idx) const override {
        return prop_.symbol(idx);
    }
    void extend_model(Model &m) override {
        prop_.extend_model(m);
    }
    bool has_lower_bound(uint32_t thread_id, size_t index) const override {
        return prop_.has_lower_bound(thread_id, index);
    }
    double lower_bound(uint32_t thread_id, size_t index) const override {
        return prop_.lower_bound(thread_id, index);
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
    std::unique_ptr<Storage> storage{nullptr};
    bool strict{false};
    bool rdl{false};
    PropagationMode mode{PropagationMode::Check};
};

extern "C" bool clingodl_create_propagator(clingodl_propagator_t **prop) {
    CLINGODL_TRY { *prop = new clingodl_propagator{}; }
    CLINGODL_CATCH;
}

extern "C" bool clingodl_register_propagator(clingodl_propagator_t *prop, clingo_control_t* ctl) {
    CLINGODL_TRY {
        if (!prop->rdl) {
            prop->storage = std::make_unique<PropagatorStorage<int>>(ctl, prop->strict, prop->mode);
        }
        else {
            prop->storage = std::make_unique<PropagatorStorage<double>>(ctl, prop->strict, prop->mode);
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

extern "C" bool clingodl_add_options(clingodl_propagator_t *prop, clingo_options_t* options) {
    CLINGODL_TRY {
        char const * group = "Clingo.DL Options";
        CLINGO_CALL(clingo_options_add(options, group, "propagate",
            "Set propagation mode [no]\n"
            "    <mode>: {no,partial,full}\n"
            "      no      : No propagation; only detect conflicts\n"
            "      inverse : Check inverse constraints\n"
            "      partial : Detect some conflicting constraints\n"
            "      partial+: Detect some more conflicting constraints\n"
            "      full    : Detect all conflicting constraints",
            &parse_mode, &prop->mode, false, "mode"));
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
        prop->storage->extend_model(m);
    }
    CLINGODL_CATCH;
}

extern "C" void clingodl_assignment_begin(clingodl_propagator_t *, uint32_t, size_t *current) {
    *current = 0;
}

extern "C" bool clingodl_assignment_next(clingodl_propagator_t *prop, uint32_t thread_id, size_t *current, clingo_symbol_t *name, double* value, bool *ret) {
    CLINGODL_TRY {
        *ret = false;
        for (++*current; *current <= prop->storage->num_vertices(); ++*current) {
            size_t i = *current - 1;
            if (prop->storage->has_lower_bound(thread_id, i)) {
                *name = prop->storage->symbol(i).to_c();
                *value = prop->storage->lower_bound(thread_id, i);
                *ret = true;
                break;
            }
        }
    }
    CLINGODL_CATCH;
}

extern "C" bool clingodl_on_statistics(clingodl_propagator_t *prop, clingo_statistics_t* step, clingo_statistics_t* accu) {
    CLINGODL_TRY {
        uint64_t root_s, root_a;
        CLINGO_CALL(clingo_statistics_root(step, &root_s));
        CLINGO_CALL(clingo_statistics_root(accu, &root_a));
        UserStatistics s(step, root_s);
        UserStatistics a(accu, root_a);
        prop->storage->on_statistics(s, a);
    }
    CLINGODL_CATCH;
}

#undef CLINGODL_TRY
#undef CLINGODL_CATCH
#undef CLINGO_CALL
