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

#include "clingo-dl.h"
#include "propagator.h"
#include <map>

#define CLINGODL_TRY try
#define CLINGODL_CATCH catch(const std::exception& ex) {\
                           clingo_set_error(clingo_error_unknown, ex.what());\
                           return false;\
                       }\
                       catch (...){\
                           return false;\
                       }\
                       return true

#define CLINGO_CALL(x) if (!x) { throw std::runtime_error(clingo_error_message()); }

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

class Storage {

public:
    Storage() : currentVar_(0) {}
    virtual ~Storage() {};
    virtual clingo_control_t *getControl() = 0;
    virtual clingo_propagator_t *getClingoPropagator() = 0;
    virtual Propagator* getPropagator() = 0;
    virtual ExtendedPropagator* getExtendedPropagator() = 0;
    virtual void onStatistics(UserStatistics& step, UserStatistics &accu) = 0;
    char const* valueToString(double value) {
        auto x = valueToString_.find(value);
        if (x != valueToString_.end()) {
            return &(x->second[0]);
        }
        valueToString_[value].resize(snprintf( nullptr, 0, "%f.10", value));
        snprintf( &(valueToString_[value][0]), 0, "%f.10", value);
        return &(valueToString_[value][0]);
    }

    uint32_t currentThread_;
    size_t currentVar_;
    std::map<double, std::vector<char> > valueToString_;
};

static std::unique_ptr<Storage> storage(nullptr);
static bool strict(false);
static bool rdl(false);
static bool prop(false);

template<typename T>
class PropagatorStorage : public Storage {
public:
    PropagatorStorage(clingo_control_t *ctl, bool strict, bool prop)
    : Storage()
    , clingoProp_ {
        init<int>,
        propagate<int>,
        undo<int>,
        check<int>,
        nullptr
    }
    , diffProp_{step_, strict, prop}
    , ctl_(ctl) { }

    clingo_control_t *getControl() override { return ctl_; }
    clingo_propagator_t *getClingoPropagator() override { return &clingoProp_; }
    Propagator* getPropagator() override { return &diffProp_; }
    ExtendedPropagator* getExtendedPropagator() override { return &diffProp_; }

    void onStatistics(UserStatistics& step, UserStatistics &accu) override {
        accu_.accu(step_);
        addStatistics(step, step_);
        addStatistics(accu, accu_);
        step_.reset();
    }

    void addStatistics(UserStatistics& root, Stats const &stats) {
        UserStatistics diff = root.add_subkey("DifferenceLogic", StatisticsType::Map);
        diff.add_subkey("Time init(s)", StatisticsType::Value).set_value(stats.time_init.count());
        UserStatistics threads = diff.add_subkey("Thread", StatisticsType::Array);
        for (DLStats& stat : step_.dl_stats) {
            auto thread = threads.push(StatisticsType::Map);
            thread.add_subkey("Propagation(s)", StatisticsType::Value).set_value(stat.time_propagate.count());
            thread.add_subkey("Dijkstra(s)", StatisticsType::Value).set_value(stat.time_dijkstra.count());
            thread.add_subkey("Undo(s)", StatisticsType::Value).set_value(stat.time_undo.count());
            thread.add_subkey("True edges", StatisticsType::Value).set_value(stat.true_edges);
            thread.add_subkey("False edges", StatisticsType::Value).set_value(stat.false_edges);
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
    clingo_propagator_t clingoProp_;
    DifferenceLogicPropagator<T> diffProp_;
    clingo_control_t* ctl_;
};

extern "C" bool theory_create_propagator(clingo_control_t* ctl) {
    CLINGODL_TRY {
        if (!rdl) {
            storage = std::make_unique<PropagatorStorage<int>>(ctl, strict, prop);
        }
        else {
            storage = std::make_unique<PropagatorStorage<double>>(ctl, strict, prop);
        }

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

        CLINGO_CALL(clingo_control_register_propagator(ctl, storage->getClingoPropagator(),
                                             static_cast<void*>(storage->getPropagator()), false));
    }
    CLINGODL_CATCH;
}

extern "C" bool theory_destroy_propagator() {
    CLINGODL_TRY {
        //TODO: currently no way to unregister a propagator
        storage.reset(nullptr);
    }
    CLINGODL_CATCH;
}

// FIXME: storing prop, strict, and rdl in static variables is stark ugly
extern "C" bool theory_add_options(clingo_options_t* options) {
    CLINGODL_TRY {
        char const * group = "DLPropagator";
        CLINGO_CALL(clingo_options_add_flag(options, group, "propagate,p", "Enable propagation.", &prop));
        CLINGO_CALL(clingo_options_add_flag(options, group, "rdl", "Enable support for real numbers.", &rdl));
        CLINGO_CALL(clingo_options_add_flag(options, group, "strict", "Enable strict mode.", &strict));
    }
    CLINGODL_CATCH;
}

extern "C" bool theory_validate_options() {
    CLINGODL_TRY {
        if (strict && rdl) {
            throw std::runtime_error("real difference logic not available with strict semantics");
        }
    }
    CLINGODL_CATCH;
}

extern "C" bool theory_on_model(clingo_model_t* model) {
     CLINGODL_TRY {
        Model m(model);
        return storage->getExtendedPropagator()->extend_model(m);
     }
     CLINGODL_CATCH;
 }

extern "C" bool theory_assignment_first(uint32_t threadId, char const** name, char const** comp, char const** value) {
    CLINGODL_TRY {
        storage->currentThread_ = threadId;
        storage->currentVar_ = 1;
        CLINGO_CALL(theory_assignment_next(name, comp, value));
    }
    CLINGODL_CATCH;
}

extern "C" bool theory_assignment_next(char const** name, char const** comp, char const** value) {
    CLINGODL_TRY {
        auto x = storage->getExtendedPropagator();
        while (! x->hasLowerBound(storage->currentThread_, storage->currentVar_)) {
            ++storage->currentVar_;
            if (storage->currentVar_ >= x->numVars()) {
                name = comp = value = 0;
                return true;
            }
        }
        auto lb = x->lowerBound(storage->currentThread_, storage->currentVar_);
        *name = x->name(storage->currentVar_);
        static char const* c = "<=";
        *comp = c;
        *value = storage->valueToString(lb);
        ++storage->currentVar_;
    }
    CLINGODL_CATCH;
}

extern "C" bool theory_on_statistics(clingo_statistics_t* step, clingo_statistics_t* accu) {
    CLINGODL_TRY {
        uint64_t root_s, root_a;
        CLINGO_CALL(clingo_statistics_root(step, &root_s));
        CLINGO_CALL(clingo_statistics_root(accu, &root_a));
        UserStatistics s(step, root_s);
        UserStatistics a(accu, root_a);
        storage->onStatistics(s, a);
    }
    CLINGODL_CATCH;
}

#undef CLINGODL_TRY
#undef CLINGODL_CATCH
#undef CLINGO_CALL
