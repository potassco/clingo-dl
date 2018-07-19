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
bool undo(clingo_propagate_control_t* i, const clingo_literal_t *changes, size_t size, void* data)
{
    CLINGODL_TRY {
        PropagateControl in(i);
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
    virtual ~Storage() {};
    virtual void create(clingo_control_t* ctl) = 0;
    virtual clingo_control_t* getControl() const = 0;
    virtual clingo_propagator* getClingoPropagator() const = 0;
    virtual Propagator* getPropagator() const = 0;
    virtual ExtendedPropagator* getExtendedPropagator() const = 0;
    virtual void onStatisticsStep(UserStatistics& step) = 0;
    virtual void onStatisticsAccu(UserStatistics& accu) = 0;


};

static std::unique_ptr<Storage> storage(nullptr);
static bool strict(false);
static bool rdl(false);
static bool prop(false);

template<typename T>
class PropagatorStorage : public Storage {
public:
    PropagatorStorage() : clingoProp_(nullptr), diffProp_(nullptr), ctl_(nullptr){
    }
    void create(clingo_control_t* ctl) override {
        stats_      = std::make_unique<Stats>();
        diffProp_   = std::make_unique<DifferenceLogicPropagator<T>>(*(stats_.get()),strict,prop);
        clingoProp_ = std::make_unique<clingo_propagator>();
        clingoProp_->init = (bool (*) (clingo_propagate_init_t *, void *))init<int>;
        clingoProp_->propagate = (bool (*) (clingo_propagate_control_t *, clingo_literal_t const *, size_t, void *))propagate<int>;
        clingoProp_->undo = (bool (*) (clingo_propagate_control_t *, clingo_literal_t const *, size_t, void *))undo<int>;
        clingoProp_->check = (bool (*) (clingo_propagate_control_t *, void *))check<int>;
        ctl_ = ctl;
        statCalls_["Time init(s)"] = [this]() { return static_cast<double>(stats_->time_init.count()); };
        statCallsThread_["Propagation(s)"] = [this](size_t thread) { return static_cast<double>(stats_->dl_stats[thread].time_propagate.count()); };
        statCallsThread_["Dijkstra(s)"] = [this](size_t thread) { return static_cast<double>(stats_->dl_stats[thread].time_dijkstra.count()); };
        statCallsThread_["True edges"] = [this](size_t thread) { return static_cast<double>(stats_->dl_stats[thread].true_edges); };
        statCallsThread_["False edges"] = [this](size_t thread) { return static_cast<double>(stats_->dl_stats[thread].false_edges); };
        statCallsThread_["Undo(s)"] = [this](size_t thread) { return static_cast<double>(stats_->dl_stats[thread].time_undo.count()); };
    }

    clingo_control_t* getControl() const override { return ctl_; }
    clingo_propagator* getClingoPropagator() const override { return clingoProp_.get(); }
    Propagator* getPropagator() const override { return diffProp_.get(); }
    ExtendedPropagator* getExtendedPropagator() const override { return diffProp_.get(); }
    void onStatisticsStep(UserStatistics& step) override {
        assert(step.type() == StatisticsType::Map);
        for (const auto& it : statCalls_) {
            if (!step.has_subkey(it.first)) {
                auto t = step.add_subkey(it.first, StatisticsType::Value);
                t.set_value(it.second());
            } else {
                auto t = step[it.first];
                t.set_value(it.second());
            }
        }
        size_t num = 0;
        for (auto& stat : stats_->dl_stats) {
            std::string threadname("Thread[");
            threadname += std::to_string(num) + "]";
            UserStatistics thread(0,0);
            if (!step.has_subkey(threadname.c_str())) {
                thread = step.add_subkey(threadname.c_str(), StatisticsType::Map);
            } else {
                thread = step[threadname.c_str()];
            }
            for (const auto& it : statCallsThread_) {
                if (!thread.has_subkey(it.first)) {
                    auto t = thread.add_subkey(it.first, StatisticsType::Value);
                    t.set_value(it.second(num));
                } else {
                    auto t = thread[it.first];
                    t.set_value(it.second(num));
                }
            }
            ++num;
        }
        stats_->reset();
    }

    void onStatisticsAccu(UserStatistics& accu) override {
        assert(accu.type() == StatisticsType::Map);
        for (const auto& it : statCalls_) {
            if (!accu.has_subkey(it.first)) {
                auto t = accu.add_subkey(it.first, StatisticsType::Value);
                t.set_value(it.second());
            } else {
                auto t = accu[it.first];
                t.set_value(t.value() + it.second());
            }
        }
        size_t num = 0;
        for (auto& stat : stats_->dl_stats) {
            std::string threadname("Thread[");
            threadname += std::to_string(num) + "]";
            UserStatistics thread(0,0);
            if (!accu.has_subkey(threadname.c_str())) {
                thread = accu.add_subkey(threadname.c_str(), StatisticsType::Map);
            } else {
                thread = accu[threadname.c_str()];
            }
            for (const auto& it : statCallsThread_) {
                if (!thread.has_subkey(it.first)) {
                    auto t = thread.add_subkey(it.first, StatisticsType::Value);
                    t.set_value(it.second(num));
                } else {
                    auto t = thread[it.first];
                    t.set_value(t.value() + it.second(num));
                }
            }
            ++num;
        }
    }


private:
    std::unique_ptr<Stats> stats_;
    std::unique_ptr<clingo_propagator> clingoProp_;
    std::unique_ptr<DifferenceLogicPropagator<T>> diffProp_;
    clingo_control_t* ctl_;
    std::map<const char*, std::function<double()>> statCalls_;
    std::map<const char*, std::function<double(size_t)>> statCallsThread_;
};

extern "C" bool theory_create_propagator(clingo_control_t* ctl) {
    CLINGODL_TRY {
        if (!rdl) {
            storage = std::make_unique<PropagatorStorage<int>>();
        }
        else {
            storage = std::make_unique<PropagatorStorage<double>>();
        }

        storage->create(ctl);

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

extern "C" bool theory_on_statistics(clingo_statistics_t* step, clingo_statistics_t* accu) {
    CLINGODL_TRY {
        uint64_t root_a = 0;
        CLINGO_CALL(clingo_statistics_root(accu, &root_a));
        bool ret = false;
        CLINGO_CALL(clingo_statistics_map_has_subkey(accu, root_a, "DifferenceLogic", &ret));
        if (!ret) {
            uint64_t new_s = 0;
            CLINGO_CALL(clingo_statistics_map_add_subkey(accu, root_a, "DifferenceLogic", clingo_statistics_type_map, &new_s));
            root_a = new_s;
        } else {
            CLINGO_CALL(clingo_statistics_map_at(accu, root_a, "DifferenceLogic", &root_a));
        }
        UserStatistics a(accu, root_a);
        storage->onStatisticsAccu(a);

        uint64_t root_s = 0;
        CLINGO_CALL(clingo_statistics_root(step, &root_s));
        ret = false;
        CLINGO_CALL(clingo_statistics_map_has_subkey(step, root_s, "DifferenceLogic", &ret));
        if (!ret) {
            uint64_t new_s = 0;
            CLINGO_CALL(clingo_statistics_map_add_subkey(step, root_s, "DifferenceLogic", clingo_statistics_type_map, &new_s));
            root_s = new_s;
        } else {
            CLINGO_CALL(clingo_statistics_map_at(step, root_s, "DifferenceLogic", &root_s));
        }
        UserStatistics s(step, root_s);
        storage->onStatisticsStep(s);
    }
    CLINGODL_CATCH;
}

#undef CLINGODL_TRY
#undef CLINGODL_CATCH
#undef CLINGO_CALL
