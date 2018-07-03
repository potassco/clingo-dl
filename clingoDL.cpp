#include "clingoDL.h"
#include "propagator.h"
#include <vector>
#include <map>
#include <memory>
#include <iostream>

#define CLINGODL_CALLBACK_TRY try
#define CLINGODL_CALLBACK_CATCH catch (...){ return false; } return true

#define CLINGODL_TRY try
#define CLINGODL_CATCH catch (...){ return false; } return true

template <typename T>
bool init(clingo_propagate_init_t* i, void* data)
{
    CLINGODL_CALLBACK_TRY {
        PropagateInit in(i);
        static_cast<DifferenceLogicPropagator<T>*>(data)->init(in);
    }
    CLINGODL_CALLBACK_CATCH;
}

template <typename T>
bool propagate(clingo_propagate_control_t* i, const clingo_literal_t *changes, size_t size, void* data)
{
    CLINGODL_CALLBACK_TRY {
        PropagateControl in(i);
        static_cast<DifferenceLogicPropagator<T>*>(data)->propagate(in, {changes, size});
    }
    CLINGODL_CALLBACK_CATCH;
}

template <typename T>
bool undo(clingo_propagate_control_t* i, const clingo_literal_t *changes, size_t size, void* data)
{
    CLINGODL_CALLBACK_TRY {
        PropagateControl in(i);
        static_cast<DifferenceLogicPropagator<T>*>(data)->undo(in, {changes, size});
    }
    CLINGODL_CALLBACK_CATCH;
}

template <typename T>
bool check(clingo_propagate_control_t* i, void* data)
{
    CLINGODL_CALLBACK_TRY {
        PropagateControl in(i);
        static_cast<DifferenceLogicPropagator<T>*>(data)->check(in);
    }
    CLINGODL_CALLBACK_CATCH;
}


class Storage {

public:
    virtual ~Storage() {};
    virtual void create() = 0;
    virtual std::pair<clingo_propagator*, void*> getPropagator() const = 0;

};

static std::unique_ptr<Storage> storage(nullptr);
static bool strict(false);
static bool prop(false);

template<typename T>
class PropagatorStorage : public Storage {
public:
    PropagatorStorage() : clingoProp_(nullptr), diffProp_(nullptr){
    }
    void create() override {
        /// TODO: do something about stats writing, can be changed later
        static Stats s;    
        diffProp_   = std::make_unique<DifferenceLogicPropagator<T>>(s,strict,prop);
        clingoProp_ = std::make_unique<clingo_propagator>();
        clingoProp_->init = (bool (*) (clingo_propagate_init_t *, void *))init<int>;
        clingoProp_->propagate = (bool (*) (clingo_propagate_control_t *, clingo_literal_t const *, size_t, void *))propagate<int>;
        clingoProp_->undo = (bool (*) (clingo_propagate_control_t *, clingo_literal_t const *, size_t, void *))undo<int>;
        clingoProp_->check = (bool (*) (clingo_propagate_control_t *, void *))check<int>;
    }

    std::pair<clingo_propagator*, void*> getPropagator() const override {
        return std::make_pair(clingoProp_.get(), diffProp_.get());
    }

private:
    std::unique_ptr<clingo_propagator> clingoProp_;
    std::unique_ptr<DifferenceLogicPropagator<T>> diffProp_;
};

extern "C" bool theory_create_propagator(clingo_control_t* ctl, char const* option) {
    CLINGODL_TRY {
        if (strcmp(option,"int")==0) {
            storage = std::make_unique<PropagatorStorage<int>>();
        }
        else if (strcmp(option,"double")==0) {
            storage = std::make_unique<PropagatorStorage<double>>();
        }
        else { return false; }

        storage->create();

        bool ret = clingo_control_add(ctl,"base", nullptr, 0, R"(#theory dl {
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
        }.)");
        if (!ret) { return ret; }

        auto x = storage->getPropagator();
        ret = clingo_control_register_propagator(ctl, x.first, x.second, false);
        if (!ret) { return ret; }
    }
    CLINGODL_CATCH;
}

extern "C" bool theory_destroy_propagator() {
    CLINGODL_TRY {
        storage.reset(nullptr);
    }
    CLINGODL_CATCH;
}

extern "C" bool theory_add_options(clingo_options_t* options) {
    auto group = "DLPropagator";
    clingo_options_add_flag(options, group, "propagate,p", "Enable propagation.", &prop);
    clingo_options_add_flag(options, group, "strict", "Enable strict mode.", &strict);
}

#undef CLINGODL_CALLBACK_TRY
#undef CLINGODL_CALLBACK_CATCH
#undef CLINGODL_TRY
#undef CLINGODL_CATCH
