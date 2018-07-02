#include "clingoDL.h"
#include "propagator.h"
#include <vector>
#include <map>
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
/*
template<typename T>
class Prop;

template<>
class Prop<int>
{
protected:
    
}
*/

Storage* create(char const* option) {
    if(strcmp(option,"int")==0) {
        return new PropagatorStorage<int>();
            }
    else if(strcmp(option,"double")==0){
        int ref = (doubleProps_.size()-1)*2+1;
        auto fs = strict_.find(ref);
        if (fs != strict_.end()) strict_[ref] = false;
        auto fp = props_.find(ref);
        if (fp != props_.end()) props_[ref] = false;
        auto c = new DifferenceLogicPropagator<double> (s,strict_[ref],props_[ref]);
        auto p = new clingo_propagator{
            (bool (*) (clingo_propagate_init_t *, void *))init<double>,
            (bool (*) (clingo_propagate_control_t *, clingo_literal_t const *, size_t, void *))propagate<double>,
            (bool (*) (clingo_propagate_control_t *, clingo_literal_t const *, size_t, void *))undo<double>,
            (bool (*) (clingo_propagate_control_t *, void *))check<double>
        };
        doubleProps_.push_back(std::make_pair(p,c));
        return ref;
    } else throw std::invalid_argument("Difference Logic only supports int and double as options");
}

class Storage {

public:
    virtual void create(char const* option) = 0;
    virtual std::pair<clingo_propagator*, void*> getPropagator() const = 0;

};

static Storage* storage(nullptr);
static bool strict(false);
static bool prop(false);

template<typename T>
class PropagatorStorage : public Storage {
public:
    PropagatorStorage() : clingoProp_(nullptr), diffProp_(nullptr), strict(false), prop_(false) {
    }
    ~PropagatorStorage() {
        delete clingoProp_;
        delete diffProp_;
    }
    void create(char const* option) override {
        /// TODO: do something about stats writing, can be changed later
        static Stats s;    
        diffProp_   = new DifferenceLogicPropagator<T> (s,strict,prop]);
        clingoProp_ = new clingo_propagator{
            (bool (*) (clingo_propagate_init_t *, void *))init<int>,
            (bool (*) (clingo_propagate_control_t *, clingo_literal_t const *, size_t, void *))propagate<int>,
            (bool (*) (clingo_propagate_control_t *, clingo_literal_t const *, size_t, void *))undo<int>,
            (bool (*) (clingo_propagate_control_t *, void *))check<int>
        };
    }

    std::pair<clingo_propagator*, void*> getPropagator() const override {
        return std::make_pair(clingoProp_, diffProp_);
    }

private:
    clingo_propagator* clingoProp_;
    DifferenceLogicPropagator<T>* diffProp_;
    bool strict_; 
    bool props_;
} propStorage;

extern "C" bool theory_create_propagator(clingo_control_t* ctl, char const* option, int* ref) {
    CLINGODL_TRY {
        storage = create(option);
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
        /// TODO: error handling function ?
        if (!ret) {
            std::cout << "control_add failed" << std::endl;
            return ret;
        }
    
        auto x = storage.getPropagator();
        ret = clingo_control_register_propagator(ctl, x.first, x.second, false);
        if (!ret) {
            std::cout << "control_add failed" << std::endl;
            return ret;
        }
        return ret;
    }
    CLINGODL_CATCH;
}

extern "C" bool theory_destroy_propagator() {
    CLINGODL_TRY {
        delete storage;
    }
    CLINGODL_CATCH;
}

extern "C" bool theory_add_options(clingo_options_t* options) {
    propStorage.addOptions(options);
    clingo_options_add_flag(options, group, "propagate,p", "Enable propagation.", &prop);
    clingo_options_add_flag(options, group, "strict", "Enable strict mode.", &strict);
}

#undef CLINGODL_CALLBACK_TRY
#undef CLINGODL_CALLBACK_CATCH
#undef CLINGODL_TRY
#undef CLINGODL_CATCH
