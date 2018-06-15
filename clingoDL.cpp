#include "clingoDL.h"
#include "propagator.h"
#include <map>

#define CLINGODL_CALLBACK_TRY try
#define CLINGODL_CALLBACK_CATCH catch (...){ return false; } return true

#define CLINGODL_TRY try
#define CLINGODL_CATCH catch (...){ return false; } return true

template <typename T>
bool init(clingo_propagate_init_t* i, void* data)
{
    CLINGODL_CALLBACK_TRY {
    PropagateInit in(i);
    reinterpret_cast<DifferenceLogicPropagator<T>*>(data)->init(in);
    }
    CLINGODL_CALLBACK_CATCH;
}

template <typename T>
bool propagate(clingo_propagate_control_t* i, const clingo_literal_t *changes, size_t size, void* data)
{
    CLINGODL_CALLBACK_TRY {
    PropagateControl in(i);
    reinterpret_cast<DifferenceLogicPropagator<T>*>(data)->propagate(in, {changes, size});
    }
    CLINGODL_CALLBACK_CATCH;
}

template <typename T>
bool undo(clingo_propagate_control_t* i, const clingo_literal_t *changes, size_t size, void* data)
{
    CLINGODL_CALLBACK_TRY {
    PropagateControl in(i);
    reinterpret_cast<DifferenceLogicPropagator<T>*>(data)->undo(in, {changes, size});
    }
    CLINGODL_CALLBACK_CATCH;
}

template <typename T>
bool check(clingo_propagate_control_t* i, void* data)
{
    CLINGODL_CALLBACK_TRY {
    PropagateControl in(i);
    reinterpret_cast<DifferenceLogicPropagator<T>*>(data)->check(in);
    }
    CLINGODL_CALLBACK_CATCH;
}


class PropagatorStorage {
public:
    std::pair<clingo_propagator*, void*> create(char const* option) {
        /// TODO: do something about stats writing, can be changed later
        static Stats s;    
        if(strcmp(option,"int")==0) {
            auto c = new DifferenceLogicPropagator<int> (s,false,false);
            auto p = new clingo_propagator{
                (bool (*) (clingo_propagate_init_t *, void *))init<int>,
                (bool (*) (clingo_propagate_control_t *, clingo_literal_t const *, size_t, void *))propagate<int>,
                (bool (*) (clingo_propagate_control_t *, clingo_literal_t const *, size_t, void *))undo<int>,
                (bool (*) (clingo_propagate_control_t *, void *))check<int>
            };
            intProps_[p] = c;
        }
        else if(strcmp(option,"double")==0){
            auto c = new DifferenceLogicPropagator<double> (s,false,false);
            auto p = new clingo_propagator{
                (bool (*) (clingo_propagate_init_t *, void *))init<double>,
                (bool (*) (clingo_propagate_control_t *, clingo_literal_t const *, size_t, void *))propagate<double>,
                (bool (*) (clingo_propagate_control_t *, clingo_literal_t const *, size_t, void *))undo<double>,
                (bool (*) (clingo_propagate_control_t *, void *))check<double>
            };
            doubleProps_[p] = c;
        } else throw std::invalid_argument("Difference Logic only supports int and double as options");
    }

    void destroy(clingo_propagator* p)
    {
        auto x = intProps_.find(p);
        if (x != intProps_.end()) {
            delete x->second;
            intProps_.erase(x);
            delete p;
            return;
        }
        auto y = doubleProps_.find(p);
        if (y != doubleProps_.end()) {
            delete y->second;
            doubleProps_.erase(y);
            delete p;
            return;
        }
    }
private:
    std::map<clingo_propagator*, DifferenceLogicPropagator<int>*>    intProps_;
    std::map<clingo_propagator*, DifferenceLogicPropagator<double>*> doubleProps_;
} propStorage;

extern "C" bool theory_create_propagator(char const* option, clingo_propagator** p, void** prop_data) {
    CLINGODL_TRY {
    auto x = propStorage.create(option);
    *p = x.first;
    *prop_data = x.second;
    }
    CLINGODL_CATCH;
}

extern bool theory_destroy_propagator(clingo_propagator* p) {
    CLINGODL_TRY {
    propStorage.destroy(p);
    }
    CLINGODL_CATCH;
}

#undef CLINGODL_CALLBACK_TRY
#undef CLINGODL_CALLBACK_CATCH
#undef CLINGODL_TRY
#undef CLINGODL_CATCH
