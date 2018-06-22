#include "clingoDL.h"
#include "propagator.h"
#include <vector>
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
    ~PropagatorStorage() {
        for (auto i : intProps_) {
            delete i.first;
            delete i.second;
        }
        for (auto i : doubleProps_) {
            delete i.first;
            delete i.second;
        }
    }
    int create(char const* option) {
        // TODO: do something about stats writing, can be changed later
        // NOTE: this is not how to write exception-safe code (hint unique_ptr)
        static Stats s;
        if(strcmp(option,"int")==0) {
            auto c = new DifferenceLogicPropagator<int> (s,false,false);
            auto p = new clingo_propagator{
                (bool (*) (clingo_propagate_init_t *, void *))init<int>,
                (bool (*) (clingo_propagate_control_t *, clingo_literal_t const *, size_t, void *))propagate<int>,
                (bool (*) (clingo_propagate_control_t *, clingo_literal_t const *, size_t, void *))undo<int>,
                (bool (*) (clingo_propagate_control_t *, void *))check<int>
            };
            intProps_.push_back(std::make_pair(p,c));
            return (intProps_.size()-1)*2;
        }
        else if(strcmp(option,"double")==0){
            auto c = new DifferenceLogicPropagator<double> (s,false,false);
            auto p = new clingo_propagator{
                (bool (*) (clingo_propagate_init_t *, void *))init<double>,
                (bool (*) (clingo_propagate_control_t *, clingo_literal_t const *, size_t, void *))propagate<double>,
                (bool (*) (clingo_propagate_control_t *, clingo_literal_t const *, size_t, void *))undo<double>,
                (bool (*) (clingo_propagate_control_t *, void *))check<double>
            };
            doubleProps_.push_back(std::make_pair(p,c));
            return (doubleProps_.size()-1)*2+1;
        } else throw std::invalid_argument("Difference Logic only supports int and double as options");
    }

    std::pair<clingo_propagator*, void*> getPropagator(int ref) const {
        if (ref % 2 == 0) {
            ref/=2;
            return intProps_[ref];
        }
        else {
            ref= (ref-1)/2;
            return doubleProps_[ref];
        }
    }

    void destroy(int ref) {
        // vector space is currently not reused
        if (ref % 2 == 0) {
            ref/=2;
            delete intProps_[ref].first;
            delete intProps_[ref].second;
            intProps_[ref].first = 0;
            intProps_[ref].second = 0;
        }
        else {
            ref= (ref-1)/2;
            delete doubleProps_[ref].first;
            delete doubleProps_[ref].second;
            doubleProps_[ref].first = 0;
            doubleProps_[ref].second = 0;
        }
    }
private:
    std::vector<std::pair<clingo_propagator*, DifferenceLogicPropagator<int>*>>    intProps_;
    std::vector<std::pair<clingo_propagator*, DifferenceLogicPropagator<double>*>> doubleProps_;
} propStorage;

extern "C" bool theory_create_propagator(clingo_control_t* ctl, char const* option, int* ref) {
    CLINGODL_TRY {
        *ref = propStorage.create(option);

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

        auto x = propStorage.getPropagator(*ref);
        ret = clingo_control_register_propagator(ctl, x.first, x.second, false);
        if (!ret) { return ret; }
    }
    CLINGODL_CATCH;
}

extern "C" bool theory_destroy_propagator(int* ref) {
    CLINGODL_TRY {
        propStorage.destroy(*ref);
    }
    CLINGODL_CATCH;
}

#undef CLINGODL_CALLBACK_TRY
#undef CLINGODL_CALLBACK_CATCH
#undef CLINGODL_TRY
#undef CLINGODL_CATCH
