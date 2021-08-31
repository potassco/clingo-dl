// {{{ MIT License

// Copyright 2018 Roland Kaminski, Philipp Wanko, Max Ostrowski

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

// }}}

#ifndef CLINGODL_APP_HH
#define CLINGODL_APP_HH

#include <clingo.hh>
#include <clingo-dl.h>
#include <optional>
#include <iostream>
#include <limits>
#include <cmath>

namespace ClingoDL {

//! Type used for integer values.
using int_value_t = int;

//! Helper class to rewrite logic programs to use with the clingo DL theory.
class Rewriter {
public:
    Rewriter(clingodl_theory_t *theory, clingo_program_builder_t *builder)
    : theory_{theory}
    , builder_{builder} {
    }

    //! Rewrite the given files.
    void rewrite(Clingo::StringSpan files) {
        Clingo::Detail::handle_error(clingo_ast_parse_files(files.begin(), files.size(), rewrite_, this, nullptr, nullptr, 0));
    }

    //! Rewrite the given program.
    void rewrite(char const *str) {
        Clingo::Detail::handle_error(clingo_ast_parse_string(str, rewrite_, this, nullptr, nullptr, 0));
    }

private:
    //! C callback to add a statement using the builder.
    static bool add_(clingo_ast_t *stm, void *data) {
        auto *self = static_cast<Rewriter*>(data);
        return clingo_program_builder_add(self->builder_, stm);
    }

    //! C callback to rewrite a statement and add it via the builder.
    static bool rewrite_(clingo_ast_t *stm, void *data) {
        auto *self = static_cast<Rewriter*>(data);
        return clingodl_rewrite_ast(self->theory_, stm, add_, self);
    }

    clingodl_theory_t *theory_;         //!< A theory handle to rewrite statements.
    clingo_program_builder_t *builder_; //!< The builder to add rewritten statements to.
};

//! The configuration of the optimization algorithm.
struct OptimizerConfig {
    Clingo::Symbol symbol;   //!< The variable to minimize.
    double factor{1.0};      //!< The factor by which to increase optimization step size.
    mutable size_t index{0}; //!< The index of the symbol to minimize.
    int_value_t initial{0};  //!< The initial value of the bound.
    bool has_initial{false}; //!< Whether the bound is set initially.
    bool active{false};      //!< Whether optimization is enabled.
};

//! Class providing a configurable optimization algorithm.
class Optimizer : private Clingo::SolveEventHandler {
private:
    //! Type used to store bounds.
    using Bound = std::optional<int_value_t>;

public:
    Optimizer(OptimizerConfig const &opt_cfg, Clingo::SolveEventHandler &handler, clingodl_theory_t *theory)
    : opt_cfg_{opt_cfg}
    , handler_{handler}
    , theory_{theory} {
    }

    //! Run the optimization algorithm.
    //!
    //! \note
    //! With an API extension to implement a custom enumerator, one could
    //! implement this more nicely. Right now, this implementation is
    //! restricted to the application.
    void solve(Clingo::Control& ctl) {
        Clingo::AST::with_builder(ctl, [&](Clingo::AST::ProgramBuilder &builder) {
            Rewriter rewriter{theory_, builder.to_c()};
            rewriter.rewrite(
                // add a fixed bound
                "#program __ub(s,b)."
                "&diff { s-0 } <= b."
                // add a retractable bound
                "#program __sb(s,b)."
                "#external __sb(b). [true]"
                "&diff { s-0 } <= b :- __sb(b).");
        });
        if (opt_cfg_.has_initial) {
            upper_bound_ = opt_cfg_.initial;
        }
        for (;;) {
            prepare_(ctl);
            auto ret = ctl.solve(Clingo::SymbolicLiteralSpan{}, this, false, false).get();
            if (ret.is_interrupted()) {
                break;
            }
            if (ret.is_unsatisfiable()) {
                if (search_bound_) {
                    lower_bound_ = *search_bound_ + 1;
                }
                search_bound_ = upper_bound_;
                adjust_ = 1;
                if (!lower_bound_ || *lower_bound_ > *upper_bound_) {
                    break;
                }
            }
        }
    }

private:
    //! Function to add DL specific statistics.
    void on_statistics(Clingo::UserStatistics step, Clingo::UserStatistics accu) override {
        add_stats(step);
        add_stats(accu);
        handler_.on_statistics(step, accu);
    }

    //! Add information about bounds to the given root statistics object and
    //! pass call the theory specific handler.
    void add_stats(Clingo::UserStatistics root) const {
        if (optimization || lower_bound_) {
            Clingo::UserStatistics diff = root.add_subkey("DifferenceLogic", Clingo::StatisticsType::Map);
            if (optimization) {
                diff.add_subkey("Optimization", Clingo::StatisticsType::Value).set_value(*optimization);
            }
            if (lower_bound_) {
                diff.add_subkey("Lower bound", Clingo::StatisticsType::Value).set_value(*lower_bound_);
            }
        }
    }

    //! Function to extract the current bound and pass the model to the theory.
    bool on_model(Clingo::Model &model) override {
        // update (upper) bound
        optimization = get_bound(model);
        upper_bound_ = *optimization - 1;

        // determine search bound
        double aux = *optimization - adjust_;
        if (lower_bound_ && aux <= *lower_bound_) {
            aux = *lower_bound_ + 1.0;
        }
        if (aux < std::numeric_limits<int_value_t>::min()) {
            aux = std::numeric_limits<int_value_t>::min();
        }
        search_bound_ = static_cast<int_value_t>(aux);

        // update (exponential) adujustment value
        adjust_ = adjust_ * opt_cfg_.factor;

        // pass model to theory
        handler_.on_model(model);
        return false;
    }

    //! Extract the bound from the given model.
    int_value_t get_bound(Clingo::Model &model) {
        // get bound
        if (opt_cfg_.index == 0) {
            if (!clingodl_lookup_symbol(theory_, opt_cfg_.symbol.to_c(), &opt_cfg_.index)) {
                throw std::runtime_error("variable to minimize not found");
            }
        }
        if (!clingodl_assignment_has_value(theory_, model.thread_id(), opt_cfg_.index)) {
            throw std::runtime_error("variable to minimize is unassigned");
        }
        clingodl_value_t value;
        clingodl_assignment_get_value(theory_, model.thread_id(), opt_cfg_.index, &value);
        // NOTE: minimizinig real values would require an epsilon
        if (value.type != clingodl_value_type_int) {
            throw std::runtime_error("only integer minimization is supported");
        }
        return value.int_number; // NOLINT
    }

    //! Prepare the program for solving.
    //!
    //! This adds constraints to enforce the current upper or search bound as
    //! well as removes no longer required bounds.
    void prepare_(Clingo::Control& ctl) {
        if (upper_bound_ && upper_bound_ != upper_bound_last_) {
            upper_bound_last_ = upper_bound_;
            ctl.ground({{"__ub", {opt_cfg_.symbol, Clingo::Number(*upper_bound_)}}});
        }
        if (search_bound_ != search_bound_last_) {
            if (search_bound_last_) {
                ctl.release_external(Clingo::Function("__sb", {Clingo::Number(*search_bound_last_)}));
            }
            if (search_bound_ && search_bound_ != upper_bound_) {
                search_bound_last_ = search_bound_;
                ctl.ground({{"__sb", {opt_cfg_.symbol, Clingo::Number(*search_bound_)}}});
            }
            else {
                search_bound_last_ = Bound{};
            }
        }
    }

    OptimizerConfig const &opt_cfg_;     //!< Configuration of the optimization algorithm.
    Clingo::SolveEventHandler &handler_; //!< Theory specific solve event handler.
    clingodl_theory_t *theory_;          //!< The underlying DL theory.
    Bound search_bound_;                 //!< The current (volatile) search bound.
    Bound search_bound_last_;            //!< The previous search bound.
    Bound lower_bound_;                  //!< The current lower bound (from UNSAT results).
    Bound upper_bound_;                  //!< The current upper bound (from SAT results).
    Bound upper_bound_last_;             //!< The last value of the upper bound.
    Bound optimization;                  //!< The current optimization value (last model found).
    double adjust_{1};                   //!< Amount by which to decrease seach bound.
};

} // namespace ClingoDL

#endif // CLINGODL_APP_HH
