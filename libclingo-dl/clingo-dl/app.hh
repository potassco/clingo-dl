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
#include <iostream>
#include <limits>
#include <cmath>

class Rewriter {
public:
    Rewriter(clingodl_theory_t *theory, clingo_program_builder_t *builder)
    : theory_{theory}
    , builder_{builder} {
    }

    void rewrite(Clingo::StringSpan files) {
        Clingo::Detail::handle_error(clingo_ast_parse_files(files.begin(), files.size(), rewrite_, this, nullptr, nullptr, 0));
    }

    void rewrite(char const *str) {
        Clingo::Detail::handle_error(clingo_ast_parse_string(str, rewrite_, this, nullptr, nullptr, 0));
    }

private:
    static bool add_(clingo_ast_t *stm, void *data) {
        auto *self = static_cast<Rewriter*>(data);
        return clingo_program_builder_add(self->builder_, stm);
    }

    static bool rewrite_(clingo_ast_t *stm, void *data) {
        auto *self = static_cast<Rewriter*>(data);
        return clingodl_rewrite_ast(self->theory_, stm, add_, self);
    }

    clingodl_theory_t *theory_;
    clingo_program_builder_t *builder_;
};

class Bound {
public:
    int operator *() const {
        return value_;
    }
    explicit operator bool() const {
        return has_value_;
    }
    Bound &operator=(int value) {
        value_ = value;
        has_value_ = true;
        return *this;
    }
    friend bool operator==(Bound const &a, Bound const &b) {
        return a.value_ == b.value_ && a.has_value_ == b.has_value_;
    }
    friend bool operator!=(Bound const &a, Bound const &b) {
        return !(a == b);
    }
    friend std::ostream &operator<<(std::ostream &out, Bound const &b) {
        if (b) {
            out << *b;
        }
        else {
            out << "unset";
        }
        return out;
    }
private:
    int value_{0};
    bool has_value_{false};
};

struct OptimizerConfig {
    Clingo::Symbol symbol;
    double factor{1.0};
    mutable size_t index{0};
    int initial{0};
    bool has_initial{false};
    bool active{false};
};

class Optimizer : private Clingo::SolveEventHandler {
public:
    Optimizer(OptimizerConfig const &opt_cfg, Clingo::SolveEventHandler &handler, clingodl_theory_t *theory)
    : opt_cfg_{opt_cfg}
    , handler_{handler}
    , theory_{theory} {
    }

    // NOTE: with an API extension to implement a custom enumerator,
    //       one could implement this more nicely
    //       right now this implementation is restricted to the application
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
    void on_statistics(Clingo::UserStatistics step, Clingo::UserStatistics accu) override {
        add_stats(step);
        add_stats(accu);
        handler_.on_statistics(step, accu);
    }

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

    bool on_model(Clingo::Model &model) override {
        // update (upper) bound
        optimization = get_bound(model);
        upper_bound_ = *optimization - 1;

        // determine search bound
        double aux = *optimization - adjust_;
        if (lower_bound_ && aux <= *lower_bound_) {
            aux = *lower_bound_ + 1.0;
        }
        if (aux < std::numeric_limits<int>::min()) {
            aux = std::numeric_limits<int>::min();
        }
        search_bound_ = static_cast<int>(aux);

        // update (exponential) adujustment value
        adjust_ = adjust_ * opt_cfg_.factor;

        // pass model to theory
        handler_.on_model(model);
        return false;
    }

    int get_bound(Clingo::Model &model) {
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
        return value.int_number;
    }

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

    OptimizerConfig const &opt_cfg_;
    Clingo::SolveEventHandler &handler_;
    clingodl_theory_t *theory_;
    Bound search_bound_;
    Bound search_bound_last_;
    Bound lower_bound_;
    Bound upper_bound_;
    Bound upper_bound_last_;
    Bound optimization;
    double adjust_{1};
};

#endif // CLINGODL_APP_HH
