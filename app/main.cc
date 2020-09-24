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

#include <clingo.hh>
#include <clingo-dl.h>
#include <sstream>
#include <fstream>
#include <cmath>

using namespace Clingo;

#define CLINGO_CALL(x) Clingo::Detail::handle_error(x)

class Rewriter {
public:
    Rewriter(clingodl_theory_t *theory, clingo_program_builder_t *builder)
    : theory_{theory}
    , builder_{builder} {
    }

    void rewrite(StringSpan files) {
        CLINGO_CALL(clingo_parse_files(files.begin(), files.size(), rewrite_, this, nullptr, nullptr, 0));
    }

    void rewrite(char const *str) {
        CLINGO_CALL(clingo_parse_program(str, rewrite_, this, nullptr, nullptr, 0));
    }

private:
    static bool add_(clingo_ast_statement_t const *stm, void *data) {
        auto *self = static_cast<Rewriter*>(data);
        return clingo_program_builder_add(self->builder_, stm);
    }

    static bool rewrite_(clingo_ast_statement_t const *stm, void *data) {
        auto *self = static_cast<Rewriter*>(data);
        return clingodl_rewrite_statement(self->theory_, stm, add_, self);
    }

    clingodl_theory_t *theory_;
    clingo_program_builder_t *builder_;
};


class ClingoDLApp : public Clingo::Application, private SolveEventHandler {
public:
    ClingoDLApp() {
        CLINGO_CALL(clingodl_create(&theory_));
    }
    ~ClingoDLApp() { clingodl_destroy(theory_); }
    char const *program_name() const noexcept override { return "clingo-dl"; }
    char const *version() const noexcept override { return CLINGODL_VERSION; }
    bool on_model(Model &model) override {
        CLINGO_CALL(clingodl_on_model(theory_, model.to_c()));
        return true;
    }

    int get_bound(Model const &m) {
        if (!bound_index_) {
            if (!clingodl_lookup_symbol(theory_, bound_symbol.to_c(), &bound_index_)) {
                throw std::runtime_error("variable to minimize not found");
            }
        }
        if (!clingodl_assignment_has_value(theory_, m.thread_id(), bound_index_)) {
            throw std::runtime_error("variable to minimize is unassigned");
        }
        clingodl_value_t value;
        clingodl_assignment_get_value(theory_, m.thread_id(), bound_index_, &value);
        // NOTE: minimizinig real values would require an epsilon
        if (value.type != clingodl_value_type_int) {
            throw std::runtime_error("only integer minimization is supported");
        }
        return value.int_number;

    }

    void add_stats(UserStatistics root) {
        if (found_bound_) {
            UserStatistics diff = root.add_subkey("DifferenceLogic", StatisticsType::Map);
            diff.add_subkey("Optimization", StatisticsType::Value).set_value(bound_value_);
            if (lower_bound_ > std::numeric_limits<int>::min()) {
                diff.add_subkey("Lowerbound", StatisticsType::Value).set_value(lower_bound_);
            }
        }
    }

    void on_statistics(UserStatistics step, UserStatistics accu) override {
        add_stats(step);
        add_stats(accu);
        CLINGO_CALL(clingodl_on_statistics(theory_, step.to_c(), accu.to_c()));
    }

    void solve_satisfiability(Control& ctl) {
        ctl.solve(Clingo::SymbolicLiteralSpan{}, this, false, false).get();
    }

    void solve_optimal_linear(Control& ctl) {
        // NOTE: with an API extension to implement a custom enumerator,
        //       one could implement this more nicely
        //       right now this implementation is restricted to the application
        ctl.with_builder([&](Clingo::ProgramBuilder &builder) {
           Rewriter rewriter{theory_, builder.to_c()};
           rewriter.rewrite("#program __bound(s,b)."
                            "&diff { s-0 } <= b."
                            "#program base.");
        });
        do
        {
            if (has_bound_) {
                ctl.ground({{"__bound", {bound_symbol, Number(bound_value_ - 1)}}});
            }
            has_bound_ = false;
            auto h = ctl.solve(Clingo::SymbolicLiteralSpan{}, this, false, true);
            for (auto &&m : h) {
                bound_value_ = get_bound(m);
                has_bound_ = true;
                found_bound_ = true;
                break;
            }
            if (h.get().is_interrupted()) {
                break;
            }
        }
        while (has_bound_);
    }

    void solve_optimal_exponential(Control& ctl) {
        ctl.with_builder([&](Clingo::ProgramBuilder &builder) {
           Rewriter rewriter{theory_, builder.to_c()};
           rewriter.rewrite("#program __ub(s,b)."
                            "#external ub(b)."
                            "&diff { s-0 } <= b :- ub(b)."
                            "#program __lb(s,b)."
                            "#external lb(b)."
                            "&diff { 0-s } <= -b-1 :- lb(b).");
        });

        int upper_bound = std::numeric_limits<int>::max();
        double adapt = 1;

        do
        {
            if (has_bound_) {
                ctl.ground({{"__ub", {bound_symbol, Number(bound_value_)}}});
                ctl.assign_external(Function("ub", {Number(bound_value_)}), TruthValue::True);
            }
            auto h = ctl.solve(Clingo::SymbolicLiteralSpan{}, this, false, true);
            has_bound_ = false;
            for (auto &&m : h) {
                upper_bound = get_bound(m);
                bound_value_ = std::max(static_cast<int>(upper_bound - adapt),lower_bound_+1);
                if (bound_value_ < upper_bound) {
                    has_bound_ = true;
                }
                found_bound_ = true;
                break;
            }
            if (h.get().is_unsatisfiable()) {
                lower_bound_ = bound_value_;
                ctl.ground({{"__lb", {bound_symbol, Number(lower_bound_)}}});
                ctl.assign_external(Function("lb", {Number(lower_bound_)}), TruthValue::True);
                adapt=1;
                bound_value_ = std::max(static_cast<int>(upper_bound - adapt), lower_bound_);
                if (lower_bound_<bound_value_) {
                    has_bound_ = true;
                    ctl.release_external(Function("ub", {Number(lower_bound_)}));
                }
            }
            adapt = std::min(static_cast<double>(std::numeric_limits<int>::max()), adapt*2);
            if (h.get().is_interrupted()) {
                break;
            }
        }
        while (has_bound_);

    }

    void main(Control &ctl, StringSpan files) override {
        CLINGO_CALL(clingodl_register(theory_, ctl.to_c()));

        ctl.with_builder([&](Clingo::ProgramBuilder &builder) {
            Rewriter rewriter{theory_, builder.to_c()};
            rewriter.rewrite(files);
        });

        ctl.ground({{"base", {}}});
        if (!minimize_) {
            solve_satisfiability(ctl);
        }
        else if (strategy_linear_) {
            solve_optimal_linear(ctl);
        }
        else {
            solve_optimal_exponential(ctl);
        }

    }

    bool parse_bound(char const *value) {
        std::ostringstream oss;
        oss << "(" << value << ",)";
        auto term = Clingo::parse_term(oss.str().c_str());
        auto args = term.arguments();
        auto size = args.size();
        if (args.empty() || size > 2 || (size > 1 && args[1].type() != Clingo::SymbolType::Number)) {
            return false;
        }
        minimize_ = true;
        bound_symbol = args[0];
        if (size > 1) {
            has_bound_ = true;
            bound_value_ = args[1].number();
        }
        return true;
    }

    bool parse_strategy(char const *value) {
        if (strcmp(value,"linear") == 0) {
            strategy_linear_ = true;
        }
        else if (strcmp(value,"exponential") == 0) {
            strategy_linear_ = false;
        }
        else {
            return false;
        }
        return true;
    }

    void register_options(ClingoOptions &options) override {
        CLINGO_CALL(clingodl_register_options(theory_, options.to_c()));
        char const * group = "Clingo.DL Options";
        options.add(group, "minimize-variable",
            "Minimize the given variable\n"
            "      <arg>   : <variable>[,<initial>]\n"
            "      <variable>: the variable to minimize\n"
            "      <initial> : upper bound for the variable",
            [this](char const *value) { return parse_bound(value); });
        options.add(group, "minimize-strategy",
            "Minimize with the given strategy\n"
            "      <arg>   : {linear (default), exponential}",
            [this](char const *value) { return parse_strategy(value); });

    }

    void validate_options() override {
        CLINGO_CALL(clingodl_validate_options(theory_));
    }

private:
    clingodl_theory_t *theory_;
    Clingo::Symbol bound_symbol;
    size_t bound_index_ = 0;
    int bound_value_ = 0;
    int lower_bound_ = std::numeric_limits<int>::min();
    bool minimize_ = false;
    bool has_bound_ = false;
    bool found_bound_ = false;
    bool strategy_linear_ = true;
};

int main(int argc, char *argv[]) {
    ClingoDLApp app;
    return Clingo::clingo_main(app, {argv + 1, static_cast<size_t>(argc - 1)});
}


#undef CLINGO_CALL
