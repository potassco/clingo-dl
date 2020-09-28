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
            if (!clingodl_lookup_symbol(theory_, bound_symbol_.to_c(), &bound_index_)) {
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
            diff.add_subkey("Optimization", StatisticsType::Value).set_value(opt_.bound_value());
            diff.add_subkey("Lowerbound", StatisticsType::Value).set_value(opt_.lower_bound());
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

    class Optimizer {
    public:

        void set_symbol(Clingo::Symbol opt) {
            bound_symbol_ = opt;
        }

        void set_bound(int bound) {
            has_bound_ = true;
            bound_value_ = bound;
        }

        void set_factor(double factor) {
            factor_ = factor;
        }

        bool unproven() const {
            return has_bound_;
        }

        int bound_value() const {
            return bound_value_;
        }

        int lower_bound() const {
            return lower_bound_;
        }

        void setup(Control& ctl, clingodl_theory_t* theory) {
           ctl.with_builder([&](Clingo::ProgramBuilder &builder) {
           Rewriter rewriter{theory, builder.to_c()};
           rewriter.rewrite("#program __ub(s,b)."
                            "#external __ub(b)."
                            "&diff { s-0 } <= b :- __ub(b)."
                            "#program __lb(s,b)."
                            "&diff { 0-s } <= -b-1.");
            });
        }

        void prepare_solve(Control& ctl) {
            release_external(ctl); // release any old externals
            ubs_ = Number(bound_value_);
            if (has_bound_) {
                ctl.ground({{"__ub", {bound_symbol_, ubs_}}});
                ctl.assign_external(Function("__ub", {ubs_}), TruthValue::True);
                assigned_ = true;
            }
            has_bound_ = false;
        }

        // Note: it would be possible to add ub(b) as facts after an upper bound has
        //       been verified. We would get rid of the external. The constraint in
        //       clingo-dl could be removed nevertheless.
        void found_model(int bound) {
            upper_bound_ = bound;
            double aux = std::max(upper_bound_ - adapt_, lower_bound_+1.0);
            if (aux > std::numeric_limits<int>::max()) {
                throw std::overflow_error("Integer overflow during optimization");
            }
            bound_value_ = aux;
            if (bound_value_ < upper_bound_) {
                has_bound_ = true;
            }
            update_adapt();
        }

        void found_unsat(Control& ctl) {
            lower_bound_ = bound_value_;
            adapt_=1;
            bound_value_ = std::max(static_cast<int>(upper_bound_ - adapt_), lower_bound_);
            if (lower_bound_<bound_value_) {
                has_bound_ = true;
                ctl.ground({{"__lb", {bound_symbol_, Number(lower_bound_)}}});
            }
            update_adapt();
        }

    public:

        void release_external(Control& ctl) const {
            if (assigned_) {
                ctl.release_external(Function("__ub", {ubs_}));
            }
        }

        void update_adapt() {
            adapt_ = std::min(static_cast<double>(std::numeric_limits<int>::max()), adapt_*factor_);
        }

        Clingo::Symbol bound_symbol_;   // which symbol to optimize
        Clingo::Symbol ubs_;            // last used upper bound symbol
        int bound_value_ = 0;
        bool has_bound_ = false;
        int lower_bound_ = std::numeric_limits<int>::min();
        int upper_bound_ = std::numeric_limits<int>::max();
        double factor_ = 1;
        double adapt_ = 1;
        bool assigned_ = false;
    };

    void solve_optimal(Control& ctl) {
        // NOTE: with an API extension to implement a custom enumerator,
        //       one could implement this more nicely
        //       right now this implementation is restricted to the application
        opt_.setup(ctl, theory_);
        do
        {
            opt_.prepare_solve(ctl);
            auto h = ctl.solve(Clingo::SymbolicLiteralSpan{}, this, false, true);
            for (auto &&m : h) {
                opt_.found_model(get_bound(m));
                found_bound_ = true;
                break;
            }
            if (h.get().is_unsatisfiable()) {
                opt_.found_unsat(ctl);
            }
            if (h.get().is_interrupted()) {
                break;
            }
        }
        while (opt_.unproven());
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
        else {
            solve_optimal(ctl);
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
        opt_.set_symbol(args[0]);
        bound_symbol_ = args[0];
        if (size > 1) {
            opt_.set_bound(args[1].number());
        }
        return true;
    }

    bool parse_factor(char const *value) {
        std::stringstream strValue;
        strValue.imbue(std::locale::classic());
        strValue << value;
        double factor;
        strValue >> factor;
        if (factor < std::numeric_limits<int>::min() || factor > std::numeric_limits<int>::max()) {
            throw std::overflow_error("minimize-factor out of bounds");
        }
        opt_.set_factor(factor);
        return strValue.rdbuf()->in_avail() == 0;
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
        options.add(group, "minimize-factor",
            "Decrease optimization steps \n"
            "      <factor>   : {multiplication factor, 1=linear (default)}",
            [this](char const *value) { return parse_factor(value); });

    }

    void validate_options() override {
        CLINGO_CALL(clingodl_validate_options(theory_));
    }

private:
    clingodl_theory_t *theory_;
    Optimizer opt_;
    Clingo::Symbol bound_symbol_;
    size_t bound_index_ = 0;
    bool minimize_ = false;
    bool found_bound_ = false;
};

int main(int argc, char *argv[]) {
    ClingoDLApp app;
    return Clingo::clingo_main(app, {argv + 1, static_cast<size_t>(argc - 1)});
}


#undef CLINGO_CALL
