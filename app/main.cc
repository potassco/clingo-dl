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
#include <clingo-dl/app.hh>
#include <sstream>
#include <fstream>
#include <limits>

using Clingo::Detail::handle_error;

class ClingoDLApp : public Clingo::Application, private Clingo::SolveEventHandler {
public:
    ClingoDLApp() {
        handle_error(clingodl_create(&theory_));
    }
    ~ClingoDLApp() { clingodl_destroy(theory_); }
    char const *program_name() const noexcept override { return "clingo-dl"; }
    char const *version() const noexcept override { return CLINGODL_VERSION; }
    bool on_model(Clingo::Model &model) override {
        handle_error(clingodl_on_model(theory_, model.to_c()));
        return true;
    }

    void on_statistics(Clingo::UserStatistics step, Clingo::UserStatistics accu) override {
        handle_error(clingodl_on_statistics(theory_, step.to_c(), accu.to_c()));
    }

    void main(Clingo::Control &ctl, Clingo::StringSpan files) override {
        handle_error(clingodl_register(theory_, ctl.to_c()));

        Clingo::AST::with_builder(ctl, [&](Clingo::AST::ProgramBuilder &builder) {
            Rewriter rewriter{theory_, builder.to_c()};
            rewriter.rewrite(files);
        });

        ctl.ground({{"base", {}}});
        if (!opt_cfg_.active) {
            ctl.solve(Clingo::SymbolicLiteralSpan{}, this, false, false).get();
        }
        else {
            Optimizer{opt_cfg_, *this, theory_}.solve(ctl);
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
        opt_cfg_.active = true;
        opt_cfg_.symbol = args[0];
        if (size > 1) {
            opt_cfg_.has_initial = true;
            opt_cfg_.initial = args[1].number();
        }
        return true;
    }

    bool parse_factor(char const *value) {
        std::stringstream strValue;
        strValue.imbue(std::locale::classic());
        strValue << value;
        double factor{0};
        strValue >> factor;
        if (factor < 1 || factor > std::numeric_limits<int>::max()) {
            throw std::overflow_error("minimize-factor out of bounds");
        }
        opt_cfg_.factor = factor;
        return strValue.rdbuf()->in_avail() == 0;
    }

    void register_options(Clingo::ClingoOptions &options) override {
        handle_error(clingodl_register_options(theory_, options.to_c()));
        char const * group = "Clingo.DL Options";
        options.add(group, "minimize-variable",
            "Minimize the given variable\n"
            "      <arg>   : <variable>[,<initial>]\n"
            "      <variable>: the variable to minimize\n"
            "      <initial> : upper bound for the variable",
            [this](char const *value) { return parse_bound(value); });
        options.add(group, "minimize-factor",
            "Factor to adjust optimization \n"
            "      <factor>   : {multiplication factor, 1=linear (default)}\n",
            [this](char const *value) { return parse_factor(value); });
    }

    void validate_options() override {
        handle_error(clingodl_validate_options(theory_));
    }

private:
    clingodl_theory_t *theory_;
    OptimizerConfig opt_cfg_;
};

int main(int argc, char *argv[]) {
    ClingoDLApp app;
    return Clingo::clingo_main(app, {argv + 1, static_cast<size_t>(argc - 1)});
}
