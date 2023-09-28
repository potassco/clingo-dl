// {{{ MIT License

// Copyright Roland Kaminski, Philipp Wanko, and Max Ostrowski

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

#include <clingo-dl-app/app.hh>
#include <clingo-dl.h>
#include <clingo.hh>
#include <fstream>
#include <limits>
#include <sstream>

#ifdef CLINGODL_PROFILE
#include <gperftools/profiler.h>
#endif

namespace ClingoDL {

using Clingo::Detail::handle_error;

//! Application class to run clingo-dl.
class App : public Clingo::Application, private Clingo::SolveEventHandler {
  public:
    App() { handle_error(clingodl_create(&theory_)); }
    App(App const &) = default;
    App(App &&) = default;
    auto operator=(App const &) -> App & = default;
    auto operator=(App &&) -> App & = default;
    ~App() override { clingodl_destroy(theory_); }
    //! Set program name to clingo-dl.
    auto program_name() const noexcept -> char const * override { return "clingo-dl"; }
    //! Set the version.
    auto version() const noexcept -> char const * override { return CLINGODL_VERSION; }
    //! Pass models to the theory.
    auto on_model(Clingo::Model &model) -> bool override {
        handle_error(clingodl_on_model(theory_, model.to_c()));
        return true;
    }
    //! Pass statistics to the theory.
    void on_statistics(Clingo::UserStatistics step, Clingo::UserStatistics accu) override {
        handle_error(clingodl_on_statistics(theory_, step.to_c(), accu.to_c()));
    }
    //! Run main solving function.
    void main(Clingo::Control &ctl, Clingo::StringSpan files) override { // NOLINT
        handle_error(clingodl_register(theory_, ctl.to_c()));

        Clingo::AST::with_builder(ctl, [&](Clingo::AST::ProgramBuilder &builder) {
            Rewriter rewriter{theory_, builder.to_c()};
            rewriter.rewrite(ctl, files);
        });

        ctl.ground({{"base", {}}});
#ifdef CLINGODL_PROFILE
        ProfilerStart("clingodl.solve.prof");
#endif
        if (!opt_cfg_.active) {
            ctl.solve(Clingo::SymbolicLiteralSpan{}, this, false, false).get();
        } else {
            Optimizer{opt_cfg_, *this, theory_}.solve(ctl);
        }
#ifdef CLINGODL_PROFILE
        ProfilerStop();
#endif
    }
    //! Parse the variable to minimize and an optional initial bound.
    auto parse_bound(char const *value) -> bool {
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
    //! Parse factor to adjust optimization step length.
    auto parse_factor(char const *value) -> bool {
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
    //! Register options of the theory and optimization related options.
    void register_options(Clingo::ClingoOptions &options) override {
        handle_error(clingodl_register_options(theory_, options.to_c()));
        char const *group = "Clingo.DL Options";
        options.add(group, "minimize-variable",
                    "Minimize the given variable\n"
                    "      <arg>     : <variable>[,<initial>]\n"
                    "      <variable>: the variable to minimize\n"
                    "      <initial> : upper bound for the variable",
                    [this](char const *value) { return parse_bound(value); });
        options.add(
            group, "minimize-factor", "Factor to adjust minimization step size [1]",
            [this](char const *value) { return parse_factor(value); }, false, "<factor>");
    }
    //! Validate options of the theory.
    void validate_options() override { handle_error(clingodl_validate_options(theory_)); }

  private:
    clingodl_theory_t *theory_{nullptr}; //!< The underlying DL theory.
    OptimizerConfig opt_cfg_;            //!< The optimization configuration.
};

} // namespace ClingoDL

//! Run the clingo-dl application.
auto main(int argc, char *argv[]) -> int { // NOLINT(bugprone-exception-escape)
    ClingoDL::App app;
    return Clingo::clingo_main(app, {argv + 1, static_cast<size_t>(argc - 1)});
}
