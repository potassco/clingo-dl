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

#include <clingo-dl/optimizer.hh>

#include <clingo/app.hh>

#include <limits>
#include <sstream>

#ifdef CLINGODL_PROFILE
#include <gperftools/profiler.h>
#endif

namespace ClingoDL {

//! Application class to run clingo-dl.
class App : public Clingo::App, private Clingo::SolveEventHandler {
  public:
    App(Clingo::Library const &lib) : lib_{lib} {}
    App(App &&other) = delete;
    //! Set program name to clingo-dl.
    auto do_program_name() noexcept -> std::string_view override { return "clingo-dl"; }
    //! Set the version.
    auto do_version() noexcept -> std::string_view override { return CLINGODL_VERSION; }
    //! Pass models to the theory.
    auto do_model(Clingo::Model model) -> bool override {
        theory_.model(model);
        return true;
    }
    //! Pass statistics to the theory.
    void do_stats(Clingo::Stats step, Clingo::Stats accu) override { theory_.stats(step, accu); }
    //! Run main solving function.
    void do_main(Clingo::Control const &ctl, Clingo::StringSpan files) override { // NOLINT
        theory_.register_theory(ctl);
        theory_.rewrite(lib_, ctl, files);
        if (ctl.mode() == Clingo::ControlMode::solve) {
            ctl.ground();
#ifdef CLINGODL_PROFILE
            ProfilerStart("clingodl.solve.prof");
#endif
            if (!opt_cfg_.active) {
                theory_.prepare(ctl);
                std::ignore = ctl.solve({}, std::ref<Clingo::SolveEventHandler>(*this));
            } else {
                Optimizer{lib_, opt_cfg_, *this, theory_}.solve(ctl);
            }
#ifdef CLINGODL_PROFILE
            ProfilerStop();
#endif
        } else {
            ctl.main();
        }
    }
    //! Parse the variable to minimize and an optional initial bound.
    auto parse_bound(std::string_view value) -> bool {
        std::ostringstream oss;
        oss << "(" << value << ",)";
        auto term = Clingo::parse_term(lib_, oss.view());
        auto args = term.arguments();
        auto size = args.size();
        if (args.empty() || size > 2 || (size > 1 && args[1].type() != Clingo::SymbolType::number)) {
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
    auto parse_factor(std::string_view value) -> bool {
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
    void do_register_options(Clingo::Options options) override {
        using namespace std::string_view_literals;
        theory_.register_options(options);
        auto group = "Clingo.DL Options"sv;
        options.add(group, "minimize-variable",
                    "Minimize the given variable\n"
                    "      <arg>     : <variable>[,<initial>]\n"
                    "      <variable>: the variable to minimize\n"
                    "      <initial> : upper bound for the variable",
                    [this](std::string_view value) { return parse_bound(value); });
        options.add(
            group, "minimize-factor", "Factor to adjust minimization step size [1]",
            [this](std::string_view value) { return parse_factor(value); }, false, "<factor>");
    }
    //! Validate options of the theory.
    void do_validate_options() override { theory_.validate_options(); }

  private:
    Clingo::Library lib_;
    Clingo::Theory theory_{lib_, clingodl_create};
    OptimizerConfig opt_cfg_; //!< The optimization configuration.
};

} // namespace ClingoDL

//! Run the clingo-dl application.
auto main(int argc, char *argv[]) -> int { // NOLINT(bugprone-exception-escape)
    Clingo::Library lib;
    ClingoDL::App app{lib};
    auto args = std::vector<std::string_view>{argv + 1, argv + argc};
    return Clingo::main(lib, args, &app);
}
