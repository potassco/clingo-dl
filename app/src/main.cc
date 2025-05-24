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
#include <clingo/app.hh>
#include <limits>
#include <sstream>

#ifdef CLINGODL_PROFILE
#include <gperftools/profiler.h>
#endif

namespace ClingoDL {

using Clingo::Detail::handle_error;

//! Application class to run clingo-dl.
class App : public Clingo::App, private Clingo::SolveEventHandler {
  public:
    App() { handle_error(clingodl_create(c_cast(lib_), theory_)); }
    App(App &&other) = delete;
    ~App() override { theory_->destroy(theory_->self); }
    //! Set program name to clingo-dl.
    auto do_program_name() noexcept -> std::string_view override { return "clingo-dl"; }
    //! Set the version.
    auto do_version() noexcept -> std::string_view override { return CLINGODL_VERSION; }
    //! Pass models to the theory.
    auto do_model(Clingo::Model &model) -> bool override {
        handle_error(theory_->on_model(theory_->self, c_cast(model)));
        return true;
    }
    //! Pass statistics to the theory.
    void do_stats(Clingo::Stats step, Clingo::Stats accu) override {
        handle_error(theory_->on_stats(theory_->self, c_cast(accu)));
    }
    //! Run main solving function.
    void do_main(Clingo::Control const &ctl, Clingo::StringSpan files) override { // NOLINT
        handle_error(theory_->register_theory(theory_->self, c_cast(ctl)));
        auto prg = Clingo::AST::Program{lib_};
        rewrite(lib_, theory_, prg, files);
        ctl.join(prg);
        ctl.ground();
#ifdef CLINGODL_PROFILE
        ProfilerStart("clingodl.solve.prof");
#endif
        if (!opt_cfg_.active) {
            std::ignore = ctl.solve(*this).get();
        } else {
            Optimizer{lib_, opt_cfg_, *this, theory_}.solve(ctl);
        }
#ifdef CLINGODL_PROFILE
        ProfilerStop();
#endif
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
        handle_error(theory_->register_options(theory_->self, c_cast(options)));
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
    void do_validate_options() override { handle_error(theory_->validate_options(theory_->self)); }

  private:
    Clingo::Library lib_;
    clingo_theory_t *theory_{nullptr}; //!< The underlying DL theory.
    OptimizerConfig opt_cfg_;          //!< The optimization configuration.
};

} // namespace ClingoDL

//! Run the clingo-dl application.
auto main(int argc, char *argv[]) -> int { // NOLINT(bugprone-exception-escape)
    Clingo::Library lib;
    ClingoDL::App app;
    auto args = std::vector<std::string_view>{argv + 1, argv + argc - 1};
    return Clingo::main(lib, args, &app);
}
