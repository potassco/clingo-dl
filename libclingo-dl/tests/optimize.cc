// {{{ MIT License
//
// Copyright Roland Kaminski, Philipp Wanko, and Max Ostrowski
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
//
// }}}

#include <catch2/catch_test_macros.hpp>
#include <clingo-dl-app/app.hh>
#include <clingo-dl.h>
#include <clingo-dl/propagator.hh>
#include <clingo.hh>

namespace ClingoDL {

namespace {

//! Helper class to extract the last model while solving.
class SolveHandler : public Clingo::SolveEventHandler {
  public:
    SolveHandler(clingodl_theory_t *theory) : theory_{theory} {}
    //! Return the symbols in the last model found while solving.
    [[nodiscard]] auto symbols() const -> Clingo::SymbolVector const & { return symbols_; }

  private:
    //! Stores the last model.
    auto on_model(Clingo::Model &model) -> bool override {
        REQUIRE(clingodl_on_model(theory_, model.to_c()));
        symbols_ = model.symbols(Clingo::ShowType::Theory);
        std::sort(symbols_.begin(), symbols_.end());
        return false;
    }
    //! Let's the theory add statistics.
    void on_statistics(Clingo::UserStatistics step, Clingo::UserStatistics accu) override {
        REQUIRE(clingodl_on_statistics(theory_, step.to_c(), accu.to_c()));
    }

    clingodl_theory_t *theory_;    //!< The DL theory.
    Clingo::SymbolVector symbols_; //!< The symbols in the last model.
};

//! Parse and rewrite the given logic program.
void parse_program(clingodl_theory_t *theory, Clingo::Control &ctl, const char *str) {
    Clingo::AST::with_builder(ctl, [&](Clingo::AST::ProgramBuilder &builder) {
        Rewriter rewriter{theory, builder.to_c()};
        rewriter.rewrite(ctl, str);
    });
}

//! Create symbols representing DL assignments.
auto assign(Clingo::Symbol name, int value) -> Clingo::Symbol {
    return Clingo::Function("dl", {name, Clingo::Number(value)});
}

//! Run the optimization algorithm minimizing the given variable.
auto optimize(Clingo::Control &ctl, Clingo::Symbol bound, double factor, char const *prg) -> Clingo::SymbolVector {
    clingodl_theory_t *theory{nullptr};
    REQUIRE(clingodl_create(&theory));
    REQUIRE(clingodl_register(theory, ctl.to_c()));
    parse_program(theory, ctl, prg);
    ctl.ground({{"base", {}}});
    REQUIRE(clingodl_prepare(theory, ctl.to_c()));
    SolveHandler handler = {theory};
    OptimizerConfig cfg;
    cfg.symbol = bound;
    cfg.factor = factor;
    Optimizer{cfg, handler, theory}.solve(ctl);
    REQUIRE(clingodl_destroy(theory));
    return handler.symbols();
}

} // namespace

TEST_CASE("optimize", "[clingo-dl]") { // NOLINT
    auto a = Clingo::Id("a");
    auto b = Clingo::Id("b");
    SECTION("sat") {
        for (auto factor : {1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0}) {
            Clingo::Control ctl{{"1"}};
            Clingo::SymbolVector ret{assign(a, 100), assign(b, -50)}; // NOLINT
            REQUIRE(optimize(ctl, b, factor,
                             "&diff { a - 0 } >=  100.\n"
                             "&diff { b - 0 } >= -100.\n"
                             "&diff { a - b } <=  150.\n") == ret);
            REQUIRE(ctl.statistics()["user_step"]["DifferenceLogic"].has_subkey("Optimization"));
            REQUIRE(ctl.statistics()["user_step"]["DifferenceLogic"]["Optimization"].value() == -50);
        }
    }
    SECTION("unsat") {
        Clingo::Control ctl{{"1"}};
        REQUIRE(optimize( // NOLINT
                    ctl, b, 1.0,
                    "&diff { a - 0 } >=  100.\n"
                    "&diff { b - 0 } >= -100.\n"
                    "&diff { a - b } <=  150.\n"
                    "&diff { b - 0 } <  -50.\n") == Clingo::SymbolVector{});
        REQUIRE(!ctl.statistics()["user_step"]["DifferenceLogic"].has_subkey("Optimization"));
    }
}

} // namespace ClingoDL
