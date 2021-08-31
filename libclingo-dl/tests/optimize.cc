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

#include <clingo.hh>
#include <clingo-dl.h>
#include <clingo-dl/propagator.hh>
#include <clingo-dl-app/app.hh>
#include "catch.hpp"

namespace {

struct SolveHandler : public Clingo::SolveEventHandler {
public:
    SolveHandler(clingodl_theory_t *theory)
    : theory_{theory} {
    }
    [[nodiscard]] Clingo::SymbolVector const &symbols() const {
        return symbols_;
    }
private:
    bool on_model(Clingo::Model &model) override {
        REQUIRE(clingodl_on_model(theory_, model.to_c()));
        symbols_ = model.symbols(Clingo::ShowType::Theory);
        std::sort(symbols_.begin(), symbols_.end());
        return false;
    }
    void on_statistics(Clingo::UserStatistics step, Clingo::UserStatistics accu) override {
        REQUIRE(clingodl_on_statistics(theory_, step.to_c(), accu.to_c()));
    }
    clingodl_theory_t *theory_;
    Clingo::SymbolVector symbols_;
};

void parse_program(clingodl_theory_t *theory, Clingo::Control &ctl, const char *str) {
    Clingo::AST::with_builder(ctl, [&](Clingo::AST::ProgramBuilder &builder) {
        ClingoDL::Rewriter rewriter{theory, builder.to_c()};
        rewriter.rewrite(str);
    });
}

Clingo::Symbol assign(Clingo::Symbol name, int value) {
    return Clingo::Function("dl", {name, Clingo::Number(value)});
}

Clingo::SymbolVector optimize(Clingo::Control &ctl, Clingo::Symbol bound, double factor, char const *prg) {
    clingodl_theory_t *theory{nullptr};
    REQUIRE(clingodl_create(&theory));
    REQUIRE(clingodl_register(theory, ctl.to_c()));
    parse_program(theory, ctl, prg);
    ctl.ground({{"base", {}}});
    REQUIRE(clingodl_prepare(theory, ctl.to_c()));
    SolveHandler handler = {theory};
    ClingoDL::OptimizerConfig cfg;
    cfg.symbol = bound;
    cfg.factor = factor;
    ClingoDL::Optimizer{cfg, handler, theory}.solve(ctl);
    REQUIRE(clingodl_destroy(theory));
    return handler.symbols();
}

} // namespace

TEST_CASE("optimize", "[clingo-dl]") {
    auto a = Clingo::Id("a");
    auto b = Clingo::Id("b");
    SECTION("sat") {
        for (auto factor : {1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0}) {
            Clingo::Control ctl{{"1"}};
            REQUIRE(optimize(
                ctl, b, factor,
                "&diff { a - 0 } >=  100.\n"
                "&diff { b - 0 } >= -100.\n"
                "&diff { a - b } <=  150.\n") == Clingo::SymbolVector{assign(a, 100), assign(b, -50)});
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
