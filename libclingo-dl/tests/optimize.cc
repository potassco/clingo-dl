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

#include <clingo-dl-app/app.hh>
#include <clingo-dl.h>
#include <clingo-dl/propagator.hh>

#include <clingo/theory.hh>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

namespace ClingoDL {

namespace {

struct Fixture : public Clingo::SolveEventHandler {
    Fixture() = default;
    Fixture(Fixture &&other) = delete;
    //! Store the last model.
    auto do_model(Clingo::Model model) -> bool override {
        theory.model(model);
        symbols = model.symbols(Clingo::ShowFlags::theory);
        std::ranges::sort(symbols);
        return false;
    }
    //! Let's the theory add statistics.
    void do_stats(Clingo::Stats step, Clingo::Stats accu) override { theory.stats(step, accu); }

    //! Create symbols representing DL assignments.
    auto assign(Clingo::Symbol const &name, int value) -> Clingo::Symbol {
        return Clingo::Function(lib, "dl", {name, Clingo::Number(value)});
    }

    //! Run the optimization algorithm minimizing the given variable.
    auto optimize(Clingo::Control const &ctl, Clingo::Symbol const &bound, double factor, std::string_view prg)
        -> Clingo::SymbolVector {
        theory.register_theory(ctl);
        theory.rewrite(lib, ctl, prg);
        ctl.ground();
        theory.prepare(ctl);
        auto cfg = OptimizerConfig{};
        cfg.symbol = bound;
        cfg.factor = factor;
        Optimizer{lib, cfg, *this, theory}.solve(ctl);
        return symbols;
    }

    Clingo::Library lib;
    Clingo::Theory theory{lib, clingodl_create}; //!< The DL theory.
    Clingo::SymbolVector symbols;                //!< The symbols in the last model.
};

} // namespace

TEST_CASE_METHOD(Fixture, "optimize sat", "[clingo-dl]") {
    auto factor = GENERATE(1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0);
    auto a = Clingo::Function(lib, "a");
    auto b = Clingo::Function(lib, "b");
    auto ctl = Clingo::Control{lib, {"1"}};
    auto ret = Clingo::SymbolVector{assign(a, 100), assign(b, -50)}; // NOLINT
    REQUIRE(optimize(ctl, b, factor,
                     "&diff { a - 0 } >=  100.\n"
                     "&diff { b - 0 } >= -100.\n"
                     "&diff { a - b } <=  150.\n") == ret);
    REQUIRE(ctl.stats()["user_step"]["DifferenceLogic"].map().contains("Optimization"));
    REQUIRE(ctl.stats()["user_step"]["DifferenceLogic"]["Optimization"].value() == -50);
}

TEST_CASE_METHOD(Fixture, "optimize unsat", "[clingo-dl]") {
    auto ctl = Clingo::Control{lib, {"1"}};
    auto b = Clingo::Function(lib, "b");
    REQUIRE(optimize( // NOLINT
                ctl, b, 1.0,
                "&diff { a - 0 } >=  100.\n"
                "&diff { b - 0 } >= -100.\n"
                "&diff { a - b } <=  150.\n"
                "&diff { b - 0 } <  -50.\n") == Clingo::SymbolVector{});
    REQUIRE(!ctl.stats()["user_step"]["DifferenceLogic"].map().contains("Optimization"));
}

} // namespace ClingoDL
