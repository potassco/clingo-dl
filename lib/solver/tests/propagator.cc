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

#include <clingo-dl/parsing.hh>
#include <clingo-dl/propagator.hh>

#include <clingo/ast.hh>
#include <clingo/control.hh>

#include <catch2/catch_test_macros.hpp>

namespace ClingoDL {

namespace {

using MV = std::vector<std::vector<std::string>>;

template <class N> class MCB : public Clingo::SolveEventHandler {
  public:
    MCB(DLPropagator<N> &prp, MV &models) : prp_{&prp}, models_{&models} { models_->clear(); }
    ~MCB() override { std::ranges::sort(*models_); }

  private:
    auto do_model(Clingo::Model model) -> bool override {
        prp_->extend_model(model);
        if (!proven && model.optimality_proven()) {
            models_->clear();
            proven = true;
        }
        models_->emplace_back();
        for (auto &sym : model.symbols(Clingo::ShowFlags::shown)) {
            models_->back().push_back(sym.to_string());
        }
        std::ranges::sort(models_->back());
        return true;
    }
    DLPropagator<N> *prp_;
    MV *models_;
    bool proven = false;
};

template <class N> auto solve(std::string_view str) -> MV {
    auto sts = Statistics{};
    auto cfg = PropagatorConfig{};
    auto lib = Clingo::Library{};
    auto ctl = Clingo::Control{lib, {"0"}};
    ctl.parse_string(THEORY);
    auto &prp = ctl.register_propagator(std::make_unique<DLPropagator<N>>(lib, sts, cfg));
    {
        auto prg = Clingo::AST::Program{lib};
        Clingo::AST::parse(lib, str, [&](Clingo::AST::Node const &stm) {
            rewrite(lib, std::move(stm), [&]<class T>(T &&stm) { prg.add(std::forward<T>(stm)); }, true);
        });
        ctl.join(prg);
    }
    ctl.ground();
    MV models;
    {
        auto mcb = MCB{prp, models};
        auto hnd = ctl.solve(mcb);
        std::ignore = hnd.get();
    }
    return models;
}

} // namespace

TEST_CASE("propagator", "[propagator]") {
    REQUIRE(solve<int>(R"(
        &diff { a - b } < 0 :- a. 
        &diff { b - a } < 0 :- b.
        1 {a; b}.
        )") == MV{
                   {"a", "dl(a,0)", "dl(b,1)"},
                   {"b", "dl(a,1)", "dl(b,0)"},
               });
    REQUIRE(solve<double>(R"(
        &diff { a - b } <= -"0.1" :- a. 
        &diff { b - a } <= -"0.1" :- b.
        1 {a; b}.
        )") == MV{
                   {"a", "dl(a,\"-0.000000\")", "dl(b,\"0.100000\")"},
                   {"b", "dl(a,\"0.100000\")", "dl(b,\"-0.000000\")"},
               });
}

} // namespace ClingoDL
