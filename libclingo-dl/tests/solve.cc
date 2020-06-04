#include <clingo.hh>
#include <clingo-dl.h>
#include <clingo-dl/propagator.hh>
#include "catch.hpp"

class Handler : public Clingo::SolveEventHandler {
public:
    Handler(clingodl_theory_t *theory) : theory_{theory} { }
    void on_statistics(Clingo::UserStatistics step, Clingo::UserStatistics accu) override {
        clingodl_on_statistics(theory_, step.to_c(), accu.to_c());
    }

private:
    clingodl_theory_t *theory_;
};

class Rewriter {
public:
    Rewriter(clingodl_theory_t *theory, clingo_program_builder_t *builder)
    : theory_{theory}
    , builder_{builder} {
    }

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

using ResultVec = std::vector<std::pair<std::vector<std::pair<Clingo::Symbol, double>>,std::vector<Clingo::Symbol>>>;
ResultVec solve(clingodl_theory_t *theory, Clingo::Control &ctl) {
    Handler h{theory};
    using namespace Clingo;
    ResultVec result;
    for (auto &&m : ctl.solve(LiteralSpan{}, &h)) {
        result.emplace_back();
        auto &sol = result.back().first;
        auto &sol_bool = result.back().second;
        auto id = m.thread_id();
        size_t index;
        for (clingodl_assignment_begin(theory, id, &index); clingodl_assignment_next(theory, id, &index); ) {
            clingodl_value_t value;
            clingodl_assignment_get_value(theory, id, index, &value);
            if (value.type == clingodl_value_type_int) {
                sol.emplace_back(Symbol{clingodl_get_symbol(theory, index)}, value.int_number);
            }
            else if (value.type == clingodl_value_type_double) {
                sol.emplace_back(Symbol{clingodl_get_symbol(theory, index)}, value.double_number);
            }
            else {
                REQUIRE(false);
            }
        }
        std::sort(sol.begin(), sol.end());
        for (auto s : m.symbols()) {
            sol_bool.emplace_back(s);
        }
        std::sort(sol_bool.begin(), sol_bool.end());
    }
    std::sort(result.begin(), result.end());
    return result;
}

void parse_program(clingodl_theory_t *theory, Clingo::Control &ctl, const char *str) {
    ctl.with_builder([&](Clingo::ProgramBuilder &builder) {
        Rewriter rewriter{theory, builder.to_c()};
        clingo_parse_program(str, Rewriter::rewrite_, &rewriter, nullptr, nullptr, 0);
    });
}

TEST_CASE("solving", "[clingo]") {
    SECTION("with control") {
        using namespace Clingo;
        auto a = Id("a"),  b = Id("b"), c = Id("c");
        Control ctl{{"0"}};
        clingodl_theory_t *theory;
        REQUIRE(clingodl_create(&theory));
        SECTION("solve") {
            REQUIRE(clingodl_register(theory, ctl.to_c()));
            parse_program(theory, ctl,
                "#program base.\n"
                "1 { a; b } 1. &diff { a - b } <= 3.\n"
                "&diff { 0 - a } <= -5 :- a.\n"
                "&diff { 0 - b } <= -7 :- b.\n"
                "&show_assignment{}.\n");
            ctl.ground({{"base", {}}});
            REQUIRE(clingodl_prepare(theory, ctl.to_c()));
            auto result = solve(theory, ctl);
            REQUIRE(result == (ResultVec{{{{a, 0}, {b, 7}},{b}}, {{{a, 5}, {b, 2}},{a}}}));

            parse_program(theory, ctl,
                "#program ext.\n"
                "&diff { a - 0 } <= 4.\n");
            ctl.ground({{"ext", {}}});
            REQUIRE(clingodl_prepare(theory, ctl.to_c()));
            result = solve(theory, ctl);
            REQUIRE(result == (ResultVec{{{{a, 0}, {b, 7}},{b}}}));
        }
        SECTION("cc") {
            REQUIRE(clingodl_register(theory, ctl.to_c()));
            parse_program(theory, ctl,
                "#program base.\n"
                "&diff { 0 - a } <= -5.\n"
                "&diff { 0 - b } <= -10.\n"
                "&show_assignment{}.\n");
            ctl.ground({{"base", {}}});
            REQUIRE(clingodl_prepare(theory, ctl.to_c()));
            auto result = solve(theory, ctl);
            REQUIRE(result == (ResultVec{{{{a, 5}, {b, 10}},{}}}));
            REQUIRE(ctl.statistics()["user_step"]["DifferenceLogic"]["CCs"] == 2);

            parse_program(theory, ctl,
                "#program ext.\n"
                "&diff { b - a } <= 3.\n");
            ctl.ground({{"ext", {}}});
            REQUIRE(clingodl_prepare(theory, ctl.to_c()));
            result = solve(theory, ctl);
            REQUIRE(result == (ResultVec{{{{a, 7}, {b, 10}},{}}}));
            REQUIRE(ctl.statistics()["user_step"]["DifferenceLogic"]["CCs"] == 1);
        }

        SECTION("configure") {
            REQUIRE(clingodl_configure(theory, "propagate", "full"));
            REQUIRE(clingodl_register(theory, ctl.to_c()));

            parse_program(theory, ctl,
                "#program base.\n"
                "&diff { a - 0 } <= 0.\n"
                "a :- &diff { a - 0 } <=  0.\n"
                "b :- &diff { 0 - a } <= -1.\n"
                "c :- &diff { a } <= -1.\n"
                );
            ctl.ground({{"base", {}}});
            REQUIRE(clingodl_prepare(theory, ctl.to_c()));

            auto result = solve(theory, ctl);
            REQUIRE(result == (ResultVec{{{{a, -1}},{a, c}},{{{a, 0}},{a}}}));

        }
        SECTION("rdl") {
            REQUIRE(clingodl_configure(theory, "rdl", "yes"));
            REQUIRE(clingodl_register(theory, ctl.to_c()));

            parse_program(theory, ctl,
                "#program base.\n"
                "&diff { a } >= \"0.5\" * 3.\n"
                );
            ctl.ground({{"base", {}}});
            REQUIRE(clingodl_prepare(theory, ctl.to_c()));

            auto result = solve(theory, ctl);
            REQUIRE(result == (ResultVec{{{{a, 1.5}},{}}}));
        }

        SECTION("parse") {
            REQUIRE(clingodl_register(theory, ctl.to_c()));

            parse_program(theory, ctl,
                "#program base.\n"
                "&diff { p( 1 + 2 ) - q( 3 * 4 - 7 ) } <= 3 - \"9.0\".\n"
                );
            ctl.ground({{"base", {}}});
            REQUIRE(clingodl_prepare(theory, ctl.to_c()));

            auto result = solve(theory, ctl);
            auto p = Clingo::parse_term("p(3)"),  q = Clingo::parse_term("q(5)");
            REQUIRE(result == (ResultVec{{{{p, 0}, {q, 6}},{}}}));
        }
        SECTION("normalize") {
            REQUIRE(clingodl_register(theory, ctl.to_c()));

            parse_program(theory, ctl,
                "#program base.\n"
                "&diff { a } = b.\n"
                "&diff { 5 } >= 0.\n"
                "&diff { b } > c.\n"
                "&diff { c } >= d + 1.\n"
                "&diff { e } != f.\n"
                );
            ctl.ground({{"base", {}}});
            REQUIRE(clingodl_prepare(theory, ctl.to_c()));

            auto result = solve(theory, ctl);
            auto b = Id("b"),  c = Id("c"), d = Id("d"),  e = Id("e"), f = Id("f");
            REQUIRE(result == (ResultVec{{{{a, 2},{b, 2},{c, 1},{d, 0},{e, 0}, {f, 1}},{}},{{{a, 2},{b, 2},{c, 1},{d, 0},{e, 1}, {f, 0}},{}}}));
        }
        SECTION("empty constraints") {
            REQUIRE(clingodl_register(theory, ctl.to_c()));

            parse_program(theory, ctl,
                "#program base.\n"
                "a :- &diff { a - a } <= 5.\n"
                "{ b }.\n"
                "&diff { 0 } < -4 :- b.\n"
                );
            ctl.ground({{"base", {}}});
            REQUIRE(clingodl_prepare(theory, ctl.to_c()));

            auto result = solve(theory, ctl);
            REQUIRE(result == (ResultVec{{{},{a}}}));
            REQUIRE(ctl.statistics()["solving"]["solvers"]["choices"] == 0);
        }
        clingodl_destroy(theory);
    }
}
