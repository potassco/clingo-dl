#include <clingo.hh>
#include <clingo-dl.h>
#include "catch.hpp"

TEST_CASE("solving", "[clingo]") {
    SECTION("with control") {
        using namespace Clingo;
        Control ctl{{"0"}};
        clingodl_propagator_t *prop;
        REQUIRE(clingodl_create_propagator(&prop));
        REQUIRE(clingodl_register_propagator(prop, ctl.to_c()));
        SECTION("solve") {
            ctl.add("base", {},
                "1 { a; b } 1. &diff { a - b } <= 3.\n"
                "&diff { 0 - a } <= -5 :- a.\n"
                "&diff { 0 - b } <= -7 :- b.\n"
                "&show_assignment{}.\n");
            ctl.ground({{"base", {}}});
            using ResultVec = std::vector<std::vector<std::pair<Clingo::Symbol, int>>>;
            ResultVec result;
            for (auto &&m : ctl.solve()) {
                result.emplace_back();
                auto &sol = result.back();
                auto id = m.thread_id();
                size_t index;
                for (clingodl_assignment_begin(prop, id, &index); clingodl_assignment_next(prop, id, &index); ) {
                    clingodl_value_t value;
                    clingodl_assignment_get_value(prop, id, index, &value);
                    REQUIRE(value.type == clingodl_value_type_int);
                    sol.emplace_back(Symbol{clingodl_get_symbol(prop, index)}, value.int_number);
                }
                std::sort(sol.begin(), sol.end());
            }
            std::sort(result.begin(), result.end());
            auto a = Id("a"),  b = Id("b");
            REQUIRE(result == (ResultVec{{{a, 0}, {b, 7}}, {{a, 5}, {b, 2}}}));
        }
    }
}
