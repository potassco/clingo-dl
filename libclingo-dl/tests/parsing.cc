#include <clingo.hh>
#include <clingo-dl/parsing.hh>
#include <clingo-dl/propagator.hh>
#include "catch.hpp"
#include <sstream>
#include <map>

namespace ClingoDL {

namespace {

using V = std::vector<std::string>;

V parse(char const *prg) {
    Clingo::Control ctl;
    {
        Clingo::AST::ProgramBuilder builder{ctl};
        Clingo::AST::parse_string(prg, [&](Clingo::AST::Node const &ast) {
            ClingoDL::transform(ast, [&](Clingo::AST::Node &&trans) {
                builder.add(trans);
            }, true);
        });
    }
    ctl.add("base", {}, THEORY);
    ctl.ground({{"base", {}}});

    std::map<Clingo::Symbol, vertex_t> vertex_map;
    std::vector<Clingo::Symbol> vertices;
    std::ostringstream oss;
    V ret;
    for (auto &&atom : ctl.theory_atoms()) {
        auto edge = ClingoDL::parse<int>(atom, [&](Clingo::Symbol sym) {
            auto [it, ins] = vertex_map.try_emplace(sym, vertices.size());
            if (ins) {
                vertices.emplace_back(it->first);
            }
            return it->second;
        });
        bool plus = false;
        for (auto [co, var] : edge.lhs) {
            if (plus) {
                oss << " + ";
            }
            oss << co << "*" << vertices[var];
            plus = true;
        }
        oss << " " << edge.rel << " " << edge.rhs << " (" << (edge.strict ? "strict" : "non-strict") << ")";
        ret.emplace_back(oss.str());
    }
    return ret;
}

} // namespace

TEST_CASE("parsing", "[parsing]") {
    SECTION("strict / non-strict") {
        REQUIRE(parse("&diff { a - b } < 0.") == V{"1*a + -1*b < 0 (non-strict)"});
        REQUIRE(parse(":- &diff { a - b } < 0.") == V{"1*a + -1*b >= 0 (non-strict)"});
        REQUIRE(parse("a :- &diff { a - b } < 0.") == V{"1*a + -1*b < 0 (strict)"});
    }
    SECTION("relations") {
        REQUIRE(parse("&diff { a - b } < 0.") == V{"1*a + -1*b < 0 (non-strict)"});
        REQUIRE(parse("&diff { a - b } <= 0.") == V{"1*a + -1*b <= 0 (non-strict)"});
        REQUIRE(parse("&diff { a - b } > 0.") == V{"1*a + -1*b > 0 (non-strict)"});
        REQUIRE(parse("&diff { a - b } >= 0.") == V{"1*a + -1*b >= 0 (non-strict)"});
        REQUIRE(parse("&diff { a - b } = 0.") == V{"1*a + -1*b = 0 (non-strict)"});
        REQUIRE(parse("&diff { a - b } != 0.") == V{"1*a + -1*b != 0 (non-strict)"});
    }
    SECTION("complex") {
        REQUIRE(parse("&diff { 2 * (a - (x + b)) + 6*c } < 2*(3*c - x) + a - b.") == V{"1*a + -1*b < 0 (non-strict)"});
    }
}

} // namespace ClingoDL
