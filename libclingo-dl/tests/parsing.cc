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
#include <clingo-dl/parsing.hh>
#include <clingo.hh>
#include <map>
#include <sstream>

namespace ClingoDL {

namespace {

//! Vector of strings for capturing the representation of theory atoms.
using V = std::vector<std::string>;

//! Parse the theory atoms in the given program and return their string
//! representations.
auto parse(char const *prg) -> V {
    Clingo::Control ctl;
    {
        Clingo::AST::ProgramBuilder builder{ctl};
        Clingo::AST::parse_string(prg, [&](Clingo::AST::Node const &ast) {
            transform(
                ast, [&](Clingo::AST::Node &&trans) { builder.add(trans); }, true);
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

TEST_CASE("parsing", "[parsing]") { // NOLINT
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
