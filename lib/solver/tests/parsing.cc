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

#include <clingo/control.hh>

#include <catch2/catch_test_macros.hpp>

#include <map>

namespace ClingoDL {

namespace {

//! Vector of strings for capturing results.
using V = std::vector<std::string>;

//! Rewrite the given statement and return the string representation of the
//! resulting statements.
auto rewrite(std::string_view prg) -> V {
    auto ret = V{};
    auto lib = Clingo::Library{};
    auto stm = Clingo::AST::parse(lib, prg);
    ClingoDL::rewrite(lib, stm, [&]<class T>(T &&node) { ret.emplace_back(node.to_string()); }, true);
    return ret;
}

//! Parse the theory atoms in the given program and return their string
//! representations.
template <class N> auto parse(std::string_view str) -> V {
    auto lib = Clingo::Library{};
    auto ctl = Clingo::Control{lib};
    auto prg = Clingo::AST::Program{lib};
    auto stm = Clingo::AST::parse(lib, str);
    ClingoDL::rewrite(lib, stm, [&]<class T>(T &&node) { prg.add(std::forward<T>(node)); }, true);
    ctl.join(prg);
    ctl.parse_string(THEORY);
    ctl.ground();

    auto vertex_map = std::map<Clingo::Symbol, vertex_t>{};
    auto vertices = std::vector<Clingo::Symbol>{};
    auto oss = std::ostringstream{};
    V ret;
    for (auto atom : ctl.base().theory()) {
        auto edge = ClingoDL::parse<N>(lib, atom, [&](Clingo::Symbol sym) {
            auto [it, ins] = vertex_map.try_emplace(sym, static_cast<vertex_t>(vertices.size()));
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
        oss << " " << relation_to_string(edge.rel) << " " << edge.rhs << " (" << (edge.strict ? "strict" : "non-strict")
            << ")";
        ret.emplace_back(oss.str());
    }
    return ret;
}

} // namespace

TEST_CASE("parsing rewrite", "[parsing][rewrite]") {
    REQUIRE(rewrite("&diff { a - b } < 0.") == V{"&__diff_h { (a - b) } < 0."});
    REQUIRE(rewrite("x :- &diff { a - b } < 0.") == V{"x :- &__diff_b { (a - b) } < 0."});
    REQUIRE(rewrite(" :- &diff { a - b } < 0.") == V{"&__diff_h { (a - b) } >= 0."});
    REQUIRE(rewrite(" :- not &diff { a - b } >= 0.") == V{"&__diff_h { (a - b) } >= 0."});
}

TEST_CASE("parsing parse int", "[parsing][parse][int]") {
    SECTION("strict / non-strict") {
        REQUIRE(parse<int>("&diff { a - b } < 0.") == V{"1*a + -1*b < 0 (non-strict)"});
        REQUIRE(parse<int>(":- &diff { a - b } < 0.") == V{"1*a + -1*b >= 0 (non-strict)"});
        REQUIRE(parse<int>("a :- &diff { a - b } < 0.") == V{"1*a + -1*b < 0 (strict)"});
    }
    SECTION("relations") {
        REQUIRE(parse<int>("&diff { a - b } < 0.") == V{"1*a + -1*b < 0 (non-strict)"});
        REQUIRE(parse<int>("&diff { a - b } <= 0.") == V{"1*a + -1*b <= 0 (non-strict)"});
        REQUIRE(parse<int>("&diff { a - b } > 0.") == V{"1*a + -1*b > 0 (non-strict)"});
        REQUIRE(parse<int>("&diff { a - b } >= 0.") == V{"1*a + -1*b >= 0 (non-strict)"});
        REQUIRE(parse<int>("&diff { a - b } = 0.") == V{"1*a + -1*b = 0 (non-strict)"});
        REQUIRE(parse<int>("&diff { a - b } != 0.") == V{"1*a + -1*b != 0 (non-strict)"});
    }
    SECTION("complex") {
        REQUIRE(parse<int>("&diff { 2 * (a - (x + b)) + 6*c } < 2*(3*c - x) + a - b.") ==
                V{"1*a + -1*b < 0 (non-strict)"});
    }
}

TEST_CASE("parsing parse double", "[parsing][parse][double]") {
    REQUIRE(parse<double>(R"p(&diff { a - b } < "10".)p") == V{"1*a + -1*b < 10 (non-strict)"});
    REQUIRE(parse<double>(R"p(&diff { a - b } < "10.5".)p") == V{"1*a + -1*b < 10.5 (non-strict)"});
    REQUIRE(parse<double>(R"p(&diff { "2.0" * (a - (x + b)) + "6"*c } < 2*("3.123"*c - x) + a - b.)p") ==
            V{"-0.246*c + 1*a + -1*b < 0 (non-strict)"});
}

} // namespace ClingoDL
