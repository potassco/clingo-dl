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

#ifndef CLINGODL_PARSING_HH
#define CLINGODL_PARSING_HH

#include <clingo-dl/theory.hh>
#include <clingo.hh>

namespace ClingoDL {

//! The difference logic theory necessary to parse difference constraints.
static constexpr char const *THEORY = R"(#theory dl {
term {
  + : 1, binary, left;
  - : 1, binary, left;
  * : 2, binary, left;
  / : 2, binary, left;
  - : 3, unary
};
&__diff_h/0 : term, {<=,>=,<,>,=,!=}, term, head;
&__diff_b/0 : term, {<=,>=,<,>,=,!=}, term, body
}.)";

//! Throw a syntax error with the given message.
template <typename T = void> inline auto throw_syntax_error(char const *message = "Invalid Syntax") -> T {
    throw std::runtime_error(message);
}

//! The a syntax error if the given condition is false.
inline void check_syntax(bool condition, char const *message = "Invalid Syntax") {
    if (!condition) {
        throw_syntax_error(message);
    }
}

//! Callback function to retrieve AST nodes.
using NodeCallback = std::function<void(Clingo::AST::Node &&ast)>;

//! Transform the given statement with dl constraints and pass it on to the
//! given callback.
//!
//! Optionally shifts constraints from rule bodies into heads of integrity
//! constraints if possible.
void transform(Clingo::AST::Node const &ast, NodeCallback const &cb, bool shift);

//! Return true if the given theory term matches the given signature.
[[nodiscard]] auto match(Clingo::TheoryTerm const &term, char const *name, size_t arity) -> bool;

//! Parse a theory atom for a difference constraint.
template <class N>
[[nodiscard]] auto parse(Clingo::TheoryAtom const &atom, std::function<int(Clingo::Symbol)> const &map_vert)
    -> EdgeAtom<N>;

} // namespace ClingoDL

#endif
