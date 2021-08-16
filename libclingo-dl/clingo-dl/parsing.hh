// {{{ MIT License
//
// // Copyright 2018 Roland Kaminski, Philipp Wanko, Max Ostrowski
//
// // Permission is hereby granted, free of charge, to any person obtaining a copy
// // of this software and associated documentation files (the "Software"), to
// // deal in the Software without restriction, including without limitation the
// // rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// // sell copies of the Software, and to permit persons to whom the Software is
// // furnished to do so, subject to the following conditions:
//
// // The above copyright notice and this permission notice shall be included in
// // all copies or substantial portions of the Software.
//
// // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// // IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// // FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// // AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// // LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// // FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// // IN THE SOFTWARE.
//
// // }}}

#ifndef CLINGODL_PARSING_HH
#define CLINGODL_PARSING_HH

#include <clingo.hh>

namespace ClingoDL {

//! Epsilon value depending on number type.
template <class N>
[[nodiscard]] N epsilon();

//! Vector of coefficients and variables.
template <class N>
using CoVarVec = std::vector<std::pair<N, int>>;

//! An edge in the difference logic graph.
template <typename N>
struct EdgeAtom {
    CoVarVec<N> lhs;
    char const *rel;
    N rhs;
    int literal;
    bool strict;
};

//! Throw a syntax error with the given message.
template <typename T=void>
inline T throw_syntax_error(char const *message = "Invalid Syntax") {
    throw std::runtime_error(message);
}

//! The a syntax error if the given condition is false.
inline void check_syntax(bool condition, char const *message = "Invalid Syntax") {
    if (!condition) {
        throw_syntax_error(message);
    }
}

//! Return true if the given theory term matches the given signature.
[[nodiscard]] bool match(Clingo::TheoryTerm const &term, char const *name, size_t arity);

//! Parse a theory atom for a difference constraint.
template <class N>
[[nodiscard]] EdgeAtom<N> parse(Clingo::TheoryAtom const &atom, std::function<int (Clingo::Symbol)> const &map_vert);

} // namespace ClingoDL

#endif
