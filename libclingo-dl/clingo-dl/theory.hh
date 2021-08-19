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

#ifndef CLINGODL_THEORY_HH
#define CLINGODL_THEORY_HH

#include <clingo.hh>

namespace ClingoDL {

//! Type for vertices/variables in the theory.
using vertex_t = uint32_t;
//! Vector for vertex indices.
using VertexIndexVec = std::vector<vertex_t>;
//! Type for edge indices in the theory.
using edge_t = uint32_t;
//! Import literal_t from Clingo namespace.
using Clingo::literal_t;
//! Type for decision levels.
using level_t = uint32_t;

//! Vector of coefficients and variables.
template <class N>
using CoVarVec = std::vector<std::pair<N, vertex_t>>;

//! An edge in the difference logic graph.
template <typename N>
struct EdgeAtom {
    CoVarVec<N> lhs;
    char const *rel;
    N rhs;
    Clingo::literal_t literal;
    bool strict;
};

//! Epsilon value depending on number type.
template <class N>
[[nodiscard]] N epsilon();

} // namespace

#endif
