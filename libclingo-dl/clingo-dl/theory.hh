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
//! Import id_t from Clingo namespace.
using Clingo::id_t;
//! Import literal_t from Clingo namespace.
using Clingo::literal_t;
//! Type for decision levels.
using level_t = uint32_t;
//! Type for array indices/sizes.
using index_t = uint32_t;

//! Vector of coefficients and variables.
template <class T> using CoVarVec = std::vector<std::pair<T, vertex_t>>;

//! An edge in the difference logic graph.
template <typename T> struct EdgeAtom {
    CoVarVec<T> lhs;           //!< The terms associated with the atom.
    char const *rel;           //!< The comparision relation of the atom.
    T rhs;                     //!< The value on the right hand side.
    Clingo::literal_t literal; //!< The literal associated with the atom.
    bool strict;               //!< Whether the atom is strict.
};

//! Epsilon value depending on number type.
template <class T> [[nodiscard]] auto epsilon() -> T;

} // namespace ClingoDL

#endif
