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

#include <clingo/core.hh>

namespace ClingoDL {

//! Type for vertices/variables in the theory.
using vertex_t = uint32_t;
//! Vector for vertex indices.
using VertexIndexVec = std::vector<vertex_t>;
//! Type for edge indices in the theory.
using edge_t = uint32_t;
//! Import id_t from Clingo namespace.
using Clingo::ProgramId;
//! Import literal_t from Clingo namespace.
using Clingo::SolverLiteral;
//! Type for decision levels.
using level_t = uint32_t;
//! Type for array indices/sizes.
using index_t = uint32_t;
//! Type for program and solver literals.
using literal_t = Clingo::ProgramLiteral;
//! Type for ids.
using id_t = Clingo::ProgramId;

enum class Relation {
    less_than,
    less_equal,
    greater_than,
    greater_equal,
    equal_to,
    no_equal_to,
};

inline auto relation_from_string(std::string_view str) -> Relation {
    if (str == "<") {
        return Relation::less_than;
    }
    if (str == ">") {
        return Relation::greater_than;
    }
    if (str == "<=") {
        return Relation::less_equal;
    }
    if (str == ">=") {
        return Relation::greater_equal;
    }
    if (str == "=" || str == "==") {
        return Relation::equal_to;
    }
    if (str == "!=") {
        return Relation::no_equal_to;
    }
    throw std::logic_error{"invalid relation"};
}

inline auto relation_to_string(Relation rel) -> std::string_view {
    switch (rel) {
        case Relation::less_than: {
            return "<";
        }
        case Relation::less_equal: {
            return "<=";
        }
        case Relation::greater_than: {
            return ">";
        }
        case Relation::greater_equal: {
            return ">=";
        }
        case Relation::equal_to: {
            return "=";
        }
        case Relation::no_equal_to: {
            return "!=";
        }
    }
}

//! Vector of coefficients and variables.
template <class T> using CoVarVec = std::vector<std::pair<T, vertex_t>>;

//! An edge in the difference logic graph.
template <typename T> struct EdgeAtom {
    CoVarVec<T> lhs;               //!< The terms associated with the atom.
    Relation rel;                  //!< The comparision relation of the atom.
    T rhs;                         //!< The value on the right hand side.
    Clingo::SolverLiteral literal; //!< The literal associated with the atom.
    bool strict;                   //!< Whether the atom is strict.
};

//! Epsilon value depending on number type.
template <class T> [[nodiscard]] auto epsilon() -> T;

} // namespace ClingoDL

#endif
