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

#include "clingo-dl/theory.hh"

namespace ClingoDL {

namespace {

//! Epsilon value for floating point arithmetics.
constexpr double DOUBLE_EPSILON = 0.00001;

//! Epsilon value for integral numbers.
template <class T, typename std::enable_if<std::is_integral_v<T>, bool>::type = true>
[[nodiscard]] auto epsilon_() -> T {
    return 1;
}

//! Epsilon value for floating point numbers.
template <class T, typename std::enable_if<std::is_floating_point_v<T>, bool>::type = true>
[[nodiscard]] auto epsilon_() -> T {
    return DOUBLE_EPSILON;
}

} // namespace

template <class T> auto epsilon() -> T { return epsilon_<T>(); }

template int epsilon<int>();
template double epsilon<double>();

} // namespace ClingoDL
