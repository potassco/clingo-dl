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

#ifndef CLINGODL_UTIL_HH
#define CLINGODL_UTIL_HH

#include <chrono>
#include <ostream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>

namespace std {

//! Specialization to compute hashes for pairs.
template<>
struct hash<std::pair<int, int>> {
    size_t operator()(std::pair<int, int> const &p) const;
};

} // namespace std

namespace ClingoDL {

//! A cast for integrals that causes an assertion if the target type cannot
//! represent the number.
template <class T, class S>
inline T numeric_cast(S s);

//! Helper to print unordered maps.
template <class K, class V>
std::ostream &operator<<(std::ostream &out, std::unordered_map<K, V> const &map);
//! Helper to print vectors.
template <class T>
std::ostream &operator<<(std::ostream &out, std::vector<T> const &vec);
//! Helper to print pairs.
template <class K, class V>
std::ostream &operator<<(std::ostream &out, std::pair<K, V> const &pair);

//! Helper to ensure that a resizable container can hold the given index.
template <class C>
void ensure_index(C &c, size_t index);

//! Duration used by the Timer class.
using Duration = std::chrono::duration<double>;

//! Simple timer class to measure durations.
class Timer {
private:
    //! Data type for time points.
    using TimePoint = std::chrono::time_point<std::chrono::steady_clock>;

public:
    //! The constructor expects a reference to a duration which is incremented
    //! by the lifetime of the timer.
    Timer(Duration &elapsed);
    Timer(Timer const &other) = delete;
    Timer(Timer &&other) = delete;
    Timer &operator=(Timer const &other) = delete;
    Timer &operator=(Timer &&other) = delete;
    void stop();
    ~Timer();

private:
    Duration &elapsed_;   //!< The duration to increment.
    TimePoint start_;     //!< The time when the timer was started.
    bool stopped_{false}; //!< Whether the timer has been stopped.
};

//! Heap class designed for DL propagation with N children per node.
//!
//! The heap stores indices and each member function takes a struct as argument
//! to get its current cost and set it's offset/index in the heap.
template <int N>
class Heap {
public:
    //! We use 32bit integers for indexing.
    using index_type = uint32_t;
    //! We use 32bit integers for sizes.
    using size_type = uint32_t;
    //! Add an element to the heap.
    template <class M>
    void push(M &m, index_type item);
    //! Remove the element with the lowest cost.
    template <class M>
    index_type pop(M &m);
    //! Update an element after it's cost has been decreased.
    template <class M>
    void decrease(M &m, index_type item);
    //! Update an element after it's cost has been increased.
    template <class M>
    void increase(M &m, index_type item);
    //! Get the number of elements in the heap.
    size_type size();
    //! Test if the heap is empty.
    bool empty();
    //! Remove all elements from the heap.
    void clear();

private:
    //! Swap to elements in the heap.
    template <class M>
    void swap_(M &m, index_type i, index_type j);
    //! Get the parent index of an element.
    index_type parent_(index_type offset);
    //! Get the index of the left most child of an element in the heap.
    index_type children_(index_type offset);
    //! Compare two elements by costs preferring relevant elements if the costs
    //! are equal.
    template <class M>
    bool less_(M &m, index_type a, index_type b);

    //! The vector storing the elements.
    std::vector<index_type> heap_;
};

// Some of the functions below could also be implemented using (much faster)
// compiler specific built-ins. For more information check the following links:
// - https://gcc.gnu.org/onlinedocs/gcc/Integer-Overflow-Builtins.html
// - https://wiki.sei.cmu.edu/confluence/display/c/INT32-C.+Ensure+that+operations+on+signed+integers+do+not+result+in+overflow

//! Safely add a and b throwing an exception in case of overflow/underflow.
template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int> = 0>
Int safe_add(Int a, Int b);

//! Add floating point numbers without specific overflow checking.
template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int> = 0>
Float safe_add(Float a, Float b);

//! Safely subtract a and b throwing an exception in case of
//! overflow/underflow.
template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int> = 0>
Int safe_sub(Int a, Int b);

//! Subtract floating point numbers without specific overflow checking.
template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int> = 0>
Float safe_sub(Float a, Float b);

//! Safely multiply a and b throwing an exception in case of
//! overflow/underflow.
template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int> = 0>
Int safe_mul(Int a, Int b);

//! Multiply floating point numbers without specific overflow checking.
template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int> = 0>
Float safe_mul(Float a, Float b);

//! Safely divide a and b throwing an exception in case of overflow/underflow.
template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int> = 0>
Int safe_div(Int a, Int b);

//! Divide floating point numbers without specific overflow checking.
template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int> = 0>
Float safe_div(Float a, Float b);

//! Safely calculate the modulo of a and b throwing an exception in case of
//! overflow/underflow.
template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int> = 0>
Int safe_mod(Int a, Int b);

//! Compute the modulo for floating point numbers without specific overflow checking.
template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int> = 0>
Float safe_mod(Float a, Float b);

//! Safely invert a throwing an exception in case of an underflow.
template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int> = 0>
Int safe_inv(Int a);

//! Invert a floating point number without specific overflow checking.
template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int> = 0>
Float safe_inv(Float a);

//! Safely exponentiate a and b throwing an exception in case of
//! overflow/underflow.
template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int> = 0>
Int safe_pow(Int a, Int b);

//! Exponentiate floating point numbers without specific overflow checking.
template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int> = 0>
Float safe_pow(Float a, Float b);

//! Restores quoted characters in the given string.
std::string unquote(char const* str);

} // namespace ClingoDL

#include <clingo-dl/impl/util.hh>

#endif // CLINGODL_UTIL_HH
