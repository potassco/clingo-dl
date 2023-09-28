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

#define CLINGODL_UTIL_IMPL
#include <clingo-dl/util.hh>

namespace std {

inline auto hash<std::pair<int, int>>::operator()(std::pair<int, int> const &p) const -> size_t {
    return static_cast<size_t>(p.first) + (static_cast<size_t>(p.second) >> 32); // NOLINT
}

} // namespace std

namespace ClingoDL {

namespace Detail {

//! Check if casting a signed to a signed integer preserves the number.
template <int X> using int_type = std::integral_constant<int, X>;
template <class T, class S> inline void nc_check(S s, int_type<0> t) { // same sign
    static_cast<void>(s);
    static_cast<void>(t);
    assert((std::is_same_v<T, S>) || (s >= std::numeric_limits<T>::min() && s <= std::numeric_limits<T>::max()));
}

//! Check if casting a signed to an unsigned integer preserves the number.
template <class T, class S> inline void nc_check(S s, int_type<-1> t) { // Signed -> Unsigned
    static_cast<void>(s);
    static_cast<void>(t);
    assert(s >= 0 && static_cast<S>(static_cast<T>(s)) == s);
}

//! Check if casting an unsigned to a signed integer preserves the number.
template <class T, class S> inline void nc_check(S s, int_type<1> t) { // Unsigned -> Signed
    static_cast<void>(s);
    static_cast<void>(t);
    assert(!(s > static_cast<std::make_unsigned_t<T>>(std::numeric_limits<T>::max())));
}

} // namespace Detail

template <class T, class S> inline auto numeric_cast(S s) -> T {
    constexpr int sv = int(std::numeric_limits<T>::is_signed) - int(std::numeric_limits<S>::is_signed);
    Detail::nc_check<T>(s, Detail::int_type<sv>());
    return static_cast<T>(s);
}

template <class T> auto operator<<(std::ostream &out, std::vector<T> const &vec) -> std::ostream & {
    out << "{";
    for (auto &x : vec) {
        out << " " << x;
    }
    out << " }";
    return out;
}

template <class K, class V> auto operator<<(std::ostream &out, std::unordered_map<K, V> const &map) -> std::ostream & {
    using T = std::pair<K, V>;
    std::vector<T> vec;
    vec.assign(map.begin(), map.end());
    std::sort(vec.begin(), vec.end(), [](T const &a, T const &b) { return a.first < b.first; });
    out << vec;
    return out;
}

template <class K, class V> auto operator<<(std::ostream &out, std::pair<K, V> const &pair) -> std::ostream & {
    out << "( " << pair.first << " " << pair.second << " )";
    return out;
}

template <class C> void ensure_index(C &c, size_t index) {
    if (index >= c.size()) {
        c.resize(index + 1);
    }
}

inline Timer::Timer(Duration &elapsed) : elapsed_(elapsed), start_(std::chrono::steady_clock::now()) {}

inline void Timer::stop() {
    if (!stopped_) {
        elapsed_ += std::chrono::steady_clock::now() - start_;
        stopped_ = true;
    }
}

inline Timer::~Timer() { stop(); }

template <int N> template <class M> void Heap<N>::push(M &m, index_type item) {
    m.offset(item) = size();
    heap_.push_back(item);
    decrease(m, item);
}

template <int N> template <class M> auto Heap<N>::pop(M &m) -> index_type {
    assert(!heap_.empty());
    auto ret = heap_[0];
    if (size() > 1) {
        heap_[0] = heap_.back();
        m.offset(heap_[0]) = 0;
        heap_.pop_back();
        increase(m, heap_[0]);
    } else {
        heap_.pop_back();
    }
    return ret;
}

template <int N> template <class M> void Heap<N>::decrease(M &m, index_type item) {
    auto i = m.offset(item);
    while (i > 0) {
        index_type p = parent_(i);
        if (less_(m, i, p)) {
            swap_(m, i, p);
            i = p;
        } else {
            break;
        }
    }
}

template <int N> template <class M> void Heap<N>::increase(M &m, index_type item) {
    auto i = m.offset(item);
    for (index_type p = i, j = children_(p), s = size(); j < s; j = children_(p)) {
        index_type min = j;
        for (index_type k = j + 1; k < j + N; ++k) {
            if (k < s && less_(m, k, min)) {
                min = k;
            }
        }
        if (less_(m, min, p)) {
            swap_(m, min, p);
            p = min;
        } else {
            return;
        }
    }
}

template <int N> auto Heap<N>::size() -> size_type { return numeric_cast<size_type>(heap_.size()); }

template <int N> auto Heap<N>::empty() -> bool { return heap_.empty(); }

template <int N> void Heap<N>::clear() { heap_.clear(); }

template <int N> template <class M> void Heap<N>::swap_(M &m, index_type i, index_type j) {
    m.offset(heap_[j]) = i;
    m.offset(heap_[i]) = j;
    std::swap(heap_[i], heap_[j]);
}

template <int N> auto Heap<N>::parent_(index_type offset) -> index_type { return (offset - 1) / N; }

template <int N> auto Heap<N>::children_(index_type offset) -> index_type { return N * offset + 1; }

template <int N> template <class M> auto Heap<N>::less_(M &m, index_type a, index_type b) -> bool {
    return m.less(heap_[a], heap_[b]);
}

template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int>> auto safe_add(Int a, Int b) -> Int {
    if (b > 0) {
        if (a > std::numeric_limits<Int>::max() - b) {
            throw std::overflow_error("integer overflow");
        }
    } else if (b < 0) {
        if (a < std::numeric_limits<Int>::min() - b) {
            throw std::underflow_error("integer underflow");
        }
    }
    return a + b;
}

template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int>>
auto safe_add(Float a, Float b) -> Float {
    return a + b;
}

template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int>> auto safe_sub(Int a, Int b) -> Int {
    if (b > 0) {
        if (a < std::numeric_limits<Int>::min() + b) {
            throw std::underflow_error("integer underflow");
        }
    } else if (b < 0) {
        if (a > std::numeric_limits<Int>::max() + b) {
            throw std::overflow_error("integer overflow");
        }
    }
    return a - b;
}

template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int>>
auto safe_sub(Float a, Float b) -> Float {
    return a - b;
}

template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int>> auto safe_mul(Int a, Int b) -> Int {
    if (a > 0) {
        if (b > 0) {
            if (a > (std::numeric_limits<Int>::max() / b)) {
                throw std::overflow_error("integer overflow");
            }
        } else if (b < (std::numeric_limits<Int>::min() / a)) {
            throw std::underflow_error("integer underflow");
        }
    } else {
        if (b > 0) {
            if (a < (std::numeric_limits<Int>::min() / b)) {
                throw std::underflow_error("integer underflow");
            }
        } else if (a < 0 && b < (std::numeric_limits<Int>::max() / a)) {
            throw std::overflow_error("integer overflow");
        }
    }
    return a * b;
}

template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int>>
auto safe_mul(Float a, Float b) -> Float {
    return a * b;
}

template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int>> auto safe_div(Int a, Int b) -> Int {
    if (a == std::numeric_limits<Int>::min() && b == -1) {
        throw std::overflow_error("integer overflow");
    }
    if (b == 0) {
        if (a < 0) {
            throw std::underflow_error("integer underflow");
        }
        throw std::overflow_error("integer overflow");
    }
    return a / b;
}

template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int>>
auto safe_div(Float a, Float b) -> Float {
    return a / b;
}

template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int>> auto safe_mod(Int a, Int b) -> Int {
    if (a == std::numeric_limits<Int>::min() && b == -1) {
        throw std::overflow_error("integer overflow");
    }
    if (b == 0) {
        if (a < 0) {
            throw std::underflow_error("integer underflow");
        }
        throw std::overflow_error("integer overflow");
    }
    return a % b;
}

template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int>>
auto safe_mod(Float a, Float b) -> Float {
    return fmod(a, b);
}

template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int>> auto safe_inv(Int a) -> Int {
    if (a == std::numeric_limits<Int>::min()) {
        throw std::overflow_error("integer overflow");
    }
    return -a;
}

template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int>> auto safe_inv(Float a) -> Float {
    return -a;
}

template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int>> auto safe_pow(Int a, Int b) -> Int {
    if (a == 0) {
        throw std::overflow_error("integer overflow");
    }
    auto ret = std::pow(static_cast<double>(a), b);
    if (ret > std::numeric_limits<Int>::max()) {
        throw std::overflow_error("integer overflow");
    }
    if (ret < std::numeric_limits<Int>::min()) {
        throw std::underflow_error("integer underflow");
    }
    return static_cast<Int>(ret);
}

template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int>>
auto safe_pow(Float a, Float b) -> Float {
    return std::pow(a, b);
}

inline auto unquote(char const *str) -> std::string {
    std::string res;
    bool slash = false;
    for (char const *it = *str == '"' ? str + 1 : str; *it != '\0'; ++it) { // NOLINT
        if (slash) {
            switch (*it) {
                case 'n': {
                    res.push_back('\n');
                    break;
                }
                case '\\': {
                    res.push_back('\\');
                    break;
                }
                case '"': {
                    res.push_back('"');
                    break;
                }
                default: {
                    assert(false);
                    break;
                }
            }
            slash = false;
        } else if (*it == '"' && *(it + 1) == '\0') {
            break;
        } // NOLINT
        else if (*it == '\\') {
            slash = true;
        } else {
            res.push_back(*it);
        }
    }
    return res;
}

} // namespace ClingoDL
