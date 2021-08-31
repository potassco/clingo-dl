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

template<>
struct hash<std::pair<int, int>> {
    size_t operator()(std::pair<int, int> const &p) const {
        return static_cast<size_t>(p.first) + (static_cast<size_t>(p.second) >> 32); // NOLINT
    }
};

} // namespace std

namespace Detail {

template <int X>
using int_type = std::integral_constant<int, X>;
template <class T, class S>
inline void nc_check(S s, int_type<0> t) { // same sign
    static_cast<void>(s);
    static_cast<void>(t);
    assert((std::is_same<T, S>::value) || (s >= std::numeric_limits<T>::min() && s <= std::numeric_limits<T>::max()));
}
template <class T, class S>
inline void nc_check(S s, int_type<-1> t) { // Signed -> Unsigned
    static_cast<void>(s);
    static_cast<void>(t);
    assert(s >= 0 && static_cast<S>(static_cast<T>(s)) == s);
}
template <class T, class S>
inline void nc_check(S s, int_type<1> t) { // Unsigned -> Signed
    static_cast<void>(s);
    static_cast<void>(t);
    assert(!(s > static_cast<std::make_unsigned_t<T>>(std::numeric_limits<T>::max())));
}

} // namespace Detail

template <class T, class S>
inline T numeric_cast(S s) {
    constexpr int sv = int(std::numeric_limits<T>::is_signed) - int(std::numeric_limits<S>::is_signed);
    ::Detail::nc_check<T>(s, ::Detail::int_type<sv>());
    return static_cast<T>(s);
}

template <class K, class V>
std::ostream &operator<<(std::ostream &out, std::unordered_map<K, V> const &map);
template <class T>
std::ostream &operator<<(std::ostream &out, std::vector<T> const &vec);
template <class K, class V>
std::ostream &operator<<(std::ostream &out, std::pair<K, V> const &pair);

template <class T>
std::ostream &operator<<(std::ostream &out, std::vector<T> const &vec) {
    out << "{";
    for (auto &x : vec) {
        out << " " << x;
    }
    out << " }";
    return out;
}

template <class K, class V>
std::ostream &operator<<(std::ostream &out, std::unordered_map<K, V> const &map) {
    using T = std::pair<K, V>;
    std::vector<T> vec;
    vec.assign(map.begin(), map.end());
    std::sort(vec.begin(), vec.end(), [](T const &a, T const &b) { return a.first < b.first; });
    out << vec;
    return out;
}

template <class K, class V>
std::ostream &operator<<(std::ostream &out, std::pair<K, V> const &pair) {
    out << "( " << pair.first << " " << pair.second << " )";
    return out;
}

template <class C>
void ensure_index(C &c, size_t index) {
    if (index >= c.size()) {
        c.resize(index + 1);
    }
}

using Duration = std::chrono::duration<double>;

class Timer {
public:
    Timer(Duration &elapsed)
        : elapsed_(elapsed)
        , start_(std::chrono::steady_clock::now()) {}
    Timer(Timer const &other) = delete;
    Timer(Timer &&other) = delete;
    Timer &operator=(Timer const &other) = delete;
    Timer &operator=(Timer &&other) = delete;
    ~Timer() { elapsed_ += std::chrono::steady_clock::now() - start_; }

private:
    Duration &elapsed_;
    std::chrono::time_point<std::chrono::steady_clock> start_;
};

template <int N>
class Heap {
public:
    using index_type = uint32_t;
    using size_type = uint32_t;
    template <class M>
    void push(M &m, index_type item) {
        auto i = m.offset(item) = size();
        heap_.push_back(item);
        decrease(m, i);
    }
    template <class M>
    index_type pop(M &m) {
        assert(!heap_.empty());
        auto ret = heap_[0];
        if (size() > 1) {
            heap_[0] = heap_.back();
            m.offset(heap_[0]) = 0;
            heap_.pop_back();
            increase(m, 0);
        }
        else {
            heap_.pop_back();
        }
        return ret;
    }

    template <class M>
    void decrease(M &m, index_type i) {
        while (i > 0) {
            index_type p = parent_(i);
            if (m.cost(heap_[p]) > m.cost(heap_[i])) {
                swap_(m, i, p);
                i = p;
            }
            else {
                break;
            }
        }
    }
    template <class M>
    void increase(M &m, index_type i) {
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
            }
            else {
                return;
            }
        }
    }
    size_type size() { return numeric_cast<size_type>(heap_.size()); }
    bool empty() { return heap_.empty(); }
    void clear() { heap_.clear(); }

private:
    template <class M>
    void swap_(M &m, index_type i, index_type j) {
        m.offset(heap_[j]) = i;
        m.offset(heap_[i]) = j;
        std::swap(heap_[i], heap_[j]);
    }
    index_type parent_(index_type offset) { return (offset - 1) / N; }
    index_type children_(index_type offset) { return N * offset + 1; }
    template <class M>
    bool less_(M &m, index_type a, index_type b) {
        a = heap_[a], b = heap_[b];
        auto ca = m.cost(a);
        auto cb = m.cost(b);
        return ca < cb || (ca == cb && m.relevant(a) < m.relevant(b));
    }

    std::vector<index_type> heap_;
};


// Some of the functions below could also be implemented using (much faster)
// compiler specific built-ins. For more information check the following links:
// - https://gcc.gnu.org/onlinedocs/gcc/Integer-Overflow-Builtins.html
// - https://wiki.sei.cmu.edu/confluence/display/c/INT32-C.+Ensure+that+operations+on+signed+integers+do+not+result+in+overflow

//! Safely add a and b throwing an exception in case of overflow/underflow.
template <typename Int, typename std::enable_if<std::is_integral<Int>::value, int>::type = 0>
Int safe_add(Int a, Int b) {
    if (b > 0) {
        if (a > std::numeric_limits<Int>::max() - b) {
            throw std::overflow_error("integer overflow");
        }
    }
    else if (b < 0) {
        if (a < std::numeric_limits<Int>::min() - b) {
            throw std::underflow_error("integer underflow");
        }
    }
    return a + b;
}

template <typename Float, typename std::enable_if<std::is_floating_point<Float>::value, int>::type = 0>
Float safe_add(Float a, Float b) {
	return a + b;
}

//! Safely subtract a and b throwing an exception in case of
//! overflow/underflow.
template <typename Int, typename std::enable_if<std::is_integral<Int>::value, int>::type = 0>
Int safe_sub(Int a, Int b) {
    if (b > 0) {
        if (a < std::numeric_limits<Int>::min() + b) {
            throw std::underflow_error("integer underflow");
        }
    }
    else if (b < 0) {
        if (a > std::numeric_limits<Int>::max() + b) {
            throw std::overflow_error("integer overflow");
        }
    }
    return a - b;
}

template <typename Float, typename std::enable_if<std::is_floating_point<Float>::value, int>::type = 0>
Float safe_sub(Float a, Float b) {
	return a - b;
}

//! Safely multiply a and b throwing an exception in case of
//! overflow/underflow.
template <typename Int, typename std::enable_if<std::is_integral<Int>::value, int>::type = 0>
Int safe_mul(Int a, Int b) {
    if (a > 0) {
        if (b > 0) {
            if (a > (std::numeric_limits<Int>::max() / b)) {
                throw std::overflow_error("integer overflow");
            }
        }
        else if (b < (std::numeric_limits<Int>::min() / a)) {
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

template <typename Float, typename std::enable_if<std::is_floating_point<Float>::value, int>::type = 0>
Float safe_mul(Float a, Float b) {
	return a * b;
}

//! Safely divide a and b throwing an exception in case of overflow/underflow.
template <typename Int, typename std::enable_if<std::is_integral<Int>::value, int>::type = 0>
Int safe_div(Int a, Int b) {
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

template <typename Float, typename std::enable_if<std::is_floating_point<Float>::value, int>::type = 0>
Float safe_div(Float a, Float b) {
	return a / b;
}


//! Safely calculate the modulo of a and b throwing an exception in case of
//! overflow/underflow.
template <typename Int, typename std::enable_if<std::is_integral<Int>::value, int>::type = 0>
Int safe_mod(Int a, Int b) {
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

template <typename Float, typename std::enable_if<std::is_floating_point<Float>::value, int>::type = 0>
Float safe_mod(Float a, Float b) {
	return fmod(a,b);
}

//! Safely invert a throwing an exception in case of an underflow.
template <typename Int, typename std::enable_if<std::is_integral<Int>::value, int>::type = 0>
Int safe_inv(Int a) {
    if (a == std::numeric_limits<Int>::min()) {
        throw std::overflow_error("integer overflow");
    }
    return -a;
}

template <typename Float, typename std::enable_if<std::is_floating_point<Float>::value, int>::type = 0>
Float safe_inv(Float a) {
	return -a;
}

template <typename Int, typename std::enable_if<std::is_integral<Int>::value, int>::type = 0>
Int safe_pow(Int a, Int b) {
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

template <typename Float, typename std::enable_if<std::is_floating_point<Float>::value, int>::type = 0>
Float safe_pow(Float a, Float b) {
	return std::pow(a,b);
}

//! Expects quoted string with protected characters and
//! removes quotes and protection
inline std::string unquote(char const* str) {
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
        }
        else if (*it == '"' && *(it + 1) == '\0') { break; } // NOLINT
        else if (*it == '\\') { slash = true; }
        else { res.push_back(*it); }
    }
    return res;
}

#endif // CLINGODL_UTIL_HH
