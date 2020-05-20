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

#ifndef CLINGODL_UTIL_HH
#define CLINGODL_UTIL_HH

#include <chrono>
#include <ostream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>

namespace std {

template<>
struct hash<std::pair<int, int>> {
    size_t operator()(std::pair<int, int> const &p) const {
        return static_cast<size_t>(p.first) + (static_cast<size_t>(p.second) >> 32);
    }
};

}

namespace Detail {

template <int X>
using int_type = std::integral_constant<int, X>;
template <class T, class S>
inline void nc_check(S s, int_type<0>) { // same sign
    (void)s;
    assert((std::is_same<T, S>::value) || (s >= std::numeric_limits<T>::min() && s <= std::numeric_limits<T>::max()));
}
template <class T, class S>
inline void nc_check(S s, int_type<-1>) { // Signed -> Unsigned
    (void)s;
    assert(s >= 0 && static_cast<S>(static_cast<T>(s)) == s);
}
template <class T, class S>
inline void nc_check(S s, int_type<1>) { // Unsigned -> Signed
    (void)s;
    assert(!(s > std::numeric_limits<T>::max()));
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
    ~Timer() { elapsed_ += std::chrono::steady_clock::now() - start_; }

private:
    Duration &elapsed_;
    std::chrono::time_point<std::chrono::steady_clock> start_;
};

template <int N>
class Heap {
public:
    template <class M>
    void push(M &m, int item) {
        auto i = m.offset(item) = static_cast<int>(heap_.size());
        heap_.push_back(item);
        decrease(m, i);
    }
    template <class M>
    int pop(M &m) {
        assert(!heap_.empty());
        auto ret = heap_[0];
        if (heap_.size() > 1) {
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
    void decrease(M &m, int i) {
        while (i > 0) {
            int p = parent_(i);
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
    void increase(M &m, int i) {
        for (int p = i, j = children_(p), s = numeric_cast<int>(heap_.size()); j < s; j = children_(p)) {
            int min = j;
            for (int k = j + 1; k < j + N; ++k) {
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
    int size() { return heap_.size(); }
    bool empty() { return heap_.empty(); }
    void clear() { heap_.clear(); }

private:
    template <class M>
    void swap_(M &m, int i, int j) {
        m.offset(heap_[j]) = i;
        m.offset(heap_[i]) = j;
        std::swap(heap_[i], heap_[j]);
    }
    int parent_(int offset) { return (offset - 1) / N; }
    int children_(int offset) { return N * offset + 1; }
    template <class M>
    bool less_(M &m, int a, int b) {
        a = heap_[a], b = heap_[b];
        auto ca = m.cost(a), cb = m.cost(b);
        return ca < cb || (ca == cb && m.relevant(a) < m.relevant(b));
    }

private:
    std::vector<int> heap_;
};


// Some of the functions below could also be implemented using (much faster)
// compiler specific built-ins. For more information check the following links:
// - https://gcc.gnu.org/onlinedocs/gcc/Integer-Overflow-Builtins.html
// - https://wiki.sei.cmu.edu/confluence/display/c/INT32-C.+Ensure+that+operations+on+signed+integers+do+not+result+in+overflow

//! Safely add a and b throwing an exception in case of overflow/underflow.
template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int> = 0> 
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

template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int> = 0>
Float safe_add(Float a, Float b) {
	return a + b; 
}

//! Safely subtract a and b throwing an exception in case of
//! overflow/underflow.
template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int> = 0> 
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

template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int> = 0>
Float safe_sub(Float a, Float b) {
	return a - b; 
}

//! Safely multiply a and b throwing an exception in case of
//! overflow/underflow.
template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int> = 0> 
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

template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int> = 0>
Float safe_mul(Float a, Float b) {
	return a * b; 
}

//! Safely divide a and b throwing an exception in case of overflow/underflow.
template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int> = 0> 
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

template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int> = 0>
Float safe_div(Float a, Float b) {
	return a / b; 
}


//! Safely calculate the modulo of a and b throwing an exception in case of
//! overflow/underflow.
template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int> = 0> 
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

template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int> = 0>
Float safe_mod(Float a, Float b) {
	return fmod(a,b); 
}

//! Safely invert a throwing an exception in case of an underflow.
template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int> = 0> 
Int safe_inv(Int a) {
    if (a == std::numeric_limits<Int>::min()) {
        throw std::overflow_error("integer overflow");
    }
    return -a;
}

template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int> = 0>
Float safe_inv(Float a) {
	return -a; 
}

template <typename Int, std::enable_if_t<std::is_integral_v<Int>, int> = 0> 
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

template <typename Float, std::enable_if_t<std::is_floating_point_v<Float>, int> = 0>
Float safe_pow(Float a, Float b) {
	return std::pow(a,b); 
}




#endif // CLINGODL_UTIL_HH
