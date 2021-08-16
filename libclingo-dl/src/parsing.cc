// {{{ MIT License
//
// Copyright 2018 Roland Kaminski, Philipp Wanko, Max Ostrowski
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

#include "clingo-dl/parsing.hh"
#include "clingo-dl/util.hh"

namespace ClingoDL {

namespace {

//! Epsilon value for floating point arithmetics.
constexpr double DOUBLE_EPSILON = 0.00001;

//! Index that represents an invalid variable.
constexpr int INVALID_VAR{std::numeric_limits<int>::max()};

//! Test whether a variable is valid.
[[nodiscard]] inline bool is_valid_var(int var) {
    return var < INVALID_VAR;
}

//! Convert a symbol to a double or integer.
template <class T>
[[nodiscard]] T to_number(Clingo::Symbol const &a) {
    if (a.type() == Clingo::SymbolType::Number) {
        return static_cast<T>(a.number());
    }
    if (a.type() == Clingo::SymbolType::String) {
        return std::stod(a.string());
    }
    return throw_syntax_error<T>();
}

//! Evaluate a theory term to a number (represented by a symbol).
template <class N>
[[nodiscard]] Clingo::Symbol evaluate(Clingo::TheoryTerm const &term);

//! Evaluate two theory terms involved involved in a binary operation to an integral number.
template <class N, class F, typename std::enable_if<std::is_integral_v<N>, bool>::type = true>
[[nodiscard]] Clingo::Symbol evaluate_binary(Clingo::TheoryTerm const &a, Clingo::TheoryTerm const &b, F &&f) {
    auto ea = evaluate<N>(a);
    check_syntax(ea.type() == Clingo::SymbolType::Number);
    auto eb = evaluate<N>(b);
    check_syntax(eb.type() == Clingo::SymbolType::Number);
    return Clingo::Number(f(to_number<N>(ea), to_number<N>(eb)));
}

//! Evaluate two theory terms involved involved in a binary operation to a floating point number.
template <class N, class F, typename std::enable_if<std::is_floating_point_v<N>, bool>::type = true>
[[nodiscard]] Clingo::Symbol evaluate_binary(Clingo::TheoryTerm const &a, Clingo::TheoryTerm const &b, F &&f) {
    auto ea = evaluate<N>(a);
    auto eb = evaluate<N>(b);
    return Clingo::String(std::to_string(f(to_number<N>(ea), to_number<N>(eb))).c_str());
}

//! Epsilon value for integral numbers.
template <class N, typename std::enable_if<std::is_integral_v<N>, bool>::type = true>
[[nodiscard]] N epsilon_() {
    return 1;
}

//! Epsilon value for floating point numbers.
template <class N, typename std::enable_if<std::is_floating_point_v<N>, bool>::type = true>
[[nodiscard]] N epsilon_() {
    return DOUBLE_EPSILON;
}

//! Parse a string to a number.
template <class N>
[[nodiscard]] std::optional<N> parse_number(char const *name) {
    static const std::string chars = "\"";
    auto len = std::strlen(name);
    if (len < 2 || name[0] != '"' || name[len - 1] != '"') { // NOLINT
        return std::nullopt;
    }
    char *parsed = nullptr;
    auto res = std::strtod(name + 1, &parsed); // NOLINT
    if (parsed != name + len - 1) { // NOLINT
        return std::nullopt;
    }
    auto ret = static_cast<N>(res);
    if (static_cast<double>(ret) != res) {
        throw std::runtime_error("could not evaluate term: for floating point numbers use option rdl");
    }
    return ret;
}

template <class N>
Clingo::Symbol evaluate(Clingo::TheoryTerm const &term) {
    if (term.type() == Clingo::TheoryTermType::Symbol) {
        auto const *name = term.name();
        if (std::strncmp(name, "\"", 1) == 0) {
            return Clingo::String(unquote(name).c_str());
        }
        return Clingo::Function(name, {});
    }

    if (term.type() == Clingo::TheoryTermType::Number) {
        return Clingo::Number(term.number());
    }

    if (match(term, "+", 2)) {
        return evaluate_binary<N>(term.arguments().front(), term.arguments().back(), safe_add<N>);
    }
    if (match(term, "-", 2)) {
        return evaluate_binary<N>(term.arguments().front(), term.arguments().back(), safe_sub<N>);
    }
    if (match(term, "*", 2)) {
        return evaluate_binary<N>(term.arguments().front(), term.arguments().back(), safe_mul<N>);
    }
    if (match(term, "/", 2)) {
        return evaluate_binary<N>(term.arguments().front(), term.arguments().back(), safe_div<N>);
    }
    if (match(term, "\\", 2)) {
        return evaluate_binary<N>(term.arguments().front(), term.arguments().back(), safe_mod<N>);
    }
    if (match(term, "**", 2)) {
        return evaluate_binary<N>(term.arguments().front(), term.arguments().back(), safe_pow<N>);
    }

    if (match(term, "-", 1)) {
        auto ea = evaluate<N>(term.arguments().front());
        if (ea.type() == Clingo::SymbolType::Number) {
            return Clingo::Number(safe_inv(ea.number()));
        }
        if (ea.type() == Clingo::SymbolType::Function && std::strlen(ea.name()) > 0) {
            return Clingo::Function(ea.name(), ea.arguments(), !ea.is_positive());
        }
        return throw_syntax_error<Clingo::Symbol>();
    }

    check_syntax(!match(term, "..", 2));

    if (term.type() == Clingo::TheoryTermType::Tuple || term.type() == Clingo::TheoryTermType::Function) {
        std::vector<Clingo::Symbol> args;
        args.reserve(term.arguments().size());
        for (auto const &arg : term.arguments()) {
            args.emplace_back(evaluate<N>(arg));
        }
        return Clingo::Function(term.type() == Clingo::TheoryTermType::Function ? term.name() : "", args);
    }
    return throw_syntax_error<Clingo::Symbol>();
}

//! Parse the given theory term for an arithmetic expression.
template <class N>
void parse_elem(Clingo::TheoryTerm const &term, std::function<int (Clingo::Symbol)> const &map_vert, CoVarVec<N> &res) {
    if (term.type() == Clingo::TheoryTermType::Number) {
        res.emplace_back(term.number(), INVALID_VAR);
    }
    else if (match(term, "+", 2)) {
        auto args = term.arguments();
        parse_elem(args.front(), map_vert, res);
        parse_elem(args.back(), map_vert, res);
    }
    else if (match(term, "-", 2)) {
        auto args = term.arguments();
        parse_elem(args.front(), map_vert, res);
        auto pos = res.size();
        parse_elem(args.back(), map_vert, res);
        for (auto it = res.begin() + pos, ie = res.end(); it != ie; ++it) {
            it->first = safe_inv(it->first);
        }
    }
    else if (match(term, "-", 1)) {
        auto pos = res.size();
        parse_elem(term.arguments().front(), map_vert, res);
        for (auto it = res.begin() + pos, ie = res.end(); it != ie; ++it) {
            it->first = safe_inv(it->first);
        }
    }
    else if (match(term, "+", 1)) {
        parse_elem(term.arguments().front(), map_vert, res);
    }
    else if (match(term, "*", 2)) {
        auto args = term.arguments();
        CoVarVec<N> lhs;
        parse_elem(args.front(), map_vert, lhs);
        CoVarVec<N> rhs;
        parse_elem(args.back(), map_vert, rhs);
        for (auto &l : lhs) {
            for (auto &r : rhs) {
                if (!is_valid_var(l.second)) {
                    res.emplace_back(safe_mul(l.first, r.first), r.second);
                }
                else if (!is_valid_var(r.second)) {
                    res.emplace_back(safe_mul(l.first, r.first), l.second);
                }
                else {
                    throw_syntax_error("Invalid Syntax: only linear difference constraints are supported");
                }
            }
        }
    }
    else if (term.type() == Clingo::TheoryTermType::Symbol) {
        if (auto val = parse_number<N>(term.name()); val) {
            res.emplace_back(*val, INVALID_VAR);
        }
        else {
            res.emplace_back(1, map_vert(evaluate<N>(term)));
        }
    }
    else if (term.type() == Clingo::TheoryTermType::Function || term.type() == Clingo::TheoryTermType::Tuple) {
        res.emplace_back(1, map_vert(evaluate<N>(term)));
    }
    else {
        throw_syntax_error("Invalid Syntax: invalid diff constraint");
    }
}

//! Simplify the given vector of terms.
template <class N>
[[nodiscard]] N simplify(CoVarVec<N> &vec) {
    static thread_local std::unordered_map<int, typename CoVarVec<N>::iterator> seen;
    N rhs = 0;
    seen.clear();

    auto jt = vec.begin();
    for (auto it = jt, ie = vec.end(); it != ie; ++it) {
        auto &[co, var] = *it;
        if (co == 0) {
            continue;
        }
        if (!is_valid_var(var)) {
            rhs = safe_sub<N>(rhs, co);
        }
        else {
            auto r = seen.emplace(var, jt);
            auto kt = r.first;
            auto ins = r.second;
            if (!ins) {
                kt->second->first = safe_add<N>(kt->second->first, co);
            }
            else {
                if (it != jt) {
                    *jt = *it;
                }
                ++jt;
            }
        }
    }

    jt = std::remove_if(vec.begin(), jt, [](auto &co_var) { return co_var.first == 0; } );
    vec.erase(jt, vec.end());
    return rhs;
}

} // namespace

template <class N>
N epsilon() {
    return epsilon_<N>();
}

bool match(Clingo::TheoryTerm const &term, char const *name, size_t arity) {
    return (term.type() == Clingo::TheoryTermType::Symbol &&
        std::strcmp(term.name(), name) == 0 &&
        arity == 0) ||
        (term.type() == Clingo::TheoryTermType::Function &&
        std::strcmp(term.name(), name) == 0 &&
        term.arguments().size() == arity);
}

template <class N>
EdgeAtom<N> parse(Clingo::TheoryAtom const &atom, std::function<int (Clingo::Symbol)> const &map_vert) {
    char const *msg = "parsing difference constraint failed: only constraints of form &diff {u - v} <= b are accepted";
    if (!atom.has_guard()) {
        throw std::runtime_error(msg);
    }
    auto term = atom.term();
    bool strict = match(term, "__diff_b", 0);
    if (strict && std::is_floating_point_v<N>) {
        throw std::runtime_error("strict semantics not available with floating point numbers");
    }
    CoVarVec<N> covec;
    parse_elem(atom.guard().second, map_vert, covec);
    for (auto &[co, var] : covec) {
        co = safe_inv<N>(co);
    }
    auto const *rel = atom.guard().first;

    auto elems = atom.elements();
    if (elems.size() > 1) {
        throw std::runtime_error(msg);
    }
    for (auto const &element : elems) {
        auto tuple = element.tuple();
        check_syntax(!tuple.empty() && element.condition().empty(), "Invalid Syntax: invalid sum constraint");
        parse_elem(element.tuple().front(), map_vert, covec);
    }

    auto rhs = simplify(covec);
    return {std::move(covec), rel, rhs, atom.literal(), strict};
}

template int epsilon<int>();
template double epsilon<double>();

template EdgeAtom<int> parse<int>(Clingo::TheoryAtom const &, std::function<int (Clingo::Symbol)> const &);
template EdgeAtom<double> parse<double>(Clingo::TheoryAtom const &, std::function<int (Clingo::Symbol)> const &);

} // namespace ClingoDL
