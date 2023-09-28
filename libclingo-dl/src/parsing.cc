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

#include "clingo-dl/parsing.hh"
#include "clingo-dl/util.hh"
#include <optional>

namespace ClingoDL {

namespace {

//! Negate a relation symbol.
auto negate_relation(char const *op) -> char const * {
    if (std::strcmp(op, "=") == 0) {
        return "!=";
    }
    if (std::strcmp(op, "!=") == 0) {
        return "=";
    }
    if (std::strcmp(op, "<") == 0) {
        return ">=";
    }
    if (std::strcmp(op, "<=") == 0) {
        return ">";
    }
    if (std::strcmp(op, ">") == 0) {
        return "<=";
    }
    if (std::strcmp(op, ">=") == 0) {
        return "<";
    }
    throw std::runtime_error("unexpected operator");
}

using ClingoDL::match;

//! Match if the given node represents a constant with the given name.
auto match_constant(Clingo::AST::Node const &ast, char const *name) -> bool {
    using namespace Clingo::AST;
    switch (ast.type()) {
        case Type::SymbolicTerm: {
            return ast.get<Clingo::Symbol>(Attribute::Term).match(name, 0);
        }
        case Type::Function: {
            if (ast.get<int>(Attribute::External) != 0) {
                return false;
            }
            if (!ast.get<NodeVector>(Attribute::Arguments).empty()) {
                return false;
            }
            return std::strcmp(ast.get<char const *>(Attribute::Name), name) == 0;
        }
        default: {
            return false;
        }
    }
}

//! Shift difference constraints in integrity constraints to the head.
auto shift_rule(Clingo::AST::Node ast) -> Clingo::AST::Node {
    using namespace Clingo::AST;
    if (ast.type() != Type::Rule) {
        return ast;
    }
    auto head = ast.get<Node>(Attribute::Head);
    if (head.type() != Type::Literal) {
        return ast;
    }
    auto atom = head.get<Node>(Attribute::Atom);
    if (atom.type() != Type::BooleanConstant) {
        return ast;
    }
    auto sign = head.get<int>(Attribute::Sign);
    auto value = atom.get<int>(Attribute::Value);
    if ((value == 0 && static_cast<Sign>(sign) == Sign::Negation) ||
        (value == 1 && static_cast<Sign>(sign) != Sign::Negation)) {
        return ast;
    }

    auto body = ast.get<NodeVector>(Attribute::Body);
    for (auto it = body.begin(), ie = body.end(); it != ie; ++it) {
        Node lit = *it;
        if (lit.type() != Type::Literal) {
            continue;
        }
        auto atom = lit.get<Node>(Attribute::Atom);
        if (atom.type() != Type::TheoryAtom || !match_constant(atom.get<Node>(Attribute::Term), "diff")) {
            continue;
        }
        auto ret = ast.copy();
        auto ret_bd = ret.get<NodeVector>(Attribute::Body);
        auto jt = ret_bd.begin() + (it - body.begin());
        auto ret_hd = atom.copy();
        auto guard = ret_hd.get<Clingo::Optional<Node>>(Attribute::Guard);
        check_syntax(guard.get() != nullptr);
        if (static_cast<Sign>(lit.get<int>(Attribute::Sign)) != Sign::Negation) {
            auto const *rel = guard->get<char const *>(Attribute::OperatorName);
            auto ret_guard = guard->copy();
            ret_guard.set(Attribute::OperatorName, negate_relation(rel));
            ret_hd.set(Attribute::Guard, Clingo::Optional<Node>{std::move(ret_guard)});
        }
        ret.set(Attribute::Head, std::move(ret_hd));
        ret_bd.erase(jt);
        return ret;
    }
    return ast;
}

//! Tag terms depending on whether they occur in heads or bodies.
auto tag_terms(Clingo::AST::Node &ast, char const *tag) -> Clingo::AST::Node {
    using namespace Clingo::AST;
    if (ast.type() == Type::SymbolicTerm) {
        auto ret = ast.copy();
        auto term = ast.get<Clingo::Symbol>(Attribute::Term);
        check_syntax(term.type() == Clingo::SymbolType::Function);
        std::string name{"__"};
        name += term.name();
        name += tag;
        ast.set(Attribute::Symbol, Clingo::Function(name.c_str(), {}));
        return ret;
    }
    if (ast.type() == Type::Function) {
        auto ret = ast.copy();
        std::string name{"__"};
        name += ret.get<char const *>(Attribute::Name);
        name += tag;
        ret.set(Attribute::Name, Clingo::add_string(name.c_str()));
        return ret;
    }
    return throw_syntax_error<Node>();
}

//! Tag head/body theory atoms.
struct TheoryRewriter {
    auto operator()(Clingo::AST::Node const &ast) -> Clingo::AST::Node {
        using namespace Clingo::AST;
        if (ast.type() == Type::Literal) {
            in_literal = true;
            auto ret = ast.transform_ast(*this);
            in_literal = false;
            return ret;
        }
        if (ast.type() == Type::TheoryAtom) {
            auto term = ast.get<Node>(Attribute::Term);
            if (match_constant(term, "diff")) {
                auto atom = ast.copy();

                auto elements = atom.get<NodeVector>(Attribute::Elements);
                check_syntax(elements.size() == 1);
                Clingo::AST::Node element = *elements.begin();
                auto tuple = element.get<NodeVector>(Attribute::Terms);
                check_syntax(tuple.size() == 1);
                auto condition = element.get<NodeVector>(Attribute::Condition);
                check_syntax(condition.empty());

                atom.set(Attribute::Term, tag_terms(term, in_literal ? "_b" : "_h"));

                return atom;
            }
        }
        return ast.transform_ast(*this);
    }
    bool in_literal{false};
};

//! Index that represents an invalid variable.
constexpr int INVALID_VAR{std::numeric_limits<int>::max()};

//! Test whether a variable is valid.
[[nodiscard]] inline auto is_valid_var(int var) -> bool { return var < INVALID_VAR; }

//! Convert a symbol to a double or integer.
template <class T> [[nodiscard]] auto to_number(Clingo::Symbol const &a) -> T {
    if (a.type() == Clingo::SymbolType::Number) {
        return static_cast<T>(a.number());
    }
    if (a.type() == Clingo::SymbolType::String) {
        return std::stod(a.string());
    }
    return throw_syntax_error<T>();
}

//! Evaluate a theory term to a number (represented by a symbol).
template <class N> [[nodiscard]] auto evaluate(Clingo::TheoryTerm const &term) -> Clingo::Symbol;

//! Evaluate two theory terms involved involved in a binary operation to an integral number.
template <class N, class F, typename std::enable_if<std::is_integral_v<N>, bool>::type = true>
[[nodiscard]] auto evaluate_binary(Clingo::TheoryTerm const &a, Clingo::TheoryTerm const &b, F &&f) -> Clingo::Symbol {
    auto ea = evaluate<N>(a);
    check_syntax(ea.type() == Clingo::SymbolType::Number);
    auto eb = evaluate<N>(b);
    check_syntax(eb.type() == Clingo::SymbolType::Number);
    return Clingo::Number(f(to_number<N>(ea), to_number<N>(eb)));
}

//! Evaluate two theory terms involved involved in a binary operation to a floating point number.
template <class N, class F, typename std::enable_if<std::is_floating_point_v<N>, bool>::type = true>
[[nodiscard]] auto evaluate_binary(Clingo::TheoryTerm const &a, Clingo::TheoryTerm const &b, F &&f) -> Clingo::Symbol {
    auto ea = evaluate<N>(a);
    auto eb = evaluate<N>(b);
    return Clingo::String(std::to_string(f(to_number<N>(ea), to_number<N>(eb))).c_str());
}

//! Parse a string to a number.
template <class N> [[nodiscard]] auto parse_number(char const *name) -> std::optional<N> {
    static const std::string chars = "\"";
    auto len = std::strlen(name);
    if (len < 2 || name[0] != '"' || name[len - 1] != '"') { // NOLINT
        return std::nullopt;
    }
    char *parsed = nullptr;                    // NOLINT
    auto res = std::strtod(name + 1, &parsed); // NOLINT
    if (parsed != name + len - 1) {            // NOLINT
        return std::nullopt;
    }
    auto ret = static_cast<N>(res);
    if (static_cast<double>(ret) != res) {
        throw std::runtime_error("could not evaluate term: for floating point numbers use option rdl");
    }
    return ret;
}

template <class N> auto evaluate(Clingo::TheoryTerm const &term) -> Clingo::Symbol {
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
void parse_elem(Clingo::TheoryTerm const &term, std::function<int(Clingo::Symbol)> const &map_vert,
                CoVarVec<N> &res) { // NOLINT
    if (term.type() == Clingo::TheoryTermType::Number) {
        res.emplace_back(term.number(), INVALID_VAR);
    } else if (match(term, "+", 2)) {
        auto args = term.arguments();
        parse_elem(args.front(), map_vert, res);
        parse_elem(args.back(), map_vert, res);
    } else if (match(term, "-", 2)) {
        auto args = term.arguments();
        parse_elem(args.front(), map_vert, res);
        auto pos = res.size();
        parse_elem(args.back(), map_vert, res);
        for (auto it = res.begin() + pos, ie = res.end(); it != ie; ++it) {
            it->first = safe_inv(it->first);
        }
    } else if (match(term, "-", 1)) {
        auto pos = res.size();
        parse_elem(term.arguments().front(), map_vert, res);
        for (auto it = res.begin() + pos, ie = res.end(); it != ie; ++it) {
            it->first = safe_inv(it->first);
        }
    } else if (match(term, "+", 1)) {
        parse_elem(term.arguments().front(), map_vert, res);
    } else if (match(term, "*", 2)) {
        auto args = term.arguments();
        CoVarVec<N> lhs;
        parse_elem(args.front(), map_vert, lhs);
        CoVarVec<N> rhs;
        parse_elem(args.back(), map_vert, rhs);
        for (auto &l : lhs) {
            for (auto &r : rhs) {
                if (!is_valid_var(l.second)) {
                    res.emplace_back(safe_mul(l.first, r.first), r.second);
                } else if (!is_valid_var(r.second)) {
                    res.emplace_back(safe_mul(l.first, r.first), l.second);
                } else {
                    throw_syntax_error("Invalid Syntax: only linear difference constraints are supported");
                }
            }
        }
    } else if (term.type() == Clingo::TheoryTermType::Symbol) {
        if (auto val = parse_number<N>(term.name()); val) {
            res.emplace_back(*val, INVALID_VAR);
        } else {
            res.emplace_back(1, map_vert(evaluate<N>(term)));
        }
    } else if (term.type() == Clingo::TheoryTermType::Function || term.type() == Clingo::TheoryTermType::Tuple) {
        res.emplace_back(1, map_vert(evaluate<N>(term)));
    } else {
        throw_syntax_error("Invalid Syntax: invalid diff constraint");
    }
}

//! Simplify the given vector of terms.
template <class N> [[nodiscard]] auto simplify(CoVarVec<N> &vec) -> N {
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
        } else {
            auto r = seen.emplace(var, jt);
            auto kt = r.first;
            auto ins = r.second;
            if (!ins) {
                kt->second->first = safe_add<N>(kt->second->first, co);
            } else {
                if (it != jt) {
                    *jt = *it;
                }
                ++jt;
            }
        }
    }

    jt = std::remove_if(vec.begin(), jt, [](auto &co_var) { return co_var.first == 0; });
    vec.erase(jt, vec.end());
    return rhs;
}

} // namespace

auto match(Clingo::TheoryTerm const &term, char const *name, size_t arity) -> bool {
    return (term.type() == Clingo::TheoryTermType::Symbol && std::strcmp(term.name(), name) == 0 && arity == 0) ||
           (term.type() == Clingo::TheoryTermType::Function && std::strcmp(term.name(), name) == 0 &&
            term.arguments().size() == arity);
}

void transform(Clingo::AST::Node const &ast, NodeCallback const &cb, bool shift) {
    for (auto &unpooled : ast.unpool()) {
        if (shift) {
            unpooled = shift_rule(unpooled);
        }
        cb(unpooled.transform_ast(TheoryRewriter{}));
    }
}

template <class N>
auto parse(Clingo::TheoryAtom const &atom, std::function<int(Clingo::Symbol)> const &map_vert) -> EdgeAtom<N> {
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

template EdgeAtom<int> parse<int>(Clingo::TheoryAtom const &, std::function<int(Clingo::Symbol)> const &);
template EdgeAtom<double> parse<double>(Clingo::TheoryAtom const &, std::function<int(Clingo::Symbol)> const &);

} // namespace ClingoDL
