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
auto negate_relation(std::string_view op) -> std::string_view {
    if (op == "=") {
        return "!=";
    }
    if (op == "!=") {
        return "=";
    }
    if (op == "<") {
        return ">=";
    }
    if (op == "<=") {
        return ">";
    }
    if (op == ">") {
        return "<=";
    }
    if (op == ">=") {
        return "<";
    }
    return throw_syntax_error<std::string_view>("unexpected operator");
}

//! Match if the given node represents a constant with the given name.
auto match_constant(Clingo::AST::Node const &ast, std::string_view name) -> bool {
    using namespace Clingo::AST;
    switch (ast.type()) {
        case NodeType::term_symbolic: {
            return ast.symbol(Attribute::symbol).match(name, 0);
        }
        case NodeType::term_function: {
            check_syntax(ast.number(Attribute::external) == 0, "theory atom names must not be external");
            if (ast.string(Attribute::name) != name) {
                return false;
            }
            auto pool = ast.nodes(Attribute::pool);
            check_syntax(pool.size() == 1, "theory atom names must not contain pools");
            return pool.front().nodes(Attribute::arguments).empty();
        }
        default: {
            return false;
        }
    }
}

//! Shift difference constraints in integrity constraints to the head.
//!
//! Moves the first difference constraint in the body of a rule into the head of an
//! integrity constraint.
auto shift_rule(Clingo::Library &lib, Clingo::AST::Node ast) -> Clingo::AST::Node {
    using namespace Clingo::AST;
    if (ast.type() != NodeType::statement_rule) {
        return ast;
    }
    auto head = ast.node(Attribute::head);
    if (head.type() != NodeType::head_simple_literal) {
        return ast;
    }
    auto lit = head.node(Attribute::literal);
    if (lit.type() != NodeType::literal_boolean) {
        return ast;
    }
    auto sign = lit.number(Attribute::sign);
    auto value = lit.number(Attribute::value);
    if ((value == 0 && sign == Sign::single) || (value == 1 && sign != Sign::single)) {
        return ast;
    }

    auto body = ast.nodes(Attribute::body);
    for (auto it = body.begin(), ie = body.end(); it != ie; ++it) {
        auto lit = *it;
        if (lit.type() != NodeType::body_theory_atom) {
            continue;
        }
        auto name = lit.node(Attribute::name);
        if (!match_constant(name, "diff")) {
            continue;
        }
        auto elems = lit.nodes(Attribute::elements);
        auto guard = lit.optional_node(Attribute::right);
        check_syntax(guard.has_value());
        if (lit.number(Attribute::sign) != Sign::single) {
            guard = guard->update<NodeType::theory_right_guard>(lib, [&]<Attribute Attr>() {
                if constexpr (Attr == Attribute::theory_operator) {
                    return negate_relation(guard->string(Attribute::theory_operator));
                }
            });
        }
        body.erase(it);
        return ast.update<NodeType::statement_rule>(lib, [&]<Attribute Attr>() {
            if constexpr (Attr == Attribute::head) {
                return Node::create<NodeType::head_theory_atom>(lib, lit.location(Attribute::location), name, elems,
                                                                guard);
            }
            if constexpr (Attr == Attribute::body) {
                return std::move(body);
            }
        });
    }
    return ast;
}

//! Tag terms depending on whether they occur in heads or bodies.
auto tag_terms(Clingo::Library &lib, Clingo::AST::Node &ast, std::string_view tag) -> Clingo::AST::Node {
    using namespace Clingo::AST;
    auto tag_string = [tag](std::string_view str) {
        auto res = std::string{"__"};
        res.reserve(res.size() + str.size() + tag.size());
        res += str;
        res += tag;
        return res;
    };
    if (ast.type() == NodeType::term_symbolic) {
        auto term = ast.symbol(Attribute::symbol);
        assert(term.match("diff", 0));
        auto sym = Clingo::Function(lib, tag_string(term.name()), {});
        return Node::create<NodeType::term_symbolic>(lib, ast.location(Attribute::location), sym);
    }
    if (ast.type() == NodeType::term_function) {
        return ast.update<NodeType::term_function>(lib, [&]<Attribute Attr>() {
            if constexpr (Attr == Attribute::name) {
                return tag_string(ast.string(Attribute::name));
            }
        });
    }
    return throw_syntax_error<Node>();
}

auto rewrite_theory(Clingo::Library &lib, Clingo::AST::Node const &ast) -> std::optional<Clingo::AST::Node> {
    using namespace Clingo::AST;
    using T = NodeType;
    using A = Attribute;
    Transformer trans = [&](Clingo::AST::Node const &ast) -> std::optional<Node> {
        auto update = [&]<T N>() -> std::optional<Node> {
            auto term = ast.node(A::name);
            if (match_constant(term, "diff")) {
                return ast.update<N>(lib, [&]<A attr>() {
                    if constexpr (attr == A::elements) {
                        auto elements = ast.nodes(attr);
                        check_syntax(elements.size() == 1);
                        Clingo::AST::Node element = *elements.begin();
                        auto tuple = element.nodes(A::terms);
                        check_syntax(tuple.size() == 1);
                        auto condition = element.nodes(Attribute::condition);
                        check_syntax(condition.empty());
                    }
                    if constexpr (attr == A::name) {
                        return tag_terms(lib, term, N == T::body_theory_atom ? "_b" : "_h");
                    }
                });
            }
            return std::nullopt;
        };
        if (ast.type() == T::body_theory_atom) {
            return update.template operator()<T::body_theory_atom>();
        }
        if (ast.type() == T::head_theory_atom) {
            return update.template operator()<T::head_theory_atom>();
        }
        return ast.accept(lib, trans);
    };
    return trans(ast);
}

//! Index that represents an invalid variable.
constexpr int INVALID_VAR{std::numeric_limits<int>::max()};

//! Test whether a variable is valid.
[[nodiscard]] inline auto is_valid_var(int var) -> bool { return var < INVALID_VAR; }

//! Parse a string to a number.
template <class T> [[nodiscard]] auto parse_number(std::string_view name) -> std::optional<T> {
    T res = 0;
    auto end = name.data() + name.size();
    auto [ptr, ec] = std::from_chars(name.data(), end, res);
    if (ec != std::errc{} || ptr != end) {
        return std::nullopt;
    }
    return res;
}

//! Convert a symbol to a double or integer.
template <class T> [[nodiscard]] auto to_number(Clingo::Symbol const &a) -> T {
    if (a.type() == Clingo::SymbolType::number) {
        return static_cast<T>(a.number());
    }
    if (a.type() == Clingo::SymbolType::string) {
        if (auto res = parse_number<T>(a.string())) {
            return *res;
        }
    }
    return throw_syntax_error<T>("failed to parse number");
}

//! Evaluate a theory term to a number (represented by a symbol).
template <class N> [[nodiscard]] auto evaluate(Clingo::Library &lib, Clingo::TheoryTerm const &term) -> Clingo::Symbol;

//! Evaluate two theory terms involved involved in a binary operation to an integral number.
template <class N, class F, typename std::enable_if<std::is_integral_v<N>, bool>::type = true>
[[nodiscard]] auto evaluate_binary(Clingo::Library &lib, Clingo::TheoryTerm const &a, Clingo::TheoryTerm const &b,
                                   F &&f) -> Clingo::Symbol {
    auto ea = evaluate<N>(lib, a);
    check_syntax(ea.type() == Clingo::SymbolType::number);
    auto eb = evaluate<N>(lib, b);
    check_syntax(eb.type() == Clingo::SymbolType::number);
    return Clingo::Number(f(to_number<N>(ea), to_number<N>(eb)));
}

//! Evaluate two theory terms involved involved in a binary operation to a floating point number.
template <class N, class F, typename std::enable_if<std::is_floating_point_v<N>, bool>::type = true>
[[nodiscard]] auto evaluate_binary(Clingo::Library &lib, Clingo::TheoryTerm const &a, Clingo::TheoryTerm const &b,
                                   F &&f) -> Clingo::Symbol {
    auto ea = evaluate<N>(lib, a);
    auto eb = evaluate<N>(lib, b);
    return Clingo::String(lib, std::to_string(f(to_number<N>(ea), to_number<N>(eb))));
}

template <class N> auto evaluate(Clingo::Library &lib, Clingo::TheoryTerm const &term) -> Clingo::Symbol {
    if (term.type() == Clingo::TheoryTermType::symbol) {
        auto name = term.name();
        if (name.starts_with('"')) {
            return Clingo::String(lib, unquote(name));
        }
        return Clingo::Function(lib, name, {});
    }

    if (term.type() == Clingo::TheoryTermType::number) {
        return Clingo::Number(term.number());
    }

    if (match(term, "+", 2)) {
        return evaluate_binary<N>(lib, term.arguments().front(), term.arguments().back(), safe_add<N>);
    }
    if (match(term, "-", 2)) {
        return evaluate_binary<N>(lib, term.arguments().front(), term.arguments().back(), safe_sub<N>);
    }
    if (match(term, "*", 2)) {
        return evaluate_binary<N>(lib, term.arguments().front(), term.arguments().back(), safe_mul<N>);
    }
    if (match(term, "/", 2)) {
        return evaluate_binary<N>(lib, term.arguments().front(), term.arguments().back(), safe_div<N>);
    }
    if (match(term, "\\", 2)) {
        return evaluate_binary<N>(lib, term.arguments().front(), term.arguments().back(), safe_mod<N>);
    }
    if (match(term, "**", 2)) {
        return evaluate_binary<N>(lib, term.arguments().front(), term.arguments().back(), safe_pow<N>);
    }

    if (match(term, "-", 1)) {
        auto ea = evaluate<N>(lib, term.arguments().front());
        if (ea.type() == Clingo::SymbolType::number) {
            return Clingo::Number(safe_inv(ea.number()));
        }
        if (ea.type() == Clingo::SymbolType::function) {
            return Clingo::Function(lib, ea.name(), ea.arguments(), !ea.is_positive());
        }
        return throw_syntax_error<Clingo::Symbol>();
    }

    check_syntax(!match(term, "..", 2));

    if (term.type() == Clingo::TheoryTermType::tuple || term.type() == Clingo::TheoryTermType::function) {
        std::vector<Clingo::Symbol> args;
        args.reserve(term.arguments().size());
        for (auto const &arg : term.arguments()) {
            args.emplace_back(evaluate<N>(lib, arg));
        }
        return Clingo::Function(lib, term.type() == Clingo::TheoryTermType::function ? term.name() : "", args);
    }
    return throw_syntax_error<Clingo::Symbol>();
}

//! Parse the given theory term for an arithmetic expression.
template <class N>
void parse_elem(Clingo::Library &lib, Clingo::TheoryTerm const &term,
                std::function<int(Clingo::Symbol)> const &map_vert,
                CoVarVec<N> &res) { // NOLINT
    if (term.type() == Clingo::TheoryTermType::number) {
        res.emplace_back(term.number(), INVALID_VAR);
    } else if (match(term, "+", 2)) {
        auto args = term.arguments();
        parse_elem(lib, args.front(), map_vert, res);
        parse_elem(lib, args.back(), map_vert, res);
    } else if (match(term, "-", 2)) {
        auto args = term.arguments();
        parse_elem(lib, args.front(), map_vert, res);
        auto pos = res.size();
        parse_elem(lib, args.back(), map_vert, res);
        for (auto it = res.begin() + pos, ie = res.end(); it != ie; ++it) {
            it->first = safe_inv(it->first);
        }
    } else if (match(term, "-", 1)) {
        auto pos = res.size();
        parse_elem(lib, term.arguments().front(), map_vert, res);
        for (auto it = res.begin() + pos, ie = res.end(); it != ie; ++it) {
            it->first = safe_inv(it->first);
        }
    } else if (match(term, "+", 1)) {
        parse_elem(lib, term.arguments().front(), map_vert, res);
    } else if (match(term, "*", 2)) {
        auto args = term.arguments();
        CoVarVec<N> lhs;
        parse_elem(lib, args.front(), map_vert, lhs);
        CoVarVec<N> rhs;
        parse_elem(lib, args.back(), map_vert, rhs);
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
    } else if (term.type() == Clingo::TheoryTermType::symbol) {
        auto name = term.name();
        if (auto val = name.starts_with('"') ? parse_number<N>(unquote(name)) : std::nullopt; val) {
            res.emplace_back(*val, INVALID_VAR);
        } else {
            res.emplace_back(1, map_vert(evaluate<N>(lib, term)));
        }
    } else if (term.type() == Clingo::TheoryTermType::function || term.type() == Clingo::TheoryTermType::tuple) {
        res.emplace_back(1, map_vert(evaluate<N>(lib, term)));
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

auto match(Clingo::TheoryTerm const &term, std::string_view name, size_t arity) -> bool {
    return (term.type() == Clingo::TheoryTermType::symbol && term.name() == name && arity == 0) ||
           (term.type() == Clingo::TheoryTermType::function && term.name() == name && term.arguments().size() == arity);
}

void rewrite(Clingo::Library &lib, Clingo::AST::Node ast, NodeCallback const &cb, bool shift) {
    if (shift) {
        ast = shift_rule(lib, ast);
    }
    if (auto res = rewrite_theory(lib, ast); res) {
        cb(*std::move(res));
    } else {
        cb(std::move(ast));
    }
}

template <class N>
auto parse(Clingo::Library &lib, Clingo::TheoryAtom const &atom,
           std::function<int(Clingo::Symbol const &)> const &map_vert) -> EdgeAtom<N> {
    char const *msg = "parsing difference constraint failed";
    auto guard = atom.guard();
    if (!guard) {
        throw_syntax_error(msg);
    }
    auto term = atom.name();
    bool strict = match(term, "__diff_b", 0);
    if (strict && std::is_floating_point_v<N>) {
        throw_syntax_error("strict semantics not available with floating point numbers");
    }
    CoVarVec<N> covec;
    parse_elem(lib, guard->second, map_vert, covec);
    for (auto &[co, var] : covec) {
        co = safe_inv<N>(co);
    }
    auto rel = guard->first;

    auto elems = atom.elements();
    if (elems.size() > 1) {
        throw std::runtime_error(msg);
    }
    for (auto const &element : elems) {
        auto tuple = element.tuple();
        check_syntax(!tuple.empty() && element.condition().empty(), "invalid diff constraint");
        parse_elem(lib, element.tuple().front(), map_vert, covec);
    }

    auto rhs = simplify(covec);
    return {std::move(covec), relation_from_string(rel), rhs, atom.literal(), strict};
}

template EdgeAtom<int> parse<int>(Clingo::Library &lib, Clingo::TheoryAtom const &,
                                  std::function<int(Clingo::Symbol const &)> const &);
template EdgeAtom<double> parse<double>(Clingo::Library &lib, Clingo::TheoryAtom const &,
                                        std::function<int(Clingo::Symbol const &)> const &);

} // namespace ClingoDL
