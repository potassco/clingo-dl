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

#include <clingo-dl/propagator.hh>

#include <unordered_set>

namespace ClingoDL {

char const *negate_relation(char const *op) {
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

template <typename CStr>
bool match(Clingo::ASTv2::AST const &ast, CStr str) {
    using namespace Clingo::ASTv2;
    if (ast.type() == Type::SymbolicTerm) {
        auto sym = ast.get<Clingo::Symbol>(Attribute::Term);
        return sym.match(str, 0);
    }
    if (ast.type() == Type::Function) {
        if (ast.get<int>(Attribute::External) != 0) {
            return false;
        }
        if (!ast.get<ASTVector>(Attribute::Arguments).empty()) {
            return false;
        }
        auto const *name = ast.get<char const *>(Attribute::Name);
        return (std::strcmp(name, str) == 0);
    }
    return false;
}

Clingo::ASTv2::AST shift_rule(Clingo::ASTv2::AST ast) {
    using namespace Clingo::ASTv2;
    if (ast.type() != Type::Rule) {
        return ast;
    }
    auto head = ast.get<AST>(Attribute::Head);
    if (head.type() != Type::Literal) {
        return ast;
    }
    auto atom = head.get<AST>(Attribute::Atom);
    if (atom.type() != Type::BooleanConstant) {
        return ast;
    }
    auto sign = head.get<int>(Attribute::Sign);
    auto value = atom.get<int>(Attribute::Value);
    if ((value == 0 && static_cast<Sign>(sign) == Sign::Negation) ||
        (value == 1 && static_cast<Sign>(sign) != Sign::Negation)) {
        return ast;
    }

    auto body = ast.get<ASTVector>(Attribute::Body);
    for (auto it = body.begin(), ie = body.end(); it != ie; ++it) {
        AST lit = *it;
        auto atom = lit.get<AST>(Attribute::Atom);
        if (atom.type() == Type::TheoryAtom) {
            if (match(atom.get<AST>(Attribute::Term), "diff")) {
                auto ret = ast.copy();
                auto ret_bd = ret.get<ASTVector>(Attribute::Body);
                auto jt = ret_bd.begin() + (it - body.begin());
                auto ret_hd = atom.copy();
                auto guard = ret_hd.get<Clingo::Optional<AST>>(Attribute::Guard);
                check_syntax(guard.get() != nullptr);
                if (static_cast<Sign>(lit.get<int>(Attribute::Sign)) != Sign::Negation) {
                    auto const *rel = guard->get<char const *>(Attribute::OperatorName);
                    auto ret_guard = guard->copy();
                    ret_guard.set(Attribute::OperatorName, negate_relation(rel));
                    ret_hd.set(Attribute::Guard, Clingo::Optional<AST>{std::move(ret_guard)});
                }
                ret.set(Attribute::Head, std::move(ret_hd));
                ret_bd.erase(jt);
                return ret;
            }
        }
    }
    return ast;
}

Clingo::ASTv2::AST tag_terms(Clingo::ASTv2::AST &ast, char const *tag) {
    using namespace Clingo::ASTv2;
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
        name += ret.get<char const*>(Attribute::Name);
        name += tag;
        ret.set(Attribute::Name, Clingo::add_string(name.c_str()));
        return ret;
    }
    return throw_syntax_error<AST>();
}

// Mark head/body theory atoms.
struct TheoryRewriterR {
    Clingo::ASTv2::AST operator()(Clingo::ASTv2::AST const &ast) {
        using namespace Clingo::ASTv2;
        if (ast.type() == Type::Literal) {
            in_literal = true;
            auto ret = ast.transform_ast(*this);
            in_literal = false;
            return ret;
        }
        if (ast.type() == Type::TheoryAtom) {
            auto term = ast.get<AST>(Attribute::Term);
            if (match(term, "diff")) {
                auto atom = ast.copy();

                auto elements = atom.get<ASTVector>(Attribute::Elements);
                check_syntax(elements.size() == 1);
                Clingo::ASTv2::AST element = *elements.begin();
                auto tuple = element.get<ASTVector>(Attribute::Terms);
                check_syntax(tuple.size() == 1);
                auto condition = element.get<ASTVector>(Attribute::Condition);
                check_syntax(condition.empty());

                atom.set(Attribute::Term, tag_terms(term, in_literal ? "_b" : "_h"));

                return atom;
            }
        }
        return ast.transform_ast(*this);
    }
    bool in_literal{false};
};

// Checks if the given theory atom is shiftable.
struct SigMatcher {
    template <typename CStr>
    static bool visit(Clingo::Symbol const &f, CStr str) {
        return f.match(str, 0);
    }

    template <typename CStr, typename... CStrs>
    static bool visit(Clingo::Symbol const &f, CStr str, CStrs... strs) {
        return (f.match(str, 0) || visit(f, strs...));
    }

    template <typename CStr>
    static bool visit(Clingo::AST::Function const &f, CStr str) {
        return !f.external && f.arguments.empty() && (std::strcmp(f.name, str) == 0);
    }

    template <typename CStr, typename... CStrs>
    static bool visit(Clingo::AST::Function const &f, CStr str, CStrs... strs) {
        return !f.external && f.arguments.empty() && ((std::strcmp(f.name, str) == 0) || visit(f, strs...));
    }

    template <class T, typename CStr>
    static bool visit(T const &x, CStr str) {
        static_cast<void>(x);
        static_cast<void>(str);
        return false;
    }

    template <class T, typename CStr, typename... CStrs>
    static bool visit(T const &x, CStr str, CStrs... strs) {
        static_cast<void>(x);
        static_cast<void>(str);
        return visit(x, strs...);
    }
};

// Shifts constraints into rule heads.
struct TheoryShifter {
    static void visit(Clingo::AST::Rule &rule) {
        if (!rule.head.data.is<Clingo::AST::Literal>()) {
            return;
        }
        auto &head = rule.head.data.get<Clingo::AST::Literal>();
        if (!head.data.is<Clingo::AST::Boolean>() || head.data.get<Clingo::AST::Boolean>().value) {
            return;
        }
        auto it = rule.body.begin();
        auto ie = rule.body.end();
        auto jt = it;
        for (; it != ie; ++it) {
            if (it->data.is<Clingo::AST::TheoryAtom>()) {
                auto &atom = it->data.get<Clingo::AST::TheoryAtom>();
                SigMatcher matcher;
                if (atom.term.data.accept(matcher, "diff")) {
                    check_syntax(atom.guard.get() != nullptr);
                    if (it->sign != Clingo::AST::Sign::Negation) {
                        auto *guard = atom.guard.get();
                        guard->operator_name = negate_relation(guard->operator_name);
                    }
                    rule.head.location = it->location;
                    rule.head.data = std::move(atom);
                    for (++it; it != ie; ++it, ++jt) {
                        if (it != jt) {
                            std::iter_swap(it, jt);
                        }
                    }
                    break;
                }
            }
            if (it != jt) {
                std::iter_swap(it, jt);
            }
            ++jt;
        }
        rule.body.erase(jt, ie);
    }

    template <class T>
    static void visit(T &value) {
        static_cast<void>(value);
    }
};

struct TermTagger {
    static void visit(Clingo::AST::Function &term, char const *tag) {
        std::string name{"__"};
        name += term.name;
        name += tag;
        term.name = Clingo::add_string(name.c_str());
    }

    static void visit(Clingo::Symbol &term, char const *tag) {
        std::string name{"__"};
        name += term.name();
        name += tag;
        term = Clingo::Function(name.c_str(), {});
    }

    template <typename T>
    static void visit(T &value, char const *tag) {
        static_cast<void>(value);
        static_cast<void>(tag);
        throw_syntax_error();
    }
};

// Tags head and body atoms and ensures multiset semantics.
struct TheoryRewriter {
    static char const *tag(Clingo::AST::HeadLiteral const &lit) {
        static_cast<void>(lit);
        return "_h";
    }

    static char const *tag(Clingo::AST::BodyLiteral const &lit) {
        static_cast<void>(lit);
        return "_b";
    }

    // Mark head/body literals and ensure multiset semantics for theory atoms.
    template <class Lit>
    static void visit(Lit &node, Clingo::AST::TheoryAtom &atom) {
        SigMatcher matcher;
        if (atom.term.data.accept(matcher, "diff")) {
            TermTagger tagger;
            char const *msg = "Diff atoms must consist of a single term without conditions.";
            check_syntax(atom.elements.size() == 1, msg);
            auto &element = atom.elements.front();
            check_syntax(element.condition.empty() && element.tuple.size() == 1, msg);
            atom.term.data.accept(tagger, tag(node));
        }
    }
};

void transform(Clingo::AST::Statement &&stm, Clingo::StatementCallback const &cb, bool shift) {
    unpool(std::move(stm), [&](Clingo::AST::Statement &&unpooled) {
        if (shift) {
            TheoryShifter shifter;
            unpooled.data.accept(shifter);
        }
        TheoryRewriter tagger;
        transform_ast(tagger, unpooled);
        cb(std::move(unpooled));
    });
}

void transform(Clingo::ASTv2::AST const &ast, ASTCallback const &cb, bool shift) {
    for (auto &unpooled : ast.unpool()) {
        if (shift) {
            unpooled = shift_rule(unpooled);
        }
        cb(unpooled.transform_ast(TheoryRewriterR{}));
    }
}

bool match(Clingo::TheoryTerm const &term, char const *name, size_t arity) {
    return (term.type() == Clingo::TheoryTermType::Symbol &&
        std::strcmp(term.name(), name) == 0 &&
        arity == 0) ||
        (term.type() == Clingo::TheoryTermType::Function &&
        std::strcmp(term.name(), name) == 0 &&
        term.arguments().size() == arity);
}

} // namespace ClingoDL
