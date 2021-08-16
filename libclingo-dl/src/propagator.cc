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
bool match(Clingo::AST::Node const &ast, CStr str) {
    using namespace Clingo::AST;
    if (ast.type() == Type::SymbolicTerm) {
        auto sym = ast.get<Clingo::Symbol>(Attribute::Term);
        return sym.match(str, 0);
    }
    if (ast.type() == Type::Function) {
        if (ast.get<int>(Attribute::External) != 0) {
            return false;
        }
        if (!ast.get<NodeVector>(Attribute::Arguments).empty()) {
            return false;
        }
        auto const *name = ast.get<char const *>(Attribute::Name);
        return (std::strcmp(name, str) == 0);
    }
    return false;
}

Clingo::AST::Node shift_rule(Clingo::AST::Node ast) {
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
        if (atom.type() != Type::TheoryAtom || !match(atom.get<Node>(Attribute::Term), "diff")) {
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

Clingo::AST::Node tag_terms(Clingo::AST::Node &ast, char const *tag) {
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
        name += ret.get<char const*>(Attribute::Name);
        name += tag;
        ret.set(Attribute::Name, Clingo::add_string(name.c_str()));
        return ret;
    }
    return throw_syntax_error<Node>();
}

// Mark head/body theory atoms.
struct TheoryRewriter {
    Clingo::AST::Node operator()(Clingo::AST::Node const &ast) {
        using namespace Clingo::AST;
        if (ast.type() == Type::Literal) {
            in_literal = true;
            auto ret = ast.transform_ast(*this);
            in_literal = false;
            return ret;
        }
        if (ast.type() == Type::TheoryAtom) {
            auto term = ast.get<Node>(Attribute::Term);
            if (match(term, "diff")) {
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

void transform(Clingo::AST::Node const &ast, NodeCallback const &cb, bool shift) {
    for (auto &unpooled : ast.unpool()) {
        if (shift) {
            unpooled = shift_rule(unpooled);
        }
        cb(unpooled.transform_ast(TheoryRewriter{}));
    }
}

} // namespace ClingoDL
