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
    // Add variables to tuple to ensure multiset semantics.
    static void rewrite_tuple(Clingo::AST::TheoryAtomElement &element, int number) {
        check_syntax(element.tuple.size() == 1);
        auto vars_condition = collect_variables(element.condition.begin(), element.condition.end());
        for (auto const &name : collect_variables(element.tuple.begin(), element.tuple.end())) {
            vars_condition.erase(name);
        }
        vars_condition.erase("_");
        if (number >= 0) {
            element.tuple.push_back({element.tuple.front().location, Clingo::Number(number)});
        }
        for (auto const &name : vars_condition) {
            element.tuple.push_back({element.tuple.front().location, Clingo::AST::Variable{name}});
        }
    }

    // Add variables to tuples of elements to ensure multiset semantics.
    static void rewrite_tuples(Clingo::AST::TheoryAtom &atom) {
        int number = atom.elements.size() > 1 ? 0 : -1;
        for (auto &element : atom.elements) {
            rewrite_tuple(element, number++);
        }
    }

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
            int number = atom.elements.size() > 1 ? 0 : -1;
            for (auto &element : atom.elements) {
                rewrite_tuple(element, number++);
            }
        }

        if (atom.term.data.accept(matcher, "diff")) {
            TermTagger tagger;
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


bool match(Clingo::TheoryTerm const &term, char const *name, size_t arity) {
    return (term.type() == Clingo::TheoryTermType::Symbol &&
        std::strcmp(term.name(), name) == 0 &&
        arity == 0) ||
        (term.type() == Clingo::TheoryTermType::Function &&
        std::strcmp(term.name(), name) == 0 &&
        term.arguments().size() == arity);
}

} // namespace ClingoDL
