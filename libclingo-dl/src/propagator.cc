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
template <typename T>
T evaluate_real(char const *);

template <>
int evaluate_real(char const *) {
    throw std::runtime_error("could not evaluate term: integer expected");
}

template <>
double evaluate_real(char const *name) {
    static const std::string chars = "\"";
    auto len = std::strlen(name);
    if (len < 2 || name[0] != '"' || name[len - 1] != '"') {
        throw std::runtime_error("could not evaluate term: real numbers have to be represented as strings");
    }
    char *parsed = nullptr;
    auto ret = std::strtod(name + 1, &parsed);
    if (parsed != name + len - 1) {
        throw std::runtime_error("could not evaluate term: not a valid real number");
    }
    return ret;
}


//template <typename T>
//T evaluate(Clingo::TheoryTerm term) {
//    switch (term.type()) {
//        case Clingo::TheoryTermType::Number: {
//            return term.number();
//        }
//        case Clingo::TheoryTermType::Symbol: {
//            return evaluate_real<T>(term.name());
//        }
//        case Clingo::TheoryTermType::Function: {
//            auto args = term.arguments();
//            if (args.size() == 2) {
//                return evaluate_binary(term.name(), evaluate<T>(args[0]), evaluate<T>(args[1]));
//            }
//            if (args.size() == 1) {
//                if (std::strcmp(term.name(), "-") == 0) {
//                    return -evaluate<T>(args[0]);
//                }
//                else {
//                    throw std::runtime_error("could not evaluate term: unknown unary operator");
//                }
//            }
//            // [[fallthrough]]
//        }
//        default: { throw std::runtime_error("could not evaluate term: only numeric terms with basic arithmetic operations are supported"); }
//    }
//}
//
//template int evaluate<int>(Clingo::TheoryTerm term);
//template double evaluate<double>(Clingo::TheoryTerm term);

int require_number(Clingo::Symbol sym) { return sym.type() == Clingo::SymbolType::Number ? sym.number() : throw std::runtime_error("could not evaluate term: artithmetic on non-integer"); }

//extern Clingo::Symbol evaluate_term(Clingo::TheoryTerm term) {
//    switch (term.type()) {
//        case Clingo::TheoryTermType::Number: {
//            return Clingo::Number(term.number());
//        }
//        case Clingo::TheoryTermType::Symbol: {
//            return Clingo::Id(term.name(), true);
//        }
//        case Clingo::TheoryTermType::Function: {
//            auto op = term.name();
//            std::vector<Clingo::Symbol> args;
//            for (auto arg : term.arguments()) {
//                args.emplace_back(evaluate_term(arg));
//            }
//            if (args.size() == 2) {
//                if (std::strcmp(op, "+") == 0) {
//                    return Clingo::Number(require_number(args[0]) + require_number(args[1]));
//                }
//                else if (std::strcmp(op, "-") == 0) {
//                    return Clingo::Number(require_number(args[0]) - require_number(args[1]));
//                }
//                else if (std::strcmp(op, "*") == 0) {
//                    return Clingo::Number(require_number(args[0]) * require_number(args[1]));
//                }
//                else if (std::strcmp(op, "/") == 0) {
//                    if (args[1] == Clingo::Number(0)) {
//                        throw std::runtime_error("could not evaluate term: division by zero");
//                    }
//                    return Clingo::Number(require_number(args[0]) / require_number(args[1]));
//                }
//            }
//            else if (args.size() == 1) {
//                if (std::strcmp(op, "-") == 0) {
//                    switch (args[0].type()) {
//                        case Clingo::SymbolType::Number: {
//                            return Clingo::Number(-args[0].number());
//                        }
//                        case Clingo::SymbolType::Function: {
//                            if (std::strcmp(args[0].name(), "") != 0) {
//                                return Clingo::Function(args[0].name(), args[0].arguments(), args[0].is_negative());
//                            }
//                            // [[fallthrough]]
//                        }
//                        default: { throw std::runtime_error("could not evaluate term: only numbers and functions can be inverted"); }
//                    }
//                }
//            }
//            return Clingo::Function(op, args);
//        }
//        case Clingo::TheoryTermType::Tuple: {
//            std::vector<Clingo::Symbol> args;
//            for (auto arg : term.arguments()) {
//                args.emplace_back(evaluate_term(arg));
//            }
//            return Clingo::Function("", args);
//        }
//        default: { throw std::runtime_error("could not evaluate term: sets and lists are not supported"); }
//    }
//}

bool match(Clingo::TheoryTerm const &term, char const *name, size_t arity) {
    return (term.type() == Clingo::TheoryTermType::Symbol &&
        std::strcmp(term.name(), name) == 0 &&
        arity == 0) ||
        (term.type() == Clingo::TheoryTermType::Function &&
        std::strcmp(term.name(), name) == 0 &&
        term.arguments().size() == arity);
}


