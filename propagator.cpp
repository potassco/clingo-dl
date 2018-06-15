#include "propagator.h"


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

template <typename T>
T evaluate(Clingo::TheoryTerm term) {
    switch (term.type()) {
        case Clingo::TheoryTermType::Number: {
            return term.number();
        }
        case Clingo::TheoryTermType::Symbol: {
            return evaluate_real<T>(term.name());
        }
        case Clingo::TheoryTermType::Function: {
            auto args = term.arguments();
            if (args.size() == 2) {
                return evaluate_binary(term.name(), evaluate<T>(args[0]), evaluate<T>(args[1]));
            }
            if (args.size() == 1) {
                if (std::strcmp(term.name(), "-") == 0) {
                    return -evaluate<T>(args[0]);
                }
                else {
                    throw std::runtime_error("could not evaluate term: unknown unary operator");
                }
            }
            // [[fallthrough]]
        }
        default: { throw std::runtime_error("could not evaluate term: only numeric terms with basic arithmetic operations are supported"); }
    }
}

template int evaluate<int>(Clingo::TheoryTerm term);
template double evaluate<double>(Clingo::TheoryTerm term);




int require_number(Clingo::Symbol sym) { return sym.type() == Clingo::SymbolType::Number ? sym.number() : throw std::runtime_error("could not evaluate term: artithmetic on non-integer"); }


template <>
int get_weight(TheoryAtom const &atom) {
    return evaluate<int>(atom.guard().second);
}
template <>
double get_weight(TheoryAtom const &atom) {
    return evaluate<double>(atom.guard().second);
}

extern Clingo::Symbol evaluate_term(Clingo::TheoryTerm term) {
    switch (term.type()) {
        case Clingo::TheoryTermType::Number: {
            return Clingo::Number(term.number());
        }
        case Clingo::TheoryTermType::Symbol: {
            return Clingo::Id(term.name(), true);
        }
        case Clingo::TheoryTermType::Function: {
            auto op = term.name();
            std::vector<Clingo::Symbol> args;
            for (auto arg : term.arguments()) {
                args.emplace_back(evaluate_term(arg));
            }
            if (args.size() == 2) {
                if (std::strcmp(op, "+") == 0) {
                    return Clingo::Number(require_number(args[0]) + require_number(args[1]));
                }
                else if (std::strcmp(op, "-") == 0) {
                    return Clingo::Number(require_number(args[0]) - require_number(args[1]));
                }
                else if (std::strcmp(op, "*") == 0) {
                    return Clingo::Number(require_number(args[0]) * require_number(args[1]));
                }
                else if (std::strcmp(op, "/") == 0) {
                    if (args[1] == Clingo::Number(0)) {
                        throw std::runtime_error("could not evaluate term: division by zero");
                    }
                    return Clingo::Number(require_number(args[0]) / require_number(args[1]));
                }
            }
            else if (args.size() == 1) {
                if (std::strcmp(op, "-") == 0) {
                    switch (args[0].type()) {
                        case Clingo::SymbolType::Number: {
                            return Clingo::Number(-args[0].number());
                        }
                        case Clingo::SymbolType::Function: {
                            if (std::strcmp(args[0].name(), "") != 0) {
                                return Clingo::Function(args[0].name(), args[0].arguments(), args[0].is_negative());
                            }
                            // [[fallthrough]]
                        }
                        default: { throw std::runtime_error("could not evaluate term: only numbers and functions can be inverted"); }
                    }
                }
            }
            return Clingo::Function(op, args);
        }
        case Clingo::TheoryTermType::Tuple: {
            std::vector<Clingo::Symbol> args;
            for (auto arg : term.arguments()) {
                args.emplace_back(evaluate_term(arg));
            }
            return Clingo::Function("", args);
        }
        default: { throw std::runtime_error("could not evaluate term: sets and lists are not supported"); }
    }
}

