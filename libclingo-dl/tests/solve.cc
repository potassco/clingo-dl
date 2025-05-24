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

#include <clingo-dl-app/app.hh>
#include <clingo-dl.h>
#include <clingo-dl/propagator.hh>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <iostream>

namespace ClingoDL {

using namespace std::string_view_literals;

namespace {

// TODO: move to cxx api
class TheoryAssignment {
  public:
    struct sentinel {};
    class iterator {
      public:
        using iterator_category = std::input_iterator_tag;
        using value_type = std::pair<Clingo::Symbol, std::variant<int, double, Clingo::Symbol>>;
        using difference_type = std::ptrdiff_t;
        using pointer = value_type *;
        using reference = value_type &;

        iterator(clingo_theory_t *theory, uint32_t thread_id) : theory_{theory}, thread_id_{thread_id} { advance(); }

        iterator() = default;

        reference operator*() const { return current_; }

        iterator &operator++() {
            advance();
            return *this;
        }

        iterator operator++(int) { return std::exchange(*this, ++iterator(*this)); }

        bool operator==([[maybe_unused]] sentinel const &other) const { return !has_value_; }

      private:
        void advance() {
            Clingo::Detail::handle_error(
                theory_->assignment_next(theory_->self, thread_id_, &init_, &index_, &has_value_));
            if (has_value_) {
                clingo_symbol_t symbol;
                clingo_theory_value_t value;
                Clingo::Detail::handle_error(
                    theory_->assignment_get_value(theory_->self, thread_id_, index_, &symbol, &value, nullptr));
                current_.first = Clingo::Symbol{symbol, false};
                switch (value.type) {
                    case clingo_theory_value_type_int: {
                        current_.second = value.int_number;
                        break;
                    }
                    case clingo_theory_value_type_double: {
                        current_.second = value.double_number;
                        break;
                    }
                    case clingo_theory_value_type_symbol: {
                        current_.second = Clingo::Symbol{value.symbol, false};
                        break;
                    }
                }
            }
        }

        mutable value_type current_;
        clingo_theory_t *theory_ = nullptr;
        size_t index_ = 0;
        uint32_t thread_id_ = 0;
        bool init_ = true;
        bool has_value_ = true;
    };
    static_assert(std::input_iterator<iterator>);
    static_assert(std::sentinel_for<sentinel, iterator>);

    explicit TheoryAssignment(clingo_theory_t *theory, uint32_t thread_id) : theory_{theory}, thread_id_{thread_id} {}

    iterator begin() const { return iterator(theory_, thread_id_); }
    sentinel end() const { return sentinel(); }

  private:
    clingo_theory_t *theory_;
    uint32_t thread_id_;
};

// TODO: move to cxx api
class Theory {
  public:
    Theory(Clingo::Library lib) { clingodl_create(c_cast(lib), &theory_); }
    Theory(Theory &&other) = delete;
    ~Theory() {
        if (theory_.destroy != nullptr) {
            theory_.destroy(theory_.self);
        }
    }

    void prepare(Clingo::Control const &ctl) {
        Clingo::Detail::handle_error(theory_.prepare(theory_.self, c_cast(ctl)));
    }

    template <class F> void rewrite(Clingo::AST::Node stm, F fun) {
        constexpr auto add = [](clingo_ast_t *stm, void *data) -> bool {
            CLINGO_TRY {
                auto *fun = static_cast<F *>(data);
                std::invoke<F &>(*fun, Clingo::AST::Node{stm, true});
            }
            CLINGO_CATCH;
        };
        Clingo::Detail::handle_error(theory_.rewrite_ast(theory_.self, c_cast(stm), +add, static_cast<void *>(&fun)));
    }

    void rewrite(Clingo::Library const &lib, Clingo::Control const &ctl, std::string_view str) {
        auto scanner = Clingo::AST::Scanner{lib, str};
        auto prg = Clingo::AST::Program{lib};
        for (auto &stm : scanner) {
            rewrite(std::move(stm), [&](Clingo::AST::Node node) { prg.add(std::move(node)); });
        }
        ctl.join(prg);
    }

    void rewrite(Clingo::Library const &lib, Clingo::Control const &ctl, Clingo::StringSpan files) {
        auto scanner = Clingo::AST::Scanner{lib, files};
        auto prg = Clingo::AST::Program{lib};
        for (auto &stm : scanner) {
            rewrite(std::move(stm), [&](Clingo::AST::Node node) { prg.add(std::move(node)); });
        }
        ctl.join(prg);
    }

    auto assignment(uint32_t thread_id) -> TheoryAssignment { return TheoryAssignment{&theory_, thread_id}; }

    void stats(Clingo::Stats step, [[maybe_unused]] Clingo::Stats accu) {
        Clingo::Detail::handle_error(theory_.on_stats(theory_.self, c_cast(step)));
    }

    void register_theory(Clingo::Control const &ctl) {
        Clingo::Detail::handle_error(theory_.register_theory(theory_.self, c_cast(ctl)));
    }

    void configure(std::string_view key, std::string_view value) {
        Clingo::Detail::handle_error(
            theory_.configure(theory_.self, key.data(), key.size(), value.data(), value.size()));
    }

  private:
    clingo_theory_t theory_;
};

struct Fixture {
    //! A DL assignment.
    using A = std::pair<Clingo::Symbol, double>;
    //! A vector of DL assignments.
    using AV = std::vector<A>;
    //! A vector of symbols.
    using SV = std::vector<Clingo::Symbol>;
    //! A solution in form of a pair of DL assignments and symbols.
    using SP = std::pair<AV, SV>;
    //! A vector solutions.
    using RV = std::vector<SP>;

    //! Encoding for the job shop problem.
    static constexpr char const *ENC = R"(
task(T):-executionTime(T,_,_).
machine(M):-executionTime(_,M,_).

% decide which operation first
{seq((T,ST1),(T,ST2),Time1)} :- assign(T,ST1,M1), assign(T,ST2,M2), ST1<ST2, executionTime(T,ST1,Time1), executionTime(T,ST2,Time2).
seq((T,ST2),(T,ST1),Time2) :- assign(T,ST1,M1), assign(T,ST2,M2), ST1<ST2, executionTime(T,ST1,Time1), executionTime(T,ST2,Time2), not seq((T,ST1),(T,ST2),Time1).

% decide which task first on machine
{seq((T1,ST1),(T2,ST2),Time1)} :- assign(T1,ST1,M), assign(T2,ST2,M), T1<T2, executionTime(T1,ST1,Time1), executionTime(T2,ST2,Time2).
seq((T2,ST2),(T1,ST1),Time2) :- assign(T1,ST1,M), assign(T2,ST2,M), T1<T2, executionTime(T1,ST1,Time1), executionTime(T2,ST2,Time2), not seq((T1,ST1),(T2,ST2),Time1).

&diff{T1-T2}<= -Time:-seq(T1,T2,Time).

&diff{0-(T,M)} <= 0 :- task(T), machine(M), bound(B).
&diff{(T,M)-0} <= B :- task(T), machine(M), bound(B).

#show seq/3.

executionTime(1,1,54).
executionTime(1,2,34).
executionTime(1,3,61).
executionTime(1,4,2).
executionTime(2,1,9).
executionTime(2,2,15).
executionTime(2,3,89).
executionTime(2,4,70).
executionTime(3,1,38).
executionTime(3,2,19).
executionTime(3,3,28).
executionTime(3,4,87).
assign(1,1,3).
assign(1,2,1).
assign(1,3,4).
assign(1,4,2).
assign(2,1,4).
assign(2,2,1).
assign(2,3,2).
assign(2,4,3).
assign(3,1,1).
assign(3,2,2).
assign(3,3,3).
assign(3,4,4).
bound(104).
)";

    //! Create a symbol for sequence atoms of task/machine pairs.
    auto seq(int a, int b, int c, int d, int e) -> Clingo::Symbol {
        return Clingo::Function(lib, "seq",
                                {Clingo::Tuple(lib, {Clingo::Number(a), Clingo::Number(b)}),
                                 Clingo::Tuple(lib, {Clingo::Number(c), Clingo::Number(d)}), Clingo::Number(e)});
    }

    //! A DL assignment for task/machine pairs.
    auto ass(int a, int b, int c) -> A { return A(Clingo::Tuple(lib, {Clingo::Number(a), Clingo::Number(b)}), c); }

    auto sols() -> RV {
        return {SP{{
                       ass(1, 1, 100), ass(1, 2, 0), ass(1, 3, 34), ass(1, 4, 95), // NOLINT
                       ass(2, 1, 95), ass(2, 2, 72), ass(2, 3, 104), ass(2, 4, 0), // NOLINT
                       ass(3, 1, 34), ass(3, 2, 0), ass(3, 3, 72), ass(3, 4, 104)  // NOLINT
                   },
                   {
                       seq(1, 2, 1, 1, 34), seq(1, 2, 1, 3, 34), seq(1, 2, 1, 4, 34), // NOLINT
                       seq(1, 2, 2, 2, 34), seq(1, 2, 3, 1, 34), seq(1, 3, 1, 1, 61), // NOLINT
                       seq(1, 3, 1, 4, 61), seq(1, 3, 2, 1, 61), seq(1, 3, 3, 4, 61), // NOLINT
                       seq(1, 4, 1, 1, 2),  seq(1, 4, 2, 3, 2),                       // NOLINT
                       seq(2, 1, 2, 3, 9),  seq(2, 1, 3, 4, 9),  seq(2, 2, 2, 1, 15), // NOLINT
                       seq(2, 2, 2, 3, 15), seq(2, 4, 1, 1, 70), seq(2, 4, 2, 1, 70), // NOLINT
                       seq(2, 4, 2, 2, 70), seq(2, 4, 2, 3, 70), seq(2, 4, 3, 3, 70), // NOLINT
                       seq(3, 1, 2, 2, 38), seq(3, 1, 3, 3, 38), seq(3, 1, 3, 4, 38), // NOLINT
                       seq(3, 2, 1, 4, 19), seq(3, 2, 2, 3, 19), seq(3, 2, 3, 1, 19), // NOLINT
                       seq(3, 2, 3, 3, 19), seq(3, 2, 3, 4, 19), seq(3, 3, 1, 1, 28), // NOLINT
                       seq(3, 3, 3, 4, 28),                                           // NOLINT
                   }},
                SP{{
                       ass(1, 1, 104), ass(1, 2, 70), ass(1, 3, 9), ass(1, 4, 0), // NOLINT
                       ass(2, 1, 0), ass(2, 2, 9), ass(2, 3, 98), ass(2, 4, 28),  // NOLINT
                       ass(3, 1, 28), ass(3, 2, 66), ass(3, 3, 0), ass(3, 4, 85)  // NOLINT
                   },
                   {
                       seq(1, 2, 1, 1, 34), seq(1, 3, 1, 1, 61), seq(1, 3, 1, 2, 61), // NOLINT
                       seq(1, 3, 3, 4, 61), seq(1, 4, 1, 1, 2),  seq(1, 4, 1, 2, 2),  // NOLINT
                       seq(1, 4, 1, 3, 2),  seq(1, 4, 2, 3, 2),  seq(1, 4, 3, 2, 2),  // NOLINT
                       seq(2, 1, 1, 3, 9),  seq(2, 1, 2, 2, 9),  seq(2, 1, 2, 3, 9),  // NOLINT
                       seq(2, 1, 2, 4, 9),  seq(2, 1, 3, 4, 9),  seq(2, 2, 1, 2, 15), // NOLINT
                       seq(2, 2, 2, 3, 15), seq(2, 2, 2, 4, 15), seq(2, 2, 3, 1, 15), // NOLINT
                       seq(2, 4, 1, 1, 70), seq(2, 4, 2, 3, 70), seq(3, 1, 1, 2, 38), // NOLINT
                       seq(3, 1, 3, 2, 38), seq(3, 1, 3, 4, 38), seq(3, 2, 2, 3, 19), // NOLINT
                       seq(3, 2, 3, 4, 19), seq(3, 3, 1, 1, 28), seq(3, 3, 2, 4, 28), // NOLINT
                       seq(3, 3, 3, 1, 28), seq(3, 3, 3, 2, 28), seq(3, 3, 3, 4, 28), // NOLINT
                   }}};
    }
    //! Solutions to the task assignment problem.
    //! A handler to gather statistics in a DL theory.
    class Handler : public Clingo::SolveEventHandler {
      public:
        Handler(Theory &theory) : theory_{&theory} {}
        //! Add theory specific statistics.
        void do_stats(Clingo::Stats step, Clingo::Stats accu) override { theory_->stats(step, accu); }

      private:
        Theory *theory_; //!< The DL theory.
    };

    //! Solve a given DL problem returning all models.
    auto solve(Clingo::Control &ctl) -> RV {
        using namespace Clingo;
        Handler h{theory};
        RV result;
        for (auto &&m : ctl.solve(h, {}, SolveFlags::yield)) {
            result.emplace_back();
            auto &sol = result.back().first;
            auto &sol_bool = result.back().second;
            for (auto &[key, value] : theory.assignment(m.thread_id())) {
                if (auto *num = std::get_if<int>(&value)) {
                    sol.emplace_back(key, *num);
                } else if (auto *num = std::get_if<double>(&value)) {
                    sol.emplace_back(key, *num);
                } else {
                    REQUIRE(false);
                }
            }
            std::ranges::sort(sol);
            for (auto s : m.symbols()) {
                sol_bool.emplace_back(s);
            }
            std::ranges::sort(sol_bool);
        }
        std::ranges::sort(result);
        return result;
    }

    void print(RV const &result) {
        for (auto const &[ass, syms] : result) {
            std::cerr << "solution:\n";
            std::cerr << "  symbols:";
            for (auto sym : syms) {
                std::cerr << " " << sym;
            }
            std::cerr << std::endl;
            std::cerr << "  assignment:";
            for (auto [sym, val] : ass) {
                std::cerr << " " << sym << "=" << val;
            }
            std::cerr << std::endl;
        }
    }

    Clingo::Library lib;
    Theory theory{lib};
    Clingo::Control ctl{lib, {"0"}};
    Clingo::Symbol sym_a = Function(lib, "a");
    Clingo::Symbol sym_b = Function(lib, "b");
    Clingo::Symbol sym_c = Function(lib, "c");
    Clingo::Symbol sym_d = Function(lib, "d");
    Clingo::Symbol sym_e = Function(lib, "e");
    Clingo::Symbol sym_f = Tuple(lib, {Function(lib, "f"), Function(lib, "f")});
};

} // namespace

TEST_CASE_METHOD(Fixture, "solving base", "[clingo]") { // NOLINT
    theory.register_theory(ctl);
    theory.rewrite(lib, ctl,
                   "#program base.\n"
                   "1 { a; b } 1. &diff { a - b } <= 3.\n"
                   "&diff { 0 - a } <= -5 :- a.\n"
                   "&diff { 0 - b } <= -7 :- b.\n");
    ctl.ground();
    theory.prepare(ctl);
    auto result = solve(ctl);
    REQUIRE(result == (RV{{{{sym_a, 0}, {sym_b, 7}}, {sym_b}}, {{{sym_a, 5}, {sym_b, 2}}, {sym_a}}}));

    theory.rewrite(lib, ctl,
                   "#program ext.\n"
                   "&diff { a - 0 } <= 4.\n");
    ctl.ground({{"ext", {}}});
    theory.prepare(ctl);
    result = solve(ctl);
    REQUIRE(result == (RV{{{{sym_a, 0}, {sym_b, 7}}, {sym_b}}}));
}

TEST_CASE_METHOD(Fixture, "solving not_equal", "[clingo]") {
    theory.register_theory(ctl);
    theory.rewrite(lib, ctl,
                   "#program base.\n"
                   "{ a }. &diff { b } != 5 :- not a.\n");
    ctl.ground();
    theory.prepare(ctl);
    auto result = solve(ctl);
    REQUIRE(result == (RV{{{}, {sym_a}}, {{{sym_b, 0}}, {}}, {{{sym_b, 6}}, {}}}));
}

TEST_CASE_METHOD(Fixture, "solving cc", "[clingo]") {
    theory.register_theory(ctl);
    theory.rewrite(lib, ctl,
                   "#program base.\n"
                   "&diff { 0 - a } <= -5.\n"
                   "&diff { 0 - b } <= -10.\n");
    ctl.ground();
    theory.prepare(ctl);
    auto result = solve(ctl);
    REQUIRE(result == (RV{{{{sym_a, 5}, {sym_b, 10}}, {}}}));
    REQUIRE(ctl.stats()["user_step"]["DifferenceLogic"]["CCs"].value() == 2);

    theory.rewrite(lib, ctl,
                   "#program ext.\n"
                   "&diff { b - a } <= 3.\n");
    ctl.ground({{"ext", {}}});
    theory.prepare(ctl);
    result = solve(ctl);
    REQUIRE(result == (RV{{{{sym_a, 7}, {sym_b, 10}}, {}}}));
    REQUIRE(ctl.stats()["user_step"]["DifferenceLogic"]["CCs"].value() == 1);
}

TEST_CASE_METHOD(Fixture, "solving configure", "[clingo]") {
    theory.configure("propagate", "full");
    theory.register_theory(ctl);
    theory.rewrite(lib, ctl,
                   "#program base.\n"
                   "&diff { a - 0 } <= 0.\n"
                   "a :- &diff { a - 0 } <=  0.\n"
                   "b :- &diff { 0 - a } <= -1.\n"
                   "c :- &diff { a } <= -1.\n");
    ctl.ground();
    theory.prepare(ctl);
    auto result = solve(ctl);
    REQUIRE(result == (RV{{{{sym_a, -1}}, {sym_a, sym_c}}, {{{sym_a, 0}}, {sym_a}}}));
}

TEST_CASE_METHOD(Fixture, "solving rdl", "[clingo]") {
    theory.configure("rdl", "yes");
    theory.register_theory(ctl);
    theory.rewrite(lib, ctl,
                   "#program base.\n"
                   "&diff { a } >= \"0.5\" * 3.\n");
    ctl.ground();
    theory.prepare(ctl);
    auto result = solve(ctl);
    REQUIRE(result == (RV{{{{sym_a, 1.5}}, {}}})); // NOLINT
}

TEST_CASE_METHOD(Fixture, "solving parse", "[clingo]") {
    theory.register_theory(ctl);
    theory.rewrite(lib, ctl,
                   "#program base.\n"
                   "&diff { p( 1 + 2 ) - q( 3 * 4 - 7 ) } <= 3 - 9.\n");
    ctl.ground();
    theory.prepare(ctl);
    auto result = solve(ctl);
    auto p = Clingo::parse_term(lib, "p(3)");
    auto q = Clingo::parse_term(lib, "q(5)");
    REQUIRE(result == (RV{{{{p, 0}, {q, 6}}, {}}}));
}

TEST_CASE_METHOD(Fixture, "solving normalize", "[clingo]") {
    theory.register_theory(ctl);
    theory.rewrite(lib, ctl,
                   "#program base.\n"
                   "&diff { a } = b.\n"
                   "&diff { 5 } >= 0.\n"
                   "&diff { b } > c.\n"
                   "&diff { c } >= d + 1.\n"
                   "&diff { e } != (f,f).\n");
    ctl.ground();
    theory.prepare(ctl);
    auto result = solve(ctl);
    REQUIRE(result == RV{
                          {{{sym_f, 0}, {sym_a, 2}, {sym_b, 2}, {sym_c, 1}, {sym_d, 0}, {sym_e, 1}}, {}},
                          {{{sym_f, 1}, {sym_a, 2}, {sym_b, 2}, {sym_c, 1}, {sym_d, 0}, {sym_e, 0}}, {}},
                      });
}

TEST_CASE_METHOD(Fixture, "solving empty", "[clingo]") {
    theory.register_theory(ctl);
    theory.rewrite(lib, ctl,
                   "#program base.\n"
                   "a :- &diff { a - a } <= 5.\n"
                   "{ b }.\n"
                   "&diff { 0 } < -4 :- b.\n");
    ctl.ground();
    theory.prepare(ctl);
    auto result = solve(ctl);
    REQUIRE(result == (RV{{{}, {sym_a}}}));
    REQUIRE(ctl.stats()["solving"]["solvers"]["choices"].value() == 0);
}

TEST_CASE_METHOD(Fixture, "solving symbols", "[clingo]") {
    theory.register_theory(ctl);
    theory.rewrite(lib, ctl,
                   "#program base.\n"
                   "&diff{ (\"foo\\\\\\nbar\\\"foo\",123) - 0 } <= 17.\n");
    ctl.ground();
    theory.prepare(ctl);
    auto result = solve(ctl);
    REQUIRE(result == (RV{{{{Tuple(lib, {String(lib, "foo\\\nbar\"foo"), Clingo::Number(123)}), 0}}, {}}}));
}

TEST_CASE_METHOD(Fixture, "solving task-assignment", "[clingo]") {
    auto mode = GENERATE("no", "inverse", "partial", "partial+", "zero", "full");
    auto mutex = GENERATE("0", "10,100");
    auto sort_edges = GENERATE("no", "weight", "potential");
    theory.configure("propagate", mode);
    theory.configure("add-mutexes", mutex);
    theory.configure("sort-edges", sort_edges);
    theory.register_theory(ctl);
    theory.rewrite(lib, ctl, ENC);
    ctl.ground();
    theory.prepare(ctl);
    REQUIRE(solve(ctl) == sols());
}

} // namespace ClingoDL
