// {{{ MIT License

// Copyright Roland Kaminski, Philipp Wanko, and Max Ostrowski

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

// }}}

#ifndef CLINGODL_APP_HH
#define CLINGODL_APP_HH

#include <clingo-dl.h>
#include <clingo.hh>
#include <optional>

namespace ClingoDL {

//! Type used for integer values.
using int_value_t = int;

//! Helper class to rewrite logic programs to use with the clingo DL theory.
class Rewriter {
  public:
    Rewriter(clingodl_theory_t *theory, clingo_program_builder_t *builder);
    //! Rewrite the given files.
    void rewrite(Clingo::Control &ctl, Clingo::StringSpan files);
    //! Rewrite the given program.
    void rewrite(Clingo::Control &ctl, char const *str);

  private:
    //! C callback to add a statement using the builder.
    static auto add_(clingo_ast_t *stm, void *data) -> bool;

    //! C callback to rewrite a statement and add it via the builder.
    static auto rewrite_(clingo_ast_t *stm, void *data) -> bool;

    clingodl_theory_t *theory_;         //!< A theory handle to rewrite statements.
    clingo_program_builder_t *builder_; //!< The builder to add rewritten statements to.
};

//! The configuration of the optimization algorithm.
struct OptimizerConfig {
    Clingo::Symbol symbol;   //!< The variable to minimize.
    double factor{1.0};      //!< The factor by which to increase optimization step size.
    mutable size_t index{0}; //!< The index of the symbol to minimize.
    int_value_t initial{0};  //!< The initial value of the bound.
    bool has_initial{false}; //!< Whether the bound is set initially.
    bool active{false};      //!< Whether optimization is enabled.
};

//! Class providing a configurable optimization algorithm.
class Optimizer : private Clingo::SolveEventHandler {
  private:
    //! Type used to store bounds.
    using Bound = std::optional<int_value_t>;
    using EventHandler = Clingo::SolveEventHandler;

  public:
    Optimizer(OptimizerConfig const &opt_cfg, EventHandler &handler, clingodl_theory_t *theory);
    //! Run the optimization algorithm.
    //!
    //! \note
    //! With an API extension to implement a custom enumerator, one could
    //! implement this more nicely. Right now, this implementation is
    //! restricted to the application.
    void solve(Clingo::Control &ctl);

  private:
    //! Function to add DL specific statistics.
    void on_statistics(Clingo::UserStatistics step, Clingo::UserStatistics accu) override;
    //! Add information about bounds to the given root statistics object and
    //! pass call the theory specific handler.
    void add_stats(Clingo::UserStatistics root) const;
    //! Function to extract the current bound and pass the model to the theory.
    auto on_model(Clingo::Model &model) -> bool override;
    //! Extract the bound from the given model.
    auto get_bound(Clingo::Model &model) -> int_value_t;
    //! Prepare the program for solving.
    //!
    //! This adds constraints to enforce the current upper or search bound as
    //! well as removes no longer required bounds.
    void prepare_(Clingo::Control &ctl);

    OptimizerConfig const &opt_cfg_; //!< Configuration of the optimization algorithm.
    EventHandler &handler_;          //!< Theory specific solve event handler.
    clingodl_theory_t *theory_;      //!< The underlying DL theory.
    Bound search_bound_;             //!< The current (volatile) search bound.
    Bound search_bound_last_;        //!< The previous search bound.
    Bound lower_bound_;              //!< The current lower bound (from UNSAT results).
    Bound upper_bound_;              //!< The current upper bound (from SAT results).
    Bound upper_bound_last_;         //!< The last value of the upper bound.
    Bound optimization;              //!< The current optimization value (last model found).
    double adjust_{1};               //!< Amount by which to decrease seach bound.
};

} // namespace ClingoDL

#endif // CLINGODL_APP_HH
