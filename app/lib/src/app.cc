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

#include <clingo-dl-app/app.hh>

#include <cmath>
#include <limits>

namespace ClingoDL {

#define CLINGODL_TRY try // NOLINT
#define CLINGODL_CATCH                                                                                                 \
    catch (...) {                                                                                                      \
        Clingo::Detail::store_error();                                                                                 \
        return false;                                                                                                  \
    }                                                                                                                  \
    return true // NOLINT

namespace {
auto add_(clingo_ast_t *stm, void *data) -> bool {
    auto *program = static_cast<clingo_program_t *>(data);
    return clingo_program_add(program, stm);
}
} // namespace

Optimizer::Optimizer(Clingo::Library const &lib, OptimizerConfig const &opt_cfg, Clingo::SolveEventHandler &handler,
                     Clingo::Theory const &theory)
    : lib_{lib}, opt_cfg_{opt_cfg}, handler_{handler}, theory_{theory} {}

void Optimizer::solve(Clingo::Control const &ctl) {
    theory_.rewrite(lib_, ctl,
                    // add a fixed bound
                    "#program __ub(s,b)."
                    "&diff { s-0 } <= b."
                    // retract previous bound
                    "#program __rb(b)."
                    "#external __sb(b). [release]"
                    // add a retractable bound
                    "#program __sb(s,b)."
                    "#external __sb(b). [true]"
                    "&diff { s-0 } <= b :- __sb(b).");
    if (opt_cfg_.has_initial) {
        upper_bound_ = opt_cfg_.initial;
    }
    for (;;) {
        prepare_(ctl);
        auto ret = ctl.solve(*this).get();
        if (ret.interrupted()) {
            break;
        }
        if (ret.unsatisfiable()) {
            if (search_bound_) {
                lower_bound_ = *search_bound_ + 1;
            }
            search_bound_ = upper_bound_;
            adjust_ = 1;
            if (!lower_bound_ || *lower_bound_ > *upper_bound_) {
                break;
            }
        }
    }
}

void Optimizer::do_stats(Clingo::Stats step, Clingo::Stats accu) {
    add_stats(step);
    add_stats(accu);
    handler_.stats(step, accu);
}

void Optimizer::add_stats(Clingo::Stats root) const {
    if (optimization || lower_bound_) {
        auto diff = root.map().insert("DifferenceLogic", Clingo::StatsType::map).map();
        if (optimization) {
            diff.insert("Optimization", Clingo::StatsType::value).value(*optimization);
        }
        if (lower_bound_) {
            diff.insert("Lower bound", Clingo::StatsType::value).value(*lower_bound_);
        }
    }
}

auto Optimizer::do_model(Clingo::Model model) -> bool {
    // update (upper) bound
    optimization = get_bound(model);
    upper_bound_ = *optimization - 1;

    // determine search bound
    double aux = *optimization - adjust_;
    if (lower_bound_ && aux <= *lower_bound_) {
        aux = *lower_bound_ + 1.0;
    }
    if (aux < std::numeric_limits<int_value_t>::min()) {
        aux = std::numeric_limits<int_value_t>::min();
    }
    search_bound_ = static_cast<int_value_t>(aux);

    // update (exponential) adujustment value
    adjust_ = adjust_ * opt_cfg_.factor;

    // pass model to theory
    handler_.model(model);
    return false;
}

auto Optimizer::get_bound(Clingo::Model model) -> int_value_t {
    // get bound
    bool found = false;
    auto ass = theory_.assignment(model.thread_id());
    if (opt_cfg_.index == 0) {
        if (auto index = ass.lookup(opt_cfg_.symbol)) {
            opt_cfg_.index = *index;
        } else {
            throw std::logic_error{"bound symbol not found"};
        }
    }
    auto value = ass.at(opt_cfg_.index);
    // NOTE: minimizinig real values would require an epsilon
    if (auto *ret = std::get_if<int>(&value.second); ret != nullptr) {
        return *ret;
    }
    throw std::runtime_error("only integer minimization is supported");
}

void Optimizer::prepare_(Clingo::Control const &ctl) {
    std::vector<Clingo::Part> parts = {};
    parts.reserve(3);
    if (upper_bound_ && upper_bound_ != upper_bound_last_) {
        upper_bound_last_ = upper_bound_;
        parts.push_back({"__ub", {opt_cfg_.symbol, Clingo::Number(*upper_bound_)}});
    }
    if (search_bound_ != search_bound_last_) {
        if (search_bound_last_) {
            parts.push_back({"__rb", {Clingo::Number(*search_bound_last_)}});
        }
        if (search_bound_ && search_bound_ != upper_bound_) {
            search_bound_last_ = search_bound_;
            parts.push_back({"__sb", {opt_cfg_.symbol, Clingo::Number(*search_bound_)}});
        } else {
            search_bound_last_ = Bound{};
        }
    }
    if (!parts.empty()) {
        ctl.ground(parts);
    }
}

} // namespace ClingoDL
