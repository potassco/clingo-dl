// {{{ MIT License

// Copyright 2017 Roland Kaminski

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

#include "propagator.h"
#include "clingoDL.h"
#include <clingo.hh>


using namespace Clingo;


class MyHandler : public SolveEventHandler
{
//    virtual bool on_model(Model const &model) override
//    {
//       auto v = model.symbols(ShowType::Theory);
//       if (std::find(v.begin(), v.end(), Function("dl",SymbolVector{Id("y",true),String("6")}, true))!=v.end()) {
//           // an answer where y is 6
//       }
//       return true;
//    }
};

template <typename T>
void solve(Stats &stats, Control &ctl, bool strict, bool propagate) {
    DifferenceLogicPropagator<T> p{stats, strict, propagate};
    ctl.register_propagator(p);
    ctl.ground({{"base", {}}});
    MyHandler h;
    for(auto m : ctl.solve(Clingo::SymbolicLiteralSpan(),&h)) {
        (void)(m);
        }
    std::cout << "\n";
}

class ClingoDLApp : public Clingo::ClingoApplication {
public:
    ClingoDLApp(Stats &stats)
        : stats_{stats} {}
    char const *program_name() const noexcept override { return "clingoDL"; }
    char const *version() const noexcept override { return CLINGODL_VERSION; }
    void main(Control &ctl, StringSpan files) override {
        ctl.add("base", {}, R"(#theory dl {
term{};
constant {
  + : 1, binary, left;
  - : 1, binary, left;
  * : 2, binary, left;
  / : 2, binary, left;
  - : 3, unary
};
diff_term {- : 1, binary, left};
&diff/0 : diff_term, {<=}, constant, any;
&show_assignment/0 : term, directive
}.)");
        for (auto &file : files) {
            ctl.load(file);
        }
        if (files.empty()) {
            ctl.load("-");
        }

        if (rdl_) {
            solve<double>(stats_, ctl, strict_, propagate_);
        }
        else {
            solve<int>(stats_, ctl, strict_, propagate_);
        }
    }

    void register_options(ClingoOptions &options) override {
        char const *group = "Clingo(DL) Options";
        options.add_flag(group, "propagate,p", "Enable propagation.", propagate_);
        options.add_flag(group, "rdl", "Enable support for real numbers.", rdl_);
        options.add_flag(group, "strict", "Enable strict mode.", strict_);
    }

    void validate_options() override {
        if (rdl_ && strict_) {
            // NOTE: could be implemented by introducing and epsilon
            throw std::runtime_error("real difference logic not available with strict semantics");
        }
    }

private:
    Stats &stats_;
    bool propagate_ = false;
    bool strict_ = false;
    bool rdl_ = false;
};

int main(int argc, char *argv[]) {
    Stats stats;
    int ret;
    {
        Timer t{stats.time_total};
        ClingoDLApp app{stats};
        ret = Clingo::clingo_main(app, {argv + 1, numeric_cast<size_t>(argc - 1)});
    }

    if (stats.dl_stats.size() > 0) {
        std::cout << "\n";
        std::cout << "propagator statistics:\n";
        std::cout << "  total: " << stats.time_total.count() << "s\n";
        std::cout << "    init: " << stats.time_init.count() << "s\n";
        int thread = 0;
        for (auto &stat : stats.dl_stats) {
            std::cout << "    total[" << thread << "]: ";
            std::cout << (stat.time_undo + stat.time_propagate).count() << "s\n";
            std::cout << "      propagate: " << stat.time_propagate.count() << "s\n";
            std::cout << "        dijkstra: " << stat.time_dijkstra.count() << "s\n";
            std::cout << "        true edges: " << stat.true_edges << "\n";
            std::cout << "        false edges: " << stat.false_edges << "\n";
            std::cout << "      undo: " << stat.time_undo.count() << "s\n";
            ++thread;
        }
    }
    return ret;
}
