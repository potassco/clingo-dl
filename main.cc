// {{{ MIT License

// Copyright 2018 Roland Kaminski, Philipp Wanko, Max Ostrowski

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

#include <clingo.hh>
#include <clingo-dl.h>

using namespace Clingo;

#define CLINGO_CALL(x) Clingo::Detail::handle_error(x)

class ClingoDLApp : public Clingo::Application, private SolveEventHandler {
public:
    ClingoDLApp() {
        CLINGO_CALL(clingodl_create_propagator(&prop_));
    }
    ~ClingoDLApp() { clingodl_destroy_propagator(prop_); }
    char const *program_name() const noexcept override { return "clingo-dl"; }
    char const *version() const noexcept override { return CLINGODL_VERSION; }
    bool on_model(Model &model) override {
        CLINGO_CALL(clingodl_on_model(prop_, model.to_c()));
        return true;
    }

    void on_statistics(UserStatistics step, UserStatistics accu) override {
        CLINGO_CALL(clingodl_on_statistics(prop_, step.to_c(), accu.to_c()));
    }
    void main(Control &ctl, StringSpan files) override {
        CLINGO_CALL(clingodl_register_propagator(prop_, ctl.to_c()));
        for (auto &file : files) {
            ctl.load(file);
        }
        if (files.empty()) {
            ctl.load("-");
        }

        ctl.ground({{"base", {}}});
        ctl.solve(Clingo::SymbolicLiteralSpan{}, this, false, false).get();
    }

    void register_options(ClingoOptions &options) override {
        CLINGO_CALL(clingodl_register_options(prop_, options.to_c()));
    }

    void validate_options() override {
        CLINGO_CALL(clingodl_validate_options(prop_));
    }
    clingodl_propagator_t *prop_;
};

int main(int argc, char *argv[]) {
    ClingoDLApp app;
    return Clingo::clingo_main(app, {argv + 1, static_cast<size_t>(argc - 1)});
}


#undef CLINGO_CALL
