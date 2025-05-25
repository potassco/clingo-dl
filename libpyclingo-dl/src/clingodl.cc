#include <pybind11/pybind11.h>

namespace PyClingo {

void register_clingo(pybind11::module &m) {
    m.doc() = R"doc(
TODO:
)doc";
}

} // namespace PyClingo

PYBIND11_MODULE(clingo, m) { PyClingo::register_clingo(m); }
