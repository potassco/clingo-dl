#include <clingo-dl.h>

#include <pybind11/pybind11.h>

namespace PyClingoDL {

namespace {

auto create_theory() -> pybind11::object {
    return pybind11::capsule{reinterpret_cast<void *>(&clingodl_create), "clingo_theory_create"};
}

} // namespace

void register_clingodl(pybind11::module &m) {
    m.doc() = R"doc(TODO)doc";
    m.def("create_theory", create_theory, R"(TODO)");
}

} // namespace PyClingoDL

PYBIND11_MODULE(clingodl, m) { PyClingoDL::register_clingodl(m); }
