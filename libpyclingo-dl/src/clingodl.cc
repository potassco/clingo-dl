#include <clingo-dl.h>

#include <pybind11/pybind11.h>

namespace PyClingoDL {

namespace {

struct Free {
    void operator()(clingo_theory_t *theory) noexcept {
        if (theory->destroy != nullptr && theory->self != nullptr) {
            theory->destroy(theory->self);
            theory->self = nullptr;
        }
        delete theory;
    }
};

void destroy(void *data) noexcept {
    std::ignore = std::unique_ptr<clingo_theory_t, Free>(static_cast<clingo_theory_t *>(data));
}

pybind11::object make_theory(pybind11::handle lib) {
    using namespace std::string_view_literals;
    auto theory_module = pybind11::module::import("clingo.theory");
    auto lib_capsule = lib.attr("_capsule")().cast<pybind11::capsule>();
    if (lib_capsule.name() != "clingo_lib_t"sv) {
        throw std::invalid_argument("clingo_lib_t pointer expected");
    }
    clingo_lib_t *lib_ptr = static_cast<clingo_lib_t *>(lib_capsule.get_pointer());
    if (lib_ptr == nullptr) {
        throw std::invalid_argument{"library must not be null"};
    }
    auto theory_dl = std::unique_ptr<clingo_theory_t, Free>(new clingo_theory_t{});
    clingodl_create(lib_ptr, theory_dl.get());
    auto theory_capsule = pybind11::capsule{theory_dl.release(), "clingo_theory_t", destroy};
    return theory_module.attr("Theory")(std::move(theory_capsule));
}

} // namespace

void register_clingodl(pybind11::module &m) {
    m.doc() = R"doc(TODO)doc";
    m.def("theory", make_theory, pybind11::arg("lib"), R"(TODO)");
}

} // namespace PyClingoDL

PYBIND11_MODULE(clingodl, m) { PyClingoDL::register_clingodl(m); }
