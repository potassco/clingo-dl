#include <Python.h>

#include <clingo-dl.h>

#include <array>
#include <utility>

namespace {

template <class T, auto Destroy> class Guard {
  public:
    constexpr explicit Guard(T v = {}) noexcept : value_{v} {}

    Guard(Guard const &other) = delete;
    auto operator=(Guard const &other) -> Guard & = delete;
    constexpr Guard(Guard &&other) noexcept : value_{std::exchange(other.value_, T{})} {}
    constexpr auto operator=(Guard &&other) noexcept -> Guard & {
        if (this != &other) {
            destroy_();
            value_ = std::exchange(other.value_, T{});
        }
        return *this;
    }

    ~Guard() noexcept { destroy_(); }

    [[nodiscard]] constexpr auto get() noexcept -> T & { return value_; }

    auto release() noexcept -> T {
        auto v = value_;
        value_ = T{};
        return v;
    }

  private:
    void destroy_() {
        if (value_ != T{}) {
            Destroy(value_);
        }
    }
    T value_{};
};

using PyGuard = Guard<PyObject *, Py_DECREF>;

auto steal(PyObject *obj) -> PyGuard { return PyGuard(obj); }

auto borrow(PyObject *obj) -> PyGuard {
    Py_XINCREF(obj);
    return PyGuard(obj);
}

auto create_theory([[maybe_unused]] PyObject *module, [[maybe_unused]] PyObject *arg) -> PyObject * {
    // NOLINTNEXTLINE
    return PyCapsule_New(reinterpret_cast<void *>(&clingodl_create), "clingo_theory_create", nullptr);
}

auto run_app([[maybe_unused]] PyObject *module, [[maybe_unused]] PyObject *arg) -> PyObject * {
    static const char *source = R"py(
import sys
from typing import Sequence

from clingo import ast
from clingo.app import App, AppOptions, clingo_main
from clingo.control import Control, ControlMode
from clingo.core import Library
from clingo.script import enable_python
from clingo.theory import Theory
from clingodl import create_theory


class ClingoDLApp(App):
    def __init__(self, lib: Library):
        theory = Theory(lib, create_theory())
        major, minor, revision = theory.version
        super().__init__(theory.name, f"{major}.{minor}.{revision}", theory.commit_hash)
        self._lib = lib
        self._theory = theory
        globals()["THEORY"] = theory

    def main(self, control: Control, files: Sequence[str]) -> None:
        self._theory.register(control)
        self._theory.rewrite_files(self._lib, control, files)
        if control.mode != ControlMode.Solve or "main" in globals():
            control.main()
        else:
            control.ground()
            self._theory.prepare(control)
            control.solve(on_model=self._theory.on_model,
                          on_stats=self._theory.on_stats)

    def register_options(self, options: AppOptions) -> None:
        self._theory.register_options(options)

    def validate_options(self) -> None:
        self._theory.validate_options()


def run():
    lib = Library()
    enable_python(lib)
    app = ClingoDLApp(lib)
    clingo_main(lib, sys.argv[1:], app)


run()
)py";
    auto main = borrow(PyImport_AddModule("__main__"));
    if (main.get() == nullptr) {
        return nullptr;
    }
    auto scope = borrow(PyModule_GetDict(main.get()));
    if (scope.get() == nullptr) {
        return nullptr;
    }
    auto src = steal(PyBytes_FromString(source));
    if (src.get() == nullptr) {
        return nullptr;
    }
    auto compiled = steal(Py_CompileString(PyBytes_AsString(src.get()), "<string>", Py_file_input));
    if (compiled.get() == nullptr) {
        return nullptr;
    }

    auto res = steal(PyEval_EvalCode(compiled.get(), scope.get(), scope.get()));
    if (res.get() == nullptr) {
        return nullptr;
    }
    Py_RETURN_NONE;
}

// NOLINTNEXTLINE
auto methods = std::to_array<PyMethodDef>({
    {
        "create_theory",
        create_theory,
        METH_NOARGS,
        "Get the theory constructor.",
    },
    {
        "_pyclingodl",
        run_app,
        METH_NOARGS,
        "Run clingo-dl.",
    },
    {nullptr, nullptr, 0, nullptr},
});

// NOLINTNEXTLINE
PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "clingodl",
    "The clingo-dl Python module.",
    -1,
    methods.data(),
    nullptr,
    nullptr,
    nullptr,
    nullptr,
};

} // namespace

// NOLINTNEXTLINE
PyMODINIT_FUNC PyInit_clingodl() { return PyModule_Create(&module); }
