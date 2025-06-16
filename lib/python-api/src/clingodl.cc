#include <clingo-dl.h>

#include <pybind11/eval.h>
#include <pybind11/pybind11.h>

namespace PyClingoDL {

namespace {

auto create_theory() -> pybind11::object {
    return pybind11::capsule{reinterpret_cast<void *>(&clingodl_create), "clingo_theory_create"};
}

auto main() {
    pybind11::exec(
        R"py(
import sys
from sys import stdout
from typing import Callable, Sequence

from clingo.app import App, AppOptions, clingo_main
from clingo.core import Library
from clingo.control import Control
from clingo.symbol import SymbolType
from clingo.theory import Theory
from clingo.solve import Model
from clingo import ast

from clingodl import create_theory


class ClingoDLApp(App):
    def __init__(self, lib: Library):
        theory = Theory(lib, create_theory())
        major, minor, revision = theory.version
        super().__init__(theory.name, f"{major}.{minor}.{revision}")
        self._lib = lib
        self._theory = theory

    def main(self, control: Control, files: Sequence[str]) -> None:
        """
        Run the main execution flow of the application.
        """
        self._theory.register(control)
        with ast.Scanner(self._lib, files) as scn:
            prg = ast.Program(self._lib)
            for stm in scn:
                self._theory.rewrite(stm, prg.add)
            control.join(prg)
        control.ground()
        with control.solve(on_model=self._on_model, on_stats=self._on_stats) as hnd:
            hnd.get()

    def register_options(self, options: AppOptions) -> None:
        """
        Register command-line options for the application.
        """
        self._theory.register_options(options)

    def validate_options(self) -> None:
        """
        Validate the options passed to the application.
        """
        self._theory.validate_options()

    def _on_model(self, model: Model):
        self._theory.on_model(model)

    def _on_stats(self, step, accu):
        self._theory.on_stats(step, accu)


def run():
    lib = Library()
    app = ClingoDLApp(lib)
    print(sys.argv[1:]);
    clingo_main(lib, sys.argv[1:], app)

run()
)py");
}

} // namespace

void register_clingodl(pybind11::module &m) {
    m.doc() = R"doc(The clingo-dl python module.)doc";
    m.def("create_theory", create_theory, R"(Get the theory constructor.)");
    m.def("main", main, R"(Run clingo-dl.)");
}

} // namespace PyClingoDL

PYBIND11_MODULE(clingodl, m) { PyClingoDL::register_clingodl(m); }
