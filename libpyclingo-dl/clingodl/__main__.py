#!/usr/bin/python

import sys
import clingo
from clingo import ast
from . import ClingoDLTheory

class Application(clingo.Application):
    def __init__(self, name):
        self.program_name = name
        self.version = "1.2.0"
        self.__theory = ClingoDLTheory()

    def register_options(self, options):
        self.__theory.register_options(options)

    def validate_options(self):
        self.__theory.validate_options()
        return True

    def print_model(self, model, default_printer):
        # print model
        symbols = model.symbols(shown=True)
        sys.stdout.write(" ".join(str(symbol) for symbol in sorted(symbols) if not self.__hidden(symbol)))
        sys.stdout.write('\n')

        # print assignment
        sys.stdout.write('Assignment:\n')
        symbols = model.symbols(theory=True)
        assignment = []
        cost = None
        for symbol in sorted(symbols):
            if symbol.match("dl", 2):
                assignment.append("{}={}".format(*symbol.arguments))
        sys.stdout.write(" ".join(assignment))
        sys.stdout.write('\n')

        sys.stdout.flush()

    def main(self, ctl, files):
        self.__theory.register(ctl)

        with ast.ProgramBuilder(ctl) as bld:
            ast.parse_files(files, lambda stm: self.__theory.rewrite_ast(stm, bld.add))

        ctl.ground([("base", [])])
        self.__theory.prepare(ctl)

        ctl.solve(on_model=self.__on_model, on_statistics=self.__on_statistics)

    def __on_model(self, model):
        self.__theory.on_model(model)

    def __on_statistics(self, step, accu):
        self.__theory.on_statistics(step, accu)

    def __hidden(self, symbol):
        return symbol.type == clingo.SymbolType.Function and symbol.name.startswith("__")

if __name__ == "__main__":
    sys.exit(int(clingo.clingo_main(Application("clingo-dl"), sys.argv[1:])))
