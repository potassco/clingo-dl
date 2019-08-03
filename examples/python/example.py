#!/usr/bin/python

import sys
import clingo
import theory

class Application:
    def __init__(self, name):
        self.program_name = name
        self.version = "1.0"
        self.__theory = theory.Theory("clingodl", "clingo-dl")

    def __on_model(self, model):
        self.__theory.on_model(model)

    def register_options(self, options):
        self.__theory.register_options(options)

    def validate_options(self):
        self.__theory.validate_options()
        return True

    def __on_statistics(self, step, accu):
        self.__theory.on_statistics(step, accu)
        pass

    def main(self, prg, files):
        self.__theory.configure_propagator("propagate", "full,1")
        self.__theory.register_propagator(prg)
        if not files:
            files.append("-")
        for f in files:
            prg.load(f)

        prg.ground([("base", [])])

        # Note: this symbol is created upon propagator creation
        #       right now it is a bit tricky to get into a state where Propagator.init has been called for all propagators
        #       one possibility would be to register a propgator after all other propagators
        #       another to lazily build a lookup table when required
        adjust = self.__theory.lookup_symbol(clingo.Number(0))

        with prg.solve(on_model=self.__on_model, on_statistics=self.__on_statistics, yield_=True) as handle:
            for model in handle:
                sys.stdout.write("assignment:")
                for name, value in self.__theory.assignment(model.thread_id):
                    sys.stdout.write(" {}={}".format(name, value))
                sys.stdout.write("\n")
                if self.__theory.has_value(model.thread_id, adjust):
                    sys.stdout.write("adjustment: {}\n".format(self.__theory.get_value(model.thread_id, adjust)))

sys.exit(int(clingo.clingo_main(Application("test"), sys.argv[1:])))
