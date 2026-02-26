import gc

from clingo.theory import Theory
from clingo.core import Library
from clingo.control import Control
from clingo import ast

from clingodl import create_theory


def test_solve():
    """
    Test theory functionality.
    """
    lib = Library()
    thy = Theory(lib, create_theory())

    ctl = Control(lib)
    thy.register(ctl)
    thy.rewrite_string(lib, ctl, "a. b. c. &diff{x - y} <= -1.")
    ctl.ground()
    thy.prepare(ctl)

    models = []

    def on_model(model):
        nonlocal models
        thy.on_model(model)
        ass = thy.assignment(model.thread_id)
        models.append(
            (
                [str(sym) for sym in sorted(model.symbols(shown=True))],
                [(str(sym), val) for sym, val in sorted(ass)],
            )
        )

    ctl.solve(on_model=on_model)

    assert models == [(["__dl(x,0)", "__dl(y,1)", "a", "b", "c"], [("x", 0), ("y", 1)])]

    gc.collect()
