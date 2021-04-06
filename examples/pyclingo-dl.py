from clingo import Control
from clingo.ast import parse_string, ProgramBuilder
from clingodl import ClingoDLTheory

prg = '&diff { x } >= 1. &diff { y } >= 3.'

thy = ClingoDLTheory()
ctl = Control(['0'])
thy.register(ctl)
with ProgramBuilder(ctl) as bld:
    parse_string(prg, lambda ast: thy.rewrite_ast(ast, bld.add))

ctl.ground([('base', [])])
thy.prepare(ctl)
with ctl.solve(yield_=True, on_model=thy.on_model) as hnd:
    for mdl in hnd:
        print([f'{key}={val}' for key, val in thy.assignment(mdl.thread_id)])

