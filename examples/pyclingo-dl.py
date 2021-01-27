from clingo import Control
from clingo.ast import parse_string, ProgramBuilder
from clingodl import ClingoDLTheory

prg = '&sum { x } >= 1. &sum { x } <= 3.'

thy = ClingoDLTheory()
ctl = Control(['0'])
thy.register(ctl)
with ProgramBuilder(ctl) as bld:
    parse_string(prg, lambda ast: thy.rewrite_ast(ast, bld.add))

ctl.ground([('base', [])])
thy.prepare(ctl)
with ctl.solve(yield_=True) as hnd:
    for mdl in hnd:
        print([f'{key}={val}' for key, val in thy.assignment(mdl.thread_id)])

