#script (python)
import clingo
import ctypes

def main(prg):
    #clingoDL = ctypes.cdll.libclingoDL # for windows
    clingoDL = ctypes.cdll.LoadLibrary("libclingoDL.so")
    clingoDL.load(ctypes.c_void_p(prg.get_c_control()))

    prg.ground([("base", [])])
    prg.solve()

#end.

{a;b}.

% both can't be true at the same time
&diff{x-y} <= 5 :- a.
&diff{y-x} <= -6 :- b.

