#const pen_init = 56.
#const pen_factor = 2.
#const pen_steps = 5.

dpenalty(1,1).
dpenalty(pen_init,pen_init-1).
dpenalty(pen_init*pen_factor**X,(pen_init*pen_factor**X-pen_init*pen_factor**(X-1))) :- X=1..pen_steps-1.
