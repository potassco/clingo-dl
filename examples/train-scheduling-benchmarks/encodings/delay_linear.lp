#const pen_init = 180.
#const pen_steps = 5.

dpenalty(1,1).
dpenalty(pen_init,pen_init-1).
dpenalty(pen_init*X,(pen_init*X-pen_init*(X-1))) :- X=2..pen_steps.
