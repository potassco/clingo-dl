Propagation Difference Constraints
==================================

This project comprises a small prototype to propagate difference constraints
using clingo's C++ API and theory language.

Usage
-----

    ./propagator [-p] [FILE]... [-- [CLINGO OPTION]...]

Option `-p` can be used to enable propagation of difference constraints.

Example
-------

    ./propagator -p encoding-dl.lp tai4_4_1.lp -- -c n=132

Development
-----------

### Compiling

The propagator needs a C++14 conforming compiler and the clingo solver. If both
are available, the propagator can be compiled. First run

    make FLAGS

to the default configuration file `FLAGS`.  The default paths there point to
locations in our development environment and probably have to be adjusted.
Note that it is not necessary to install clingo but just to compile it. Point
the include path to the folder where the `clingo.hh` file is, and the library
paths to the place where the `libclingo.so` file is.  If clingo is installed
globally the CPPFLAGS and LDFLAGS can be left empty. Then run

    make

to build the propagator.

### Code Formatting

To format the code with clang-format simply run:

    make format

Check the Makefile to fine-tune the style.

### TODO

There are still some things that could be done to improve the implementation.

- For summing up paths, 64 bit integers could be used to avoid all possible
  overflows in practice.
- The shortest path search could still be sped up by using a state-of-the-art
  algorthim.
- There might be further opportunities to terminate the shortest path search
  early, e.g., if no more candidate nodes are reachable. Maybe it is possible
  to implement this efficiently.
- By introducing an epsilon, real valued constraints can be made to work with
  strict semantics, too.