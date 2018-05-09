ASP modulo Difference Constraints
=================================

This project comprises a small prototype to extend clingo with difference constraints
using its C++ API and theory language.

Usage
-----

    clingoDL [CLINGO OPTION]... [-p] [--strict] [--rdl] [FILE]...

Option `-p` can be used to enable propagation of difference constraints,
`--strict` to turn on strict mode, and `--rdl` to use real numbers.

Example
-------

    clingoDL -c n=132 -p examples/taskassignment/encoding-dl.lp examples/taskassignment/tai4_4_1.lp

Development
-----------

### Compiling

The clingoDL system needs a C++14 conforming compiler and at least clingo version 5.3.
First run

    cmake -H. -Bbuild

to create the default configuration for building in `./build`.

    cmake --build build

to build clingoDL.

### TODO

There are still some things that could be done to improve the implementation.

- For summing up paths, 64 bit integers could be used to avoid all possible
  overflows in practice.
- By introducing an epsilon, real valued constraints can be made to work with
  strict semantics, too.
- Some kind of preprocessing on the difference constraints might be
  interesting.
- Operators and operations in difference constraints should be made more
  flexible.
