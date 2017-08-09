ASP modulo Difference Constraints
=================================

This project comprises a small prototype to extend clingo with difference constraints
using its C++ API and theory language.

Usage
-----

    ./clingoDL [CLINGO OPTION]... [-p] [--strict] [--rdl] [FILE]...

Option `-p` can be used to enable propagation of difference constraints,
`--strict` to turn on strict mode, and `--rdl` to use real numbers.

Example
-------

    ./clingoDL -c n=132 -p encoding-dl.lp tai4_4_1.lp

Development
-----------

### Compiling

The clingoDL system needs a C++14 conforming compiler and at least clingo version 5.3. If
both are available, it can be compiled. First run

    make FLAGS

to create the default configuration file `FLAGS`.  The default paths there
point to locations in our development environment and probably have to be
adjusted.  Note that it is not necessary to install clingo but just to compile
it. Point the include path to the folder where the `clingo.hh` file is, and the
library paths to the place where the `libclingo.so` file is.  If clingo is
installed globally the CPPFLAGS and LDFLAGS can be left empty. Then run

    make

to build clingoDL.

### Code Formatting

To format the code with clang-format simply run:

    make format

Check the Makefile to fine-tune the style.

### TODO

There are still some things that could be done to improve the implementation.

- For summing up paths, 64 bit integers could be used to avoid all possible
  overflows in practice.
- By introducing an epsilon, real valued constraints can be made to work with
  strict semantics, too.
- Some kind of preprocessing on the difference constraints might be
  interesting.
