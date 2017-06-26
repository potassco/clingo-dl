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

Currently, the prototype only compiles out of the box in the haiti environment.
Just run:

    make

### Code Formatting

To format the code with clang-format simply run:

    make format

Check the Makefile to fine-tune the style.
