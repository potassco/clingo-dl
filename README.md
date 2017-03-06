Propagation Difference Constraints
==================================

This project comprises a small prototype to propagate difference constraints
using clingo's C++ API and theory language.

Usage
-----

    ./propagator [FILE]... [-- [CLINGO OPTION]...]

Example
-------

    ./propagator open_shop_dl.lp tai4_4_1.lp -- -c n=132

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
