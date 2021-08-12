# Clingo-dl: A grounder and solver for solving ASP modulo Difference Constraints

![tests](https://github.com/potassco/clingoDL/workflows/tests/badge.svg)

Clingo-dl is part of the [Potassco](https://potassco.org) project for *Answer Set
Programming* (ASP). 
It extends ASP with constraints over difference logic
and extends the ASP grounder and solver [clingo](https://potassco.org/clingo/).


Please consult the following resources for further information:

  - [**Downloading source and binary releases**](https://github.com/potassco/clingoDL/releases)
  - [**Installation and software requirements**](INSTALL.md)
  - [Changes between releases](CHANGES.md)
  - [Potassco clingo-dl page](https://potassco.org/labs/clingodl/)

Clingo-dl is distributed under the [MIT License](LICENSE.md).

### Usage

    clingo-dl [CLINGO OPTION]... [--propagate MODE] [--strict] [--rdl] [FILE]...

Option `--propagate MODE` can be used to enable propagation of difference
constraints, `--strict` to turn on strict mode, and `--rdl` to use real
numbers.

### Example

    clingo-dl -c n=132 --propagate=full examples/taskassignment/encoding-dl.lp examples/taskassignment/tai4_4_1.lp

For an example of how to use the clingoDL library, see [here](https://github.com/potassco/clingo-dl/blob/master/examples/pyclingo-dl.py)
