# Clingo-dl: A grounder and solver for solving ASP modulo Difference Constraints

![tests](https://github.com/potassco/clingoDL/workflows/tests/badge.svg)

Clingo-dl is part of the [Potassco] project for *Answer Set Programming* (ASP).
It extends ASP with constraints over difference logic and extends the ASP grounder and solver [clingo].

Please consult the following resources for further information:

  - [**Downloading source and binary releases**][download]
  - [**Installation and software requirements**](INSTALL.md)
  - [Changes between releases](CHANGES.md)
  - [Potassco clingo-dl page][home]

Clingo-dl is distributed under the [MIT License](LICENSE.md).

### Usage

    clingo-dl [OPTIONS]... [FILE]...

The system accepts all of clingo's options as well as options specific to difference constraints.
Use option `--help` to see the available options.

Furthermore, clingo-dl provides the Python module `clingodl`.

### Example

    clingo-dl -c n=132 --propagate=full examples/taskassignment/encoding-dl.lp examples/taskassignment/tai4_4_1.lp

[clingo]: https://potassco.org/clingo/
[Potassco]: https://potassco.org/
[home]: https://potassco.org/labs/clingoDL/
[download]: https://github.com/potassco/clingoDL/releases/
