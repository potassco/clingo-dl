# Table of Contents

- [Installation using conda](#installation-using-conda)
- [Installation using pip](#installation-using-pip)
- [Requirements](#requirements)
- [Build, Install, and Test](#build-install-and-test)
  - [Build Options](#build-options)

# Installation using conda

The latest clingo-dl release is available using the conda-forge channel:

    conda install -c conda-forge clingo-dl

Furthermore, releases and development versions can also be installed from the `potassco` and `potassco/label/dev` channels.

# Installation using pip

For the latest release use:

    pip install --upgrade clingo-dl

or the latest development version:

    pip install --upgrade --extra-index-url https://test.pypi.org/simple/ clingo-dl

# Requirements

- a C++17 conforming compiler
  - *at least* [gcc] version 7.0
  - *at least* [clang] version 4.0
  - *at least* msvc++ 14.11 ([vs][Visual Studio] 2017 15.3)
  - other compilers might work
- the [cmake] build system
  - at least version 3.16 is recommended
  - at least version 3.1 is *required*
- the [clingo] ASP solver
  - *at least* version 5.5
- optionally, the [python] programming language
  - *at least* version 3.6

# Build, Install, and Test

To build clingo-dl in its default configurations in release mode, run

    cmake -H<SOURCE_DIR> -B<BUILD_DIR> -DCMAKE_BUILD_TYPE=Release
    cmake --build <BUILD_DIR>

The resulting binaries and shared libraries will be in `<BUILD_DIR>/bin` and are ready to use.

To install all binaries and development files under cmake's install prefix (see the [build options](#build-options)), run

    cmake --build <BUILD_DIR> --target install

## Build Options

The most important options to control the build are

- Variable `CMAKE_BUILD_TYPE` should be set to `Release`. (Default: unset)
- Variable `CMAKE_INSTALL_PREFIX` controls where to install clingo-dl. (Default: `/usr/local/bin`)

Cmake's `-L` option can be used to get an overview over the variables that can be set for building clingo-dl.
To get clingo-dl specific options, run

    cmake -H<SOURCE_DIR> -B<BUILD_DIR> -DCMAKE_BUILD_TYPE=Release -LH

or, to also print important cmake specific configuration variables, run

    cmake -H<SOURCE_DIR> -B<BUILD_DIR> -DCMAKE_BUILD_TYPE=Release -LAH

Options and variables can be passed to cmake on the command line using `-D<VARIABLE>=<VALUE>`
or by editing `<BUILD_DIR>/CMakeCache.txt` after running cmake.

[gcc]: https://gcc.gnu.org/
[clang]: http://clang.llvm.org/
[msvc]: https://www.visualstudio.com/
[clingo]: https://github.com/potassco/clingo/
[cmake]: https://www.cmake.org/
[python]: https://www.python.org/
