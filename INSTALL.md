# Table of Contents

- [Requirements](#requirements)
- [Build, Install, and Test](#build-install-and-test)
  - [Build Options](#build-options)
    - [Generic Options](#generic-options)
    - [Python Support](#python-support)
    - [Lua Support](#lua-support)
  - [Compilation to JavaScript](#compilation-to-javascript)
- [Troubleshooting](#troubleshooting)
  - [Notes for Windows Users](#notes-for-windows-users)

# Requirements

- a c++14 conforming compiler
  - *at least* [gcc](https://gcc.gnu.org/) version 4.9
  - [clang](http://clang.llvm.org/) version 3.1 (using either libstdc++
    provided by gcc 4.9 or libc++)
  - *at least* msvc++ 14.0 ([Visual Studio](https://www.visualstudio.com/) 2015
    Update 3)
  - other compilers might work
- the [cmake](https://www.cmake.org/) build system
  - at least version 3.3 is recommended
  - at least version 3.1 is *required*

# Build, Install, and Test

When cloning the git repository, do not forget to update the submodules (with
source releases, you can skip this step):

    git submodule update --init --recursive

To build clingo-dl in its default configurations in release
mode, run:

    cmake -H<SOURCE_DIR> -B<BUILD_DIR> -DCMAKE_BUILD_TYPE=Release
    cmake --build <BUILD_DIR>

The resulting binaries and shared libraries will be in `<BUILD_DIR>/bin` and
are ready to use.

To install all binaries and development files under cmake's install
prefix (see the [build options](#build-options)), run:

    cmake --build <BUILD_DIR> --target install

## Build Options

Cmake's `-L` option can be used to get an overview over the variables that can
be set for building clingo-dl. To get clingo-dl specific options, run

    cmake -H<SOURCE_DIR> -B<BUILD_DIR> -DCMAKE_BUILD_TYPE=Release -LH
    
or, to also print important cmake specific configuration variables

    cmake -H<SOURCE_DIR> -B<BUILD_DIR> -DCMAKE_BUILD_TYPE=Release -LAH

Options and variables can be passed to
cmake on the command line using `-D<VARIABLE>=<VALUE>` or by editing
`<BUILD_DIR>/CMakeCache.txt` after running cmake.


Clingo-dl uses [clingo](https://github.com/potassco/clingo).
It has its own set of configuration variables:
- [building clingo](https://github.com/potassco/clingo/blob/master/INSTALL.md#build-options)

In the following, the most important options to control the build are listed.

### Generic Options

- Variable `CMAKE_BUILD_TYPE` should be set to `Release`.
- Variable `CMAKE_INSTALL_PREFIX` controls where to install clingo-dl.
- Option `CLINGODL_BUILD_WITH_SYSTEM_CLINGO` to build with the already installed clingo version. Otherwise the local source copy will be used
- Option `CLINGODL_MANAGE_RPATH` controls how to find libraries on platforms
  where this is supported, like Linux, macOS, or BSD but not Windows. This
  option should be enabled if clingo-dl is installed in a non-default location,
  like the users home directory; otherwise it has no effect.
  (Default: `ON`)

### Python and Lua Support

Python and Lua support is enabled if it is enabled in the used clingo version.
See [Python Support](https://github.com/potassco/clingo/blob/master/INSTALL.md#python-support) or [Lua Support](https://github.com/potassco/clingo/blob/master/INSTALL.md#lua-support) for clingo.
