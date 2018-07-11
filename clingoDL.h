// {{{ MIT License
//
// // Copyright 2017 Roland Kaminski
//
// // Permission is hereby granted, free of charge, to any person obtaining a copy
// // of this software and associated documentation files (the "Software"), to
// // deal in the Software without restriction, including without limitation the
// // rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// // sell copies of the Software, and to permit persons to whom the Software is
// // furnished to do so, subject to the following conditions:
//
// // The above copyright notice and this permission notice shall be included in
// // all copies or substantial portions of the Software.
//
// // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// // IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// // FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// // AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// // LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// // FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// // IN THE SOFTWARE.
//
// // }}}


#ifndef CLINGODL_H
#define CLINGODL_H

//! Major version number.
#define CLINGODL_VERSION_MAJOR 1
//! Minor version number.
#define CLINGODL_VERSION_MINOR 0
//! Revision number.
#define CLINGODL_VERSION_REVISION 0
//! String representation of version.
#define CLINGODL_VERSION "1.0.0"

#ifdef __cplusplus
extern "C" {
#endif

#if defined _WIN32 || defined __CYGWIN__
#   define CLINGODL_WIN
#endif
#ifdef CLINGODL_NO_VISIBILITY
#   define CLINGODL_VISIBILITY_DEFAULT
#   define CLINGODL_VISIBILITY_PRIVATE
#else
#   ifdef CLINGODL_WIN
#       ifdef CLINGODL_BUILD_LIBRARY
#           define CLINGODL_VISIBILITY_DEFAULT __declspec (dllexport)
#       else
#           define CLINGODL_VISIBILITY_DEFAULT __declspec (dllimport)
#       endif
#       define CLINGODL_VISIBILITY_PRIVATE
#   else
#       if __GNUC__ >= 4
#           define CLINGODL_VISIBILITY_DEFAULT  __attribute__ ((visibility ("default")))
#           define CLINGODL_VISIBILITY_PRIVATE __attribute__ ((visibility ("hidden")))
#       else
#           define CLINGODL_VISIBILITY_DEFAULT
#           define CLINGODL_VISIBILITY_PRIVATE
#       endif
#   endif
#endif

#if defined __GNUC__
#define CLINGODL_DEPRECATED __attribute__((deprecated))
#elif defined _MSC_VER
#define CLINGODL_DEPRECATED __declspec(deprecated)
#else
#define CLINGODL_DEPRECATED
#endif

#include <clingo.h>

//! creates and registers the propagator
CLINGODL_VISIBILITY_DEFAULT bool theory_create_propagator(clingo_control_t* ctl);

//! destroys the propagator, currently no way to unregister a propagator
CLINGODL_VISIBILITY_DEFAULT bool theory_destroy_propagator();

//! add options for your theory
CLINGODL_VISIBILITY_DEFAULT bool theory_add_options(clingo_options_t* options);

//! validate options for your theory
CLINGODL_VISIBILITY_DEFAULT bool theory_validate_options();

//! callback on every model 
CLINGODL_VISIBILITY_DEFAULT bool theory_on_model(clingo_model_t* model);

//! callback on statistic updates
/// please add a subkey with the name of your propagator
CLINGO_VISIBILITY_DEFAULT bool theory_on_statistics(clingo_statistics_t* step, clingo_statistics_t* accu);

#ifdef __cplusplus
}
#endif

#endif
