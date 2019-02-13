// {{{ MIT License
//
// // Copyright 2018 Roland Kaminski, Philipp Wanko, Max Ostrowski
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

#include <clingo.h>

typedef struct clingodl_propagator clingodl_propagator_t;

//! creates the propagator
CLINGODL_VISIBILITY_DEFAULT bool clingodl_create_propagator(clingodl_propagator_t **prop);

//! registers the propagator with the control
CLINGODL_VISIBILITY_DEFAULT bool clingodl_register_propagator(clingodl_propagator_t *prop, clingo_control_t* ctl);

//! destroys the propagator, currently no way to unregister a propagator
CLINGODL_VISIBILITY_DEFAULT bool clingodl_destroy_propagator(clingodl_propagator_t *prop);

//! add options for your theory
CLINGODL_VISIBILITY_DEFAULT bool clingodl_add_options(clingodl_propagator_t *prop, clingo_options_t* options);

//! validate options for your theory
CLINGODL_VISIBILITY_DEFAULT bool clingodl_validate_options(clingodl_propagator_t *prop);

//! callback on every model
CLINGODL_VISIBILITY_DEFAULT bool clingodl_on_model(clingodl_propagator_t *prop, clingo_model_t* model);

//! return pointer to first element of the current (partial) assignment, 0 if not existent
CLINGODL_VISIBILITY_DEFAULT void clingodl_assignment_begin(clingodl_propagator_t *prop, uint32_t threadId, size_t *current);

//! return pointer to next element of the current (partial) assignment 0 if not existent
CLINGODL_VISIBILITY_DEFAULT bool clingodl_assignment_next(clingodl_propagator_t *prop, uint32_t threadId, size_t *current, clingo_symbol_t *name, double* value, bool *ret);

//! callback on statistic updates
/// please add a subkey with the name of your propagator
CLINGODL_VISIBILITY_DEFAULT bool clingodl_on_statistics(clingodl_propagator_t *prop, clingo_statistics_t* step, clingo_statistics_t* accu);

#ifdef __cplusplus
}
#endif

#endif
