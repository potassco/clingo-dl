"""
This module defines a single Theory class for using a C propagator with
clingo's python library.
"""


import ctypes
import ctypes.util
import clingo

from typing import Optional, Union, Iterator, Tuple
from ctypes import c_bool, c_void_p, c_int, c_double, c_uint, c_uint64, c_size_t, c_char_p, Structure, POINTER, byref

class _c_value(ctypes.Union):
    _fields_ = [ ("integer", c_int)
               , ("double", c_double)
               , ("symbol", c_uint64)
               ]

class _c_variant(Structure):
    _fields_ = [ ("type", c_int)
               , ("value", _c_value)
               ]


class Theory:
    """
    Interface to call functions from a C-library extending clingo's C/Python
    library.

    The functions in here are designed to be used with a `clingo.Application`
    object but can also be used with a standalone `clingo.Control` object.

    Notes
    -----
    The C library must implement the following functions:

    - `bool create_propagator(propagator_t **propagator)`
    - `bool destroy_propagator(propagator_t *propagator)`
    - `bool register_propagator(propagator_t *propagator, clingo_control_t* control)`
    - `bool register_options(propagator_t *propagator, clingo_options_t* options)`
    - `bool validate_options(propagator_t *propagator)`
    - `bool on_model(propagator_t *propagator, clingo_model_t* model)`
    - `bool on_statistics(propagator_t *propagator, clingo_statistics_t* step, clingo_statistics_t* accu)`
    - `bool lookup_symbol(propagator_t *propagator, clingo_symbol_t symbol, size_t *index)`
    - `clingo_symbol_t get_symbol(propagator_t *propagator, size_t index)`
    - `void assignment_begin(propagator_t *propagator, uint32_t thread_id, size_t *index)`
    - `bool assignment_next(propagator_t *propagator, uint32_t thread_id, size_t *index)`
    - `void assignment_has_value(propagator_t *propagator, uint32_t thread_id, size_t index)`
    - `void assignment_get_value(propagator_t *propagator, uint32_t thread_id, size_t index, value_t *value)`
    - `bool clingodl_configure_propagator(clingodl_propagator_t *prop, char const *key, char const *value)`
    """

    ValueType = Union[int, float, clingo.Symbol]

    def __init__(self, prefix: str, lib: str):
        """
        Loads a given library.

        Arguments
        ---------
        prefix: str
            Prefix of functions in the library.
        lib: str
            Name of the library to load.
        """
        self.__c_propagator = None

        # load library
        self.__theory = ctypes.cdll.LoadLibrary(ctypes.util.find_library(lib))

        # bool create_propagator(propagator_t **propagator);
        self.__create_propagator = self.__fun(prefix, "create_propagator", c_bool, [POINTER(c_void_p)])

        # bool destroy_propagator(propagator_t *propagator);
        self.__destroy_propagator = self.__fun(prefix, "destroy_propagator", c_bool, [c_void_p])

        # bool register_propagator(propagator_t *propagator, clingo_control_t* control);
        self.__register_propagator = self.__fun(prefix, "register_propagator", c_bool, [c_void_p, c_void_p])

        # bool register_options(propagator_t *propagator, clingo_options_t* options);
        self.__register_options = self.__fun(prefix, "register_options", c_bool, [c_void_p, c_void_p])

        # bool validate_options(propagator_t *propagator);
        self.__validate_options = self.__fun(prefix, "validate_options", c_bool, [c_void_p])

        # bool on_model(propagator_t *propagator, clingo_model_t* model);
        self.__on_model = self.__fun(prefix, "on_model", c_bool, [c_void_p, c_void_p])

        # bool on_statistics(propagator_t *propagator, clingo_statistics_t* step, clingo_statistics_t* accu);
        self.__on_statistics = self.__fun(prefix, "on_statistics", c_bool, [c_void_p, c_void_p, c_void_p])

        # bool lookup_symbol(propagator_t *propagator, clingo_symbol_t symbol, size_t *index);
        self.__lookup_symbol = self.__fun(prefix, "lookup_symbol", c_bool, [c_void_p, c_uint64, POINTER(c_size_t)], False)

        # clingo_symbol_t get_symbol(propagator_t *propagator, size_t index);
        self.__get_symbol = self.__fun(prefix, "get_symbol", c_uint64, [c_void_p, c_size_t], False)

        # void assignment_begin(propagator_t *propagator, uint32_t thread_id, size_t *index);
        self.__assignment_begin = self.__fun(prefix, "assignment_begin", None, [c_void_p, c_uint, POINTER(c_size_t)], False)

        # bool assignment_next(propagator_t *propagator, uint32_t thread_id, size_t *index);
        self.__assignment_next = self.__fun(prefix, "assignment_next", c_bool, [c_void_p, c_uint, POINTER(c_size_t)], False)

        # void assignment_has_value(propagator_t *propagator, uint32_t thread_id, size_t index);
        self.__assignment_has_value = self.__fun(prefix, "assignment_has_value", c_bool, [c_void_p, c_uint, c_size_t], False)

        # void assignment_get_value(propagator_t *propagator, uint32_t thread_id, size_t index, value_t *value);
        self.__assignment_get_value = self.__fun(prefix, "assignment_get_value", None, [c_void_p, c_uint, c_size_t, POINTER(_c_variant)], False)

        # bool clingodl_configure_propagator(clingodl_propagator_t *prop, char const *key, char const *value);
        self.__configure_propagator = self.__fun(prefix, "configure_propagator", c_bool, [c_void_p, c_char_p, c_char_p])

        # create propagator
        c_propagator = c_void_p()
        self.__create_propagator(byref(c_propagator))
        self.__c_propagator = c_propagator

    def __del__(self):
        if self.__c_propagator is not None:
            self.__destroy_propagator(self.__c_propagator)
            self.__c_propagator = None

    def configure_propagator(self, key: str, value: str) -> None:
        """
        Allows for configuring a propagator via key/value pairs similar to
        command line options.

        This function must be called before the propagator is registered.

        Arguments
        ---------
        key: str
            The name of the option.
        value: str
            The value of the option.
        """
        self.__configure_propagator(self.__c_propagator, key.encode(), value.encode())

    def register_propagator(self, control: clingo.Control) -> None:
        """
        Register the propagator with the given control object.

        Arguments
        ---------
        control: clingo.Control
            Target to register with.
        """
        self.__register_propagator(self.__c_propagator, control._to_c)

    def register_options(self, options: clingo.ApplicationOptions) -> None:
        """
        Register the propagator's options with the given application options
        object.

        Arguments
        ---------
        options: clingo.ApplicationOptions
            Target to register with.
        """

        self.__register_options(self.__c_propagator, options._to_c)

    def validate_options(self) -> None:
        """
        Validate the options of the propagator.
        """
        self.__validate_options(self.__c_propagator)

    def on_model(self, model: clingo.Model) -> None:
        """
        Inform the propagator that a model has been found.

        Arguments
        ---------
        model: clingo.Model
            The current model.
        """
        self.__on_model(self.__c_propagator, model._to_c)

    def on_statistics(self, step: clingo.StatisticsMap, accu: clingo.StatisticsMap) -> None:
        """
        Add the propagator's statistics to the given maps.

        Arguments
        ---------
        step: clingo.StatisticsMap
            Map for per step statistics.
        accu: clingo.StatisticsMap
            Map for accumulated statistics.
        """
        self.__on_statistics(self.__c_propagator, step._to_c, accu._to_c)

    def lookup_symbol(self, symbol: clingo.Symbol) -> Optional[int]:
        """
        Get the integer index of a symbol assigned by the propagator when a
        model is found.

        Using indices allows for efficent retreival of values.

        Arguments
        ---------
        symbol: clingo.Symbol
            The symbol to look up.

        Returns
        -------
        Optional[int]
            The index of the value if found.
        """
        c_index = c_size_t()
        if self.__lookup_symbol(self.__c_propagator, symbol._to_c, byref(c_index)):
            return c_index.value
        else:
            return None

    def get_symbol(self, index: int) -> clingo.Symbol:
        """
        Get the symbol associated with an index.

        The index must be valid.

        Arguments
        ---------
        index: int
            Index to retreive.

        Returns
        -------
        clingo.Symbol
            The associated symbol.
        """
        return clingo._Symbol(self.__get_symbol(self.__c_propagator, index))

    def has_value(self, thread_id: int, index: int) -> bool:
        """
        Check if the given symbol index has a value in the current model.

        Arguments
        ---------
        thread_id: int
            The index of the solving thread that found the model.
        index: int
            Index to retreive.

        Returns
        -------
        bool
            Whether the given index has a value.
        """
        return self.__assignment_has_value(self.__c_propagator, thread_id, index)

    def get_value(self, thread_id: int, index: int) -> ValueType:
        """
        Get the value of the symbol index in the current model.

        Arguments
        ---------
        thread_id: int
            The index of the solving thread that found the model.
        index: int
            Index to retreive.

        Returns
        -------
        ValueType
            The value of the index in form of an int, float, or clingo.Symbol.
        """
        c_value = _c_variant()
        self.__assignment_get_value(self.__c_propagator, thread_id, index, byref(c_value))
        if c_value.type == 0:
            return c_value.value.integer
        elif c_value.type == 1:
            return c_value.value.double
        elif c_value.type == 2:
            return clingo._Symbol(c_value.value.symbol)
        else:
            return None

    def assignment(self, thread_id: int) -> Iterator[Tuple[clingo.Symbol,ValueType]]:
        """
        Get all values symbol/value pairs assigned by the propagator in the
        current model.

        Arguments
        ---------
        thread_id: int
            The index of the solving thread that found the model.

        Returns
        -------
        Iterator[Tuple[clingo.Symbol,ValueType]]
            An iterator over symbol/value pairs.
        """
        c_index = c_size_t()
        self.__assignment_begin(self.__c_propagator, thread_id, byref(c_index))
        while self.__assignment_next(self.__c_propagator, thread_id, byref(c_index)):
            yield (self.get_symbol(c_index), self.get_value(thread_id, c_index))

    def __fun(self, prefix, name, res, args, error=True):
        ret = self.__theory["{}_{}".format(prefix, name)]
        ret.restype = res
        ret.argtypes = args
        ret.errcheck = self.__handle_error if error else self.__skip_error
        return ret

    def __skip_error(self, ret, func, arguments):
        return ret

    def __handle_error(self, success, func, arguments):
        if not success:
            msg = clingo._error_message()
            code = clingo._error_code()
            if msg is None:
                msg = "no message"
            if code in (1, 2, 4):
                raise RuntimeError(msg)
            if code == 3:
                raise MemoryError(msg)
            raise RuntimeError("unknow error")

