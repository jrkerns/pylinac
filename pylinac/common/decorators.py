
"""Yes I know it's not Pythonic to type check, but for pylinac, I don't see many alternatives. Note that these are for Classes."""
from __future__ import unicode_literals, print_function, division, absolute_import
from builtins import zip
from builtins import str

from future import standard_library


standard_library.install_aliases()

# The following is adapted from: http://code.activestate.com/recipes/578809-decorator-to-check-method-param-types/
from functools import wraps

def type_accept(**types):
    """Function decorator enforcing input types.

    :param types: expected types (e.g. str, int, float).
    """
    def type_decor(func):  # func == function to decorate
        @wraps(func)
        def new_func(*args, **kwargs):  # the decorated function's arguments'
            # check positional arguments
            # if self is first argument, i.e. a class method, remove it. Otherwise the arg values and arg names are not the same length.
            if func.__code__.co_varnames[0] == 'self':
                arg_names = func.__code__.co_varnames[1:]
                arg_vals = args[1:]
            else:
                arg_names = func.__code__.co_varnames
                arg_vals = args
            for arg_name, arg_val in zip(func.__code__.co_varnames, arg_vals):  # [1:] syntax because of "self" parameter when used on classes.
                if arg_name in types:
                    if not isinstance(arg_val, types[arg_name]):
                        raise TypeError("Argument '{}' was not of type {}".format(arg_val, str(types[arg_name]).split("'")[1]))
            # check keyword arguments
            for arg_name, arg_val in kwargs.items():
                if arg_name in types:
                    if not isinstance(arg_val, types[arg_name]):
                        raise TypeError("Argument '{}' was not of type {}".format(arg_val, str(types[arg_name]).split("'")[1]))
            return func(*args, **kwargs)
        return new_func
    return type_decor

def value_accept(**val_accept):
    """Function decorator enforcing input values.

    :param val_accept: expected values. These can be in (lower, upper) for int or float, and in ('str1', 'str2') for string matching. E.g.
    (1, 10) means accept a value between 1 and 10. ('str1', 'str2') means the value must be one of the listed strings.
    """
    def decorator(func):  # func == function to decorate
        @wraps(func)
        def new_func(*args, **kwargs):  # the decorated function's arguments'
            # check positional arguments
            # if self is first argument, i.e. a class method, remove it. Otherwise the arg values and arg names are not the same length.
            if func.__code__.co_varnames[0] == 'self':
                arg_names = func.__code__.co_varnames[1:]
                arg_vals = args[1:]
            else:
                arg_names = func.__code__.co_varnames
                arg_vals = args
            for arg_name, arg_val in zip(arg_names, arg_vals):
                if arg_name in val_accept:
                    if type(arg_val) in (float, int):
                        if not val_accept[arg_name][0] <= arg_val <= val_accept[arg_name][1]:
                            raise ValueError("Argument '{:f}' needs to be between {:f} and {:f}".format(arg_name, val_accept[arg_name][
                                0], val_accept[arg_name][1]))
                    elif type(arg_val) in (str,):
                        if not arg_val in val_accept[arg_name]:
                            raise ValueError("Argument '{}' passed to '{}' needs to be one of these: {}".format(arg_name, func.__name__,
                                                                                                                    val_accept[arg_name]))
            # check keyword arguments
            for arg_name, arg_val in kwargs.items():
                if arg_name in val_accept:
                    if type(val_accept[arg_name][0]) in (float, int):
                        if not val_accept[arg_name][0] <= arg_val <= val_accept[arg_name][1]:
                            raise ValueError("Argument '{:f}' passed to '{}' needs to be between {:f} and {:f}".format(arg_name,
                                                                                                                      func.__name__,
                                                                                                           val_accept[arg_name][0],
                                                                                                        val_accept[
                                                                                                            arg_name][1]))
                    elif type(val_accept[arg_name][0]) in (str,):
                        if not arg_val in val_accept[arg_name]:
                            raise ValueError("Argument '{}' passed to '{}' needs to be one of these: {}".format(arg_name, func.__name__,
                                                                                                             val_accept[arg_name]))
            return func(*args, **kwargs)
        return new_func
    return decorator