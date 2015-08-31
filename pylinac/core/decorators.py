
# The following is adapted from: http://code.activestate.com/recipes/578809-decorator-to-check-method-param-types/
# Another type checking decorator: http://code.activestate.com/recipes/454322-type-checking-decorator/
from abc import ABCMeta
from functools import wraps
from inspect import signature
import time

import numpy as np


def timethis(func):
    """Report execution time of function."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print(func.__name__, 'took {:3.2f}s'.format(end-start))
        return result
    return wrapper

def type_accept(*type_args, **type_kwargs):
    """Decorator to check function/method input types. Based on Python Cookbook 3rd ed. #9.7."""
    def decorate(func):

        # Map function argument names to supplied types
        sig = signature(func)
        bound_types = sig.bind_partial(*type_args, **type_kwargs).arguments

        @wraps(func)
        def wrapper(*args, **kwargs):
            bound_values = sig.bind(*args, **kwargs)
            # Enforce type assertions across supplied arguments
            for name, value in bound_values.arguments.items():
                if name in bound_types:
                    if type(bound_types[name]) in (type, ABCMeta):  # Single-type comparisons
                        if not isinstance(value, bound_types[name]):
                            raise TypeError("Argument '{}' must be {}".format(name, bound_types[name]))
                    else:
                        if type(value) not in bound_types[name]:
                            if value not in bound_types[name]:
                                raise TypeError("Argument '{}' must be {}".format(name, bound_types[name]))
            return func(*args, **kwargs)
        return wrapper
    return decorate

def value_accept(*value_args, **value_kwargs):
    """Decorator to check function/method input types. Based on Python Cookbook 3rd ed. #9.7."""
    def decorate(func):
        sig = signature(func)
        # convert any dictionary value acceptances to tuples
        vkw = convert_dictvals2tuple(value_kwargs)
        # Map function argument names to supplied types
        bound_values = sig.bind_partial(*value_args, **vkw).arguments

        @wraps(func)
        def wrapper(*args, **kwargs):
            passed_values = sig.bind(*args, **kwargs)
            # Enforce value assertions across supplied arguments
            for name, value in passed_values.arguments.items():
                if name in bound_values:
                    if type(value) in (float, int, np.float64):
                        # value must be within a number range
                        if not bound_values[name][0] <= value <= bound_values[name][1]:
                            raise ValueError("Argument '{}' needs to be between {:f} and {:f}".format(name,
                                                                                                    bound_values[name][0],
                                                                                                    bound_values[name][1]))
                    else:
                        # value is a str and must be one of the accepted str values
                        if value not in bound_values[name]:
                            raise ValueError("Argument '{}' must be one of {}".format(name, bound_values[name]))
            return func(*args, **kwargs)
        return wrapper
    return decorate

def convert_dictvals2tuple(args):
    """Convert from dictionary to tuple of dictionary values."""
    for arg in args:
        if type(args[arg]) == dict:
            args[arg] = tuple(args[arg].values())
    return args

def unwrap_func(wrapped_func, unwraps=0):
    """Return an unwrapped function from the wrapped function.

    Functions or methods can use one or more decorators. If the decorator needs access to the base function (e.g. its arguments), it
    won't get it if more than 1 decorator is used. This function allows one to get to the base function, or a function wrapped a certain
    number of times. The function must be wrapped with the stdlib functools.wraps for proper unwrapping.

    Discussion of this issue: http://bugs.python.org/issue17482

    :param unwraps: The number of times to unwrap a function, if wrapped with multiple decorators. If 0 or negative, will return the base
        function.
    :type unwraps: int
    """
    func = wrapped_func
    # if we want the base function
    if unwraps <= 0:
        at_base = False
        while not at_base:
            try:
                func = func.__wrapped__
            except AttributeError:
                at_base = True

    # else we want the function unwrapped n times
    else:
        for unwrap in range(unwraps):
            try:
                func = func.__wrapped__
            except AttributeError:
                break

    return func