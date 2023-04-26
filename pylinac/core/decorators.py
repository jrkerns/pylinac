import functools
import inspect
import weakref
from typing import Iterable


def lru_cache(*lru_args, **lru_kwargs):
    """Method-safe LRU cache; https://stackoverflow.com/a/33672499"""

    def decorator(func):
        @functools.wraps(func)
        def wrapped_func(self, *args, **kwargs):
            # We're storing the wrapped method inside the instance. If we had
            # a strong reference to self the instance would never die.
            self_weak = weakref.ref(self)

            @functools.wraps(func)
            @functools.lru_cache(*lru_args, **lru_kwargs)
            def cached_method(*args, **kwargs):
                return func(self_weak(), *args, **kwargs)

            setattr(self, func.__name__, cached_method)
            return cached_method(*args, **kwargs)

        return wrapped_func

    return decorator


def validate(**validate_kwargs):
    """Validate arguments to a function with validator-like functions.

    def is_float(value):
        if not isinstance(value, float):
            raise ValueError

    @validate(a=is_float)
    def add(a, b):
        return a + b

    # this will fail
    add(3, 4.4)

    # this is fine
    add(1.1, 2)  # b is not checked
    """

    def decorator(func):
        sig = inspect.signature(func)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            passed_values = sig.bind(*args, **kwargs).arguments
            for arg, value in passed_values.items():
                if arg in validate_kwargs.keys():
                    if isinstance(validate_kwargs[arg], Iterable):
                        for validator in validate_kwargs[arg]:
                            validator(value)
                    else:
                        validate_kwargs[arg](value)
            res = func(*args, **kwargs)
            return res

        return wrapper

    return decorator


def between(**validate_kwargs):
    """Validate arguments to a function with validator-like functions.

    def is_float(value):
        if not isinstance(value, float):
            raise ValueError

    @validate(a=is_float)
    def add(a, b):
        return a + b

    # this will fail
    add(3, 4.4)

    # this is fine
    add(1.1, 2)  # b is not checked
    """

    def decorator(func):
        sig = inspect.signature(func)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            passed_values = sig.bind(*args, **kwargs).arguments
            for arg, value in passed_values.items():
                if arg in validate_kwargs.keys():
                    bounds = validate_kwargs[arg]
                    if value > bounds[1] or value < bounds[0]:
                        raise ValueError(
                            f"The value {value} must be between {bounds[0]} and {bounds[1]}"
                        )
            res = func(*args, **kwargs)
            return res

        return wrapper

    return decorator


def bounds(**bounds_kwargs):
    message = bounds_kwargs.pop(
        "message", "Argument '{}' passed to '{}' was not in the range: ({}, {})"
    )
    error_type = bounds_kwargs.pop("error_type", ValueError)

    for key, values in bounds_kwargs.items():
        if len(values) != 2:
            raise ValueError(
                "There must only be two items per argument to define a bound range"
            )
        for value in values:
            if not isinstance(value, (int, float)):
                raise TypeError(
                    f"The range values passed ({values}) are not ints or floats"
                )

    def decorator(func):

        sig = inspect.signature(func)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            start_stop = sig.bind_partial(**bounds_kwargs).arguments
            passed_values = sig.bind(*args, **kwargs).arguments
            for name, value in passed_values.items():
                if name in start_stop:
                    start = start_stop[name][0]
                    stop = start_stop[name][1]
                    if start > value or value > stop:
                        raise error_type(
                            message.format(value, func.__name__, start, stop)
                        )
            res = func(*args, **kwargs)
            return res

        return wrapper

    return decorator
