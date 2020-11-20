from functools import wraps
import time


def timethis(func):
    """Report execution time of function."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print(func.__name__, f'took {end-start:3.2f}s')
        return result
    return wrapper
