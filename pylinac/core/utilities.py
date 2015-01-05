"""Utility functions for pylinac."""
import decimal
import inspect
import os.path as osp
from tkinter import Tk
from collections import Iterable

import numpy as np

from pylinac.core.decorators import type_accept



# class Numeric(metaclass=ABCMeta):
#     """An abstract class that encompasses many numeric types.
#
#     Usage is for testing isinstance() for many numeric types: int, float, numpy datatypes, etc."""
#     pass
# Numeric.register(int)
# Numeric.register(float)
# Numeric.register(np.number)
# Numeric.register(decimal.Decimal)


def isnumeric(object):
    """Check whether the passed object is numeric in any sense."""
    return isinstance(object, (int, float, decimal.Decimal, np.number))

def is_iterable(object):
    """Determine if an object is iterable."""
    if isinstance(object, Iterable):
        return True
    else:
        return False

@type_accept(array=np.ndarray)
def array2logical(array, threshold):
    """Return a 'logical' (binary) version of the input array based on a threshold.

    :param array: numpy array to be analyzed
    :param threshold_value: int or float specifying the threshold value. If an array value is below the
        threshold value, it is converted to 0, otherwise to 1.
    """
    return np.where(array >= threshold, 1, 0)

def go_up_dirlevel(levels=0):
    """Go up directory levels from where the caller file is located.

    :param levels: Specifies how many levels to go up. 0 goes to the current directory.
    :type levels: int
    """
    calling_file = inspect.stack()[1][1]
    calling_dir = osp.dirname(calling_file)
    new_dir = calling_dir
    while levels > 0:
        old_dir = new_dir
        new_dir = osp.dirname(old_dir)
        levels -= 1
    return new_dir

def withdraw_tkinter():
    """Opens and withdraws a Tk window. Necessary so a base window doesn't open."""
    Tk().withdraw()