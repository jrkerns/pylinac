"""Utility functions for pylinac."""
import decimal
import inspect
import os.path as osp
from tkinter import Tk
from collections import Iterable

import numpy


def isnumeric(object):
    """Check whether the passed object is numeric in any sense."""
    return isinstance(object, (int, float, decimal.Decimal, numpy.number))

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

def is_iterable(object):
    """Determine if an object is iterable."""
    if isinstance(object, Iterable):
        return True
    else:
        return False

def withdraw_tkinter():
    """Opens and withdraws a Tk window. Necessary so a base window doesn't open."""
    Tk().withdraw()