from __future__ import annotations

import numpy as np


def array_not_empty(array: np.ndarray) -> None:
    """Check an array isn't empty"""
    if not array.size:
        raise ValueError("Array must not be empty")


def single_dimension(array: np.ndarray) -> None:
    """Check an array is a single dimension"""
    if array.ndim > 1:
        raise ValueError(
            f"Array was multidimensional. Must pass 1D array; found {array.ndim}"
        )


def double_dimension(array: np.ndarray) -> None:
    """Check an array is a double dimension"""
    if array.ndim != 2:
        raise ValueError(f"Array was not 2D. Must pass 2D array; found {array.ndim}")
