"""Module for processing "masked" arrays, i.e. binary images."""
from typing import Tuple

import numpy as np


def bounding_box(array: np.array) -> Tuple[float, ...]:
    """Get the bounding box values of an ROI in a 2D array."""
    binary_arr = np.argwhere(array)
    (ymin, xmin), (ymax, xmax) = binary_arr.min(0), binary_arr.max(0) + 1
    return ymin, ymax, xmin, xmax


def filled_area_ratio(array: np.array) -> float:
    """Return the ratio of filled pixels to empty pixels in the ROI bounding box.

    For example a solid square would be 1.0, while a sold circle would be ~0.785.
    """
    ymin, ymax, xmin, xmax = bounding_box(array)
    box_area = (ymax - ymin) * (xmax - xmin)
    filled_area = np.sum(array)
    return float(filled_area / box_area)


def square_ratio(array: np.array) -> float:
    """Determine the width/height ratio of the ROI"""
    ymin, ymax, xmin, xmax = bounding_box(array)
    y = abs(ymax - ymin)
    x = abs(xmax - xmin)
    return y/x
