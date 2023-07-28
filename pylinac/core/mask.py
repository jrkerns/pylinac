"""Module for processing "masked" arrays, i.e. binary images."""
from __future__ import annotations

import numpy as np


def bounding_box(array: np.ndarray) -> (float, ...):
    # TODO: replace with regionprops
    """Get the bounding box values of an ROI in a 2D array."""
    binary_arr = np.argwhere(array)
    (ymin, xmin), (ymax, xmax) = binary_arr.min(0), binary_arr.max(0) + 1
    return ymin, ymax, xmin, xmax
