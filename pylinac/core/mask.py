"""Module for processing "masked" arrays, i.e. binary images."""
import numpy as np


def bounding_box(array):
    """Get the bounding box values of an ROI in a 2D array."""
    binary_arr = np.argwhere(array)
    (ymin, xmin), (ymax, xmax) = binary_arr.min(0), binary_arr.max(0) + 1
    return ymin, ymax, xmin, xmax


def filled_area_ratio(array):
    """Return the ratio of filled pixels to empty pixels in the ROI bounding box.

    For example a solid square would be 1.0, while a sold circle would be ~0.785.
    """
    ymin, ymax, xmin, xmax = bounding_box(array)
    box_area = (ymax - ymin) * (xmax - xmin)
    filled_area = np.sum(array)
    return filled_area / box_area
