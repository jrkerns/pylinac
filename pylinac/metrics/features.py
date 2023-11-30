from __future__ import annotations

import numpy as np
from skimage.measure._regionprops import RegionProperties


def is_symmetric(region: RegionProperties, *args, **kwargs) -> bool:
    """Whether the binary object's dimensions are symmetric, i.e. a perfect circle. Used to find the BB."""
    ymin, xmin, ymax, xmax = region.bbox
    y = abs(ymax - ymin)
    x = abs(xmax - xmin)
    if x > max(y * 1.05, y + 3) or x < min(y * 0.95, y - 3):
        return False
    return True


def is_near_center(region: RegionProperties, *args, **kwargs) -> bool:
    """Whether the bb is <2cm from the center of the field"""
    dpmm = kwargs["dpmm"]
    shape = kwargs["shape"]
    extent_limit_mm = 20
    bottom, left, top, right = region.bbox
    bb_center_x = left + (right - left) / 2
    bb_center_y = bottom + (top - bottom) / 2
    x_lo_limit = shape[1] / 2 - dpmm * extent_limit_mm
    x_hi_limit = shape[1] / 2 + dpmm * extent_limit_mm
    is_bb_x_centered = x_lo_limit < bb_center_x < x_hi_limit
    y_lo_limit = shape[0] / 2 - dpmm * extent_limit_mm
    y_hi_limit = shape[0] / 2 + dpmm * extent_limit_mm
    is_bb_y_centered = y_lo_limit < bb_center_y < y_hi_limit
    return is_bb_x_centered and is_bb_y_centered


def is_right_size_bb(region: RegionProperties, *args, **kwargs) -> bool:
    """Decide whether the ROI is roughly the size of a BB; not noise and not an artifact. Used to find the BB."""
    bb_area = region.area_filled / (kwargs["dpmm"] ** 2)
    bb_size = kwargs["bb_size"]  # radius in mm
    tolerance = kwargs["tolerance"]  # radius tolerance in mm
    # A = pi * r^2
    larger_bb_area = np.pi * (bb_size + tolerance) ** 2
    smaller_bb_area = max(
        (np.pi * (bb_size - tolerance) ** 2, 2)
    )  # set a min of 2 to avoid a lower bound of 0 when radius<=2. Having a small lower bound is much more likely to find noise in a block.
    # this is actually really important. A lower bound of 1 will catch SIGNIFICANT noise and produce erroneous results.
    return smaller_bb_area < bb_area < larger_bb_area


def is_solid(region: RegionProperties, *args, **kwargs) -> bool:
    """Whether the ROI is spiculated. We want nice, round ROIs,
    and this will drop such ROIs. Generally, these spiculations are noise or a BB rod.
    """
    return region.solidity > 0.9


def is_round(region: RegionProperties, *args, **kwargs) -> bool:
    """Decide if the ROI is circular in nature by testing the filled area vs bounding box. Used to find the BB."""
    expected_fill_ratio = np.pi / 4  # area of a circle inside a square
    actual_fill_ratio = region.filled_area / region.bbox_area
    return expected_fill_ratio * 1.2 > actual_fill_ratio > expected_fill_ratio * 0.8


def is_right_circumference(region: RegionProperties, *args, **kwargs) -> bool:
    """Test the regionprop's perimeter attr to see if it matches
    that of an equivalent circle"""
    upper_circumference = 2 * np.pi * (kwargs["bb_size"] + kwargs["tolerance"])
    lower_circumference = 2 * np.pi * (kwargs["bb_size"] - kwargs["tolerance"])
    actual_perimeter = region.perimeter / kwargs["dpmm"]
    return upper_circumference > actual_perimeter > lower_circumference


def is_right_square_perimeter(region: RegionProperties, *args, **kwargs) -> bool:
    """Test the regionprop's perimeter attr to see if it matches
    that of an equivalent square. In reality, edges aren't perfectly straight, so
    the real perimeter is always going to be higher than the theoretical perimeter.
    We thus add a larger tolerance (20%) to the upper perimeter"""
    actual_perimeter = region.perimeter / kwargs["dpmm"]
    upper_perimeter = 1.20 * 2 * (
        kwargs["field_width_mm"] + kwargs["field_tolerance_mm"]
    ) + 2 * (kwargs["field_height_mm"] + kwargs["field_tolerance_mm"])
    lower_perimeter = 2 * (
        kwargs["field_width_mm"] - kwargs["field_tolerance_mm"]
    ) + 2 * (kwargs["field_height_mm"] - kwargs["field_tolerance_mm"])
    return upper_perimeter > actual_perimeter > lower_perimeter


def is_square(region: RegionProperties, *args, **kwargs) -> bool:
    """Decide if the ROI is square in nature by testing the filled area vs bounding box. Used to find the BB."""
    actual_fill_ratio = region.filled_area / region.bbox_area
    return actual_fill_ratio > 0.8


def is_right_area_square(region: RegionProperties, *args, **kwargs) -> bool:
    """Decide if the ROI is square in nature by testing the filled area vs bounding box. Used to find the BB."""
    field_area = region.area_filled / (kwargs["dpmm"] ** 2)
    low_bound_expected_area = (
        kwargs["field_width_mm"] - kwargs["field_tolerance_mm"]
    ) * (kwargs["field_height_mm"] - kwargs["field_tolerance_mm"])
    high_bound_expected_area = (
        kwargs["field_width_mm"] + kwargs["field_tolerance_mm"]
    ) * (kwargs["field_height_mm"] + kwargs["field_tolerance_mm"])
    return low_bound_expected_area < field_area < high_bound_expected_area
