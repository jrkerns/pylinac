from __future__ import annotations

import copy
import math
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
from skimage import measure
from skimage.measure._regionprops import RegionProperties

from pylinac.core.array_utils import bit_invert, convert_to_dtype
from pylinac.core.geometry import Point

if TYPE_CHECKING:
    from pylinac.core.image import BaseImage


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
    bb_size = kwargs["bb_size"]
    larger_bb_area = np.pi * ((bb_size + 2) / 2) ** 2
    smaller_bb_area = max(
        (np.pi * ((bb_size - 2) / 2) ** 2, 2)
    )  # set a min of 2 to avoid a lower bound of 0 when radius=2. This is much more likely to find noise in a block.
    return smaller_bb_area < bb_area < larger_bb_area


def is_round(region: RegionProperties, *args, **kwargs) -> bool:
    """Decide if the ROI is circular in nature by testing the filled area vs bounding box. Used to find the BB."""
    expected_fill_ratio = np.pi / 4  # area of a circle inside a square
    actual_fill_ratio = region.filled_area / region.bbox_area
    return expected_fill_ratio * 1.2 > actual_fill_ratio > expected_fill_ratio * 0.8


def is_square(region: RegionProperties, *args, **kwargs) -> bool:
    """Decide if the ROI is square in nature by testing the filled area vs bounding box. Used to find the BB."""
    actual_fill_ratio = region.filled_area / region.bbox_area
    return actual_fill_ratio > 0.8


def is_right_size_square(region: RegionProperties, *args, **kwargs) -> bool:
    """Decide if the ROI is square in nature by testing the filled area vs bounding box. Used to find the BB."""
    field_area = region.area_filled / (kwargs["dpmm"] ** 2)
    rad_size = max((kwargs["rad_size"], 5))
    larger_bb_area = (rad_size + 5) ** 2
    smaller_bb_area = (rad_size - 5) ** 2
    return smaller_bb_area < field_area < larger_bb_area


class MetricBase(ABC):
    """Base class for any 2D metric. This class is abstract and should not be instantiated.

    The subclass should implement the ``calculate`` method and the ``name`` attribute.

    As a best practice, the ``image_compatibility`` attribute should be set to a list of image classes that the metric
    is compatible with. Image types that are not in the list will raise an error. This allows
    compatibility to be explicit. However, by default this is None and no compatibility checking is done.
    """

    unit: str = ""
    image: BaseImage
    image_compatibility: list[BaseImage] | None = None

    def inject_image(self, image: BaseImage):
        """Inject the image into the metric."""
        if self.image_compatibility is not None and not isinstance(
            image, self.image_compatibility
        ):
            raise TypeError(f"Image must be one of {self.image_compatibility}")
        self.image = image

    def context_calculate(self) -> Any:
        """Calculate the metric, passing in an image copy so that
        modifications to the image don't affect the original.

        This is also **kinda** memory efficient since the original
        image is a reference. The copy here will get destroyed
        after the call returns vs keeping a copy around.

        So at any given time, only 2x the memory is required instead of
        Nx. This is important when computing multiple metrics.
        """
        image_copy = copy.deepcopy(self.image)
        return self.calculate(image_copy)

    @abstractmethod
    def calculate(self, image: BaseImage) -> Any:
        """Calculate the metric. Can return anything"""
        pass

    def plot(self, axis: plt.Axes):
        """Plot the metric"""
        pass

    def additional_plots(self) -> list[plt.figure]:
        """Plot additional information on a separate figure as needed.

        This should NOT show the figure. The figure will be shown
        via the ``metric_plots`` method. Calling show here would
        block other metrics from plotting their own separate metrics.
        """
        pass

    @property
    def name(self) -> str:
        """The name of the metric"""
        raise NotImplementedError("The metric must have a name")


class DiskRegion(MetricBase):
    """A metric to find a disk/BB in an image where the BB is near an expected position and size.
    This will calculate the scikit-image regionprops of the BB."""

    x_offset: float
    y_offset: float

    def __init__(
        self,
        expected_position_mm: Point | tuple[float, float],
        search_window_mm: tuple[float, float],
        radius_mm: float,
        radius_tolerance_mm: float = 0.25,
        detection_conditions: list[callable] = (is_round, is_right_size_bb),
        name="Disk Region",
    ):
        self.expected_position = Point(expected_position_mm)
        self.radius = radius_mm
        self.tolerance = radius_tolerance_mm
        self.search_window = search_window_mm
        self.detection_conditions = detection_conditions
        self.name = name

    def calculate(self, image: BaseImage) -> RegionProperties:
        arr = image.array
        arr_uint16 = convert_to_dtype(arr, np.uint16)
        arr_inverted = bit_invert(arr_uint16)
        # sample the image in the search window; need to convert to mm
        left = math.floor(self.expected_position.x - self.search_window[0] / 2)
        right = math.ceil(self.expected_position.x + self.search_window[0] / 2)
        top = math.floor(self.expected_position.y - self.search_window[1] / 2)
        bottom = math.ceil(self.expected_position.y + self.search_window[1] / 2)
        sample = arr_inverted[top:bottom, left:right]
        # search for the BB by iteratively lowering the low-pass threshold value until the BB is found.
        found = False
        threshold_percentile = 5
        while not found:
            try:
                binary_array = sample > np.percentile(sample, threshold_percentile)
                labeled_arr = measure.label(binary_array)
                regions = measure.regionprops(labeled_arr)
                conditions_met = [
                    all(
                        condition(
                            region,
                            dpmm=image.dpmm,
                            bb_size=self.radius,
                            shape=binary_array.shape,
                        )
                        for condition in self.detection_conditions
                    )
                    for region in regions
                ]
                if not any(conditions_met):
                    raise ValueError
                else:
                    region_idx = [
                        idx for idx, value in enumerate(conditions_met) if value
                    ][0]
                    found = True
            except (IndexError, ValueError):
                threshold_percentile += 2
                if threshold_percentile >= 100:
                    raise ValueError("Couldn't find a disk")

        # determine the center of mass of the BB
        # we invert so BB intensity increases w/ attenuation
        inv_img = bit_invert(sample)
        bb_rprops = measure.regionprops(labeled_arr, intensity_image=inv_img)[
            region_idx
        ]
        self.x_offset = left
        self.y_offset = top
        return bb_rprops


class DiskLocator(DiskRegion):
    point: Point

    def calculate(self, image) -> Point:
        """Get the weighted centroid of the region prop of the BB."""
        region = super().calculate(image)
        self.point = Point(
            region.weighted_centroid[1] + self.x_offset,
            region.weighted_centroid[0] + self.y_offset,
        )
        return self.point

    def plot(self, axis: plt.Axes):
        """Plot the BB center"""
        axis.plot(self.point.x, self.point.y, "ro")


class GlobalFieldLocator(MetricBase):
    def __init__(
        self,
        low_threshold_percentile: float = 5,
        high_threshold_percentile: float = 99.9,
        name: str = "Field Finder",
    ):
        self.low_threshold_percentile = low_threshold_percentile
        self.high_threshold_percentile = high_threshold_percentile
        self.name = name

    def calculate(self, image: BaseImage) -> Point:
        min, max = np.percentile(
            image.array, [self.low_threshold_percentile, self.high_threshold_percentile]
        )
        threshold_img = image.as_binary((max - min) / 2 + min)
        filled_img = ndimage.binary_fill_holes(threshold_img)
        coords = ndimage.center_of_mass(filled_img)
        return Point(x=coords[-1], y=coords[0])


class GlobalSizedFieldRegion(MetricBase):
    def __init__(
        self,
        low_threshold_percentile: float = 5,
        high_threshold_percentile: float = 99.9,
        field_size: float = 10,
        name: str = "Field Finder",
    ):
        self.low_threshold_percentile = low_threshold_percentile
        self.high_threshold_percentile = high_threshold_percentile
        self.name = name

    def calculate(self, image: BaseImage) -> RegionProperties:
        min, max = np.percentile(
            image.array, [self.low_threshold_percentile, self.high_threshold_percentile]
        )
        threshold_img = image.as_binary((max - min) / 2 + min)
        filled_img = ndimage.binary_fill_holes(threshold_img)
        return RegionProperties(filled_img)
