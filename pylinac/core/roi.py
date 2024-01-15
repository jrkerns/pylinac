from __future__ import annotations

from functools import cached_property

import matplotlib.pyplot as plt
import numpy as np
from skimage import draw
from skimage.measure._regionprops import _RegionProperties

from .contrast import Contrast, contrast, michelson, ratio, rms, visibility, weber
from .decorators import lru_cache
from .geometry import Circle, Point, Rectangle
from .image import ArrayImage


def bbox_center(region: _RegionProperties) -> Point:
    """Return the center of the bounding box of an scikit-image region.

    Parameters
    ----------
    region
        A scikit-image region as calculated by skimage.measure.regionprops().

    Returns
    -------
    point : :class:`~pylinac.core.geometry.Point`
    """
    bbox = region.bbox
    y = abs(bbox[0] - bbox[2]) / 2 + min(bbox[0], bbox[2])
    x = abs(bbox[1] - bbox[3]) / 2 + min(bbox[1], bbox[3])
    return Point(x, y)


class DiskROI(Circle):
    """An class representing a disk-shaped Region of Interest."""

    def __init__(
        self,
        array: np.ndarray,
        angle: float,
        roi_radius: float,
        dist_from_center: float,
        phantom_center: tuple | Point,
    ):
        """
        Parameters
        ----------
        array : ndarray
            The 2D array representing the image the disk is on.
        angle : int, float
            The angle of the ROI in degrees from the phantom center.
        roi_radius : int, float
            The radius of the ROI from the center of the phantom.
        dist_from_center : int, float
            The distance of the ROI from the phantom center.
        phantom_center : tuple
            The location of the phantom center.
        """
        center = self._get_shifted_center(angle, dist_from_center, phantom_center)
        super().__init__(center_point=center, radius=roi_radius)
        self._array = array

    @staticmethod
    def _get_shifted_center(
        angle: float,
        dist_from_center: float,
        phantom_center: Point,
    ) -> Point:
        """The center of the ROI; corrects for phantom dislocation and roll."""
        y_shift = np.sin(np.deg2rad(angle)) * dist_from_center
        x_shift = np.cos(np.deg2rad(angle)) * dist_from_center
        return Point(phantom_center.x + x_shift, phantom_center.y + y_shift)

    @cached_property
    def pixel_values(self) -> np.ndarray:
        return self.circle_mask()

    @cached_property
    def pixel_value(self) -> float:
        """The median pixel value of the ROI."""
        masked_img = self.circle_mask()
        return float(np.median(masked_img))

    @cached_property
    def std(self) -> float:
        """The standard deviation of the pixel values."""
        masked_img = self.circle_mask()
        return float(np.std(masked_img))

    @lru_cache()
    def circle_mask(self) -> np.ndarray:
        """Return a mask of the image, only showing the circular ROI."""
        rr, cc = draw.disk(center=(self.center.y, self.center.x), radius=self.radius)
        return self._array[rr, cc]

    def plot2axes(
        self,
        axes: plt.Axes | None = None,
        edgecolor: str = "black",
        fill: bool = False,
        text: str = "",
        fontsize: str = "medium",
        **kwargs,
    ) -> None:
        """Plot the Circle on the axes.

        Parameters
        ----------
        axes : matplotlib.axes.Axes
            An MPL axes to plot to.
        edgecolor : str
            The color of the circle.
        fill : bool
            Whether to fill the circle with color or leave hollow.
        text: str
            If provided, plots the given text at the center. Useful for differentiating ROIs on a plotted image.
        fontsize: str
            The size of the text, if provided. See https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
            for options.
        """
        if axes is None:
            fig, axes = plt.subplots()
            axes.imshow(self._array)
        super().plot2axes(
            axes, edgecolor=edgecolor, text=str(text), fontsize=fontsize, **kwargs
        )

    def as_dict(self) -> dict:
        """Convert to dict. Useful for dataclasses/Result"""
        data = super().as_dict()
        data.update({"median": self.pixel_value, "std": self.std})
        return data


class LowContrastDiskROI(DiskROI):
    """A class for analyzing the low-contrast disks."""

    contrast_threshold: float | None
    cnr_threshold: float | None
    contrast_reference: float | None

    def __init__(
        self,
        array: np.ndarray | ArrayImage,
        angle: float,
        roi_radius: float,
        dist_from_center: float,
        phantom_center: tuple | Point,
        contrast_threshold: float | None = None,
        contrast_reference: float | None = None,
        cnr_threshold: float | None = None,
        contrast_method: str = Contrast.MICHELSON,
        visibility_threshold: float | None = 0.1,
    ):
        """
        Parameters
        ----------
        contrast_threshold : float, int
            The threshold for considering a bubble to be "seen".
        """
        super().__init__(array, angle, roi_radius, dist_from_center, phantom_center)
        self.contrast_threshold = contrast_threshold
        self.cnr_threshold = cnr_threshold
        self.contrast_reference = contrast_reference
        self.contrast_method = contrast_method
        self.visibility_threshold = visibility_threshold

    @property
    def _contrast_array(self) -> np.ndarray:
        return np.array((self.pixel_value, self.contrast_reference))

    @property
    def signal_to_noise(self) -> float:
        """The signal-to-noise ratio. Cast to numpy first to use numpy overflow handling."""
        return float(np.array(self.pixel_value) / self.std)

    @property
    def contrast_to_noise(self) -> float:
        """The contrast to noise ratio of the ROI. Cast to numpy first to use numpy overflow handling."""
        return float(np.array(self.contrast) / self.std)

    @property
    def michelson(self) -> float:
        """The Michelson contrast"""
        return michelson(self._contrast_array)

    @property
    def weber(self) -> float:
        """The Weber contrast"""
        return weber(feature=self.pixel_value, background=self.contrast_reference)

    @property
    def rms(self) -> float:
        """The root-mean-square contrast"""
        return rms(self._contrast_array)

    @property
    def ratio(self) -> float:
        """The ratio contrast"""
        return ratio(self._contrast_array)

    @property
    def contrast(self) -> float:
        """The contrast of the bubble. Uses the contrast method passed in the constructor. See https://en.wikipedia.org/wiki/Contrast_(vision)."""
        return contrast(self._contrast_array, self.contrast_method)

    @property
    def cnr_constant(self) -> float:
        """The contrast-to-noise value times the bubble diameter."""
        DeprecationWarning(
            "The 'cnr_constant' property will be deprecated in a future release. Use .visibility instead."
        )
        return self.contrast_to_noise * self.diameter

    @property
    def visibility(self) -> float:
        """The visual perception of CNR. Uses the model from A Rose: https://www.osapublishing.org/josa/abstract.cfm?uri=josa-38-2-196.
        See also here: https://howradiologyworks.com/x-ray-cnr/.
        Finally, a review paper here: http://xrm.phys.northwestern.edu/research/pdf_papers/1999/burgess_josaa_1999.pdf
        Importantly, the Rose model is not applicable for high-contrast use cases."""
        return visibility(
            array=self._contrast_array,
            radius=self.radius,
            std=self.std,
            algorithm=self.contrast_method,
        )

    @property
    def contrast_constant(self) -> float:
        """The contrast value times the bubble diameter."""
        DeprecationWarning(
            "The 'contrast_constant' property will be deprecated in a future release. Use .visibility instead."
        )
        return self.contrast * self.diameter

    @property
    def passed(self) -> bool:
        """Whether the disk ROI contrast passed."""
        return self.contrast > self.contrast_threshold

    @property
    def passed_visibility(self) -> bool:
        """Whether the disk ROI's visibility passed."""
        return self.visibility > self.visibility_threshold

    @property
    def passed_contrast_constant(self) -> bool:
        """Boolean specifying if ROI pixel value was within tolerance of the nominal value."""
        return self.contrast_constant > self.contrast_threshold

    @property
    def passed_cnr_constant(self) -> bool:
        """Boolean specifying if ROI pixel value was within tolerance of the nominal value."""
        return self.cnr_constant > self.cnr_threshold

    @property
    def plot_color(self) -> str:
        """Return one of two colors depending on if ROI passed."""
        return "green" if self.passed_visibility else "red"

    @property
    def plot_color_constant(self) -> str:
        """Return one of two colors depending on if ROI passed."""
        return "green" if self.passed_contrast_constant else "red"

    @property
    def plot_color_cnr(self) -> str:
        """Return one of two colors depending on if ROI passed."""
        return "green" if self.passed_cnr_constant else "red"

    def as_dict(self) -> dict:
        """Dump important data as a dictionary. Useful when exporting a `results_data` output"""
        return {
            "contrast method": self.contrast_method,
            "visibility": self.visibility,
            "visibility threshold": self.visibility_threshold,
            "passed visibility": self.passed_visibility,
            "contrast": self.contrast,
            "cnr": self.contrast_to_noise,
            "signal to noise": self.signal_to_noise,
        }

    def percentile(self, percentile: float) -> float:
        """Return the pixel value at the given percentile."""
        return np.percentile(self.circle_mask(), percentile)

    @cached_property
    def std(self) -> float:
        """The std within the ROI."""
        return float(np.std(self.circle_mask()))

    @cached_property
    def max(self) -> float:
        """The max pixel value of the ROI."""
        masked_img = self.circle_mask()
        return np.max(masked_img)

    @cached_property
    def min(self) -> float:
        """The min pixel value of the ROI."""
        masked_img = self.circle_mask()
        return np.min(masked_img)


class HighContrastDiskROI(DiskROI):
    """A class for analyzing the high-contrast disks."""

    contrast_threshold: float | None

    def __init__(
        self,
        array: np.ndarray,
        angle: float,
        roi_radius: float,
        dist_from_center: float,
        phantom_center: tuple | Point,
        contrast_threshold: float,
    ):
        """
        Parameters
        ----------
        contrast_threshold : float, int
            The threshold for considering a bubble to be "seen".
        """
        super().__init__(array, angle, roi_radius, dist_from_center, phantom_center)
        self.contrast_threshold = contrast_threshold

    def __repr__(self):
        return f"High-Contrast Disk; max pixel: {self.max}, min pixel: {self.min}"

    @cached_property
    def max(self) -> float:
        """The max pixel value of the ROI."""
        masked_img = self.circle_mask()
        return float(np.max(masked_img))

    @cached_property
    def min(self) -> float:
        """The min pixel value of the ROI."""
        masked_img = self.circle_mask()
        return float(np.min(masked_img))

    @cached_property
    def mean(self) -> float:
        """The mean pixel value of the ROI."""
        masked_img = self.circle_mask()
        return float(np.mean(masked_img))


class RectangleROI(Rectangle):
    """Class that represents a rectangular ROI."""

    def __init__(self, array, width, height, angle, dist_from_center, phantom_center):
        y_shift = np.sin(np.deg2rad(angle)) * dist_from_center
        x_shift = np.cos(np.deg2rad(angle)) * dist_from_center
        center = Point(phantom_center.x + x_shift, phantom_center.y + y_shift)
        # the ROI must be 'real', i.e. >= 2x2 matrix
        if width < 2:
            raise ValueError(f"The width must be >= 2. Given {width}")
        if height < 2:
            raise ValueError(f"The height must be >= 2. Given {height}")
        super().__init__(width, height, center, as_int=True)
        self._array = array

    def __repr__(self):
        return f"Rectangle ROI @ {self.center}; mean pixel: {self.pixel_value}"

    # TODO: See if I could use this somewhere
    # @classmethod
    # def from_regionprop(cls, regionprop: _RegionProperties, phan_center: Point):
    #     width = regionprop.bbox[3] - regionprop.bbox[1]
    #     height = regionprop.bbox[2] - regionprop.bbox[0]
    #     angle = np.rad2deg(np.arctan2((regionprop.centroid[0] - phan_center.y), (regionprop.centroid[1] - phan_center.x)))
    #     distance = phan_center.distance_to(Point(regionprop.centroid[1], regionprop.centroid[0]))
    #     return cls(regionprop.intensity_image, width=width, height=height,
    #                angle=angle, dist_from_center=distance, phantom_center=phan_center)

    @cached_property
    def pixel_array(self) -> np.ndarray:
        """The pixel array within the ROI."""
        return self._array[
            self.bl_corner.y : self.tr_corner.y, self.bl_corner.x : self.tr_corner.x
        ]

    @cached_property
    def pixel_value(self) -> float:
        """The pixel array within the ROI."""
        return float(np.mean(self.pixel_array))

    @cached_property
    def mean(self) -> float:
        """The mean value within the ROI."""
        return float(np.mean(self.pixel_array))

    @cached_property
    def std(self) -> float:
        """The std within the ROI."""
        return float(np.std(self.pixel_array))

    @cached_property
    def min(self) -> float:
        """The min value within the ROI."""
        return float(np.min(self.pixel_array))

    @cached_property
    def max(self) -> float:
        """The max value within the ROI."""
        return float(np.max(self.pixel_array))
