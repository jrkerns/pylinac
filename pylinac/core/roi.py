
from functools import lru_cache
from typing import Union, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle as mpl_Circle

from skimage.measure._regionprops import _RegionProperties

from .geometry import Circle, Point, Rectangle


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
    def __init__(self, array: np.ndarray, angle: Union[float, int], roi_radius: Union[float, int],
                 dist_from_center: Union[float, int], phantom_center: Union[Tuple, Point]):
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
    def _get_shifted_center(angle: Union[float, int], dist_from_center: Union[float, int], phantom_center: Point):
        """The center of the ROI; corrects for phantom dislocation and roll."""
        y_shift = np.sin(np.deg2rad(angle)) * dist_from_center
        x_shift = np.cos(np.deg2rad(angle)) * dist_from_center
        return Point(phantom_center.x + x_shift, phantom_center.y + y_shift)

    @property
    @lru_cache()
    def pixel_value(self) -> np.ndarray:
        """The median pixel value of the ROI."""
        masked_img = self.circle_mask()
        return np.nanmedian(masked_img)

    @property
    def std(self) -> np.ndarray:
        """The standard deviation of the pixel values."""
        masked_img = self.circle_mask()
        return np.nanstd(masked_img)

    @lru_cache()
    def circle_mask(self) -> np.ndarray:
        """Return a mask of the image, only showing the circular ROI."""
        # http://scikit-image.org/docs/dev/auto_examples/plot_camera_numpy.html
        masked_array = np.copy(self._array).astype(np.float)
        l_x, l_y = self._array.shape[0], self._array.shape[1]
        X, Y = np.ogrid[:l_x, :l_y]
        outer_disk_mask = (X - self.center.y) ** 2 + (Y - self.center.x) ** 2 > self.radius ** 2
        masked_array[outer_disk_mask] = np.NaN
        return masked_array

    def plot2axes(self, axes=None, edgecolor: str='black', fill: bool=False):
        """Plot the Circle on the axes.

        Parameters
        ----------
        axes : matplotlib.axes.Axes
            An MPL axes to plot to.
        edgecolor : str
            The color of the circle.
        fill : bool
            Whether to fill the circle with color or leave hollow.
        """
        if axes is None:
            fig, axes = plt.subplots()
            axes.imshow(self._array)
        axes.add_patch(mpl_Circle((self.center.x, self.center.y), edgecolor=edgecolor, radius=self.radius, fill=fill))


class LowContrastDiskROI(DiskROI):
    """A class for analyzing the low-contrast disks."""
    contrast_threshold: Optional[float]
    cnr_threshold: Optional[float]
    background: Optional[float]

    def __init__(self, array, angle, roi_radius, dist_from_center, phantom_center, contrast_threshold=None, background=None,
                 cnr_threshold=None):
        """
        Parameters
        ----------
        contrast_threshold : float, int
            The threshold for considering a bubble to be "seen".
        """
        super().__init__(array, angle, roi_radius, dist_from_center, phantom_center)
        self.contrast_threshold = contrast_threshold
        self.cnr_threshold = cnr_threshold
        self.background = background

    @property
    def contrast_to_noise(self) -> float:
        """The contrast to noise ratio of the bubble: (Signal - Background)/Stdev."""
        return abs(self.pixel_value - self.background) / self.std

    @property
    def contrast(self) -> float:
        """The contrast of the bubble compared to background: (ROI - backg) / (ROI + backg)."""
        return abs((self.pixel_value - self.background) / (self.pixel_value + self.background))

    @property
    def cnr_constant(self) -> float:
        """The contrast-to-noise value times the bubble diameter."""
        return self.contrast_to_noise * self.diameter

    @property
    def contrast_constant(self) -> float:
        """The contrast value times the bubble diameter."""
        return self.contrast * self.diameter

    @property
    def passed(self) -> bool:
        """Whether the disk ROI contrast passed."""
        return self.contrast > self.contrast_threshold

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
        return 'blue' if self.passed else 'red'

    @property
    def plot_color_constant(self) -> str:
        """Return one of two colors depending on if ROI passed."""
        return 'blue' if self.passed_contrast_constant else 'red'

    @property
    def plot_color_cnr(self) -> str:
        """Return one of two colors depending on if ROI passed."""
        return 'blue' if self.passed_cnr_constant else 'red'


class HighContrastDiskROI(DiskROI):
    """A class for analyzing the high-contrast disks."""
    contrast_threshold: Optional[float]

    def __init__(self, array, angle, roi_radius, dist_from_center, phantom_center, contrast_threshold):
        """
        Parameters
        ----------
        contrast_threshold : float, int
            The threshold for considering a bubble to be "seen".
        """
        super().__init__(array, angle, roi_radius, dist_from_center, phantom_center)
        self.contrast_threshold = contrast_threshold

    @property
    def max(self) -> np.ndarray:
        """The max pixel value of the ROI."""
        masked_img = self.circle_mask()
        return np.nanmax(masked_img)

    @property
    def min(self) -> np.ndarray:
        """The min pixel value of the ROI."""
        masked_img = self.circle_mask()
        return np.nanmin(masked_img)


class RectangleROI(Rectangle):
    """Class that represents a rectangular ROI."""

    def __init__(self, array, width, height, angle, dist_from_center, phantom_center):
        y_shift = np.sin(np.deg2rad(angle)) * dist_from_center
        x_shift = np.cos(np.deg2rad(angle)) * dist_from_center
        center = Point(phantom_center.x + x_shift, phantom_center.y + y_shift)
        super().__init__(width, height, center, as_int=True)
        self._array = array

    @property
    def pixel_array(self) -> np.ndarray:
        """The pixel array within the ROI."""
        return self._array[self.bl_corner.x:self.tr_corner.x, self.bl_corner.y:self.tr_corner.y]