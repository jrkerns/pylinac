from __future__ import annotations

import typing
from functools import cached_property

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from skimage import draw
from skimage.draw import polygon, polygon2mask
from skimage.measure._regionprops import _RegionProperties

from .contrast import Contrast, contrast, michelson, ratio, rms, visibility, weber
from .decorators import lru_cache
from .geometry import Circle, Point, Rectangle

if typing.TYPE_CHECKING:
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
    """A class representing a disk-shaped Region of Interest."""

    @classmethod
    def from_phantom_center(
        cls,
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

        Notes
        -----
        Parameter names are different from the regular class constructor
        due to historical reasons and semi-backwards compatibility.
        """
        center = cls._get_shifted_center(angle, dist_from_center, phantom_center)
        return cls(array=array, center=center, radius=roi_radius)

    def __init__(
        self,
        array: np.ndarray,
        radius: float,
        center: Point,
    ):
        """
        Parameters
        ----------
        array : ndarray
            The 2D array representing the image the disk is on.
        radius : float
            The radius of the ROI in pixels
        center : Point
            The center of the Disk ROI.
        """
        super().__init__(center_point=center, radius=radius)
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
    def mean(self) -> float:
        """The mean value within the ROI."""
        return float(np.mean(self.circle_mask()))

    @cached_property
    def std(self) -> float:
        """The standard deviation of the pixel values."""
        masked_img = self.circle_mask()
        return float(np.std(masked_img))

    @cached_property
    def min(self) -> float:
        """The min value within the ROI."""
        return float(np.min(self.circle_mask()))

    @cached_property
    def max(self) -> float:
        """The max value within the ROI."""
        return float(np.max(self.circle_mask()))

    @lru_cache()
    def circle_mask(self) -> np.ndarray:
        """Return a mask of the image, only showing the circular ROI."""
        rr, cc = draw.disk(center=(self.center.y, self.center.x), radius=self.radius)
        return self._array[rr, cc]

    def masked_array(self) -> np.ndarray:
        """A 2D array the same shape as the underlying image array, with the pixels
        within the ROI set to their pixel values, and the rest set to nan.
        """
        shape = self._array.shape
        img = np.full(shape, np.nan, dtype=self._array.dtype)
        rr, cc = draw.disk(
            center=(self.center.y, self.center.x), radius=self.radius, shape=shape
        )
        img[rr, cc] = self._array[rr, cc]
        return img

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

    @classmethod
    def from_phantom_center(
        cls,
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
        contrast_threshold : float, int
            The threshold for considering a bubble to be "seen".
        contrast_reference : float, int
            The reference contrast value.
        cnr_threshold : float, int
            The threshold for the CNR constant.
        contrast_method : str
            The method to calculate contrast.
        visibility_threshold : float, int
            The threshold for the visibility. If not provided, no visibility is calculated.

        Notes
        -----
        Parameter names are different from the regular class constructor
        due to historical reasons and semi-backwards compatibility.
        """
        center = cls._get_shifted_center(angle, dist_from_center, phantom_center)
        return cls(
            array=array,
            radius=roi_radius,
            center=center,
            contrast_threshold=contrast_threshold,
            contrast_reference=contrast_reference,
            cnr_threshold=cnr_threshold,
            contrast_method=contrast_method,
            visibility_threshold=visibility_threshold,
        )

    def __init__(
        self,
        array: np.ndarray | ArrayImage,
        radius: float,
        center: Point,
        contrast_threshold: float | None = None,
        contrast_reference: float | None = None,
        cnr_threshold: float | None = None,
        contrast_method: str = Contrast.MICHELSON,
        visibility_threshold: float = 0.1,
    ):
        """
        Parameters
        ----------
        array : ndarray
            The 2D array representing the image the disk is on.
        radius : int, float
            The radius of the ROI.
        center : Point
            The ROI center location.
        contrast_threshold : float, int
            The threshold for considering a bubble to be "seen".
        contrast_reference : float, int
            The reference contrast value.
        cnr_threshold : float, int
            The threshold for the CNR constant.
        contrast_method : str
            The method to calculate contrast.
        visibility_threshold : float, int
            The threshold for the visibility to be considered passing.
        """
        super().__init__(array, radius, center=center)
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
            "passed visibility": bool(self.passed_visibility),
            "contrast": self.contrast,
            "cnr": self.contrast_to_noise,
            "signal to noise": self.signal_to_noise,
        }

    def percentile(self, percentile: float) -> float:
        """Return the pixel value at the given percentile."""
        return float(np.percentile(self.circle_mask(), percentile))


class HighContrastDiskROI(DiskROI):
    """A class for analyzing the high-contrast disks."""

    contrast_threshold: float | None

    @classmethod
    def from_phantom_center(
        cls,
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
        contrast_threshold : float, int
            The threshold for considering a bubble to be "seen".

        Notes
        -----
        Parameter names are different from the regular class constructor
        due to historical reasons and semi-backwards compatibility.
        """
        center = cls._get_shifted_center(angle, dist_from_center, phantom_center)
        return cls(
            array=array,
            radius=roi_radius,
            center=center,
            contrast_threshold=contrast_threshold,
        )

    def __init__(
        self,
        array: np.ndarray,
        radius: float,
        center: Point,
        contrast_threshold: float,
    ):
        """
        Parameters
        ----------
        array : ndarray
            The 2D array representing the image the disk is on.
        radius : int, float
            The radius of the ROI.
        center : Point
            The ROI center location.
        contrast_threshold : float, int
            The threshold for considering a bubble to be "seen".
        """
        super().__init__(array=array, radius=radius, center=center)
        self.contrast_threshold = contrast_threshold

    def __repr__(self):
        return f"High-Contrast Disk; max pixel: {self.max}, min pixel: {self.min}"


class RectangleROI(Rectangle):
    """Class that represents a rectangular ROI."""

    @classmethod
    def from_phantom_center(
        cls,
        array: np.ndarray,
        width: float,
        height: float,
        angle: float,
        dist_from_center: float,
        phantom_center: Point,
        rotation: float = 0.0,
    ):
        """
        Parameters
        ----------
        array : ndarray
            The 2D array representing the image the disk is on.
        width : float
            The width of the ROI in pixels.
        height : float
            The height of the ROI in pixels.
        angle : float
            The angle of the ROI in degrees from the phantom center.
        dist_from_center : float
            The distance of the ROI from the phantom center.
        phantom_center : tuple
            The location of the phantom center.
        rotation : float
            The rotation of the ROI itself in degrees. Defaults to 0.0.

            .. warning::

                This is separate from the angle parameter, which is the angle of the ROI from the phantom center.

        Notes
        -----
        Parameter names are different from the regular class constructor
        due to historical reasons and semi-backwards compatibility.
        """
        y_shift = np.sin(np.deg2rad(angle)) * dist_from_center
        x_shift = np.cos(np.deg2rad(angle)) * dist_from_center
        center = Point(phantom_center.x + x_shift, phantom_center.y + y_shift)
        return cls(
            array=array,
            width=width,
            height=height,
            center=center,
            rotation=rotation,
        )

    def __init__(
        self,
        array: np.ndarray,
        width: float,
        height: float,
        center: Point,
        rotation: float = 0.0,
    ):
        """
        Parameters
        ----------
        array : ndarray
            The 2D array representing the image the disk is on.
        width : float
            The width of the ROI in pixels.
        height : float
            The height of the ROI in pixels.
        center : Point
            The location of the ROI center.
        rotation : float
            The rotation of the ROI itself in degrees. Defaults to 0.0.

            .. warning::

                This is separate from the angle parameter, which is the angle of the ROI from the phantom center.

        """
        # the ROI must be 'real', i.e. >= 2x2 matrix
        if width < 2:
            raise ValueError(f"The width must be >= 2. Given {width}")
        if height < 2:
            raise ValueError(f"The height must be >= 2. Given {height}")
        super().__init__(width, height, center, rotation=rotation)
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

    def plotly_debug(self):
        """Plot the ROI against the image array along with highlighting of the pixels within the ROI.

        Useful for debugging and visualizing the ROI."""
        # largely copied from Image.plotly()
        fig = go.Figure()
        fig.add_heatmap(
            z=self._array.array,
            colorscale="gray",
            name="Image",
            showlegend=True,
            showscale=False,
        )
        fig.update_traces(showscale=True)
        fig.add_heatmap(
            z=self.masked_array,
            colorscale="Viridis",
            name="ROI pixels",
            showlegend=True,
            showscale=False,
        )
        fig.update_layout(
            xaxis_showticklabels=False,
            yaxis_showticklabels=False,
            # this inverts the y axis so 0 is at the top
            # note that this will cause later `range=(...)` calls to fail;
            # appears to be bug in plotly.
            yaxis_autorange="reversed",
            yaxis_scaleanchor="x",
            yaxis_constrain="domain",
            xaxis_scaleanchor="y",
            xaxis_constrain="domain",
            legend={"x": 0},
            showlegend=True,
        )
        self.plotly(fig=fig, name="ROI Outline", showlegend=True)
        fig.show()

    @cached_property
    def masked_array(self) -> np.ndarray:
        """A 2D array the same shape as the underlying image array,
        with the pixels within the ROI set to their pixel values, and the rest set to NaN.

        This can be useful for plotting. We can also calculate non-spatial statistics on this array,
        but requires using special functions like `np.nanmean()` or `np.nanstd()`.

        It's better to calculate statistics on the `pixels_flat` property, which is always available.
        """
        vertices_array = np.array(
            # polygon2mask expects the vertices in (row, col) format
            [v.as_array(("y", "x")) for v in self.vertices]
        )
        mask = polygon2mask(image_shape=self._array.shape, polygon=vertices_array)
        mask = mask.astype(
            self._array.dtype
        )  # convert from boolean  to the same dtype as the image
        mask[mask == 0] = np.nan  # set the mask to NaN where it is not part of the ROI
        image_mask = mask * self._array
        return image_mask

    @cached_property
    def pixels_flat(self) -> np.ndarray:
        """A flattened array of the pixel values within the ROI.

        This is always available, even for rotated ROIs. However,
        it is flattened, so spatial statistics cannot be calculated.
        """
        corners = np.array(
            [
                # note the -1; this is due to the numpy exclusive indexing in `.pixel_array`
                # to keep these two properties consistent.
                (self.bl_corner.x, self.bl_corner.y - 1),  # bottom-left
                (self.br_corner.x - 1, self.br_corner.y - 1),  # bottom-right
                (self.tr_corner.x - 1, self.tr_corner.y),  # top-right
                (self.tl_corner.x, self.tl_corner.y),  # top-left
            ]
        )
        row_coords = corners[:, 1]
        col_coords = corners[:, 0]
        rr, cc = polygon(r=row_coords, c=col_coords, shape=self._array.shape)
        return self._array[rr, cc]

    @cached_property
    def pixel_array(self) -> np.ndarray:
        """A 2D array of the pixel values within the ROI. Only available for non-rotated ROIs.

        Raises
        ------
        ValueError
            If the rotation is != 0, the array cannot be reshaped into a 2D array.
        """
        if self.rotation != 0:
            raise ValueError(
                "The pixel array cannot be reshaped into a 2D array when the rotation is not 0."
            )
        # note that numpy indexing is exclusive of the end index! This might be considered a bug,
        # but for historical compatibility, we keep it this way.
        return self._array[
            int(np.round(self.tl_corner.y)) : int(np.round(self.bl_corner.y)),
            int(np.round(self.bl_corner.x)) : int(np.round(self.br_corner.x)),
        ]

    @cached_property
    def pixel_value(self) -> float:
        """The pixel array within the ROI."""
        return float(np.mean(self.pixels_flat))

    @cached_property
    def mean(self) -> float:
        """The mean value within the ROI."""
        return float(np.mean(self.pixels_flat))

    @cached_property
    def std(self) -> float:
        """The std within the ROI."""
        return float(np.std(self.pixels_flat))

    @cached_property
    def min(self) -> float:
        """The min value within the ROI."""
        return float(np.min(self.pixels_flat))

    @cached_property
    def max(self) -> float:
        """The max value within the ROI."""
        return float(np.max(self.pixels_flat))
