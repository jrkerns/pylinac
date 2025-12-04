"""The CT module automatically analyzes DICOM images of a CatPhan 503, 504, 600, 604, 700, Quart DVT, or ACR phantoms acquired when doing CBCT or CT quality assurance.
It can load a folder or zip file that the images are in and automatically correct for translational and rotational errors.
It can analyze the HU regions and image scaling (CTP404), the high-contrast line pairs (CTP528) to calculate the modulation transfer function (MTF),
the HU uniformity (CTP486), and Low Contrast (CTP515) on the corresponding slices.

For ACR and Quart phantoms, the equivalent sections are analyzed where applicable even though each module does not have an explicit name.
Where intuitive similarities between the phantoms exist, the library usage is the same.

Features:

* **Automatic phantom registration** - Your phantom can be tilted, rotated, or translated--pylinac will automatically register the phantom.
* **Automatic testing of all major modules** - Major modules are automatically registered and analyzed.
* **Any scan protocol** - Scan your CatPhan with any protocol; even scan it in a regular CT scanner.
  Any field size or field extent is allowed.
"""

from __future__ import annotations

import abc
import io
import itertools
import os
import textwrap
import webbrowser
import zipfile
from collections.abc import Sequence
from functools import cached_property
from io import BytesIO
from os import path as osp
from pathlib import Path
from typing import BinaryIO, Callable

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from plotly import graph_objects as go
from py_linq import Enumerable
from pydantic import BaseModel, Field
from scipy import ndimage
from skimage import draw, filters, measure, segmentation
from skimage.measure._regionprops import RegionProperties
from skimage.transform import EuclideanTransform

from .core import image, pdf
from .core.contrast import Contrast
from .core.geometry import Line, Point
from .core.image import ArrayImage, DicomImageStack, ImageLike, z_position
from .core.io import get_url, retrieve_demo_file
from .core.mtf import MTF
from .core.nps import (
    average_power,
    max_frequency,
    noise_power_spectrum_1d,
    noise_power_spectrum_2d,
)
from .core.plotly_utils import add_title, add_vertical_line
from .core.profile import CollapsedCircleProfile, FWXMProfile
from .core.roi import DiskROI, LowContrastDiskROI, RectangleROI
from .core.utilities import QuaacDatum, QuaacMixin, ResultBase, ResultsDataMixin
from .core.warnings import capture_warnings

# The ramp angle ratio is from the Catphan manual ("Scan slice geometry" section)
# and represents the fact that the wire is at an oblique angle (23°), making it appear
# longer than it is if it were normal or perpendicular to the z (imaging) axis. This ratio
# fixes the length to represent it as if it were perpendicular to the imaging axis.
RAMP_ANGLE_RATIO = 0.42


# The values in comment are the HU ranges taken from the Catphan manuals
# The nominal values are approximately the mid point.
AIR = -1000  # -1046 : -986
LUNG_7112 = -868  #  -925 : -810
PMP = -196  #  -220 : -172
LDPE = -104  #  -121 :  -87
POLY = -47  #   -65 :  -29
WATER = 0  #    -7 :    7
ACRYLIC = 115  #    92 :  137
BONE_20 = 237  #   211 :  263
DELRIN = 365  #   344 :  387
BONE_50 = 725  #   667 :  783
TEFLON = 1000  #   941 : 1060


class ROIResult(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    name: str = Field(description="The region the ROI was sampled from.")
    value: float = Field(description="The measured HU value.")
    stdev: float = Field(description="The pixel value standard deviation of the ROI.")
    difference: float | None = Field(
        description="The difference between the measured and nominal values."
    )
    nominal_value: float | None = Field(description="The nominal HU value.")
    passed: bool | None = Field(description=" Whether the ROI passed.")


class CTP404Result(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    offset: int = Field(
        description="The offset of the module from the origin slice in mm."
    )
    low_contrast_visibility: float = Field(
        description="The low contrast visibility score.",
        title="Low Contrast Visibility",
    )
    thickness_passed: bool = Field(description="Whether the slice thickness passed.")
    measured_slice_thickness_mm: float = Field(
        description="The measured slice thickness in mm.",
        title="Measured Slice Thickness (mm)",
    )
    thickness_num_slices_combined: int = Field(
        description="The number of slices combined when measuring slice thickness."
    )
    geometry_passed: bool = Field(
        description="Whether the geometry test passed, using the 4 nodes."
    )
    avg_line_distance_mm: float = Field(
        description="The average distance between the 4 nodes in mm.",
        title="Average Line Distance (mm)",
    )
    line_distances_mm: list[float] = Field(
        description="A list of the individual distances between nodes."
    )
    hu_linearity_passed: bool = Field(
        description="Whether all the HU ROIs were within tolerance."
    )
    hu_tolerance: float = Field(
        description="The tolerance used for the HU linearity test."
    )
    hu_rois: dict[str, ROIResult] = Field(
        description="A dictionary of the HU ROIs and their values. The keys will be the material name such as ``Acrylic``."
    )


class CTP486Result(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    uniformity_index: float = Field(
        description="The uniformity index as defined in Equation 2 of Elstrom et al",
        title="Uniformity Index",
    )
    integral_non_uniformity: float = Field(
        description="The integral non-uniformity as defined in Equation 1 of Elstrom et al",
        title="Integral Non-Uniformity",
    )
    nps_avg_power: float = Field(
        description="The average power of the noise power spectrum."
    )
    nps_max_freq: float = Field(
        description="The most populous frequency of the noise power spectrum."
    )
    passed: bool = Field(description="Whether the uniformity test passed.")
    rois: dict[str, ROIResult] = Field(
        description="A dictionary of the uniformity ROIs and their values. The keys will be the region name such as ``Top``."
    )


class CTP515Result(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    cnr_threshold: float = Field(
        description="The contrast-to-noise ratio threshold used to determine if a low contrast ROI was 'seen'."
    )
    num_rois_seen: int = Field(
        description="The number of ROIs that were above the threshold; ie. 'seen'.",
        title="Number of low-contrast ROIs detected",
    )
    roi_settings: dict = Field(
        description="The settings used for each low contrast ROI. The key names are ``n`` where ``<n>`` is the size of the ROI."
    )
    roi_results: dict = Field(
        description="Results of the low contrast ROIs and their values. The keys will be the size of the ROI such as ``'9'``."
    )


class CTP528Result(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    start_angle_radians: float | None = Field(
        description="The angle where the circular profile started."
    )
    mtf_lp_mm: dict = Field(
        description="A dictionary from 10% to 90% resolution in steps of 10 of the MTF in lp/mm. E.g. ``'20': 0.748``."
    )
    roi_settings: dict[str, dict[str, int | float]] = Field(
        description="A dictionary of the settings used for each MTF ROI. The key names are ``region_<n>`` where ``<n>`` is the region number."
    )


class CatphanResult(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    catphan_model: str = Field(description="The model of CatPhan that was analyzed.")
    catphan_roll_deg: float = Field(
        description="The roll of the phantom in degrees.",
        title="Phantom Roll (\N{DEGREE SIGN})",
    )
    origin_slice: int = Field(
        description="The 'origin' slice number. For CatPhan, this is the center of the HU module."
    )
    num_images: int = Field(description="The number of images in the passed dataset.")
    ctp404: CTP404Result = Field(
        description="The results of the CTP404 (HU linearity, spacing) module.",
        title="CTP404 (HU Linearity)",
    )
    ctp486: CTP486Result | None = Field(
        description="The results of the CTP486 module (HU uniformity) module.",
        title="CTP486 (HU Uniformity)",
        default=None,
    )
    ctp528: CTP528Result | None = Field(
        description="The results of the CTP528 (spatial resolution) module.",
        title="CTP528 (Spatial Resolution)",
        default=None,
    )
    ctp515: CTP515Result | None = Field(
        description="The results of the CTP515 (Low contrast) module.",
        default=None,
        title="CTP515 (Low Contrast)",
    )


class SpatialResolutionROI(RectangleROI):
    """Rectangle ROI specially for CatPhan."""

    @property
    def plot_color(self):
        return "blue"


class HUDiskROI(DiskROI):
    """An HU ROI object. Represents a circular area measuring either HU sample (Air, Poly, ...)
    or HU uniformity (bottom, left, ...).
    """

    def __init__(
        self,
        array: np.ndarray | ArrayImage,
        angle: float,
        roi_radius: float,
        dist_from_center: float,
        phantom_center: tuple | Point,
        nominal_value: float | None = None,
        tolerance: float | None = None,
        background_mean: float | None = None,
        background_std: float | None = None,
    ):
        """
        Parameters
        ----------
        nominal_value
            The nominal pixel value of the HU ROI.
        tolerance
            The roi pixel value tolerance.
        """
        new_center = self._get_shifted_center(angle, dist_from_center, phantom_center)
        super().__init__(array, roi_radius, new_center)
        self.nominal_val = nominal_value
        self.tolerance = tolerance

    @property
    def value_diff(self) -> float:
        """The difference in HU between measured and nominal."""
        return self.pixel_value - self.nominal_val

    @property
    def passed(self) -> bool:
        """Boolean specifying if ROI pixel value was within tolerance of the nominal value."""
        if self.tolerance:
            return abs(self.value_diff) <= self.tolerance
        else:
            return True

    @property
    def plot_color(self) -> str:
        """Return one of two colors depending on if ROI passed."""
        return "green" if self.passed else "red"


class ThicknessROI(RectangleROI):
    """A rectangular ROI that measures the angled wire rod in the HU linearity slice which determines slice thickness."""

    @cached_property
    def long_profile(self) -> FWXMProfile:
        """The profile along the axis perpendicular to ramped wire."""
        img = image.load(self.pixel_array)
        img.filter(size=1, kind="gaussian")
        return FWXMProfile(values=img.array.max(axis=np.argmin(img.shape)))

    @cached_property
    def wire_fwhm(self) -> float:
        """The FWHM of the wire in pixels."""
        return self.long_profile.field_width_px

    @property
    def plot_color(self) -> str:
        """The plot color."""
        return "blue"


class Slice:
    """Base class for analyzing specific slices of a CBCT dicom set."""

    def __init__(
        self,
        catphan,
        slice_num: int | None = None,
        combine: bool = True,
        combine_method: str = "mean",
        num_slices: int = 0,
        clear_borders: bool = True,
        original_image: ImageLike | None = None,
    ):
        """
        Parameters
        ----------

        catphan : :class:`~pylinac.cbct.CatPhanBase` instance.
            The catphan instance.
        slice_num : int
            The slice number of the DICOM array desired. If None, will use the ``slice_num`` property of subclass.
        combine : bool
            If True, combines the slices +/- ``num_slices`` around the slice of interest to improve signal/noise.
        combine_method : {'mean', 'max'}
            How to combine the slices if ``combine`` is True.
        num_slices : int
            The number of slices on either side of the nominal slice to combine to improve signal/noise; only
            applicable if ``combine`` is True.
        clear_borders : bool
            If True, clears the borders of the image to remove any ROIs that may be present.
        original_image : :class:`~pylinac.core.image.Image` or None
            The array of the slice. This is a bolt-on parameter for optimization.
            Leaving as None is fine, but can increase analysis speed if 1) this image is passed and
            2) there is no combination of slices happening, which is most of the time.
        """
        if slice_num is not None:
            self.slice_num = slice_num
        if combine and num_slices > 0:
            array = combine_surrounding_slices(
                catphan.dicom_stack,
                self.slice_num,
                mode=combine_method,
                slices_plusminus=num_slices,
            )
        elif original_image is not None:
            array = original_image
        else:
            array = catphan.dicom_stack[self.slice_num].array
        self.image = image.load(array)
        self.catphan_size = catphan.catphan_size
        self.mm_per_pixel = catphan.mm_per_pixel
        self.clear_borders = clear_borders
        self.clip_in_localization = catphan.clip_in_localization
        if catphan._phantom_center_func:
            self._phantom_center_func = catphan._phantom_center_func

    @property
    def __getitem__(self, item):
        return self.image.array[item]

    @cached_property
    def phantom_roi(self) -> RegionProperties:
        """Get the Scikit-Image ROI of the phantom

        The image is analyzed to see if:
        1) the CatPhan is even in the image (if there were any ROIs detected)
        2) an ROI is within the size criteria of the catphan
        3) the ROI area that is filled compared to the bounding box area is close to that of a circle
        """
        # convert the slice to binary and label ROIs
        edges = filters.scharr(self.image.as_type(float))
        if np.max(edges) < 0.1:
            raise ValueError(
                "No edges were found in the image that look like the phantom"
            )
        # we clip the image to avoid issues where
        # very high or very low HU values cause
        # thresholding problems. E.g. a very high HU bb
        # can make the region detection not see the phantom
        # I can see this causing problems in the future if the
        # HU values are insanely off. This also causes issues
        # with MRI images, which aren't HU values, hence the flag.
        if self.clip_in_localization:
            clipped_arr = np.clip(self.image.array, a_min=-1000, a_max=1000)
        else:
            clipped_arr = self
        larr, regionprops, num_roi = get_regions(
            clipped_arr,
            fill_holes=True,
            threshold="otsu",
            clear_borders=self.clear_borders,
        )
        # check that there is at least 1 ROI
        if num_roi < 1 or num_roi is None:
            raise ValueError(
                f"The number of ROIs detected {num_roi} was not the number expected (1)"
            )
        catphan_region = sorted(
            regionprops, key=lambda x: np.abs(x.filled_area - self.catphan_size)
        )[0]
        is_too_large = self.catphan_size * 1.3 < catphan_region.filled_area
        is_too_small = catphan_region.filled_area < self.catphan_size / 1.3
        if is_too_small or is_too_large:
            raise ValueError("Unable to find ROI of expected size of the phantom")
        return catphan_region

    def is_phantom_in_view(self) -> bool:
        """Whether the phantom appears to be within the slice."""
        try:
            self.phantom_roi
            return True
        except ValueError:
            return False

    @property
    def phan_center(self) -> Point:
        """Determine the location of the center of the phantom."""
        x = self._phantom_center_func[0](self.slice_num)
        y = self._phantom_center_func[1](self.slice_num)
        return Point(x=x, y=y)


class CatPhanModule(Slice):
    """Base class for a CTP module."""

    common_name: str = ""
    combine_method: str = "mean"
    num_slices: int = 0
    roi_settings: dict = {}
    background_roi_settings: dict = {}
    roi_dist_mm = float
    roi_radius_mm = float
    rois: dict = {}  # dicts of HUDiskROIs
    background_rois: dict = {}  # dict of HUDiskROIs; possibly empty
    window_min: int | None = None  # plt visualization
    window_max: int | None = None  # plt visualization

    def __init__(
        self,
        catphan,
        tolerance: float | None = None,
        offset: int = 0,
        clear_borders: bool = True,
    ):
        self.model = ""
        self._offset = offset
        self.origin_slice = catphan.origin_slice
        self.tolerance = tolerance
        self.slice_thickness = catphan.dicom_stack.metadata.SliceThickness
        self.slice_spacing = catphan.dicom_stack.slice_spacing
        self.catphan_roll = catphan.catphan_roll
        self.roi_size_factor = catphan.roi_size_factor
        self.scaling_factor = catphan.scaling_factor
        self.roll_slice_offset = catphan.roll_slice_offset
        self.mm_per_pixel = catphan.mm_per_pixel
        self.rois: dict[str, HUDiskROI] = {}
        self.background_rois: dict[str, HUDiskROI] = {}
        Slice.__init__(
            self,
            catphan,
            combine_method=self.combine_method,
            num_slices=self.num_slices,
            clear_borders=clear_borders,
        )
        self._convert_units_in_settings()
        self.preprocess(catphan)
        self._setup_rois()

    def _convert_units_in_settings(self) -> None:
        setting_groups = [
            getattr(self, attr) for attr in dir(self) if attr.endswith("roi_settings")
        ]
        for roi_settings in setting_groups:
            for roi, settings in roi_settings.items():
                if isinstance(settings, dict):
                    if settings.get("distance") is not None:
                        settings["distance_pixels"] = (
                            settings["distance"]
                            * self.scaling_factor
                            / self.mm_per_pixel
                        )
                    if settings.get("radial_distance") is not None:
                        settings["radial_distance_pixels"] = (
                            settings["radial_distance"]
                            * self.scaling_factor
                            / self.mm_per_pixel
                        )
                    if settings.get("transversal_distance") is not None:
                        settings["transversal_distance_pixels"] = (
                            settings["transversal_distance"]
                            * self.scaling_factor
                            / self.mm_per_pixel
                        )
                    if settings.get("angle") is not None:
                        settings["angle_corrected"] = (
                            settings["angle"] + self.catphan_roll
                        )
                    if settings.get("radius") is not None:
                        settings["radius_pixels"] = (
                            settings["radius"]
                            * self.roi_size_factor
                            / self.mm_per_pixel
                        )
                    if settings.get("width") is not None:
                        settings["width_pixels"] = (
                            settings["width"] * self.roi_size_factor / self.mm_per_pixel
                        )
                    if settings.get("height") is not None:
                        settings["height_pixels"] = (
                            settings["height"]
                            * self.roi_size_factor
                            / self.mm_per_pixel
                        )

    def preprocess(self, catphan):
        """A preprocessing step before analyzing the CTP module.

        Parameters
        ----------
        catphan : `~pylinac.cbct.CatPhanBase` instance.
        """
        pass

    @property
    def slice_num(self) -> int:
        """The slice number of the spatial resolution module.

        Returns
        -------
        float
        """
        return int(self.origin_slice + round(self._offset / self.slice_spacing))

    def _setup_rois(self) -> None:
        for name, setting in self.background_roi_settings.items():
            self.background_rois[name] = HUDiskROI(
                self.image,
                setting["angle_corrected"],
                setting["radius_pixels"],
                setting["distance_pixels"],
                self.phan_center,
            )
        if self.background_rois:
            background_mean = np.mean(
                [roi.pixel_value for roi in self.background_rois.values()]
            )
            background_std = np.std(
                [roi.pixel_value for roi in self.background_rois.values()]
            )
        else:
            background_mean = None
            background_std = None

        for name, setting in self.roi_settings.items():
            nominal_value = setting.get("value", 0)
            self.rois[name] = HUDiskROI(
                self.image,
                setting["angle_corrected"],
                setting["radius_pixels"],
                setting["distance_pixels"],
                self.phan_center,
                nominal_value,
                self.tolerance,
                background_mean=background_mean,
                background_std=background_std,
            )

    def plot_rois(self, axis: plt.Axes) -> None:
        """Plot the ROIs to the axis."""
        for roi in self.rois.values():
            roi.plot2axes(axis, edgecolor=roi.plot_color)
        for roi in self.background_rois.values():
            roi.plot2axes(axis, edgecolor="blue")

    def plotly_rois(self, fig: go.Figure) -> None:
        for name, roi in self.rois.items():
            roi.plotly(fig, line_color=roi.plot_color, name=name)
        for name, roi in self.background_rois.items():
            roi.plotly(fig, line_color="blue", name=f"{name} Background")

    def plot(self, axis: plt.Axes):
        """Plot the image along with ROIs to an axis"""
        self.image.plot(ax=axis, show=False, vmin=self.window_min, vmax=self.window_max)
        self.plot_rois(axis)
        axis.autoscale(tight=True)
        axis.set_title(f"{self.common_name} ({self.slice_num+1})")
        axis.axis("off")

    def plotly(self, **kwargs) -> go.Figure:
        """Plot the image along with the ROIs to a plotly figure."""
        fig = go.Figure()
        self.image.plotly(
            fig, show=False, zmin=self.window_min, zmax=self.window_max, **kwargs
        )
        self.plotly_rois(fig)
        add_title(fig, f"{self.common_name} ({self.slice_num+1})")
        return fig

    @property
    def roi_vals_as_str(self) -> str:
        return ", ".join(
            f"{name}: {roi.pixel_value}" for name, roi in self.rois.items()
        )


class CTP404CP504(CatPhanModule):
    """Class for analysis of the HU linearity, geometry, and slice thickness regions of the CTP404."""

    attr_name = "ctp404"
    common_name = "HU Linearity"
    roi_dist_mm = 58.7
    roi_radius_mm = 5
    roi_settings = {
        "Air": {
            "value": AIR,
            "angle": -90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "PMP": {
            "value": PMP,
            "angle": -120,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "LDPE": {
            "value": LDPE,
            "angle": 180,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Poly": {
            "value": POLY,
            "angle": 120,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Acrylic": {
            "value": ACRYLIC,
            "angle": 60,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Delrin": {
            "value": DELRIN,
            "angle": 0,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Teflon": {
            "value": TEFLON,
            "angle": -60,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
    }
    background_roi_settings = {
        "1": {"angle": -30, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "2": {"angle": -150, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "3": {"angle": -210, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "4": {"angle": 30, "distance": roi_dist_mm, "radius": roi_radius_mm},
    }
    # thickness
    thickness_roi_height = 40
    thickness_roi_width = 10
    thickness_roi_distance_mm = 38
    thickness_roi_settings = {
        "Left": {
            "angle": 180,
            "width": thickness_roi_width,
            "height": thickness_roi_height,
            "distance": thickness_roi_distance_mm,
        },
        "Bottom": {
            "angle": 90,
            "width": thickness_roi_height,
            "height": thickness_roi_width,
            "distance": thickness_roi_distance_mm,
        },
        "Right": {
            "angle": 0,
            "width": thickness_roi_width,
            "height": thickness_roi_height,
            "distance": thickness_roi_distance_mm,
        },
        "Top": {
            "angle": -90,
            "width": thickness_roi_height,
            "height": thickness_roi_width,
            "distance": thickness_roi_distance_mm,
        },
    }
    # geometry
    geometry_roi_size_mm = 35
    geometry_roi_settings = {
        "Top-Horizontal": (0, 1),
        "Bottom-Horizontal": (2, 3),
        "Left-Vertical": (0, 2),
        "Right-Vertical": (1, 3),
    }
    pad: str | int
    thickness_image: Slice

    def __init__(
        self,
        catphan,
        offset: int,
        hu_tolerance: float,
        thickness_tolerance: float,
        scaling_tolerance: float,
        clear_borders: bool = True,
        thickness_slice_straddle: str | int = "auto",
        expected_hu_values: dict[str, float | int] | None = None,
    ):
        """
        Parameters
        ----------
        catphan : `~pylinac.cbct.CatPhanBase` instance.
        offset : int
        hu_tolerance : float
        thickness_tolerance : float
        scaling_tolerance : float
        clear_borders : bool
        """
        self.mm_per_pixel = catphan.mm_per_pixel
        self.hu_tolerance = hu_tolerance
        self.thickness_tolerance = thickness_tolerance
        self.scaling_tolerance = scaling_tolerance
        self.thickness_rois = {}
        self.lines = {}
        self.thickness_slice_straddle = thickness_slice_straddle
        self.expected_hu_values = expected_hu_values
        super().__init__(
            catphan, tolerance=hu_tolerance, offset=offset, clear_borders=clear_borders
        )

    def preprocess(self, catphan) -> None:
        # for the thickness analysis image, combine thin slices or just use one slice if slices are thick
        if (
            isinstance(self.thickness_slice_straddle, str)
            and self.thickness_slice_straddle.lower() == "auto"
        ):
            if float(catphan.dicom_stack.metadata.SliceThickness) < 3.5:
                self.pad = 1
            else:
                self.pad = 0
        else:
            self.pad = self.thickness_slice_straddle
        self.thickness_image = Slice(
            catphan,
            combine_method="mean",
            num_slices=self.num_slices + self.pad,
            slice_num=self.slice_num,
            clear_borders=self.clear_borders,
        ).image

    def _replace_hu_values(self):
        """Possibly replace the HU values in the ROI settings with the expected values if the key is present."""
        if self.expected_hu_values is not None:
            for name, value in self.expected_hu_values.items():
                if name in self.roi_settings:
                    self.roi_settings[name]["value"] = value

    def _setup_rois(self) -> None:
        self._replace_hu_values()
        super()._setup_rois()
        self._setup_thickness_rois()
        if len(self.geometry_roi_settings) > 0:
            self._setup_geometry_rois()

    def _setup_thickness_rois(self) -> None:
        for name, setting in self.thickness_roi_settings.items():
            self.thickness_rois[name] = ThicknessROI.from_phantom_center(
                self.thickness_image,
                setting["width_pixels"],
                setting["height_pixels"],
                setting["angle_corrected"],
                setting["distance_pixels"],
                self.phan_center,
            )

    def _setup_geometry_rois(self) -> None:
        boxsize = self.geometry_roi_size_mm / self.mm_per_pixel
        xbounds = (int(self.phan_center.x - boxsize), int(self.phan_center.x + boxsize))
        ybounds = (int(self.phan_center.y - boxsize), int(self.phan_center.y + boxsize))
        geo_img = self.image[ybounds[0] : ybounds[1], xbounds[0] : xbounds[1]].copy()
        # clip to the nearest of the two extremes
        # this can arise from direct density scans. In that case the
        # 1 teflon node will not get detected as the edge intensity is much less than the other nodes (unlike normal)
        # So, we clip the sub-image to the nearest extreme to the median.
        # This does very little to normal scans. RAM-4056
        geo_img -= np.median(geo_img)
        nearest_extreme = min(abs(geo_img.max()), abs(geo_img.min()))
        geo_clipped = np.clip(geo_img, a_min=-nearest_extreme, a_max=nearest_extreme)
        geo_clipped_abs = np.abs(geo_clipped)
        larr, regionprops, num_roi = get_regions(
            geo_clipped_abs, fill_holes=True, clear_borders=False
        )
        # check that there is at least 1 ROI
        if num_roi < 4:
            raise ValueError("Unable to locate the Geometric nodes")
        elif num_roi > 4:
            regionprops = sorted(
                regionprops, key=lambda x: x.filled_area, reverse=True
            )[:4]
        sorted_regions = sorted(
            regionprops, key=lambda x: (2 * x.centroid[0] + x.centroid[1])
        )
        centers = [
            Point(
                r.weighted_centroid[1] + xbounds[0], r.weighted_centroid[0] + ybounds[0]
            )
            for r in sorted_regions
        ]
        #  setup the geometric lines
        for name, order in self.geometry_roi_settings.items():
            self.lines[name] = GeometricLine(
                centers[order[0]],
                centers[order[1]],
                self.mm_per_pixel,
                self.scaling_tolerance,
            )

    @property
    def lcv(self) -> float:
        """The low-contrast visibility"""
        return (
            2
            * abs(self.rois["LDPE"].pixel_value - self.rois["Poly"].pixel_value)
            / (self.rois["LDPE"].std + self.rois["Poly"].std)
        )

    def plotly_linearity(
        self, plot_delta: bool = True, show_legend: bool = True
    ) -> go.Figure:
        fig = go.Figure()
        nominal_x_values = [roi.nominal_val for roi in self.rois.values()]
        if plot_delta:
            values = [roi.value_diff for roi in self.rois.values()]
            nominal_measurements = [0] * len(values)
            ylabel = "HU Delta"
        else:
            values = [roi.pixel_value for roi in self.rois.values()]
            nominal_measurements = nominal_x_values
            ylabel = "Measured Values"
        fig.add_scatter(
            x=nominal_x_values,
            y=values,
            name="Measured values",
            mode="markers",
            marker=dict(color="green", size=10, symbol="cross", line=dict(width=1)),
        )
        fig.add_scatter(
            x=nominal_x_values,
            y=nominal_measurements,
            mode="lines",
            name="Nominal Values",
            marker_color="green",
        )
        fig.add_scatter(
            x=nominal_x_values,
            y=np.array(nominal_measurements) + self.hu_tolerance,
            mode="lines",
            name="Upper Tolerance",
            line=dict(dash="dash", color="red"),
        )
        fig.add_scatter(
            x=nominal_x_values,
            y=np.array(nominal_measurements) - self.hu_tolerance,
            mode="lines",
            name="Lower Tolerance",
            line=dict(dash="dash", color="red"),
        )
        fig.update_layout(
            xaxis_title="Nominal Values", yaxis_title=ylabel, showlegend=show_legend
        )
        add_title(fig, "HU Linearity")
        return fig

    def plot_linearity(
        self, axis: plt.Axes | None = None, plot_delta: bool = True
    ) -> tuple:
        """Plot the HU linearity values to an axis.

        Parameters
        ----------
        axis : None, matplotlib.Axes
            The axis to plot the values on. If None, will create a new figure.
        plot_delta : bool
            Whether to plot the actual measured HU values (False), or the difference from nominal (True).
        """
        nominal_x_values = [roi.nominal_val for roi in self.rois.values()]
        if axis is None:
            fig, axis = plt.subplots()
        if plot_delta:
            values = [roi.value_diff for roi in self.rois.values()]
            nominal_measurements = [0] * len(values)
            ylabel = "HU Delta"
        else:
            values = [roi.pixel_value for roi in self.rois.values()]
            nominal_measurements = nominal_x_values
            ylabel = "Measured Values"
        points = axis.plot(nominal_x_values, values, "g+", markersize=15, mew=2)
        axis.plot(nominal_x_values, nominal_measurements)
        axis.plot(
            nominal_x_values, np.array(nominal_measurements) + self.hu_tolerance, "r--"
        )
        axis.plot(
            nominal_x_values, np.array(nominal_measurements) - self.hu_tolerance, "r--"
        )
        axis.margins(0.05)
        axis.grid(True)
        axis.set_xlabel("Nominal Values")
        axis.set_ylabel(ylabel)
        axis.set_title("HU linearity")
        return points

    @property
    def passed_hu(self) -> bool:
        """Boolean specifying whether all the ROIs passed within tolerance."""
        return all(roi.passed for roi in self.rois.values())

    def plotly_rois(self, fig: go.Figure) -> None:
        super().plotly_rois(fig)
        # plot slice thickness / ramp ROIs
        for name, roi in self.thickness_rois.items():
            roi.plotly(fig, line_color="blue", name=f"Ramp {name}")
        # plot geometry lines
        for name, line in self.lines.items():
            line.plotly(fig, color=line.pass_fail_color, name=f"Geometry {name}")

    def plot_rois(self, axis: plt.Axes) -> None:
        """Plot the ROIs onto the image, as well as the background ROIs"""
        # plot HU linearity ROIs
        super().plot_rois(axis)
        # plot thickness ROIs
        for roi in self.thickness_rois.values():
            roi.plot2axes(axis, edgecolor="blue")
        # plot geometry lines
        for line in self.lines.values():
            line.plot2axes(axis, color=line.pass_fail_color)

    @property
    def passed_thickness(self) -> bool:
        """Whether the slice thickness was within tolerance from nominal."""
        return (
            self.slice_thickness - self.thickness_tolerance
            < self.meas_slice_thickness
            < self.slice_thickness + self.thickness_tolerance
        )

    @property
    def meas_slice_thickness(self) -> float:
        """The average slice thickness for the 4 wire measurements in mm."""
        return np.mean(
            sorted(
                roi.wire_fwhm * self.mm_per_pixel * RAMP_ANGLE_RATIO
                for roi in self.thickness_rois.values()
            )
        ) / (1 + 2 * self.pad)

    @property
    def avg_line_length(self) -> float:
        return float(np.mean([line.length_mm for line in self.lines.values()]))

    @property
    def passed_geometry(self) -> bool:
        """Returns whether all the line lengths were within tolerance."""
        return all(line.passed for line in self.lines.values())


class CTP404CP503(CTP404CP504):
    """Alias for namespace consistency"""

    pass


class CTP404CP600(CTP404CP504):
    roi_dist_mm = 58.7
    roi_radius_mm = 5
    roi_settings = {
        "Air": {
            "value": AIR,
            "angle": 90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "PMP": {
            "value": PMP,
            "angle": 60,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "LDPE": {
            "value": LDPE,
            "angle": 0,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Poly": {
            "value": POLY,
            "angle": -60,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Acrylic": {
            "value": ACRYLIC,
            "angle": -120,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Delrin": {
            "value": DELRIN,
            "angle": -180,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Teflon": {
            "value": TEFLON,
            "angle": 120,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Vial": {
            "value": WATER,
            "angle": -90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm
            - 1,  # the vial sits inside the ROI and needs some clearance
        },
    }

    def _setup_rois(self) -> None:
        """For the 600, the top ROI is an optional water vial slot. If the HU is near water we leave it, otherwise we remove it so as not to flag false failures"""
        super()._setup_rois()
        if self.rois["Vial"].pixel_value < -500:  # closer to air than water
            self.rois.pop("Vial")


class CTP404CP604(CTP404CP504):
    roi_dist_mm = 58.7
    roi_radius_mm = 5
    roi_settings = {
        "Air": {
            "value": AIR,
            "angle": -90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "PMP": {
            "value": PMP,
            "angle": -120,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "50% Bone": {
            "value": BONE_50,
            "angle": -150,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "LDPE": {
            "value": LDPE,
            "angle": 180,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Poly": {
            "value": POLY,
            "angle": 120,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Acrylic": {
            "value": ACRYLIC,
            "angle": 60,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "20% Bone": {
            "value": BONE_20,
            "angle": 30,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Delrin": {
            "value": DELRIN,
            "angle": 0,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Teflon": {
            "value": TEFLON,
            "angle": -60,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
    }
    background_roi_settings = {
        "1": {"angle": -30, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "2": {"angle": -210, "distance": roi_dist_mm, "radius": roi_radius_mm},
    }


class CTP404CP700(CTP404CP504):
    roi_dist_mm = 58.7
    roi_radius_mm = 5
    # The angles in the available scans are 180 - drawings.
    # If it is found that the slices need flipping some additional logic is needed, but for now it's hardcoded
    roi_settings = {
        "Air": {
            "value": AIR,
            "angle": 180 - -90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "PMP": {
            "value": PMP,
            "angle": 180 - -120,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Lung #7112": {
            "value": LUNG_7112,
            "angle": 180 - -165,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Delrin": {
            "value": DELRIN,
            "angle": 180 - 165,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Poly": {
            "value": POLY,
            "angle": 180 - 120,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Teflon": {
            "value": TEFLON,
            "angle": 180 - 90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Bone 20%": {
            "value": BONE_20,
            "angle": 180 - 60,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "LDPE": {
            "value": LDPE,
            "angle": 180 - 15,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Bone 50%": {
            "value": BONE_50,
            "angle": 180 - -15,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Acrylic": {
            "value": ACRYLIC,
            "angle": 180 - -60,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Vial": {
            "value": WATER,
            "angle": 180 - -135,
            "distance": 28,
            "radius": roi_radius_mm,
        },
    }

    background_roi_settings = {
        "1": {
            "angle": -37.5,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },  # NEE Between BONE_50 and ACRYLIC
        "2": {
            "angle": -142.5,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },  # NWW Between PMP and Lung_7112
        "3": {
            "angle": 142.5,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },  # SWW Between DELRIN and POLY
        "4": {
            "angle": 37.5,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },  # SEE Between BONE_20 and LDPE
    }
    # thickness
    # The CatPhan 700 has wire ramps on top and bottom only.
    # On the left and right there are bead ramps which require a different evaluation
    thickness_roi_height = 40
    thickness_roi_width = 10
    thickness_roi_distance_mm = 40
    thickness_roi_settings = {
        "Bottom": {
            "angle": 90,
            "width": thickness_roi_height,
            "height": thickness_roi_width,
            "distance": thickness_roi_distance_mm,
        },
        "Top": {
            "angle": -90,
            "width": thickness_roi_height,
            "height": thickness_roi_width,
            "distance": thickness_roi_distance_mm,
        },
    }
    # geometry
    geometry_roi_settings = {}

    def _setup_rois(self) -> None:
        """For the 700, the water vial is an optional slot. If the HU is near water we leave it, otherwise we remove it so as not to flag false failures"""
        super()._setup_rois()
        if self.rois["Vial"].pixel_value < -500:  # closer to air than water
            self.rois.pop("Vial")


class CTP486(CatPhanModule):
    """Class for analysis of the Uniformity slice of the CTP module. Measures 5 ROIs around the slice that
    should all be close to the same value.
    """

    attr_name = "ctp486"
    common_name = "HU Uniformity"
    roi_dist_mm = 53
    roi_radius_mm = 10
    nominal_value = 0
    nps_rois: dict[str, RectangleROI]
    roi_settings = {
        "Top": {
            "value": nominal_value,
            "angle": -90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Right": {
            "value": nominal_value,
            "angle": 0,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Bottom": {
            "value": nominal_value,
            "angle": 90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Left": {
            "value": nominal_value,
            "angle": 180,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Center": {
            "value": nominal_value,
            "angle": 0,
            "distance": 0,
            "radius": roi_radius_mm,
        },
    }

    def plot_profiles(self, axis: plt.Axes | None = None) -> None:
        """Plot the horizontal and vertical profiles of the Uniformity slice.

        Parameters
        ----------
        axis : None, matplotlib.Axes
            The axis to plot on; if None, will create a new figure.
        """
        if axis is None:
            fig, axis = plt.subplots()
        horiz_data = self.image[int(self.phan_center.y), :]
        vert_data = self.image[:, int(self.phan_center.x)]
        axis.plot(horiz_data, "g", label="Horizontal")
        axis.plot(vert_data, "b", label="Vertical")
        axis.autoscale(tight=True)
        axis.axhline(self.nominal_value + self.tolerance, color="r", linewidth=3)
        axis.axhline(self.nominal_value - self.tolerance, color="r", linewidth=3)
        axis.grid(True)
        axis.set_ylabel("HU")
        axis.legend(loc=8, fontsize="small", title="")
        axis.set_title("Uniformity Profiles")

    def _setup_rois(self) -> None:
        """Generate our NPS ROIs. They are just square versions of the existing ROIs."""
        super()._setup_rois()
        self.nps_rois = {}
        for name, setting in self.roi_settings.items():
            self.nps_rois[name] = RectangleROI.from_phantom_center(
                array=self.image,
                width=setting["radius_pixels"] * 2,
                height=setting["radius_pixels"] * 2,
                angle=setting["angle_corrected"],
                dist_from_center=setting["distance_pixels"],
                phantom_center=self.phan_center,
            )

    def plot(self, axis: plt.Axes):
        """Plot the ROIs but also the noise power spectrum ROIs"""
        for nps_roi in self.nps_rois.values():
            nps_roi.plot2axes(axis, edgecolor="green", linestyle="-.")
        super().plot(axis)

    def plotly(self, **kwargs) -> go.Figure:
        fig = super().plotly(**kwargs)
        for name, nps_roi in self.nps_rois.items():
            nps_roi.plotly(
                fig, line_color="green", line_dash="dash", name=f"NPS {name}"
            )
        return fig

    @property
    def overall_passed(self) -> bool:
        """Boolean specifying whether all the ROIs passed within tolerance."""
        return all(roi.passed for roi in self.rois.values())

    @property
    def uniformity_index(self) -> float:
        """The Uniformity Index. Elstrom et al equation 2. https://www.tandfonline.com/doi/pdf/10.3109/0284186X.2011.590525"""
        center = self.rois["Center"]
        uis = [
            100 * ((roi.pixel_value - center.pixel_value) / (center.pixel_value + 1000))
            for roi in self.rois.values()
        ]
        abs_uis = np.abs(uis)
        return uis[np.argmax(abs_uis)]

    @property
    def integral_non_uniformity(self) -> float:
        """The Integral Non-Uniformity. Elstrom et al equation 1. https://www.tandfonline.com/doi/pdf/10.3109/0284186X.2011.590525"""
        maxhu = max(roi.pixel_value for roi in self.rois.values())
        minhu = min(roi.pixel_value for roi in self.rois.values())
        return (maxhu - minhu) / (maxhu + minhu + 2000)

    @cached_property
    def power_spectrum_2d(self) -> np.ndarray:
        """The power spectrum of the uniformity ROI."""
        return noise_power_spectrum_2d(
            pixel_size=self.mm_per_pixel,
            rois=[r.pixel_array for r in self.nps_rois.values()],
        )

    @cached_property
    def power_spectrum_1d(self) -> np.ndarray:
        """The 1D power spectrum of the uniformity ROI."""
        return noise_power_spectrum_1d(self.power_spectrum_2d)

    @property
    def avg_noise_power(self) -> float:
        """The average noise power of the uniformity ROI."""
        return average_power(self.power_spectrum_1d)

    @property
    def max_noise_power_frequency(self) -> float:
        """The frequency of the maximum noise power. 0 means no pattern."""
        return max_frequency(self.power_spectrum_1d)


class CTP528(abc.ABC, CatPhanModule):
    """ABC for module detection purposes."""

    mtf: MTF


class CTP528CP504(CTP528):
    """Class for analysis of the Spatial Resolution slice of the CBCT dicom data set.

    A collapsed circle profile is taken of the line-pair region. This profile is search for
    peaks and valleys. The MTF is calculated from those peaks & valleys.

    Attributes
    ----------

    radius2linepairs_mm : float
        The radius in mm to the line pairs.

    """

    attr_name: str = "ctp528"
    common_name: str = "Spatial Resolution"
    radius2linepairs_mm = 47
    combine_method: str = "max"
    num_slices: int = 3
    boundaries: tuple[float, ...] = (
        0,
        0.107,
        0.173,
        0.236,
        0.286,
        0.335,
        0.387,
        0.434,
        0.479,
    )
    start_angle: float = np.pi
    ccw: bool = True
    roi_settings = {
        "region 1": {
            "start": boundaries[0],
            "end": boundaries[1],
            "num peaks": 2,
            "num valleys": 1,
            "peak spacing": 0.021,
            "gap size (cm)": 0.5,
            "lp/mm": 0.1,
        },
        "region 2": {
            "start": boundaries[1],
            "end": boundaries[2],
            "num peaks": 3,
            "num valleys": 2,
            "peak spacing": 0.01,
            "gap size (cm)": 0.25,
            "lp/mm": 0.2,
        },
        "region 3": {
            "start": boundaries[2],
            "end": boundaries[3],
            "num peaks": 4,
            "num valleys": 3,
            "peak spacing": 0.006,
            "gap size (cm)": 0.167,
            "lp/mm": 0.3,
        },
        "region 4": {
            "start": boundaries[3],
            "end": boundaries[4],
            "num peaks": 4,
            "num valleys": 3,
            "peak spacing": 0.00557,
            "gap size (cm)": 0.125,
            "lp/mm": 0.4,
        },
        "region 5": {
            "start": boundaries[4],
            "end": boundaries[5],
            "num peaks": 4,
            "num valleys": 3,
            "peak spacing": 0.004777,
            "gap size (cm)": 0.1,
            "lp/mm": 0.5,
        },
        "region 6": {
            "start": boundaries[5],
            "end": boundaries[6],
            "num peaks": 5,
            "num valleys": 4,
            "peak spacing": 0.00398,
            "gap size (cm)": 0.083,
            "lp/mm": 0.6,
        },
        "region 7": {
            "start": boundaries[6],
            "end": boundaries[7],
            "num peaks": 5,
            "num valleys": 4,
            "peak spacing": 0.00358,
            "gap size (cm)": 0.071,
            "lp/mm": 0.7,
        },
        "region 8": {
            "start": boundaries[7],
            "end": boundaries[8],
            "num peaks": 5,
            "num valleys": 4,
            "peak spacing": 0.0027866,
            "gap size (cm)": 0.063,
            "lp/mm": 0.8,
        },
    }

    def _setup_rois(self):
        pass

    def _convert_units_in_settings(self):
        pass

    @cached_property
    def mtf(self) -> MTF:
        """The Relative MTF of the line pairs, normalized to the first region.

        Returns
        -------
        dict
        """
        maxs = list()
        mins = list()
        for key, value in self.roi_settings.items():
            max_indices, max_values = self.circle_profile.find_peaks(
                min_distance=value["peak spacing"],
                max_number=value["num peaks"],
                search_region=(value["start"], value["end"]),
            )
            # check that the right number of peaks were found before continuing, otherwise stop searching for regions
            if len(max_values) != value["num peaks"]:
                break
            maxs.append(max_values.mean())
            _, min_values = self.circle_profile.find_valleys(
                min_distance=value["peak spacing"],
                max_number=value["num valleys"],
                search_region=(min(max_indices), max(max_indices)),
            )
            mins.append(min_values.mean())
        if not maxs:
            raise ValueError(
                "Did not find any spatial resolution pairs to analyze. File an issue on github (https://github.com/jrkerns/pylinac/issues) if this is a valid dataset."
            )

        spacings = [roi["lp/mm"] for roi in self.roi_settings.values()]
        mtf = MTF(lp_spacings=spacings, lp_maximums=maxs, lp_minimums=mins)
        return mtf

    @property
    def radius2linepairs(self) -> float:
        """Radius from the phantom center to the line-pair region, corrected for pixel spacing."""
        return self.radius2linepairs_mm * self.scaling_factor / self.mm_per_pixel

    def plotly_rois(self, fig: go.Figure) -> None:
        self.circle_profile.plotly(fig, color="blue", plot_peaks=False)
        fig.update_layout(
            showlegend=False,
        )

    def plot_rois(self, axis: plt.Axes) -> None:
        """Plot the circles where the profile was taken within."""
        self.circle_profile.plot2axes(axis, edgecolor="blue", plot_peaks=False)

    @cached_property
    def circle_profile(self) -> CollapsedCircleProfile:
        """Calculate the median profile of the Line Pair region.

        Returns
        -------
        :class:`pylinac.core.profile.CollapsedCircleProfile` : A 1D profile of the Line Pair region.
        """
        circle_profile = CollapsedCircleProfile(
            self.phan_center,
            self.radius2linepairs,
            image_array=self.image,
            start_angle=self.start_angle + np.deg2rad(self.catphan_roll),
            width_ratio=0.04 * self.roi_size_factor,
            sampling_ratio=2,
            ccw=self.ccw,
        )
        circle_profile.filter(0.001, kind="gaussian")
        circle_profile.ground()
        return circle_profile


class CTP528CP604(CTP528CP504):
    """Alias for namespace consistency."""

    pass


class CTP528CP600(CTP528CP504):
    start_angle = np.pi - 0.1
    ccw = False
    boundaries = (0, 0.127, 0.195, 0.255, 0.304, 0.354, 0.405, 0.453, 0.496)
    roi_settings = {
        "region 1": {
            "start": boundaries[0],
            "end": boundaries[1],
            "num peaks": 2,
            "num valleys": 1,
            "peak spacing": 0.021,
            "gap size (cm)": 0.5,
            "lp/mm": 0.1,
        },
        "region 2": {
            "start": boundaries[1],
            "end": boundaries[2],
            "num peaks": 3,
            "num valleys": 2,
            "peak spacing": 0.01,
            "gap size (cm)": 0.25,
            "lp/mm": 0.2,
        },
        "region 3": {
            "start": boundaries[2],
            "end": boundaries[3],
            "num peaks": 4,
            "num valleys": 3,
            "peak spacing": 0.006,
            "gap size (cm)": 0.167,
            "lp/mm": 0.3,
        },
        "region 4": {
            "start": boundaries[3],
            "end": boundaries[4],
            "num peaks": 4,
            "num valleys": 3,
            "peak spacing": 0.00557,
            "gap size (cm)": 0.125,
            "lp/mm": 0.4,
        },
        "region 5": {
            "start": boundaries[4],
            "end": boundaries[5],
            "num peaks": 4,
            "num valleys": 3,
            "peak spacing": 0.004777,
            "gap size (cm)": 0.1,
            "lp/mm": 0.5,
        },
        "region 6": {
            "start": boundaries[5],
            "end": boundaries[6],
            "num peaks": 5,
            "num valleys": 4,
            "peak spacing": 0.00398,
            "gap size (cm)": 0.083,
            "lp/mm": 0.6,
        },
        "region 7": {
            "start": boundaries[6],
            "end": boundaries[7],
            "num peaks": 5,
            "num valleys": 4,
            "peak spacing": 0.00358,
            "gap size (cm)": 0.071,
            "lp/mm": 0.7,
        },
        "region 8": {
            "start": boundaries[7],
            "end": boundaries[8],
            "num peaks": 5,
            "num valleys": 4,
            "peak spacing": 0.0027866,
            "gap size (cm)": 0.063,
            "lp/mm": 0.8,
        },
    }


class CTP528CP503(CTP528CP504):
    start_angle = 0
    ccw = False
    boundaries = (0, 0.111, 0.176, 0.240, 0.289, 0.339, 0.390, 0.436, 0.481)


class CTP528CP700(CTP528):
    """Class for analysis of the Spatial Resolution slice of the CatPhan 700 phantom

    Whereas the other modules use a collapsed circle profile, the CatPhan 700 uses a
    set of Rectangle ROIs. This is partly out of necessity as the ROIs are not
    radially-aligned. Because of this, we also don't inherit from other CTP528 module classes.
    """

    rois: dict[str, RectangleROI]
    attr_name: str = "ctp528"
    common_name: str = "Spatial Resolution"
    combine_method: str = "max"
    num_slices: int = 3
    start_angle = None  # this is for compatibility with other CTP528 classes where we use a circle profile.
    # We don't use a circle profile here but when serializing results, it's expected. We use None since it doesn't apply here.
    roi_settings = {
        # regions are from the top-left going clockwise on CT
        "region 1": {
            "lp/mm": 0.1,  # from the CatPhan 700 manual
            "radial_distance": 50,
            "transversal_distance": -7,
            "rotation": -90,
            "width": 3,
            "height": 11,
        },
        "region 2": {
            "lp/mm": 0.2,
            "radial_distance": 50,
            "transversal_distance": 11,
            "rotation": -90,
            "width": 3,
            "height": 11,
        },
        "region 3": {
            "lp/mm": 0.3,
            "radial_distance": 50,
            "transversal_distance": -5.5,
            "rotation": -45,
            "width": 3,
            "height": 10,
        },
        "region 4": {
            "lp/mm": 0.4,
            "radial_distance": 50,
            "transversal_distance": 9.5,
            "rotation": -45,
            "width": 3,
            "height": 8.5,
        },
        "region 5": {
            "lp/mm": 0.5,
            "radial_distance": 50,
            "transversal_distance": -9,
            "rotation": 0,
            "width": 3,
            "height": 8,
        },
        "region 6": {
            "lp/mm": 0.6,
            "radial_distance": 50,
            "transversal_distance": 2,
            "rotation": 0,
            "width": 3,
            "height": 7,
        },
        "region 7": {
            "lp/mm": 0.7,
            "radial_distance": 50,
            "transversal_distance": 12,
            "rotation": 0,
            "width": 3,
            "height": 6,
        },
        "region 8": {
            "lp/mm": 0.8,
            "radial_distance": 50,
            "transversal_distance": -10.5,
            "rotation": 45,
            "width": 3,
            "height": 4,
        },
    }

    def _setup_rois(self) -> None:
        tform_catphan_global = EuclideanTransform(
            rotation=np.deg2rad(self.catphan_roll),
            translation=[self.phan_center.x, self.phan_center.y],
        )
        for name, setting in self.roi_settings.items():
            # The roi's are placed in polar coordinates, i.e. rotation then intrinsic translation.
            # That's the same as extrinsic translation and then rotation which is implemented by the + operator
            tform_roi_catphan = EuclideanTransform(
                translation=[
                    setting["radial_distance_pixels"],
                    setting["transversal_distance_pixels"],
                ]
            ) + EuclideanTransform(rotation=np.deg2rad(setting["rotation"]))
            # The roi in global coordinates, is the placement of the roi wrt phantom and then the phantom wrt global.
            # This is implemented by the + operator
            tform_roi_global = tform_roi_catphan + tform_catphan_global
            self.rois[name] = SpatialResolutionROI(
                array=self.image.array,
                width=setting["width_pixels"],
                height=setting["height_pixels"],
                center=tform_roi_global.translation,
                rotation=np.rad2deg(tform_roi_global.rotation),
            )

    @cached_property
    def mtf(self) -> MTF:
        """The Relative MTF of the line pairs, normalized to the first region."""
        return MTF.from_high_contrast_diskset(
            spacings=[(r["lp/mm"]) for r in self.roi_settings.values()],
            diskset=self.rois.values(),
        )


class GeometricLine(Line):
    """Represents a line connecting two nodes/ROIs on the Geometry Slice.

    Attributes
    ----------
    nominal_length_mm : int, float
        The nominal distance between the geometric nodes, in mm.
    """

    nominal_length_mm: float | int = 50

    def __init__(
        self,
        geo_roi1: Point,
        geo_roi2: Point,
        mm_per_pixel: float,
        tolerance: int | float,
    ):
        """
        Parameters
        ----------
        geo_roi1 : GEO_ROI
            One of two ROIs representing one end of the line.
        geo_roi2 : GEO_ROI
            The other ROI which is the other end of the line.
        mm_per_pixel : float
            The mm/pixel value.
        tolerance : int, float
            The tolerance of the geometric line, in mm.
        """
        super().__init__(geo_roi1, geo_roi2)
        self.mm_per_pixel = mm_per_pixel
        self.tolerance = tolerance

    @property
    def passed(self) -> bool:
        """Whether the line passed tolerance."""
        return (
            self.nominal_length_mm - self.tolerance
            < self.length_mm
            < self.nominal_length_mm + self.tolerance
        )

    @property
    def pass_fail_color(self) -> str:
        """Plot color for the line, based on pass/fail status."""
        return "blue" if self.passed else "red"

    @property
    def length_mm(self) -> float:
        """Return the length of the line in mm."""
        return self.length * self.mm_per_pixel


class CTP515(CatPhanModule):
    """Class for analysis of the low contrast slice of the CTP module. Low contrast is measured by obtaining
    the average pixel value of the contrast ROIs and comparing that value to the average background value. To obtain
    a more "human" detection level, the contrast (which is largely the same across different-sized ROIs) is multiplied
    by the diameter. This value is compared to the contrast threshold to decide if it can be "seen".
    """

    attr_name = "ctp515"
    common_name = "Low Contrast"
    num_slices = 1
    roi_dist_mm = 50
    roi_radius_mm = [6, 3.5, 3, 2.5, 2, 1.5]
    roi_angles = [-87.4, -69.1, -52.7, -38.5, -25.1, -12.9]
    roi_settings = {
        "15": {
            "angle": roi_angles[0],
            "distance": roi_dist_mm,
            "radius": roi_radius_mm[0],
        },
        "9": {
            "angle": roi_angles[1],
            "distance": roi_dist_mm,
            "radius": roi_radius_mm[1],
        },
        "8": {
            "angle": roi_angles[2],
            "distance": roi_dist_mm,
            "radius": roi_radius_mm[2],
        },
        "7": {
            "angle": roi_angles[3],
            "distance": roi_dist_mm,
            "radius": roi_radius_mm[3],
        },
        "6": {
            "angle": roi_angles[4],
            "distance": roi_dist_mm,
            "radius": roi_radius_mm[4],
        },
        "5": {
            "angle": roi_angles[5],
            "distance": roi_dist_mm,
            "radius": roi_radius_mm[5],
        },
    }
    background_roi_dist_ratio = 0.75
    background_roi_radius_mm = 4
    WINDOW_SIZE = 50

    def __init__(
        self,
        catphan,
        tolerance: float,
        cnr_threshold: float,
        offset: int,
        contrast_method: str,
        visibility_threshold: float,
        clear_borders: bool = True,
    ):
        self.cnr_threshold = cnr_threshold
        self.contrast_method = contrast_method
        self.visibility_threshold = visibility_threshold
        super().__init__(
            catphan, tolerance=tolerance, offset=offset, clear_borders=clear_borders
        )

    def _setup_rois(self):
        # create both background rois dynamically, then create the actual sample ROI as normal
        for name, setting in self.roi_settings.items():
            self.background_rois[name + "-outer"] = (
                LowContrastDiskROI.from_phantom_center(
                    self.image,
                    setting["angle_corrected"],
                    self.background_roi_radius_mm / self.mm_per_pixel,
                    setting["distance_pixels"] * (2 - self.background_roi_dist_ratio),
                    self.phan_center,
                )
            )
            self.background_rois[name + "-inner"] = (
                LowContrastDiskROI.from_phantom_center(
                    self.image,
                    setting["angle_corrected"],
                    self.background_roi_radius_mm / self.mm_per_pixel,
                    setting["distance_pixels"] * self.background_roi_dist_ratio,
                    self.phan_center,
                )
            )
            background_val = float(
                np.mean(
                    [
                        self.background_rois[name + "-outer"].pixel_value,
                        self.background_rois[name + "-inner"].pixel_value,
                    ]
                )
            )

            self.rois[name] = LowContrastDiskROI.from_phantom_center(
                self.image,
                setting["angle_corrected"],
                setting["radius_pixels"],
                setting["distance_pixels"],
                self.phan_center,
                contrast_reference=background_val,
                cnr_threshold=self.cnr_threshold,
                contrast_method=self.contrast_method,
                visibility_threshold=self.visibility_threshold,
            )

    @property
    def rois_visible(self) -> int:
        """The number of ROIs "visible"."""
        return sum(roi.passed_visibility for roi in self.rois.values())

    @property
    def window_min(self) -> float:
        """Lower bound of CT window/leveling to show on the plotted image. Improves apparent contrast."""
        return (
            Enumerable(self.background_rois.values()).min(lambda r: r.pixel_value)
            - self.WINDOW_SIZE
        )

    @property
    def window_max(self) -> float:
        """Upper bound of CT window/leveling to show on the plotted image. Improves apparent contrast"""
        return (
            Enumerable(self.rois.values()).max(lambda r: r.pixel_value)
            + self.WINDOW_SIZE
        )


class CTP515CP600(CTP515):
    roi_angles = [
        -87.4 + 180,
        -69.1 + 180,
        -52.7 + 180,
        -38.5 + 180,
        -25.1 + 180,
        -12.9 + 180,
    ]
    roi_dist_mm = 50
    roi_radius_mm = [6, 3.5, 3, 2.5, 2, 1.5]
    roi_settings = {
        "15": {
            "angle": roi_angles[0],
            "distance": roi_dist_mm,
            "radius": roi_radius_mm[0],
        },
        "9": {
            "angle": roi_angles[1],
            "distance": roi_dist_mm,
            "radius": roi_radius_mm[1],
        },
        "8": {
            "angle": roi_angles[2],
            "distance": roi_dist_mm,
            "radius": roi_radius_mm[2],
        },
        "7": {
            "angle": roi_angles[3],
            "distance": roi_dist_mm,
            "radius": roi_radius_mm[3],
        },
        "6": {
            "angle": roi_angles[4],
            "distance": roi_dist_mm,
            "radius": roi_radius_mm[4],
        },
        "5": {
            "angle": roi_angles[5],
            "distance": roi_dist_mm,
            "radius": roi_radius_mm[5],
        },
    }


class CTP515CP700(CTP515CP600):
    """Alias for namespace consistency."""

    pass


class CatPhanBase(ResultsDataMixin[CatphanResult], QuaacMixin):
    """A class for loading and analyzing CT DICOM files of a CatPhan 504 & CatPhan 503. Can be from a CBCT or CT scanner
    Analyzes: Uniformity (CTP486), High-Contrast Spatial Resolution (CTP528), Image Scaling & HU Linearity (CTP404).
    """

    _demo_url: str = ""
    _model: str = ""
    air_bubble_radius_mm: int | float = 7
    localization_radius: int | float = 59
    was_from_zip: bool = False
    min_num_images = 39
    clear_borders: bool = True
    hu_origin_slice_variance = 400  # the HU variance required on the origin slice
    _phantom_center_func: tuple[Callable, Callable] | None = None
    modules: dict[CatPhanModule, dict[str, int]]
    dicom_stack: image.DicomImageStack | image.LazyDicomImageStack
    clip_in_localization: bool = False
    x_adjustment: float
    y_adjustment: float
    roi_size_factor: float
    scaling_factor: float
    angle_adjustment: float
    roll_slice_offset: float = 0

    def __init__(
        self,
        folderpath: str | Sequence[str] | Path | Sequence[Path] | Sequence[BytesIO],
        check_uid: bool = True,
        memory_efficient_mode: bool = False,
        is_zip: bool = False,
    ):
        """
        Parameters
        ----------
        folderpath : str, list of strings, or Path to folder
            String that points to the CBCT image folder location.
        check_uid : bool
            Whether to enforce raising an error if more than one UID is found in the dataset.
        memory_efficient_mode : bool
            Whether to use a memory efficient mode. If True, the DICOM stack will be loaded on demand rather than all at once.
            This will reduce the memory footprint but will be slower by ~25%. Default is False.

        Raises
        ------
        NotADirectoryError
            If folder str passed is not a valid directory.
        FileNotFoundError
            If no CT images are found in the folder
        """
        super().__init__()
        self.origin_slice = 0
        self.catphan_roll = 0
        if isinstance(folderpath, (str, Path)) and not is_zip:
            if not osp.isdir(folderpath):
                raise NotADirectoryError("Path given was not a Directory/Folder")
        if not memory_efficient_mode:
            stack = image.DicomImageStack
        elif memory_efficient_mode and is_zip:
            stack = image.LazyZipDicomImageStack
        else:
            stack = image.LazyDicomImageStack
        if is_zip:
            self.dicom_stack = stack.from_zip(
                folderpath, check_uid=check_uid, min_number=self.min_num_images
            )
        else:
            self.dicom_stack = stack(
                folderpath, check_uid=check_uid, min_number=self.min_num_images
            )

    @classmethod
    def from_demo_images(cls):
        """Construct a CBCT object from the demo images."""
        demo_file = retrieve_demo_file(name=cls._demo_url)
        return cls.from_zip(demo_file)

    @classmethod
    def from_url(cls, url: str, check_uid: bool = True):
        """Instantiate a CBCT object from a URL pointing to a .zip object.

        Parameters
        ----------
        url : str
            URL pointing to a zip archive of CBCT images.
        check_uid : bool
            Whether to enforce raising an error if more than one UID is found in the dataset.
        """
        filename = get_url(url)
        return cls.from_zip(filename, check_uid=check_uid)

    @classmethod
    def from_zip(
        cls,
        zip_file: str | zipfile.ZipFile | BinaryIO,
        check_uid: bool = True,
        memory_efficient_mode: bool = False,
    ):
        """Construct a CBCT object and pass the zip file.

        Parameters
        ----------
        zip_file : str, ZipFile
            Path to the zip file or a ZipFile object.
        check_uid : bool
            Whether to enforce raising an error if more than one UID is found in the dataset.
        memory_efficient_mode : bool
            Whether to use a memory efficient mode. If True, the DICOM stack will be loaded on demand rather than all at once.
            This will reduce the memory footprint but will be slower by ~25%. Default is False.

        Raises
        ------
        FileExistsError : If zip_file passed was not a legitimate zip file.
        FileNotFoundError : If no CT images are found in the folder
        """
        return cls(
            folderpath=zip_file,
            check_uid=check_uid,
            memory_efficient_mode=memory_efficient_mode,
            is_zip=True,
        )

    def plotly_analyzed_images(
        self,
        show: bool = True,
        show_colorbar: bool = True,
        show_legend: bool = True,
        **kwargs,
    ) -> dict[str, go.Figure]:
        """Plot the analyzed set of images to Plotly figures.


        Parameters
        ----------
        show : bool
            Whether to show the plot.
        show_colorbar : bool
            Whether to show the colorbar on the plot.
        show_legend : bool
            Whether to show the legend on the plot.
        kwargs
            Additional keyword arguments to pass to the plot.

        Returns
        -------
        dict
            A dictionary of the Plotly figures where the key is the name of the
            image and the value is the figure.
        """
        figs = {}
        figs["CTP404"] = self.ctp404.plotly(
            show_legend=show_legend, show_colorbar=show_colorbar
        )
        figs["HU Linearity"] = self.ctp404.plotly_linearity(show_legend=show_legend)
        figs["Side View"] = self.plotly_side_view(show_legend=show_legend)
        if self._has_module(CTP486):
            figs["CTP486"] = self.ctp486.plotly(
                show_legend=show_legend, show_colorbar=show_colorbar
            )
        if self._has_module(CTP528):
            figs["CTP528"] = self.ctp528.plotly(
                show_legend=show_legend, show_colorbar=show_colorbar
            )
            figs["MTF"] = self.ctp528.mtf.plotly(show_legend=show_legend)
        if self._has_module(CTP515):
            figs["CTP515"] = self.ctp515.plotly(
                show_legend=show_legend, show_colorbar=show_colorbar
            )

        if show:
            for fig in figs.values():
                fig.show()
        return figs

    def plot_analyzed_image(self, show: bool = True, **plt_kwargs) -> None:
        """Plot the images used in the calculation and summary data.

        Parameters
        ----------
        show : bool
            Whether to plot the image or not.
        plt_kwargs : dict
            Keyword args passed to the plt.figure() method. Allows one to set things like figure size.
        """

        # set up grid and axes
        plt.figure(**plt_kwargs)
        grid_size = (2, 4)
        hu_ax = plt.subplot2grid(grid_size, (0, 1))
        self.ctp404.plot(hu_ax)
        hu_lin_ax = plt.subplot2grid(grid_size, (0, 2))
        self.ctp404.plot_linearity(hu_lin_ax)
        # plot side view w/ module locations
        side_ax = plt.subplot2grid(grid_size, (1, 2))
        self.plot_side_view(side_ax)
        # plot individual modules
        if self._has_module(CTP486):
            unif_ax = plt.subplot2grid(grid_size, (0, 0))
            self.ctp486.plot(unif_ax)
            unif_prof_ax = plt.subplot2grid(grid_size, (1, 3))
            self.ctp486.plot_profiles(unif_prof_ax)
        if self._has_module(CTP528):
            sr_ax = plt.subplot2grid(grid_size, (1, 0))
            self.ctp528.plot(sr_ax)
            mtf_ax = plt.subplot2grid(grid_size, (0, 3))
            self.ctp528.mtf.plot(mtf_ax)
        if self._has_module(CTP515):
            locon_ax = plt.subplot2grid(grid_size, (1, 1))
            self.ctp515.plot(locon_ax)

        # finish up
        plt.tight_layout()
        if show:
            plt.show()

    def save_analyzed_image(self, filename: str | Path | BinaryIO, **kwargs) -> None:
        """Save the analyzed summary plot.

        Parameters
        ----------
        filename : str, file object
            The name of the file to save the image to.
        kwargs :
            Any valid matplotlib kwargs.
        """
        self.plot_analyzed_image(show=False)
        plt.savefig(filename, **kwargs)

    def plot_analyzed_subimage(
        self,
        subimage: str = "hu",
        delta: bool = True,
        show: bool = True,
    ) -> plt.Figure | None:
        """Plot a specific component of the CBCT analysis.

        Parameters
        ----------
        subimage : {'hu', 'un', 'sp', 'lc', 'mtf', 'lin', 'prof', 'side'}
            The subcomponent to plot. Values must contain one of the following letter combinations.
            E.g. ``linearity``, ``linear``, and ``lin`` will all draw the HU linearity values.

            * ``hu`` draws the HU linearity image.
            * ``un`` draws the HU uniformity image.
            * ``sp`` draws the Spatial Resolution image.
            * ``lc`` draws the Low Contrast image (if applicable).
            * ``mtf`` draws the RMTF plot.
            * ``lin`` draws the HU linearity values. Used with ``delta``.
            * ``prof`` draws the HU uniformity profiles.
            * ``side`` draws the side view of the phantom with lines of the module locations.
        delta : bool
            Only for use with ``lin``. Whether to plot the HU delta or actual values.
        show : bool
            Whether to actually show the plot.
        """
        subimage = subimage.lower()
        fig, ax = plt.subplots()
        plt.axis("off")

        if "hu" in subimage:  # HU, GEO & thickness objects
            self.ctp404.plot(ax)
            plt.autoscale(tight=True)
        elif "un" in subimage:  # uniformity
            self.ctp486.plot(ax)
            plt.autoscale(tight=True)
        elif "sp" in subimage:  # SR objects
            self.ctp528.plot(ax)
            plt.autoscale(tight=True)
        elif "mtf" in subimage:
            plt.axis("on")
            self.ctp528.mtf.plot(ax)
        elif "lc" in subimage:
            if self._has_module(CTP515):
                self.ctp515.plot(ax)
                plt.autoscale(tight=True)
            else:
                return
        elif "lin" in subimage:
            plt.axis("on")
            self.ctp404.plot_linearity(ax, delta)
        elif "prof" in subimage:
            plt.axis("on")
            self.ctp486.plot_profiles(ax)
        elif "side" in subimage:
            ax = plt.gca()
            self.plot_side_view(ax)
        else:
            raise ValueError(f"Subimage parameter {subimage} not understood")

        if show:
            plt.show()
        return fig

    def save_analyzed_subimage(
        self,
        filename: str | BinaryIO,
        subimage: str = "hu",
        delta: bool = True,
        **kwargs,
    ) -> plt.Figure | None:
        """Save a component image to file.

        Parameters
        ----------
        filename : str, file object
            The file to write the image to.
        subimage : str
            See :meth:`~pylinac.ct.CatPhanBase.plot_analyzed_subimage` for parameter info.
        delta : bool
            Only for use with ``lin``. Whether to plot the HU delta or actual values.
        """
        fig = self.plot_analyzed_subimage(subimage, delta=delta, show=False)
        if fig:  # no fig if we plot low contrast
            plt.savefig(filename, **kwargs)
            if isinstance(filename, str):
                print(f"CatPhan subimage figure saved to {osp.abspath(filename)}")
            return fig

    def _results(self) -> None:
        """Helper function to spit out values that will be tested."""
        print(self.results())
        print(f"Phantom roll: {self.catphan_roll}")
        print(f"Origin slice: {self.origin_slice}")
        mtfs = {}
        for mtf in (95, 90, 80, 50, 30):
            mtfval = self.ctp528.mtf.relative_resolution(mtf)
            mtfs[mtf] = mtfval
        print(f"MTFs: {mtfs}")

    def localize(self, origin_slice: int | None) -> None:
        """Find the slice number of the catphan's HU linearity module and roll angle"""
        self._phantom_center_func = self.find_phantom_axis()
        if origin_slice is not None:
            self.origin_slice = origin_slice
        else:
            self.origin_slice = self.find_origin_slice()
        self.catphan_roll = self.find_phantom_roll() + self.angle_adjustment
        if origin_slice is None:
            self.origin_slice = self.refine_origin_slice(
                initial_slice_num=self.origin_slice
            )
        # now that we have the origin slice, ensure we have scanned all linked modules
        if not self._ensure_physical_scan_extent():
            raise ValueError(
                "The physical scan extent does not match the module configuration. "
                "This means not all modules were included in the scan. Rescan the phantom to include all"
                "relevant modules, or remove modules from the analysis."
            )

    def _module_offsets(self) -> list[float]:
        """A list of the module offsets. Used to confirm scan extent"""
        absolute_origin_position = self.dicom_stack[self.origin_slice].z_position
        return [
            absolute_origin_position + config["offset"]
            for config in self.modules.values()
        ]

    def _ensure_physical_scan_extent(self) -> bool:
        """Ensure that all the modules of the phantom have been scanned. If a CBCT isn't
        positioned correctly, some modules might not be included.

        It appears there can be rounding errors between the DICOM tag and the actual slice position. See RAM-2897.
        """
        z_positions = [z_position(m) for m in self.dicom_stack.metadatas]
        min_scan_extent_slice = round(min(z_positions), 1)
        max_scan_extent_slice = round(max(z_positions), 1)
        min_config_extent_slice = round(min(self._module_offsets()), 1)
        max_config_extent_slice = round(max(self._module_offsets()), 1)
        return (min_config_extent_slice >= min_scan_extent_slice) and (
            max_config_extent_slice <= max_scan_extent_slice
        )

    def find_phantom_axis(self) -> (Callable, Callable):
        """We fit all the center locations of the phantom across all slices to a 1D poly function instead of finding them individually for robustness.

        Normally, each slice would be evaluated individually, but the RadMachine jig gets in the way of
        detecting the HU module (🤦‍♂️). To work around that in a backwards-compatible way we instead
        look at all the slices and if the phantom was detected, capture the phantom center.
        ALL the centers are then fitted to a 1D poly function and passed to the individual slices.
        This way, even if one slice is messed up (such as because of the phantom jig), the poly function
        is robust to give the real center based on all the other properly-located positions on the other slices.
        """
        z = []
        center_x = []
        center_y = []
        for idx, img in enumerate(self.dicom_stack):
            slice = Slice(
                self,
                slice_num=idx,
                clear_borders=self.clear_borders,
                original_image=img,
            )
            if slice.is_phantom_in_view():
                roi = slice.phantom_roi
                z.append(idx)
                center_y.append(roi.centroid[0])
                center_x.append(roi.centroid[1])
        # clip to exclude any crazy values
        zs = np.array(z)
        center_xs = np.array(center_x) + self.x_adjustment
        center_ys = np.array(center_y) + self.y_adjustment
        # gives an absolute and relative range so tight ranges are all included
        # but extreme values are excluded. Sometimes the range is very tight
        # and thus percentiles are not a sure thing
        x_idxs = np.argwhere(
            np.isclose(np.median(center_xs), center_xs, atol=3, rtol=0.01)
        )
        y_idxs = np.argwhere(
            np.isclose(np.median(center_ys), center_ys, atol=3, rtol=0.01)
        )
        common_idxs = np.intersect1d(x_idxs, y_idxs)
        # fit to 1D polynomials; inspiration: https://stackoverflow.com/a/45351484
        # rcond should be explicitly passed. Started randomly failing in the pipe. v1.14.0 numpy release notes
        # say it should be explicitly passed. Value is arbitrary but small and tests pass.
        fit_zx = np.poly1d(
            np.polyfit(zs[common_idxs], center_xs[common_idxs], deg=1, rcond=0.00001)
        )
        fit_zy = np.poly1d(
            np.polyfit(zs[common_idxs], center_ys[common_idxs], deg=1, rcond=0.00001)
        )
        return fit_zx, fit_zy

    @property
    def mm_per_pixel(self) -> float:
        """The millimeters per pixel of the DICOM images."""
        return self.dicom_stack.metadata.PixelSpacing[0]

    def find_origin_slice(self) -> int:
        """Using a brute force search of the images, find the median HU linearity slice.

        This method walks through all the images and takes a collapsed circle profile where the HU
        linearity ROIs are. If the profile contains both low (<800) and high (>800) HU values and most values are the same
        (i.e. it's not an artifact), then
        it can be assumed it is an HU linearity slice. The median of all applicable slices is the
        center of the HU slice.

        Returns
        -------
        int
            The middle slice of the HU linearity module.
        """
        hu_slices = []
        for image_number in range(0, self.num_images, 2):
            slice = Slice(
                self, image_number, combine=False, clear_borders=self.clear_borders
            )
            # print(image_number)
            # slice.image.plot()
            if slice.is_phantom_in_view():
                circle_prof = CollapsedCircleProfile(
                    slice.phan_center,
                    radius=self.localization_radius / self.mm_per_pixel,
                    image_array=slice.image,
                    width_ratio=0.05,
                    num_profiles=5,
                )
                prof = circle_prof.values
                # determine if the profile contains both low and high values and that most values are the same
                low_end, high_end = np.percentile(prof, [2, 98])
                median = np.median(prof)
                middle_variation = np.percentile(prof, 80) - np.percentile(prof, 20)
                variation_limit = max(
                    100, self.dicom_stack.metadata.SliceThickness * -100 + 300
                )
                if (
                    (low_end < median - self.hu_origin_slice_variance)
                    and (high_end > median + self.hu_origin_slice_variance)
                    and (middle_variation < variation_limit)
                ):
                    hu_slices.append(image_number)

        if not hu_slices:
            raise ValueError(
                "No slices were found that resembled the HU linearity module"
            )
        hu_slices = np.array(hu_slices)
        c = int(round(float(np.median(hu_slices))))
        ln = len(hu_slices)
        # drop slices that are way far from median
        hu_slices = hu_slices[((c + ln / 2) >= hu_slices) & (hu_slices >= (c - ln / 2))]
        center_hu_slice = int(round(float(np.median(hu_slices))))
        if self._is_within_image_extent(center_hu_slice):
            return center_hu_slice

    def refine_origin_slice(self, initial_slice_num: int) -> int:
        """Apply a refinement to the origin slice. This was added to handle
        the catphan 604 at least due to variations in the length of the HU plugs."""
        return initial_slice_num

    def _is_right_area(self, region: RegionProperties):
        thresh = np.pi * ((self.air_bubble_radius_mm / self.mm_per_pixel) ** 2)
        return thresh * 2 > region.filled_area > thresh / 2

    def _is_right_eccentricity(self, region: RegionProperties):
        return region.eccentricity < 0.5

    def find_phantom_roll(self, func: Callable | None = None) -> float:
        """Determine the "roll" of the phantom.

        This algorithm uses the two air bubbles in the HU slice and the resulting angle between them.

        Parameters
        ----------
        func
            A callable to sort the air ROIs.

        Returns
        -------
        float : the angle of the phantom in **degrees**.
        """
        # get slice for phantom roll detection
        slice_offset = round(self.roll_slice_offset / self.dicom_stack.slice_spacing)
        slice_num = self.origin_slice + slice_offset
        slice = Slice(self, slice_num, clear_borders=self.clear_borders)
        # get edges and make ROIs from it
        larr, regions, _ = get_regions(slice)
        # find appropriate ROIs and grab the two most centrally positioned ones
        hu_bubbles = [
            r
            for r in regions
            if (self._is_right_area(r) and self._is_right_eccentricity(r))
        ]
        func = func or (lambda x: abs(x.centroid[1] - slice.phan_center.x))
        central_bubbles = sorted(hu_bubbles, key=func)[:2]
        sorted_bubbles = sorted(
            central_bubbles, key=lambda x: x.centroid[0]
        )  # top, bottom
        if not sorted_bubbles:
            raise ValueError(
                "No air bubbles were found in the HU slice. The origin slice algorithm likely failed or the origin slice was passed and is incorrect."
            )
        y_dist = sorted_bubbles[1].centroid[0] - sorted_bubbles[0].centroid[0]
        x_dist = sorted_bubbles[1].centroid[1] - sorted_bubbles[0].centroid[1]
        phan_roll = np.arctan2(y_dist, x_dist)
        anglroll = np.rad2deg(phan_roll) - 90
        return anglroll

    @property
    def num_images(self) -> int:
        """The number of images loaded."""
        return len(self.dicom_stack)

    def _is_within_image_extent(self, image_num: int) -> bool:
        """Determine if the image number is beyond the edges of the images (negative or past last image)."""
        if self.num_images - 1 > image_num > 1:
            return True
        else:
            raise ValueError(
                "The determined image number is beyond the image extent. Either the entire dataset "
                "wasn't loaded or the entire phantom wasn't scanned."
            )

    @property
    def catphan_size(self) -> float:
        """The expected size of the phantom in pixels, based on a 20cm wide phantom."""
        phan_area = np.pi * (self.catphan_radius_mm**2)
        return phan_area / (self.mm_per_pixel**2)

    def publish_pdf(
        self,
        filename: str | Path,
        notes: str | None = None,
        open_file: bool = False,
        metadata: dict | None = None,
        logo: Path | str | None = None,
    ) -> None:
        """Publish (print) a PDF containing the analysis and quantitative results.

        Parameters
        ----------
        filename : (str, file-like object}
            The file to write the results to.
        notes : str, list of strings
            Text; if str, prints single line.
            If list of strings, each list item is printed on its own line.
        open_file : bool
            Whether to open the file using the default program after creation.
        metadata : dict
            Extra data to be passed and shown in the PDF. The key and value will be shown with a colon.
            E.g. passing {'Author': 'James', 'Unit': 'TrueBeam'} would result in text in the PDF like:
            --------------
            Author: James
            Unit: TrueBeam
            --------------
        logo: Path, str
            A custom logo to use in the PDF report. If nothing is passed, the default pylinac logo is used.
        """
        analysis_title = f"CatPhan {self._model} Analysis"
        module_images = [("hu", "lin")]
        if self._has_module(CTP528):
            module_images.append(("sp", "mtf"))
        if self._has_module(CTP486):
            module_images.append(("un", "prof"))
        if self._has_module(CTP515):
            module_images.append(("lc", None))
        module_images.append(("side", None))

        self._publish_pdf(
            filename,
            metadata,
            notes,
            analysis_title,
            [*self.results(as_list=True), ""],
            module_images,
            logo,
        )
        if open_file:
            webbrowser.open(filename)

    def _publish_pdf(
        self,
        filename: str,
        metadata: dict | None,
        notes: str,
        analysis_title: str,
        texts: Sequence[str],
        imgs: Sequence[tuple[str, str]],
        logo: Path | str | None = None,
    ):
        canvas = pdf.PylinacCanvas(
            filename, page_title=analysis_title, metadata=metadata, logo=logo
        )
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 4.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 4))

        for page, ((img1, img2), text) in enumerate(zip(imgs, texts)):
            for img, offset in zip((img1, img2), (12, 2)):
                if img is not None:
                    data = io.BytesIO()
                    self.save_analyzed_subimage(data, img)
                    canvas.add_image(data, location=(4, offset), dimensions=(15, 10))
            canvas.add_text(text=text, location=(1.5, 23))
            canvas.add_new_page()
        canvas.finish()

    def _zip_images(self) -> None:
        """Compress the raw images into a ZIP archive and remove the uncompressed images."""
        zip_name = rf'{osp.dirname(self.dicom_stack[0].path)}\CBCT - {self.dicom_stack[0].date_created(format="%A, %I-%M-%S, %B %d, %Y")}.zip'
        with zipfile.ZipFile(zip_name, "w", compression=zipfile.ZIP_DEFLATED) as zfile:
            for img in self.dicom_stack:
                zfile.write(img.path, arcname=osp.basename(img.path))
        for img in self.dicom_stack:
            try:
                os.remove(img.path)
            except Exception:
                pass

    def plotly_side_view(self, show_legend: bool) -> go.Figure:
        fig = go.Figure()
        side_array = self.dicom_stack.side_view(axis=1)
        add_title(fig, "Side View")
        fig.add_heatmap(z=side_array, colorscale="gray", showscale=False)

        for module in self._detected_modules():
            add_vertical_line(
                fig,
                module.slice_num,
                width=3,
                color="blue",
                name=module.common_name,
            )
        fig.update_layout(showlegend=show_legend)
        return fig

    def plot_side_view(self, axis: Axes) -> None:
        """Plot a view of the scan from the side with lines showing detected module positions"""
        side_array = self.dicom_stack.side_view(axis=1)
        axis.set_yticks([])
        axis.set_title("Side View")
        axis.imshow(side_array, aspect="auto", cmap="gray", interpolation="none")
        for module in self._detected_modules():
            axis.axvline(module.slice_num)

    def _detected_modules(self) -> list[CatPhanModule]:
        """A list of the modules detected. Unlike _get_module, this returns the instances"""
        modules = [self.ctp404]
        if self._has_module(CTP515):
            modules.append(self.ctp515)
        if self._has_module(CTP486):
            modules.append(self.ctp486)
        if self._has_module(CTP528):
            modules.append(self.ctp528)
        return modules

    def analyze(
        self,
        hu_tolerance: int | float = 40,
        scaling_tolerance: int | float = 1,
        thickness_tolerance: int | float = 0.2,
        low_contrast_tolerance: int | float = 1,
        cnr_threshold: int | float = 15,
        zip_after: bool = False,
        contrast_method: str = Contrast.MICHELSON,
        visibility_threshold: float = 0.15,
        thickness_slice_straddle: str | int = "auto",
        expected_hu_values: dict[str, int | float] | None = None,
        x_adjustment: float = 0,
        y_adjustment: float = 0,
        angle_adjustment: float = 0,
        roi_size_factor: float = 1,
        scaling_factor: float = 1,
        origin_slice: int | None = None,
        roll_slice_offset: float = 0,
    ):
        """Single-method full analysis of CBCT DICOM files.

        Parameters
        ----------
        hu_tolerance : int
            The HU tolerance value for both HU uniformity and linearity.
        scaling_tolerance : float, int
            The scaling tolerance in mm of the geometric nodes on the HU linearity slice (CTP404 module).
        thickness_tolerance : float, int
            The tolerance of the thickness calculation in mm, based on the wire ramps in the CTP404 module.

            .. warning:: Thickness accuracy degrades with image noise; i.e. low mAs images are less accurate.

        low_contrast_tolerance : int
            The number of low-contrast bubbles needed to be "seen" to pass.
        cnr_threshold : float, int
            The threshold for "detecting" low-contrast image. See RTD for calculation info.

            .. deprecated:: 3.0

                Use visibility parameter instead.

        zip_after : bool
            If the CT images were not compressed before analysis and this is set to true, pylinac will compress
            the analyzed images into a ZIP archive.
        contrast_method
            The contrast equation to use. See :ref:`low_contrast_topic`.
        visibility_threshold
            The threshold for detecting low-contrast ROIs. Use instead of ``cnr_threshold``. Follows the Rose equation.
            See :ref:`visibility`.
        thickness_slice_straddle
            The number of extra slices **on each side** of the HU module slice to use for slice thickness determination.
            The rationale is that for thin slices the ramp FWHM can be very noisy. I.e. a 1mm slice might have a 100%
            variation with a low-mAs protocol. To account for this, slice thicknesses < 3.5mm have 1 slice added
            on either side of the HU module (so 3 total slices) and then averaged. The default is 'auto',
            which follows the above logic. Set to an integer to explicitly use a certain amount of padding. Typical
            values are 0, 1, and 2.

            .. warning:: This is the padding **on either side**. So a value of 1 => 3 slices, 2 => 5 slices, 3 => 7 slices, etc.

        expected_hu_values
            An optional dictionary of the expected HU values for the HU linearity module. The keys are the ROI names and the values
            are the expected HU values. If a key is not present or the parameter is None, the default values will be used.
        x_adjustment: float
            A fine-tuning adjustment to the detected x-coordinate of the phantom center. This will move the
            detected phantom position by this amount in the x-direction in mm. Positive values move the phantom to the right.
        y_adjustment: float
            A fine-tuning adjustment to the detected y-coordinate of the phantom center. This will move the
            detected phantom position by this amount in the y-direction in mm. Positive values move the phantom down.
        angle_adjustment: float
            A fine-tuning adjustment to the detected angle of the phantom. This will rotate the phantom by this amount in degrees.
            Positive values rotate the phantom clockwise.
        roi_size_factor: float
            A fine-tuning adjustment to the ROI sizes of the phantom. This will scale the ROIs by this amount.
            Positive values increase the ROI sizes. In contrast to the scaling adjustment, this
            adjustment effectively makes the ROIs bigger or smaller, but does not adjust their position.
        scaling_factor: float
            A fine-tuning adjustment to the detected magnification of the phantom. This will zoom the ROIs and phantom outline (if applicable) by this amount.
            In contrast to the roi size adjustment, the scaling adjustment effectively moves the phantom and ROIs
            closer or further from the phantom center. I.e. this zooms the outline and ROI positions, but not ROI size.
        origin_slice : int, None
            The slice number of the HU linearity module. If None, the slice will be determined automatically. This is
            a fallback method if the automatic localization algorithm fails.
        roll_slice_offset : float
            The offset in mm from ``origin_slice`` used to select the slice for phantom
            roll detection. The phantom roll is determined based on the two
            inserts in the central vertical axis of the HU module, but this
            detection can be influenced by other structures. Adjusting this offset to
            select a different slice can enhance the accuracy of phantom roll detection.
        """
        self.x_adjustment = x_adjustment
        self.y_adjustment = y_adjustment
        self.angle_adjustment = angle_adjustment
        self.roi_size_factor = roi_size_factor
        self.scaling_factor = scaling_factor
        self.roll_slice_offset = roll_slice_offset
        self.localize(origin_slice)
        ctp404, offset = self._get_module(CTP404CP504, raise_empty=True)
        self.ctp404 = ctp404(
            self,
            offset=offset,
            hu_tolerance=hu_tolerance,
            thickness_tolerance=thickness_tolerance,
            scaling_tolerance=scaling_tolerance,
            clear_borders=self.clear_borders,
            thickness_slice_straddle=thickness_slice_straddle,
            expected_hu_values=expected_hu_values,
        )
        if self._has_module(CTP486):
            ctp486, offset = self._get_module(CTP486)
            self.ctp486 = ctp486(
                self,
                offset=offset,
                tolerance=hu_tolerance,
                clear_borders=self.clear_borders,
            )
        if self._has_module(CTP528):
            ctp528, offset = self._get_module(CTP528)
            self.ctp528 = ctp528(
                self, offset=offset, tolerance=None, clear_borders=self.clear_borders
            )
        if self._has_module(CTP515):
            ctp515, offset = self._get_module(CTP515)
            self.ctp515 = ctp515(
                self,
                tolerance=low_contrast_tolerance,
                cnr_threshold=cnr_threshold,
                offset=offset,
                contrast_method=contrast_method,
                visibility_threshold=visibility_threshold,
                clear_borders=self.clear_borders,
            )
        if zip_after and not self.was_from_zip:
            self._zip_images()

    def _has_module(self, module_of_interest: type[CatPhanModule]) -> bool:
        return any(
            issubclass(module, module_of_interest) for module in self.modules.keys()
        )

    def _get_module(
        self, module_of_interest: type[CatPhanModule], raise_empty: bool = False
    ) -> tuple[type[CatPhanModule], int]:
        """Grab the module that is, or is a subclass of, the module of interest. This allows users to subclass a CTP module and pass that in."""
        for module, values in self.modules.items():
            if issubclass(module, module_of_interest):
                return module, values.get("offset")
        if raise_empty:
            raise ValueError(
                f"Tried to find the {module_of_interest} or a subclass of it. Did you override `modules` and not pass this module in?"
            )

    def results(self, as_list: bool = False) -> str | list[list[str]]:
        """Return the results of the analysis as a string. Use with print().

        Parameters
        ----------
        as_list : bool
            Whether to return as a list of list of strings vs single string. Pretty much for internal usage.
        """
        results = []
        result = [
            f" - CBCT/CT {self._model} QA Test - ",
            " - CTP 404 Results - ",
            f"HU Linearity tolerance: {self.ctp404.hu_tolerance}",
            "HU Linearity ROIs:",
            # wrap so it doesn't fall off the page in PDFs
            *textwrap.wrap(self.ctp404.roi_vals_as_str, width=50),
            f"HU Passed?: {self.ctp404.passed_hu}",
            f"Low contrast visibility: {self.ctp404.lcv:2.2f}",
            f"Geometric Line Average (mm): {self.ctp404.avg_line_length:2.2f}",
            f"Geometry Passed?: {self.ctp404.passed_geometry}",
            f"Measured Slice Thickness (mm): {self.ctp404.meas_slice_thickness:2.3f}",
            f"Slice Thickness Passed? {self.ctp404.passed_thickness}",
        ]
        results.append(result)
        if self._has_module(CTP528):
            ctp528_result = [
                " - CTP528 Results - ",
                f"MTF 80% (lp/mm): {self.ctp528.mtf.relative_resolution(80):2.2f}",
                f"MTF 50% (lp/mm): {self.ctp528.mtf.relative_resolution(50):2.2f}",
                f"MTF 30% (lp/mm): {self.ctp528.mtf.relative_resolution(30):2.2f}",
            ]
            results.append(ctp528_result)
        if self._has_module(CTP486):
            ctp486_result = [
                " - CTP486 Results - ",
                f"Uniformity tolerance: {self.ctp486.tolerance}",
                f"Uniformity ROIs: {self.ctp486.roi_vals_as_str}",
                f"Uniformity index: {self.ctp486.uniformity_index:2.3f}",
                f"Integral non-uniformity: {self.ctp486.integral_non_uniformity:2.4f}",
                f"Uniformity Passed?: {self.ctp486.overall_passed}",
                f"Max Noise Power frequency: {self.ctp486.max_noise_power_frequency}",
                f"Average Noise Power: {self.ctp486.avg_noise_power}",
            ]
            results.append(ctp486_result)
        if self._has_module(CTP515):
            ctp515_result = [
                " - CTP515 Results - ",
                f"CNR threshold: {self.ctp515.cnr_threshold}",
                f'Low contrast ROIs "seen": {self.ctp515.rois_visible}',
            ]
            results.append(ctp515_result)
        if not as_list:
            result = "\n".join(itertools.chain(*results))
        else:
            result = results
        return result

    def _quaac_datapoints(self) -> dict[str, QuaacDatum]:
        results_data = self.results_data(as_dict=True)
        data = {
            "Phantom Roll": QuaacDatum(
                value=results_data["catphan_roll_deg"],
                unit="degrees",
            )
        }
        for material, roi in results_data["ctp404"]["hu_rois"].items():
            data[f"{material} HU"] = QuaacDatum(
                value=roi["value"], unit="HU", reference_value=roi["nominal_value"]
            )
        ctp404_keys = (
            ("avg_line_distance_mm", "Geometric Line Average", "mm"),
            ("measured_slice_thickness_mm", "Measured Slice Thickness", "mm"),
            ("low_contrast_visibility", "Low Contrast Visibility", ""),
        )
        for key, desc, unit in ctp404_keys:
            data[desc] = QuaacDatum(
                value=results_data["ctp404"][key],
                unit=unit,
            )
        if results_data["ctp486"] is not None:
            ctp486_keys = (
                ("uniformity_index", "Uniformity Index", ""),
                ("integral_non_uniformity", "Integral Non-Uniformity", ""),
            )
            for key, desc, unit in ctp486_keys:
                data[desc] = QuaacDatum(
                    value=results_data["ctp486"][key],
                    unit=unit,
                )
            for location, roi in results_data["ctp486"]["rois"].items():
                data[f"{location} Uniformity"] = QuaacDatum(
                    value=roi["value"], unit="HU", reference_value=roi["nominal_value"]
                )
        if results_data["ctp528"] is not None:
            for percent, mtf in results_data["ctp528"]["mtf_lp_mm"].items():
                data[f"MTF {percent}%"] = QuaacDatum(
                    value=mtf,
                    unit="lp/mm",
                )
        if results_data["ctp515"] is not None:
            for diameter, roi in results_data["ctp515"]["roi_results"].items():
                data[f"Low Contrast {diameter}mm visibility"] = QuaacDatum(
                    value=roi["visibility"],
                    unit="",
                )
                data[f"Low Contrast {diameter}mm CNR"] = QuaacDatum(
                    value=roi["cnr"],
                    unit="",
                )
                data[f"Low Contrast {diameter}mm SNR"] = QuaacDatum(
                    value=roi["signal to noise"],
                    unit="",
                )
                data[f"Low Contrast {diameter}mm Contrast"] = QuaacDatum(
                    value=roi["contrast"],
                    unit="",
                )
        return data

    def _generate_results_data(self) -> CatphanResult:
        """Present the results data and metadata as a dataclass or dict.
        The default return type is a dataclass."""
        ctp404_result = CTP404Result(
            offset=self.ctp404._offset,
            low_contrast_visibility=self.ctp404.lcv,
            thickness_passed=self.ctp404.passed_thickness,
            measured_slice_thickness_mm=self.ctp404.meas_slice_thickness,
            thickness_num_slices_combined=self.ctp404.num_slices + self.ctp404.pad,
            geometry_passed=self.ctp404.passed_geometry,
            avg_line_distance_mm=self.ctp404.avg_line_length,
            line_distances_mm=[
                line.length_mm for name, line in self.ctp404.lines.items()
            ],
            hu_linearity_passed=self.ctp404.passed_hu,
            hu_tolerance=self.ctp404.hu_tolerance,
            hu_rois=rois_to_results(self.ctp404.rois),
        )
        data = CatphanResult(
            catphan_model=self._model,
            catphan_roll_deg=self.catphan_roll,
            origin_slice=self.origin_slice,
            num_images=self.num_images,
            ctp404=ctp404_result,
        )

        # CTP 486 Uniformity stuff
        if self._has_module(CTP486):
            data.ctp486 = CTP486Result(
                passed=self.ctp486.overall_passed,
                uniformity_index=self.ctp486.uniformity_index,
                integral_non_uniformity=self.ctp486.integral_non_uniformity,
                rois=rois_to_results(self.ctp486.rois),
                nps_avg_power=self.ctp486.avg_noise_power,
                nps_max_freq=self.ctp486.max_noise_power_frequency,
            )

        # CTP 528 stuff
        if self._has_module(CTP528):
            data.ctp528 = CTP528Result(
                roi_settings=self.ctp528.roi_settings,
                start_angle_radians=self.ctp528.start_angle,
                mtf_lp_mm={
                    p: self.ctp528.mtf.relative_resolution(p) for p in range(10, 91, 10)
                },
            )

        # CTP 515 stuff
        if self._has_module(CTP515):
            data.ctp515 = CTP515Result(
                cnr_threshold=self.ctp515.cnr_threshold,
                num_rois_seen=self.ctp515.rois_visible,
                roi_settings=self.ctp515.roi_settings,
                roi_results={
                    key: roi.as_dict() for key, roi in self.ctp515.rois.items()
                },
            )
        return data


@capture_warnings
class CatPhan503(CatPhanBase):
    """A class for loading and analyzing CT DICOM files of a CatPhan 503.
    Analyzes: Uniformity (CTP486), High-Contrast Spatial Resolution (CTP528), Image Scaling & HU Linearity (CTP404).
    """

    _demo_url = "CatPhan503.zip"
    _model = "503"
    catphan_radius_mm = 97
    modules = {
        CTP404CP503: {"offset": 0},
        CTP486: {"offset": -110},
        CTP528CP503: {"offset": -30},
    }

    @staticmethod
    def run_demo(show: bool = True):
        """Run the CBCT demo using high-quality head protocol images."""
        cbct = CatPhan503.from_demo_images()
        cbct.analyze()
        print(cbct.results())
        cbct.plot_analyzed_image(show)


@capture_warnings
class CatPhan504(CatPhanBase):
    """A class for loading and analyzing CT DICOM files of a CatPhan 504. Can be from a CBCT or CT scanner
    Analyzes: Uniformity (CTP486), High-Contrast Spatial Resolution (CTP528),
    Image Scaling & HU Linearity (CTP404), and Low contrast (CTP515).
    """

    _demo_url = "CatPhan504.zip"
    _model = "504"
    catphan_radius_mm = 101
    modules = {
        CTP404CP504: {"offset": 0},
        CTP486: {"offset": -65},
        CTP528CP504: {"offset": 30},
        CTP515: {"offset": -30},
    }

    @staticmethod
    def run_demo(show: bool = True):
        """Run the CBCT demo using high-quality head protocol images."""
        cbct = CatPhan504.from_demo_images()
        cbct.analyze()
        print(cbct.results())
        cbct.plot_analyzed_image(show)


@capture_warnings
class CatPhan604(CatPhanBase):
    """A class for loading and analyzing CT DICOM files of a CatPhan 604. Can be from a CBCT or CT scanner
    Analyzes: Uniformity (CTP486), High-Contrast Spatial Resolution (CTP528),
    Image Scaling & HU Linearity (CTP404), and Low contrast (CTP515).
    """

    _demo_url = "CatPhan604.zip"
    _model = "604"
    catphan_radius_mm = 101
    modules = {
        CTP404CP604: {"offset": 0},
        CTP486: {"offset": -80},
        CTP528CP604: {"offset": 40},
        CTP515: {"offset": -40},
    }

    @staticmethod
    def run_demo(show: bool = True):
        """Run the CBCT demo using high-quality head protocol images."""
        cbct = CatPhan604.from_demo_images()
        cbct.analyze()
        print(cbct.results())
        cbct.plot_analyzed_image(show)

    def refine_origin_slice(self, initial_slice_num: int) -> int:
        """The HU plugs are longer than the 'wire section'. This applies a refinement to find the
        slice that has the least angle between the centers of the left and right wires.

        Under normal conditions, we would simply apply an offset to the
        initial slice. The rods extend 4-5mm past the wire section.
        Unfortunately, we also sometimes need to account for users with the
        RM R1-4 jig which will cause localization issues due to the base of the
        jig touching the phantom.
        This solution is robust to the jig being present, but can suffer
        from images where the angle of the wire, due to noise,
        appears small but doesn't actually represent the wire ramp.

        Starting with the initial slice, go +/- 5 slices to find the slice with the least angle
        between the left and right wires.

        This suffers from a weakness that the roll is not yet determined.
        This will thus return the slice that has the least ABSOLUTE
        roll. If the phantom has an inherent roll, this will not be detected and may be off by a slice or so.
        Given the angle of the wire, the error would be small and likely only 1-2 slices max.
        """
        angles = []
        # make a CTP module for the purpose of easily extracting the thickness ROIs
        ctp404, offset = self._get_module(CTP404CP604, raise_empty=True)
        # we don't want to set up our other ROIs (they sometimes fail) so we temporarily override the method
        original_setup = ctp404._setup_rois
        ctp404._setup_rois = lambda x: x
        ctp = ctp404(
            self,
            offset=offset,
            clear_borders=self.clear_borders,
            hu_tolerance=0,
            scaling_tolerance=0,
            thickness_tolerance=0,
        )
        # we have to reset the method after we're done for future calls
        ctp404._setup_rois = original_setup
        for slice_num in range(initial_slice_num - 5, initial_slice_num + 5):
            slice = Slice(self, slice_num, clear_borders=self.clear_borders)
            # make a slice and add the wire ROIs to it.
            troi = {}
            for name, setting in ctp.thickness_roi_settings.items():
                troi[name] = ThicknessROI.from_phantom_center(
                    slice.image,
                    setting["width_pixels"],
                    setting["height_pixels"],
                    setting["angle_corrected"],
                    setting["distance_pixels"],
                    slice.phan_center,
                )
            # now find the angle between the left and right and top and bottom wires via the long profile
            left_wire = troi["Left"].long_profile.center_idx
            right_wire = troi["Right"].long_profile.center_idx
            h_angle = abs(left_wire - right_wire)
            top_wire = troi["Top"].long_profile.center_idx
            bottom_wire = troi["Bottom"].long_profile.center_idx
            v_angle = abs(top_wire - bottom_wire)
            angle = (h_angle + v_angle) / 2

            angles.append(
                {
                    "slice": slice_num,
                    "angle": angle,
                    "left width": troi["Left"].long_profile.field_width_px,
                    "right width": troi["Right"].long_profile.field_width_px,
                    # the values AT the FWXM
                    "left center": troi["Left"].long_profile.y_at_x(left_wire),
                    "right center": troi["Right"].long_profile.y_at_x(right_wire),
                    "left profile": troi["Left"].long_profile.values,
                    "right profile": troi["Right"].long_profile.values,
                }
            )

        # some slices might not include the wire
        # we need to drop those; we do so by dropping pairs that have a field width well below the median
        # or by fields who don't appear to have the wire in them (max value is near median)
        median_width_l = np.median([angle["left width"] for angle in angles])
        median_width_r = np.median([angle["right width"] for angle in angles])
        median_width = (median_width_l + median_width_r) / 2
        # get median and max pixel values of all the profiles
        median_left_pixel_val = np.median(
            np.concatenate([a["left profile"] for a in angles])
        )
        median_right_pixel_val = np.median(
            np.concatenate([a["right profile"] for a in angles])
        )
        median_pixel_val = (median_left_pixel_val + median_right_pixel_val) / 2
        max_left_pixel_val = np.max(np.concatenate([a["left profile"] for a in angles]))
        max_right_pixel_val = np.max(
            np.concatenate([a["right profile"] for a in angles])
        )
        max_pixel_val = (max_left_pixel_val + max_right_pixel_val) / 2

        for angle_set in angles.copy():
            # field width is well below the median; probably not in the slice; drop it
            if (
                angle_set["left width"] < median_width * 0.7
                or angle_set["right width"] < median_width * 0.7
            ):
                angles.remove(angle_set)
                continue
            # if the max pixel value of the angle set is closer to the overall median than the max
            # it means the wire isn't in the slice; drop it
            fwxm_pixel = np.mean((angle_set["left center"], angle_set["right center"]))
            delta_median = abs(median_pixel_val - fwxm_pixel)
            delta_max = abs(max_pixel_val - fwxm_pixel)
            if delta_median < delta_max:
                angles.remove(angle_set)

        # now find the slice with the least angle, accounting for the phantom roll
        m_slice_num = np.argsort([a["angle"] - self.catphan_roll for a in angles])
        refined_slice = angles[m_slice_num[0]]["slice"]
        return refined_slice


@capture_warnings
class CatPhan600(CatPhanBase):
    """A class for loading and analyzing CT DICOM files of a CatPhan 600.
    Analyzes: Uniformity (CTP486), High-Contrast Spatial Resolution (CTP528),
    Image Scaling & HU Linearity (CTP404), and Low contrast (CTP515).
    """

    _demo_url = "CatPhan600.zip"
    _model = "600"
    catphan_radius_mm = 101
    modules = {
        CTP404CP600: {"offset": 0},
        CTP486: {"offset": -160},
        CTP515CP600: {"offset": -110},
        CTP528CP600: {"offset": -70},
    }

    @staticmethod
    def run_demo(show: bool = True):
        """Run the CatPhan 600 demo."""
        cbct = CatPhan600.from_demo_images()
        cbct.analyze()
        print(cbct.results())
        cbct.plot_analyzed_image(show)

    def find_phantom_roll(self, func: Callable | None = None) -> float:
        """With the CatPhan 600, we have to consider that the top air ROI
        has a water vial in it (see pg 12 of the manual). If so, the top air ROI won't be detected.
        Rather, the default algorithm will find the bottom air ROI and teflon to the left.
        It may also find the top air ROI if the water vial isn't there.
        We use the below lambda to select the bottom air and teflon ROIs consistently.
        These two ROIs are at 75 degrees from cardinal. We thus offset the default outcome by 75.

        HOWEVER, for direct density scans, the Teflon ROI might not register because of the reduced
        HU. 🤦‍♂️. We make a best guess depending on the detected roll. If it's ~75 degrees,
        we have caught the bottom Air and Teflon. If it's near zero, we have
        caught the top and bottom Air ROIs.
        """
        angle = super().find_phantom_roll(lambda x: -x.centroid[0])
        if abs(angle) < 10:
            return angle
        return angle + 75


@capture_warnings
class CatPhan700(CatPhanBase):
    """A class for loading and analyzing CT DICOM files of a CatPhan 700.
    Analyzes: Uniformity (CTP486/CTP712), High-Contrast Spatial Resolution (CTP714),
    Image Scaling & HU Linearity (CTP404CP700/CTP682), Low contrast (CTP515) and slice geometry (CTP721)
    """

    # _demo_url = "CatPhan600.zip"
    _model = "700"
    catphan_radius_mm = 101
    modules = {
        CTP404CP700: {"offset": 0},
        CTP515CP700: {"offset": -80},
        CTP486: {"offset": -160},
        CTP528CP700: {"offset": -40},
    }


def get_regions(
    slice_or_arr: Slice | np.array,
    fill_holes: bool = False,
    clear_borders: bool = True,
    threshold: str = "otsu",
) -> tuple[np.array, list, int]:
    """Get the skimage regions of a black & white image."""
    if threshold == "otsu":
        thresmeth = filters.threshold_otsu
    elif threshold == "mean":
        thresmeth = np.mean
    if isinstance(slice_or_arr, Slice):
        edges = filters.scharr(slice_or_arr.image.array.astype(float))
        center = slice_or_arr.image.center
    elif isinstance(slice_or_arr, np.ndarray):
        edges = filters.scharr(slice_or_arr.astype(float))
        center = (int(edges.shape[1] / 2), int(edges.shape[0] / 2))
    edges = filters.gaussian(edges, sigma=1)
    if isinstance(slice_or_arr, Slice):
        radius = 110 / slice_or_arr.mm_per_pixel
        rr, cc = draw.disk(
            center=(center.y, center.x), radius=radius, shape=edges.shape
        )
        thres = thresmeth(edges[rr, cc]) * 0.8
    else:
        thres = thresmeth(edges)
    bw = edges > thres
    if clear_borders:
        bw = segmentation.clear_border(bw, buffer_size=min(int(max(bw.shape) / 100), 3))
    if fill_holes:
        bw = ndimage.binary_fill_holes(bw)
    labeled_arr, num_roi = measure.label(bw, return_num=True)
    regionprops = measure.regionprops(labeled_arr, edges)
    return labeled_arr, regionprops, num_roi


def combine_surrounding_slices(
    dicomstack: DicomImageStack,
    nominal_slice_num: int,
    slices_plusminus: int = 1,
    mode: str = "mean",
) -> np.array:
    """Return an array that is the combination of a given slice and a number of slices surrounding it.

    Parameters
    ----------
    dicomstack : `~pylinac.core.image.DicomImageStack`
        The CBCT DICOM stack.
    nominal_slice_num : int
        The slice of interest (along 3rd dim).
    slices_plusminus: int
        How many slices plus and minus to combine (also along 3rd dim).
    mode : {'mean', 'median', 'max}
        Specifies the method of combination.

    Returns
    -------
    combined_array : numpy.array
        The combined array of the DICOM stack slices.
    """
    slices = range(
        nominal_slice_num - slices_plusminus, nominal_slice_num + slices_plusminus + 1
    )
    arrays = tuple(dicomstack[s].array for s in slices)
    array_stack = np.dstack(arrays)
    if mode == "mean":
        combined_array = np.mean(array_stack, 2)
    elif mode == "median":
        combined_array = np.median(array_stack, 2)
    else:
        combined_array = np.max(array_stack, 2)
    return combined_array


def rois_to_results(dict_mapping: dict[str, DiskROI]) -> dict[str, ROIResult]:
    """Converts a dict of HUDiskROIs to a dict of ROIResults. This is for dumping to simple data formats for results_data and RadMachine"""
    flat_dict = {}
    for name, roi in dict_mapping.items():
        flat_dict[name] = ROIResult(
            name=name,
            value=roi.pixel_value,
            stdev=roi.std,
            difference=getattr(roi, "value_diff", None),
            nominal_value=getattr(roi, "nominal_val", None),
            passed=getattr(roi, "passed", None),
        )
    return flat_dict
