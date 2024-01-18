"""The Winston-Lutz module loads and processes EPID images that have acquired Winston-Lutz type images.

Features:

* **Couch shift instructions** - After running a WL test, get immediate feedback on how to shift the couch.
  Couch values can also be passed in and the new couch values will be presented so you don't have to do that pesky conversion.
  "Do I subtract that number or add it?"
* **Automatic field & BB positioning** - When an image or directory is loaded, the field CAX and the BB
  are automatically found, along with the vector and scalar distance between them.
* **Isocenter size determination** - Using backprojections of the EPID images, the 3D gantry isocenter size
  and position can be determined *independent of the BB position*. Additionally, the 2D planar isocenter size
  of the collimator and couch can also be determined.
* **Image plotting** - WL images can be plotted separately or together, each of which shows the field CAX, BB and
  scalar distance from BB to CAX.
* **Axis deviation plots** - Plot the variation of the gantry, collimator, couch, and EPID in each plane
  as well as RMS variation.
* **File name interpretation** - Rename DICOM filenames to include axis information for linacs that don't include
  such information in the DICOM tags. E.g. "myWL_gantry45_coll0_couch315.dcm".
"""
from __future__ import annotations

import copy
import dataclasses
import enum
import io
import math
import os.path as osp
import statistics
import tempfile
import webbrowser
from dataclasses import dataclass
from functools import cached_property
from itertools import zip_longest
from pathlib import Path
from textwrap import wrap
from typing import BinaryIO, Iterable, Sequence, TypedDict

import argue
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import art3d
from scipy import linalg, ndimage, optimize
from scipy.ndimage import zoom
from scipy.spatial.transform import Rotation
from skimage import measure
from skimage.measure._regionprops import RegionProperties
from tabulate import tabulate

from .core import image, pdf
from .core.array_utils import array_to_dicom
from .core.decorators import lru_cache
from .core.geometry import Line, Point, Vector, cos, sin
from .core.image import DicomImageStack, LinacDicomImage, is_image, tiff_to_dicom
from .core.io import TemporaryZipDirectory, get_url, retrieve_demo_file
from .core.mask import bounding_box
from .core.scale import MachineScale, convert
from .core.utilities import ResultBase, convert_to_enum, is_close
from .metrics.features import (
    is_right_circumference,
    is_right_size_bb,
    is_round,
    is_solid,
    is_symmetric,
)
from .metrics.image import SizedDiskLocator

BB_ERROR_MESSAGE = (
    "Unable to locate the BB. Make sure the field edges do not obscure the BB, that there are no artifacts in the images, that the 'bb_size' parameter is close to reality, "
    "and that the BB is near the center (within 2cm). If this is a large-field image or kV image try setting 'low_density_bb' to True."
)


class BBArrangement:
    """Presets for multi-target phantoms."""

    # locations: https://www.postersessiononline.eu/173580348_eu/congresos/ESTRO2020/aula/-PO_1320_ESTRO2020.pdf
    SNC_MULTIMET = (
        {
            "name": "Iso",
            "offset_left_mm": 0,
            "offset_up_mm": 0,
            "offset_in_mm": 0,
            "bb_size_mm": 5,
            "rad_size_mm": 20,
        },
        {
            "name": "1",
            "offset_left_mm": 0,
            "offset_up_mm": 0,
            "offset_in_mm": 30,
            "bb_size_mm": 5,
            "rad_size_mm": 20,
        },
        {
            "name": "2",
            "offset_left_mm": -30,
            "offset_up_mm": 0,
            "offset_in_mm": 15,
            "bb_size_mm": 5,
            "rad_size_mm": 20,
        },
        {
            "name": "3",
            "offset_left_mm": 0,
            "offset_up_mm": 0,
            "offset_in_mm": -30,
            "bb_size_mm": 5,
            "rad_size_mm": 20,
        },
        {
            "name": "4",
            "offset_left_mm": 30,
            "offset_up_mm": 0,
            "offset_in_mm": -50,
            "bb_size_mm": 5,
            "rad_size_mm": 20,
        },
        {
            "name": "5",
            "offset_left_mm": 0,
            "offset_up_mm": 0,
            "offset_in_mm": -70,
            "bb_size_mm": 5,
            "rad_size_mm": 20,
        },
    )
    DEMO = (
        {
            "name": "Iso",
            "offset_left_mm": 0,
            "offset_up_mm": 0,
            "offset_in_mm": 0,
            "bb_size_mm": 5,
            "rad_size_mm": 20,
        },
        {
            "name": "Left,Down,In",
            "offset_left_mm": 20,
            "offset_up_mm": -20,
            "offset_in_mm": 60,
            "bb_size_mm": 5,
            "rad_size_mm": 20,
        },
    )

    @staticmethod
    def to_human(arrangement: dict) -> str:
        """Convert one BB location to a human-readable str"""
        a = arrangement
        lr = "Left" if a["offset_left_mm"] >= 0 else "Right"
        ud = "Up" if a["offset_up_mm"] >= 0 else "Down"
        io = "In" if a["offset_in_mm"] >= 0 else "Out"
        return f"'{a['name']}': {lr} {abs(a['offset_left_mm'])}mm, {ud} {abs(a['offset_up_mm'])}mm, {io} {abs(a['offset_in_mm'])}mm"


class NominalBB(TypedDict):
    """Input for BB location"""

    offset_left_mm: float
    offset_up_mm: float
    offset_in_mm: float
    bb_diameter_mm: float


class BB:
    """A representation of a BB in 3D space"""

    def __repr__(self):
        return self.nominal_position

    def __init__(self, nominal_bb: NominalBB, ray_lines: list[Line]):
        self.nominal_bb = nominal_bb
        self.ray_lines = ray_lines

    @cached_property
    def measured_position(self) -> Point:
        """The 3D measured position of the BB based on the ray-tracing lines in MM"""
        initial_guess = self.nominal_position.as_array()
        bounds = [(-200, 200), (-200, 200), (-200, 200)]
        result = optimize.minimize(
            max_distance_to_lines, initial_guess, args=self.ray_lines, bounds=bounds
        )
        return Point(result.x)

    @cached_property
    def nominal_position(self) -> Point:
        """The nominal location of the BB in MM"""
        return Point(
            x=-self.nominal_bb["offset_left_mm"],
            y=-self.nominal_bb["offset_in_mm"],
            z=self.nominal_bb["offset_up_mm"],
        )

    @cached_property
    def delta_vector(self) -> Vector:
        """The shift from measured BB location to nominal as a vector in MM"""
        return self.measured_position - self.nominal_position

    @cached_property
    def delta_distance(self):
        """The scalar distance between the measured BB location and nominal in MM"""
        return self.measured_position.distance_to(self.nominal_position)

    def plot_nominal(self, axes: plt.Axes, color: str, **kwargs):
        """Plot the BB nominal position"""
        x, y, z = create_sphere_surface(
            radius=self.nominal_bb["bb_diameter_mm"] / 2, center=self.nominal_position
        )
        axes.plot_surface(x, y, z, color=color, **kwargs)

    def plot_measured(self, axes: plt.Axes, color: str, **kwargs):
        """Plot the BB measured position"""
        x, y, z = create_sphere_surface(
            radius=self.nominal_bb["bb_diameter_mm"] / 2, center=self.measured_position
        )
        axes.plot_surface(x, y, z, color=color, **kwargs)


def create_sphere_surface(
    radius: float, center: Point
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create a sphere surface for plotting"""
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = radius * np.outer(np.cos(u), np.sin(v)) + center.x
    y = radius * np.outer(np.sin(u), np.sin(v)) + center.y
    z = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + center.z
    return x, y, z


class Axis(enum.Enum):
    GANTRY = "Gantry"  #:
    COLLIMATOR = "Collimator"  #:
    COUCH = "Couch"  #:
    GB_COMBO = "GB Combo"  #:
    GBP_COMBO = "GBP Combo"  #:
    EPID = "Epid"  #:
    REFERENCE = "Reference"  #:


@dataclass
class WinstonLutz2DResult(ResultBase):
    variable_axis: str  #:
    cax2epid_vector: Vector  #:
    cax2epid_distance: float  #:
    cax2bb_distance: float  #:
    cax2bb_vector: Vector  #:
    bb_location: Point  #:
    field_cax: Point  #:


@dataclass
class WinstonLutzResult(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    num_gantry_images: int  #:
    num_gantry_coll_images: int  #:
    num_coll_images: int  #:
    num_couch_images: int  #:
    num_total_images: int  #:
    max_2d_cax_to_bb_mm: float  #:
    median_2d_cax_to_bb_mm: float  #:
    mean_2d_cax_to_bb_mm: float  #:
    max_2d_cax_to_epid_mm: float  #:
    median_2d_cax_to_epid_mm: float  #:
    mean_2d_cax_to_epid_mm: float  #:
    gantry_3d_iso_diameter_mm: float  #:
    max_gantry_rms_deviation_mm: float  #:
    max_epid_rms_deviation_mm: float  #:
    gantry_coll_3d_iso_diameter_mm: float  #:
    coll_2d_iso_diameter_mm: float  #:
    max_coll_rms_deviation_mm: float  #:
    couch_2d_iso_diameter_mm: float  #:
    max_couch_rms_deviation_mm: float  #:
    image_details: list[WinstonLutz2DResult]  #:
    keyed_image_details: dict[str, WinstonLutz2DResult]  #:


@dataclass
class WinstonLutzMultiTargetMultiFieldResult(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    num_total_images: int  #:
    max_2d_field_to_bb_mm: float  #:
    median_2d_field_to_bb_mm: float  #:
    mean_2d_field_to_bb_mm: float  #:
    bb_arrangement: Iterable[dict]  #:
    bb_maxes: dict  #:


def plot_image(img: WinstonLutz2D | None, axis: plt.Axes) -> None:
    """Helper function to plot a WLImage to an axis."""
    if img is None:
        axis.set_frame_on(False)
        axis.axis("off")
    else:
        img.plot(ax=axis, show=False)


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


def is_modest_size(region: RegionProperties, *args, **kwargs) -> bool:
    """Decide whether the ROI is roughly the size of a BB; not noise and not an artifact. Used to find the BB."""
    bb_area = region.area_filled / (kwargs["dpmm"] ** 2)
    bb_size = kwargs["bb_size"]
    larger_bb_area = np.pi * ((bb_size + 2) / 2) ** 2
    smaller_bb_area = max(
        (np.pi * ((bb_size - 2) / 2) ** 2, 2)
    )  # set a min of 2 to avoid a lower bound of 0 when radius=2. This is much more likely to find noise in a block.
    return smaller_bb_area < bb_area < larger_bb_area


def is_square(region: RegionProperties, *args, **kwargs) -> bool:
    """Decide if the ROI is square in nature by testing the filled area vs bounding box. Used to find the BB."""
    actual_fill_ratio = region.filled_area / region.bbox_area
    return actual_fill_ratio > 0.8


def is_right_square_size(region: RegionProperties, *args, **kwargs) -> bool:
    """Decide if the ROI is square in nature by testing the filled area vs bounding box. Used to find the BB."""
    field_area = region.area_filled / (kwargs["dpmm"] ** 2)
    rad_size = max((kwargs["rad_size"], 5))
    larger_bb_area = (rad_size + 5) ** 2
    smaller_bb_area = (rad_size - 5) ** 2
    return smaller_bb_area < field_area < larger_bb_area


class WinstonLutz2D(image.LinacDicomImage):
    """Holds individual Winston-Lutz EPID images, image properties, and automatically finds the field CAX and BB."""

    bb: Point
    field_cax: Point
    _rad_field_bounding_box: list
    is_from_tiff: bool = False
    detection_conditions: list[callable] = [
        is_right_size_bb,
        is_round,
        is_right_circumference,
        is_symmetric,
        is_solid,
    ]

    def __init__(
        self,
        file: str | BinaryIO | Path,
        use_filenames: bool = False,
        **kwargs,
    ):
        """
        Parameters
        ----------
        file : str
            Path to the image file.
        use_filenames: bool
            Whether to try to use the file name to determine axis values.
            Useful for Elekta machines that do not include that info in the DICOM data.
        """
        # override detection conditions if passed
        if conditions := kwargs.pop("detection_conditions", False):
            self.detection_conditions = conditions
        super().__init__(file, use_filenames=use_filenames, **kwargs)
        self._is_analyzed = False
        self.flipud()

    def analyze(
        self,
        bb_size_mm: float = 5,
        low_density_bb: bool = False,
        open_field: bool = False,
    ) -> None:
        """Analyze the image. See WinstonLutz.analyze for parameter details."""
        self.check_inversion_by_histogram(percentiles=(0.01, 50, 99.99))
        self._clean_edges()
        self.ground()
        self.normalize()
        self.field_cax, self._rad_field_bounding_box = self._find_field_centroid(
            open_field
        )
        self.bb = self._find_bb(bb_size_mm, low_density_bb)
        self._is_analyzed = True

    def __repr__(self):
        return f"WLImage(gantry={self.gantry_angle:.1f}, coll={self.collimator_angle:.1f}, couch={self.couch_angle:.1f})"

    def to_axes(self) -> str:
        """Give just the axes values as a human-readable string"""
        return f"Gantry={self.gantry_angle:.1f}, Coll={self.collimator_angle:.1f}, Couch={self.couch_angle:.1f}"

    def _clean_edges(self, window_size: int = 2) -> None:
        """Clean the edges of the image to be near the background level."""

        def has_noise(self, window_size):
            """Helper method to determine if there is spurious signal at any of the image edges.

            Determines if the min or max of an edge is within 10% of the baseline value and trims if not.
            """
            near_min, near_max = np.percentile(self.array, [5, 99.5])
            img_range = near_max - near_min
            top = self[:window_size, :]
            left = self[:, :window_size]
            bottom = self[-window_size:, :]
            right = self[:, -window_size:]
            edge_array = np.concatenate(
                (top.flatten(), left.flatten(), bottom.flatten(), right.flatten())
            )
            edge_too_low = edge_array.min() < (near_min - img_range / 10)
            edge_too_high = edge_array.max() > (near_max + img_range / 10)
            return edge_too_low or edge_too_high

        safety_stop = np.min(self.shape) / 10
        while has_noise(self, window_size) and safety_stop > 0:
            self.crop(window_size)
            safety_stop -= 1

    def _find_field_centroid(self, is_open_field: bool) -> tuple[Point, list]:
        """Find the centroid of the radiation field.

        Parameters
        ----------
        is_open_field
            If True, simply uses the image/EPID center as the field center.
            If False, finds the radiation field based on a 50% height threshold.

        Returns
        -------
        p
            The CAX point location.
        edges
            The bounding box of the field, plus a small margin.
        """
        if is_open_field:
            p = self.center
            edges = [0, self.shape[0], 0, self.shape[1]]
        else:
            min, max = np.percentile(self.array, [5, 99.9])
            threshold_img = self.as_binary((max - min) / 2 + min)
            filled_img = ndimage.binary_fill_holes(threshold_img)
            # clean single-pixel noise from outside field
            cleaned_img = ndimage.binary_erosion(threshold_img)
            [*edges] = bounding_box(cleaned_img)
            edges[0] -= 10
            edges[1] += 10
            edges[2] -= 10
            edges[3] += 10
            coords = ndimage.center_of_mass(filled_img)
            p = Point(x=coords[-1], y=coords[0])
        return p, edges

    def _calculate_bb_tolerance(self, bb_diameter: float) -> int:
        """Calculate the BB tolerance based on the BB diameter.
        Min will be 2 for 1.5mm and under. Will be 4 for diameters at or above 30mm."""
        y = (2, 4)
        x = (1.5, 30)
        return np.interp(bb_diameter, x, y)

    def _find_bb(self, bb_diameter: float, low_density: bool) -> Point:
        """Find the BB within the radiation field. Iteratively searches for a circle-like object
        by lowering a low-pass threshold value until found.

        Returns
        -------
        Point
            The weighted-pixel value location of the BB.
        """
        bb_tolerance_mm = self._calculate_bb_tolerance(bb_diameter)
        center = self.compute(
            metrics=SizedDiskLocator.from_center_physical(
                expected_position_mm=(0, 0),
                search_window_mm=(40 + bb_diameter, 40 + bb_diameter),
                radius_mm=bb_diameter / 2,
                radius_tolerance_mm=bb_tolerance_mm,
                invert=not low_density,
                detection_conditions=self.detection_conditions,
            )
        )
        return center

    @property
    def epid(self) -> Point:
        """Center of the EPID panel"""
        return self.center

    @property
    def cax_line_projection(self) -> Line:
        """The projection of the field CAX through space around the area of the BB.
        Used for determining gantry isocenter size.

        Returns
        -------
        Line
            The virtual line in space made by the beam CAX.
        """
        p1 = Point()
        p2 = Point()
        # point 1 - ray origin
        p1.x = self.cax2bb_vector.x * cos(self.gantry_angle) + 20 * sin(
            self.gantry_angle
        )
        p1.z = self.cax2bb_vector.x * -sin(self.gantry_angle) + 20 * cos(
            self.gantry_angle
        )
        p1.y = self.cax2bb_vector.y
        # point 2 - ray destination
        p2.x = self.cax2bb_vector.x * cos(self.gantry_angle) - 20 * sin(
            self.gantry_angle
        )
        p2.z = self.cax2bb_vector.x * -sin(self.gantry_angle) - 20 * cos(
            self.gantry_angle
        )
        p2.y = self.cax2bb_vector.y
        line = Line(p1, p2)
        return line

    @property
    def cax2bb_vector(self) -> Vector:
        """The vector in mm from the CAX to the BB."""
        dist = (self.bb - self.field_cax) / self.dpmm
        return Vector(dist.x, dist.y, dist.z)

    @property
    def cax2bb_distance(self) -> float:
        """The scalar distance in mm from the CAX to the BB."""
        dist = self.field_cax.distance_to(self.bb)
        return dist / self.dpmm

    @property
    def cax2epid_vector(self) -> Vector:
        """The vector in mm from the CAX to the EPID center pixel"""
        dist = (self.epid - self.field_cax) / self.dpmm
        return Vector(dist.x, dist.y, dist.z)

    @property
    def cax2epid_distance(self) -> float:
        """The scalar distance in mm from the CAX to the EPID center pixel"""
        return self.field_cax.distance_to(self.epid) / self.dpmm

    def plot(
        self, ax: plt.Axes | None = None, show: bool = True, clear_fig: bool = False
    ) -> plt.Axes:
        """Plot the image, zoomed-in on the radiation field, along with the detected
        BB location and field CAX location.

        Parameters
        ----------
        ax : None, matplotlib Axes instance
            The axis to plot to. If None, will create a new figure.
        show : bool
            Whether to actually show the image.
        clear_fig : bool
            Whether to clear the figure first before drawing.
        """
        ax = super().plot(ax=ax, show=False, clear_fig=clear_fig)
        ax.plot(self.field_cax.x, self.field_cax.y, "gs", ms=8)
        ax.axvline(x=self.epid.x, color="b")
        ax.axhline(y=self.epid.y, color="b")
        ax.set_ylim([self._rad_field_bounding_box[0], self._rad_field_bounding_box[1]])
        ax.set_xlim([self._rad_field_bounding_box[2], self._rad_field_bounding_box[3]])
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_title("\n".join(wrap(str(self.path), 30)), fontsize=10)
        ax.set_xlabel(
            f"G={self.gantry_angle:.0f}, B={self.collimator_angle:.0f}, P={self.couch_angle:.0f}"
        )
        ax.set_ylabel(
            f"CAX to BB: {self.cax2bb_distance:3.2f}mm\nX: {self.cax2bb_vector.x:3.2f}mm; Y: {self.cax2bb_vector.y:3.2f}mm"
        )
        if show:
            plt.show()
        return ax

    def save_plot(self, filename: str, **kwargs):
        """Save the image plot to file."""
        self.plot(show=False)
        plt.tight_layout()
        plt.savefig(filename, **kwargs)

    @property
    def variable_axis(self) -> Axis:
        """The axis that is varying.

        There are five types of images:

        * Reference : All axes are at 0.
        * Gantry: All axes but gantry at 0.
        * Collimator : All axes but collimator at 0.
        * Couch : All axes but couch at 0.
        * Combo : More than one axis is not at 0.
        """
        G0 = is_close(self.gantry_angle, [0, 360])
        B0 = is_close(self.collimator_angle, [0, 360])
        P0 = is_close(self.couch_angle, [0, 360])
        if G0 and B0 and not P0:
            return Axis.COUCH
        elif G0 and P0 and not B0:
            return Axis.COLLIMATOR
        elif P0 and B0 and not G0:
            return Axis.GANTRY
        elif P0 and B0 and G0:
            return Axis.REFERENCE
        elif P0:
            return Axis.GB_COMBO
        else:
            return Axis.GBP_COMBO

    def results_data(self, as_dict: bool = False) -> WinstonLutz2DResult | dict:
        """Present the results data and metadata as a dataclass or dict.
        The default return type is a dataclass."""
        if not self._is_analyzed:
            raise ValueError("The image is not analyzed. Use .analyze() first.")

        data = WinstonLutz2DResult(
            variable_axis=self.variable_axis.value,
            cax2bb_vector=self.cax2bb_vector,
            cax2epid_vector=self.cax2epid_vector,
            cax2bb_distance=self.cax2bb_distance,
            cax2epid_distance=self.cax2epid_distance,
            bb_location=self.bb,
            field_cax=self.field_cax,
        )
        if as_dict:
            return dataclasses.asdict(data)
        return data


class WinstonLutz:
    """Class for performing a Winston-Lutz test of the radiation isocenter."""

    images: list[WinstonLutz2D]  #:
    machine_scale: MachineScale  #:
    image_type = WinstonLutz2D
    is_from_cbct: bool = False
    _bb_diameter = float
    detection_conditions: list[callable] = [
        is_right_size_bb,
        is_round,
        is_right_circumference,
        is_symmetric,
        is_solid,
    ]

    def __init__(
        self,
        directory: str | list[str] | Path,
        use_filenames: bool = False,
        axis_mapping: dict[str, tuple[int, int, int]] | None = None,
        axes_precision: int | None = None,
        dpi: float | None = None,
        sid: float | None = None,
    ):
        """
        Parameters
        ----------
        directory : str, list[str]
            Path to the directory of the Winston-Lutz EPID images or a list of the image paths
        use_filenames: bool
            Whether to try to use the file name to determine axis values.
            Useful for Elekta machines that do not include that info in the DICOM data.
            This is mutually exclusive to axis_mapping. If True, axis_mapping is ignored.
        axis_mapping: dict
            An optional way of instantiating by passing each file along with the axis values.
            Structure should be <filename>: (<gantry>, <coll>, <couch>).
        axes_precision: int | None
            How many significant digits to represent the axes values. If None, no precision is set and the input/DICOM values are used raw.
            If set to an integer, rounds the axes values (gantry, coll, couch) to that many values. E.g. gantry=0.1234 => 0.1 with precision=1.
            This is mostly useful for plotting/rounding (359.9=>0) and if using the ``keyed_image_details`` with ``results_data``.
        dpi
            The dots-per-inch setting. Only needed if using TIFF images and the images do not contain the resolution tag.
            An error will raise if dpi is not passed and the TIFF resolution cannot be determined.
        sid
            The Source-to-Image distance in mm. Only needed when using TIFF images.
        """
        self.images = []
        if axis_mapping and not use_filenames:
            for filename, (gantry, coll, couch) in axis_mapping.items():
                self.images.append(
                    self._load_image(
                        Path(directory) / filename,
                        sid=sid,
                        dpi=dpi,
                        gantry=gantry,
                        coll=coll,
                        couch=couch,
                        axes_precision=axes_precision,
                    )
                )
        elif isinstance(directory, list):
            for file in directory:
                if is_image(file):
                    self.images.append(
                        self._load_image(
                            file,
                            dpi=dpi,
                            sid=sid,
                            use_filenames=use_filenames,
                            axes_precision=axes_precision,
                        )
                    )
        elif not osp.isdir(directory):
            raise ValueError(
                "Invalid directory passed. Check the correct method and file was used."
            )
        else:
            image_files = image.retrieve_image_files(directory)
            for file in image_files:
                self.images.append(
                    self._load_image(
                        file,
                        dpi=dpi,
                        sid=sid,
                        use_filenames=use_filenames,
                        axes_precision=axes_precision,
                    )
                )
        if len(self.images) < 2:
            raise ValueError(
                "<2 valid WL images were found in the folder/file or passed. Ensure you chose the correct folder/file for analysis."
            )
        self.images.sort(
            key=lambda i: (i.gantry_angle, i.collimator_angle, i.couch_angle)
        )
        self._is_analyzed = False

    def _load_image(
        self,
        file: str | Path,
        sid: float | None,
        dpi: float | None,
        **kwargs,
    ) -> WinstonLutz2D:
        """A helper method to load either DICOM or TIFF files appropriately."""
        try:
            return self.image_type(
                file, detection_conditions=self.detection_conditions, **kwargs
            )
        except AttributeError:
            if kwargs.get("gantry") is None:
                raise ValueError(
                    "TIFF images detected. Must pass `axis_mapping` parameter."
                )
            if sid is None:
                raise ValueError("TIFF images detected. Must pass `sid` parameter")
            with io.BytesIO() as stream:
                ds = tiff_to_dicom(
                    tiff_file=file,
                    sid=sid,
                    dpi=dpi,
                    gantry=kwargs.pop("gantry"),
                    coll=kwargs.pop("coll"),
                    couch=kwargs.pop("couch"),
                )
                ds.save_as(stream, write_like_original=False)
                img = self.image_type(
                    stream, detection_conditions=self.detection_conditions, **kwargs
                )
                if not img.dpmm:
                    raise ValueError(
                        "TIFF images were detected but the dpi tag was not available. Pass the `dpi` parameter manually."
                    )
                img.filter(size=0.01, kind="median")
                return img

    @classmethod
    def from_demo_images(cls, **kwargs):
        """Instantiate using the demo images.

        Parameters
        ----------
        kwargs
            See parameters of the __init__ method for details.
        """
        demo_file = retrieve_demo_file(name="winston_lutz.zip")
        return cls.from_zip(demo_file, **kwargs)

    @classmethod
    def from_zip(cls, zfile: str | BinaryIO, **kwargs):
        """Instantiate from a zip file rather than a directory.

        Parameters
        ----------
        zfile
            Path to the archive file.
        kwargs
            See parameters of the __init__ method for details.
        """
        with TemporaryZipDirectory(zfile) as tmpz:
            obj = cls(tmpz, **kwargs)
        return obj

    @classmethod
    def from_url(cls, url: str, **kwargs):
        """Instantiate from a URL.

        Parameters
        ----------
        url : str
            URL that points to a zip archive of the DICOM images.
        kwargs
            See parameters of the __init__ method for details.
        """
        zfile = get_url(url)
        return cls.from_zip(zfile, **kwargs)

    @classmethod
    def from_cbct_zip(cls, file: Path | str, raw_pixels: bool = False, **kwargs):
        """Instantiate from a zip file containing CBCT images.

        Parameters
        ----------
        file
            Path to the archive file.
        raw_pixels
            If True, uses the raw pixel values of the DICOM files. If False, uses the rescaled Hounsfield units.
            Generally, this should be true.
        kwargs
            See parameters of the __init__ method for details.
        """
        with TemporaryZipDirectory(file) as tmpz:
            obj = cls.from_cbct(tmpz, raw_pixels=raw_pixels, **kwargs)
        return obj

    @classmethod
    def from_cbct(cls, directory: Path | str, raw_pixels: bool = False, **kwargs):
        """Create a 4-angle WL test from a CBCT dataset.

        The dataset is loaded and the array is "viewed" from top, bottom, left, and right to create the 4 angles.
        The dataset has to be rescaled so that the z-axis spacing is equal to the x/y axis. This is because the
        typical slice thickness is much larger than the in-plane resolution.

        Parameters
        ----------
        directory
            The directory containing the CBCT DICOM files.
        raw_pixels
            If True, uses the raw pixel values of the DICOM files. If False, uses the rescaled Hounsfield units.
            Generally, this should be true.
        kwargs
            See parameters of the __init__ method for details.
        """
        dicom_stack = DicomImageStack(
            folder=directory, min_number=10, raw_pixels=raw_pixels
        )
        np_stack = np.stack(dicom_stack.images, axis=-1)
        zoom_ratio = (
            1,
            dicom_stack.metadata.SliceThickness / dicom_stack.metadata.PixelSpacing[0],
        )
        left_arr = np.rot90(
            zoom(
                np_stack.max(axis=0),
                zoom=zoom_ratio,
                grid_mode=True,
                mode="nearest",
                order=1,
            ),
            k=1,
        )
        top_arr = np.rot90(
            zoom(
                np_stack.max(axis=1),
                zoom=zoom_ratio,
                grid_mode=True,
                mode="nearest",
                order=1,
            ),
            k=1,
        )
        right_arr = np.fliplr(left_arr)
        bottom_arr = np.fliplr(top_arr)
        dicom_dir = Path(tempfile.mkdtemp())
        dpi = 25.4 / dicom_stack.metadata.PixelSpacing[0]
        for array, gantry in zip(
            (left_arr, top_arr, right_arr, bottom_arr), (270, 0, 90, 180)
        ):
            ds = array_to_dicom(
                array=np.ascontiguousarray(array),  # pydicom complains due to np.rot90
                sid=1000,
                gantry=gantry,
                coll=0,
                couch=0,
                dpi=dpi,
            )
            ds.save_as(dicom_dir / f"G={gantry}", write_like_original=False)
        # now we load these as normal images into the WL algorithm
        instance = cls(dicom_dir, **kwargs)
        instance.is_from_cbct = True
        return instance

    @staticmethod
    def run_demo() -> None:
        """Run the Winston-Lutz demo, which loads the demo files, prints results, and plots a summary image."""
        wl = WinstonLutz.from_demo_images()
        wl.analyze(machine_scale=MachineScale.VARIAN_IEC)
        print(wl.results())
        wl.plot_summary()

    def analyze(
        self,
        bb_size_mm: float = 5,
        machine_scale: MachineScale = MachineScale.IEC61217,
        low_density_bb: bool = False,
        open_field: bool = False,
    ) -> None:
        """Analyze the WL images.

        Parameters
        ----------
        bb_size_mm
            The expected diameter of the BB in mm. The actual size of the BB can be +/-2mm from the passed value.
        machine_scale
            The scale of the machine. Shift vectors depend on this value.
        low_density_bb
            Set this flag to True if the BB is lower density than the material surrounding it.
        open_field
            If True, sets the field center to the EPID center under the assumption the field is not the focus of interest or is too wide to be calculated.
            This is often helpful for kV WL analysis where the blades are wide open and even then the blade edge is of
            less interest than simply the imaging iso vs the BB.
        """
        self.machine_scale = machine_scale
        if self.is_from_cbct:
            low_density_bb = True
            open_field = True
        for img in self.images:
            img.analyze(bb_size_mm, low_density_bb, open_field)
        self._is_analyzed = True
        self._bb_diameter = bb_size_mm

    @lru_cache()
    def _minimize_axis(self, axes: Axis | tuple[Axis, ...] = (Axis.GANTRY,)):
        """Return the minimization result of the given axis."""
        if isinstance(axes, Axis):
            axes = (axes,)

        things = [
            image.cax_line_projection
            for image in self.images
            if image.variable_axis in (axes + (Axis.REFERENCE,))
        ]
        if len(things) <= 1:
            raise ValueError(
                "Not enough images of the given type to identify the axis isocenter"
            )
        initial_guess = np.array([0, 0, 0])
        bounds = [(-20, 20), (-20, 20), (-20, 20)]
        result = optimize.minimize(
            max_distance_to_lines, initial_guess, args=things, bounds=bounds
        )
        return result

    @property
    def gantry_iso_size(self) -> float:
        """The diameter of the 3D gantry isocenter size in mm. Only images where the collimator
        and couch were at 0 are used to determine this value."""
        num_gantry_like_images = self._get_images((Axis.GANTRY, Axis.REFERENCE))[0]
        if num_gantry_like_images > 1:
            return self._minimize_axis(Axis.GANTRY).fun * 2
        else:
            return 0

    @property
    def gantry_coll_iso_size(self) -> float:
        """The diameter of the 3D gantry isocenter size in mm *including collimator and gantry/coll combo images*.
        Images where the couch!=0 are excluded."""
        num_gantry_like_images = self._get_images(
            (Axis.GANTRY, Axis.COLLIMATOR, Axis.GB_COMBO, Axis.REFERENCE)
        )[0]
        if num_gantry_like_images > 1:
            return (
                self._minimize_axis((Axis.GANTRY, Axis.COLLIMATOR, Axis.GB_COMBO)).fun
                * 2
            )
        else:
            return 0

    @staticmethod
    def _find_max_distance_between_points(images) -> float:
        """Find the maximum distance between a set of points. Used for 2D images like collimator and couch."""
        points = [
            Point(image.cax2bb_vector.x, image.cax2bb_vector.y) for image in images
        ]
        dists = []
        for point1 in points:
            for point2 in points:
                p = point1.distance_to(point2)
                dists.append(p)
        return max(dists)

    @property
    def collimator_iso_size(self) -> float:
        """The 2D collimator isocenter size (diameter) in mm. The iso size is in the plane
        normal to the gantry."""
        num_collimator_like_images, images = self._get_images(
            (Axis.COLLIMATOR, Axis.REFERENCE)
        )
        if num_collimator_like_images > 1:
            return self._find_max_distance_between_points(images)
        else:
            return 0

    @property
    def couch_iso_size(self) -> float:
        """The diameter of the 2D couch isocenter size in mm. Only images where
        the gantry and collimator were at zero are used to determine this value."""
        num_couch_like_images, images = self._get_images((Axis.COUCH, Axis.REFERENCE))
        if num_couch_like_images > 1:
            return self._find_max_distance_between_points(images)
        else:
            return 0

    @property
    def bb_shift_vector(self) -> Vector:
        """The shift necessary to place the BB at the radiation isocenter.
        The values are in the coordinates defined in the documentation.

        The shift is based on the paper by Low et al. See online documentation for more.
        """
        A = np.empty([2 * len(self.images), 3])
        epsilon = np.empty([2 * len(self.images), 1])
        for idx, img in enumerate(self.images):
            # convert from input scale to Varian scale
            # Low's paper assumes Varian scale input and this is easier than changing the actual signs in the equation which have a non-intuitive relationship
            gantry, _, couch = convert(
                input_scale=self.machine_scale,
                output_scale=MachineScale.VARIAN_STANDARD,
                gantry=img.gantry_angle,
                collimator=img.collimator_angle,
                rotation=img.couch_angle,
            )
            A[2 * idx : 2 * idx + 2, :] = np.array(
                [
                    [-cos(couch), -sin(couch), 0],
                    [-cos(gantry) * sin(couch), cos(gantry) * cos(couch), -sin(gantry)],
                ]
            )  # equation 6 (minus delta)
            epsilon[2 * idx : 2 * idx + 2] = np.array(
                [[img.cax2bb_vector.y], [img.cax2bb_vector.x]]
            )  # equation 7

        B = linalg.pinv(A)
        delta = B.dot(epsilon)  # equation 9
        # we use the negative for all values because it's from the iso POV -> linac not the linac -> iso POV
        return Vector(x=-delta[1][0], y=-delta[0][0], z=-delta[2][0])

    def bb_shift_instructions(
        self,
        couch_vrt: float | None = None,
        couch_lng: float | None = None,
        couch_lat: float | None = None,
    ) -> str:
        """Returns a string describing how to shift the BB to the radiation isocenter looking from the foot of the couch.
        Optionally, the current couch values can be passed in to get the new couch values. If passing the current
        couch position all values must be passed.

        Parameters
        ----------
        couch_vrt : float
            The current couch vertical position in cm.
        couch_lng : float
            The current couch longitudinal position in cm.
        couch_lat : float
            The current couch lateral position in cm.
        """
        sv = self.bb_shift_vector
        x_dir = "LEFT" if sv.x < 0 else "RIGHT"
        y_dir = "IN" if sv.y > 0 else "OUT"
        z_dir = "UP" if sv.z > 0 else "DOWN"
        move = f"{x_dir} {abs(sv.x):2.2f}mm; {y_dir} {abs(sv.y):2.2f}mm; {z_dir} {abs(sv.z):2.2f}mm"
        if all(val is not None for val in [couch_vrt, couch_lat, couch_lng]):
            new_lat = round(couch_lat + sv.x / 10, 2)
            new_vrt = round(couch_vrt + sv.z / 10, 2)
            new_lng = round(couch_lng + sv.y / 10, 2)
            move += f"\nNew couch coordinates (mm): VRT: {new_vrt:3.2f}; LNG: {new_lng:3.2f}; LAT: {new_lat:3.2f}"
        return move

    @argue.options(value=("all", "range"))
    def axis_rms_deviation(
        self, axis: Axis | tuple[Axis, ...] = Axis.GANTRY, value: str = "all"
    ) -> Iterable | float:
        """The RMS deviations of a given axis/axes.

        Parameters
        ----------
        axis : ('Gantry', 'Collimator', 'Couch', 'Epid', 'GB Combo',  'GBP Combo')
            The axis desired.
        value : {'all', 'range'}
            Whether to return all the RMS values from all images for that axis, or only return the maximum range of
            values, i.e. the 'sag'.
        """
        if isinstance(axis, Iterable):
            axis = [convert_to_enum(ax, Axis) for ax in axis]
        else:
            axis = convert_to_enum(axis, Axis)
        if axis != Axis.EPID:
            attr = "cax2bb_vector"
        else:
            attr = "cax2epid_vector"
            axis = (Axis.GANTRY, Axis.COLLIMATOR, Axis.REFERENCE)
        imgs = self._get_images(axis=axis)[1]
        if len(imgs) <= 1:
            return (0,)
        rms = [getattr(img, attr).as_scalar() for img in imgs]
        if value == "range":
            rms = max(rms) - min(rms)
        return rms

    @argue.options(metric=("max", "median", "mean"))
    def cax2bb_distance(self, metric: str = "max") -> float:
        """The distance in mm between the CAX and BB for all images according to the given metric.

        Parameters
        ----------
        metric : {'max', 'median', 'mean'}
            The metric of distance to use.
        """
        if metric == "max":
            return max(image.cax2bb_distance for image in self.images)
        elif metric == "median":
            return float(np.median([image.cax2bb_distance for image in self.images]))
        elif metric == "mean":
            return float(np.mean([image.cax2bb_distance for image in self.images]))

    @argue.options(metric=("max", "median", "mean"))
    def cax2epid_distance(self, metric: str = "max") -> float:
        """The distance in mm between the CAX and EPID center pixel for all images according to the given metric.

        Parameters
        ----------
        metric : {'max', 'median', 'mean'}
            The metric of distance to use.
        """
        if metric == "max":
            return max(image.cax2epid_distance for image in self.images)
        elif metric == "median":
            return float(np.median([image.cax2epid_distance for image in self.images]))
        elif metric == "mean":
            return float(np.mean([image.cax2epid_distance for image in self.images]))

    def _plot_deviation(
        self, axis: Axis, ax: plt.Axes | None = None, show: bool = True
    ) -> None:
        """Helper function: Plot the sag in Cartesian coordinates.

        Parameters
        ----------
        axis : {'gantry', 'epid', 'collimator', 'couch'}
            The axis to plot.
        ax : None, matplotlib.Axes
            The axis to plot to. If None, creates a new plot.
        show : bool
            Whether to show the image.
        """
        title = f"In-plane {axis.value} displacement"
        if axis == Axis.EPID:
            attr = "cax2epid_vector"
            axis = Axis.GANTRY
        else:
            attr = "cax2bb_vector"
        # get axis images, angles, and shifts
        imgs = [
            image
            for image in self.images
            if image.variable_axis in (axis, Axis.REFERENCE)
        ]
        angles = [getattr(image, f"{axis.value.lower()}_angle") for image in imgs]
        xz_sag = np.array([getattr(img, attr).x for img in imgs])
        y_sag = np.array([getattr(img, attr).y for img in imgs])
        rms = np.sqrt(xz_sag**2 + y_sag**2)

        # plot the axis deviation
        if ax is None:
            ax = plt.subplot(111)
        ax.plot(angles, y_sag, "bo", label="Y-axis", ls="-.")
        ax.plot(angles, xz_sag, "m^", label="X/Z-axis", ls="-.")
        ax.plot(angles, rms, "g+", label="RMS", ls="-")
        ax.set_title(title)
        ax.set_ylabel("mm")
        ax.set_xlabel(f"{axis.value} angle")
        ax.set_xticks(np.arange(0, 361, 45))
        ax.set_xlim(-15, 375)
        ax.grid(True)
        ax.legend(numpoints=1)
        if show:
            plt.show()

    def _get_images(
        self, axis: Axis | tuple[Axis, ...] = (Axis.GANTRY,)
    ) -> tuple[float, list]:
        if isinstance(axis, Axis):
            axis = (axis,)
        images = [image for image in self.images if image.variable_axis in axis]
        return len(images), images

    def plot_axis_images(
        self, axis: Axis = Axis.GANTRY, show: bool = True, ax: plt.Axes | None = None
    ) -> None:
        """Plot all CAX/BB/EPID positions for the images of a given axis.

        For example, axis='Couch' plots a reference image, and all the BB points of the other
        images where the couch was moving.

        Parameters
        ----------
        axis : {'Gantry', 'Collimator', 'Couch', 'GB Combo',  'GBP Combo'}
            The images/markers from which accelerator axis to plot.
        show : bool
            Whether to actually show the images.
        ax : None, matplotlib.Axes
            The axis to plot to. If None, creates a new plot.
        """
        axis = convert_to_enum(axis, Axis)
        images = [
            image
            for image in self.images
            if image.variable_axis in (axis, Axis.REFERENCE)
        ]
        ax = images[0].plot(
            show=False, ax=ax
        )  # plots the first marker; plot the rest of the markers below
        if axis != Axis.COUCH:
            # plot EPID
            epid_xs = [img.epid.x for img in images[1:]]
            epid_ys = [img.epid.y for img in images[1:]]
            ax.plot(epid_xs, epid_ys, "b+", ms=8)
            # get CAX positions
            xs = [img.field_cax.x for img in images[1:]]
            ys = [img.field_cax.y for img in images[1:]]
            marker = "gs"
        else:
            # get BB positions
            xs = [img.bb.x for img in images[1:]]
            ys = [img.bb.y for img in images[1:]]
            marker = "ro"
        ax.plot(xs, ys, marker, ms=8)
        # set labels
        ax.set_title(axis.value + " wobble")
        ax.set_xlabel(axis.value + " positions superimposed")
        ax.set_ylabel(
            axis.value
            + f" iso size: {getattr(self, axis.value.lower() + '_iso_size'):3.2f}mm"
        )
        if show:
            plt.show()

    def plot_location(
        self,
        show: bool = True,
        viewbox_mm: float | None = None,
        plot_bb: bool = True,
        plot_isocenter_sphere: bool = True,
        plot_couch_iso: bool = True,
        plot_coll_iso: bool = True,
        show_legend: bool = True,
    ):
        """Plot the isocenter and size as a sphere in 3D space relative to the BB. The
        iso is at the origin.

        Only images where the couch was at zero are considered.

        Parameters
        ----------
        show : bool
            Whether to plot the image.
        viewbox_mm : float
            The default size of the 3D space to plot in mm in each axis.
        plot_bb : bool
            Whether to plot the BB location; the size is also considered.
        plot_isocenter_sphere : bool
            Whether to plot the gantry + collimator isocenter size.
        plot_couch_iso : bool
            Whether to plot the couch-plane-only isocenter size.
            This will be zero if there are no images where the couch rotated.
        plot_coll_iso : bool
            Whether to plot the collimator-plane-only isocenter size.
            This is shown along the Z/Y plane only to differentiate from the couch iso visualization.
            The collimator plane is always normal to the gantry angle.
            This will be zero if there are no images where the collimator rotated.
        show_legend : bool
            Whether to show the legend.
        """
        limit = (
            viewbox_mm
            or max(
                np.abs(
                    (
                        self.bb_shift_vector.x,
                        self.bb_shift_vector.y,
                        self.bb_shift_vector.z,
                    )
                )
            )
            + self._bb_diameter
        )
        ax = plt.axes(projection="3d")
        _, relevant_images = self._get_images(
            axis=(Axis.REFERENCE, Axis.GB_COMBO, Axis.COLLIMATOR, Axis.GANTRY)
        )
        # we can represent the iso sphere as a BB object; the nominal object isn't used, just the BB size
        # the ray lines are what we want to plot as a sphere
        iso_sphere = BB(
            nominal_bb={
                "offset_left_mm": 0,
                "offset_in_mm": 0,
                "offset_up_mm": 0,
                "bb_diameter_mm": self._bb_diameter,
            },
            ray_lines=[image.cax_line_projection for image in relevant_images],
        )
        # plot the x,y,z origin lines
        x_line = Line(Point(-limit, 0, 0), Point(limit, 0, 0))
        x_line.plot2axes(ax, color="green", alpha=0.5)
        y_line = Line(Point(0, -limit, 0), Point(0, limit, 0))
        y_line.plot2axes(ax, color="green", alpha=0.5)
        z_line = Line(Point(0, 0, -limit), Point(0, 0, limit))
        z_line.plot2axes(
            ax, color="green", alpha=0.5, label="Determined isocenter (x,y,z)"
        )
        if plot_bb:
            iso_sphere.plot_measured(ax, color="cyan", alpha=0.6)
            # create an empty, fake line so we can add a label for the legend
            fake_line = Line(Point(0, 0, 0), Point(0, 0, 0))
            fake_line.plot2axes(ax, color="cyan", label=f"BB ({self._bb_diameter}mm)")
        if plot_isocenter_sphere:
            x, y, z = create_sphere_surface(
                radius=self.gantry_coll_iso_size / 2, center=Point(0, 0, 0)
            )
            ax.plot_surface(x, y, z, alpha=0.3, color="magenta")
            # create an empty, fake line so we can add a label for the legend
            fake_line = Line(Point(0, 0, 0), Point(0, 0, 0))
            fake_line.plot2axes(
                ax,
                color="magenta",
                label=f"Gantry + Coll Isosphere ({self.gantry_coll_iso_size:3.2f}mm)",
            )
        if plot_couch_iso:
            circle = plt.Circle(
                (0, 0),
                radius=self.couch_iso_size / 2,
                fill=True,
                color="yellow",
                alpha=0.4,
                label=f"Couch-only iso ({self.couch_iso_size:3.2f}mm)",
            )
            ax.add_patch(circle)
            art3d.pathpatch_2d_to_3d(circle, z=0, zdir="z")
        if plot_coll_iso:
            circle = plt.Circle(
                (0, 0),
                radius=self.collimator_iso_size / 2,
                fill=True,
                color="blue",
                alpha=0.4,
                label=f"Collimator-only iso ({self.collimator_iso_size:3.2f}mm)",
            )
            ax.add_patch(circle)
            art3d.pathpatch_2d_to_3d(circle, z=0, zdir="x")
        if show_legend:
            ax.legend()
        # set the limits of the 3D plot; they must be the same in all axes for equal aspect ratio
        ax.set(
            xlabel="X (mm), Right (+)",
            ylabel="Y (mm), In (+)",
            zlabel="Z (mm), Up (+)",
            title="Isocenter Visualization",
            ylim=[-limit, limit],
            xlim=[-limit, limit],
            zlim=[-limit, limit],
        )

        if show:
            plt.show()

    def plot_images(
        self,
        axis: Axis = Axis.GANTRY,
        show: bool = True,
        split: bool = False,
        **kwargs,
    ) -> (list[plt.Figure], list[str]):
        """Plot a grid of all the images acquired.

        Four columns are plotted with the titles showing which axis that column represents.

        Parameters
        ----------
        axis : {'Gantry', 'Collimator', 'Couch', 'GB Combo', 'GBP Combo', 'All'}
            The axis to plot.
        show : bool
            Whether to show the image.
        split : bool
            Whether to show/plot the images individually or as one large figure.
        """
        axis = convert_to_enum(axis, Axis)
        if not self._is_analyzed:
            raise ValueError("The set is not analyzed. Use .analyze() first.")

        # get axis images
        if axis == Axis.GANTRY:
            images = [
                image
                for image in self.images
                if image.variable_axis in (Axis.GANTRY, Axis.REFERENCE)
            ]
        elif axis == Axis.COLLIMATOR:
            images = [
                image
                for image in self.images
                if image.variable_axis in (Axis.COLLIMATOR, Axis.REFERENCE)
            ]
        elif axis == Axis.COUCH:
            images = [
                image
                for image in self.images
                if image.variable_axis in (Axis.COUCH, Axis.REFERENCE)
            ]
        elif axis == Axis.GB_COMBO:
            images = [
                image
                for image in self.images
                if image.variable_axis
                in (Axis.GB_COMBO, Axis.GANTRY, Axis.COLLIMATOR, Axis.REFERENCE)
            ]
        elif axis == Axis.GBP_COMBO:
            images = self.images

        # set the figsize if it wasn't passed
        if not kwargs.get("figsize"):
            dpi = 72
            width_px = 1080
            width_in = width_px / dpi
            if not split:
                max_num_images = math.ceil(len(images) / 4)
                height_in = (width_in / 4) * max_num_images
            else:
                height_in = width_in = 3
            kwargs["figsize"] = (width_in, height_in)

        figs = []
        names = []
        # create plots
        if not split:
            fig, axes = plt.subplots(nrows=max_num_images, ncols=4, **kwargs)
            for mpl_axis, wl_image in zip_longest(axes.flatten(), images):
                plot_image(wl_image, mpl_axis)

            # set titles
            fig.suptitle(f"{axis.value} images", fontsize=14, y=1)
            fig.tight_layout()
            figs.append(fig)
            names.append("image")
        else:
            for wl_image in images:
                fig, axes = plt.subplots(**kwargs)
                plot_image(wl_image, axes)
                figs.append(fig)
                names.append(str(wl_image))

        if show:
            plt.show()

        return figs, names

    def save_images(
        self, filename: str | BinaryIO, axis: Axis = Axis.GANTRY, **kwargs
    ) -> None:
        """Save the figure of `plot_images()` to file. Keyword arguments are passed to `matplotlib.pyplot.savefig()`.

        Parameters
        ----------
        filename : str
            The name of the file to save to.
        axis
            The axis to save.
        """
        self.plot_images(axis=axis, show=False)
        plt.savefig(filename, **kwargs)

    def save_images_to_stream(self, **kwargs) -> dict[str, io.BytesIO]:
        """Save the individual image plots to stream"""
        figs, names = self.plot_images(
            axis=Axis.GBP_COMBO, show=False, split=True
        )  # all images
        streams = [io.BytesIO() for _ in figs]
        for fig, stream in zip(figs, streams):
            fig.savefig(stream, **kwargs)
        return {name: stream for name, stream in zip(names, streams)}

    def plot_summary(self, show: bool = True, fig_size: tuple | None = None) -> None:
        """Plot a summary figure showing the gantry sag and wobble plots of the three axes."""
        if not self._is_analyzed:
            raise ValueError("The set is not analyzed. Use .analyze() first.")
        figsize = (11, 9) if fig_size is None else fig_size
        plt.figure(figsize=figsize)
        grid = (3, 6)
        gantry_sag_ax = plt.subplot2grid(grid, (0, 0), colspan=3)
        self._plot_deviation(Axis.GANTRY, gantry_sag_ax, show=False)
        epid_sag_ax = plt.subplot2grid(grid, (0, 3), colspan=3)
        self._plot_deviation(Axis.EPID, epid_sag_ax, show=False)
        if self._get_images((Axis.COLLIMATOR, Axis.REFERENCE))[0] > 1:
            coll_sag_ax = plt.subplot2grid(grid, (1, 0), colspan=3)
            self._plot_deviation(Axis.COLLIMATOR, coll_sag_ax, show=False)
        if self._get_images((Axis.COUCH, Axis.REFERENCE))[0] > 1:
            couch_sag_ax = plt.subplot2grid(grid, (1, 3), colspan=3)
            self._plot_deviation(Axis.COUCH, couch_sag_ax, show=False)

        for axis, axnum in zip((Axis.GANTRY, Axis.COLLIMATOR, Axis.COUCH), (0, 2, 4)):
            if self._get_images((axis, Axis.REFERENCE))[0] > 1:
                ax = plt.subplot2grid(grid, (2, axnum), colspan=2)
                self.plot_axis_images(axis=axis, ax=ax, show=False)
        if show:
            plt.tight_layout()
            plt.show()

    def save_summary(self, filename: str | BinaryIO, **kwargs) -> None:
        """Save the summary image."""
        self.plot_summary(show=False, fig_size=kwargs.pop("fig_size", None))
        plt.tight_layout()
        plt.savefig(filename, **kwargs)

    def results(self, as_list: bool = False) -> str:
        """Return the analysis results summary.

        Parameters
        ----------
        as_list : bool
            Whether to return as a list of strings vs single string. Pretty much for internal usage.
        """
        if not self._is_analyzed:
            raise ValueError("The set is not analyzed. Use .analyze() first.")
        num_gantry_imgs = self._get_images(axis=(Axis.GANTRY, Axis.REFERENCE))[0]
        num_gantry_coll_imgs = self._get_images(
            axis=(Axis.GANTRY, Axis.COLLIMATOR, Axis.GB_COMBO, Axis.REFERENCE)
        )[0]
        num_coll_imgs = self._get_images(axis=(Axis.COLLIMATOR, Axis.REFERENCE))[0]
        num_couch_imgs = self._get_images(axis=(Axis.COUCH, Axis.REFERENCE))[0]
        num_imgs = len(self.images)
        result = [
            "Winston-Lutz Analysis",
            "=================================",
            f"Number of images: {num_imgs}",
            f"Maximum 2D CAX->BB distance: {self.cax2bb_distance('max'):.2f}mm",
            f"Median 2D CAX->BB distance: {self.cax2bb_distance('median'):.2f}mm",
            f"Mean 2D CAX->BB distance: {self.cax2bb_distance('mean'):.2f}mm",
            f"Shift to iso: facing gantry, move BB: {self.bb_shift_instructions()}",
            f"Gantry 3D isocenter diameter: {self.gantry_iso_size:.2f}mm ({num_gantry_imgs}/{num_imgs} images considered)",
            f"Maximum Gantry RMS deviation (mm): {max(self.axis_rms_deviation((Axis.GANTRY, Axis.REFERENCE))):.2f}mm",
            f"Maximum EPID RMS deviation (mm): {max(self.axis_rms_deviation(Axis.EPID)):.2f}mm",
            f"Gantry+Collimator 3D isocenter diameter: {self.gantry_coll_iso_size:.2f}mm ({num_gantry_coll_imgs}/{num_imgs} images considered)",
            f"Collimator 2D isocenter diameter: {self.collimator_iso_size:.2f}mm ({num_coll_imgs}/{num_imgs} images considered)",
            f"Maximum Collimator RMS deviation (mm): {max(self.axis_rms_deviation((Axis.COLLIMATOR, Axis.REFERENCE))):.2f}",
            f"Couch 2D isocenter diameter: {self.couch_iso_size:.2f}mm ({num_couch_imgs}/{num_imgs} images considered)",
            f"Maximum Couch RMS deviation (mm): {max(self.axis_rms_deviation((Axis.COUCH, Axis.REFERENCE))):.2f}",
        ]
        if not as_list:
            result = "\n".join(result)
        return result

    def results_data(self, as_dict: bool = False) -> WinstonLutzResult | dict:
        """Present the results data and metadata as a dataclass or dict.
        The default return type is a dataclass."""
        if not self._is_analyzed:
            raise ValueError("The set is not analyzed. Use .analyze() first.")
        num_gantry_imgs = self._get_images(axis=(Axis.GANTRY, Axis.REFERENCE))[0]
        num_gantry_coll_imgs = self._get_images(
            axis=(Axis.GANTRY, Axis.COLLIMATOR, Axis.GB_COMBO, Axis.REFERENCE)
        )[0]
        num_coll_imgs = self._get_images(axis=(Axis.COLLIMATOR, Axis.REFERENCE))[0]
        num_couch_imgs = self._get_images(axis=(Axis.COUCH, Axis.REFERENCE))[0]

        individual_image_data = [i.results_data(as_dict=as_dict) for i in self.images]
        if as_dict:
            # convert classes to dicts; little wonky but we have to get it through
            # to radmachine and we want to dynamically convert classes to dicts
            for img in individual_image_data:
                for key, value in img.items():
                    try:
                        img[key] = value.__dict__
                    except AttributeError:
                        pass

        data = WinstonLutzResult(
            num_total_images=len(self.images),
            num_gantry_images=num_gantry_imgs,
            num_coll_images=num_coll_imgs,
            num_gantry_coll_images=num_gantry_coll_imgs,
            num_couch_images=num_couch_imgs,
            max_2d_cax_to_bb_mm=self.cax2bb_distance("max"),
            median_2d_cax_to_bb_mm=self.cax2bb_distance("median"),
            mean_2d_cax_to_bb_mm=self.cax2bb_distance("mean"),
            max_2d_cax_to_epid_mm=self.cax2epid_distance("max"),
            median_2d_cax_to_epid_mm=self.cax2epid_distance("median"),
            mean_2d_cax_to_epid_mm=self.cax2epid_distance("mean"),
            coll_2d_iso_diameter_mm=self.collimator_iso_size,
            couch_2d_iso_diameter_mm=self.couch_iso_size,
            gantry_3d_iso_diameter_mm=self.gantry_iso_size,
            gantry_coll_3d_iso_diameter_mm=self.gantry_coll_iso_size,
            max_gantry_rms_deviation_mm=max(
                self.axis_rms_deviation(axis=(Axis.GANTRY, Axis.REFERENCE))
            ),
            max_coll_rms_deviation_mm=max(
                self.axis_rms_deviation(axis=(Axis.COLLIMATOR, Axis.REFERENCE))
            ),
            max_couch_rms_deviation_mm=max(
                self.axis_rms_deviation(axis=(Axis.COUCH, Axis.REFERENCE))
            ),
            max_epid_rms_deviation_mm=max(self.axis_rms_deviation(axis=Axis.EPID)),
            image_details=individual_image_data,
            keyed_image_details=self._generate_keyed_images(individual_image_data),
        )
        if as_dict:
            return dataclasses.asdict(data)
        return data

    def _generate_keyed_images(
        self, individual_image_data: list[WinstonLutz2D] | dict
    ) -> dict:
        """Generate a dict where each key is based on the axes values and the key is an image. Used in the results_data method.
        We can't do a simple dict comprehension because we may have duplicate axes sets. We pass individual data
        because we may have already converted to a dict; we don't want to do that again.
        """
        data = {}
        for img_idx, img in enumerate(self.images):
            key = f"G{img.gantry_angle}B{img.collimator_angle}P{img.couch_angle}"
            suffix = ""
            idx = 1
            while key + suffix in data.keys():
                suffix = f"_{idx}"
                idx += 1
            data[key + suffix] = individual_image_data[img_idx]
        return data

    def publish_pdf(
        self,
        filename: str,
        notes: str | list[str] | None = None,
        open_file: bool = False,
        metadata: dict | None = None,
        logo: Path | str | None = None,
    ):
        """Publish (print) a PDF containing the analysis, images, and quantitative results.

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
        if not self._is_analyzed:
            raise ValueError("The set is not analyzed. Use .analyze() first.")
        plt.ioff()
        title = "Winston-Lutz Analysis"
        canvas = pdf.PylinacCanvas(
            filename, page_title=title, metadata=metadata, logo=logo
        )
        text = self.results(as_list=True)
        canvas.add_text(text=text, location=(7, 25.5))
        # draw summary image on 1st page
        data = io.BytesIO()
        self.save_summary(data, fig_size=(8, 8))
        canvas.add_image(image_data=data, location=(2, 3), dimensions=(16, 16))
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 4.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 4))
        # add more pages showing individual axis images
        for ax in (
            Axis.GANTRY,
            Axis.COLLIMATOR,
            Axis.COUCH,
            Axis.GB_COMBO,
            Axis.GBP_COMBO,
        ):
            if self._contains_axis_images(ax):
                canvas.add_new_page()
                data = io.BytesIO()
                self.save_images(data, axis=ax)
                canvas.add_image(data, location=(2, 7), dimensions=(18, 18))

        canvas.finish()

        if open_file:
            webbrowser.open(filename)

    def _contains_axis_images(self, axis: Axis = Axis.GANTRY) -> bool:
        """Return whether or not the set of WL images contains images pertaining to a given axis"""
        return any(True for image in self.images if image.variable_axis in (axis,))


class WinstonLutz2DMultiTarget(WinstonLutz2D):
    """A 2D image of a WL delivery, but where multiple BBs are in use."""

    detection_conditions = [is_round, is_symmetric, is_modest_size]
    field_conditions = [is_square, is_right_square_size]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.flipud()  # restore to original view; vanilla WL may need to revert

    def as_analyzed(self, bb_location: dict) -> WinstonLutz2DMultiTarget:
        """Analyze the image of the multi-BB setup. We return a copy of the
        WL image because we analyze images more than once. Each "analyzed" image
        is really the analysis of a BB/image combo.

        Parameters
        ----------
        bb_location
            An iterable of dictionaries. Each dict contains keys for the offsets and size of the BB in mm.
            Use the ``BBArrangement`` class as a guide.
        """
        self.field_cax, self._rad_field_bounding_box = self._find_field_centroid(
            bb_location
        )
        self.bb = self._find_bb(bb_location)
        self._is_analyzed = True
        return copy.deepcopy(self)

    def plot(
        self, ax: plt.Axes | None = None, show: bool = True, clear_fig: bool = False
    ) -> plt.Axes:
        ax = super(LinacDicomImage, self).plot(ax=ax, show=False, clear_fig=clear_fig)
        ax.plot(self.field_cax.x, self.field_cax.y, "gs", ms=8)
        ax.set_ylim([self._rad_field_bounding_box[0], self._rad_field_bounding_box[2]])
        ax.set_xlim([self._rad_field_bounding_box[1], self._rad_field_bounding_box[3]])
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_title("\n".join(wrap(Path(self.path).name, 30)), fontsize=10)
        ax.set_xlabel(
            f"G={self.gantry_angle:.0f}, B={self.collimator_angle:.0f}, P={self.couch_angle:.0f}"
        )
        ax.set_ylabel(f"CAX to BB: {self.cax2bb_distance:3.2f}mm")
        if show:
            plt.show()
        return ax

    def _nominal_point(self, bb: dict) -> Point:
        """Calculate the expected point position in 2D"""
        shift_y_mm = bb_projection_long(
            offset_in=bb["offset_in_mm"],
            offset_up=bb["offset_up_mm"],
            offset_left=bb["offset_left_mm"],
            sad=self.sad,
            gantry=self.gantry_angle,
        )
        shift_x_mm = bb_projection_gantry_plane(
            offset_left=bb["offset_left_mm"],
            offset_up=bb["offset_up_mm"],
            sad=self.sad,
            gantry=self.gantry_angle,
        )
        # unlike vanilla WL, the field can be asymmetric, so use center of image
        expected_y = self.epid.y - shift_y_mm * self.dpmm
        expected_x = self.epid.x + shift_x_mm * self.dpmm
        return Point(x=expected_x, y=expected_y)

    def _find_field_centroid(self, location: dict) -> tuple[Point, list]:
        """Find the centroid of the radiation field based on a 50% height threshold.
        This applies the field detection conditions and also a nearness condition.

        Returns
        -------
        p
            The CAX point location.
        edges
            The bounding box of the field, plus a small margin.
        """
        min, max = np.percentile(self.array, [5, 99.9])
        threshold_img = self.as_binary((max - min) / 2 + min)
        filled_img = ndimage.binary_fill_holes(threshold_img)
        labeled_arr, num_roi = ndimage.label(filled_img)
        regions = measure.regionprops(labeled_arr)
        field_candidates = [
            all(
                condition(
                    region,
                    dpmm=self.dpmm,
                    rad_size=location["rad_size_mm"],
                    shape=labeled_arr.shape,
                )
                for condition in self.field_conditions
            )
            for region in regions
        ]
        if not any(field_candidates):
            raise ValueError(
                "Did not find an ROI that looked like the radiation field. Ensure the rad_size parameter is correct."
            )
        field_localized = [
            self.location_near_nominal(region, location) for region in regions
        ]
        if not any(field_localized):
            raise ValueError("Did not find the radiation field where it was expected.")
        else:
            region_idx = [
                idx
                for idx, _ in enumerate(regions)
                if field_candidates[idx] and field_localized[idx]
            ][0]
            field_roi = regions[region_idx]
        bbox = list(field_roi.bbox)
        bbox[0] -= 40
        bbox[1] -= 40
        bbox[2] += 40
        bbox[3] += 40
        return Point(field_roi.centroid[1], field_roi.centroid[0]), bbox

    def _find_bb(self, bb_of_interest: dict) -> Point:
        """Find the specific BB based on the arrangement rather than a single one. This is in local pixel coordinates"""
        # get initial starting conditions
        bb_diameter = bb_of_interest["bb_size_mm"]
        window = bb_of_interest["rad_size_mm"]
        expected_position = self._nominal_point(bb_of_interest)
        expected_position_mm = expected_position / self.dpmm
        bb_tolerance_mm = self._calculate_bb_tolerance(bb_diameter)
        center = self.compute(
            metrics=SizedDiskLocator.from_physical(
                expected_position_mm=expected_position_mm,
                search_window_mm=(window, window),
                radius_mm=bb_diameter / 2,
                radius_tolerance_mm=bb_tolerance_mm,
                invert=True,  # only MV images are involved for MTWL AFAICT
                detection_conditions=self.detection_conditions,
            )
        )
        return center

    def location_near_nominal(self, region: RegionProperties, location: dict) -> bool:
        """Determine whether the given BB ROI is near where the BB is expected to be"""
        # since we are dealing with images at the isoplane we have to calculate the expected position
        # of the BB at that plane from the 3D coordinates
        if region.area < 5:
            return False  # skip single or very small pixel regions
        expected = self._nominal_point(location)
        near_y = math.isclose(expected.y, region.centroid[0], abs_tol=5 * self.dpmm)
        near_x = math.isclose(expected.x, region.centroid[1], abs_tol=5 * self.dpmm)
        return near_y and near_x

    def results_data(self, as_dict: bool = False) -> WinstonLutz2DResult | dict:
        raise NotImplementedError(
            "Results data is not available for a multi-bb 2D WL image"
        )


class WinstonLutzMultiTargetMultiField(WinstonLutz):
    images: Sequence[WinstonLutz2DMultiTarget]  #:
    analyzed_images: dict[str, list[WinstonLutz2DMultiTarget]]  #:
    image_type = WinstonLutz2DMultiTarget
    bb_arrangement: Sequence[dict]  #:

    def __init__(self, *args, **kwargs):
        """We cannot yet handle non-0 couch angles so we drop them. Analysis fails otherwise"""
        super().__init__(*args, **kwargs)
        orig_length = len(self.images)
        self.images = [
            i for i in self.images if is_close(i.couch_angle, [0, 360], delta=5)
        ]
        new_length = len(self.images)
        if new_length != orig_length:
            print(
                f"Non-zero couch angles not yet allowed. Dropped {orig_length-new_length} images"
            )

    @classmethod
    def from_demo_images(cls):
        """Instantiate using the demo images."""
        demo_file = retrieve_demo_file(name="mt_mf_wl.zip")
        return cls.from_zip(demo_file)

    @staticmethod
    def run_demo():
        """Run the Winston-Lutz MT MF demo, which loads the demo files, prints results, and plots a summary image."""
        arrangement = (
            {
                "name": "Iso",
                "offset_left_mm": 0,
                "offset_up_mm": 0,
                "offset_in_mm": 0,
                "bb_size_mm": 5,
                "rad_size_mm": 20,
            },
            {
                "name": "Left,Down,In",
                "offset_left_mm": 20,
                "offset_up_mm": -20,
                "offset_in_mm": 60,
                "bb_size_mm": 5,
                "rad_size_mm": 20,
            },
        )
        wl = WinstonLutzMultiTargetMultiField.from_demo_images()
        wl.analyze(bb_arrangement=arrangement)
        print(wl.results())
        wl.plot_images()

    def analyze(self, bb_arrangement: Sequence[dict]):
        """Analyze the WL images.

        Parameters
        ----------
        bb_arrangement
            The arrangement of the BBs in the phantom. A dict with offset and BB size keys. See the ``BBArrangement`` class for
            keys and syntax.
        """
        self.analyzed_images = {}
        self.bb_arrangement = bb_arrangement
        for idx, bb in enumerate(bb_arrangement):
            image_set = []
            for img in self.images:
                try:
                    image_set.append(img.as_analyzed(bb))
                except ValueError:
                    pass
            if not image_set:
                raise ValueError(f"Did not find any field/bb pairs for bb: {bb}")
            self.analyzed_images[BBArrangement.to_human(bb)] = image_set
        self._is_analyzed = True

    def plot_images(self, show: bool = True, **kwargs) -> (list[plt.Figure], list[str]):
        """Make a plot for each BB. Each plot contains the analysis of that BB on each image
        it was found."""
        figs, names = [], []
        figsize = kwargs.pop("figsize", None) or (8, 8)
        for bb, img_set in self.analyzed_images.items():
            rows = len(img_set) // 3 + 1
            fig, axes = plt.subplots(nrows=rows, ncols=3, figsize=figsize, **kwargs)
            for mpl_axis, wl_image in zip_longest(axes.flatten(), img_set):
                plot_image(wl_image, mpl_axis)

            # set titles
            fig.suptitle(f"BB {bb}", fontsize=14, y=1)
            fig.tight_layout()
            figs.append(fig)
            names.append(f"BB_{bb}")
        if show:
            plt.show()
        return figs, names

    def save_images(self, prefix: str = "", **kwargs):
        """Save the figure of `plot_images()` to file as PNG. Keyword arguments are passed to `matplotlib.pyplot.savefig()`.

        Parameters
        ----------
        prefix : str
            The prefix name of the file to save to. The BB name is appended to the prefix.
        """
        figs, names = self.plot_images(show=False, **kwargs)
        for fig, name in zip(figs, names):
            fig.savefig(prefix + "_" + str(name) + ".png", **kwargs)

    def save_images_to_stream(self, **kwargs) -> dict[str, io.BytesIO]:
        """Save the individual image plots to stream"""
        figs, names = self.plot_images(show=False, **kwargs)
        streams = [io.BytesIO() for _ in figs]
        for fig, stream in zip(figs, streams):
            fig.savefig(stream, **kwargs)
        return {name: stream for name, stream in zip(names, streams)}

    def cax2bb_distance(self, bb: str, metric: str = "max") -> float:
        """The distance in mm between the CAX and BB for all images according to the given metric.

        Parameters
        ----------
        metric : {'max', 'median', 'mean'}
            The metric of distance to use.
        bb : str
            The BB to analyze
        """
        if metric == "max":
            return max(image.cax2bb_distance for image in self.analyzed_images[bb])
        elif metric == "median":
            return statistics.median(
                [image.cax2bb_distance for image in self.analyzed_images[bb]]
            )
        elif metric == "mean":
            return statistics.mean(
                [image.cax2bb_distance for image in self.analyzed_images[bb]]
            )

    def results_data(
        self, as_dict: bool = False
    ) -> WinstonLutzMultiTargetMultiFieldResult | dict:
        """Present the results data and metadata as a dataclass or dict.
        The default return type is a dataclass."""
        if not self._is_analyzed:
            raise ValueError("The set is not analyzed. Use .analyze() first.")

        data = WinstonLutzMultiTargetMultiFieldResult(
            num_total_images=len(self.images),
            max_2d_field_to_bb_mm=self.max_bb_deviation_2d,
            mean_2d_field_to_bb_mm=self.mean_bb_deviation_2d,
            median_2d_field_to_bb_mm=self.median_bb_deviation_2d,
            bb_maxes={
                bb: self.cax2bb_distance(bb) for bb in self.analyzed_images.keys()
            },
            bb_arrangement=self.bb_arrangement,
        )
        if as_dict:
            return dataclasses.asdict(data)
        return data

    def plot_summary(self, show: bool = True, fig_size: tuple | None = None):
        raise NotImplementedError("Not yet implemented")

    def plot_axis_images(
        self, axis: Axis = Axis.GANTRY, show: bool = True, ax: plt.Axes | None = None
    ):
        raise NotImplementedError("Not yet implemented")

    @property
    def max_bb_deviation_2d(self) -> float:
        """The maximum distance from any measured BB to its nominal position"""
        dists = []
        for bb, images in self.analyzed_images.items():
            dists.append(self.cax2bb_distance(bb))
        return max(dists)

    @property
    def mean_bb_deviation_2d(self) -> float:
        """The mean distance from any measured BB to its nominal position"""
        dists = []
        for bb, images in self.analyzed_images.items():
            dists.append(self.cax2bb_distance(bb))
        return statistics.mean(dists)

    @property
    def median_bb_deviation_2d(self) -> float:
        """The median distance from any measured BB to its nominal position"""
        dists = []
        for bb, images in self.analyzed_images.items():
            dists.append(self.cax2bb_distance(bb))
        return statistics.median(dists)

    def results(self, as_list: bool = False) -> str:
        """Return the analysis results summary.

        Parameters
        ----------
        as_list : bool
            Whether to return as a list of strings vs single string. Pretty much for internal usage.
        """
        if not self._is_analyzed:
            raise ValueError("The set is not analyzed. Use .analyze() first.")
        num_imgs = len(self.images)
        result = [
            "Winston-Lutz Multi-Target Multi-Field Analysis",
            "==============================================",
            f"Number of images: {num_imgs}",
            "",
            "2D distances",
            "============",
            f"Max 2D distance of any BB: {self.max_bb_deviation_2d:.2f} mm",
            f"Mean 2D distance of any BB: {self.mean_bb_deviation_2d:.2f} mm",
            f"Median 2D distance of any BB: {self.median_bb_deviation_2d:.2f} mm",
            "",
        ]
        bb_descriptions = [
            [idx, bb] for idx, bb in enumerate(self.analyzed_images.keys())
        ]
        result += tabulate(bb_descriptions, headers=["BB #", "Description"]).split("\n")
        result += [
            "",
        ]

        # construct the image -> bb table
        # we use abbreviations and truncations so that the table will more likely fit the PDF
        data = {
            "Image": [],
            "G": [],
            "Co": [],
            "Ch": [],
        }
        bbs = {f"BB #{idx}": [] for idx in range(len(self.analyzed_images.keys()))}
        data.update(bbs)
        for img in self.images:
            data["Image"].append(
                img.base_path[-20:]
            )  # textwrap doesn't work here because files may be all one word
            data["G"].append(f"{img.gantry_angle:.1f}")
            data["Co"].append(f"{img.collimator_angle:.1f}")
            data["Ch"].append(f"{img.couch_angle:.1f}")
            for bb_idx, (bb, img_set) in enumerate(self.analyzed_images.items()):
                has_value = False
                for sub_img in img_set:
                    if img.base_path == sub_img.base_path:
                        data[f"BB #{bb_idx}"].append(f"{sub_img.cax2bb_distance:.2f}")
                        has_value = True
                if not has_value:
                    data[f"BB #{bb_idx}"].append("---")
        result += tabulate(data, headers="keys").split("\n")
        if not as_list:
            result = "\n".join(result)
        return result

    def publish_pdf(
        self,
        filename: str,
        notes: str | list[str] | None = None,
        open_file: bool = False,
        metadata: dict | None = None,
        logo: Path | str | None = None,
    ):
        """Publish (print) a PDF containing the analysis, images, and quantitative results.

        Parameters
        ----------
        filename : (str, file-like object)
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
        if not self._is_analyzed:
            raise ValueError("The set is not analyzed. Use .analyze() first.")
        plt.ioff()
        title = "Winston-Lutz Multi-BB Analysis"
        canvas = pdf.PylinacCanvas(
            filename,
            page_title=title,
            metadata=metadata,
            logo=logo,
            metadata_location=(15, 25.5),
        )
        text = self.results(as_list=True)
        canvas.add_text(text=text, location=(1, 25.5), font="Courier")
        # draw summary image on 1st page
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 4.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 4))
        # plot each BB's images
        bb_streams = self.save_images_to_stream()
        for stream in bb_streams.values():
            canvas.add_new_page()
            canvas.add_image(stream, location=(2, 7), dimensions=(18, 18))
        canvas.finish()

        if open_file:
            webbrowser.open(filename)


def max_distance_to_lines(p, lines: Iterable[Line]) -> float:
    """Calculate the maximum distance to any line from the given point."""
    point = Point(p[0], p[1], p[2])
    return max(line.distance_to(point) for line in lines)


def bb_projection_long(
    offset_in: float, offset_up: float, offset_left: float, sad: float, gantry: float
) -> float:
    """Calculate the isoplane projection in the sup/inf/longitudinal direction in mm"""
    # the divergence of the beam causes the BB to be closer or further depending on the
    # up/down position, left/right position and gantry angle
    addtl_long_shift_cos = (
        offset_up * offset_in / (sad - cos(gantry) * offset_up) * cos(gantry)
    )
    addtl_left_shift_sin = (
        offset_left * offset_in / (sad + sin(gantry) * offset_left) * -sin(gantry)
    )
    return offset_in + addtl_long_shift_cos + addtl_left_shift_sin


def bb_projection_gantry_plane(
    offset_left: float, offset_up: float, sad: float, gantry: float
) -> float:
    """Calculate the isoplane projection in the plane of gantry rotation (X/Z)"""
    addtl_left_shift = (
        -offset_up * offset_left / (sad + cos(gantry) * offset_up) * abs(cos(gantry))
    )
    addtl_up_shift = (
        offset_left * offset_up / (sad + sin(gantry) * offset_left) * abs(sin(gantry))
    )
    return (
        offset_up * -sin(gantry)
        + addtl_up_shift
        + offset_left * -cos(gantry)
        + addtl_left_shift
    )


def _bb_projection_with_rotation(
    offset_left: float,
    offset_up: float,
    offset_in: float,
    gantry: float,
    couch: float = 0,
    sad: float = 1000,
) -> (float, float):
    """Calculate the isoplane projection onto the panel at the given SSD.

    This function applies a rotation around the gantry plane (X/Z) to the
    ball bearing (BB) position and calculates its projection onto the isocentre plane in the beam's eye view.

    Could be used to calculate couch rotations, but not validated yet.

    Args
    ----
    offset_left (float): The BB position in the left/right direction.
    offset_up (float): The BB position in the superior/inferior direction.
    offset_in (float): The BB position in the anterior/posterior direction.
    gantry (float): The gantry angle in degrees.
    couch (float, optional): The couch angle in degrees. Defaults to 0.
    sad (float, optional): The source-to-axis distance in mm. Defaults to 1000.

    Returns
    -------
    np.ndarray: The projection of the BB onto the panel at the given SSD.
        The array has shape (2,) where the first element is the projection in the
        left/right direction and the second element is the projection in the
        superior/inferior direction.
    """
    # Define the BB positions in the patient coordinate system (ap, lr, si)
    bb_positions = np.array([offset_up, offset_left, offset_in])

    # Apply the rotation matrix to the BB positions
    collimator = 0  # Collimator doesn't change positional projection onto panel
    rotation_matrix = Rotation.from_euler(
        "xyz", [couch, collimator, gantry], degrees=True
    )
    rotated_positions = rotation_matrix.apply(bb_positions)

    # Calculate the projection onto the panel at the given SSD
    bb_magnification = sad / (
        sad - rotated_positions[0]
    )  # Distance from source to panel
    imager_projection = (
        np.array([rotated_positions[1], rotated_positions[2]]) * bb_magnification
    )
    return imager_projection[0], imager_projection[1]
