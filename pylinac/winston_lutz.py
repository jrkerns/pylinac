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

import dataclasses
import enum
import io
import math
import os.path as osp
import statistics
import tempfile
import webbrowser
from collections.abc import Iterable, Sequence
from functools import cached_property
from itertools import zip_longest
from pathlib import Path
from textwrap import wrap
from typing import BinaryIO

import argue
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import art3d
from plotly import graph_objects as go
from py_linq import Enumerable
from pydantic import BaseModel, Field
from scipy import ndimage, optimize
from scipy.ndimage import zoom
from scipy.spatial.transform import Rotation
from skimage.measure._regionprops import RegionProperties
from tabulate import tabulate

from .core import image, pdf
from .core.array_utils import array_to_dicom
from .core.decorators import lru_cache
from .core.geometry import (
    Line,
    Point,
    PointSerialized,
    Vector,
    VectorSerialized,
    cos,
    sin,
)
from .core.image import DicomImageStack, is_image, tiff_to_dicom
from .core.io import TemporaryZipDirectory, get_url, retrieve_demo_file
from .core.plotly_utils import add_horizontal_line, add_title, add_vertical_line
from .core.scale import MachineScale, convert
from .core.utilities import (
    QuaacDatum,
    QuaacMixin,
    ResultBase,
    ResultsDataMixin,
    convert_to_enum,
    is_close_degrees,
)
from .core.warnings import capture_warnings
from .metrics.features import (
    is_right_circumference,
    is_right_size_bb,
    is_round,
    is_solid,
    is_symmetric,
)
from .metrics.image import GlobalSizedFieldLocator, SizedDiskLocator

BB_ERROR_MESSAGE = (
    "Unable to locate the BB. Make sure the field edges do not obscure the BB, that there are no artifacts in the images, that the 'bb_size' parameter is close to reality, "
    "and that the BB is near the center (within 2cm). If this is a large-field image or kV image try setting 'low_density_bb' to True."
)


class BBConfig(BaseModel):
    name: str
    offset_left_mm: float
    offset_up_mm: float
    offset_in_mm: float
    bb_size_mm: float
    rad_size_mm: float

    def to_human(self) -> str:
        """Convert one BB location to a human-readable str"""
        lr = "Left" if self.offset_left_mm >= 0 else "Right"
        ud = "Up" if self.offset_up_mm >= 0 else "Down"
        io = "In" if self.offset_in_mm >= 0 else "Out"
        return f"{lr} {abs(self.offset_left_mm)}mm, {ud} {abs(self.offset_up_mm)}mm, {io} {abs(self.offset_in_mm)}mm"


class BBArrangement:
    """Presets for multi-target phantoms."""

    # a BB at iso; represents the simplest case
    ISO = (
        BBConfig(
            name="Iso",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=0,
            bb_size_mm=5,  # overridden later dynamically
            rad_size_mm=20,  # not used;
        ),
    )
    # ISOCAL = (
    #     BBConfig(
    #         name="1",
    #         offset_left_mm=0,
    #         offset_up_mm=-170,
    #         offset_in_mm=-30,
    #         bb_size_mm=5,
    #         rad_size_mm=20,
    #     ),
    # {
    #     'name': '2',
    #     'offset_left_mm': -170,
    #     'offset_up_mm': 0,
    #     'offset_in_mm': -45,
    #     "bb_size_mm": 5,
    #     'rad_size_mm': 20,
    # }
    # )
    # locations: https://www.postersessiononline.eu/173580348_eu/congresos/ESTRO2020/aula/-PO_1320_ESTRO2020.pdf
    SNC_MULTIMET = (
        BBConfig(
            name="Iso",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=0,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
        BBConfig(
            name="1",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=30,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
        BBConfig(
            name="2",
            offset_left_mm=-30,
            offset_up_mm=0,
            offset_in_mm=15,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
        BBConfig(
            name="3",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=-30,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
        BBConfig(
            name="4",
            offset_left_mm=30,
            offset_up_mm=0,
            offset_in_mm=-50,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
        BBConfig(
            name="5",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=-70,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
    )
    DEMO = (
        BBConfig(
            name="Iso",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=0,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
        BBConfig(
            name="1",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=30,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
        BBConfig(
            name="2",
            offset_left_mm=-30,
            offset_up_mm=0,
            offset_in_mm=15,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
        BBConfig(
            name="3",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=-30,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
        BBConfig(
            name="4",
            offset_left_mm=30,
            offset_up_mm=0,
            offset_in_mm=-50,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
        BBConfig(
            name="5",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=-70,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
    )

    @staticmethod
    def to_human(arrangement: dict) -> str:
        """Convert one BB location to a human-readable str"""
        a = arrangement
        lr = "Left" if a["offset_left_mm"] >= 0 else "Right"
        ud = "Up" if a["offset_up_mm"] >= 0 else "Down"
        io = "In" if a["offset_in_mm"] >= 0 else "Out"
        return f"'{a['name']}': {lr} {abs(a['offset_left_mm'])}mm, {ud} {abs(a['offset_up_mm'])}mm, {io} {abs(a['offset_in_mm'])}mm"


@dataclasses.dataclass
class BBFieldMatch:
    """A match of a BB and field to an expected arrangement position. I.e. the nominal BB, measured BB, and field.
    This can calculate distances, create backprojections, etc for a single BB/Field."""

    epid: Point
    field: Point
    bb: Point
    dpmm: float
    gantry_angle: float
    couch_angle: float
    sad: float

    @property
    def field_epid_vector_mm(self) -> Vector:
        """The vector from the field CAX to the EPID center *IN COORDINATE SPACE*."""
        v = (self.field - self.epid) / self.dpmm
        v.y = -v.y  # invert the y-axis; positive is down in image space but negative in coordinate space
        return v

    @property
    def bb_field_vector_mm(self) -> Vector:
        """The vector from the BB to the field CAX *IN COORDINATE SPACE*."""
        v = (self.bb - self.field) / self.dpmm
        v.y = -v.y  # invert the y-axis; positive is down in image space but negative in coordinate space
        return v

    @property
    def bb_epid_vector_mm(self) -> Vector:
        """The vector from the BB to the field CAX *IN COORDINATE SPACE*."""
        v = (self.bb - self.epid) / self.dpmm
        v.y = -v.y  # invert the y-axis; positive is down in image space but negative in coordinate space
        return v

    @property
    def bb_field_distance_mm(self) -> float:
        """The distance from the BB to the field CAX in mm."""
        return self.field.distance_to(self.bb) / self.dpmm

    @property
    def bb_epid_distance_mm(self) -> float:
        """The distance from the BB to the EPID center in mm."""
        return self.epid.distance_to(self.bb) / self.dpmm

    @property
    def field_epid_distance_mm(self) -> float:
        """The distance from the field CAX to the EPID center in mm."""
        return self.epid.distance_to(self.field) / self.dpmm

    @property
    def bb_to_field_projection(self):
        """The projection from the BB to the field. Used by vanilla WL
        to determine the gantry, coll, couch isocenter size because there the BB is the reference point.

        Returns
        -------
        Line
            The virtual line in space made by the BB.
        """
        return straight_ray(self.bb_field_vector_mm, self.gantry_angle)


class BB3D:
    """A representation of a BB in 3D space"""

    def __repr__(self):
        return self.nominal_bb_position

    def __init__(
        self,
        bb_config: BBConfig,
        bb_matches: Sequence[BBFieldMatch, ...],
        scale: MachineScale,
    ):
        self.bb_config = bb_config
        self.matches = bb_matches
        self.scale = scale

    @cached_property
    def measured_bb_position(self) -> Point:
        """The 3D measured position of the BB based on rotation matrices and gantry/couch angles."""
        xs = [m.bb_epid_vector_mm.x for m in self.matches]
        ys = [m.bb_epid_vector_mm.y for m in self.matches]
        thetas = [m.gantry_angle for m in self.matches]
        phis = [m.couch_angle for m in self.matches]
        vector = solve_3d_position_from_2d_planes(
            xs=xs, ys=ys, thetas=thetas, phis=phis, scale=self.scale
        )
        # vectors and points are effectively the same thing here but we convert to a point for clarity
        return Point(x=vector.x, y=vector.y, z=vector.z)

    @cached_property
    def nominal_bb_position(self) -> Point:
        """The nominal location of the BB in MM in coordinate space"""
        return Point(
            x=-self.bb_config.offset_left_mm,
            y=self.bb_config.offset_in_mm,
            z=self.bb_config.offset_up_mm,
        )

    @cached_property
    def measured_field_position(self) -> Point:
        """The position of the field CAXs in 3D space"""
        xs = [m.field_epid_vector_mm.x for m in self.matches]
        ys = [m.field_epid_vector_mm.y for m in self.matches]
        thetas = [m.gantry_angle for m in self.matches]
        phis = [m.couch_angle for m in self.matches]
        vector = solve_3d_position_from_2d_planes(
            xs=xs, ys=ys, thetas=thetas, phis=phis, scale=self.scale
        )
        # vectors and points are effectively the same thing here but we convert to a point for clarity
        return Point(x=vector.x, y=vector.y, z=vector.z)

    def plotly_nominal(self, fig: go.Figure, color: str, **kwargs) -> None:
        x, y, z = create_sphere_surface(
            radius=self.bb_config.bb_size_mm / 2, center=self.nominal_bb_position
        )
        fig.add_surface(
            x=x,
            y=y,
            z=z,
            name=f"Nominal BB - {self.bb_config.name}",
            showscale=False,
            colorscale=[[0, color], [1, color]],
            showlegend=True,
            **kwargs,
        )

    def plot_nominal(self, axes: plt.Axes, color: str, **kwargs):
        """Plot the BB nominal position"""
        x, y, z = create_sphere_surface(
            radius=self.bb_config.bb_size_mm / 2, center=self.nominal_bb_position
        )
        axes.plot_surface(x, y, z, color=color, **kwargs)

    def plotly_measured(self, fig: go.Figure, color: str, **kwargs):
        """Plot the BB measured position"""
        x, y, z = create_sphere_surface(
            radius=self.bb_config.bb_size_mm / 2, center=self.measured_bb_position
        )
        fig.add_surface(
            x=x,
            y=y,
            z=z,
            name=f"Measured BB - {self.bb_config.name}",
            showscale=False,
            colorscale=[[0, color], [1, color]],
            showlegend=True,
            **kwargs,
        )

    def plot_measured(self, axes: plt.Axes, color: str, **kwargs):
        """Plot the BB measured position"""
        x, y, z = create_sphere_surface(
            radius=self.bb_config.bb_size_mm / 2, center=self.measured_bb_position
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


class WinstonLutz2DResult(ResultBase):
    variable_axis: str = Field(description="The axis that varied in the image.")
    bb_location: PointSerialized = Field(
        description="The location of the BB in the image as a Point in pixels."
    )
    cax2epid_vector: VectorSerialized = Field(
        description="The vector (in Cartesian coordinates) from the field CAX to the EPID center in mm."
    )
    cax2epid_distance: float = Field(
        description="The distance from the field CAX to the EPID center in mm.",
        title="Scalar distance from CAX to EPID center (mm)",
    )
    cax2bb_vector: VectorSerialized = Field(
        description="The vector (in Cartesian coordinates) from the field CAX to the BB in mm."
    )
    cax2bb_distance: float = Field(
        description="The scalar distance from the field CAX to the BB in mm.",
        title="Scalar distance from CAX to BB (mm)",
    )
    field_cax: PointSerialized = Field(
        description="The location of the field CAX in the image as a Point in pixels."
    )


class WinstonLutzResult(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    max_2d_cax_to_bb_mm: float = Field(
        description="The maximum 2D distance from the field CAX to the BB across all images analyzed in mm.",
        title="Max scalar in-plane distance from BB to CAX (mm)",
    )
    median_2d_cax_to_bb_mm: float = Field(
        description="The median 2D distance from the field CAX to the BB across all images analyzed in mm.",
        title="Median absolute scalar in-plane distance from BB to CAX (mm)",
    )
    mean_2d_cax_to_bb_mm: float = Field(
        description="The mean 2D distance from the field CAX to the BB across all images analyzed in mm.",
        title="Mean absolute scalar in-plane distance from BB to CAX (mm)",
    )
    max_2d_cax_to_epid_mm: float = Field(
        description="The maximum 2D distance from the field CAX to the EPID center across all images analyzed in mm.",
        title="Max scalar in-plane distance from EPID to CAX (mm)",
    )
    median_2d_cax_to_epid_mm: float = Field(
        description="The median 2D distance from the field CAX to the EPID center across all images analyzed in mm.",
        title="Median absolute scalar in-plane distance from EPID to CAX (mm)",
    )
    mean_2d_cax_to_epid_mm: float = Field(
        description="The mean 2D distance from the field CAX to the EPID center across all images analyzed in mm.",
        title="Mean absolute scalar in-plane distance from EPID to CAX (mm)",
    )
    gantry_3d_iso_diameter_mm: float = Field(
        description="The 3D isocenter diameter **of the gantry axis only** as determined by the gantry images in mm. This uses backprojection lines of the field center to the source and minimizes a sphere that touches all the backprojection lines.",
        title="Gantry-isolated 3D isocenter diameter (mm)",
    )
    coll_2d_iso_diameter_mm: float = Field(
        description="The 2D isocenter diameter **of the collimator axis only** as determined by the collimator images in mm.",
        title="Collimator-isolated 2D isocenter diameter (mm)",
    )
    couch_2d_iso_diameter_mm: float = Field(
        description="The 2D isocenter diameter **of the couch axis only** as determined by the couch images in mm.",
        title="Couch-isolated 2D isocenter diameter (mm)",
    )
    gantry_coll_3d_iso_diameter_mm: float = Field(
        description="The 3D isocenter diameter **of the gantry and collimator axes** as determined by the gantry and collimator images in mm.",
        title="Gantry & Collimator combined 3D isocenter diameter (mm)",
    )
    num_total_images: int = Field(
        description="The total number of images analyzed.", title="Number of images"
    )
    num_gantry_images: int = Field(
        description="The number of images that were taken at different gantry angles and all other axes were at reference.",
        title="Number of gantry-axis images",
    )
    num_coll_images: int = Field(
        description="The number of images that were taken at different collimator angles and all other axes were at reference.",
        title="Number of collimator-axis images",
    )
    num_couch_images: int = Field(
        description="The number of images that were taken at different couch angles and all other axes were at reference.",
        title="Number of couch-axis images",
    )
    num_gantry_coll_images: int = Field(
        description="The number of images that were taken at different gantry and collimator angles and the couch was at reference.",
        title="Number of gantry+collimator axis images",
    )
    max_gantry_rms_deviation_mm: float = Field(
        description="The maximum RMS value of the field CAX to BB for the gantry axis images in mm. This is an alternative to the max/mean/median calculations."
    )
    max_epid_rms_deviation_mm: float = Field(
        description="The maximum RMS value of the field CAX to EPID center for the EPID images in mm. This is an alternative to the max/mean/median calculations."
    )
    max_coll_rms_deviation_mm: float = Field(
        description="The maximum RMS deviation of the field CAX to BB for the collimator axis images in mm. This is an alternative to the max/mean/median calculations."
    )
    max_couch_rms_deviation_mm: float = Field(
        description="The maximum RMS value of the field CAX to BB for the couch axis images in mm. This is an alternative to the max/mean/median calculations. This uses backprojection lines of the field center to the source and minimizes a sphere that touches all the backprojection lines."
    )
    bb_shift_vector: VectorSerialized = Field(
        description="The Cartesian vector that would move the BB to the radiation isocenter. Each value is in mm."
    )
    image_details: list[WinstonLutz2DResult] = Field(
        description="A list of the individual image results.",
    )
    keyed_image_details: dict[str, WinstonLutz2DResult] = Field(
        description="A **dictionary** of the individual image results. This is the same as ``image_details`` but keyed by the images using the axes values as the key. E.g. ``G0B45P0``. This can be used to identify individual images vs those in ``image_details``."
    )


class WinstonLutzMultiTargetMultiFieldResult(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    num_total_images: int = Field(
        description="The total number of images analyzed.",
        title="Number of images considered",
    )
    max_2d_field_to_bb_mm: float = Field(
        description="The maximum 2D distance from any BB to its field center.",
        title="Max field center -> BB distance (mm)",
    )
    mean_2d_field_to_bb_mm: float = Field(
        description="The mean 2D distance from any BB to its field center.",
        title="Mean field center -> BB distance (mm)",
    )
    median_2d_field_to_bb_mm: float = Field(
        description="The median 2D distance from any BB to its field center.",
        title="Median field center -> BB distance (mm)",
    )
    bb_arrangement: tuple[BBConfig, ...] = Field(
        description="A list of expected arrangements of the BBs"
    )
    bb_maxes: dict[str, float] = Field(
        description="A dictionary of the maximum 2D distances of each BB to its field center. The key is the BB name as defined in the arrangement."
    )
    bb_shift_vector: VectorSerialized = Field(
        description="The vector (in 3D cartesian space) to move the phantom to align with the isocenter in mm."
    )
    bb_shift_yaw: float = Field(
        description="The yaw rotation in degrees needed to align the phantom with the radiation isocenter."
    )
    bb_shift_pitch: float = Field(
        description="The pitch rotation needed in degrees to align the phantom with the radiation isocenter."
    )
    bb_shift_roll: float = Field(
        description="The roll rotation needed in degrees to align the phantom with the radiation isocenter."
    )


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


class WLBaseImage(image.LinacDicomImage):
    """Base class for a WL Image. Represents a single image with N fields and M BBs.

    Methods are provided to find the field CAXs and BBs and matching to the expected locations.
    """

    field_caxs: list[Point]
    bb_positions: list[Point]
    bb_arrangement: tuple[BBConfig]
    arrangement_matches: dict[
        str, BBFieldMatch
    ]  # a field CAX and BB matched to their respective nominal locations
    _gantry_reference: float
    _collimator_reference: float
    _couch_reference: float
    _snap_tolerance: float

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

    def analyze(
        self,
        bb_arrangement: tuple[BBConfig],
        is_open_field: bool = False,
        is_low_density: bool = False,
        shift_vector: Vector | None = None,
        snap_tolerance: float = 3,
        gantry_reference: float = 0,
        collimator_reference: float = 0,
        couch_reference: float = 0,
        bb_proximity_mm: float = 20,
        machine_scale: MachineScale = MachineScale.IEC61217,
    ) -> (tuple[Point], tuple[Point]):
        """Analyze the image for BBs and field CAXs.

        Parameters
        ----------
        bb_arrangement : tuple[BBConfig]
            The expected BB locations.
        is_open_field : bool
            Whether the field is open or not. If open, only one CAX is expected.
        is_low_density : bool
            Whether the BBs are low density (e.g. kV images).
        shift_vector : Vector, optional
            A vector to shift the detected BBs by. Useful for images that are not perfectly aligned.
        bb_proximity_mm : float
            The maximum distance a detected BB can be from the expected BB location in mm.

        See Also
        --------

        :meth:`~pylinac.winston_lutz.WinstonLutz.analyze`

        """
        if snap_tolerance < 0:
            raise ValueError("Snap tolerance must be >= 0")
        self._snap_tolerance = snap_tolerance
        self._gantry_reference = gantry_reference
        self._collimator_reference = collimator_reference
        self._couch_reference = couch_reference
        self.machine_scale = machine_scale
        self.check_inversion_by_histogram(percentiles=(0.01, 50, 99.99))
        self._clean_edges()
        self.ground()
        self.normalize()
        self.bb_arrangement = bb_arrangement
        field_caxs = self.find_field_centroids(is_open_field=is_open_field)
        field_matches = self.find_field_matches(
            field_caxs, bb_proximity_mm=bb_proximity_mm
        )
        detected_bb_points = self.find_bb_centroids(
            bb_diameter_mm=bb_arrangement[0].bb_size_mm,
            low_density=is_low_density,
        )
        if shift_vector:
            # apply shift to detected BB points
            lat, sup_inf = bb_projection_with_rotation(
                offset_left=-shift_vector.x,  # negative because left is negative x
                offset_up=shift_vector.z,
                offset_in=shift_vector.y,
                sad=self.sad,
                gantry=self.gantry_angle,
                couch=self.couch_angle,
                machine_scale=machine_scale,
            )
            # convert from mm to pixels and add to the detected points
            for p in detected_bb_points:
                p.x += lat * self.dpmm
                p.y -= (
                    sup_inf * self.dpmm
                )  # we subtract because the detected point is in image space, not coordinate space so we convert the shift from coordinate to image space
        bb_matches = self.find_bb_matches(
            detected_points=detected_bb_points, bb_proximity_mm=bb_proximity_mm
        )
        if len(bb_matches) != len(field_matches):
            raise ValueError("The number of detected fields and BBs do not match")
        if not field_matches:
            raise ValueError("No fields were detected")
        if not bb_matches:
            raise ValueError(BB_ERROR_MESSAGE)
        # we now have field CAXs and BBs matched to their respective nominal locations
        # merge the field and BBs per arrangement position
        combined_matches = {}
        for bb_name, bb_match in bb_matches.items():
            combined_matches[bb_name] = BBFieldMatch(
                epid=self.cax,
                field=field_matches[bb_name],
                bb=bb_match,
                dpmm=self.dpmm,
                gantry_angle=self.gantry_angle,
                couch_angle=self.couch_angle,
                sad=self.sad,
            )
        self._is_analyzed = True
        self.arrangement_matches = combined_matches

    def find_field_centroids(self, is_open_field: bool) -> list[Point]:
        """Find the field CAX(s) in the image. If the field is open or this is a vanilla WL, only one CAX is found."""
        if is_open_field:
            p = self.cax
        else:
            # TODO: Use metrics field finder
            # can't use it out of the box because the
            # analyze method doesn't pass the field size
            # using the global field locator without size will
            # find several other unrelated fields and show up in the image
            # the metric algorithm should have a post-processing function parameter
            min, max = np.percentile(self.array, [5, 99.9])
            threshold_img = self.as_binary((max - min) / 2 + min)
            filled_img = ndimage.binary_fill_holes(threshold_img)
            coords = ndimage.center_of_mass(filled_img)
            p = Point(x=coords[-1], y=coords[0])
        return [p]

    def find_field_matches(
        self, detected_points: list[Point], bb_proximity_mm: float
    ) -> dict[str, Point]:
        """Find matches between detected field points and the arrangement. See ``find_bb_matches`` for more info."""
        return self.find_bb_matches(detected_points, bb_proximity_mm=bb_proximity_mm)

    def find_bb_centroids(
        self, bb_diameter_mm: float, low_density: bool
    ) -> list[Point]:
        """Find BBs in the image. This method can return MORE than the desired number of BBs. Matching
        of the detected BBs vs the expected BBs is done in the ``find_bb_matches`` method.
        """
        bb_tolerance_mm = self._calculate_bb_tolerance(bb_diameter_mm)
        centers = self.compute(
            metrics=SizedDiskLocator.from_center_physical(
                expected_position_mm=(0, 0),
                search_window_mm=(40 + bb_diameter_mm, 40 + bb_diameter_mm),
                radius_mm=bb_diameter_mm / 2,
                radius_tolerance_mm=bb_tolerance_mm,
                invert=not low_density,
                detection_conditions=self.detection_conditions,
                name="BB",
            )
        )
        return centers

    def find_bb_matches(
        self, detected_points: list[Point], bb_proximity_mm: float
    ) -> dict[str, Point]:
        """Given an arrangement and detected BB positions, find the bbs that are closest to the expected positions.

        This is to prevent false positives from being detected as BBs (e.g. noise, couch, etc).
        The detected BBs are matched to the expected BBs based on proximity.

        These matches are linked to the individual BB arrangement by arrangement name.
        """
        bbs = {}
        for bb_arng in self.bb_arrangement:
            nominal_point = self.nominal_bb_position(bb_arng)
            distances = [
                nominal_point.distance_to(found_point)
                for found_point in detected_points
            ]
            min_distance = min(distances)
            min_distance_idx = distances.index(min_distance)
            if min_distance < bb_proximity_mm * self.dpmm:
                bbs[bb_arng.name] = detected_points[min_distance_idx]
        return bbs

    def field_to_bb_distances(self) -> list[float]:
        """The distances from the field CAXs to the BBs in mm. Useful for metrics as this is only
        the resulting floats vs a dict of points."""
        return [
            match.bb_field_distance_mm for match in self.arrangement_matches.values()
        ]

    def epid_to_bb_distances(self) -> list[float]:
        """The distances from the EPID center to the BBs in mm. Useful for metrics as this is only
        the resulting floats vs a dict of points."""
        return [
            match.bb_epid_distance_mm for match in self.arrangement_matches.values()
        ]

    def plotly(
        self,
        fig: go.Figure | None = None,
        show: bool = True,
        zoomed: bool = True,
        show_legend: bool = True,
        show_colorbar: bool = True,
    ) -> go.Figure:
        """
        Plot the image with the detected BB, outlines, and field CAX.

        Parameters
        ----------
        fig: go.Figure, None
            The Plotly figure
        show
            Whether to show the plot
        zoomed
            Whether to zoom in on the BBs. If False, no zooming is done and the entire image is shown.
        show_legend
            Whether to show the legend.
        show_colorbar
            Whether to show the colorbar.

        Returns
        -------
        go.Figure

        """
        if zoomed:
            min_x = int(
                round(
                    min([match.bb.x for match in self.arrangement_matches.values()])
                    - 20 * self.dpmm
                )
            )
            min_y = int(
                round(
                    min([match.bb.y for match in self.arrangement_matches.values()])
                    - 20 * self.dpmm
                )
            )
            max_x = int(
                round(
                    max([match.bb.x for match in self.arrangement_matches.values()])
                    + 20 * self.dpmm
                )
            )
            max_y = int(
                round(
                    max([match.bb.y for match in self.arrangement_matches.values()])
                    + 20 * self.dpmm
                )
            )
            x_indices = np.arange(min_x, max_x, 1)
            y_indices = np.arange(min_y, max_y, 1)
            window_array = self.array[int(min_y) : int(max_y), int(min_x) : int(max_x)]
        else:
            min_x, max_x = 0, self.shape[1]
            min_y, max_y = 0, self.shape[0]
            x_indices = np.arange(min_x, max_x, 1)
            y_indices = np.arange(min_y, max_y, 1)
            window_array = self.array
        fig = super().plotly(
            fig=fig,
            show=show,
            show_metrics=True,
            show_colorbar=show_colorbar,
            x=x_indices,
            y=y_indices,
            z=window_array,
        )
        # show EPID center; we use a custom line (vs `add_vertical_line`) because we might be zoomed and thus we need to shift
        # to a relative spot
        fig.add_scatter(
            x=[self.epid.x, self.epid.x],
            y=[y_indices[0], y_indices[-1]],
            line_color="blue",
            mode="lines",
            name="EPID Center (V)",
        )
        fig.add_scatter(
            x=[x_indices[0], x_indices[-1]],
            y=[self.epid.y, self.epid.y],
            line_color="blue",
            mode="lines",
            name="EPID Center (H)",
        )
        # show the field CAXs
        for match in self.arrangement_matches.values():
            fig.add_scatter(
                x=[match.field.x],
                y=[match.field.y],
                line_color="green",
                name="Field Center",
                mode="markers",
                marker_size=8,
                marker_symbol="square",
            )
            fig.add_scatter(
                x=[match.bb.x],
                y=[match.bb.y],
                line_color="cyan",
                name="Detected BB",
                mode="markers",
                marker_size=10,
                marker_symbol="circle",
            )
        fig.update_xaxes(range=[min_x, max_x])
        # bug in plotly; can't have autorange reversed and set this.
        fig.update_yaxes(range=[max_y, min_y], autorange=None)
        fig.update_layout(
            xaxis_title=f"Gantry={self.gantry_angle:.0f}, Coll={self.collimator_angle:.0f}, Couch={self.couch_angle:.0f}",
            yaxis_title=f"Max Nominal to BB: {max(self.field_to_bb_distances()):3.2f}mm",
        )
        fig.update_layout(
            showlegend=show_legend,
            title_text="\n".join(wrap(Path(self.path).name, 30)),
        )
        if show:
            fig.show()
        return fig

    def plot(
        self,
        ax: plt.Axes | None = None,
        show: bool = True,
        clear_fig: bool = False,
        zoom: bool = True,
        legend: bool = True,
        **kwargs,
    ) -> plt.Axes:
        """Plot an individual WL image.

        Parameters
        ----------
        ax : None, plt.Axes
            The axis to plot to. If None, a new figure is created.
        show : bool
            Whether to show the plot.
        clear_fig : bool
            Whether to clear the figure before plotting.
        zoom : bool
            Whether to zoom in on the BBs. If False, no zooming is done and the entire image is shown.
        legend : bool
            Whether to show the legend.
        """
        ax = super().plot(
            ax=ax, show=False, clear_fig=clear_fig, show_metrics=True, **kwargs
        )
        # show EPID center
        ax.axvline(x=self.epid.x, color="b")
        epid_handle = ax.axhline(y=self.epid.y, color="b")
        # show the field CAXs
        for match in self.arrangement_matches.values():
            (field_handle,) = ax.plot(match.field.x, match.field.y, "gs", ms=8)
            (bb_handle,) = ax.plot(match.bb.x, match.bb.y, "co", ms=10)
        if legend:
            ax.legend(
                (field_handle, bb_handle, epid_handle),
                ("Field CAX", "Detected BB", "EPID Center"),
                loc="upper right",
            )

        if zoom:
            # find the x and y limits based on the detected BB positions
            # and add a margin of 20mm
            min_x = (
                min([match.bb.x for match in self.arrangement_matches.values()])
                - 20 * self.dpmm
            )
            min_y = (
                min([match.bb.y for match in self.arrangement_matches.values()])
                - 20 * self.dpmm
            )
            max_x = (
                max([match.bb.x for match in self.arrangement_matches.values()])
                + 20 * self.dpmm
            )
            max_y = (
                max([match.bb.y for match in self.arrangement_matches.values()])
                + 20 * self.dpmm
            )
            ax.set_ylim([max_y, min_y])
            ax.set_xlim([min_x, max_x])
        # ax.set_yticklabels([])
        # ax.set_xticklabels([])
        ax.set_title("\n".join(wrap(Path(self.path).name, 30)), fontsize=10)
        ax.set_xlabel(
            f"G={self.gantry_angle:.0f}, B={self.collimator_angle:.0f}, P={self.couch_angle:.0f}"
        )
        ax.set_ylabel(f"Max Nominal to BB: {max(self.field_to_bb_distances()):3.2f}mm")
        if show:
            plt.show()
        return ax

    def nominal_bb_position(self, bb_config: BBConfig) -> Point:
        """Calculate the expected point position in 2D"""
        shift_x_mm, shift_y_mm = bb_projection_with_rotation(
            offset_left=bb_config.offset_left_mm,
            offset_up=bb_config.offset_up_mm,
            offset_in=bb_config.offset_in_mm,
            sad=self.sad,
            gantry=self.gantry_angle,
            couch=self.couch_angle,
            machine_scale=self.machine_scale,
        )
        # the field can be asymmetric, so use center of image
        expected_y = self.epid.y - shift_y_mm * self.dpmm
        expected_x = self.epid.x + shift_x_mm * self.dpmm
        return Point(x=expected_x, y=expected_y)

    @property
    def epid(self) -> Point:
        """Center of the EPID panel"""
        return self.cax

    def _calculate_bb_tolerance(self, bb_diameter: float) -> int:
        """Calculate the BB tolerance based on the BB diameter.
        Min will be 2 for 1.5mm and under. Will be 4 for diameters at or above 30mm."""
        y = (2, 4)
        x = (1.5, 30)
        return np.interp(bb_diameter, x, y)

    def to_axes(self) -> str:
        """Give just the axes values as a human-readable string"""
        return f"Gantry={self.gantry_angle:.1f}, Coll={self.collimator_angle:.1f}, Couch={self.couch_angle:.1f}"

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
        G0 = is_close_degrees(
            self.gantry_angle, self._gantry_reference, delta=self._snap_tolerance
        )
        B0 = is_close_degrees(
            self.collimator_angle,
            self._collimator_reference,
            delta=self._snap_tolerance,
        )
        P0 = is_close_degrees(
            self.couch_angle, self._couch_reference, delta=self._snap_tolerance
        )
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


@capture_warnings
class WinstonLutz2D(WLBaseImage, ResultsDataMixin[WinstonLutz2DResult]):
    """Holds individual Winston-Lutz EPID images, image properties, and automatically finds the field CAX and BB."""

    bb: Point
    field_cax: Point
    bb_arrangement: tuple[BBConfig]
    is_from_tiff: bool = False
    detection_conditions: list[callable] = [
        is_right_size_bb,
        is_round,
        is_right_circumference,
        is_symmetric,
        is_solid,
    ]

    def analyze(
        self,
        bb_size_mm: float = 5,
        low_density_bb: bool = False,
        open_field: bool = False,
        shift_vector: Vector | None = None,
        snap_tolerance: float = 3,
        gantry_reference: float = 0,
        collimator_reference: float = 0,
        couch_reference: float = 0,
        bb_proximity_mm: float = 20,
        machine_scale: MachineScale = MachineScale.IEC61217,
    ) -> None:
        """Analyze the image. See WinstonLutz.analyze for parameter details."""
        bb_config = BBArrangement.ISO
        bb_config[0].bb_size_mm = bb_size_mm
        super().analyze(
            bb_arrangement=bb_config,
            is_open_field=open_field,
            is_low_density=low_density_bb,
            shift_vector=shift_vector,
            snap_tolerance=snap_tolerance,
            gantry_reference=gantry_reference,
            collimator_reference=collimator_reference,
            couch_reference=couch_reference,
            bb_proximity_mm=bb_proximity_mm,
            machine_scale=machine_scale,
        )
        self.bb_arrangement = bb_config
        # these are set for the deprecated properties of the 2D analysis specifically where 1 field and 1 bb are expected.
        self.field_cax = self.arrangement_matches["Iso"].field
        self.bb = self.arrangement_matches["Iso"].bb

    def __repr__(self):
        return f"WLImage(gantry={self.gantry_angle:.1f}, coll={self.collimator_angle:.1f}, couch={self.couch_angle:.1f})"

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

    def save_plot(self, filename: str, **kwargs):
        """Save the image plot to file."""
        self.plot(show=False)
        plt.tight_layout()
        plt.savefig(filename, **kwargs)

    def _generate_results_data(self) -> WinstonLutz2DResult:
        """Present the results data and metadata as a dataclass or dict.
        The default return type is a dataclass."""
        if not self._is_analyzed:
            raise ValueError("The image is not analyzed. Use .analyze() first.")

        return WinstonLutz2DResult(
            variable_axis=self.variable_axis.value,
            cax2bb_vector=self.cax2bb_vector,
            cax2epid_vector=self.cax2epid_vector,
            cax2bb_distance=self.cax2bb_distance,
            cax2epid_distance=self.cax2epid_distance,
            bb_location=self.bb,
            field_cax=self.field_cax,
        )


@capture_warnings
class WinstonLutz(ResultsDataMixin[WinstonLutzResult], QuaacMixin):
    """Class for performing a Winston-Lutz test of the radiation isocenter."""

    images: list[WinstonLutz2D]  #:
    machine_scale: MachineScale  #:
    image_type = WinstonLutz2D
    bb: BB3D  # 3D representation of the BB; there is a .bb object for 2D images but is a 2D representation
    is_from_cbct: bool = False
    _bb_diameter: float
    _virtual_shift: str | None = None
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
        super().__init__()
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
    def from_zip(cls, zfile: str | BinaryIO | Path, **kwargs):
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
        apply_virtual_shift: bool = False,
        snap_tolerance: float = 3,
        gantry_reference: float = 0,
        collimator_reference: float = 0,
        couch_reference: float = 0,
        bb_proximity_mm: float = 20,
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
        apply_virtual_shift
            If True, applies a virtual shift to the BBs based on the shift necessary to place the BB at the radiation isocenter.
        snap_tolerance
            The tolerance of the axes values that will "snap" to the reference values. I.e. if the snap tolerance is 3 and the gantry is within 3 degrees of 0, it will snap to 0.
            This is helpful, e.g., when you've forgotten to reset the couch to 0 after a CBCT.
        gantry_reference
            The reference value for the gantry. This is when pylinac will consider the image to be a reference image. E.g.
            some customers take all images with collimator=45 and want that to be considered the reference. This is used in
            combination with the snap_tolerance. I.e. a gantry of 43 with snap tolerance of 3 and reference of 45 will snap to 45.
        collimator_reference
            The reference value for the collimator. See `gantry_reference`.
        couch_reference
            The reference value for the couch. See `gantry_reference`.
        bb_proximity_mm
            The maximum distance in mm that a detected BB can be from the expected BB position.
            For single-BB WL datasets, the expected BB position is isocenter.
        """
        self.machine_scale = machine_scale
        if self.is_from_cbct:
            low_density_bb = True
            open_field = True
        for img in self.images:
            img.analyze(
                bb_size_mm=bb_size_mm,
                low_density_bb=low_density_bb,
                open_field=open_field,
                snap_tolerance=snap_tolerance,
                gantry_reference=gantry_reference,
                collimator_reference=collimator_reference,
                couch_reference=couch_reference,
                bb_proximity_mm=bb_proximity_mm,
                machine_scale=machine_scale,
            )
        # we need to construct the BB representation to get the shift vector
        bb_config = BBArrangement.ISO[0]
        bb_config.bb_size_mm = bb_size_mm
        self.bb = BB3D(
            bb_config=bb_config,
            bb_matches=[img.arrangement_matches["Iso"] for img in self.images],
            scale=self.machine_scale,
        )
        if apply_virtual_shift:
            shift = self.bb_shift_vector
            self._virtual_shift = self.bb_shift_instructions()
            for img in self.images:
                img.analyze(
                    bb_size_mm=bb_size_mm,
                    low_density_bb=low_density_bb,
                    open_field=open_field,
                    shift_vector=shift,
                    snap_tolerance=snap_tolerance,
                    gantry_reference=gantry_reference,
                    collimator_reference=collimator_reference,
                    couch_reference=couch_reference,
                    machine_scale=machine_scale,
                )

        # in the vanilla WL case, the BB can only be represented by non-couch-kick images
        # the ray trace cannot handle the kick currently
        self.bb = BB3D(
            bb_config=bb_config,
            bb_matches=[img.arrangement_matches["Iso"] for img in self.images],
            scale=self.machine_scale,
        )
        self._is_analyzed = True
        self._bb_diameter = bb_size_mm

    @lru_cache()
    def _minimize_axis(self, axes: Axis | tuple[Axis, ...] = (Axis.GANTRY,)):
        """Return the minimization result of the given axis."""
        if isinstance(axes, Axis):
            axes = (axes,)

        things = [
            # we want the bb<->field because the BB is our reference point for vanilla WL
            # whether it's BB->field or field->BB is irrelevant for this case; we are not determining shift here.
            image.arrangement_matches["Iso"].bb_to_field_projection
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
            max_distance_to_lines,
            initial_guess,
            args=things,
            bounds=bounds,
            # eps default seems to cause test differences on different CPU platforms
            options={"eps": 1e-7},
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

        The shift is based on the paper by Low et al. See online documentation and the ``solve_3d_shift_vector_from_2d_planes`` function for more.,
        which is how the measured bb and field positions are determined.
        """
        # field minus BB will give the shift vector to RETURN TO ISO which is what we want. BB minus field would give the vector from field to the BB.
        return self.bb.measured_field_position - self.bb.measured_bb_position

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
            move += f"\nNew couch coordinates (cm): VRT: {new_vrt:3.2f}; LNG: {new_lng:3.2f}; LAT: {new_lat:3.2f}"
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
        distances = []
        for img in self.images:
            distances.extend(img.field_to_bb_distances())
        if metric == "max":
            return max(distances)
        elif metric == "median":
            return statistics.median(distances)
        elif metric == "mean":
            return statistics.mean(distances)

    @argue.options(metric=("max", "median", "mean"))
    def cax2epid_distance(self, metric: str = "max") -> float:
        """The distance in mm between the CAX and EPID center pixel for all images according to the given metric.

        Parameters
        ----------
        metric : {'max', 'median', 'mean'}
            The metric of distance to use.
        """
        distances = []
        for img in self.images:
            distances.extend(img.epid_to_bb_distances())
        if metric == "max":
            return max(distances)
        elif metric == "median":
            return statistics.median(distances)
        elif metric == "mean":
            return statistics.mean(distances)

    def plotly_analyzed_images(
        self,
        zoomed: bool = True,
        show_legend: bool = True,
        show: bool = True,
        show_colorbar: bool = True,
        **kwargs,
    ) -> dict[str, go.Figure]:
        """Plot the analyzed images in a Plotly figure.

        Parameters
        ----------
        zoomed : bool
            Whether to zoom in on the BBs of the 2D images.
        show_legend : bool
            Whether to show the legend on the plot.
        show : bool
            Whether to show the plot.
        show_colorbar : bool
            Whether to show the colorbar on the plot.
        kwargs
            Additional keyword arguments to pass to the plot.

        Returns
        -------
        dict
            A dictionary of the Plotly figures where the key is the name of the
            image and the value is the figure.
        """
        figs = {}
        for idx, wl_image in enumerate(self.images):
            fig = wl_image.plotly(
                show=False,
                show_legend=show_legend,
                zoomed=zoomed,
                show_colorbar=show_colorbar,
                **kwargs,
            )
            # we add a enumerator in case there are multiple images with the same axis values
            figs[f"{idx} - {wl_image.to_axes()}"] = fig

        # 3d iso visualization
        iso_fig = go.Figure()
        # origin lines
        limit = (
            max(
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
        for x, y, z in (
            ((-limit, limit), (0, 0), (0, 0)),
            ((0, 0), (-limit, limit), (0, 0)),
            ((0, 0), (0, 0), (-limit, limit)),
        ):
            iso_fig.add_scatter3d(
                mode="lines", x=x, y=y, z=z, name="Isocenter Axis", marker_color="blue"
            )
        # isosphere
        x, y, z = create_sphere_surface(
            radius=self.cax2bb_distance("max"),
            center=Point(),
        )
        iso_fig.add_surface(
            x=x,
            y=y,
            z=z,
            opacity=0.2,
            name="Isosphere",
            showscale=False,
            colorscale=[[0, "blue"], [1, "blue"]],
            showlegend=True,
        )
        # bb
        x, y, z = create_sphere_surface(
            radius=self._bb_diameter / 2,
            center=Point(
                self.bb.measured_bb_position.x,
                self.bb.measured_bb_position.y,
                self.bb.measured_bb_position.z,
            ),
        )
        iso_fig.add_surface(
            x=x,
            y=y,
            z=z,
            opacity=0.1,
            name="BB",
            showscale=False,
            colorscale=[[0, "red"], [1, "red"]],
            showlegend=True,
        )
        # coll iso size
        theta = np.linspace(0, 2 * np.pi, 100)
        circle_y = self.collimator_iso_size / 2 * np.cos(theta)  # Radius of the circle
        circle_z = self.collimator_iso_size / 2 * np.sin(theta)  # Radius of the circle
        circle_x = np.zeros_like(theta) + limit  # Fixed z-coordinate
        iso_fig.add_scatter3d(
            x=circle_x,
            y=circle_y,
            z=circle_z,
            mode="lines",
            line=dict(color="green", width=2),
            name="Collimator axis isosize projection",
            hovertext=f"Collimator isocenter size: {self.collimator_iso_size:.2f}mm",
            hoverinfo="text",
        )
        # gantry iso size
        circle_x = self.gantry_iso_size / 2 * np.cos(theta)  # Radius of the circle
        circle_z = self.gantry_iso_size / 2 * np.sin(theta)  # Radius of the circle
        circle_y = np.zeros_like(theta) - limit  # Fixed z-coordinate
        iso_fig.add_scatter3d(
            x=circle_x,
            y=circle_y,
            z=circle_z,
            mode="lines",
            line=dict(color="green", width=2),
            name="Gantry axis isosize projection",
            hoverinfo="text",
            hovertext=f"Gantry isocenter size: {self.gantry_iso_size:.2f}mm",
        )

        # couch isosize
        circle_x = self.couch_iso_size / 2 * np.cos(theta)  # Radius of the circle
        circle_y = self.couch_iso_size / 2 * np.sin(theta)  # Radius of the circle
        circle_z = np.zeros_like(theta) - limit  # Fixed z-coordinate
        iso_fig.add_scatter3d(
            x=circle_x,
            y=circle_y,
            z=circle_z,
            mode="lines",
            line=dict(color="green", width=2),
            name="Couch axis isosize projection",
            hoverinfo="text",
            hovertext=f"Couch isocenter size: {self.couch_iso_size:.2f}mm",
        )

        iso_fig.update_layout(
            scene=dict(
                xaxis_range=[-limit, limit],
                yaxis_range=[-limit, limit],
                zaxis_range=[-limit, limit],
                aspectmode="cube",
                xaxis_title="X (mm), Right (+)",
                yaxis_title="Y (mm), In (+)",
                zaxis_title="Z (mm), Up (+)",
            ),
            # set the camera so x axis is on the lower left; makes for more natural visualization
            scene_camera_eye=dict(x=-1, y=1, z=1),
            showlegend=show_legend,
        )
        add_title(iso_fig, "3D Isocenter visualization")
        figs["Isocenter Visualization"] = iso_fig

        # polar plot and POV plots
        for axis, start_angle, clock, marker_name in zip(
            (Axis.GANTRY, Axis.COLLIMATOR, Axis.COUCH, Axis.EPID),
            (90, 270, 270, 90),
            ("clockwise", "counterclockwise", "counterclockwise", "clockwise"),
            ("BB", "BB", "BB", "EPID"),
        ):
            if axis == Axis.EPID:
                attr = "cax2epid_vector"
                variable_axis = Axis.GANTRY
            else:
                attr = "cax2bb_vector"
                variable_axis = axis
            # get axis images, angles, and shifts
            imgs = [
                image
                for image in self.images
                if image.variable_axis in (variable_axis, Axis.REFERENCE)
            ]
            if not imgs:
                continue
            angles = [
                getattr(image, f"{variable_axis.value.lower()}_angle") for image in imgs
            ]
            xz_sag = np.array([getattr(img, attr).x for img in imgs])
            y_sag = np.array([getattr(img, attr).y for img in imgs])
            rms = np.sqrt(xz_sag**2 + y_sag**2)
            # append the first point to the end to close the loop
            angles = np.append(angles, angles[0])
            xz_sag = np.append(xz_sag, xz_sag[0])
            y_sag = np.append(y_sag, y_sag[0])
            rms = np.append(rms, rms[0])

            # X/Y POV plots
            fig = go.Figure()
            title = f"{axis.value} POV displacement"
            fig.add_scatter(
                x=xz_sag,
                y=y_sag,
                hovertext=[
                    f"Angle: {angle}\N{DEGREE SIGN}; Total: {r:.3f}mm"
                    for angle, r in zip(angles, rms)
                ],
                hoverinfo="text+x+y",
                mode="lines+markers",
                name=f"{marker_name} positions",
            )
            fig.add_scatter(
                x=[0],
                y=[0],
                name="Field Center",
                mode="markers",
            )
            fig.add_scatter(
                x=[xz_sag.mean()],
                y=[y_sag.mean()],
                hoverinfo="text+x+y",
                hovertext=f"Displacement: {math.hypot(xz_sag.mean(), y_sag.mean()):.3f}mm",
                name=f"{marker_name} Centroid",
                mode="markers",
            )
            add_title(fig, title)
            add_vertical_line(fig, 0, "black", name="y=0")
            add_horizontal_line(fig, 0, "black", name="x=0")
            fig.update_layout(
                showlegend=show_legend,
                xaxis_title="X (+Left) (mm)",
                yaxis_title="Y (+In) (mm)",
                xaxis_scaleanchor="y",
            )
            figs[title] = fig

            # polar plots
            fig = go.Figure()
            title = f"In-plane {axis.value} displacement"
            for name, data in zip(
                ["Y-axis (In/Out)", "X/Z-axis (Gantry plane)", "RMS"],
                [y_sag, xz_sag, rms],
            ):
                fig.add_scatterpolar(
                    r=data,
                    theta=angles,
                    name=name,
                )
                add_title(fig, f"{axis.value} Error Plot")
                fig.update_layout(
                    title=title,
                    showlegend=show_legend,
                    polar=dict(
                        angularaxis=dict(
                            rotation=start_angle,
                            direction=clock,  # Change the direction to clockwise (can also be 'counterclockwise')
                        ),
                        radialaxis=dict(
                            title="Displacement from Field (mm)", visible=True
                        ),
                    ),
                )
            # add 0-value highlight
            fig.add_scatterpolar(
                r=[0] * 100,
                theta=np.linspace(0, 360, 100),
                mode="lines",
                name="0-line",
            )
            figs[title] = fig

        if show:
            for f in figs.values():
                f.show()
        return figs

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
            marker = "co"
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
        # we can represent the iso sphere as a BB object; the nominal object isn't used, just the BB size
        # plot the field isocenter as x,y,z lines
        x_line = Line(
            Point(
                -limit,
                self.bb.measured_field_position.y,
                self.bb.measured_field_position.z,
            ),
            Point(
                limit, self.bb.measured_field_position.y, self.bb.measured_bb_position.z
            ),
        )
        x_line.plot2axes(ax, color="green", alpha=0.5)
        y_line = Line(
            Point(
                self.bb.measured_field_position.x,
                -limit,
                self.bb.measured_field_position.z,
            ),
            Point(
                self.bb.measured_field_position.x, limit, self.bb.measured_bb_position.z
            ),
        )
        y_line.plot2axes(ax, color="green", alpha=0.5)
        z_line = Line(
            Point(
                self.bb.measured_field_position.x,
                self.bb.measured_field_position.y,
                -limit,
            ),
            Point(
                self.bb.measured_field_position.x,
                self.bb.measured_field_position.y,
                limit,
            ),
        )
        z_line.plot2axes(ax, color="green", alpha=0.5, label="Field isocenter (x,y,z)")
        if plot_bb:
            self.bb.plot_measured(ax, color="cyan", alpha=0.6)
            # create an empty, fake line so we can add a label for the legend
            fake_line = Line(Point(0, 0, 0), Point(0, 0, 0))
            fake_line.plot2axes(ax, color="cyan", label=f"BB ({self._bb_diameter}mm)")
        if plot_isocenter_sphere:
            x, y, z = create_sphere_surface(
                radius=self.gantry_coll_iso_size / 2,
                center=Point(
                    self.bb.measured_bb_position.x,
                    self.bb.measured_bb_position.y,
                    self.bb.measured_bb_position.z,
                ),
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
                (self.bb.measured_field_position.x, self.bb.measured_field_position.y),
                radius=self.couch_iso_size / 2,
                fill=True,
                color="yellow",
                alpha=0.4,
                label=f"Couch-only iso ({self.couch_iso_size:3.2f}mm)",
            )
            ax.add_patch(circle)
            art3d.pathpatch_2d_to_3d(
                circle, z=self.bb.measured_field_position.z, zdir="z"
            )
        if plot_coll_iso:
            circle = plt.Circle(
                (self.bb.measured_field_position.y, self.bb.measured_field_position.z),
                radius=self.collimator_iso_size / 2,
                fill=True,
                color="blue",
                alpha=0.4,
                label=f"Collimator-only iso ({self.collimator_iso_size:3.2f}mm)",
            )
            ax.add_patch(circle)
            art3d.pathpatch_2d_to_3d(
                circle, z=self.bb.measured_field_position.x, zdir="x"
            )
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
        zoom: bool = True,
        legend: bool = True,
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
        zoom : bool
            Whether to zoom in around the BB.
        legend : bool
            Whether to show the legend.
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
                # plot the images and turn off extra axes
                if wl_image:
                    wl_image.plot(ax=mpl_axis, show=False, zoom=zoom, legend=legend)
                else:
                    mpl_axis.set_frame_on(False)
                    mpl_axis.axis("off")

            # set titles
            fig.suptitle(f"{axis.value} images", fontsize=14, y=1)
            fig.tight_layout()
            figs.append(fig)
            names.append("image")
        else:
            for wl_image in images:
                fig, axes = plt.subplots(**kwargs)
                wl_image.plot(ax=axes, show=False, zoom=zoom, legend=legend)
                # plot_image(wl_image, axes)
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
        ]
        if self._virtual_shift:
            result.append(
                f"Virtual shift applied to BB to place at isocenter: {self._virtual_shift}"
            )
        else:
            result.append(
                f"Shift to iso: facing gantry, move BB: {self.bb_shift_instructions()}"
            )
        result += [
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

    def _generate_results_data(self) -> WinstonLutzResult:
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

        individual_image_data = [i.results_data() for i in self.images]

        return WinstonLutzResult(
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
            bb_shift_vector=self.bb_shift_vector,
            image_details=individual_image_data,
            keyed_image_details=self._generate_keyed_images(individual_image_data),
        )

    def _generate_keyed_images(
        self, individual_image_data: list[WinstonLutz2DResult]
    ) -> dict[str, WinstonLutz2DResult]:
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

    def _quaac_datapoints(self) -> dict[str, QuaacDatum]:
        if not self._is_analyzed:
            raise ValueError("The set is not analyzed. Use .analyze() first.")
        result_data = self.results_data()
        dataset = {
            "Max 2D CAX->BB": QuaacDatum(
                value=result_data.max_2d_cax_to_bb_mm,
                unit="mm",
                description="The maximum 2D distance of any image from the CAX to the BB.",
            ),
            "Median 2D CAX->BB": QuaacDatum(
                value=result_data.median_2d_cax_to_bb_mm,
                unit="mm",
                description="The median 2D distance of any image from the CAX to the BB.",
            ),
            "Max 2D CAX->EPID": QuaacDatum(
                value=result_data.max_2d_cax_to_epid_mm,
                unit="mm",
                description="The maximum 2D distance of any image from the CAX to the EPID.",
            ),
            "Median 2D CAX->EPID": QuaacDatum(
                value=result_data.median_2d_cax_to_epid_mm,
                unit="mm",
                description="The median 2D distance of any image from the CAX to the EPID.",
            ),
            "Gantry-only 3D Isocenter Diameter": QuaacDatum(
                value=result_data.gantry_3d_iso_diameter_mm,
                unit="mm",
                description="The diameter of the 3D isocenter sphere when considering the gantry-only images.",
            ),
            "Gantry+Collimator 3D Isocenter Diameter": QuaacDatum(
                value=result_data.gantry_coll_3d_iso_diameter_mm,
                unit="mm",
                description="The diameter of the 3D isocenter sphere when considering the gantry and collimator images.",
            ),
            "Collimator 2D Isocenter Diameter": QuaacDatum(
                value=result_data.coll_2d_iso_diameter_mm,
                unit="mm",
                description="The diameter of the 2D isocenter circle when considering the collimator images.",
            ),
            "Couch 2D Isocenter Diameter": QuaacDatum(
                value=result_data.couch_2d_iso_diameter_mm,
                unit="mm",
                description="The diameter of the 2D isocenter circle when considering the couch images.",
            ),
        }
        return dataset


class WinstonLutzMultiTargetMultiFieldImage(WLBaseImage):
    """A 2D image of a WL delivery, but where multiple BBs are in use."""

    detection_conditions = [is_round, is_symmetric, is_modest_size]
    field_conditions = [is_square, is_right_square_size]

    def find_field_centroids(self, is_open_field: bool) -> list[Point]:
        """Find the centroid of the radiation field based on a 50% height threshold.
        This applies the field detection conditions and also a nearness condition.

        Returns
        -------
        points
            The CAX point locations.
        """
        if is_open_field:
            return [self.cax]

        # find all the fields by setting the field to the mean rad size and tolerance
        # to max-min field sizes across the arrangements
        max_field_size = max(
            self.bb_arrangement, key=lambda x: x.rad_size_mm
        ).rad_size_mm
        min_field_size = min(
            self.bb_arrangement, key=lambda x: x.rad_size_mm
        ).rad_size_mm
        mean_field_size = (max_field_size + min_field_size) / 2
        tolerance_field_size = max(
            (max_field_size - min_field_size) * 1.2, 0.1 * mean_field_size
        )
        points = self.compute(
            metrics=GlobalSizedFieldLocator.from_physical(
                max_number=len(self.bb_arrangement),
                field_height_mm=mean_field_size,
                field_width_mm=mean_field_size,
                field_tolerance_mm=tolerance_field_size,
            )
        )
        return points

    def find_bb_centroids(
        self, bb_diameter_mm: float, low_density: bool
    ) -> list[Point]:
        """Find the specific BB based on the arrangement rather than a single one. This is in local pixel coordinates"""
        centers = []
        for bb in self.bb_arrangement:
            bb_diameter_mm = bb.bb_size_mm
            bb_tolerance_mm = self._calculate_bb_tolerance(bb_diameter_mm)
            left, sup = bb_projection_with_rotation(
                offset_left=bb.offset_left_mm,
                offset_up=bb.offset_up_mm,
                offset_in=bb.offset_in_mm,
                gantry=self.gantry_angle,
                couch=self.couch_angle,
                sad=self.sad,
            )
            try:
                new_centers = self.compute(
                    metrics=SizedDiskLocator.from_center_physical(
                        expected_position_mm=Point(
                            x=left, y=-sup
                        ),  # we do -sup because sup is in WL coordinate space but the disk locator is in image coordinates.
                        search_window_mm=(40 + bb_diameter_mm, 40 + bb_diameter_mm),
                        radius_mm=bb_diameter_mm / 2,
                        radius_tolerance_mm=bb_tolerance_mm / 2,
                        invert=not low_density,
                        detection_conditions=self.detection_conditions,
                    )
                )
                centers.extend(new_centers)
            except ValueError:
                pass
        return centers


@capture_warnings
class WinstonLutzMultiTargetMultiField(WinstonLutz):
    machine_scale: MachineScale  #:
    images: Sequence[WinstonLutzMultiTargetMultiFieldImage]  #:
    image_type = WinstonLutzMultiTargetMultiFieldImage
    bb_arrangement: tuple[BBConfig]  #:
    bbs: list[BB3D]  #:  3D representation of the BBs

    @classmethod
    def from_demo_images(cls):
        """Instantiate using the demo images."""
        demo_file = retrieve_demo_file(name="SNC_MTWL_demo.zip")
        return cls.from_zip(demo_file)

    @staticmethod
    def run_demo():
        """Run the Winston-Lutz MT MF demo, which loads the demo files, prints results, and plots a summary image."""
        wl = WinstonLutzMultiTargetMultiField.from_demo_images()
        wl.analyze(bb_arrangement=BBArrangement.DEMO)
        print(wl.results())
        wl.plot_images()

    def analyze(
        self,
        bb_arrangement: tuple[BBConfig, ...],
        is_open_field: bool = False,
        is_low_density: bool = False,
        machine_scale: MachineScale = MachineScale.IEC61217,
        bb_proximity_mm: float = 10,
    ):
        """Analyze the WL images.

        Parameters
        ----------
        bb_arrangement
            The arrangement of the BBs in the phantom. A dict with offset and BB size keys. See the ``BBArrangement`` class for
            keys and syntax.

        See Also
        --------
        WinstonLutz.analyze for other parameter info.
        """
        self.machine_scale = machine_scale
        self.bb_arrangement = bb_arrangement
        for img in self.images:
            img.analyze(
                bb_arrangement=bb_arrangement,
                is_open_field=is_open_field,
                is_low_density=is_low_density,
                bb_proximity_mm=bb_proximity_mm,
                machine_scale=machine_scale,
            )

        self.bbs = []
        for arrangement in self.bb_arrangement:
            # add bbs to the matches if the match is in the given image
            matches = []
            for img in self.images:
                if arrangement.name in img.arrangement_matches:
                    match = img.arrangement_matches[arrangement.name]
                    matches.append(match)
            # ray lines are used for plotting
            bb = BB3D(
                bb_config=arrangement,
                bb_matches=matches,
                scale=self.machine_scale,
            )
            self.bbs.append(bb)
        self._is_analyzed = True

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
        x_lim = max(
            max([np.abs(bb.measured_bb_position.x) for bb in self.bbs]) * 1.3, 10
        )
        y_lim = max(
            max([np.abs(bb.measured_bb_position.y) for bb in self.bbs]) * 1.3, 10
        )
        z_lim = max(
            max([np.abs(bb.measured_bb_position.z) for bb in self.bbs]) * 1.3, 10
        )
        limit = viewbox_mm or max(x_lim, y_lim, z_lim)
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        _, relevant_images = self._get_images(
            axis=(Axis.REFERENCE, Axis.GB_COMBO, Axis.COLLIMATOR, Axis.GANTRY)
        )
        # we can represent the iso sphere as a BB object; the nominal object isn't used, just the BB size
        # the ray lines are what we want to plot as a sphere
        # plot the x,y,z origin lines
        x_line = Line(Point(-100, 0, 0), Point(100, 0, 0))
        x_line.plot2axes(ax, color="green", alpha=0.5)
        y_line = Line(Point(0, -100, 0), Point(0, 100, 0))
        y_line.plot2axes(ax, color="green", alpha=0.5)
        z_line = Line(Point(0, 0, -100), Point(0, 0, 100))
        z_line.plot2axes(
            ax, color="green", alpha=0.5, label="Nominal isocenter (x,y,z)"
        )
        if plot_bb:
            for bb in self.bbs:
                bb.plot_measured(ax, color="cyan", alpha=0.6)
                bb.plot_nominal(ax, color="green", alpha=0.6)

            # create an empty, fake line so we can add a label for the legend
            fake_line = Line(Point(0, 0, 0), Point(0, 0, 0))
            fake_line.plot2axes(ax, color="cyan", label="Measured BB")
            fake_line = Line(Point(0, 0, 0), Point(0, 0, 0))
            fake_line.plot2axes(ax, color="green", label="Nominal BB")

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
        return fig, ax

    @property
    def bb_shift_vector(self) -> (Vector, float, float, float):
        """Calculate the ideal shift in 6 degrees of freedom to place the BB at the isocenter.

        Returns
        -------
        Vector
            The ideal shift vector in mm for the cartesian coordinates. X,Y, and Z follow the pylinac coordinate convention.
        float
            Yaw; The ideal rotation about the Z-axis in degrees.
        float
            Pitch; The ideal rotation about the X-axis in degrees.
        float
            Roll; The ideal rotation about the Y-axis in degrees.

        See Also
        --------
        Euler Angles: https://en.wikipedia.org/wiki/Euler_angles
        Gimbal Lock: https://en.wikipedia.org/wiki/Gimbal_lock
        """
        return align_points(
            measured_points=[bb.measured_bb_position for bb in self.bbs],
            ideal_points=[bb.measured_field_position for bb in self.bbs],
        )

    def bb_shift_instructions(self) -> str:
        """Return a string that provides instructions on how to shift the BB to the isocenter."""
        translation, yaw, pitch, roll = self.bb_shift_vector
        x_dir = "LEFT" if translation.x < 0 else "RIGHT"
        y_dir = "IN" if translation.y > 0 else "OUT"
        z_dir = "UP" if translation.z > 0 else "DOWN"
        move = f"{x_dir} {abs(translation.x):2.2f}mm; {y_dir} {abs(translation.y):2.2f}mm; {z_dir} {abs(translation.z):2.2f}mm; Rotation {yaw:2.2f}°; Pitch {pitch:2.2f}°; Roll {roll:2.2f}°"
        return move

    def _quaac_datapoints(self) -> dict[str, QuaacDatum]:
        """Generate the Quaac datapoints for MTMF Winston-Lutz analysis"""
        if not self._is_analyzed:
            raise ValueError("The set is not analyzed. Use .analyze() first.")
        result_data = self.results_data()
        dataset = {
            "Max 2D CAX->BB": QuaacDatum(
                value=result_data.max_2d_field_to_bb_mm,
                unit="mm",
                description="The maximum 2D distance of any image from the CAX to the BB.",
            ),
            "Median 2D CAX->BB": QuaacDatum(
                value=result_data.median_2d_field_to_bb_mm,
                unit="mm",
                description="The median 2D distance of any image from the CAX to the BB.",
            ),
            "Mean 2D CAX->BB": QuaacDatum(
                value=result_data.mean_2d_field_to_bb_mm,
                unit="mm",
                description="The mean 2D distance of any image from the CAX to the BB.",
            ),
            "BB Shift (Yaw)": QuaacDatum(
                value=result_data.bb_shift_yaw,
                unit="degrees",
                description="The ideal yaw rotation to place the BB at the isocenter.",
            ),
            "BB Shift (Pitch)": QuaacDatum(
                value=result_data.bb_shift_pitch,
                unit="degrees",
                description="The ideal pitch rotation to place the BB at the isocenter.",
            ),
            "BB Shift (Roll)": QuaacDatum(
                value=result_data.bb_shift_roll,
                unit="degrees",
                description="The ideal roll rotation to place the BB at the isocenter.",
            ),
        }
        return dataset

    def _couch_rotation_error(self) -> dict[str, dict[str, float]]:
        """Calculate the couch rotation error in degrees for reference and couch-kicked images.
        This just for feature parity with SNC 🤦; the BB shift vector is more important.

        Returns
        -------
        dict
            A dictionary where the keys are the image paths and the values are a dictionary with the yaw error and nominal couch angle.
        """
        couch_results = {}
        couch_images = [
            img
            for img in self.images
            if img.variable_axis in (Axis.COUCH, Axis.REFERENCE)
        ]
        for img in couch_images:
            measured_points = [m.bb for m in img.arrangement_matches.values()]
            ideal_points = [m.field for m in img.arrangement_matches.values()]
            _, yaw, _, _ = align_points(measured_points, ideal_points)
            couch_results[img.base_path] = {
                "yaw error": yaw,
                "couch angle": img.couch_angle,
            }
        return couch_results

    @property
    def gantry_coll_iso_size(self) -> float:
        raise NotImplementedError("Not yet implemented")

    @property
    def collimator_iso_size(self) -> float:
        raise NotImplementedError("Not yet implemented")

    @property
    def couch_iso_size(self) -> float:
        raise NotImplementedError("Not yet implemented")

    @property
    def gantry_iso_size(self) -> float:
        raise NotImplementedError("Not yet implemented")

    def plotly_analyzed_images(
        self,
        zoomed: bool = True,
        show_legend: bool = True,
        show: bool = True,
        show_colorbar: bool = True,
        **kwargs,
    ) -> dict[str, go.Figure]:
        """Plot the analyzed set of images to Plotly figures.


        Parameters
        ----------
        zoomed : bool
            Whether to zoom in on the 2D image plots.
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
        for idx, wl_image in enumerate(self.images):
            fig = wl_image.plotly(
                show=False,
                show_legend=show_legend,
                zoomed=zoomed,
                show_colorbar=show_colorbar,
                **kwargs,
            )
            # we add a enumerator in case there are multiple images with the same axis values
            figs[f"{idx} - {wl_image.to_axes()}"] = fig

        x_lim = max(
            max([np.abs(bb.nominal_bb_position.x) for bb in self.bbs]) * 1.3, 10
        )
        y_lim = max(
            max([np.abs(bb.nominal_bb_position.y) for bb in self.bbs]) * 1.3, 10
        )
        z_lim = max(
            max([np.abs(bb.nominal_bb_position.z) for bb in self.bbs]) * 1.3, 10
        )
        limit = max(x_lim, y_lim, z_lim)
        # 3d iso visualization
        iso_fig = go.Figure()
        figs["Isocenter Visualization"] = iso_fig
        for x, y, z in (
            ((-limit, limit), (0, 0), (0, 0)),
            ((0, 0), (-limit, limit), (0, 0)),
            ((0, 0), (0, 0), (-limit, limit)),
        ):
            iso_fig.add_scatter3d(
                mode="lines", x=x, y=y, z=z, name="Isocenter Axis", marker_color="blue"
            )
        # bbs
        for bb in self.bbs:
            bb.plotly_measured(iso_fig, color="cyan", opacity=0.6)
            bb.plotly_nominal(iso_fig, color="green", opacity=0.6)

        # isosphere
        x, y, z = create_sphere_surface(
            radius=self.cax2bb_distance("max"),
            center=Point(),
        )
        iso_fig.add_surface(
            x=x,
            y=y,
            z=z,
            opacity=0.2,
            name="Isosphere",
            showscale=False,
            colorscale=[[0, "blue"], [1, "blue"]],
            showlegend=True,
        )
        iso_fig.update_layout(
            scene=dict(
                xaxis_range=[-limit, limit],
                yaxis_range=[-limit, limit],
                zaxis_range=[-limit, limit],
                aspectmode="cube",
                xaxis_title="X (mm), Right (+)",
                yaxis_title="Y (mm), In (+)",
                zaxis_title="Z (mm), Up (+)",
            ),
            # set the camera so x axis is on the lower left; makes for more natural visualization
            scene_camera_eye=dict(x=-1, y=1, z=1),
            showlegend=show_legend,
        )
        add_title(iso_fig, "3D Isocenter visualization")
        if show:
            for f in figs.values():
                f.show()
        return figs

    def plot_images(
        self, show: bool = True, zoomed: bool = True, legend: bool = True, **kwargs
    ) -> (list[plt.Figure], list[str]):
        """Make a plot for each BB. Each plot contains the analysis of that BB on each image
        it was found."""
        figs, names = [], []
        figsize = kwargs.pop("figsize", None) or (8, 8)
        for img in self.images:
            fig, axes = plt.subplots(figsize=figsize, **kwargs)
            img.plot(ax=axes, show=False, zoom=zoomed, legend=legend)
            fig.tight_layout()
            figs.append(fig)
            names.append(img.base_path)
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

    def _generate_results_data(self) -> WinstonLutzMultiTargetMultiFieldResult:
        """Present the results data and metadata as a dataclass or dict.
        The default return type is a dataclass."""
        if not self._is_analyzed:
            raise ValueError("The set is not analyzed. Use .analyze() first.")

        # for backward-compatibility, we have to find the max distance for each BB
        # across the images.
        bb_maxes = {}
        for bb in self.bb_arrangement:
            max_d = 0.0
            for img in self.images:
                if bb.name in img.arrangement_matches:
                    max_d = max(
                        max_d, img.arrangement_matches[bb.name].bb_field_distance_mm
                    )
            bb_maxes[bb.name] = max_d

        translation, yaw, pitch, roll = self.bb_shift_vector
        return WinstonLutzMultiTargetMultiFieldResult(
            num_total_images=len(self.images),
            max_2d_field_to_bb_mm=self.max_bb_deviation_2d,
            mean_2d_field_to_bb_mm=self.mean_bb_deviation_2d,
            median_2d_field_to_bb_mm=self.median_bb_deviation_2d,
            bb_maxes=bb_maxes,
            bb_arrangement=self.bb_arrangement,
            bb_shift_vector=translation,
            bb_shift_yaw=yaw,
            bb_shift_pitch=pitch,
            bb_shift_roll=roll,
        )

    def plot_summary(self, show: bool = True, fig_size: tuple | None = None):
        raise NotImplementedError("Not yet implemented")

    def plot_axis_images(
        self, axis: Axis = Axis.GANTRY, show: bool = True, ax: plt.Axes | None = None
    ):
        raise NotImplementedError("Not yet implemented")

    @property
    def max_bb_deviation_2d(self) -> float:
        """The maximum distance from any measured BB to its nominal position"""
        return self.cax2bb_distance(metric="max")

    @property
    def mean_bb_deviation_2d(self) -> float:
        """The mean distance from any measured BB to its nominal position"""
        return self.cax2bb_distance(metric="mean")

    @property
    def median_bb_deviation_2d(self) -> float:
        """The median distance from any measured BB to its nominal position"""
        return self.cax2bb_distance(metric="median")

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
            f"Max 2D distance of any BB->Field: {self.max_bb_deviation_2d:.2f} mm",
            f"Mean 2D distance of any BB->Field: {self.mean_bb_deviation_2d:.2f} mm",
            f"Median 2D distance of any BB->Field: {self.median_bb_deviation_2d:.2f} mm",
            "",
        ]
        bb_descriptions = [[bb.name, bb.to_human()] for bb in self.bb_arrangement]
        bb_names = [bb[0] for bb in bb_descriptions]
        result += tabulate(bb_descriptions, headers=["BB #", "Description"]).split("\n")
        result += [
            "",
        ]

        data = []
        for img in self.images:
            img_name = img.base_path[-20:]
            gantry = f"{img.gantry_angle:.1f}"
            collimator = f"{img.collimator_angle:.1f}"
            couch = f"{img.couch_angle:.1f}"
            deviations = []
            # loop through the expected BBs.
            # not every BB may be in every image, so we have to loop through the
            # BBs and find the match if there is one.
            for bb in self.bb_arrangement:
                match = Enumerable(img.arrangement_matches.items()).single_or_default(
                    lambda x: x[0] == bb.name
                )
                if match:
                    deviations.append(f"{match[1].bb_field_distance_mm:.2f}")
                else:
                    deviations.append("---")
            data.append([img_name, gantry, collimator, couch, *deviations])
        result += tabulate(data, headers=["Image", "G", "C", "P", *bb_names]).split(
            "\n"
        )
        result += [""]

        # calculate couch-kick errors
        couch_results = self._couch_rotation_error()
        couch_data = [
            [name[-20:], v["couch angle"], f"{v['yaw error']:.2f}"]
            for name, v in couch_results.items()
        ]
        result += tabulate(
            couch_data, headers=["Image", "Couch Angle", "Yaw Error (\N{DEGREE SIGN})"]
        ).split("\n")
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


# class WinstonLutzMultiTargetSingleFieldImage(WinstonLutzMultiTargetMultiFieldImage):
#     def find_field_centroids(self, is_open_field: bool) -> list[Point]:
#         """For the single field case, the field centroid is the same as the CAX."""
#         # TODO: use the global field finder
#         return [self.center]
#
#     def find_field_matches(self, detected_points: list[Point]) -> dict[str, Point]:
#         """For the single field case, the field centroid is the same as the CAX for every BB."""
#         return {bb.name: detected_points[0] for bb in self.bb_arrangement}
#
#     def field_to_bb_distances(self) -> list[float]:
#         """For a single-field, multi-BB setup we shift the BBs to the isocenter by shifting by the nominal offset.
#         The delta between the nominal and actual BB position is the same distance as if we were at iso.
#         """
#         distances = []
#         for match in self.arrangement_matches.values():
#             distances.append(match["bb nominal"].distance_to(match["bb"]) / self.dpmm)
#         return distances
#
#
# class WinstonLutzMultiTargetSingleField(WinstonLutzMultiTargetMultiField):
#     image_type = WinstonLutzMultiTargetSingleFieldImage
#     images: Sequence[WinstonLutzMultiTargetSingleFieldImage]


def max_distance_to_lines(p, lines: Iterable[Line]) -> float:
    """Calculate the maximum distance to any line from the given point."""
    point = Point(p[0], p[1], p[2])
    return max(line.distance_to(point) for line in lines)


def bb_projection_with_rotation(
    offset_left: float,
    offset_up: float,
    offset_in: float,
    gantry: float,
    couch: float,
    sad: float = 1000,
    machine_scale: MachineScale = MachineScale.IEC61217,
) -> (float, float):
    """Calculate the isoplane projection onto the panel at the given SSD.

    This function applies a rotation around the gantry plane (X/Z) to the
    ball bearing (BB) position and calculates its projection onto the isocenter plane in the beam's eye view.

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
    left_right_projection (float): The projection of the BB onto the panel in the left/right direction.
        Left is negative, right is positive. This is always in the plane normal to the CAX.
    superior_inferior_projection (float): The projection of the BB onto the panel in the superior/inferior direction.
        Superior is positive, inferior is negative. This is always in the plane normal to the CAX which
        in this case is also the absolute Y coordinate system axis.
    """
    # Define the BB positions in the patient coordinate system (ap, lr, si)
    bb_positions = np.array([offset_up, offset_left, offset_in])

    gantry_rot, _, couch_rot = convert(
        input_scale=machine_scale,
        output_scale=MachineScale.IEC61217,
        gantry=gantry,
        collimator=0,
        rotation=couch,
    )
    # Apply the rotation matrix to the BB positions
    collimator = 0  # Collimator doesn't change positional projection onto panel
    rotation_matrix = Rotation.from_euler(
        "xyz",
        [-couch_rot, collimator, gantry_rot],
        degrees=True,  # negative couch due to origin shift vs coordinate space
    )
    rotated_positions = rotation_matrix.apply(bb_positions)

    # Calculate the projection onto the panel at the given SAD
    bb_magnification = sad / (
        sad - rotated_positions[0]
    )  # Distance from source to panel
    imager_projection = (
        np.array([rotated_positions[1], rotated_positions[2]]) * bb_magnification
    )
    return -imager_projection[0], imager_projection[1]


def straight_ray(vector: Vector, gantry_angle: float) -> Line:
    """A ray projecting from the gantry source point in the given direction.

    Notes
    -----
    This does NOT account for couch rotation. It also generates a straight line
    that goes through the vector. I.e. if the vector is 5mm left, the line will be straight
    up and down 5mm left. It does not account for 2D projection from the gantry source.
    The reason is that in vanilla WL, the assumption is that the field is about isocenter (i.e. not purposely offset)
    and thus creating straight lines is a good approximation.
    In the Multi-BB case, our fields are purposely off-center,
    so this assumption does not hold for that case. See the ``ray`` function for that.

    TL;DR: this is useful for vanilla WL but not Multi-BB.
    """
    p1 = Point()
    p2 = Point()
    # point 1 - ray origin
    p1.x = vector.x * cos(gantry_angle) + 20 * sin(gantry_angle)
    p1.z = vector.x * -sin(gantry_angle) + 20 * cos(gantry_angle)
    p1.y = vector.y
    # point 2 - ray destination
    p2.x = vector.x * cos(gantry_angle) - 20 * sin(gantry_angle)
    p2.z = vector.x * -sin(gantry_angle) - 20 * cos(gantry_angle)
    p2.y = vector.y
    line = Line(p1, p2)
    return line


def solve_3d_shift_vector_from_2d_planes(
    xs: Sequence[float],
    ys: Sequence[float],
    thetas: Sequence[float],
    phis: Sequence[float],
    scale: MachineScale,
) -> Vector:
    """Solve for long, lat, and vert given arrays of x's, y's, theta's, and phi's.
    This is Y, X, and Z in pylinac coordinate conventions.

    This is a generalization of the Low et al. equation 6 and 7. This is used for both vanilla WL and Multi-BB WL.
    This is intended to be a relatively strict interpretation of the Low et al. equations.
    E.g. the gantry and couch angles are converted to Varian Standard scale.
    The idea is this is meant for high readability when comparing to the Low et al. equations.
    However, this does mean we have to invert the LONG/Y axis. See the conversion table in the docs

    Parameters
    ----------
    xs : array_like
        Array of x coordinates. Positive is to the right in the local plane.
    ys : array_like
        Array of y coordinates. Positive is up in the local plane.
    thetas : array_like
        Array of gantry angles.
    phis : array_like
        Array of couch angles.
    scale: MachineScale
        The scale of the machine. IEC61217 is the default. This converts the gantry and couch angles to Varian Standard scale
        which is to match Low's equations.

    Returns
    -------
    Vector
        The coordinates are in the pylinac coordinate system. See docs.
    """
    if not (len(xs) == len(ys) == len(thetas) == len(phis)):
        raise ValueError("The x, y, theta, and phi arrays must all be the same length.")
    n = len(xs)

    # convert the angles to Varian Standard scale
    f_thetas, f_phis = [], []
    for theta, phi in zip(thetas, phis):
        g, _, c = convert(
            scale,
            MachineScale.VARIAN_STANDARD,
            gantry=theta,
            collimator=0,
            rotation=phi,
        )
        f_thetas.append(g)
        f_phis.append(c)

    # Initialize A matrix and xi matrices
    A = np.zeros((2 * n, 3))
    xi = np.zeros(2 * n)

    for i in range(n):
        # A general rotation matrix and also Low eqn 6a
        A[2 * i, :] = [-cos(f_phis[i]), -sin(f_phis[i]), 0]
        A[2 * i + 1, :] = [
            -cos(f_thetas[i]) * sin(f_phis[i]),
            cos(f_thetas[i]) * cos(f_phis[i]),
            -sin(f_thetas[i]),
        ]
        # Low eqn 7
        xi[2 * i] = ys[
            i
        ]  # usually this would be (x, y) but Figure 1 of Low is rotated 90 degrees.
        xi[2 * i + 1] = -xs[i]

    # equation 9a; B is the pseudo-inverse of A
    B = np.linalg.pinv(A)
    # full equation 9
    long, lat, vert = B.dot(xi).squeeze()

    # At this point we have the **Low et al shift vector**. Conversion to pylinac coordinates
    # requires flipping the LONG axis; see the conversion table in the docs.
    # we use the X/Y/Z terms for the shift to be consistent with our own convention
    # and for clarity of the conversion
    # Finally, this is the **SHIFT VECTOR**. The 3D position in space is the inverse of this
    y = -long
    x = lat
    z = vert
    return Vector(x=x, y=y, z=z)


def solve_3d_position_from_2d_planes(
    xs: Sequence[float],
    ys: Sequence[float],
    thetas: Sequence[float],
    phis: Sequence[float],
    scale: MachineScale,
) -> Vector:
    """Solve for the 3D position of a BB in space given 2D images/vectors and the gantry/couch angles.

    The good news is that the position in space is the inverse of the shift vector!
    """
    return -solve_3d_shift_vector_from_2d_planes(xs, ys, thetas, phis, scale)


def conventional_to_euler_notation(axes_resolution: str) -> str:
    """Convert conventional understandings of 6DOF rotations into Euler notation.

    Ensures we don't mix up x, y, z with the pylinac coordinate system.
    """
    EULER = {
        # FROM THE COUCH PERSPECTIVE
        "pitch": "x",  # positive pitch goes up
        "yaw": "z",
        "roll": "y",  # positive angle rolls to the right
    }
    axes = axes_resolution.split(",")
    euler = "".join([EULER[a.strip()] for a in axes])
    return euler


def align_points(
    measured_points: Sequence[Point],
    ideal_points: Sequence[Point],
    axes_order: str = "roll,pitch,yaw",
) -> (Vector, float, float, float):
    """
    Aligns a set of measured points to a set of ideal points in 3D space, returning the
    translation and yaw rotation needed.

    Parameters
    ----------
    measured_points : np.ndarray
        The measured points as an Nx3 numpy array.
    ideal_points : np.ndarray
        The ideal points as an Nx3 numpy array.
    axes_order : str
        The order in which to resolve the axes.
        Resolution is **not** independent of the axes order. I.e. doing 'yaw,pitch,roll' may result
        in a different outcome.

    Returns
    -------
    Vector, float, float, float
        The vector is the cartesian translation (dx, dy, dz) and yaw, pitch, and roll angle in degrees required to align the measured points
        to the ideal points.
    """
    # convert from Point to stacked array (x, y, z) x N
    measured_array = [[p.x, p.y, p.z] for p in measured_points]
    ideal_array = [[p.x, p.y, p.z] for p in ideal_points]
    # Ensure the points are centered at their centroids
    measured_centroid = np.mean(measured_array, axis=0)
    ideal_centroid = np.mean(ideal_array, axis=0)
    measured_centered = measured_array - measured_centroid
    ideal_centered = ideal_array - ideal_centroid

    # Compute the covariance matrix
    H = measured_centered.T @ ideal_centered

    # Singular Value Decomposition
    U, _, Vt = np.linalg.svd(H)
    rotation_matrix = Vt.T @ U.T

    # Ensure a right-handed coordinate system
    if np.linalg.det(rotation_matrix) < 0:
        Vt[2, :] *= -1
        rotation_matrix = Vt.T @ U.T

    # Compute the euler angles
    rotation = Rotation.from_matrix(rotation_matrix)
    euler = conventional_to_euler_notation(axes_order)
    roll, pitch, yaw = rotation.as_euler(euler, degrees=True)

    rotated_measured_centroid = rotation.apply(measured_centroid)
    translation = ideal_centroid - rotated_measured_centroid

    return Vector(*translation), yaw, pitch, roll
