"""The VMAT module consists of the class VMAT, which is capable of loading an
EPID DICOM Open field image and MLC field image. The analysis is based on recommendations
from the Clif-Ling paper, Varian RapidArc QA tests and procedures, and
Varian RapidArc Dynamic QA Test Procedures for TrueBeam, covering:

* **Dose-Rate & Gantry-Speed (DRGS)** (aka T2 test)
* **Dose-Rate & MLC speed (DRMLC)** (aka T3 test)
* **Dose-Rate & Collimator speed (DRCS)** (aka T4 test / RapidArc Dynamic)

Features:

* **Do all tests** - Pylinac can handle DRGS, DRMLC or DRCS tests.
* **Automatic open/DMLC identification** - Pass in both images--don't worry about naming.
  Pylinac will automatically identify the right images.
* **Automatic offset correction** - Older VMAT tests had the ROIs offset, newer ones are centered.
  No worries, pylinac finds the ROIs automatically, with DRCS assuming a centered image.
"""

from __future__ import annotations

import copy
import enum
import math
import webbrowser
from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass
from io import BytesIO
from pathlib import Path
from typing import BinaryIO

import matplotlib.pyplot as plt
import numpy as np
from plotly import graph_objects as go
from pydantic import BaseModel, ConfigDict, Field
from scipy.ndimage import median_filter
from skimage.transform import EuclideanTransform

from . import Normalization
from .core import image
from .core.array_utils import normalize
from .core.geometry import Point, PointSerialized
from .core.image import DicomImage, ImageLike
from .core.io import TemporaryZipDirectory, get_url, retrieve_demo_file
from .core.pdf import PylinacCanvas
from .core.profile import CircleProfile, FWXMProfile
from .core.roi import RectangleROI
from .core.scale import wrap180
from .core.utilities import QuaacDatum, QuaacMixin, ResultBase, ResultsDataMixin
from .core.warnings import capture_warnings


class ImageType(enum.Enum):
    """The image type options"""

    DMLC = "dmlc"  #:
    OPEN = "open"  #:
    PROFILE = "profile"  #:


class SegmentResult(BaseModel):
    """An individual segment/ROI result"""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    passed: bool = Field(
        description="A boolean indicating if the segment passed or failed."
    )
    x_position_mm: float = Field(
        description="The position of the segment ROI in mm from CAX (lateral offset if DRGS/DRMLC, radial distance if DRCS)."
    )
    angular_position_deg: float = Field(
        description="The angle of the segment ROI in degrees."
    )
    r_corr: float = Field(
        description="R corrected (ratio)", title="R corrected (ratio)"
    )
    r_dev: float = Field(description="R deviation (%)", title="R deviation (%)")
    center_x_y: PointSerialized = Field(
        description="The center of the segment in pixel coordinates."
    )
    stdev: float = Field(
        description="The standard deviation of the segment of the ratioed images (DMLC / Open)"
    )


class CollimatorResult(BaseModel):
    """An individual Collimator line result"""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    angle_deviation: float = Field(  # measured - ideal (215.2 - 215.0 = +0.2)
        description="Collimator Deviation at angle"
    )
    angle_nominal: float = Field(
        description="The nominal angle of the collimator", title="Nominal Angle (deg)"
    )


class VMATResult(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    test_type: str = Field(
        description="The type of test that was performed as a string."
    )
    tolerance_percent: float = Field(
        description=" The tolerance used to determine if the test passed or failed."
    )
    max_deviation_percent: float = Field(
        description="The maximum deviation of any segment.", title="Max Deviation (%)"
    )
    abs_mean_deviation: float = Field(
        description="The average absolute deviation of all segments.",
        title="Absolute Mean Deviation (%)",
    )
    passed: bool = Field(
        description="A boolean indicating if the test passed or failed."
    )
    segment_data: list[SegmentResult] = Field(
        description="List of individual segment data."
    )
    named_segment_data: dict[str, SegmentResult] = Field(
        description="Named individual segment data."
    )


class DRCSResult(VMATResult):
    # this is implicitly a named_collimator_data field
    collimator_data: dict[str, CollimatorResult] = Field(
        description="List of individual collimator deviation data"
    )


class Segment(RectangleROI):
    """A class for holding and analyzing segment data of VMAT tests.

    For VMAT tests, there are either 4 or 7 'segments', which represents a section of the image that received
    radiation under the same conditions.

    Attributes
    ----------
    r_dev : float
            The reading deviation (R_dev) from the average readings of all the segments. See documentation for equation info.
    """

    r_dev: float

    def __init__(
        self,
        center_point: Point,
        width: float,
        height: float,
        ratio_image: np.ndarray,
        tolerance: float | int,
        rotation: float = 0,
    ):
        self.r_dev: float = 0.0  # is assigned after all segments constructed
        self._tolerance = tolerance
        self._ratio_image = ratio_image
        super().__init__(ratio_image, width, height, center_point, rotation)

    @property
    def r_corr(self) -> float:
        """Return the ratio of the mean pixel values of DMLC/OPEN images."""
        return self.pixels_flat.mean() * 100

    @property
    def stdev(self) -> float:
        """Return the standard deviation of the segment."""
        return self.pixels_flat.std()

    @property
    def passed(self) -> bool:
        """Return whether the segment passed or failed."""
        return abs(self.r_dev) < self._tolerance * 100

    def get_bg_color(self) -> str:
        """Get the background color of the segment when plotted, based on the pass/fail status."""
        return "blue" if self.passed else "red"


@dataclass
class CollimatorDeviation:
    """A class for holding collimator deviations of DRCS tests.

    Attributes
    ----------
    name : str
        The name of the collimator deviation line.
    angle_nominal : float
        The nominal angle of the line in degrees (IEC).
    points: tuple[Point, Point]
        The two points that make the line.
    """

    name: str
    angle_nominal: float
    points: tuple[Point, Point]

    @property
    def angle_measured(self) -> float:
        dy = self.points[1].y - self.points[0].y
        dx = self.points[1].x - self.points[0].x
        angle_im = np.arctan2(dy, dx)
        angle_iec = -(np.rad2deg(angle_im) + 90) % 360
        return angle_iec

    @property
    def angle_deviation(self) -> float:
        return wrap180(self.angle_measured - self.angle_nominal)


class VMATBase(ABC, ResultsDataMixin[VMATResult], QuaacMixin):
    _url_suffix: str
    _result_header: str
    _result_short_header: str
    roi_config: dict
    default_roi_config: dict
    dmlc_image: image.DicomImage
    open_image: image.DicomImage
    segments: list[Segment]
    _tolerance: float
    ratio_image: np.ndarray
    text_rotation: float | int

    @property
    @abstractmethod
    def default_segment_size_mm(self) -> tuple[float, float]:
        pass

    @property
    @abstractmethod
    def default_roi_config(self) -> dict:
        pass

    def __init__(
        self,
        image_paths: Sequence[str | BinaryIO | Path],
        ground=True,
        check_inversion=True,
        **kwargs,
    ):
        """
        Parameters
        ----------
        image_paths : iterable (list, tuple, etc)
            A sequence of paths to the image files.
        kwargs
            Passed to the image loading function. See :func:`~pylinac.core.image.load`.
        """
        super().__init__()
        ground = kwargs.get("ground", False) or ground
        check_inversion = kwargs.get("check_inversion", False) or check_inversion
        if len(image_paths) != 2:
            raise ValueError("Exactly 2 images (open, DMLC) must be passed")
        image1, image2 = self._load_images(image_paths, ground=ground, **kwargs)
        if check_inversion:
            image1, image2 = self._check_inversion(image1, image2)
        self._identify_images(image1, image2)
        self.segments = []
        self._tolerance = 0

    @classmethod
    def from_url(cls, url: str):
        """Load a ZIP archive from a URL.  Must follow the naming convention.

        Parameters
        ----------
        url : str
            Must point to a valid URL that is a ZIP archive of two VMAT images.
        """
        zfile = get_url(url)
        return cls.from_zip(zfile)

    @classmethod
    def from_zip(cls, path: str | Path, **kwargs):
        """Load VMAT images from a ZIP file that contains both images. Must follow the naming convention.

        Parameters
        ----------
        path : str
            Path to the ZIP archive which holds the VMAT image files.
        kwargs
            Passed to the constructor.
        """
        with TemporaryZipDirectory(path) as tmpzip:
            image_files = image.retrieve_image_files(tmpzip)
            return cls(image_paths=image_files, **kwargs)

    @classmethod
    def from_demo_images(cls, **kwargs):
        """Construct a VMAT instance using the demo images."""
        demo_file = retrieve_demo_file(name=cls._url_suffix)
        return cls.from_zip(demo_file, **kwargs)

    def analyze(
        self,
        tolerance: float | int = 1.5,
        segment_size_mm: tuple | None = None,
        roi_config: dict | None = None,
    ):
        """Analyze the open and DMLC field VMAT images, according to 1 of 2 possible tests.

        Parameters
        ----------
        tolerance : float, int, optional
            The tolerance of the sample deviations in percent. Default is 1.5.
            Must be between 0 and 8.
        segment_size_mm : tuple(int, int)
            The (width, height) of the ROI segments in mm.
        roi_config : dict
            A dict of the ROI settings. The keys are the names of the ROIs and each value is a dict containing the offset in mm 'offset_mm'.
        """
        if segment_size_mm is None:
            segment_size_mm = self.default_segment_size_mm
        if roi_config is None:
            roi_config = self.default_roi_config

        self._tolerance = tolerance / 100
        self.roi_config = roi_config
        self.ratio_image = self.dmlc_image.array / self.open_image.array

        """Analysis"""
        self._calculate_segments(segment_size_mm)
        self._update_r_corrs()

    @staticmethod
    def _load_images(
        image_paths: Sequence[str | BytesIO], ground, **kwargs
    ) -> tuple[ImageLike, ImageLike]:
        image1 = image.load(image_paths[0], **kwargs)
        image2 = image.load(image_paths[1], **kwargs)
        if ground:
            image1.ground()
            image2.ground()
        return image1, image2

    @abstractmethod
    def _identify_images(self, image1: DicomImage, image2: DicomImage):
        """Identify which image is the DMLC and which is the open field."""
        pass

    @abstractmethod
    def _calculate_segments(self, segment_size_mm: tuple[float, float]):
        """Construct the center points of the segments based on the roi_config."""
        pass

    def results(self) -> str:
        """A string of the summary of the analysis results.

        Returns
        -------
        str
            The results string showing the overall result and deviation statistics by segment.
        """
        if self.passed:
            passfail_str = "PASS"
        else:
            passfail_str = "FAIL"

        string = f"{self._result_header}\nTest Results (Tol. +/-{self._tolerance*100:2.2}%): {passfail_str}\n"

        string += f"Max Deviation: {self.max_r_deviation:2.3}%\nAbsolute Mean Deviation: {self.avg_abs_r_deviation:2.3}%"
        return string

    def _quaac_datapoints(self) -> dict[str, QuaacDatum]:
        results_data = self.results_data(as_dict=True)
        data = {
            "Max Deviation": QuaacDatum(
                value=results_data["max_deviation_percent"],
                unit="%",
            ),
            "Absolute Mean Deviation": QuaacDatum(
                value=results_data["abs_mean_deviation"],
                unit="%",
            ),
        }
        for segment, segment_data in results_data["named_segment_data"].items():
            data[f"{segment} Rcorr"] = QuaacDatum(
                value=segment_data["r_corr"],
            )
            data[f"{segment} Rdev"] = QuaacDatum(
                value=segment_data["r_dev"],
                unit="%",
            )
        return data

    def _update_r_corrs(self):
        """After the Segment constructions, the R_corr must be set for each segment."""
        avg_r_corr = np.array([segment.r_corr for segment in self.segments]).mean()
        for segment in self.segments:
            segment.r_dev = ((segment.r_corr / avg_r_corr) * 100) - 100

    @property
    def passed(self) -> bool:
        return all(segment.passed for segment in self.segments)

    @property
    def r_devs(self) -> np.ndarray:
        """Return the deviations of all segments as an array."""
        return np.array([segment.r_dev for segment in self.segments])

    @property
    def avg_abs_r_deviation(self) -> float:
        """Return the average of the absolute R_deviation values."""
        return np.abs(self.r_devs).mean()

    @property
    def avg_r_deviation(self) -> float:
        """Return the average of the R_deviation values, including the sign."""
        return self.r_devs.mean()

    @property
    def max_r_deviation(self) -> float:
        """Return the value of the maximum R_deviation segment."""
        return np.max(np.abs(self.r_devs))

    @abstractmethod
    def _roi_profiles(
        self, image1: DicomImage, image2: DicomImage
    ) -> tuple[FWXMProfile, FWXMProfile]:
        """Return two profiles from the open and DMLC image. Used for qualitative visualization and image identification."""
        pass

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

        # images
        fig_open = self.open_image.plotly(
            show=False,
            title="Open Image",
            show_colorbar=show_colorbar,
            show_legend=show_legend,
            **kwargs,
        )
        self._draw_plotly_segments(fig=fig_open)
        fig_dmlc = self.dmlc_image.plotly(
            show=False,
            title="DMLC Image",
            show_colorbar=show_colorbar,
            show_legend=show_legend,
            **kwargs,
        )
        self._draw_plotly_segments(fig=fig_dmlc)

        # ROI profiles
        dmlc_prof, open_prof = self._roi_profiles(self.dmlc_image, self.open_image)
        fig_profile = go.Figure()
        dmlc_prof.plotly(fig_profile, name="DMLC", show=False)
        open_prof.plotly(fig_profile, name="Open", show=False)
        fig_profile.update_layout(
            title={
                "text": "Median Profiles",
                "x": 0.5,
            },
            xaxis_title="Pixel",
            yaxis_title="Normalized Response",
            coloraxis_showscale=show_colorbar,
            showlegend=show_legend,
        )

        if show:
            fig_open.show()
            fig_dmlc.show()
            fig_profile.show()

        return {"Open": fig_open, "DMLC": fig_dmlc, "Profile": fig_profile}

    def plot_analyzed_image(
        self, show: bool = True, show_text: bool = True, **plt_kwargs: dict
    ):
        """Plot the analyzed images. Shows the open and dmlc images with the segments drawn; also plots the median
        profiles of the two images for visual comparison.

        Parameters
        ----------
        show : bool
            Whether to actually show the image.
        show_text : bool
            Whether to show the ROI names on the image.
        plt_kwargs : dict
            Keyword args passed to the plt.subplots() method. Allows one to set things like figure size.
        """
        fig, axes = plt.subplots(ncols=3, sharex=True, **plt_kwargs)
        subimages = (ImageType.OPEN, ImageType.DMLC, ImageType.PROFILE)
        titles = ("Open", "DMLC", "Median Profiles")
        for subimage, axis, title in zip(subimages, axes, titles):
            self._plot_analyzed_subimage(
                subimage=subimage, ax=axis, show=False, show_text=show_text
            )
            axis.set_title(title)
        axis.set_ylabel("Normalized Response")
        axis.legend(loc="lower center")

        if show:
            plt.tight_layout(h_pad=1.5)
            plt.show()

    def _save_analyzed_subimage(
        self,
        filename: str | BytesIO,
        subimage: ImageType,
        show_text: bool,
        **kwargs,
    ):
        """Save the analyzed images as a png file.

        Parameters
        ----------
        filename : str, file-object
            Where to save the file to.
        kwargs
            Passed to matplotlib.
        """
        self._plot_analyzed_subimage(subimage=subimage, show=False, show_text=show_text)
        plt.savefig(filename, **kwargs)

    def _plot_analyzed_subimage(
        self,
        subimage: ImageType,
        show: bool = True,
        ax: plt.Axes | None = None,
        show_text: bool = True,
    ):
        """Plot an individual piece of the VMAT analysis.

        Parameters
        ----------
        subimage : str
            Specifies which image to plot.
        show : bool
            Whether to actually plot the image.
        ax : matplotlib Axes, None
            If None (default), creates a new figure to plot to, otherwise plots to the given axes.
        show_text : bool
            Whether to show the ROI names on the image.
        """
        plt.ioff()
        if ax is None:
            fig, ax = plt.subplots()

        # plot DMLC or OPEN image
        if subimage in (ImageType.DMLC, ImageType.OPEN):
            if subimage == ImageType.DMLC:
                img = self.dmlc_image
            elif subimage == ImageType.OPEN:
                img = self.open_image
            img.plot(ax=ax, show=False)
            self._draw_segments(ax, show_text)
            plt.sca(ax)
            plt.axis("off")
            plt.tight_layout()

        # plot profile
        elif subimage == ImageType.PROFILE:
            dmlc_prof, open_prof = self._roi_profiles(self.dmlc_image, self.open_image)
            ax.plot(dmlc_prof.values, label="DMLC")
            ax.plot(open_prof.values, label="Open")
            ax.autoscale(axis="x", tight=True)
            ax.legend(loc=8, fontsize="large")
            ax.grid()

        if show:
            plt.show()

    def _draw_plotly_segments(self, fig: go.Figure) -> None:
        """Draw the segments onto a plotly figure.

        Parameters
        ----------
        fig : go.Figure
            The figure to draw the objects on.
        """
        for segment, roi_name in zip(self.segments, self.roi_config.keys()):
            segment.plotly(
                fig,
                line_color=segment.get_bg_color(),
                name=f"{roi_name} ({segment.r_dev:2.2f}%)",
            )

    def _draw_segments(self, axis: plt.Axes, show_text: bool):
        """Draw the segments onto a plot.

        Parameters
        ----------
        axis : matplotlib.axes.Axes
            The plot to draw the objects on.
        show_text : bool
            Whether to show the ROI name on the image
        """
        for segment, (roi_name, roi_config) in zip(
            self.segments, self.roi_config.items()
        ):
            color = segment.get_bg_color()
            if show_text:
                text = f"{roi_name} : {segment.r_dev:2.2f}%"
            else:
                text = ""
            segment.plot2axes(
                axis,
                edgecolor=color,
                text=text,
                text_rotation=self.text_rotation,
                fontsize="small",
            )

    def publish_pdf(
        self,
        filename: str,
        notes: str = None,
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
        canvas = PylinacCanvas(
            filename=filename,
            page_title=f"{self._result_short_header} VMAT Analysis",
            metadata=metadata,
            logo=logo,
        )
        for y, x, width, img in zip(
            (9, 9, -2),
            (1, 11, 3),
            (9, 9, 14),
            (ImageType.OPEN, ImageType.DMLC, ImageType.PROFILE),
        ):
            data = BytesIO()
            self._save_analyzed_subimage(data, subimage=img, show_text=True)
            canvas.add_image(data, location=(x, y), dimensions=(width, 18))
            # canvas.add_text(text=f"{img} Image", location=(x + 2, y + 10), font_size=18)
        canvas.add_text(text="Open Image", location=(4, 22), font_size=18)
        canvas.add_text(text=f"{self.open_image.base_path}", location=(4, 21.5))
        canvas.add_text(text="DMLC Image", location=(14, 22), font_size=18)
        canvas.add_text(text=f"{self.dmlc_image.base_path}", location=(14, 21.5))
        canvas.add_text(text="Median profiles", location=(8, 12), font_size=18)
        text = [
            f"{self._result_header} VMAT results:",
            f"Source-to-Image Distance (mm): {self.open_image.sid:2.0f}",
            f"Tolerance (%): {self._tolerance*100:2.1f}",
            f"Absolute mean deviation (%): {self.avg_abs_r_deviation:2.2f}",
            f"Maximum deviation (%): {self.max_r_deviation:2.2f}",
        ]
        canvas.add_text(text=text, location=(10, 25.5))
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 5.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 5))

        canvas.finish()

        if open_file:
            webbrowser.open(filename)

    @staticmethod
    def _check_inversion(image1, image2):
        for img in (image1, image2):
            img.check_inversion()
        return image1, image2


class VMATLinearBase(VMATBase, ABC):
    """Class representing linear VMAT tests:
    - DRGS: Dose-Rate vs Gantry-speed
    - DRMLC: Dose-Rate vs MLC-speed
    Will accept, analyze, and return the results."""

    text_rotation = 90  # rotation of text on image

    @property
    def default_segment_size_mm(self) -> tuple[float, float]:
        return 5, 100

    def _identify_images(self, image1: DicomImage, image2: DicomImage):
        """Identify which image is the DMLC and which is the open field."""
        profile1, profile2 = self._roi_profiles(image1=image1, image2=image2)
        field_profile1 = profile1.field_values()
        field_profile2 = profile2.field_values()
        # first check if the profiles have a very different length
        # if so, the longer one is the open field
        # this leverages the shortcoming in FWXMProfile where the field might be very small because
        # it "caught" on one of the first dips of the DMLC image
        # catches most often with Halcyon images
        if abs(len(field_profile1) - len(field_profile2)) > min(
            len(field_profile1), len(field_profile2)
        ):
            if len(field_profile1) > len(field_profile2):
                self.open_image = image1
                self.dmlc_image = image2
            else:
                self.open_image = image2
                self.dmlc_image = image1
        # normal check of the STD compared; for flat-ish beams this works well.
        elif np.std(field_profile1) > np.std(field_profile2):
            self.dmlc_image = image1
            self.open_image = image2
        else:
            self.dmlc_image = image2
            self.open_image = image1

    def _roi_profiles(
        self, image1: DicomImage, image2: DicomImage
    ) -> list[FWXMProfile]:
        profiles: list[FWXMProfile] = []
        for orig_img in (image1, image2):
            img = copy.deepcopy(orig_img)
            img.ground()
            img.check_inversion()
            profile = FWXMProfile(
                np.mean(img.array, axis=0),
                ground=True,
                normalization=Normalization.BEAM_CENTER,
            )
            profile.stretch()
            norm_val = np.percentile(profile.values, 90)
            profile.normalize(norm_val)
            profiles.append(profile)
        return profiles

    def _generate_results_data(self) -> VMATResult:
        """Present the results data and metadata as a dataclass or dict.
        The default return type is a dataclass."""
        segment_data = []
        named_segment_data = {}
        for segment, (roi_name, roi_data) in zip(
            self.segments, self.roi_config.items()
        ):
            segment = SegmentResult(
                passed=segment.passed,
                r_corr=segment.r_corr,
                r_dev=segment.r_dev,
                center_x_y=segment.center,
                x_position_mm=roi_data["offset_mm"],
                stdev=segment.stdev,
                angular_position_deg=0,
            )
            segment_data.append(segment)
            named_segment_data[roi_name] = segment
        return VMATResult(
            test_type=self._result_header,
            tolerance_percent=self._tolerance * 100,
            max_deviation_percent=self.max_r_deviation,
            abs_mean_deviation=self.avg_abs_r_deviation,
            passed=self.passed,
            segment_data=segment_data,
            named_segment_data=named_segment_data,
        )

    def _calculate_segments(self, segment_size_mm: tuple[float, float]):
        """Construct the center points of the segments based on the field center and known x-offsets."""
        y = self.open_image.center.y
        _, open_prof = self._roi_profiles(self.dmlc_image, self.open_image)
        x_field_center = round(open_prof.center_idx)
        dpmm = self.open_image.dpmm
        for roi_data in self.roi_config.values():
            x_offset_mm = roi_data["offset_mm"]
            x_offset_pixels = x_offset_mm * dpmm
            x = x_field_center + x_offset_pixels
            segment = Segment(
                Point(x, y),
                width=segment_size_mm[0] * dpmm,
                height=segment_size_mm[1] * dpmm,
                ratio_image=self.ratio_image,
                tolerance=self._tolerance,
            )
            self.segments.append(segment)


@capture_warnings
class DRGS(VMATLinearBase):
    """Class representing a Dose-Rate, Gantry-speed VMAT test. Will accept, analyze, and return the results."""

    _url_suffix = "drgs.zip"
    _result_header = "Dose Rate & Gantry Speed"
    _result_short_header = "DR/GS"

    @property
    def default_roi_config(self) -> dict:
        return {
            "ROI 1": {"offset_mm": -60},
            "ROI 2": {"offset_mm": -40},
            "ROI 3": {"offset_mm": -20},
            "ROI 4": {"offset_mm": 0},
            "ROI 5": {"offset_mm": 20},
            "ROI 6": {"offset_mm": 40},
            "ROI 7": {"offset_mm": 60},
        }

    @staticmethod
    def run_demo():
        """Run the demo for the Dose Rate & Gantry Speed test."""
        vmat = DRGS.from_demo_images()
        vmat.analyze()  # old images (rev1, not new rev2's), which are offset
        print(vmat.results())
        vmat.plot_analyzed_image()


@capture_warnings
class DRMLC(VMATLinearBase):
    """Class representing a Dose-Rate, MLC speed VMAT test. Will accept, analyze, and return the results."""

    _url_suffix = "drmlc.zip"
    _result_header = "Dose Rate & MLC Speed"
    _result_short_header = "DR/MLCS"

    @property
    def default_roi_config(self) -> dict:
        return {
            "ROI 1": {"offset_mm": -45},
            "ROI 2": {"offset_mm": -15},
            "ROI 3": {"offset_mm": 15},
            "ROI 4": {"offset_mm": 45},
        }

    @staticmethod
    def run_demo():
        """Run the demo for the MLC leaf speed test."""
        vmat = DRMLC.from_demo_images()
        vmat.analyze()
        print(vmat.results())
        vmat.plot_analyzed_image()


@capture_warnings
class DRCS(VMATBase):
    """Class representing a Dose-Rate, Collimator speed VMAT test.
    Will accept, analyze, and return the results."""

    collimator_deviations = list[float]
    text_rotation = 0  # rotation of text on image

    _url_suffix = "drcs.zip"
    _result_header = "Dose Rate & Collimator Speed"
    _result_short_header = "DR/CS"
    _default_radial_distance = 50  # in mm

    @property
    def default_segment_size_mm(self) -> tuple[float, float]:
        return 40, 10

    @property
    def default_roi_config(self) -> dict:
        return {
            "ROI 1": {"radial_distance": self._default_radial_distance, "angle": -120},
            "ROI 2": {"radial_distance": self._default_radial_distance, "angle": -60},
            "ROI 3": {"radial_distance": self._default_radial_distance, "angle": 0},
            "ROI 4": {"radial_distance": self._default_radial_distance, "angle": 60},
            "ROI 5": {"radial_distance": self._default_radial_distance, "angle": 120},
        }

    @property
    def default_collimator_config(self) -> dict[str, float]:
        return {"A": 210, "B": 270, "C": 330, "D": 30, "E": 90, "F": 150}  # IEC

    @property
    def default_collimator_radial_distances(self) -> tuple[float, float]:
        return 30, 70  # mm

    def analyze(
        self,
        tolerance: float | int = 1.5,  # Segments, in %
        segment_size_mm: tuple | None = None,
        roi_config: dict | None = None,
        collimator_radial_distances: dict | None = None,
        collimator_config: dict | None = None,
    ):
        super().analyze(tolerance, segment_size_mm, roi_config)
        cc = collimator_config or self.default_collimator_config
        crd = collimator_radial_distances or self.default_collimator_radial_distances
        self._calculate_collimator_deviations(cc, crd)

    def _identify_images(self, image1: DicomImage, image2: DicomImage):
        """Identify which image is the DMLC and which is the open field.

        Notes
        -----

        In DRCS, the boundaries between regions on the DMLC image have higher intensity.
        Normalizing to the max will produce a lower mean value. Furthermore, for the
        example images provided, the open image is a full circle, whereas the DMLC image
        has a pie slice missing, lowering the mean value even further.
        """
        filter_size = 10
        sum1 = normalize(median_filter(image1.array, filter_size)).sum()
        sum2 = normalize(median_filter(image2.array, filter_size)).sum()

        if sum1 > sum2:
            self.open_image = image1
            self.dmlc_image = image2
        else:
            self.open_image = image2
            self.dmlc_image = image1

    def _generate_results_data(self) -> VMATResult:
        """Present the results data and metadata as a dataclass or dict.
        The default return type is a dataclass."""
        segment_data = []
        named_segment_data = {}
        for segment, (roi_name, roi_data) in zip(
            self.segments, self.roi_config.items()
        ):
            segment = SegmentResult(
                passed=segment.passed,
                r_corr=segment.r_corr,
                r_dev=segment.r_dev,
                center_x_y=segment.center,
                x_position_mm=roi_data["radial_distance"],
                stdev=segment.stdev,
                angular_position_deg=roi_data["angle"],
            )
            segment_data.append(segment)
            named_segment_data[roi_name] = segment

        coll_data = {}
        for cd in self.collimator_deviations:
            coll_data[cd.name] = CollimatorResult(
                angle_deviation=cd.angle_deviation,
                angle_nominal=cd.angle_nominal,
            )

        return DRCSResult(
            test_type=self._result_header,
            tolerance_percent=self._tolerance * 100,
            max_deviation_percent=self.max_r_deviation,
            abs_mean_deviation=self.avg_abs_r_deviation,
            passed=self.passed,
            segment_data=segment_data,
            named_segment_data=named_segment_data,
            collimator_data=coll_data,
        )

    def _calculate_segments(self, segment_size_mm: tuple[float, float]):
        """Calculate the segments based on ROI config.

        Notes
        -----
        This assumes that the ratio image is centered to the central pixel.
        Once we have more data we can verify this assumption and change if necessary.
        """
        dpmm = self.open_image.dpmm
        image_center = self.open_image.center.as_array(("x", "y"))
        im_translation = EuclideanTransform(translation=image_center)
        for roi_data in self.roi_config.values():
            radial_distance = roi_data["radial_distance"]
            radial_distance_px = radial_distance * dpmm
            roi_translation = EuclideanTransform(translation=(radial_distance_px, 0))

            coll_angle = roi_data["angle"]
            im_angle = -coll_angle - 90
            angle_rad = np.deg2rad(im_angle)
            roi_rotation = EuclideanTransform(rotation=angle_rad)

            roi_tf = roi_translation + roi_rotation  # extrinsic translation -> rotation
            tf = roi_tf + im_translation  # roi in image coordinates

            segment = Segment(
                center_point=Point(tf.translation),
                width=segment_size_mm[0] * dpmm,
                height=segment_size_mm[1] * dpmm,
                ratio_image=self.ratio_image,
                tolerance=self._tolerance,
                rotation=np.rad2deg(tf.rotation),
            )
            self.segments.append(segment)

    def _calculate_collimator_deviations(
        self,
        collimator_config: dict[str, float],
        collimator_radial_distances: tuple[float, float],
    ):
        num_angles = len(collimator_config)
        nominal_angles = np.fromiter(collimator_config.values(), dtype=float)
        max_diff_angle = max(np.abs(wrap180(np.diff(nominal_angles))))
        crd_px = np.array(collimator_radial_distances) * self.dmlc_image.dpmm
        peaks = list()
        for crd in crd_px:
            circle_profile = CircleProfile(
                center=self.dmlc_image.center,
                radius=crd,
                image_array=self.ratio_image,
                start_angle=math.pi / 2,  # start in the "empty" region and go CCW
            )
            min_distance = 2 * np.pi * crd / 360 * 0.9 * max_diff_angle
            circle_profile.find_peaks(min_distance=min_distance, max_number=num_angles)
            if len(circle_profile.peaks) != num_angles:
                raise ValueError("Could not detect collimator lines.")
            peaks.append(circle_profile.peaks)

        self.collimator_deviations = list[CollimatorDeviation]()
        for config, points in zip(collimator_config.items(), np.array(peaks).T):
            cd = CollimatorDeviation(config[0], config[1], (points[0], points[1]))
            self.collimator_deviations.append(cd)

    @staticmethod
    def run_demo():
        """Run the demo for the Dose Rate & Collimator Speed test."""
        vmat = DRCS.from_demo_images()
        vmat.analyze()
        print(vmat.results())
        vmat.plot_analyzed_image()

    def _roi_profiles(
        self, image1: DicomImage, image2: DicomImage
    ) -> list[FWXMProfile]:
        """Return two median profiles from the open and DMLC image. Compared
        to the linear VMAT tests, we first extract a circular profile, then convert it to a linear FWXM profile."""
        profiles: list[FWXMProfile] = []
        for orig_img in (image1, image2):
            img = copy.deepcopy(orig_img)
            # we need one single radius; just take the first of the ROIs
            # It's possible the user has overloaded the ROIs to be different distances, but that's not our problem
            radius_px = (
                list(self.roi_config.values())[0]["radial_distance"]
                * self.dmlc_image.dpmm
            )
            circle_profile = CircleProfile(
                center=img.center,
                radius=radius_px,
                image_array=img.array,
                start_angle=math.pi / 2,  # start in the "empty" region and go CCW
            )
            # due to signature and implementation differences in later calls, it's easiest to directly create a FWXMProfile from the circle profile values
            profile = FWXMProfile(
                values=circle_profile.values,
                ground=False,
                normalization=Normalization.NONE,
            )
            # normalize so the profiles are visually comparable
            profile.normalize(np.percentile(profile.values, 50))
            profiles.append(profile)
        return profiles
