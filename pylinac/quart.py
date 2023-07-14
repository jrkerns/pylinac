from __future__ import annotations

import dataclasses
import io
import textwrap
import webbrowser
from io import BytesIO
from pathlib import Path

import numpy as np
import scipy.ndimage
from matplotlib import pyplot as plt

from .core import pdf
from .core.geometry import Line, Point
from .core.profile import Interpolation, SingleProfile
from .core.utilities import ResultBase
from .ct import (
    AIR,
    CTP404CP504,
    CTP486,
    CatPhanBase,
    CatPhanModule,
    ThicknessROI,
    rois_to_results,
)

UNIFORMITY_OFFSET_MM = -45
GEOMETRY_OFFSET_MM = 45
ACRYLIC = 120
POLY = -35
TEFLON = 990


@dataclasses.dataclass
class QuartHUModuleOutput:
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    offset: int
    roi_settings: dict
    rois: dict
    measured_slice_thickness_mm: float
    signal_to_noise: float
    contrast_to_noise: float


@dataclasses.dataclass
class QuartGeometryModuleOutput:
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    offset: int
    roi_settings: dict
    rois: dict
    distances: dict


@dataclasses.dataclass
class QuartUniformityModuleOutput:
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    offset: int
    roi_settings: dict
    rois: dict
    passed: bool


@dataclasses.dataclass
class QuartDVTResult(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    phantom_model: str  #:
    phantom_roll_deg: float  #:
    origin_slice: int  #:
    num_images: int  #:
    hu_module: QuartHUModuleOutput  #:
    uniformity_module: QuartUniformityModuleOutput  #:
    geometric_module: QuartGeometryModuleOutput  #:


class QuartHUModule(CTP404CP504):
    roi_dist_mm = 52.5
    roi_radius_mm = 6
    roi_settings = {
        "Air": {
            "value": AIR,
            "angle": -90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Poly": {
            "value": POLY,
            "angle": 0,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Acrylic": {
            "value": ACRYLIC,
            "angle": 45,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Teflon": {
            "value": TEFLON,
            "angle": 180,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
    }
    background_roi_settings = {}
    thickness_roi_height = 25
    thickness_roi_width = 15
    thickness_roi_distance_mm = 32
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

    def _setup_geometry_rois(self) -> None:
        # no geometry ROIs
        pass

    def _setup_thickness_rois(self) -> None:
        """We invert the thickness ROIs because they are air gaps, not high-density wires"""
        self.thickness_image.invert()
        for name, setting in self.thickness_roi_settings.items():
            self.thickness_rois[name] = ThicknessROI(
                self.thickness_image,
                setting["width_pixels"],
                setting["height_pixels"],
                setting["angle_corrected"],
                setting["distance_pixels"],
                self.phan_center,
            )

    @property
    def meas_slice_thickness(self) -> float:
        """The average slice thickness for the 4 wire measurements in mm."""
        INCLINATION_CORRECTION = 0.577  # per manual; tan(30)
        return np.mean(
            sorted(
                roi.wire_fwhm * self.mm_per_pixel * INCLINATION_CORRECTION
                for roi in self.thickness_rois.values()
            )
        ) / (1 + 2 * self.pad)

    @property
    def signal_to_noise(self) -> float:
        """Calculate the SNR based on the suggested procedure in the manual:
        SNR = (HU + 1000) / sigma,
        where HU is the mean HU of a chosen insert and sigma is the stdev of the HU insert.
        We choose to use the Polystyrene as the target HU insert"""
        return (self.rois["Poly"].pixel_value + 1000) / self.rois["Poly"].std

    @property
    def contrast_to_noise(self) -> float:
        """Calculate the CNR based on the suggested procedure in the manual:
        CNR = abs(HU_target - HU_background) / sigma,
        where HU_target is the mean HU of a chosen insert, HU_background is the mean HU of the background insert
        and sigma is the stdev of the HU background.
        We choose to use the Polystyrene as the target HU insert and Acrylic (base phantom material) as the background
        """
        return (
            abs(self.rois["Poly"].pixel_value - self.rois["Acrylic"].pixel_value)
            / self.rois["Acrylic"].std
        )


class QuartUniformityModule(CTP486):
    """Class for analysis of the Uniformity slice of the CTP module. Measures 5 ROIs around the slice that
    should all be close to the same value.
    """

    common_name = "HU Uniformity"
    roi_dist_mm = 53
    roi_radius_mm = 10
    nominal_value = 120
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


class QuartGeometryModule(CatPhanModule):
    """Class for analysis of the Uniformity slice of the CTP module. Measures 5 ROIs around the slice that
    should all be close to the same value.
    """

    attr_name = "geometry_module"
    common_name = "Geometric Distortion"
    profiles: dict

    def _setup_rois(self) -> None:
        self.profiles = {}
        img = (
            self.image.array.copy()
        )  # we copy so we don't overwrite the existing image pixels
        img = scipy.ndimage.median_filter(img, size=3)
        img = img - img.min()  # ground the profile
        # calculate horizontal
        data = img[int(self.phan_center.y), :]
        prof = SingleProfile(
            data, interpolation=Interpolation.NONE, dpmm=1 / self.mm_per_pixel
        )
        fwhm = prof.fwxm_data()
        line = Line(
            Point(fwhm["left index (rounded)"], self.phan_center.y),
            Point(fwhm["right index (rounded)"], self.phan_center.y),
        )
        self.profiles["horizontal"] = {
            "width (mm)": fwhm["width (exact) mm"],
            "line": line,
        }
        # calculate vertical
        data = img[:, int(self.phan_center.x)]
        prof = SingleProfile(
            data, interpolation=Interpolation.NONE, dpmm=1 / self.mm_per_pixel
        )
        fwhm = prof.fwxm_data()
        line = Line(
            Point(self.phan_center.x, fwhm["left index (rounded)"]),
            Point(self.phan_center.x, fwhm["right index (rounded)"]),
        )
        self.profiles["vertical"] = {
            "width (mm)": fwhm["width (exact) mm"],
            "line": line,
        }

    def plot_rois(self, axis: plt.Axes):
        for name, profile_data in self.profiles.items():
            profile_data["line"].plot2axes(axis, width=2, color="blue")

    def distances(self) -> dict[str, float]:
        """The measurements of the phantom size for the two lines in mm"""
        return {f"{name} mm": p["width (mm)"] for name, p in self.profiles.items()}


class QuartDVT(CatPhanBase):
    """A class for loading and analyzing CT DICOM files of a Quart phantom that comes with the Halcyon.
    Analyzes: HU Uniformity, Image Scaling & HU Linearity.
    """

    _demo_url = "quart.zip"
    _model = "Quart DVT"
    hu_origin_slice_variance = 300
    catphan_radius_mm = 80
    hu_module: QuartHUModule
    uniformity_module: QuartUniformityModule
    geometry_module: QuartGeometryModule

    @staticmethod
    def run_demo(show: bool = True):
        """Run the Quart algorithm with a head dataset."""
        quart = QuartDVT.from_demo_images()
        quart.analyze()
        print(quart.results())
        quart.plot_analyzed_image(show)

    def analyze(
        self,
        hu_tolerance: int | float = 40,
        scaling_tolerance: int | float = 1,
        thickness_tolerance: int | float = 0.2,
        cnr_threshold: int | float = 5,
    ):
        self.localize()
        self.hu_module = QuartHUModule(
            self,
            offset=0,
            hu_tolerance=hu_tolerance,
            thickness_tolerance=thickness_tolerance,
            scaling_tolerance=scaling_tolerance,
        )
        self.uniformity_module = QuartUniformityModule(
            self, offset=UNIFORMITY_OFFSET_MM, tolerance=hu_tolerance
        )
        self.geometry_module = QuartGeometryModule(
            self, tolerance=3, offset=GEOMETRY_OFFSET_MM
        )

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
        grid_size = (2, 3)
        hu_ax = plt.subplot2grid(grid_size, (0, 1))
        self.hu_module.plot(hu_ax)
        hu_lin_ax = plt.subplot2grid(grid_size, (0, 2))
        self.hu_module.plot_linearity(hu_lin_ax)
        unif_ax = plt.subplot2grid(grid_size, (1, 0))
        self.uniformity_module.plot(unif_ax)
        unif_prof_ax = plt.subplot2grid(grid_size, (1, 2))
        self.uniformity_module.plot_profiles(unif_prof_ax)
        geometry_ax = plt.subplot2grid(grid_size, (0, 0))
        self.geometry_module.plot(geometry_ax)
        side_view_ax = plt.subplot2grid(grid_size, (1, 1))
        self.plot_side_view(side_view_ax)

        # finish up
        plt.tight_layout()
        if show:
            plt.show()

    def plot_analyzed_subimage(self, *args, **kwargs) -> None:
        raise NotImplementedError()

    def results(self, as_str: bool = True) -> str | tuple[str, ...]:
        """Return the results of the analysis as a string. Use with print()."""
        items = (
            f"\n - {self._model} QA Test - \n",
            f"HU Linearity ROIs: {self.hu_module.roi_vals_as_str}\n",
            f"HU Passed?: {self.hu_module.passed_hu}\n",
            f"Measured Slice Thickness (mm): {self.hu_module.meas_slice_thickness:2.3f}\n",
            f"Slice Thickness Passed? {self.hu_module.passed_thickness}\n",
            f"Uniformity ROIs: {self.uniformity_module.roi_vals_as_str}\n",
            f"Uniformity Passed?: {self.uniformity_module.overall_passed}\n",
            f"Geometric width: {self.geometry_module.distances()}",
        )
        if as_str:
            return "\n".join(items)
        else:
            return items

    def results_data(self, as_dict: bool = False) -> QuartDVTResult | dict:
        """Return results in a data structure for more programmatic use."""
        data = QuartDVTResult(
            phantom_model=self._model,
            phantom_roll_deg=self.catphan_roll,
            origin_slice=self.origin_slice,
            num_images=self.num_images,
            uniformity_module=QuartUniformityModuleOutput(
                offset=UNIFORMITY_OFFSET_MM,
                roi_settings=self.uniformity_module.roi_settings,
                rois=rois_to_results(self.uniformity_module.rois),
                passed=self.uniformity_module.overall_passed,
            ),
            geometric_module=QuartGeometryModuleOutput(
                offset=GEOMETRY_OFFSET_MM,
                roi_settings=self.geometry_module.roi_settings,
                rois=rois_to_results(self.geometry_module.rois),
                distances=self.geometry_module.distances(),
            ),
            hu_module=QuartHUModuleOutput(
                offset=0,
                roi_settings=self.hu_module.roi_settings,
                rois=rois_to_results(self.hu_module.rois),
                measured_slice_thickness_mm=self.hu_module.meas_slice_thickness,
                signal_to_noise=self.hu_module.signal_to_noise,
                contrast_to_noise=self.hu_module.contrast_to_noise,
            ),
        )
        if as_dict:
            return dataclasses.asdict(data)
        else:
            return data

    def plot_images(self, show: bool = True, **plt_kwargs) -> dict[str, plt.Figure]:
        """Plot all the individual images separately.

        Parameters
        ----------
        show
            Whether to show the images.
        plt_kwargs
            Keywords to pass to matplotlib for figure customization.
        """
        figs = {}
        # plot the images
        modules = {
            "HU linearity": self.hu_module,
            "HU uniformity": self.uniformity_module,
            "Geometry": self.geometry_module,
        }
        for key, module in modules.items():
            fig, ax = plt.subplots(**plt_kwargs)
            module.plot(ax)
            figs[key] = fig
        # add side-view
        fig, ax = plt.subplots(**plt_kwargs)
        self.plot_side_view(ax)
        figs["side"] = fig

        if show:
            plt.show()
        return figs

    def save_images(
        self,
        directory: Path | str | None = None,
        to_stream: bool = False,
        **plt_kwargs,
    ) -> list[Path] | dict[str, BytesIO]:
        """Save separate images to disk or stream.

        Parameters
        ----------
        directory
            The directory to write the images to. If None, will use current working directory
        to_stream
            Whether to write to stream or disk. If True, will return streams. Directory is ignored in that scenario.
        plt_kwargs
            Keywords to pass to matplotlib for figure customization.
        """
        figs = self.plot_images(show=False, **plt_kwargs)
        paths = []
        streams = {}
        for name, fig in figs.items():
            if to_stream:
                path = io.BytesIO()
            else:
                destination = Path(directory) or Path.cwd()
                path = (destination / name).with_suffix(".png").absolute()
            fig.savefig(path)
            paths.append(path)
            streams[name] = path
        if to_stream:
            return streams
        else:
            return paths

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
        analysis_title = f"{self._model} Analysis"
        analysis_images = self.save_images(to_stream=True)

        canvas = pdf.PylinacCanvas(
            filename, page_title=analysis_title, metadata=metadata, logo=logo
        )
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 4.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 4))

        shortened_texts = [
            textwrap.wrap(r, width=110) for r in self.results(as_str=False)
        ]
        idx = 0
        for items in enumerate(shortened_texts):
            for text in items:
                canvas.add_text(text=text, location=(1.5, 25 - idx * 0.5))
                idx += 1
        for page, img in enumerate(analysis_images.values()):
            canvas.add_new_page()
            canvas.add_image(img, location=(1, 5), dimensions=(18, 18))
        canvas.finish()

        if open_file:
            webbrowser.open(filename)

    def _detected_modules(self) -> list[CatPhanModule]:
        return [self.uniformity_module, self.hu_module, self.geometry_module]
