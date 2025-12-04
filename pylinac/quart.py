from __future__ import annotations

import io
import textwrap
import warnings
import webbrowser
from io import BytesIO
from pathlib import Path

import numpy as np
import scipy.ndimage
from matplotlib import pyplot as plt
from plotly import graph_objects as go
from pydantic import BaseModel, Field
from scipy.interpolate import interp1d

from .core import pdf
from .core.geometry import Line, Point
from .core.profile import FWXMProfilePhysical
from .core.utilities import ResultBase, ResultsDataMixin
from .core.warnings import capture_warnings
from .ct import (
    AIR,
    CTP404CP504,
    CTP486,
    WATER,
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


class QuartHUModuleOutput(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    offset: int = Field(
        description="The offset of the module slice in mm from the origin slice."
    )
    roi_settings: dict = Field(description="A dictionary of the ROI settings.")
    rois: dict = Field(description="A dictionary of ROI results.")
    measured_slice_thickness_mm: float = Field(
        description="The measured slice thickness in mm.",
        title="Measured Slice Thickness (mm)",
    )
    signal_to_noise: float = Field(
        description="The signal to noise ratio.", title="SNR (Poly)"
    )
    contrast_to_noise: float = Field(
        description="The contrast to noise ratio.", title="CNR (Poly/Acrylic)"
    )


class QuartGeometryModuleOutput(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    offset: int = Field(
        description="The offset of the module slice in mm from the origin slice."
    )
    roi_settings: dict = Field(description="A dictionary of the ROI settings.")
    rois: dict = Field(description="A dictionary of ROI results.")
    distances: dict = Field(
        description="A dictionary of the phantom size itself in horizontal and vertical dimensions in mm."
    )
    high_contrast_distances: dict = Field(
        description="A dictionary of the high contrast distances in mm. The key is the region of the line and the value is the distance in mm."
    )
    mean_high_contrast_distance: float = Field(
        description="The mean of the high contrast distances in mm. Four edges are measured and averaged. The absolute distance from -700HU to -200HU is measured.",
        title="Mean Distance -700->-200HU (mm)",
    )


class QuartUniformityModuleOutput(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    offset: int = Field(
        description="The offset of the module slice in mm from the origin slice."
    )
    roi_settings: dict = Field(description="A dictionary of the ROI settings.")
    rois: dict = Field(description="A dictionary of ROI results.")
    passed: bool = Field(description="A boolean indicating if the module passed.")


class QuartDVTResult(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    phantom_model: str = Field(
        description="The model of the phantom, e.g. 'Quart DVT'."
    )
    phantom_roll_deg: float = Field(
        description="The roll of the phantom in degrees.",
        title="Quart roll (\N{DEGREE SIGN})",
    )
    origin_slice: int = Field(description="The slice number of the origin image.")
    num_images: int = Field(description="The number of images given in the dataset.")
    hu_module: QuartHUModuleOutput = Field(
        description="The HU module output.", title="HU module"
    )
    uniformity_module: QuartUniformityModuleOutput = Field(
        description="The Uniformity module output.", title="Uniformity module"
    )
    geometric_module: QuartGeometryModuleOutput = Field(
        description="The Geometric module output.", title="Geometry module"
    )


class QuartHUModule(CTP404CP504):
    roi_dist_mm = 52.5
    roi_radius_mm = 6
    vial_radius_mm = 12
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
        "Water": {
            "value": WATER,
            "angle": -45,
            "distance": roi_dist_mm,
            "radius": vial_radius_mm,
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

    def _setup_rois(self) -> None:
        """On the Quart v2 phantom, there is an optional water slot.
        We need to distinguish between:
        v1, housing material: acrylic ~ 100 HU
        v2, with full vial: water ~ 0 HU
        v2, with empty vial: air ~ -1000 HU
        """
        super()._setup_rois()
        nominal_water_hu = 0  # HU
        threshold = 50  # HU
        measured_water_hu = self.rois["Water"].pixel_value
        if np.abs(measured_water_hu - nominal_water_hu) > threshold:  # not water
            self.rois.pop("Water")

    def _setup_geometry_rois(self) -> None:
        # no geometry ROIs
        pass

    def _setup_thickness_rois(self) -> None:
        """We invert the thickness ROIs because they are air gaps, not high-density wires"""
        self.thickness_image.invert()
        for name, setting in self.thickness_roi_settings.items():
            self.thickness_rois[name] = ThicknessROI.from_phantom_center(
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


class HypersightQuartHUModule(QuartHUModule):
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
        "Water": {
            "value": WATER,
            "angle": -45,
            "distance": roi_dist_mm,
            "radius": 12,
        },
    }

    def _setup_rois(self) -> None:
        CTP404CP504._setup_rois(self)  # call grandparent method


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
    horiz_array: np.ndarray
    vert_array: np.ndarray

    def _setup_rois(self) -> None:
        self.profiles = {}
        img = (
            self.image.array.copy()
        )  # we copy so we don't overwrite the existing image pixels
        img = scipy.ndimage.median_filter(img, size=3)
        img = img - img.min()  # ground the profile
        # calculate horizontal
        self.horiz_array = img[int(self.phan_center.y), :]
        prof = FWXMProfilePhysical(
            values=self.horiz_array,
            dpmm=1 / self.mm_per_pixel,
        )
        line = Line(
            Point(round(prof.field_edge_idx("left")), self.phan_center.y),
            Point(round(prof.field_edge_idx("right")), self.phan_center.y),
        )
        self.profiles["horizontal"] = {
            "width (mm)": prof.field_width_mm,
            "line": line,
        }
        # calculate vertical
        self.vert_array = img[:, int(self.phan_center.x)]
        prof = FWXMProfilePhysical(
            values=self.vert_array,
            dpmm=1 / self.mm_per_pixel,
        )
        line = Line(
            Point(self.phan_center.x, round(prof.field_edge_idx("left"))),
            Point(self.phan_center.x, round(prof.field_edge_idx("right"))),
        )
        self.profiles["vertical"] = {
            "width (mm)": prof.field_width_mm,
            "line": line,
        }

    def plot_rois(self, axis: plt.Axes):
        for name, profile_data in self.profiles.items():
            profile_data["line"].plot2axes(axis, width=2, color="blue")

    def plotly_rois(self, fig: go.Figure) -> None:
        for name, profile_data in self.profiles.items():
            profile_data["line"].plotly(fig, line_width=2, color="blue", name=name)

    def distances(self) -> dict[str, float]:
        """The measurements of the phantom size for the two lines in mm"""
        return {f"{name} mm": p["width (mm)"] for name, p in self.profiles.items()}

    def high_contrast_resolutions(self) -> dict:
        """The distance in mm from the -700 HU index to the -200 HU index.

        This calculates the distance on each edge of the horizontal and vertical
        geometric profiles for a total of 4 measurements. The result is the
        average of the 4 values. The DICOM data is already HU-corrected so
        -1000 => 0. This means we will search for 300 HU (-1000 + 700) and 800 HU (-1000 + 200) respectively.

        This cuts the profile in half, searches for the highest-gradient index (where the phantom edge is),
        then further cuts it down to +/-10 pixels. The 300/800 HU are then found from linear interpolation.
        It was found that artifacts in the image could drastically influence these values, so hence the +/-10
        subset.

        Assumptions:
        -The phantom does not cross the halfway point of the image FOV (i.e. not offset by an obscene amount).
        -10 pixels about the phantom edge is adequate to capture the full dropoff.
        -300 and 800 HU values will be in the profile"""
        dists = {"Top": np.nan, "Bottom": np.nan, "Left": np.nan, "Right": np.nan}
        edge_5mm = int(5 / self.mm_per_pixel)  # physical 5mm distance
        keys = (key for key in dists)
        for array in (self.horiz_array, self.vert_array):
            split_idx = len(array) // 2  # we need not be exact, just close
            left_data, right_data = array[:split_idx], array[split_idx:][::-1]
            for profile_data in (left_data, right_data):
                # find the phantom edge and chop about it
                edge_idx = np.argmax(np.diff(profile_data))
                edge_data = profile_data[edge_idx - edge_5mm : edge_idx + edge_5mm]
                interp_func = interp1d(edge_data, np.arange(len(edge_data)))
                idx_300, idx_800 = interp_func([300, 800])
                dists[next(keys)] = abs(idx_800 - idx_300) * self.mm_per_pixel
        return dists

    def mean_high_contrast_resolution(self) -> float:
        """Mean high-contrast resolution"""
        return float(np.mean(list(self.high_contrast_resolutions().values())))


@capture_warnings
class QuartDVT(CatPhanBase, ResultsDataMixin[QuartDVTResult]):
    """A class for loading and analyzing CT DICOM files of a Quart phantom that comes with the Halcyon.
    Analyzes: HU Uniformity, Image Scaling & HU Linearity.
    """

    _demo_url = "quart.zip"
    _model = "Quart DVT"
    hu_origin_slice_variance = 300
    catphan_radius_mm = 80
    hu_module: QuartHUModule
    hu_module_class = QuartHUModule
    uniformity_module: QuartUniformityModule
    uniformity_module_class = QuartUniformityModule
    geometry_module: QuartGeometryModule
    geometry_module_class = QuartGeometryModule

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
        x_adjustment: float = 0,
        y_adjustment: float = 0,
        angle_adjustment: float = 0,
        roi_size_factor: float = 1,
        scaling_factor: float = 1,
        origin_slice: int | None = None,
        roll_slice_offset: float = -8,
    ):
        """Single-method full analysis of Quart DICOM files.

        Parameters
        ----------
        hu_tolerance : int
            The HU tolerance value for both HU uniformity and linearity.
        scaling_tolerance : float, int
            The scaling tolerance in mm of the geometric nodes on the HU linearity slice (CTP404 module).
        thickness_tolerance : float, int
            The tolerance of the thickness calculation in mm, based on the wire ramps in the CTP404 module.

            .. warning:: Thickness accuracy degrades with image noise; i.e. low mAs images are less accurate.

        cnr_threshold : float, int
            The threshold for "detecting" low-contrast image. See RTD for calculation info.

            .. deprecated:: 3.0

                Use visibility parameter instead.

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
            The slice number of the HU linearity slice. If None, will be automatically determined. This is a
            fallback method in case the automatic method fails.
        roll_slice_offset : float
            The offset in mm from ``origin_slice`` used to select the slice
            for phantom roll detection. The phantom roll is determined based on the two
            inserts in the central vertical axis of the HU module, but this
            detection can be influenced by the inserts used for slice thickness.
            Adjusting this offset to select a different slice can enhance the
            accuracy of phantom roll detection.
        """
        self.x_adjustment = x_adjustment
        self.y_adjustment = y_adjustment
        self.angle_adjustment = angle_adjustment
        self.roi_size_factor = roi_size_factor
        self.scaling_factor = scaling_factor
        self.roll_slice_offset = roll_slice_offset
        self.localize(origin_slice=origin_slice)
        self.hu_module = self.hu_module_class(
            self,
            offset=0,
            hu_tolerance=hu_tolerance,
            thickness_tolerance=thickness_tolerance,
            scaling_tolerance=scaling_tolerance,
        )
        self.uniformity_module = self.uniformity_module_class(
            self, offset=UNIFORMITY_OFFSET_MM, tolerance=hu_tolerance
        )
        self.geometry_module = self.geometry_module_class(
            self, tolerance=3, offset=GEOMETRY_OFFSET_MM
        )

    def plotly_analyzed_images(
        self,
        show: bool = True,
        show_legend: bool = True,
        show_colorbar: bool = True,
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
        figs[self.hu_module.common_name] = self.hu_module.plotly(
            show_colorbar=show_colorbar, show_legend=show_legend, **kwargs
        )
        figs["HU Linearity plot"] = self.hu_module.plotly_linearity(
            show_legend=show_legend
        )
        figs[self.uniformity_module.common_name] = self.uniformity_module.plotly(
            show_colorbar=show_colorbar, show_legend=show_legend, **kwargs
        )
        figs[self.geometry_module.common_name] = self.geometry_module.plotly(
            show_colorbar=show_colorbar, show_legend=show_legend, **kwargs
        )
        figs["Side View"] = self.plotly_side_view(show_legend=show_legend)

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
            f"High-Contrast distance (mm): {self.geometry_module.mean_high_contrast_resolution():2.3f}",
        )
        if as_str:
            return "\n".join(items)
        else:
            return items

    def _generate_results_data(self) -> QuartDVTResult:
        """Return results in a data structure for more programmatic use."""
        return QuartDVTResult(
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
                high_contrast_distances=self.geometry_module.high_contrast_resolutions(),
                mean_high_contrast_distance=self.geometry_module.mean_high_contrast_resolution(),
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

    def _module_offsets(self) -> list[float]:
        absolute_origin_position = self.dicom_stack[self.origin_slice].z_position
        relative_offsets_mm = [0, UNIFORMITY_OFFSET_MM, GEOMETRY_OFFSET_MM]
        return [
            absolute_origin_position + offset_mm for offset_mm in relative_offsets_mm
        ]

    def _detected_modules(self) -> list[CatPhanModule]:
        return [self.uniformity_module, self.hu_module, self.geometry_module]


@capture_warnings
class HypersightQuartDVT(QuartDVT):
    """A class for loading and analyzing CT DICOM files of a Quart phantom that comes with the Halcyon, specifically
    for the Hypersight version, which includes a water ROI.
    Analyzes: HU Uniformity, Image Scaling & HU Linearity.
    """

    def __init__(self, **kwargs):
        warnings.warn(
            "This class is now deprecated. Please use the QuartDVT class instead as it now handles the water vial that differentiated this class",
            DeprecationWarning,
        )
        super().__init__(**kwargs)

    _model = "Hypersight Quart DVT"
    hu_module = HypersightQuartHUModule
    hu_module_class = HypersightQuartHUModule
