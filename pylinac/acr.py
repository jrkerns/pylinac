from __future__ import annotations

import io
import math
import textwrap
import warnings
import webbrowser
from io import BytesIO
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from plotly import graph_objects as go
from pydantic import BaseModel, ConfigDict, Field
from scipy import ndimage
from skimage.filters import threshold_li, threshold_otsu

from .core import pdf
from .core.array_utils import fill_middle_zeros, find_nearest_idx
from .core.geometry import Line, LineSerialized, Point
from .core.image import DicomImage
from .core.mtf import MTF
from .core.plotly_utils import add_title
from .core.profile import FWXMProfile
from .core.roi import HighContrastDiskROI, RectangleROI
from .core.utilities import QuaacDatum, ResultBase, ResultsDataMixin
from .core.warnings import capture_warnings
from .ct import (
    CatPhanBase,
    CatPhanModule,
    Slice,
    ThicknessROI,
    get_regions,
    rois_to_results,
)

# CT
CT_UNIFORMITY_MODULE_OFFSET_MM = 70
CT_SPATIAL_RESOLUTION_MODULE_OFFSET_MM = 100
CT_LOW_CONTRAST_MODULE_OFFSET_MM = 30

# MR
MR_SLICE11_MODULE_OFFSET_MM = 100
MR_GEOMETRIC_DISTORTION_MODULE_OFFSET_MM = 40
MR_UNIFORMITY_MODULE_OFFSET_MM = 60


class CTModule(CatPhanModule):
    common_name = "HU Linearity"
    attr_name = "ct_calibration_module"
    roi_dist_mm = 63
    roi_radius_mm = 10
    roi_settings = {
        "Air": {"angle": 45, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "Poly": {"angle": 225, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "Acrylic": {"angle": 135, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "Bone": {"angle": -45, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "Water": {"angle": 180, "distance": roi_dist_mm, "radius": roi_radius_mm},
    }
    window_min = -200
    window_max = 200


class CTModuleOutput(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.

    Use the following attributes as normal class attributes."""

    offset: float = Field(
        description="The offset of the module slice in mm from the origin slice (z-direction)."
    )
    roi_distance_from_center_mm: float = Field(
        description="The distance of the ROIs from the center of the phantom in mm in the image plane."
    )
    roi_radius_mm: float = Field(description="The radius of the ROIs in mm.")
    roi_settings: dict = Field(
        description="The ROI settings. The keys are the material names."
    )
    rois: dict[str, float] = Field(
        description="The analyzed ROIs. The key is the name of the material and the value is the mean HU value. E.g. ``'Air': -987.1``."
    )


class UniformityModule(CatPhanModule):
    """Class for analysis of the Uniformity slice of the CTP module. Measures 5 ROIs around the slice that
    should all be close to the same value.
    """

    attr_name = "uniformity_module"
    common_name = "HU Uniformity"
    roi_dist_mm = 66
    roi_radius_mm = 11
    roi_settings = {
        "Top": {"angle": -90, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "Right": {"angle": 0, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "Bottom": {"angle": 90, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "Left": {"angle": 180, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "Center": {"angle": 0, "distance": 0, "radius": roi_radius_mm},
    }
    window_min = -50
    window_max = 50


class UniformityModuleOutput(CTModuleOutput):
    """This class should not be called directly. It is returned by the ``results_data()`` method.

    Use the following attributes as normal class attributes."""

    center_roi_stdev: float = Field(
        description="The standard deviation of the center ROI.",
        title="Center ROI Standard Deviation",
    )


class SpatialResolutionModule(CatPhanModule):
    """Class for analysis of the Uniformity slice of the CTP module. Measures 5 ROIs around the slice that
    should all be close to the same value.
    """

    attr_name = "spatial_resolution_module"
    common_name = "Spatial Resolution"
    rois: dict[str, HighContrastDiskROI]
    roi_dist_mm = 70
    roi_radius_mm = 6
    roi_settings = {
        "10oclock": {
            "angle": -135,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
            "lp/mm": 0.4,
        },
        "9oclock": {
            "angle": -180,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
            "lp/mm": 0.5,
        },
        "7oclock": {
            "angle": 135,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
            "lp/mm": 0.6,
        },
        "6oclock": {
            "angle": 90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
            "lp/mm": 0.7,
        },
        "4oclock": {
            "angle": 45,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
            "lp/mm": 0.8,
        },
        "3oclock": {
            "angle": 0,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
            "lp/mm": 0.9,
        },
        "2oclock": {
            "angle": -45,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
            "lp/mm": 1.0,
        },
        "12oclock": {
            "angle": -90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
            "lp/mm": 1.2,
        },
    }

    def _setup_rois(self) -> None:
        for name, setting in self.roi_settings.items():
            self.rois[name] = HighContrastDiskROI.from_phantom_center(
                self.image,
                setting["angle_corrected"],
                setting["radius_pixels"],
                setting["distance_pixels"],
                self.phan_center,
                contrast_threshold=1.0,  # fixed to 1 so everything passes. We aren't evaluating pass/fail here
            )

    @property
    def mtf(self) -> MTF:
        spacings = [roi["lp/mm"] for roi in self.roi_settings.values()]
        return MTF.from_high_contrast_diskset(
            spacings=spacings, diskset=list(self.rois.values())
        )

    def plotly_rois(self, fig: go.Figure) -> None:
        for name, roi in self.rois.items():
            roi.plotly(fig, line_color="green", name=name)

    def plot_rois(self, axis: plt.Axes) -> None:
        """Plot the ROIs to the axis. Override to set the color"""
        for roi, mtf in zip(self.rois.values(), self.mtf.norm_mtfs.values()):
            roi.plot2axes(axis, edgecolor="g")


class SpatialResolutionModuleOutput(CTModuleOutput):
    """This class should not be called directly. It is returned by the ``results_data()`` method.

    Use the following attributes as normal class attributes."""

    lpmm_to_rmtf: dict = Field(
        description="Line pair to relative modulation transfer mapping. The keys are the line pair values and the values are the relative modulation transfer values.",
        title="Line Pair to Relative MTF",
    )


class LowContrastModule(CatPhanModule):
    """Class for analysis of the Uniformity slice of the CTP module. Measures 5 ROIs around the slice that
    should all be close to the same value.
    """

    attr_name = "low_contrast_module"
    common_name = "Low Contrast"
    roi_dist_mm = 60
    roi_radius_mm = 6
    nominal_value = 0
    roi_settings = {
        "ROI": {"angle": -90, "distance": roi_dist_mm, "radius": roi_radius_mm},
    }
    background_roi_settings = {
        "ROI": {"angle": -115, "distance": roi_dist_mm, "radius": roi_radius_mm},
    }
    window_min = 50
    window_max = 150

    def cnr(self) -> float:
        """Given in the guidance doc as |A-B|/SD where A is the contrast ROI, B is the background, and SD is stdev of B"""
        return (
            abs(self.rois["ROI"].pixel_value - self.background_rois["ROI"].pixel_value)
            / self.background_rois["ROI"].std
        )


class LowContrastModuleOutput(CTModuleOutput):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    cnr: float = Field(
        description="The contrast-to-noise ratio.", title="Contrast to Noise Ratio"
    )


class ACRCTResult(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.

    Use the following attributes as normal class attributes."""

    phantom_model: str = Field(description="The model of the phantom used.")
    phantom_roll_deg: float = Field(
        description="The roll of the phantom in degrees.",
        title="Phantom roll (\N{DEGREE SIGN})",
    )
    origin_slice: int = Field(
        description="The slice number of the 'origin' slice; for ACR this is Module 1."
    )
    num_images: int = Field(description="The number of images in the passed dataset.")
    ct_module: CTModuleOutput = Field(
        description="The results of the CT module.", title="CT Module"
    )
    uniformity_module: UniformityModuleOutput = Field(
        description="The results of the Uniformity module.",
        title="HU Uniformity",
    )
    low_contrast_module: LowContrastModuleOutput = Field(
        description="The results of the Low Contrast module.",
        title="Low Contrast Resolution",
    )
    spatial_resolution_module: SpatialResolutionModuleOutput = Field(
        description="The results of the Spatial Resolution module.",
        title="Spatial Resolution",
    )


@capture_warnings
class ACRCT(CatPhanBase, ResultsDataMixin[ACRCTResult]):
    _model = "ACR CT 464"
    catphan_radius_mm = 100
    air_bubble_radius_mm = 14
    min_num_images = 4
    localization_radius = 70
    ct_calibration_module = CTModule
    low_contrast_module = LowContrastModule
    spatial_resolution_module = SpatialResolutionModule
    uniformity_module = UniformityModule
    clear_borders = False

    def _detected_modules(self) -> list[CatPhanModule]:
        return [
            self.ct_calibration_module,
            self.low_contrast_module,
            self.spatial_resolution_module,
            self.uniformity_module,
        ]

    def plot_analyzed_subimage(self, *args, **kwargs):
        raise NotImplementedError("Use `plot_images`")

    def save_analyzed_subimage(self, *args, **kwargs):
        raise NotImplementedError("Use `save_images`")

    def analyze(
        self,
        x_adjustment: float = 0,
        y_adjustment: float = 0,
        angle_adjustment: float = 0,
        roi_size_factor: float = 1,
        scaling_factor: float = 1,
        origin_slice: int | None = None,
    ) -> None:
        """Analyze the ACR CT phantom

        Parameters
        ----------
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
        origin_slice: int, None
            The slice number of the HU linearity module. If None, will be automatically determined. This is a fallback
            method in case the automatic determination fails.
        """
        self.x_adjustment = x_adjustment
        self.y_adjustment = y_adjustment
        self.angle_adjustment = angle_adjustment
        self.roi_size_factor = roi_size_factor
        self.scaling_factor = scaling_factor
        self.localize(origin_slice=origin_slice)
        self.ct_calibration_module = self.ct_calibration_module(
            self, offset=0, clear_borders=self.clear_borders
        )
        self.uniformity_module = self.uniformity_module(
            self,
            offset=CT_UNIFORMITY_MODULE_OFFSET_MM,
            clear_borders=self.clear_borders,
        )
        self.spatial_resolution_module = self.spatial_resolution_module(
            self,
            offset=CT_SPATIAL_RESOLUTION_MODULE_OFFSET_MM,
            clear_borders=self.clear_borders,
        )
        self.low_contrast_module = self.low_contrast_module(
            self,
            offset=CT_LOW_CONTRAST_MODULE_OFFSET_MM,
            clear_borders=self.clear_borders,
        )

    def plotly_analyzed_images(
        self,
        show: bool = True,
        show_colorbar: bool = True,
        show_legend: bool = True,
        **kwargs,
    ) -> dict[str, go.Figure]:
        """Plot the analyzed image using Plotly. Will create multiple figures.

        Parameters
        ----------
        show : bool
            Whether to show the images. Set to False if doing further processing of the figure.
        show_colorbar : bool
            Whether to show the colorbar on the images.
        show_legend : bool
            Whether to show the legend on the images.
        kwargs
            Additional keyword arguments to pass to the figure.
        """
        figs = {}
        for module in (
            self.ct_calibration_module,
            self.uniformity_module,
            self.spatial_resolution_module,
            self.low_contrast_module,
        ):
            figs[module.common_name] = module.plotly(
                show_colorbar=show_colorbar, show_legend=show_legend, **kwargs
            )
        figs["MTF"] = self.spatial_resolution_module.mtf.plotly(
            show_legend=show_legend, **kwargs
        )
        figs["Side View"] = self.plotly_side_view(show_legend=show_legend)

        if show:
            for fig in figs.values():
                fig.show()
        return figs

    def plot_analyzed_image(self, show: bool = True, **plt_kwargs) -> plt.Figure:
        """Plot the analyzed image

        Parameters
        ----------
        show
            Whether to show the image.
        plt_kwargs
            Keywords to pass to matplotlib for figure customization.
        """
        # set up grid and axes
        fig = plt.figure(**plt_kwargs)
        grid_size = (2, 3)
        hu_ax = plt.subplot2grid(grid_size, (0, 0))
        self.ct_calibration_module.plot(hu_ax)
        unif_ax = plt.subplot2grid(grid_size, (0, 1))
        self.uniformity_module.plot(unif_ax)
        sr_ax = plt.subplot2grid(grid_size, (0, 2))
        self.spatial_resolution_module.plot(sr_ax)
        locon_ax = plt.subplot2grid(grid_size, (1, 0))
        self.low_contrast_module.plot(locon_ax)
        spatial_res_graph = plt.subplot2grid(grid_size, (1, 2))
        self.spatial_resolution_module.mtf.plot(spatial_res_graph)
        side_ax = plt.subplot2grid(grid_size, (1, 1))
        self.plot_side_view(side_ax)

        # finish up
        plt.tight_layout()
        if show:
            plt.show()
        return fig

    def save_analyzed_image(self, filename: str | Path | BytesIO, **plt_kwargs) -> None:
        """Save the analyzed image to disk or stream

        Parameters
        ----------
        filename
            Where to save the image to
        plt_kwargs
            Keywords to pass to matplotlib for figure customization.
        """
        fig = self.plot_analyzed_image(show=False, **plt_kwargs)
        fig.savefig(filename)

    def plot_images(self, show: bool = True, **plt_kwargs) -> dict[str, plt.Figure]:
        """Plot all the individual images separately

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
            "hu": self.ct_calibration_module,
            "uniformity": self.uniformity_module,
            "spatial resolution": self.spatial_resolution_module,
            "low contrast": self.low_contrast_module,
        }
        for key, module in modules.items():
            fig, ax = plt.subplots(**plt_kwargs)
            module.plot(ax)
            figs[key] = fig
        # plot the one-off MTF image
        fig, ax = plt.subplots(**plt_kwargs)
        figs["mtf"] = fig
        self.spatial_resolution_module.mtf.plot(ax)
        # plot the side view
        fig, ax = plt.subplots(**plt_kwargs)
        figs["side"] = fig
        self.plot_side_view(ax)

        plt.tight_layout()

        if show:
            plt.show()
        return figs

    def save_images(
        self,
        directory: Path | str | None = None,
        to_stream: bool = False,
        **plt_kwargs,
    ) -> list[Path | BytesIO]:
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
        for name, fig in figs.items():
            if to_stream:
                path = io.BytesIO()
            else:
                destination = Path(directory) or Path.cwd()
                path = (destination / name).with_suffix(".png").absolute()
            fig.savefig(path)
            paths.append(path)
        return paths

    def find_phantom_roll(self, func=lambda roi: roi.bbox_area) -> float:
        """Determine the "roll" of the phantom.

        Only difference of base method is that we sort the ROIs by size,
        not by being in the center since the two we're looking for are both right-sided.
        """
        return super().find_phantom_roll(func)

    def results(self) -> str:
        """Return the results of the analysis as a string. Use with print()."""
        string = (
            f"\n - ACR CT 464 QA Test - \n"
            f"HU ROIs: {self.ct_calibration_module.roi_vals_as_str}\n"
            f"Contrast to Noise Ratio: {self.low_contrast_module.cnr():2.2f}\n"
            f"Uniformity ROIs: {self.uniformity_module.roi_vals_as_str}\n"
            f'Uniformity Center ROI standard deviation: {self.uniformity_module.rois["Center"].std:2.2f}\n'
            f"MTF 50% (lp/mm): {self.spatial_resolution_module.mtf.relative_resolution(50):2.2f}\n"
        )
        return string

    def _generate_results_data(self) -> ACRCTResult:
        return ACRCTResult(
            phantom_model="ACR CT 464",
            phantom_roll_deg=self.catphan_roll,
            origin_slice=self.origin_slice,
            num_images=self.num_images,
            ct_module=CTModuleOutput(
                offset=0,
                roi_distance_from_center_mm=self.ct_calibration_module.roi_dist_mm,
                roi_radius_mm=self.ct_calibration_module.roi_radius_mm,
                roi_settings=self.ct_calibration_module.roi_settings,
                rois={
                    name: roi.pixel_value
                    for name, roi in self.ct_calibration_module.rois.items()
                },
            ),
            uniformity_module=UniformityModuleOutput(
                offset=CT_UNIFORMITY_MODULE_OFFSET_MM,
                roi_distance_from_center_mm=self.uniformity_module.roi_dist_mm,
                roi_radius_mm=self.uniformity_module.roi_radius_mm,
                roi_settings=self.uniformity_module.roi_settings,
                rois={
                    name: roi.pixel_value
                    for name, roi in self.uniformity_module.rois.items()
                },
                center_roi_stdev=self.uniformity_module.rois["Center"].std,
            ),
            spatial_resolution_module=SpatialResolutionModuleOutput(
                offset=CT_SPATIAL_RESOLUTION_MODULE_OFFSET_MM,
                roi_distance_from_center_mm=self.spatial_resolution_module.roi_dist_mm,
                roi_radius_mm=self.spatial_resolution_module.roi_radius_mm,
                roi_settings=self.spatial_resolution_module.roi_settings,
                rois={
                    name: roi.pixel_value
                    for name, roi in self.spatial_resolution_module.rois.items()
                },
                lpmm_to_rmtf=self.spatial_resolution_module.mtf.norm_mtfs,
            ),
            low_contrast_module=LowContrastModuleOutput(
                offset=CT_LOW_CONTRAST_MODULE_OFFSET_MM,
                roi_distance_from_center_mm=self.low_contrast_module.roi_dist_mm,
                roi_radius_mm=self.low_contrast_module.roi_radius_mm,
                roi_settings=self.low_contrast_module.roi_settings,
                rois={
                    name: roi.pixel_value
                    for name, roi in self.low_contrast_module.rois.items()
                },
                cnr=self.low_contrast_module.cnr(),
            ),
        )

    def _quaac_datapoints(self) -> dict[str, QuaacDatum]:
        results_data = self.results_data(as_dict=True)
        data = {}
        data["Phantom Roll"] = QuaacDatum(
            value=results_data["phantom_roll_deg"],
            unit="degrees",
            description="The roll of the phantom in the image",
        )
        for name, value in results_data["ct_module"]["rois"].items():
            data[f"{name} HU"] = QuaacDatum(
                value=value,
                unit="HU",
                description=f"The HU value of the {name} ROI",
            )
        for name, value in results_data["uniformity_module"]["rois"].items():
            data[f"{name} Uniformity HU"] = QuaacDatum(
                value=value,
                unit="HU",
                description=f"The HU value of the {name} Uniformity ROI",
            )
        for name, value in results_data["spatial_resolution_module"][
            "lpmm_to_rmtf"
        ].items():
            data[f"{name} lp/mm"] = QuaacDatum(
                value=value,
                unit="rMTF",
            )
        for name, value in results_data["low_contrast_module"]["rois"].items():
            data[f"{name} CNR"] = QuaacDatum(
                value=value,
                unit="CNR",
                description=f"The CNR value of the {name} ROI",
            )
        return data

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
        texts = [
            " - ACR CT 464 Results - ",
            f"HU Linearity ROIs: {self.ct_calibration_module.roi_vals_as_str}",
            f"Low contrast visibility: {self.low_contrast_module.cnr():2.2f}",
            f"Uniformity ROIs: {self.uniformity_module.roi_vals_as_str}",
        ]
        analysis_images = self.save_images(to_stream=True)

        canvas = pdf.PylinacCanvas(
            filename, page_title=analysis_title, metadata=metadata, logo=logo
        )
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 4.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 4))

        for idx, text in enumerate(texts):
            canvas.add_text(text=text, location=(1.5, 23 - idx * 0.5))
        for page, img in enumerate(analysis_images):
            canvas.add_new_page()
            canvas.add_image(img, location=(1, 5), dimensions=(18, 18))
        canvas.finish()

        if open_file:
            webbrowser.open(filename)

    def _module_offsets(self) -> list[float]:
        absolute_origin_position = self.dicom_stack[self.origin_slice].z_position
        relative_offsets_mm = [
            0,
            CT_UNIFORMITY_MODULE_OFFSET_MM,
            CT_LOW_CONTRAST_MODULE_OFFSET_MM,
            CT_SPATIAL_RESOLUTION_MODULE_OFFSET_MM,
        ]
        return [
            absolute_origin_position + offset_mm for offset_mm in relative_offsets_mm
        ]


class MRSlice11PositionModule(CatPhanModule):
    common_name = "Slice Position, Slice 11"
    roi_settings = {
        "Left": {"width": 2, "height": 25, "distance": 65, "angle": 2.5},
        "Right": {"width": 2, "height": 25, "distance": 65, "angle": -2.5},
    }
    rois: dict = {}

    def _setup_rois(self) -> None:
        for name, setting in self.roi_settings.items():
            # angle is +90 because pointing right is 0, and these rois move downward, not rightward
            self.rois[name] = RectangleROI.from_phantom_center(
                self.image,
                setting["width_pixels"],
                setting["height_pixels"],
                self.catphan_roll - 90 + setting["angle"],
                setting["distance_pixels"],
                self.phan_center,
            )

    @property
    def bar_difference_mm(self) -> float:
        """The difference in height between the two angled bars"""
        idxs = []
        for roi in (self.rois["Right"], self.rois["Left"]):
            prof = roi.pixel_array.max(axis=np.argmin(roi.pixel_array.shape))
            mid_height = (prof.max() - prof.min()) / 2 + prof.min()
            idx = find_nearest_idx(prof, mid_height)
            idxs.append(idx)
        return (idxs[0] - idxs[1]) * self.mm_per_pixel

    @property
    def slice_shift_mm(self) -> float:
        """The effective shift in phantom position in the S/I direction. Because bars are at 45 degrees, the shift is half the bar difference"""
        return self.bar_difference_mm / 2

    def plot_rois(self, axis: plt.Axes) -> None:
        """Plot the ROIs to the axis.

        We overload because simple rectangle ROIs don't have a pass/fail color.
        """
        for roi in self.rois.values():
            roi.plot2axes(axis, edgecolor="blue")

    def plotly_rois(self, fig: go.Figure) -> None:
        for name, roi in self.rois.items():
            roi.plotly(fig, line_color="blue", name=name)


class MRSlice11ModuleOutput(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.

    Use the following attributes as normal class attributes."""

    offset: int = Field(
        description="The offset of the phantom in mm from the origin slice."
    )
    roi_settings: dict = Field(
        description="The ROI settings. The keys are the ROI names."
    )
    rois: dict = Field(
        description="The results of the left and right bar ROIs. The key is the name of the bar"
    )
    bar_difference_mm: float = Field(
        description="The difference in bar positions in mm.",
        title="Bar Difference (mm)",
    )
    slice_shift_mm: float = Field(
        description="The measure shift in slice position compared to nominal.",
        title="Slice Shift (mm)",
    )


class MRSlice1Module(CatPhanModule):
    common_name = "Slice 1 (Thickness, Offset, Resolution)"
    slice_lines: dict[str, Line]
    thickness_rois: dict[str, ThicknessROI] = {}
    thickness_roi_settings = {
        "Top": {"width": 100, "height": 3, "distance": -3},
        "Bottom": {"width": 100, "height": 3, "distance": 2.5},
    }
    roi_settings = {
        "Row Reference": {"radius": 9, "distance": 58, "angle": 135, "lp/mm": 0},
        "Col Reference": {"radius": 9, "distance": 58, "angle": 135, "lp/mm": 0},
        "Row 1.1": {"radius": 3, "distance": 40, "angle": 116, "lp/mm": 1 / 1.1},
        "Col 1.1": {"radius": 3, "distance": 44, "angle": 104, "lp/mm": 1 / 1.1},
        "Row 1.0": {"radius": 3, "distance": 36, "angle": 81, "lp/mm": 1.0},
        "Col 1.0": {"radius": 3, "distance": 44, "angle": 74, "lp/mm": 1.0},
        "Row 0.9": {"radius": 2, "distance": 46, "angle": 52, "lp/mm": 1 / 0.9},
        "Col 0.9": {"radius": 2, "distance": 55, "angle": 51, "lp/mm": 1 / 0.9},
    }
    position_roi_settings = {
        "Left": {"width": 2, "height": 25, "distance": 65, "angle": 2.5},
        "Right": {"width": 2, "height": 25, "distance": 65, "angle": -2.5},
    }
    position_rois: dict = {}
    rois: dict[str, HighContrastDiskROI]
    spacings = [0, 1 / 1.1, 1, 1 / 0.9]

    def _setup_rois(self) -> None:
        # thickness
        for name, setting in self.thickness_roi_settings.items():
            # angle is +90 because pointing right is 0, and these rois move downward, not rightward
            self.thickness_rois[name] = ThicknessROI.from_phantom_center(
                self.image,
                setting["width_pixels"],
                setting["height_pixels"],
                self.catphan_roll + 90,
                setting["distance_pixels"],
                self.phan_center,
            )
        # spatial res
        for name, setting in self.roi_settings.items():
            self.rois[name] = HighContrastDiskROI.from_phantom_center(
                self.image,
                setting["angle_corrected"],
                setting["radius_pixels"],
                setting["distance_pixels"],
                self.phan_center,
                contrast_threshold=1.0,  # fixed to 1 so everything passes. We aren't evaluating pass/fail here
            )
        # slice position
        for name, setting in self.position_roi_settings.items():
            # angle is +90 because pointing right is 0, and these rois move downward, not rightward
            self.position_rois[name] = ThicknessROI.from_phantom_center(
                self.image,
                setting["width_pixels"],
                setting["height_pixels"],
                self.catphan_roll - 90 + setting["angle"],
                setting["distance_pixels"],
                self.phan_center,
            )

    def plot_rois(self, axis: plt.Axes) -> None:
        for roi in self.position_rois.values():
            roi.plot2axes(axis, edgecolor="blue")
        for roi in self.thickness_rois.values():
            roi.plot2axes(axis, edgecolor="blue")
        for roi, mtf in zip(self.rois.values(), self.rois.values()):
            roi.plot2axes(axis, edgecolor="g")

    def plotly_rois(self, fig: go.Figure) -> None:
        for name, roi in self.position_rois.items():
            roi.plotly(fig, line_color="blue", name=name)
        for name, roi in self.thickness_rois.items():
            roi.plotly(fig, line_color="blue", name=name)
        for name, roi in self.rois.items():
            roi.plotly(fig, line_color="green", name=name)

    @property
    def bar_difference_mm(self) -> float:
        """The difference in height between the two angled bars"""
        left_array = self.position_rois["Left"].long_profile.values
        left_mid_height = (left_array.max() - left_array.min()) / 2 + left_array.min()
        left_idx = find_nearest_idx(left_array, left_mid_height)
        right_array = self.position_rois["Right"].long_profile.values
        right_mid_height = (
            right_array.max() - right_array.min()
        ) / 2 + right_array.min()
        right_idx = find_nearest_idx(right_array, right_mid_height)
        return (right_idx - left_idx) * self.mm_per_pixel

    @property
    def slice_shift_mm(self) -> float:
        """The effective shift in phantom position in the S/I direction. Because bars are at 45 degrees, the shift is half the bar difference"""
        return self.bar_difference_mm / 2

    @property
    def measured_slice_thickness_mm(self) -> float:
        """The slice thickness as determined by the two angled ROIs in the center of Slice 1"""
        top = self.thickness_rois["Top"].wire_fwhm * self.mm_per_pixel
        bottom = self.thickness_rois["Bottom"].wire_fwhm * self.mm_per_pixel
        return 0.2 * (top * bottom) / (top + bottom)

    @property
    def row_mtf(self) -> MTF:
        """The MTF of the spatial resolution module looking at the row-wise ROIs"""
        return MTF.from_high_contrast_diskset(
            spacings=self.spacings,
            diskset=list(roi for name, roi in self.rois.items() if "Row" in name),
        )

    @property
    def col_mtf(self) -> MTF:
        """The MTF of the spatial resolution module looking at the column-wise ROIs"""
        return MTF.from_high_contrast_diskset(
            spacings=self.spacings,
            diskset=list(roi for name, roi in self.rois.items() if "Col" in name),
        )


class MRSlice1ModuleOutput(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.

    Use the following attributes as normal class attributes."""

    offset: int = Field(
        description="The offset of the phantom in mm from the origin slice."
    )
    roi_settings: dict = Field(
        description="A dictionary of the ROI settings. The keys are the ROI names."
    )
    rois: dict = Field(
        description=" A dictionary of the analyzed MTF ROIs. The key is the name of the ROI; e.g. ``Row 1.1``."
    )
    bar_difference_mm: float = Field(
        description="The difference in bar positions in mm.",
        title="Bar Difference (mm)",
    )
    slice_shift_mm: float = Field(
        description="The measured shift in slice position compared to nominal.",
        title="Slice Shift (mm)",
    )
    measured_slice_thickness_mm: float = Field(
        description="The measured slice thickness in mm.",
        title="Measured Slice Thickness (mm)",
    )
    row_mtf_50: float = Field(
        description="The MTF at 50% for the row-based ROIs.",
        title="Row-wise 50% MTF (lp/mm)",
    )
    col_mtf_50: float = Field(
        description="The MTF at 50% for the column-based ROIs.",
        title="Column-wise 50% MTF (lp/mm)",
    )
    row_mtf_lp_mm: dict[int, float] = Field(
        description="A key-value pair of the MTF. The key is the relative resolution in % and the value is the lp/mm at that resolution",
        title="MTF (lp/mm)",
    )
    col_mtf_lp_mm: dict[int, float] = Field(
        description="A key-value pair of the MTF. The key is the relative resolution in % and the value is the lp/mm at that resolution",
        title="MTF (lp/mm)",
    )


class MRUniformityModule(CatPhanModule):
    """Class for analysis of the Uniformity slice of the CTP module. Measures 5 ROIs around the slice that
    should all be close to the same value.
    """

    common_name = "Signal Uniformity"
    roi_settings = {
        "Center": {
            "angle": 90,
            "distance": 5,
            "radius": 80,
        },  # 80 radius ~= 200cm2, per the manual
    }
    ghost_roi_settings = {
        # size of ~900mm2 per the manual
        "Top": {"angle": -90, "distance": 110, "width": 60, "height": 15},
        "Bottom": {"angle": 90, "distance": 110, "width": 60, "height": 15},
        "Left": {"angle": 180, "distance": 110, "width": 15, "height": 60},
        "Right": {"angle": 0, "distance": 110, "width": 15, "height": 60},
    }
    ghost_rois: dict = {}

    def __init__(self, catphan, offset):
        self.tesla = float(catphan.dicom_stack.metadata.MagneticFieldStrength)
        super().__init__(catphan, tolerance=None, offset=offset)

    def _setup_rois(self) -> None:
        super()._setup_rois()
        for name, roi in self.ghost_roi_settings.items():
            self.ghost_rois[name] = RectangleROI.from_phantom_center(
                self.image,
                roi["width_pixels"],
                roi["height_pixels"],
                roi["angle"] + self.catphan_roll,
                roi["distance_pixels"],
                self.phan_center,
            )

    def plot_rois(self, axis: plt.Axes) -> None:
        super().plot_rois(axis)
        for roi in self.ghost_rois.values():
            roi.plot2axes(axis, edgecolor="yellow")

    def plotly_rois(self, fig: go.Figure) -> None:
        super().plotly_rois(fig)
        for name, roi in self.ghost_rois.items():
            roi.plotly(fig, line_color="yellow", name=name)

    @property
    def percent_image_uniformity(self) -> float:
        """PIU value calculated via section 5.3 of the manual"""
        piu_high = np.percentile(self.rois["Center"].pixel_values, 99)
        piu_low = np.percentile(self.rois["Center"].pixel_values, 1)
        return 100 * (1 - ((piu_high - piu_low) / (piu_high + piu_low)))

    @property
    def piu_passed(self) -> bool:
        """Section 5.4"""
        if self.tesla < 3:
            return self.percent_image_uniformity > 85
        else:
            return self.percent_image_uniformity > 80

    @property
    def ghosting_ratio(self) -> float:
        """Ghosting ratio of section 6.3 of the manual"""
        top = self.ghost_rois["Top"].pixel_value
        bottom = self.ghost_rois["Bottom"].pixel_value
        left = self.ghost_rois["Left"].pixel_value
        right = self.ghost_rois["Right"].pixel_value
        return abs(
            ((top + bottom) - (left + right)) / (2 * self.rois["Center"].pixel_value)
        )

    @property
    def psg(self) -> float:
        """Percent Signal Ghosting"""
        return self.ghosting_ratio * 100

    @property
    def psg_passed(self) -> bool:
        """Whether the PSG is within tolerance"""
        return self.psg < 3.0


class MRUniformityModuleOutput(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.

    Use the following attributes as normal class attributes."""

    offset: int = Field(
        description="The offset of the phantom in mm from the origin slice."
    )
    roi_settings: dict = Field(
        description="A dictionary of the ROI settings. The keys are the ROI names."
    )
    rois: dict = Field(description="A dictionary of the analyzed ROIs.")
    ghost_roi_settings: dict = Field(
        description="A dictionary of the ghost ROI settings. The keys are the ROI names."
    )
    ghost_rois: dict = Field(description="A dictionary of the ghost ROIs.")
    psg: float = Field(
        description="The percent signal ghosting.", title="Percent Signal Ghosting"
    )
    ghosting_ratio: float = Field(
        description="The ghosting ratio.", title="Ghosting Ratio"
    )
    piu_passed: bool = Field(
        description="Whether the percent integral uniformity passed the test."
    )
    piu: float = Field(
        description="The percent integral uniformity.",
        title="Percent Integral Uniformity",
    )


class GeometricDistortionModule(CatPhanModule):
    """Class for analysis of the Uniformity slice of the CTP module. Measures 5 ROIs around the slice that
    should all be close to the same value.
    """

    common_name = "Geometric Distortion"
    profiles: dict

    def _setup_rois(self) -> None:
        """This is mostly for plotting purposes. This is why we use FWXMProfile
        instead of FWXMProfilePhysical. The lines to plot should be in pixel coordinates, not physical.
        We convert to physical just for the field width calculation."""
        px_to_cut_off = int(round(5 / self.mm_per_pixel))
        self.profiles = {}
        threshold = threshold_otsu(image=self.image.array)
        bin_image = self.image.as_binary(threshold=threshold)
        bin_image = ndimage.binary_fill_holes(bin_image).astype(float)
        # calculate horizontal
        data = bin_image[int(self.phan_center.y), :]
        # cutoff 3mm from the search area
        f_data = fill_middle_zeros(data, cutoff_px=px_to_cut_off)
        prof = FWXMProfile(values=f_data)
        line = Line(
            Point(prof.field_edge_idx(side="left"), self.phan_center.y),
            Point(prof.field_edge_idx(side="right"), self.phan_center.y),
        )

        self.profiles["horizontal"] = {
            "width (mm)": prof.field_width_px * self.mm_per_pixel,
            "line": line,
        }
        # calculate vertical
        data = bin_image[:, int(self.phan_center.x)]
        f_data = fill_middle_zeros(data, cutoff_px=px_to_cut_off)
        prof = FWXMProfile(values=f_data)
        line = Line(
            Point(self.phan_center.x, prof.field_edge_idx(side="left")),
            Point(self.phan_center.x, prof.field_edge_idx(side="right")),
        )
        self.profiles["vertical"] = {
            "width (mm)": prof.field_width_px * self.mm_per_pixel,
            "line": line,
        }
        # calculate negative diagonal
        # calculate slope equation intercept
        # b = y - (+1)x
        b = self.phan_center.y - self.phan_center.x
        xs = np.arange(0, self.image.shape[1])
        ys = xs + b
        coords = ndimage.map_coordinates(bin_image, [ys, xs], order=1, mode="mirror")
        f_data = fill_middle_zeros(coords, cutoff_px=px_to_cut_off)
        prof = FWXMProfile(values=f_data)
        line = Line(
            Point(
                xs[int(round(prof.field_edge_idx(side="left")))],
                ys[int(round(prof.field_edge_idx(side="left")))],
            ),
            Point(
                xs[int(round(prof.field_edge_idx(side="right")))],
                ys[int(round(prof.field_edge_idx(side="right")))],
            ),
        )
        # pixels are now diagonal and thus spacing between pixels is now the hypotenuse
        # We don't have to fix the line above because that's in pixels. The issue is the
        # geometric distance and thus we only have to correct here.
        self.profiles["negative diagonal"] = {
            "width (mm)": prof.field_width_px * self.mm_per_pixel * math.sqrt(2),
            "line": line,
        }
        # calculate positive diagonal
        # calculate slope equation intercept
        # b = y - (-1)x
        b = self.phan_center.y + self.phan_center.x
        ys = -xs + b
        coords = ndimage.map_coordinates(bin_image, [ys, xs], order=1, mode="mirror")
        f_data = fill_middle_zeros(coords, cutoff_px=px_to_cut_off)
        prof = FWXMProfile(values=f_data)
        line = Line(
            Point(
                xs[int(round(prof.field_edge_idx(side="left")))],
                ys[int(round(prof.field_edge_idx(side="left")))],
            ),
            Point(
                xs[int(round(prof.field_edge_idx(side="right")))],
                ys[int(round(prof.field_edge_idx(side="right")))],
            ),
        )
        self.profiles["positive diagonal"] = {
            "width (mm)": prof.field_width_px * self.mm_per_pixel * math.sqrt(2),
            "line": line,
        }

    def plotly_rois(self, fig: go.Figure) -> None:
        for name, profile_data in self.profiles.items():
            profile_data["line"].plotly(fig, line_width=2, color="blue", name=name)

    def plot_rois(self, axis: plt.Axes):
        for name, profile_data in self.profiles.items():
            profile_data["line"].plot2axes(axis, width=2, color="blue")

    def distances(self) -> dict:
        """The measurements of the phantom size for all 4 lines in mm"""
        return {name: f"{p['width (mm)']:2.2f}mm" for name, p in self.profiles.items()}


class MRGeometricDistortionModuleOutput(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.

    Use the following attributes as normal class attributes."""

    model_config = ConfigDict(arbitrary_types_allowed=True)
    offset: int = Field(
        description="The offset of the phantom in mm from the origin slice."
    )
    profiles: dict[str, dict[str, float | LineSerialized]] = Field(
        description="A dictionary of the profiles used to measure the geometric distortion. The key is the name of the profile.",
        title="Profile widths (mm)",
    )
    distances: dict = Field(
        description="The lines measuring the ROI size. The key is the name of the line direction and the value is a string of the line length.",
        title="Distance measurements (mm)",
    )


class SagittalLocalizationModule:
    """Class for analysis of the Sagittal slice of the ACR Phantom."""

    """This class shares some similarities with CatphanModule but also has key differences.
    The options considered were: inherit to override, refactor into a common abstraction,
    or duplicate plot code for a standalone implementation.
    By design choice, the latter was taken."""

    common_name = "Sagittal Distortion"
    roi_settings: dict[str, dict[str, float]] = {
        "ROI1": {"offset": -75},
        "ROI2": {"offset": -25},
        "ROI3": {"offset": 25},
        "ROI4": {"offset": 75},
    }  # offsets in mm left/right from phantom centroid
    rois: dict[str, Line] = {}
    profiles: dict = {}
    image: DicomImage
    window_min: int | None = None  # plt visualization
    window_max: int | None = None  # plt visualization

    def __init__(self, image: DicomImage | None):
        if image is None:
            return

        self.image = image

        # The processing is consistent with the geometric distortion module
        threshold = round(threshold_li(image.array))
        bin_image = image.as_binary(threshold=threshold)
        bin_image = ndimage.binary_fill_holes(bin_image).astype(float)

        centroid = np.argwhere(bin_image).mean(axis=0)
        pixel_size = 1 / image.dpmm
        for key, val in self.roi_settings.items():
            offset_px = val["offset"] * pixel_size
            col = round(centroid[1] + offset_px)
            data = bin_image[:, col]
            prof = FWXMProfile(values=data)
            line = Line(
                Point(col, prof.field_edge_idx(side="left")),
                Point(col, prof.field_edge_idx(side="right")),
            )
            self.profiles[key] = {
                "width (mm)": prof.field_width_px * pixel_size,
                "line": line,
            }
            self.rois[key] = line

    def distances(self) -> dict:
        """The measurements of the phantom size in mm"""
        return {name: f"{p['width (mm)']:2.2f}mm" for name, p in self.profiles.items()}

    def plot(self, axis: plt.Axes):
        """Plot the image along with ROIs to an axis"""
        self.image.plot(ax=axis, show=False, vmin=self.window_min, vmax=self.window_max)
        self.plot_rois(axis)
        axis.autoscale(tight=True)
        axis.set_title(f"{self.common_name}")
        axis.axis("off")

    def plotly(self, **kwargs) -> go.Figure:
        """Plot the image along with the ROIs to a plotly figure."""
        fig = go.Figure()
        self.image.plotly(
            fig, show=False, zmin=self.window_min, zmax=self.window_max, **kwargs
        )
        self.plotly_rois(fig)
        add_title(fig, f"{self.common_name}")
        return fig

    def plotly_rois(self, fig: go.Figure) -> None:
        for name, profile_data in self.profiles.items():
            profile_data["line"].plotly(fig, line_width=2, color="blue", name=name)

    def plot_rois(self, axis: plt.Axes):
        for name, profile_data in self.profiles.items():
            profile_data["line"].plot2axes(axis, width=2, color="blue")


class MRSagittalLocalizationModuleOutput(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.

    Use the following attributes as normal class attributes."""

    model_config = ConfigDict(arbitrary_types_allowed=True)
    profiles: dict[str, dict[str, float | LineSerialized]] = Field(
        description="A dictionary of the profiles used to measure the geometric distortion. The key is the name of the profile.",
        title="Profile widths (mm)",
    )
    distances: dict = Field(
        description="The lines measuring the ROI size. The key is the name of the line direction and the value is a string of the line length.",
        title="Distance measurements (mm)",
    )


class ACRMRIResult(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.

    Use the following attributes as normal class attributes."""

    phantom_model: str = Field(description="The model of the phantom used.")
    phantom_roll_deg: float = Field(description="The roll of the phantom in degrees.")
    origin_slice: int = Field(
        description="The slice number of the 'origin' slice; for ACR this is Slice 1."
    )
    num_images: int = Field(description="The number of images in the passed dataset.")
    slice1: MRSlice1ModuleOutput = Field(
        description="The results for the 'Slice 1' module", title="Slice 1 Module"
    )
    slice11: MRSlice11ModuleOutput = Field(
        description="The results for the 'Slice 11' module", title="Slice 11 Module"
    )
    uniformity_module: MRUniformityModuleOutput = Field(
        description="Results from the uniformity module", title="Uniformity Module"
    )
    geometric_distortion_module: MRGeometricDistortionModuleOutput = Field(
        description="Results from the geometric distortion module",
        title="Geometric Distortion Module",
    )
    sagittal_localizer_module: MRSagittalLocalizationModuleOutput = Field(
        description="Results from the sagittal localizer module",
        title="Sagittal Localization Module",
    )


@capture_warnings
class ACRMRILarge(CatPhanBase, ResultsDataMixin[ACRMRIResult]):
    _model = "ACR MRI Large"
    catphan_radius_mm = 100
    min_num_images = 4
    air_bubble_radius_mm = 20
    slice1 = MRSlice1Module
    geometric_distortion = GeometricDistortionModule
    uniformity_module = MRUniformityModule
    slice11 = MRSlice11PositionModule
    sagittal_localization = SagittalLocalizationModule
    has_sagittal_module: bool = False
    clip_in_localization = False

    def plot_analyzed_subimage(self, *args, **kwargs):
        raise NotImplementedError("Use `plot_images`")

    def save_analyzed_subimage(self, *args, **kwargs):
        raise NotImplementedError("Use `save_images`")

    def localize(self) -> None:
        self._phantom_center_func = self.find_phantom_axis()
        self.catphan_roll = self.find_phantom_roll() + self.angle_adjustment
        # now that we have the origin slice, ensure we have scanned all linked modules
        if not self._ensure_physical_scan_extent():
            raise ValueError(
                "The physical scan extent does not cover the extent of module configuration. "
                "This means not all modules were included in the scan. Rescan the phantom to include all "
                "relevant modules, or change the offset values."
            )

    def _module_offsets(self) -> list[float]:
        absolute_origin_position = self.dicom_stack[self.origin_slice].z_position
        relative_offsets_mm = [
            0,
            MR_GEOMETRIC_DISTORTION_MODULE_OFFSET_MM,
            MR_UNIFORMITY_MODULE_OFFSET_MM,
            MR_SLICE11_MODULE_OFFSET_MM,
        ]
        return [
            absolute_origin_position + offset_mm for offset_mm in relative_offsets_mm
        ]

    def find_phantom_roll(self) -> float:
        """Determine the "roll" of the phantom. This algorithm uses the circular left-upper hole on slice 1 as the reference

        Returns
        -------
        float : the angle of the phantom in **degrees**.
        """
        # get edges and make ROIs from it
        slice = Slice(self, self.origin_slice)
        larr, regions, _ = get_regions(slice)
        try:
            # find appropriate ROIs and grab the two most centrally positioned ones
            circle_bubbles = [
                r
                for r in regions
                if (self._is_right_area(r) and self._is_right_eccentricity(r))
            ]
            exact_size = np.pi * ((self.air_bubble_radius_mm / self.mm_per_pixel) ** 2)
            most_similar_bubble = sorted(
                circle_bubbles, key=lambda r: abs(r.filled_area - exact_size)
            )[0]
            y_dist = most_similar_bubble.centroid[0] - slice.phan_center.y
            x_dist = most_similar_bubble.centroid[1] - slice.phan_center.x
            phan_roll = np.arctan2(y_dist, x_dist)
            corrected_roll = (
                np.rad2deg(phan_roll) + 135
            )  # bubble is at top-left. perfect placement is -135
            return corrected_roll
        except Exception:
            raise RuntimeError(
                "Could not determine the roll of the phantom. Ensure the 20mm top-left circle is visible on Slice 1"
            )

    def analyze(
        self,
        echo_number: int | None = None,
        x_adjustment: float = 0,
        y_adjustment: float = 0,
        angle_adjustment: float = 0,
        roi_size_factor: float = 1,
        scaling_factor: float = 1,
    ) -> None:
        """Analyze the ACR CT phantom

        Parameters
        ----------
        echo_number:
            The echo to analyze. If not passed, uses the minimum echo number found.
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
        """
        self.x_adjustment = x_adjustment
        self.y_adjustment = y_adjustment
        self.angle_adjustment = angle_adjustment
        self.roi_size_factor = roi_size_factor
        self.scaling_factor = scaling_factor
        self._select_echo_images(echo_number)
        sagittal_image = self._select_sagittal_image()
        self.has_sagittal_module = sagittal_image is not None
        self.localize()
        self.slice1 = self.slice1(self, offset=0)
        self.geometric_distortion = self.geometric_distortion(
            self, offset=MR_GEOMETRIC_DISTORTION_MODULE_OFFSET_MM
        )
        self.uniformity_module = self.uniformity_module(
            self, offset=MR_UNIFORMITY_MODULE_OFFSET_MM
        )
        self.slice11 = self.slice11(self, offset=MR_SLICE11_MODULE_OFFSET_MM)
        self.sagittal_localization = self.sagittal_localization(sagittal_image)

    def _select_echo_images(self, echo_number: int | None) -> None:
        """Select out the images that match the given echo number"""
        # we check for multiple echos. We only pick the first echo found.
        # this is probably not the best logic but we somehow have to pick
        # Echo Numbers is an int; https://dicom.innolitics.com/ciods/mr-image/mr-image/00180086

        # in case EchoNumbers isn't there, use all
        try:
            all_echos = {int(i.metadata.EchoNumbers) for i in self.dicom_stack}
        except AttributeError:
            # no manipulation; use all images
            return
        if echo_number is None:
            echo_number = min(all_echos)
            if len(all_echos) > 1:
                warnings.warn(
                    f"Multiple echoes found ({all_echos}) and no echo number was passed. Using echo # {echo_number}"
                )
        if echo_number not in all_echos:
            raise ValueError(
                f"Echo number {echo_number} was passed but not found in the dataset. Found echo numbers: {all_echos}. Remove the echo_number parameter or pick a valid echo number."
            )
        # drop images that don't have the same echo number
        to_pop = []
        for idx, img in enumerate([i for i in self.dicom_stack].copy()):
            if int(img.metadata.EchoNumbers) != echo_number:
                to_pop.append(idx)
        for idx in sorted(to_pop, reverse=True):
            del self.dicom_stack[idx]
            del self.dicom_stack.metadatas[idx]

    def _select_sagittal_image(self, max_dist: float = 0.01) -> DicomImage | None:
        """This function return the sagittal image from the dicom stack.

        Parameters
        ----------
        max_dist: float
            Max vectorial distance from the nominal to the actual image orientation.
        """
        # It uses the table in this article to find the sagittal image
        # https://www.kaggle.com/code/rickandjoe/mri-orientation-axial-sagittal-or-coronal

        nominal_sagittal_image_orientation = np.array([0, 1, 0, 0, 0, -1])
        metadatas = self.dicom_stack.metadatas
        image_orientation = [m.ImageOrientationPatient for m in metadatas]
        diff = np.array(image_orientation) - nominal_sagittal_image_orientation
        dist = np.linalg.norm(diff, axis=1)
        if np.sum(dist < max_dist) > 1:
            raise ValueError("There are too many sagittal images in the dataset.")
        min_value = dist.min()
        min_index = dist.argmin()
        if min_value >= max_dist:
            return None

        image = self.dicom_stack[min_index]
        # Remove from the stack since localize() assumes only axial slices are present
        del self.dicom_stack[min_index]
        del self.dicom_stack.metadatas[min_index]
        return image

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
        # plot the images
        modules: list = [
            self.slice1,
            self.geometric_distortion,
            self.uniformity_module,
            self.slice11,
        ]
        if self.has_sagittal_module:
            modules.append(self.sagittal_localization)
        for module in modules:
            figs[module.common_name] = module.plotly(
                show_colorbar=show_colorbar, show_legend=show_legend, **kwargs
            )
        # side view
        figs["Side View"] = self.plotly_side_view(show_legend=show_legend)
        # mtf
        fig = go.Figure()
        self.slice1.row_mtf.plotly(fig=fig, name="Row-wise rMTF")
        figs["MTF"] = self.slice1.col_mtf.plotly(
            show_legend=show_legend,
            fig=fig,
            name="Column-wise rMTF",
            marker_color="orange",
        )

        if show:
            for fig in figs.values():
                fig.show()
        return figs

    def plot_analyzed_image(self, show: bool = True, **plt_kwargs) -> plt.Figure:
        """Plot the analyzed image

        Parameters
        ----------
        show
            Whether to show the image.
        plt_kwargs
            Keywords to pass to matplotlib for figure customization.
        """
        # set up grid and axes
        fig = plt.figure(**plt_kwargs)
        if self.has_sagittal_module:
            grid_size = (2, 4)
        else:
            grid_size = (2, 3)
        slice1_ax = plt.subplot2grid(grid_size, (0, 0))
        self.slice1.plot(slice1_ax)
        geom_ax = plt.subplot2grid(grid_size, (0, 1))
        self.geometric_distortion.plot(geom_ax)
        unif_ax = plt.subplot2grid(grid_size, (0, 2))
        self.uniformity_module.plot(unif_ax)
        position_ax = plt.subplot2grid(grid_size, (1, 0))
        self.slice11.plot(position_ax)

        side_view_ax = plt.subplot2grid(grid_size, (1, 1))
        self.plot_side_view(side_view_ax)
        spatial_res_graph = plt.subplot2grid(grid_size, (1, 2))
        self.slice1.row_mtf.plot(spatial_res_graph, label="Row-wise rMTF")
        self.slice1.col_mtf.plot(spatial_res_graph, label="Column-wise rMTF")
        spatial_res_graph.legend()

        if self.has_sagittal_module:
            sag_ax = plt.subplot2grid(grid_size, (0, 3))
            self.sagittal_localization.plot(sag_ax)

        # finish up
        plt.tight_layout()
        if show:
            plt.show()
        return fig

    def plot_images(self, show: bool = True, **plt_kwargs) -> dict[str, plt.Figure]:
        """Plot all the individual images separately

        Parameters
        ----------
        show
            Whether to show the images.
        plt_kwargs
            Keywords to pass to matplotlib for figure customization.
        """
        figs = {}
        # plot the images
        modules: dict = {
            "geometric": self.geometric_distortion,
            "slice 1": self.slice1,
            "signal uniformity": self.uniformity_module,
            "slice 11": self.slice11,
        }
        if self.has_sagittal_module:
            modules["sagittal"] = self.sagittal_localization
        for key, module in modules.items():
            fig, ax = plt.subplots(**plt_kwargs)
            module.plot(ax)
            figs[key] = fig
        # plot rMTF
        fig, ax = plt.subplots(**plt_kwargs)
        self.slice1.row_mtf.plot(ax, label="Row-wise rMTF")
        self.slice1.col_mtf.plot(ax, label="Column-wise rMTF")
        ax.legend()
        figs["rMTF"] = fig
        # plot the side view
        fig, ax = plt.subplots(**plt_kwargs)
        figs["side"] = fig
        self.plot_side_view(ax)

        if show:
            plt.show()
        return figs

    def _detected_modules(self) -> list[CatPhanModule]:
        return [
            self.slice1,
            self.slice11,
            self.uniformity_module,
            self.geometric_distortion,
        ]

    def save_images(
        self,
        directory: Path | str | None = None,
        to_stream: bool = False,
        **plt_kwargs,
    ) -> list[Path | BytesIO]:
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
        for name, fig in figs.items():
            if to_stream:
                path = io.BytesIO()
            else:
                destination = Path(directory) or Path.cwd()
                path = (destination / name).with_suffix(".png").absolute()
            fig.savefig(path)
            paths.append(path)
        return paths

    def _quaac_datapoints(self) -> dict[str, QuaacDatum]:
        results_data = self.results_data(as_dict=True)
        data = {}

        data["Phantom Roll"] = QuaacDatum(
            value=results_data["phantom_roll_deg"],
            unit="degrees",
            description="The roll of the phantom in the image",
        )
        slice1_keys = (
            ("bar_difference_mm", "Bar Difference", "mm"),
            ("slice_shift_mm", "Slice Shift", "mm"),
            ("measured_slice_thickness_mm", "Measured Slice Thickness", "mm"),
            ("row_mtf_50", "Row-wise MTF 50%", "lp/mm"),
            ("col_mtf_50", "Column-wise MTF 50%", "lp/mm"),
        )
        for key, name, unit in slice1_keys:
            data[name] = QuaacDatum(
                value=results_data["slice1"][key],
                unit=unit,
            )
        for name, roi in results_data["slice11"]["rois"].items():
            data[f"Slice 11 {name} ROI"] = QuaacDatum(
                value=roi["value"],
                unit="HU",
            )
        data["Slice 11 Bar Difference"] = QuaacDatum(
            value=results_data["slice11"]["bar_difference_mm"],
            unit="mm",
        )
        data["Slice 11 Slice Shift"] = QuaacDatum(
            value=results_data["slice11"]["slice_shift_mm"],
            unit="mm",
        )
        for name, roi in results_data["uniformity_module"]["rois"].items():
            data[f"Uniformity {name} ROI"] = QuaacDatum(
                value=roi["value"],
                unit="HU",
            )
        for name, roi in results_data["uniformity_module"]["ghost_rois"].items():
            data[f"Uniformity {name} Ghost ROI"] = QuaacDatum(
                value=roi["value"],
                unit="HU",
            )
        data["Percent Signal Ghosting"] = QuaacDatum(
            value=results_data["uniformity_module"]["psg"],
            unit="%",
        )
        data["Ghosting Ratio"] = QuaacDatum(
            value=results_data["uniformity_module"]["ghosting_ratio"],
            unit="",
        )
        data["Percent Integral Uniformity"] = QuaacDatum(
            value=results_data["uniformity_module"]["piu"],
            unit="%",
        )
        for name, line in results_data["geometric_distortion_module"][
            "profiles"
        ].items():
            data[f"Geometric Distortion {name} line length"] = QuaacDatum(
                value=line["width (mm)"],
                unit="mm",
            )
        for name, line in results_data["sagittal_localizer_module"]["profiles"].items():
            data[f"Localizer Distortion {name} line length"] = QuaacDatum(
                value=line["width (mm)"],
                unit="mm",
            )
        return data

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
            filename, page_title=analysis_title, metadata=metadata
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
        for page, img in enumerate(analysis_images):
            canvas.add_new_page()
            canvas.add_image(img, location=(1, 5), dimensions=(18, 18))
        canvas.finish()

        if open_file:
            webbrowser.open(filename)

    def results(self, as_str: bool = True) -> str | tuple:
        """Return the results of the analysis as a string. Use with print()."""
        string = (
            f" - {self._model} Results - ",
            f"Geometric Distortions: {self.geometric_distortion.distances()}",
            f"Slice Thickness: {self.slice1.measured_slice_thickness_mm:2.2f}mm",
            f"Slice 1 S/I Position shift: {self.slice1.slice_shift_mm:2.2f}mm",
            f"Slice 11 S/I Position shift: {self.slice11.slice_shift_mm:2.2f}mm",
            f"Uniformity PIU: {self.uniformity_module.percent_image_uniformity:2.2f}",
            f"Percent-signal ghosting: {self.uniformity_module.psg:2.2f}%",
            f'Uniformity Center ROI standard deviation: {self.uniformity_module.rois["Center"].std:2.2f}',
            f"Row-wise MTF 50% (lp/mm): {self.slice1.row_mtf.relative_resolution(50):2.2f}",
            f"Column-wise MTF 50% (lp/mm): {self.slice1.col_mtf.relative_resolution(50):2.2f}",
            f"Sagittal Distortions: {self.sagittal_localization.distances()}",
        )
        if as_str:
            return "\n".join(string)
        else:
            return string

    def _generate_results_data(self) -> ACRMRIResult:
        resolutions = range(10, 91, 10)  # 10-90% in 10% increments
        row_mtfs = {
            resolution: self.slice1.row_mtf.relative_resolution(resolution)
            for resolution in resolutions
        }
        col_mtfs = {
            resolution: self.slice1.col_mtf.relative_resolution(resolution)
            for resolution in resolutions
        }
        return ACRMRIResult(
            phantom_model=self._model,
            phantom_roll_deg=self.catphan_roll,
            origin_slice=self.origin_slice,
            num_images=self.num_images,
            slice1=MRSlice1ModuleOutput(
                offset=0,
                roi_settings=self.slice1.roi_settings,
                rois=rois_to_results(self.slice1.rois),
                bar_difference_mm=self.slice1.bar_difference_mm,
                slice_shift_mm=self.slice1.slice_shift_mm,
                measured_slice_thickness_mm=self.slice1.measured_slice_thickness_mm,
                row_mtf_50=self.slice1.row_mtf.relative_resolution(50),
                col_mtf_50=self.slice1.col_mtf.relative_resolution(50),
                row_mtf_lp_mm=row_mtfs,
                col_mtf_lp_mm=col_mtfs,
            ),
            slice11=MRSlice11ModuleOutput(
                offset=MR_SLICE11_MODULE_OFFSET_MM,
                bar_difference_mm=self.slice11.bar_difference_mm,
                slice_shift_mm=self.slice11.slice_shift_mm,
                rois=rois_to_results(self.slice11.rois),
                roi_settings=self.slice11.roi_settings,
            ),
            geometric_distortion_module=MRGeometricDistortionModuleOutput(
                offset=MR_GEOMETRIC_DISTORTION_MODULE_OFFSET_MM,
                profiles=self.geometric_distortion.profiles,
                distances=self.geometric_distortion.distances(),
            ),
            uniformity_module=MRUniformityModuleOutput(
                offset=0,
                roi_settings=self.uniformity_module.roi_settings,
                rois=rois_to_results(self.uniformity_module.rois),
                ghost_roi_settings=self.uniformity_module.ghost_roi_settings,
                ghost_rois=rois_to_results(self.uniformity_module.ghost_rois),
                psg=self.uniformity_module.psg,
                ghosting_ratio=self.uniformity_module.ghosting_ratio,
                piu=self.uniformity_module.percent_image_uniformity,
                piu_passed=self.uniformity_module.piu_passed,
            ),
            sagittal_localizer_module=MRSagittalLocalizationModuleOutput(
                profiles=self.sagittal_localization.profiles,
                distances=self.sagittal_localization.distances(),
            ),
        )
