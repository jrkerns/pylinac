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
from pydantic import BaseModel
from scipy import ndimage

from .core import pdf
from .core.array_utils import find_nearest_idx
from .core.geometry import Line, Point
from .core.mtf import MTF
from .core.profile import FWXMProfilePhysical
from .core.roi import HighContrastDiskROI, RectangleROI
from .core.utilities import ResultBase, ResultsDataMixin
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
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    offset: int
    roi_distance_from_center_mm: int
    roi_radius_mm: int
    roi_settings: dict
    rois: dict


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
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    center_roi_stdev: float  #:


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
            self.rois[name] = HighContrastDiskROI(
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

    def plot_rois(self, axis: plt.Axes) -> None:
        """Plot the ROIs to the axis. Override to set the color"""
        for roi, mtf in zip(self.rois.values(), self.mtf.norm_mtfs.values()):
            roi.plot2axes(axis, edgecolor="g")


class SpatialResolutionModuleOutput(CTModuleOutput):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    lpmm_to_rmtf: dict  #:


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

    cnr: float  #:


class ACRCTResult(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    phantom_model: str  #:
    phantom_roll_deg: float  #:
    origin_slice: int  #:
    num_images: int  #:
    ct_module: CTModuleOutput  #:
    uniformity_module: UniformityModuleOutput  #:
    low_contrast_module: LowContrastModuleOutput  #:
    spatial_resolution_module: SpatialResolutionModuleOutput  #:


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

    def analyze(self) -> None:
        """Analyze the ACR CT phantom"""
        self.localize()
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
            self.rois[name] = RectangleROI(
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


class MRSlice11ModuleOutput(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    offset: int  #:
    roi_settings: dict  #:
    rois: dict  #:
    bar_difference_mm: float  #:
    slice_shift_mm: float  #:


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
            self.thickness_rois[name] = ThicknessROI(
                self.image,
                setting["width_pixels"],
                setting["height_pixels"],
                self.catphan_roll + 90,
                setting["distance_pixels"],
                self.phan_center,
            )
        # spatial res
        for name, setting in self.roi_settings.items():
            self.rois[name] = HighContrastDiskROI(
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
            self.position_rois[name] = ThicknessROI(
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
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    offset: int  #:
    roi_settings: dict  #:
    rois: dict  #:
    bar_difference_mm: float  #:
    slice_shift_mm: float  #:
    measured_slice_thickness_mm: float  #:
    row_mtf_50: float  #:
    col_mtf_50: float  #:


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
            self.ghost_rois[name] = RectangleROI(
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
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    offset: int  #:
    roi_settings: dict  #:
    rois: dict  #:
    ghost_roi_settings: dict  #:
    ghost_rois: dict  #:
    psg: float  #:
    ghosting_ratio: float  #:
    piu_passed: bool  #:
    piu: float  #:


class GeometricDistortionModule(CatPhanModule):
    """Class for analysis of the Uniformity slice of the CTP module. Measures 5 ROIs around the slice that
    should all be close to the same value.
    """

    common_name = "Geometric Distortion"
    profiles: dict

    def _setup_rois(self) -> None:
        self.profiles = {}
        bin_image = self.image.as_binary(threshold=np.percentile(self.image, 60))
        bin_image = ndimage.binary_fill_holes(bin_image).astype(float)
        # calculate horizontal
        data = bin_image[int(self.phan_center.y), :]
        prof = FWXMProfilePhysical(values=data, dpmm=1 / self.mm_per_pixel)
        line = Line(
            Point(prof.field_edge_idx(side="left"), self.phan_center.y),
            Point(prof.field_edge_idx(side="right"), self.phan_center.y),
        )

        self.profiles["horizontal"] = {
            "width (mm)": prof.field_width_mm,
            "line": line,
        }
        # calculate vertical
        data = bin_image[:, int(self.phan_center.x)]
        prof = FWXMProfilePhysical(values=data, dpmm=1 / self.mm_per_pixel)
        line = Line(
            Point(self.phan_center.x, prof.field_edge_idx(side="left")),
            Point(self.phan_center.x, prof.field_edge_idx(side="right")),
        )
        self.profiles["vertical"] = {
            "width (mm)": prof.field_width_mm,
            "line": line,
        }
        # calculate negative diagonal
        # calculate slope equation intercept
        # b = y - (+1)x
        b = self.phan_center.y - self.phan_center.x
        xs = np.arange(0, self.image.shape[1])
        ys = xs + b
        coords = ndimage.map_coordinates(bin_image, [ys, xs], order=1, mode="mirror")
        # pixels are now diagonal and thus spacing between pixels is now the hypotenuse
        prof = FWXMProfilePhysical(
            values=coords, dpmm=1 / (self.mm_per_pixel * math.sqrt(2))
        )
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
        self.profiles["negative diagonal"] = {
            "width (mm)": prof.field_width_mm,
            "line": line,
        }
        # calculate positive diagonal
        # calculate slope equation intercept
        # b = y - (-1)x
        b = self.phan_center.y + self.phan_center.x
        ys = -xs + b
        coords = ndimage.map_coordinates(bin_image, [ys, xs], order=1, mode="mirror")
        prof = FWXMProfilePhysical(
            values=coords, dpmm=1 / (self.mm_per_pixel * math.sqrt(2))
        )
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
            "width (mm)": prof.field_width_mm,
            "line": line,
        }

    def plot_rois(self, axis: plt.Axes):
        for name, profile_data in self.profiles.items():
            profile_data["line"].plot2axes(axis, width=2, color="blue")

    def distances(self) -> dict:
        """The measurements of the phantom size for all 4 lines in mm"""
        return {name: f"{p['width (mm)']:2.2f}mm" for name, p in self.profiles.items()}


class MRGeometricDistortionModuleOutput(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    offset: int  #:
    profiles: dict  #:
    distances: dict  #:


class ACRMRIResult(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    phantom_model: str  #:
    phantom_roll_deg: float  #:
    origin_slice: int  #:
    num_images: int  #:
    slice1: MRSlice1ModuleOutput  #:
    slice11: MRSlice11ModuleOutput  #:
    uniformity_module: MRUniformityModuleOutput  #:
    geometric_distortion_module: MRGeometricDistortionModuleOutput  #:


class ACRMRILarge(CatPhanBase, ResultsDataMixin[ACRMRIResult]):
    _model = "ACR MRI Large"
    catphan_radius_mm = 100
    min_num_images = 4
    air_bubble_radius_mm = 20
    slice1 = MRSlice1Module
    geometric_distortion = GeometricDistortionModule
    uniformity_module = MRUniformityModule
    slice11 = MRSlice11PositionModule

    def plot_analyzed_subimage(self, *args, **kwargs):
        raise NotImplementedError("Use `plot_images`")

    def save_analyzed_subimage(self, *args, **kwargs):
        raise NotImplementedError("Use `save_images`")

    def localize(self) -> None:
        self._phantom_center_func = self.find_phantom_axis()
        self.catphan_roll = self.find_phantom_roll()
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

    def analyze(self, echo_number: int | None = None) -> None:
        """Analyze the ACR CT phantom

        Parameters
        ----------
        echo_number:
            The echo to analyze. If not passed, uses the minimum echo number found.
        """
        self._select_echo_images(echo_number)
        self.localize()
        self.slice1 = self.slice1(self, offset=0)
        self.geometric_distortion = self.geometric_distortion(
            self, offset=MR_GEOMETRIC_DISTORTION_MODULE_OFFSET_MM
        )
        self.uniformity_module = self.uniformity_module(
            self, offset=MR_UNIFORMITY_MODULE_OFFSET_MM
        )
        self.slice11 = self.slice11(self, offset=MR_SLICE11_MODULE_OFFSET_MM)

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
        modules = {
            "geometric": self.geometric_distortion,
            "slice 1": self.slice1,
            "signal uniformity": self.uniformity_module,
            "slice 11": self.slice11,
        }
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
        )
        if as_str:
            return "\n".join(string)
        else:
            return string

    def _generate_results_data(self) -> ACRMRIResult:
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
        )
