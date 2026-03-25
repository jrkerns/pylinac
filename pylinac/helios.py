from __future__ import annotations

import io
import textwrap
import webbrowser
from collections.abc import Callable
from io import BytesIO
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from plotly import graph_objects as go
from pydantic import BaseModel, Field
from skimage import draw

from .core import pdf
from .core.geometry import Point
from .core.mtf import MTF
from .core.roi import RectangleROI
from .core.utilities import QuaacDatum, ResultBase, ResultsDataMixin
from .core.warnings import capture_warnings
from .ct import CatPhanBase, CatPhanModule, Slice

SECTION_3_OFFSET_MM = 60
HELIOS_LOW_CONTRAST_SLICE_OFFSETS_INDEX = {
    "slice_1": 0,
    "slice_2": -1,
    "slice_3": -2,
}


class HeliosContrastScaleModule(CatPhanModule):
    """Class for analysis of the Contrast Scale."""

    common_name = "Contrast Scale"
    attr_name = "contrast_scale_module"
    roi_settings = {
        "Plexiglass": {"width": 10, "height": 10, "distance": 35, "angle": -135},
        "Water": {"width": 10, "height": 10, "distance": 75, "angle": -90},
    }

    def _setup_rois(self) -> None:
        self.rois = {}
        for name, setting in self.roi_settings.items():
            self.rois[name] = RectangleROI.from_phantom_center(
                array=self.image,
                width=setting["width_pixels"],
                height=setting["height_pixels"],
                angle=setting["angle_corrected"],
                dist_from_center=setting["distance_pixels"],
                phantom_center=self.phan_center,
            )

    @property
    def contrast_difference(self) -> float:
        """Difference in mean HU: Plexiglass - Water."""
        return self.rois["Plexiglass"].mean - self.rois["Water"].mean

    def as_dict(self) -> dict:
        """Dump important data as a dictionary."""
        return {
            "data": {
                "mean_hu": {name: roi.mean for name, roi in self.rois.items()},
                "std": {name: roi.std for name, roi in self.rois.items()},
            }
        }

    def plot_rois(self, axis: plt.Axes) -> None:
        """Plot the ROIs to the axis."""
        for roi in self.rois.values():
            roi.plot2axes(axis, edgecolor="blue")

    def plotly_rois(self, fig: go.Figure) -> None:
        """Plot the ROIs to a Plotly figure."""
        for name, roi in self.rois.items():
            roi.plotly(fig, line_color="blue", name=name)


class HeliosContrastScaleModuleOutput(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.

    Use the following attributes as normal class attributes."""

    offset: float = Field(
        description="The offset of this module slice from the origin slice in mm."
    )
    roi_settings: dict = Field(
        description="The ROI settings. The keys are the material names."
    )
    rois: dict = Field(description="The analyzed ROIs.")


class HeliosHighContrastModule(CatPhanModule):
    """Class for analysis of the High Contrast Spatial Resolution."""

    common_name = "High Contrast"
    attr_name = "high_contrast_module"
    roi_settings = {
        "1.6mm": {
            "width": 8,
            "height": 8,
            "distance": 42,
            "angle": -53,
            "bar_size": 1.6,
        },
        "1.3mm": {
            "width": 7,
            "height": 7,
            "distance": 21,
            "angle": -62,
            "bar_size": 1.3,
        },
        "1.0mm": {
            "width": 6,
            "height": 6,
            "distance": 5,
            "angle": -120,
            "bar_size": 1.0,
        },
        "0.8mm": {
            "width": 5,
            "height": 5,
            "distance": 16,
            "angle": 146,
            "bar_size": 0.8,
        },
    }
    rois: dict[str, RectangleROI]

    def _setup_rois(self) -> None:
        self.rois = {}
        for name, setting in self.roi_settings.items():
            self.rois[name] = RectangleROI.from_phantom_center(
                array=self.image,
                width=setting["width_pixels"],
                height=setting["height_pixels"],
                angle=setting["angle_corrected"],
                dist_from_center=setting["distance_pixels"],
                phantom_center=self.phan_center,
            )

    @property
    def mtf(self) -> MTF:
        spacings = [1 / (2 * roi["bar_size"]) for roi in self.roi_settings.values()]
        diskset = list(self.rois.values())
        return MTF.from_high_contrast_diskset(spacings=spacings, diskset=diskset)

    def as_dict(self) -> dict:
        """Dump important data as a dictionary."""
        return {name: roi.std for name, roi in self.rois.items()}

    def plot_rois(self, axis: plt.Axes) -> None:
        """Plot the ROIs to the axis."""
        for roi in self.rois.values():
            roi.plot2axes(axis, edgecolor="blue")

    def plotly_rois(self, fig: go.Figure) -> None:
        """Plot the ROIs to a Plotly figure."""
        for name, roi in self.rois.items():
            roi.plotly(fig, line_color="blue", name=name)


class HeliosHighContrastModuleOutput(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.

    Use the following attributes as normal class attributes."""

    offset: float = Field(
        description="The offset of this module slice from the origin slice in mm."
    )
    rois: dict = Field(description="The analyzed ROIs.")
    mtf_lp_mm: dict[int, float] = Field(
        description="A key-value pair of the MTF. The key is the relative resolution in % and the value is the lp/mm at that resolution",
        title="MTF (lp/mm)",
    )


class HeliosLowContrastModule(CatPhanModule):
    """Class for analysis of the Low Contrast Detectability."""

    common_name = "Low Contrast Detectability"
    attr_name = "low_contrast_module"
    cell_size: float = 5.0
    num_cells: int = 15

    def _setup_rois(self) -> None:
        self.common_name = f"Low Contrast - {self.slice_num + 1}"

        roi_size_px = self.cell_size / self.mm_per_pixel
        total_size_px = roi_size_px * self.num_cells
        half_grid = total_size_px / 2
        half_roi = roi_size_px / 2

        self.rois: list[RectangleROI] = []
        for row in range(self.num_cells):
            for col in range(self.num_cells):
                center = Point(
                    self.phan_center.x - half_grid + col * roi_size_px + half_roi,
                    self.phan_center.y - half_grid + row * roi_size_px + half_roi,
                )
                self.rois.append(
                    RectangleROI(
                        array=self.image,
                        width=roi_size_px,
                        height=roi_size_px,
                        center=center,
                    )
                )

    @property
    def mean(self) -> float:
        """The mean value of the ROIs."""
        return float(np.mean([roi.mean for roi in self.rois]))

    @property
    def std(self) -> float:
        """The std value of the ROIs."""
        return float(np.std([roi.mean for roi in self.rois]))

    def plot_rois(self, axis: plt.Axes) -> None:
        """Plot the grid of ROIs on the axis."""
        for roi in self.rois:
            roi.plot2axes(axis, edgecolor="orange")

    def plotly_rois(self, fig: go.Figure) -> None:
        """Plot the grid of ROIs on a Plotly figure."""
        for roi in self.rois:
            roi.plotly(fig, line_color="orange")


class HeliosLowContrastModuleOutput(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.

    Use the following attributes as normal class attributes."""

    offset: float = Field(
        description="The offset of this module slice from the origin slice in mm."
    )
    settings: dict = Field()
    mean: float = Field(description="Mean HU values of the ROIs")
    std: float = Field(description="Standard deviation of the ROIs.")


class HeliosLowContrastMultiSliceModule:
    """Class for managing Low Contrast analysis across multiple slices."""

    roi_settings = {
        "slice_1": {"offset": HELIOS_LOW_CONTRAST_SLICE_OFFSETS_INDEX["slice_1"]},
        "slice_2": {"offset": HELIOS_LOW_CONTRAST_SLICE_OFFSETS_INDEX["slice_2"]},
        "slice_3": {"offset": HELIOS_LOW_CONTRAST_SLICE_OFFSETS_INDEX["slice_3"]},
    }

    def __init__(self, catphan) -> None:
        """Initialize the multi-slice Low Contrast Detectability module.

        Parameters
        ----------
        catphan : CatPhanBase
            The parent phantom instance.
        """
        self.slices: dict[str, HeliosLowContrastModule] = {}
        slice_spacing = catphan.dicom_stack.slice_spacing
        for key, value in self.roi_settings.items():
            self.slices[key] = HeliosLowContrastModule(
                catphan,
                offset=int(value["offset"] * slice_spacing + SECTION_3_OFFSET_MM),
            )

    @property
    def mean(self) -> float:
        """The mean HU value across all slices."""
        return float(np.mean([s.mean for s in self.slices.values()]))

    @property
    def std(self) -> float:
        """The average standard deviation across all slices."""
        return float(np.mean([s.std for s in self.slices.values()]))


class HeliosLowContrastMultiSliceModuleOutput(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.

    Use the following attributes as normal class attributes."""

    slices: dict[str, HeliosLowContrastModuleOutput] = Field(
        description="Per-slice low contrast results keyed by slice name."
    )
    mean: float = Field(
        description="Mean HU value across all slices.",
    )
    std: float = Field(
        description="Average standard deviation across all slices.",
    )


class HeliosNoiseUniformityModule(CatPhanModule):
    """Class for analysis of the Noise & Uniformity."""

    common_name = "Noise & Uniformity"
    attr_name = "noise_uniformity_module"
    roi_settings = {
        "Center": {"width": 15, "height": 15, "distance": 0, "angle": 0},
        "12 o'clock": {"width": 15, "height": 15, "distance": 75, "angle": -90},
        "3 o'clock": {"width": 15, "height": 15, "distance": 75, "angle": 0},
    }
    noise_roi_settings = {
        "Center": {"width": 25, "height": 25, "distance": 0, "angle": 0},
    }
    rois: dict
    noise_rois: dict

    def _setup_rois(self) -> None:
        self.rois = {}
        self.noise_rois = {}
        for name, setting in self.roi_settings.items():
            self.rois[name] = RectangleROI.from_phantom_center(
                array=self.image,
                width=setting["width_pixels"],
                height=setting["height_pixels"],
                angle=setting["angle_corrected"],
                dist_from_center=setting["distance_pixels"],
                phantom_center=self.phan_center,
            )
        for name, setting in self.noise_roi_settings.items():
            self.noise_rois[name] = RectangleROI.from_phantom_center(
                array=self.image,
                width=setting["width_pixels"],
                height=setting["height_pixels"],
                angle=setting["angle_corrected"],
                dist_from_center=setting["distance_pixels"],
                phantom_center=self.phan_center,
            )

    @property
    def noise_center_std(self) -> float:
        """Std of the central ROI."""
        return self.noise_rois["Center"].std

    @property
    def mean_outer(self) -> float:
        """Mean HU of the outer ROIs."""
        roi_keys = ["12 o'clock", "3 o'clock"]
        return float(np.mean([self.rois[key].mean for key in roi_keys]))

    @property
    def uniformity_difference(self) -> float:
        """Difference between the center mean and the average edge mean."""
        return float(self.rois["Center"].mean - self.mean_outer)

    def as_dict(self) -> dict:
        """Dump important data as a dictionary."""
        return {
            "mean_hu": {name: roi.mean for name, roi in self.rois.items()},
            "std": {name: roi.std for name, roi in self.rois.items()},
        }

    def plot_rois(self, axis: plt.Axes) -> None:
        """Plot the ROIs to the axis."""
        for roi in self.rois.values():
            roi.plot2axes(axis, edgecolor="blue")
        for roi in self.noise_rois.values():
            roi.plot2axes(axis, edgecolor="blue")

    def plotly_rois(self, fig: go.Figure) -> None:
        """Plot the ROIs to a Plotly figure."""
        for name, roi in self.rois.items():
            roi.plotly(fig, line_color="blue", name=name)
        for name, roi in self.noise_rois.items():
            roi.plotly(fig, line_color="blue", name=name)


class HeliosNoiseUniformityModuleOutput(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.

    Use the following attributes as normal class attributes."""

    offset: float = Field(
        description="The offset of this module slice from the origin slice in mm."
    )
    roi_settings: dict = Field(
        description="The ROI settings. The keys are the ROI locations."
    )
    rois: dict = Field(description="The analyzed ROIs.")
    noise_center_std: float = Field(description="The noise in the central ROI")
    mean_outer: float = Field(description="Mean HU values of the outer ROIs.")
    means_diff: float = Field(
        description="Difference between the center ROI mean and the average of the edge ROIs.",
        title="Uniformity Difference (HU)",
    )


class GEHeliosResult(ResultBase):
    """Top-level results for the GE Helios CT Daily phantom analysis.

    This class should not be called directly. It is returned by the
    ``results_data()`` method.
    """

    phantom_model: str = Field(description="The phantom model name.")
    phantom_roll_deg: float = Field(
        description="The roll of the phantom in degrees.",
        title="Phantom Roll (deg)",
    )
    origin_slice: int = Field(
        description="The slice index of the origin (Section 1) slice."
    )
    num_images: int = Field(description="The number of images in the dataset.")
    contrast_scale: HeliosContrastScaleModuleOutput = Field(
        description="The results of the Contrast Scale test."
    )
    high_contrast: HeliosHighContrastModuleOutput = Field(
        description="The results of the High Contrast Spatial Resolution test."
    )
    low_contrast: HeliosLowContrastMultiSliceModuleOutput = Field(
        description="The results of the Low Contrast Detectability multi-slice test."
    )
    noise_uniformity: HeliosNoiseUniformityModuleOutput = Field(
        description="The results of the Noise & Uniformity test."
    )


@capture_warnings
class GEHeliosCTDaily(CatPhanBase, ResultsDataMixin[GEHeliosResult]):
    _model = "GE Helios CT Daily"
    catphan_radius_mm = 107.5
    min_num_images = 8
    clear_borders = False

    contrast_scale_module = HeliosContrastScaleModule
    high_contrast_module = HeliosHighContrastModule
    low_contrast_multi_slice = HeliosLowContrastMultiSliceModule
    noise_uniformity_module = HeliosNoiseUniformityModule

    @classmethod
    def from_demo_image(cls):
        raise NotImplementedError("There is no demo file for this analysis")

    def plot_analyzed_subimage(self, *args, **kwargs):
        raise NotImplementedError("Use `plot_images`")

    def save_analyzed_subimage(self, *args, **kwargs):
        raise NotImplementedError("Use `save_images`")

    def analyze(
        self,
        x_adjustment: float | int = 0,
        y_adjustment: float | int = 0,
        angle_adjustment: float | int = 0,
        roi_size_factor: float | int = 1,
        scaling_factor: float | int = 1,
        origin_slice: int | None = None,
    ) -> None:
        """Analyze the GE Helios CT Daily phantom.

        Parameters
        ----------
        x_adjustment : float
            Shift the phantom center in the x-direction (pixels).
        y_adjustment : float
            Shift the phantom center in the y-direction (pixels).
        angle_adjustment : float
            Add a rotational offset to the phantom roll (degrees).
        roi_size_factor : float
            Factor to scale ROI sizes.  1.0 = nominal.
        scaling_factor : float
            Factor to scale ROI distances from center.  1.0 = nominal.
        origin_slice : int or None
            If given, use this slice index as the origin instead of
            auto-detecting it.
        """
        self.x_adjustment = x_adjustment
        self.y_adjustment = y_adjustment
        self.angle_adjustment = angle_adjustment
        self.roi_size_factor = roi_size_factor
        self.scaling_factor = scaling_factor
        self.localize(origin_slice=origin_slice)
        self.contrast_scale_module = self.contrast_scale_module(
            self, offset=0, clear_borders=self.clear_borders
        )
        self.high_contrast_module = self.high_contrast_module(
            self, offset=0, clear_borders=self.clear_borders
        )
        self.low_contrast_multi_slice = self.low_contrast_multi_slice(self)
        self.noise_uniformity_module = self.noise_uniformity_module(
            self, offset=SECTION_3_OFFSET_MM, clear_borders=self.clear_borders
        )

    def localize(self, origin_slice: int | None = None) -> None:
        """Locate the phantom in the dataset.

        This finds the phantom axis across all slices, identifies the
        Section 1 (Plexiglass block) origin slice, and sets the roll
        angle to zero (the phantom is bracket-mounted).

        Parameters
        ----------
        origin_slice : int or None
            If given, skip auto-detection and use this index directly.
        """
        self._phantom_center_func = self.find_phantom_axis()
        if origin_slice is not None:
            self.origin_slice = origin_slice
        else:
            self.origin_slice = self.find_origin_slice()
        self.catphan_roll = self.find_phantom_roll() + self.angle_adjustment
        if not self._ensure_physical_scan_extent():
            raise ValueError(
                "The physical scan extent does not cover the extent of "
                "module configuration. This means not all modules were "
                "included in the scan. Rescan the phantom to include all "
                "relevant modules, or change the offset values."
            )

    def find_origin_slice(self) -> int:
        """Find the Section 1 slice containing the Plexiglass block.

        The algorithm searches all slices for the one with the highest
        pixel-value variance within the phantom boundary. Section 1
        contains a Plexiglass block (~120 HU) embedded in water (~0 HU)
        which produces distinctly higher variance than the uniform water
        of Section 3.

        Returns
        -------
        int
            The slice index of the origin (Section 1) slice.
        """
        num_slices = len(self.dicom_stack)
        variances = np.zeros(num_slices)
        for idx in range(num_slices):
            slice_obj = Slice(
                self,
                slice_num=idx,
                combine=False,
                clear_borders=self.clear_borders,
            )
            if not slice_obj.is_phantom_in_view():
                continue
            center = slice_obj.phan_center
            radius_px = self.catphan_radius_mm * 0.8 / self.mm_per_pixel
            arr = np.asarray(slice_obj.image)
            rr, cc = draw.disk(
                center=(center.y, center.x),
                radius=radius_px,
                shape=arr.shape,
            )
            pixels = arr[rr, cc]
            variances[idx] = float(np.var(pixels))

        threshold = variances > variances.max() / 2
        best_slice = int(np.mean(np.argwhere(threshold)))
        return best_slice

    def find_phantom_roll(self, func: Callable | None = None) -> float:
        """Return the phantom roll angle."""
        return 0.0

    def _module_offsets(self) -> list[float]:
        """Return the absolute z-positions of the two active sections."""
        absolute_origin_position = self.dicom_stack[self.origin_slice].z_position
        return [
            absolute_origin_position,
            absolute_origin_position + SECTION_3_OFFSET_MM,
        ]

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
            self.contrast_scale_module,
            self.high_contrast_module,
            self.noise_uniformity_module,
        ]
        modules.extend(list(self.low_contrast_multi_slice.slices.values()))
        for module in modules:
            figs[module.common_name] = module.plotly(
                show_colorbar=show_colorbar, show_legend=show_legend, **kwargs
            )
        # side view
        figs["Side View"] = self.plotly_side_view(show_legend=show_legend)
        # mtf
        fig = go.Figure()
        figs["MTF"] = self.high_contrast_module.mtf.plotly(
            show_legend=show_legend, fig=fig, name="rMTF"
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
        modules: list = [
            self.contrast_scale_module,
            self.high_contrast_module,
            self.noise_uniformity_module,
        ]
        modules.extend(list(self.low_contrast_multi_slice.slices.values()))

        # set up grid and axes
        fig, axs = plt.subplots(2, 4, **plt_kwargs)
        axes = axs.ravel()
        ax_idx = -1
        for module in modules:
            ax_idx += 1
            module.plot(axes[ax_idx])

        ax_idx += 1
        self.plot_side_view(axes[ax_idx])
        ax_idx += 1
        self.high_contrast_module.mtf.plot(axes[ax_idx], label="rMTF")
        axes[ax_idx].legend()

        for i in range(ax_idx + 1, len(axes)):
            axes[i].set_visible(False)

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
            "contrast scale": self.contrast_scale_module,
            "high contrast": self.high_contrast_module,
            "noise and uniformity": self.noise_uniformity_module,
        }
        modules |= self.low_contrast_multi_slice.slices
        for key, module in modules.items():
            fig, ax = plt.subplots(**plt_kwargs)
            module.plot(ax)
            figs[key] = fig
        # plot rMTF
        fig, ax = plt.subplots(**plt_kwargs)
        self.high_contrast_module.mtf.plot(ax, label="rMTF")
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
        modules = [
            self.contrast_scale_module,
            self.high_contrast_module,
            self.noise_uniformity_module,
        ]
        modules.extend(list(self.low_contrast_multi_slice.slices.values()))
        return modules

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
                destination = Path(directory) if directory is not None else Path.cwd()
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
        for name, roi in results_data["contrast_scale"]["rois"]["data"][
            "mean_hu"
        ].items():
            data[f"Contrast scale {name} ROI"] = QuaacDatum(
                value=roi,
                unit="HU",
            )
        for resolution, lp_mm in results_data["high_contrast"]["mtf_lp_mm"].items():
            data[f"High contrast MTF {resolution}%"] = QuaacDatum(
                value=lp_mm,
                unit="lp/mm",
            )
        data["Low contrast Mean"] = QuaacDatum(
            value=results_data["low_contrast"]["mean"],
            unit="HU",
        )
        data["Low contrast Std"] = QuaacDatum(
            value=results_data["low_contrast"]["std"],
            unit="HU",
        )
        data["Noise Std"] = QuaacDatum(
            value=results_data["noise_uniformity"]["noise_center_std"],
            unit="HU",
        )
        data["Uniformity Difference"] = QuaacDatum(
            value=results_data["noise_uniformity"]["means_diff"],
            unit="HU",
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
            filename, page_title=analysis_title, metadata=metadata, logo=logo
        )
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 4.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 4))

        shortened_texts = [
            textwrap.wrap(r, width=110) for r in self.results(as_str=False)
        ]
        idx = 0
        for wrapped_lines in shortened_texts:
            for text in wrapped_lines:
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
            f"Contrast Difference: {self.contrast_scale_module.contrast_difference}",
            f"MTF 50% (lp/mm): {self.high_contrast_module.mtf.relative_resolution(50):2.2f}",
            f"Low Contrast Mean: {self.low_contrast_multi_slice.mean:2.2f}",
            f"Low Contrast Standard Deviation: {self.low_contrast_multi_slice.std:2.2f}",
            f"Noise Std: {self.noise_uniformity_module.noise_center_std:2.2f}",
            f"Uniformity Difference: {self.noise_uniformity_module.uniformity_difference:2.2f}",
        )
        if as_str:
            return "\n".join(string)
        else:
            return string

    def _generate_results_data(self) -> GEHeliosResult:
        resolutions = range(10, 91, 10)  # 10-90% in 10% increments
        mtfs = {
            resolution: self.high_contrast_module.mtf.relative_resolution(resolution)
            for resolution in resolutions
        }

        return GEHeliosResult(
            phantom_model=self._model,
            phantom_roll_deg=self.catphan_roll,
            origin_slice=self.origin_slice,
            num_images=self.num_images,
            contrast_scale=HeliosContrastScaleModuleOutput(
                offset=0,
                roi_settings=self.contrast_scale_module.roi_settings,
                rois=self.contrast_scale_module.as_dict(),
            ),
            high_contrast=HeliosHighContrastModuleOutput(
                offset=0,
                rois=self.high_contrast_module.as_dict(),
                mtf_lp_mm=mtfs,
            ),
            low_contrast=HeliosLowContrastMultiSliceModuleOutput(
                slices={
                    k: HeliosLowContrastModuleOutput(
                        offset=self.low_contrast_multi_slice.roi_settings[k]["offset"],
                        settings={
                            "cell_size": v.cell_size,
                            "num_cells": v.num_cells,
                        },
                        mean=v.mean,
                        std=v.std,
                    )
                    for k, v in self.low_contrast_multi_slice.slices.items()
                },
                mean=self.low_contrast_multi_slice.mean,
                std=self.low_contrast_multi_slice.std,
            ),
            noise_uniformity=HeliosNoiseUniformityModuleOutput(
                offset=SECTION_3_OFFSET_MM,
                roi_settings=self.noise_uniformity_module.roi_settings,
                rois=self.noise_uniformity_module.as_dict(),
                noise_center_std=self.noise_uniformity_module.noise_center_std,
                mean_outer=self.noise_uniformity_module.mean_outer,
                means_diff=self.noise_uniformity_module.uniformity_difference,
            ),
        )
