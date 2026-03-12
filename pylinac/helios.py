from __future__ import annotations

import io
import webbrowser
from io import BytesIO
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from plotly import graph_objects as go
from pydantic import BaseModel, Field
from skimage import draw

from .core import pdf
from .core.mtf import MTF
from .core.plotly_utils import add_title
from .core.roi import RectangleROI
from .core.utilities import QuaacDatum, ResultBase, ResultsDataMixin
from .core.warnings import capture_warnings
from .ct import CatPhanBase, CatPhanModule, Slice

SECTION_3_OFFSET_MM = 60


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
                "Std_dev": {name: roi.std for name, roi in self.rois.items()},
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
        for name, setting in self.noise_rois.items():
            self.noise_rois[name] = RectangleROI.from_phantom_center(
                array=self.image,
                width=setting["width_pixels"],
                height=setting["height_pixels"],
                angle=setting["angle_corrected"],
                dist_from_center=setting["distance_pixels"],
                phantom_center=self.phan_center,
            )

    @property
    def noise_center(self) -> float:
        """Std of the central ROI."""
        return self.rois["Center"].std

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
            "data": {
                "mean_hu": {name: roi.mean for name, roi in self.rois.items()},
                "Std_dev": {name: roi.std for name, roi in self.rois.items()},
            }
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
    noise_center: float = Field("The noise in the central ROI")
    mean_outer: float = Field("Mean HU values of the outer ROIs.")
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
    noise_uniformity_module = HeliosNoiseUniformityModule

    def _detected_modules(self) -> list[CatPhanModule]:
        return [
            self.contrast_scale_module,
            self.high_contrast_module,
            self.noise_uniformity_module,
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

    def find_phantom_roll(self, func=None) -> float:
        """Return the phantom roll angle."""
        return 0.0

    def _module_offsets(self) -> list[float]:
        """Return the absolute z-positions of the two active sections."""
        absolute_origin_position = self.dicom_stack[self.origin_slice].z_position
        return [
            absolute_origin_position,
            absolute_origin_position + SECTION_3_OFFSET_MM,
        ]

    # -----------------------------------------------------------------
    # Plotting
    # -----------------------------------------------------------------

    def plot_analyzed_image(self, show: bool = True, **plt_kwargs) -> plt.Figure:
        """Plot all analysis results in a single figure.

        Parameters
        ----------
        show : bool
            Whether to call ``plt.show()``.
        plt_kwargs
            Keyword arguments forwarded to ``plt.figure()``.

        Returns
        -------
        matplotlib.figure.Figure
        """
        fig = plt.figure(**plt_kwargs)
        grid_size = (2, 3)
        cs_ax = plt.subplot2grid(grid_size, (0, 0))
        self.contrast_scale_module.plot(cs_ax)
        hc_ax = plt.subplot2grid(grid_size, (0, 1))
        self.high_contrast_module.plot(hc_ax)
        nu_ax = plt.subplot2grid(grid_size, (0, 2))
        self.noise_uniformity_module.plot(nu_ax)
        mtf_ax = plt.subplot2grid(grid_size, (1, 0))
        self._plot_mtf(mtf_ax)
        side_ax = plt.subplot2grid(grid_size, (1, 1))
        self.plot_side_view(side_ax)
        plt.tight_layout()
        if show:
            plt.show()
        return fig

    def save_analyzed_image(self, filename: str | Path | BytesIO, **plt_kwargs) -> None:
        """Save the combined analysis image to disk or stream.

        Parameters
        ----------
        filename : str, Path, or BytesIO
            Destination file path or stream.
        plt_kwargs
            Keyword arguments forwarded to ``plt.figure()``.
        """
        fig = self.plot_analyzed_image(show=False, **plt_kwargs)
        fig.savefig(filename)

    def _plot_mtf(self, axis: plt.Axes) -> None:
        """Plot the bar-pattern relative MTF curve."""
        mtf = self.high_contrast_module.mtf
        axis.plot(list(mtf.keys()), list(mtf.values()), "bo-")
        axis.set_xlabel("Spatial Frequency (lp/mm)")
        axis.set_ylabel("Relative MTF")
        axis.set_title("MTF")
        axis.set_ylim(0, 1.1)
        axis.grid(True)

    def plotly_analyzed_images(
        self,
        show: bool = True,
        show_colorbar: bool = True,
        show_legend: bool = True,
        **kwargs,
    ) -> dict[str, go.Figure]:
        """Plot the analysis results as Plotly figures.

        Parameters
        ----------
        show : bool
            Whether to call ``fig.show()`` on each figure.
        show_colorbar : bool
            Whether to display the colour-bar on image plots.
        show_legend : bool
            Whether to display the legend.

        Returns
        -------
        dict[str, go.Figure]
            A dictionary mapping descriptive names to Plotly figures.
        """
        figs = {}
        for module in (
            self.contrast_scale_module,
            self.high_contrast_module,
            self.noise_uniformity_module,
        ):
            figs[module.common_name] = module.plotly(
                show_colorbar=show_colorbar,
                show_legend=show_legend,
                **kwargs,
            )
        figs["MTF"] = self.high_contrast_module.mtf.plotly(
            show_legend=show_legend, **kwargs
        )
        if show:
            for fig in figs.values():
                fig.show()
        return figs

    def _plotly_mtf(self, show_legend: bool = True) -> go.Figure:
        """Create a Plotly figure for the MTF curve."""
        fig = go.Figure()
        mtf = self.high_contrast_module.mtf
        fig.add_scatter(
            x=list(mtf.keys()),
            y=list(mtf.values()),
            mode="lines+markers",
            name="MTF",
        )
        add_title(fig, "MTF")
        fig.update_xaxes(title_text="Spatial Frequency (lp/mm)")
        fig.update_yaxes(title_text="Relative MTF", range=[0, 1.1])
        fig.update_layout(showlegend=show_legend)
        return fig

    def plot_images(self, show: bool = True, **plt_kwargs) -> dict[str, plt.Figure]:
        """Plot each analysis image as a separate figure.

        Parameters
        ----------
        show : bool
            Whether to call ``plt.show()``.
        plt_kwargs
            Keyword arguments forwarded to ``plt.subplots()``.

        Returns
        -------
        dict[str, matplotlib.figure.Figure]
        """
        figs = {}
        modules = {
            "contrast_scale": self.contrast_scale_module,
            "high_contrast": self.high_contrast_module,
            "noise_uniformity": self.noise_uniformity_module,
        }
        for key, module in modules.items():
            fig, ax = plt.subplots(**plt_kwargs)
            module.plot(ax)
            figs[key] = fig
        fig, ax = plt.subplots(**plt_kwargs)
        self._plot_mtf(ax)
        figs["mtf"] = fig
        fig, ax = plt.subplots(**plt_kwargs)
        self.plot_side_view(ax)
        figs["side"] = fig
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
        """Save each analysis image to disk or to in-memory streams.

        Parameters
        ----------
        directory : Path or str or None
            Output directory. Defaults to the current working directory.
        to_stream : bool
            If True, return ``BytesIO`` streams instead of writing files.
        plt_kwargs
            Keyword arguments forwarded to ``plt.subplots()``.

        Returns
        -------
        list[Path | BytesIO]
        """
        figs = self.plot_images(show=False, **plt_kwargs)
        paths: list[Path | BytesIO] = []
        for name, fig in figs.items():
            if to_stream:
                path = io.BytesIO()
            else:
                destination = Path(directory) if directory else Path.cwd()
                path = (destination / name).with_suffix(".png").absolute()
            fig.savefig(path)
            paths.append(path)
        return paths

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
            noise_uniformity=HeliosNoiseUniformityModuleOutput(
                offset=SECTION_3_OFFSET_MM,
                roi_settings=self.noise_uniformity_module.roi_settings,
                rois=self.noise_uniformity_module.as_dict(),
                noise_center=self.noise_uniformity_module.noise_center,
                mean_outer=self.noise_uniformity_module.mean_outer,
                means_diff=self.noise_uniformity_module.uniformity_difference,
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
        """Publish a PDF report containing the analysis results.

        Parameters
        ----------
        filename : str or Path
            The file to write the PDF to.
        notes : str or None
            Optional text to include in the report.
        open_file : bool
            Whether to open the PDF after creation.
        metadata : dict or None
            Extra metadata shown in the report header (e.g.
            ``{'Author': 'Jane', 'Unit': 'CT-1'}``).
        logo : Path, str, or None
            Custom logo for the report. Uses the pylinac logo if None.
        """
        analysis_title = f"{self._model} Analysis"
        cs = self.contrast_scale_module
        nu = self.noise_uniformity_module
        texts = [
            f" - {self._model} Results - ",
            f"Contrast Scale: Difference={cs.contrast_difference:.1f} HU "
            f"Noise (StdDev): {nu.center_stdev:.2f}",
            f"Uniformity Difference: {nu.uniformity_difference:.1f} HU",
        ]
        analysis_images = self.save_images(to_stream=True)
        canvas = pdf.PylinacCanvas(
            filename,
            page_title=analysis_title,
            metadata=metadata,
            logo=logo,
        )
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 4.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 4))
        for idx, text in enumerate(texts):
            canvas.add_text(text=text, location=(1.5, 23 - idx * 0.5))
        for img in analysis_images:
            canvas.add_new_page()
            canvas.add_image(img, location=(1, 5), dimensions=(18, 18))
        canvas.finish()
        if open_file:
            webbrowser.open(filename)

    def _quaac_datapoints(self) -> dict[str, QuaacDatum]:
        results = self.results_data(as_dict=True)
        data: dict[str, QuaacDatum] = {}
        data["Phantom Roll"] = QuaacDatum(
            value=results["phantom_roll_deg"],
            unit="degrees",
            description="The roll of the phantom in degrees",
        )
        data["Contrast Difference"] = QuaacDatum(
            value=results["contrast_scale"]["contrast_difference"],
            unit="HU",
            description="Plexiglass minus Water mean HU",
        )
        for bar_name, stdev in results["high_contrast"]["stdev_per_bar"].items():
            data[f"{bar_name} StdDev"] = QuaacDatum(
                value=stdev,
                unit="HU",
                description=f"Standard deviation of {bar_name} bar pattern ROI",
            )
        data["Noise StdDev"] = QuaacDatum(
            value=results["noise_uniformity"]["center_stdev"],
            unit="HU",
            description="Center ROI standard deviation (noise)",
        )
        data["Uniformity Difference"] = QuaacDatum(
            value=results["noise_uniformity"]["uniformity_difference"],
            unit="HU",
            description="Center minus edge mean HU",
        )
        return data
