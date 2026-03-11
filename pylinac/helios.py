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
    plexiglass_hu: float = Field(
        description="Mean HU of the Plexiglass ROI.",
        title="Plexiglass Mean HU",
    )
    plexiglass_stdev: float = Field(
        description="Standard deviation of the Plexiglass ROI.",
        title="Plexiglass standard deviation",
    )
    water_hu: float = Field(
        description="Mean HU of the Water ROI.",
        title="Water Mean HU",
    )
    water_stdev: float = Field(
        description="Standard deviation of the Water ROI.",
        title="Water standard deviation",
    )
    contrast_difference: float = Field(
        description="Difference in mean HU between Plexiglass and Water.",
        title="Contrast Difference (HU)",
    )


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
    def stdev_per_bar(self) -> dict[str, float]:
        """Standard deviation measured within each bar-pattern ROI."""
        return {name: float(roi.std) for name, roi in self.rois.items()}

    @property
    def mtf(self) -> MTF:
        spacings = [1 / (2 * roi["bar_size"]) for roi in self.roi_settings.values()]
        diskset = list(self.rois.values())
        return MTF.from_high_contrast_diskset(spacings=spacings, diskset=diskset)

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
    stdev_per_bar: dict[str, float] = Field(
        description="Standard deviation measured in each bar pattern ROI. Keys are the bar sizes.",
        title="Standard deviation per Bar Pattern",
    )
    mtf_lp_mm: dict[float, float] = Field(
        description="Relative MTF at each spatial frequency (lp/mm), normalized to the coarsest bar pattern.",
        title="MTF (lp/mm)",
    )


class HeliosNoiseUniformityModule(CatPhanModule):
    """Class for analysis of the Noise & Uniformity."""

    common_name = "Noise & Uniformity"
    attr_name = "noise_uniformity_module"
    center_roi_settings = {
        "Noise": {"width": 25, "height": 25, "distance": 0, "angle": 0},
        "Uniformity": {"width": 15, "height": 15, "distance": 0, "angle": 0},
    }
    edge_roi_settings = {
        "12 o'clock": {"width": 15, "height": 15, "distance": 75, "angle": -90},
        "3 o'clock": {"width": 15, "height": 15, "distance": 75, "angle": 0},
    }
    center_rois: dict
    edge_rois: dict

    def _setup_rois(self) -> None:
        self.center_rois = {}
        self.edge_rois = {}
        for name, setting in self.center_roi_settings.items():
            self.center_rois[name] = RectangleROI.from_phantom_center(
                array=self.image,
                width=setting["width_pixels"],
                height=setting["height_pixels"],
                angle=setting["angle_corrected"],
                dist_from_center=setting["distance_pixels"],
                phantom_center=self.phan_center,
            )
        for name, setting in self.edge_roi_settings.items():
            self.edge_rois[name] = RectangleROI.from_phantom_center(
                array=self.image,
                width=setting["width_pixels"],
                height=setting["height_pixels"],
                angle=setting["angle_corrected"],
                dist_from_center=setting["distance_pixels"],
                phantom_center=self.phan_center,
            )

    @property
    def center_mean_hu(self) -> float:
        """Mean HU of the 15 mm center ROI."""
        return float(self.center_rois["Uniformity"].mean)

    @property
    def center_stdev(self) -> float:
        """Standard deviation (noise) of the 25 mm center ROI."""
        return self.center_rois["Noise"].std

    @property
    def uniformity_difference(self) -> float:
        """Difference between the center mean and the average edge mean."""
        edge_mean = np.mean([roi.mean for roi in self.edge_rois.values()])
        return float(self.center_mean_hu - edge_mean)

    def plot_rois(self, axis: plt.Axes) -> None:
        """Plot the ROIs to the axis."""
        for roi in self.center_rois.values():
            roi.plot2axes(axis, edgecolor="blue")
        for roi in self.edge_rois.values():
            roi.plot2axes(axis, edgecolor="blue")

    def plotly_rois(self, fig: go.Figure) -> None:
        """Plot the ROIs to a Plotly figure."""
        for name, roi in self.center_rois.items():
            roi.plotly(fig, line_color="blue", name=name)
        for name, roi in self.edge_rois.items():
            roi.plotly(fig, line_color="blue", name=name)


class HeliosNoiseUniformityModuleOutput(BaseModel):
    """This class should not be called directly. It is returned by the ``results_data()`` method.

    Use the following attributes as normal class attributes."""

    offset: float = Field(
        description="The offset of this module slice from the origin slice in mm."
    )
    center_mean_hu: float = Field(
        description="Mean HU of the center uniformity ROI.",
        title="Center Mean HU",
    )
    center_stdev: float = Field(
        description="Standard deviation of the center noise ROI (25 mm box).",
        title="Center StdDev (Noise)",
    )
    edge_rois: dict = Field(
        description="Mean HU values of each edge ROI, keyed by position name.",
        title="Edge ROI Mean HU Values",
    )
    uniformity_difference: float = Field(
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
        return GEHeliosResult(
            phantom_model=self._model,
            phantom_roll_deg=self.catphan_roll,
            origin_slice=self.origin_slice,
            num_images=self.num_images,
            contrast_scale=HeliosContrastScaleModuleOutput(
                offset=0,
                plexiglass_hu=self.contrast_scale_module.rois["Plexiglass"].mean,
                plexiglass_stdev=self.contrast_scale_module.rois["Plexiglass"].std,
                water_hu=self.contrast_scale_module.rois["Water"].mean,
                water_stdev=self.contrast_scale_module.rois["Water"].std,
                contrast_difference=self.contrast_scale_module.contrast_difference,
            ),
            high_contrast=HeliosHighContrastModuleOutput(
                offset=0,
                stdev_per_bar={
                    name: roi.std
                    for name, roi in self.high_contrast_module.rois.items()
                },
                mtf_lp_mm=self.high_contrast_module.mtf.norm_mtfs,
            ),
            noise_uniformity=HeliosNoiseUniformityModuleOutput(
                offset=SECTION_3_OFFSET_MM,
                center_mean_hu=self.noise_uniformity_module.center_mean_hu,
                center_stdev=self.noise_uniformity_module.center_stdev,
                edge_rois={
                    name: float(roi.mean)
                    for name, roi in self.noise_uniformity_module.edge_rois.items()
                },
                uniformity_difference=self.noise_uniformity_module.uniformity_difference,
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
