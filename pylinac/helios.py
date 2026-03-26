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
        """Dump important data as a dictionary.

        Returns
        -------
        dict
            A nested dictionary with keys ``"data"`` → ``"mean_hu"`` and
            ``"std"``, each mapping ROI name to the corresponding value.
        """
        mean_hu: dict[str, float] = {}
        std: dict[str, float] = {}
        for name, roi in self.rois.items():
            mean_hu[name] = roi.mean
            std[name] = roi.std
        return {
            "data": {
                "mean_hu": mean_hu,
                "std": std,
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
        """Compute the relative Modulation Transfer Function (rMTF) from the bar-pattern ROIs.

        Each ROI targets a specific bar-pattern frequency; the spatial
        frequency (in lp/mm) is ``1 / (2 * bar_size_mm)``.  The MTF is
        derived from the min/max of the ROI (standard MTF definition)
        via :meth:`~pylinac.core.mtf.MTF.from_high_contrast_diskset`.

        Returns
        -------
        MTF
            The computed rMTF object, which can be queried for relative
            resolution at any percentage (e.g. ``mtf.relative_resolution(50)``
            returns the lp/mm at 50 % MTF).
        """
        spacings = [1 / (2 * roi["bar_size"]) for roi in self.roi_settings.values()]
        diskset = list(self.rois.values())
        return MTF.from_high_contrast_diskset(spacings=spacings, diskset=diskset)

    def as_dict(self) -> dict:
        """Dump important data as a dictionary.

        Returns
        -------
        dict
            A dictionary mapping each ROI name to the standard deviation of
            pixel values within that ROI, used as a proxy for bar-pattern
            contrast when computing the MTF.
        """
        roi_stds: dict[str, float] = {}
        for name, roi in self.rois.items():
            roi_stds[name] = roi.std
        return roi_stds

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
        """The mean HU value across all low-contrast ROIs.

        Returns
        -------
        float
            Mean of the per-ROI mean HU values.
        """
        roi_means: list[float] = []
        for roi in self.rois:
            roi_means.append(roi.mean)
        return float(np.mean(roi_means))

    @property
    def std(self) -> float:
        """The standard deviation of mean HU values across all low-contrast ROIs.

        Returns
        -------
        float
            Standard deviation of the per-ROI mean HU values.
        """
        roi_means: list[float] = []
        for roi in self.rois:
            roi_means.append(roi.mean)
        return float(np.std(roi_means))

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
    """Manages Low Contrast Detectability analysis across three adjacent slices.

    The GE Helios phantom Section 3 contains a uniform water region used to
    assess low-contrast detectability.  Rather than analysing a single slice,
    this class instantiates three :class:`HeliosLowContrastModule` objects at
    consecutive slice positions offset from :data:`SECTION_3_OFFSET_MM` by
    the inter-slice spacing defined in
    :data:`HELIOS_LOW_CONTRAST_SLICE_OFFSETS_INDEX`.

    Attributes
    ----------
    slices : dict[str, HeliosLowContrastModule]
        Mapping of slice name (``"slice_1"``, ``"slice_2"``, ``"slice_3"``)
        to the corresponding analysed module instance.
    """

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
        """The mean HU value across all slices.

        Returns
        -------
        float
            Mean of the per-slice mean HU values.
        """
        slice_means: list[float] = []
        for s in self.slices.values():
            slice_means.append(s.mean)
        mean = np.mean(slice_means)
        return float(mean)

    @property
    def std(self) -> float:
        """The average standard deviation across all slices.

        Returns
        -------
        float
            Mean of the per-slice standard deviations.
        """
        slice_stds: list[float] = []
        for s in self.slices.values():
            slice_stds.append(s.std)
        std = np.mean(slice_stds)
        return float(std)


class HeliosLowContrastMultiSliceModuleOutput(BaseModel):
    """Results for the multi-slice Low Contrast Detectability analysis.

    This class should not be called directly. It is returned by the
    ``results_data()`` method.  Use the following attributes as normal
    class attributes.

    Attributes
    ----------
    slices : dict[str, HeliosLowContrastModuleOutput]
        Per-slice results keyed by slice name (``"slice_1"``,
        ``"slice_2"``, ``"slice_3"``).  Each value contains the
        ``mean`` and ``std`` HU statistics for that individual slice.
    mean : float
        Mean HU value averaged across all three slices.
    std : float
        Average standard deviation across all three slices.
    """

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
        """Mean HU of the outer (non-center) ROIs.

        Returns
        -------
        float
            Average mean HU across the ``"12 o'clock"`` and ``"3 o'clock"``
            ROIs, representing the peripheral image uniformity.
        """
        outer_roi_keys = ["12 o'clock", "3 o'clock"]
        outer_means: list[float] = []
        for key in outer_roi_keys:
            outer_means.append(self.rois[key].mean)
        return float(np.mean(outer_means))

    @property
    def uniformity_difference(self) -> float:
        """Difference between the center mean and the average edge mean."""
        return float(self.rois["Center"].mean - self.mean_outer)

    def as_dict(self) -> dict:
        """Dump important data as a dictionary.

        Returns
        -------
        dict
            A dictionary with keys ``"mean_hu"`` and ``"std"``, each
            mapping ROI name (``"Center"``, ``"12 o'clock"``,
            ``"3 o'clock"``) to its corresponding value.
        """
        mean_hu: dict[str, float] = {}
        std: dict[str, float] = {}
        for name, roi in self.rois.items():
            mean_hu[name] = roi.mean
            std[name] = roi.std
        return {
            "mean_hu": mean_hu,
            "std": std,
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

        max_variance = variances.max()
        threshold = variances > max_variance / 2
        candidate_indices = np.argwhere(threshold)
        best_slice = int(np.mean(candidate_indices))
        return best_slice

    def find_phantom_roll(self, func: Callable | None = None) -> float:
        """Return the phantom roll angle.

        The GE Helios phantom is bracket-mounted and does not rotate; its
        roll is always zero.  The ``func`` parameter is accepted for
        interface compatibility with :class:`~pylinac.ct.CatPhanBase` but
        is not used.

        Parameters
        ----------
        func : callable or None
            Unused.  Present for API compatibility only.

        Returns
        -------
        float
            Always ``0.0``.
        """
        return 0.0

    def _module_offsets(self) -> list[float]:
        """Return the absolute z-positions (mm) of the two active phantom sections.

        Section 1 (Contrast Scale and High Contrast) is at the origin slice.
        Section 3 (Low Contrast and Noise/Uniformity) is offset by
        :data:`SECTION_3_OFFSET_MM` millimetres in the superior direction.

        Returns
        -------
        list of float
            Two-element list ``[origin_z, origin_z + SECTION_3_OFFSET_MM]``.
        """
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
        for slice_module in self.low_contrast_multi_slice.slices.values():
            modules.append(slice_module)
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
        for slice_module in self.low_contrast_multi_slice.slices.values():
            modules.append(slice_module)

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
        for key, slice_module in self.low_contrast_multi_slice.slices.items():
            modules[key] = slice_module
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
        """Return the list of all analysed phantom modules.

        Used internally by the base class for operations that iterate
        over every module (e.g. building the side-view plot).

        Returns
        -------
        list of CatPhanModule
            The six analysed modules in display order:

            * **Contrast Scale** – HU linearity / contrast scale section.
            * **High Contrast** – Spatial resolution (MTF) section.
            * **Noise & Uniformity** – Image noise and field uniformity section.
            * **Low Contrast - 1** – Low-contrast detectability, slice 1
              (``slice_1``).
            * **Low Contrast - 2** – Low-contrast detectability, slice 2
              (``slice_2``).
            * **Low Contrast - 3** – Low-contrast detectability, slice 3
              (``slice_3``).
        """
        modules = [
            self.contrast_scale_module,
            self.high_contrast_module,
            self.noise_uniformity_module,
        ]
        for slice_module in self.low_contrast_multi_slice.slices.values():
            modules.append(slice_module)
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
        """Build the dictionary of QuAAC data points from the analysis results.

        Each entry maps a human-readable label to a
        :class:`~pylinac.core.utilities.QuaacDatum` containing the numeric
        value and its physical unit.  This method is called automatically
        by the QuAAC export machinery.

        Returns
        -------
        dict[str, QuaacDatum]
            QuAAC-compatible data points covering all outputs produced by
            :meth:`_generate_results_data`: phantom roll; per-ROI
            contrast-scale mean HU and standard deviation; per-bar-pattern
            high-contrast ROI standard deviations and MTF at each 10 %
            resolution step; per-slice low-contrast mean and standard
            deviation for each of the three analysis slices; aggregate
            low-contrast mean and standard deviation; per-position
            noise/uniformity ROI mean HU and standard deviation; noise
            center standard deviation; mean outer ROI HU; and uniformity
            difference.
        """
        results_data = self.results_data(as_dict=True)

        phantom_roll = QuaacDatum(
            value=results_data["phantom_roll_deg"],
            unit="degrees",
            description="The roll of the phantom in the image",
        )

        contrast_scale_mean_hu_data = results_data["contrast_scale"]["rois"]["data"][
            "mean_hu"
        ]
        contrast_scale_mean_hu: dict[str, QuaacDatum] = {}
        for name, hu_value in contrast_scale_mean_hu_data.items():
            label = f"Contrast scale {name} mean HU"
            contrast_scale_mean_hu[label] = QuaacDatum(
                value=hu_value,
                unit="HU",
            )

        contrast_scale_std_data = results_data["contrast_scale"]["rois"]["data"]["std"]
        contrast_scale_std: dict[str, QuaacDatum] = {}
        for name, std_value in contrast_scale_std_data.items():
            label = f"Contrast scale {name} std"
            contrast_scale_std[label] = QuaacDatum(
                value=std_value,
                unit="HU",
            )

        high_contrast_roi_stds: dict[str, QuaacDatum] = {}
        for name, std_value in results_data["high_contrast"]["rois"].items():
            label = f"High contrast {name} ROI std"
            high_contrast_roi_stds[label] = QuaacDatum(
                value=std_value,
                unit="HU",
            )

        mtf_datapoints: dict[str, QuaacDatum] = {}
        for resolution, lp_mm in results_data["high_contrast"]["mtf_lp_mm"].items():
            label = f"High contrast MTF {resolution}%"
            mtf_datapoints[label] = QuaacDatum(
                value=lp_mm,
                unit="lp/mm",
            )

        low_contrast_slice_means: dict[str, QuaacDatum] = {}
        low_contrast_slice_stds: dict[str, QuaacDatum] = {}
        for slice_name, slice_data in results_data["low_contrast"]["slices"].items():
            mean_label = f"Low contrast {slice_name} mean"
            low_contrast_slice_means[mean_label] = QuaacDatum(
                value=slice_data["mean"],
                unit="HU",
            )
            std_label = f"Low contrast {slice_name} std"
            low_contrast_slice_stds[std_label] = QuaacDatum(
                value=slice_data["std"],
                unit="HU",
            )

        low_contrast_mean = QuaacDatum(
            value=results_data["low_contrast"]["mean"],
            unit="HU",
        )
        low_contrast_std = QuaacDatum(
            value=results_data["low_contrast"]["std"],
            unit="HU",
        )

        noise_uniformity_rois_mean_hu: dict[str, QuaacDatum] = {}
        for name, hu_value in results_data["noise_uniformity"]["rois"][
            "mean_hu"
        ].items():
            label = f"Noise uniformity {name} mean HU"
            noise_uniformity_rois_mean_hu[label] = QuaacDatum(
                value=hu_value,
                unit="HU",
            )

        noise_uniformity_rois_std: dict[str, QuaacDatum] = {}
        for name, std_value in results_data["noise_uniformity"]["rois"]["std"].items():
            label = f"Noise uniformity {name} std"
            noise_uniformity_rois_std[label] = QuaacDatum(
                value=std_value,
                unit="HU",
            )

        noise_center_std = QuaacDatum(
            value=results_data["noise_uniformity"]["noise_center_std"],
            unit="HU",
        )
        mean_outer = QuaacDatum(
            value=results_data["noise_uniformity"]["mean_outer"],
            unit="HU",
        )
        uniformity_difference = QuaacDatum(
            value=results_data["noise_uniformity"]["means_diff"],
            unit="HU",
        )

        return {
            "Phantom Roll": phantom_roll,
            **contrast_scale_mean_hu,
            **contrast_scale_std,
            **high_contrast_roi_stds,
            **mtf_datapoints,
            **low_contrast_slice_means,
            **low_contrast_slice_stds,
            "Low contrast Mean": low_contrast_mean,
            "Low contrast Std": low_contrast_std,
            **noise_uniformity_rois_mean_hu,
            **noise_uniformity_rois_std,
            "Noise Center Std": noise_center_std,
            "Mean Outer": mean_outer,
            "Uniformity Difference": uniformity_difference,
        }

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
                canvas.add_text(text=text, location=(2.5, 24 - idx * 0.5))
                idx += 1
        for page, img in enumerate(analysis_images):
            canvas.add_new_page()
            canvas.add_image(img, location=(1, 5), dimensions=(18, 18))
        canvas.finish()

        if open_file:
            webbrowser.open(filename)

    def results(self, as_str: bool = True) -> str | tuple:
        """Return the results of the analysis as a string. Use with print().

        All scalar and per-ROI outputs produced by
        :meth:`_generate_results_data` are represented: phantom roll,
        contrast-scale per-ROI mean HU and standard deviation, contrast
        difference, high-contrast per-bar-pattern ROI standard deviation,
        MTF at every 10 % resolution step, per-slice low-contrast mean and
        standard deviation for each of the three analysis slices followed by
        the aggregate low-contrast mean and standard deviation,
        noise/uniformity per-position mean HU and standard deviation, noise
        center standard deviation, mean outer HU, and uniformity difference.

        Parameters
        ----------
        as_str : bool
            If ``True`` (default) return a single newline-joined string
            suitable for ``print()``.  If ``False`` return a tuple of
            individual result strings.

        Returns
        -------
        str or tuple of str
            The formatted results.
        """
        lines: list[str] = []

        title = f" - {self._model} Results - "
        lines.append(title)

        phantom_roll = f"Phantom Roll: {self.catphan_roll:2.2f} deg"
        lines.append(phantom_roll)

        for name, roi in self.contrast_scale_module.rois.items():
            contrast_roi_mean = f"Contrast Scale {name} Mean HU: {roi.mean:2.2f}"
            lines.append(contrast_roi_mean)
            contrast_roi_std = f"Contrast Scale {name} Std: {roi.std:2.2f}"
            lines.append(contrast_roi_std)

        contrast_diff = f"Contrast Difference: {self.contrast_scale_module.contrast_difference:2.2f}"
        lines.append(contrast_diff)

        for name, roi in self.high_contrast_module.rois.items():
            hc_roi_std = f"High Contrast {name} ROI Std: {roi.std:2.2f}"
            lines.append(hc_roi_std)

        for resolution in range(10, 91, 10):
            lp_mm = self.high_contrast_module.mtf.relative_resolution(resolution)
            mtf_line = f"MTF {resolution}% (lp/mm): {lp_mm:2.2f}"
            lines.append(mtf_line)

        for slice_name, slice_module in self.low_contrast_multi_slice.slices.items():
            lc_slice_mean = f"Low Contrast {slice_name} Mean: {slice_module.mean:2.2f}"
            lines.append(lc_slice_mean)
            lc_slice_std = f"Low Contrast {slice_name} Std: {slice_module.std:2.2f}"
            lines.append(lc_slice_std)

        low_contrast_mean = (
            f"Low Contrast Mean: {self.low_contrast_multi_slice.mean:2.2f}"
        )
        lines.append(low_contrast_mean)

        low_contrast_std = (
            f"Low Contrast Standard Deviation: {self.low_contrast_multi_slice.std:2.2f}"
        )
        lines.append(low_contrast_std)

        for name, roi in self.noise_uniformity_module.rois.items():
            nu_roi_mean = f"Noise Uniformity {name} Mean HU: {roi.mean:2.2f}"
            lines.append(nu_roi_mean)
            nu_roi_std = f"Noise Uniformity {name} Std: {roi.std:2.2f}"
            lines.append(nu_roi_std)

        noise_center_std = (
            f"Noise Center Std: {self.noise_uniformity_module.noise_center_std:2.2f}"
        )
        lines.append(noise_center_std)

        mean_outer = f"Mean Outer HU: {self.noise_uniformity_module.mean_outer:2.2f}"
        lines.append(mean_outer)

        uniformity_diff = f"Uniformity Difference: {self.noise_uniformity_module.uniformity_difference:2.2f}"
        lines.append(uniformity_diff)

        string = tuple(lines)
        if as_str:
            return "\n".join(string)
        else:
            return string

    def _generate_results_data(self) -> GEHeliosResult:
        resolutions = range(10, 91, 10)  # 10-90% in 10% increments
        mtfs: dict[int, float] = {}
        for resolution in resolutions:
            mtfs[resolution] = self.high_contrast_module.mtf.relative_resolution(
                resolution
            )

        slice_outputs: dict[str, HeliosLowContrastModuleOutput] = {}
        for k, v in self.low_contrast_multi_slice.slices.items():
            slice_outputs[k] = HeliosLowContrastModuleOutput(
                offset=self.low_contrast_multi_slice.roi_settings[k]["offset"],
                settings={
                    "cell_size": v.cell_size,
                    "num_cells": v.num_cells,
                },
                mean=v.mean,
                std=v.std,
            )

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
                slices=slice_outputs,
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
