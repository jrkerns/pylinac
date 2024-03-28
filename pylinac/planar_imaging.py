"""The planar imaging module analyzes phantom images taken with the kV or MV imager in 2D.
The following phantoms are supported:

* Leeds TOR 18
* Standard Imaging QC-3
* Standard Imaging QC-kV
* Las Vegas
* Elekta Las Vegas
* Doselab MC2 MV
* Doselab MC2 kV
* SNC kV
* SNC MV
* PTW EPID QC

Features:

* **Automatic phantom localization** - Set up your phantom any way you like; automatic positioning,
  angle, and inversion correction mean you can set up how you like, nor will setup variations give you headache.
* **High and low contrast determination** - Analyze both low and high contrast ROIs. Set thresholds
  as you see fit.
"""
from __future__ import annotations

import io
import math
import os.path as osp
import warnings
import webbrowser
from functools import cached_property
from pathlib import Path
from typing import BinaryIO, Callable, Literal

import matplotlib.pyplot as plt
import numpy as np
from py_linq import Enumerable
from scipy.ndimage import median_filter
from skimage import exposure, feature, measure
from skimage.measure._regionprops import RegionProperties

from . import Normalization
from .core import geometry, image, pdf
from .core.contrast import Contrast
from .core.decorators import lru_cache
from .core.geometry import Circle, Point, Rectangle, Vector
from .core.io import get_url, retrieve_demo_file
from .core.mtf import MTF
from .core.profile import CollapsedCircleProfile, FWXMProfilePhysical
from .core.roi import DiskROI, HighContrastDiskROI, LowContrastDiskROI, bbox_center
from .core.utilities import ResultBase, ResultsDataMixin
from .metrics.image import SizedDiskLocator


class PlanarResult(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    analysis_type: str  #:
    median_contrast: float  #:
    median_cnr: float  #:
    num_contrast_rois_seen: int  #:
    phantom_center_x_y: tuple[float, float]  #:
    low_contrast_rois: list[dict]  #:
    phantom_area: float  #: The area of the phantom in pixels^2
    mtf_lp_mm: tuple[float, float, float] | None = None  #:
    percent_integral_uniformity: float | None = None  #:


def _middle_of_bbox_region(region: RegionProperties) -> tuple:
    return (
        (region.bbox[2] - region.bbox[0]) / 2 + region.bbox[0],
        (region.bbox[3] - region.bbox[1]) / 2 + region.bbox[1],
    )


def is_square(region: RegionProperties, instance: object, rtol=0.2) -> bool:
    """Whether the region has symmetric height and width"""
    height = region.bbox[2] - region.bbox[0]
    width = region.bbox[3] - region.bbox[1]
    return math.isclose(height / width, 1, rel_tol=rtol)


def is_centered(region: RegionProperties, instance: object, rtol=0.3) -> bool:
    """Whether the region is centered on the image"""
    img_center = (instance.image.center.y, instance.image.center.x)
    # we don't want centroid because that could be offset by missing lengths of the outline. Center of bbox is more robust
    return np.allclose(_middle_of_bbox_region(region), img_center, rtol=rtol)


def is_right_size(region: RegionProperties, instance: object, rtol=0.1) -> bool:
    """Whether the region is close to the expected size of the phantom, given the SSD and physical phantom size."""
    return bool(
        np.isclose(
            region.bbox_area,
            instance.phantom_bbox_size_px,
            rtol=rtol,
        )
    )


def percent_integral_uniformity(max: float, min: float) -> float:
    """Calculate the percent integral uniformity. A small constant is
    added to avoid possible division by zero."""
    return 100 * (1 - (max - min + 1e-6) / (max + min + 1e-6))


class ImagePhantomBase(ResultsDataMixin[PlanarResult]):
    """Base class for planar phantom classes.

    Attributes
    ----------
    common_name : str
        The human-readable name of the phantom. Used in plots and PDF report.
    phantom_outline_object : {None, 'Circle', 'Rectangle'}
        What type of outline to display on the plotted image. Helps to visually determine the accuracy of the
        phantom size, position, and scale.
    high_contrast_rois : list
        :class:`~pylinac.core.roi.HighContrastDiskROI` instances of the
        high contrast line pair regions.
    high_contrast_roi_settings : dict
        Settings of the placement of the high-contrast ROIs.
    low_contrast_rois : list
        :class:`~pylinac.core.roi.LowContrastDiskROI` instances of the low
        contrast ROIs, other than the reference ROI (below).
    low_contrast_roi_settings : dict
        Settings of the placement of the low-contrast ROIs.
    low_contrast_background_rois : list
        :class:`~pylinac.core.roi.LowContrastDiskROI` instances of the low
        contrast background ROIs.
    low_contrast_background_roi_settings : dict
        Settings of the placement of the background low-contrast ROIs.
    low_contrast_background_value : float
        The average pixel value of all the low-contrast background ROIs.
    detection_conditions: list of callables
        This should be a list of functions that return a boolean. It is used for finding the phantom outline in the image.
        E.g. is_at_center().
    phantom_bbox_size_mm2: float
        This is the expected size of the **BOUNDING BOX** of the phantom. Additionally, it is usually smaller than the
        physical bounding box because we sometimes detect an inner ring/square. Typically, x0.9-1.0 of the physical size.
    """

    _demo_filename: str
    common_name: str
    high_contrast_roi_settings = {}
    high_contrast_rois = []
    low_contrast_roi_settings = {}
    low_contrast_rois = []
    low_contrast_background_roi_settings = {}
    low_contrast_background_rois = []
    low_contrast_background_value = None
    phantom_outline_object = None
    detection_conditions: list[Callable] = [is_centered, is_right_size]
    detection_canny_settings = {"sigma": 2, "percentiles": (0.001, 0.01)}
    phantom_bbox_size_mm2: float
    roi_match_condition: Literal["max", "closest"] = "max"
    mtf: MTF
    _ssd: float

    def __init__(
        self,
        filepath: str | BinaryIO | Path,
        normalize: bool = True,
        image_kwargs: dict | None = None,
    ):
        """
        Parameters
        ----------
        filepath : str
            Path to the image file.
        normalize: bool
            Whether to "ground" and normalize the image. This can affect contrast measurements, but for
            backwards compatibility this is True. You may want to set this to False if trying to compare with other software.
        image_kwargs : dict
            Keywords passed to the image load function; this would include things like DPI or SID if applicable
        """
        img_kwargs = image_kwargs or {}
        self.image = image.load(filepath, **img_kwargs)
        if normalize:
            self.image.ground()
            self.image.normalize()
        self._angle_override = None
        self._size_override = None
        self._center_override = None
        self._high_contrast_threshold = None
        self._low_contrast_threshold = None

    @classmethod
    def from_demo_image(cls):
        """Instantiate and load the demo image."""
        demo_file = retrieve_demo_file(name=cls._demo_filename)
        return cls(demo_file)

    @classmethod
    def from_url(cls, url: str):
        """
        Parameters
        ----------
        url : str
            The URL to the image.
        """
        image_file = get_url(url)
        return cls(image_file)

    def _preprocess(self):
        pass

    def _check_inversion(self):
        pass

    def window_floor(self) -> float | None:
        """The value to use as the minimum when displaying the image (see https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.imshow.html)
        Helps show contrast of images, specifically if there is an open background"""
        return

    def window_ceiling(self) -> float | None:
        """The value to use as the maximum when displaying the image. Helps show contrast of images, specifically if there is an open background"""
        return

    @property
    def magnification_factor(self) -> float:
        """The mag factor of the image based on SSD vs SAD"""
        return self.image.sad / self._ssd

    @property
    def phantom_bbox_size_px(self) -> float:
        """The phantom bounding box size in pixels^2 at the isoplane."""
        return (
            self.phantom_bbox_size_mm2
            * (self.image.dpmm**2)
            * (self.magnification_factor**2)
        )

    @cached_property
    def phantom_ski_region(self) -> RegionProperties:
        """The skimage region of the phantom outline."""
        regions = self._get_canny_regions()
        sorted_regions = (
            Enumerable(regions)
            .where(lambda r: r.area_bbox > 100)
            .order_by_descending(lambda r: r.area_bbox)
            .to_list()
        )
        # sorted_regions = sorted(regions, key=lambda r: r.area_bbox, reverse=True)
        # search through all the canny ROIs to see which ones pass the detection conditions
        blobs = []
        for phantom_idx, region in enumerate(sorted_regions):
            conditions_met = [
                condition(region, self) for condition in self.detection_conditions
            ]
            if all(conditions_met):
                blobs.append(phantom_idx)

        if not blobs:
            raise ValueError(
                "Unable to find the phantom in the image. Potential solutions: check the SSD was passed correctly, check that the phantom isn't at the edge of the field, check that the phantom is centered along the CAX."
            )

        if self.roi_match_condition == "max":
            # take the biggest ROI and call that the phantom outline
            best_roi_idx = np.argsort(
                [sorted_regions[phan].bbox_area for phan in blobs]
            )[-1]
        elif (
            self.roi_match_condition == "closest"
        ):  # take the one most similar in size to known size
            best_roi_idx = np.argsort(
                [
                    abs(sorted_regions[phan].bbox_area - self.phantom_bbox_size_px)
                    for phan in blobs
                ]
            )[0]
        phantom_idx = blobs[best_roi_idx]

        return sorted_regions[phantom_idx]

    def analyze(
        self,
        low_contrast_threshold: float = 0.05,
        high_contrast_threshold: float = 0.5,
        invert: bool = False,
        angle_override: float | None = None,
        center_override: tuple | None = None,
        size_override: float | None = None,
        ssd: float | Literal["auto"] = "auto",
        low_contrast_method: str = Contrast.MICHELSON,
        visibility_threshold: float = 100,
    ) -> None:
        """Analyze the phantom using the provided thresholds and settings.

        Parameters
        ----------
        low_contrast_threshold : float
            This is the contrast threshold value which defines any low-contrast ROI as passing or failing.
        high_contrast_threshold : float
            This is the contrast threshold value which defines any high-contrast ROI as passing or failing.
        invert : bool
            Whether to force an inversion of the image. This is useful if pylinac's automatic inversion algorithm fails
            to properly invert the image.
        angle_override : None, float
            A manual override of the angle of the phantom. If None, pylinac will automatically determine the angle. If
            a value is passed, this value will override the automatic detection.

            .. Note::

                0 is pointing from the center toward the right and positive values go counterclockwise.

        center_override : None, 2-element tuple
            A manual override of the center point of the phantom. If None, pylinac will automatically determine the center. If
            a value is passed, this value will override the automatic detection. Format is (x, y)/(col, row).
        size_override : None, float
            A manual override of the relative size of the phantom. This size value is used to scale the positions of
            the ROIs from the center. If None, pylinac will automatically determine the size.
            If a value is passed, this value will override the automatic sizing.

            .. Note::

                 This value is not necessarily the physical size of the phantom. It is an arbitrary value.
        ssd
            The SSD of the phantom itself in mm. If set to "auto", will first search for the phantom at the SAD, then at 5cm above the SID.
        low_contrast_method
            The equation to use for calculating low contrast.
        visibility_threshold
            The threshold for whether an ROI is "seen".
        """
        self._angle_override = angle_override
        self._center_override = center_override
        self._size_override = size_override
        self._high_contrast_threshold = high_contrast_threshold
        self._low_contrast_threshold = low_contrast_threshold
        self._low_contrast_method = low_contrast_method
        self.visibility_threshold = visibility_threshold
        self._ssd = ssd
        self._find_ssd()
        self._check_inversion()
        if invert:
            self.image.invert()
        self._preprocess()
        if self.high_contrast_roi_settings:
            self.high_contrast_rois = self._sample_high_contrast_rois()
            # generate rMTF
            spacings = [
                roi["lp/mm"] for roi in self.high_contrast_roi_settings.values()
            ]
            self.mtf = MTF.from_high_contrast_diskset(
                diskset=self.high_contrast_rois, spacings=spacings
            )
        if self.low_contrast_background_roi_settings:
            (
                self.low_contrast_background_rois,
                self.low_contrast_background_value,
            ) = self._sample_low_contrast_background_rois()
        if self.low_contrast_roi_settings:
            self.low_contrast_rois = self._sample_low_contrast_rois()

    def _sample_low_contrast_rois(self) -> list[LowContrastDiskROI]:
        """Sample the low-contrast sample regions for calculating contrast values."""
        lc_rois = []
        for stng in self.low_contrast_roi_settings.values():
            roi = LowContrastDiskROI(
                self.image,
                self.phantom_angle + stng["angle"],
                self.phantom_radius * stng["roi radius"],
                self.phantom_radius * stng["distance from center"],
                self.phantom_center,
                self._low_contrast_threshold,
                self.low_contrast_background_value,
                contrast_method=self._low_contrast_method,
                visibility_threshold=self.visibility_threshold,
            )
            lc_rois.append(roi)
        return lc_rois

    def _sample_low_contrast_background_rois(
        self,
    ) -> tuple[list[LowContrastDiskROI], float]:
        """Sample the low-contrast background regions for calculating contrast values."""
        bg_rois = []
        for stng in self.low_contrast_background_roi_settings.values():
            roi = LowContrastDiskROI(
                self.image,
                self.phantom_angle + stng["angle"],
                self.phantom_radius * stng["roi radius"],
                self.phantom_radius * stng["distance from center"],
                self.phantom_center,
                self._low_contrast_threshold,
            )
            bg_rois.append(roi)
        avg_bg = np.mean([roi.pixel_value for roi in bg_rois])
        return bg_rois, avg_bg

    def _sample_high_contrast_rois(self) -> list[HighContrastDiskROI]:
        """Sample the high-contrast line pair regions."""
        hc_rois = []
        for stng in self.high_contrast_roi_settings.values():
            roi = HighContrastDiskROI(
                self.image,
                self.phantom_angle + stng["angle"],
                self.phantom_radius * stng["roi radius"],
                self.phantom_radius * stng["distance from center"],
                self.phantom_center,
                self._high_contrast_threshold,
            )
            hc_rois.append(roi)
        return hc_rois

    def save_analyzed_image(
        self,
        filename: None | str | BinaryIO = None,
        split_plots: bool = False,
        to_streams: bool = False,
        **kwargs,
    ) -> dict[str, BinaryIO] | list[str] | None:
        """Save the analyzed image to disk or to stream. Kwargs are passed to plt.savefig()

        Parameters
        ----------
        filename : None, str, stream
            A string representing where to save the file to. If split_plots and to_streams are both true, leave as None as newly-created streams are returned.
        split_plots: bool
            If split_plots is True, multiple files will be created that append a name. E.g. `my_file.png` will become `my_file_image.png`, `my_file_mtf.png`, etc.
            If to_streams is False, a list of new filenames will be returned
        to_streams: bool
            This only matters if split_plots is True. If both of these are true, multiple streams will be created and returned as a dict.
        """
        if filename is None and to_streams is False:
            raise ValueError("Must pass in a filename unless saving to streams.")
        figs, names = self.plot_analyzed_image(
            show=False, split_plots=split_plots, **kwargs
        )
        # remove plot keywords as savefig complains about extra kwargs
        for key in (
            "image",
            "low_contrast",
            "high_contrast",
            "show",
        ):
            kwargs.pop(key, None)
        if not split_plots:
            plt.savefig(filename, **kwargs)
        else:
            # append names to filename if it's file-like
            if not to_streams:
                filenames = []
                f, ext = osp.splitext(filename)
                for name in names:
                    filenames.append(f + "_" + name + ext)
            else:  # it's a stream buffer
                filenames = [io.BytesIO() for _ in names]
            for fig, name in zip(figs, filenames):
                fig.savefig(name, **kwargs)
            if to_streams:
                return {name: stream for name, stream in zip(names, filenames)}
            if split_plots:
                return filenames

    def _get_canny_regions(self) -> list[RegionProperties]:
        """Compute the canny edges of the image and return the connected regions found."""
        # compute the canny edges with very low thresholds (detects nearly everything)
        canny_img = feature.canny(
            self.image.array,
            low_threshold=self.detection_canny_settings["percentiles"][0],
            high_threshold=self.detection_canny_settings["percentiles"][1],
            use_quantiles=True,
            sigma=self.detection_canny_settings["sigma"],
        )

        # label the canny edge regions
        labeled = measure.label(canny_img)
        regions = measure.regionprops(labeled, intensity_image=self.image.array)
        return regions

    def _create_phantom_outline_object(self) -> tuple[Rectangle | Circle, dict]:
        """Construct the phantom outline object which will be plotted on the image for visual inspection."""
        outline_type = list(self.phantom_outline_object)[0]
        outline_settings = list(self.phantom_outline_object.values())[0]
        settings = {}
        if outline_type == "Rectangle":
            side_a = self.phantom_radius * outline_settings["width ratio"]
            side_b = self.phantom_radius * outline_settings["height ratio"]
            half_hyp = np.sqrt(side_a**2 + side_b**2) / 2
            internal_angle = np.rad2deg(np.arctan(side_b / side_a))
            new_x = self.phantom_center.x + half_hyp * (
                geometry.cos(internal_angle)
                - geometry.cos(internal_angle + self.phantom_angle)
            )
            new_y = self.phantom_center.y + half_hyp * (
                geometry.sin(internal_angle)
                - geometry.sin(internal_angle + self.phantom_angle)
            )
            obj = Rectangle(
                width=self.phantom_radius * outline_settings["width ratio"],
                height=self.phantom_radius * outline_settings["height ratio"],
                center=Point(new_x, new_y),
            )
            settings["angle"] = self.phantom_angle
        elif outline_type == "Circle":
            obj = Circle(
                center_point=self.phantom_center,
                radius=self.phantom_radius * outline_settings["radius ratio"],
            )
        else:
            raise ValueError(
                "An outline object was passed but was not a Circle or Rectangle."
            )
        return obj, settings

    def percent_integral_uniformity(
        self, percentiles: tuple[float, float] = (1, 99)
    ) -> float | None:
        """Calculate and return the percent integral uniformity (PIU). This uses
        a similar equation as ACR does for CT protocols. The PIU is calculated
        over all the low contrast ROIs and the lowest (worst) PIU is returned.

        If the phantom does not contain low-contrast ROIs, None is returned."""
        if not self.low_contrast_rois:
            return
        pius = []
        for roi in self.low_contrast_rois:
            low = roi.percentile(percentiles[0])
            high = roi.percentile(percentiles[1])
            pius.append(percent_integral_uniformity(max=high, min=low))
        return min(pius)

    def plot_analyzed_image(
        self,
        image: bool = True,
        low_contrast: bool = True,
        high_contrast: bool = True,
        show: bool = True,
        split_plots: bool = False,
        **plt_kwargs: dict,
    ) -> tuple[list[plt.Figure], list[str]]:
        """Plot the analyzed image.

        Parameters
        ----------
        image : bool
            Show the image.
        low_contrast : bool
            Show the low contrast values plot.
        high_contrast : bool
            Show the high contrast values plot.
        show : bool
            Whether to actually show the image when called.
        split_plots : bool
            Whether to split the resulting image into individual plots. Useful for saving images into individual files.
        plt_kwargs : dict
            Keyword args passed to the plt.figure() method. Allows one to set things like figure size.
        """
        plot_low_contrast = low_contrast and any(self.low_contrast_rois)
        plot_high_contrast = high_contrast and any(self.high_contrast_rois)
        num_plots = sum((image, plot_low_contrast, plot_high_contrast))
        if num_plots < 1:
            warnings.warn(
                "Nothing was plotted because either all parameters were false or there were no actual high/low ROIs"
            )
            return
        # set up axes and make axes iterable
        figs = []
        names = []
        if split_plots:
            axes = []
            for n in range(num_plots):
                fig, axis = plt.subplots(1, **plt_kwargs)
                figs.append(fig)
                axes.append(axis)
        else:
            fig, axes = plt.subplots(1, num_plots, **plt_kwargs)
            fig.subplots_adjust(wspace=0.4)
        if num_plots < 2:
            axes = (axes,)
        axes = iter(axes)

        # plot the marked image
        if image:
            img_ax = next(axes)
            names.append("image")
            self.image.plot(
                ax=img_ax,
                show=False,
                vmin=self.window_floor(),
                vmax=self.window_ceiling(),
            )
            img_ax.axis("off")
            img_ax.set_title(f"{self.common_name} Phantom Analysis")

            # plot the outline image
            if self.phantom_outline_object is not None:
                outline_obj, settings = self._create_phantom_outline_object()
                outline_obj.plot2axes(img_ax, edgecolor="b", **settings)
            # plot the low contrast background ROIs
            for roi in self.low_contrast_background_rois:
                roi.plot2axes(img_ax, edgecolor="b")
            # plot the low contrast ROIs
            for roi in self.low_contrast_rois:
                roi.plot2axes(img_ax, edgecolor=roi.plot_color)
            # plot the high-contrast ROIs along w/ pass/fail coloration
            if self.high_contrast_rois:
                for roi, mtf in zip(
                    self.high_contrast_rois, self.mtf.norm_mtfs.values()
                ):
                    color = "b" if mtf > self._high_contrast_threshold else "r"
                    roi.plot2axes(img_ax, edgecolor=color)
            # plot the center of the detected ROI; used for qualitative eval of detection algorithm
            img_ax.scatter(x=self.phantom_center.x, y=self.phantom_center.y, marker="x")

        # plot the low contrast value graph
        if plot_low_contrast:
            lowcon_ax = next(axes)
            names.append("low_contrast")
            self._plot_lowcontrast_graph(lowcon_ax)

        # plot the high contrast MTF graph
        if plot_high_contrast:
            hicon_ax = next(axes)
            names.append("high_contrast")
            self._plot_highcontrast_graph(hicon_ax)

        plt.tight_layout()
        if show:
            plt.show()
        return figs, names

    def _plot_lowcontrast_graph(self, axes: plt.Axes):
        """Plot the low contrast ROIs to an axes."""
        (line1,) = axes.plot(
            [roi.contrast for roi in self.low_contrast_rois],
            marker="o",
            color="m",
            label="Contrast",
        )
        axes.axhline(self._low_contrast_threshold, color="m")
        axes.grid(True)
        axes.set_title("Low-frequency Contrast")
        axes.set_xlabel("ROI #")
        axes.set_ylabel("Contrast")
        axes2 = axes.twinx()
        (line2,) = axes2.plot(
            [roi.contrast_to_noise for roi in self.low_contrast_rois],
            marker="^",
            label="CNR",
        )
        axes2.set_ylabel("CNR")
        axes.legend(handles=[line1, line2])

    def _plot_highcontrast_graph(self, axes: plt.Axes):
        """Plot the high contrast ROIs to an axes."""
        axes.plot(self.mtf.spacings, list(self.mtf.norm_mtfs.values()), marker="*")
        axes.axhline(self._high_contrast_threshold, color="k")
        axes.grid(True)
        axes.set_title("High-frequency rMTF")
        axes.set_xlabel("Line pairs / mm")
        axes.set_ylabel("relative MTF")

    def results(self, as_list: bool = False) -> str | list[str]:
        """Return the results of the analysis.

        Parameters
        ----------
        as_list : bool
            Whether to return as a list of strings vs single string. Pretty much for internal usage.
        """
        text = [f"{self.common_name} results:", f"File: {self.image.truncated_path}"]
        if self.low_contrast_rois:
            text += [
                f"Median Contrast: {np.median([roi.contrast for roi in self.low_contrast_rois]):2.2f}",
                f"Median CNR: {np.median([roi.contrast_to_noise for roi in self.low_contrast_rois]):2.1f}",
                f'# Low contrast ROIs "seen": {sum(roi.passed_visibility for roi in self.low_contrast_rois):2.0f} of {len(self.low_contrast_rois)}',
                f"Area: {self.phantom_area:2.2f} mm^2",
            ]
        if self.high_contrast_rois:
            text += [
                f"MTF 80% (lp/mm): {self.mtf.relative_resolution(80):2.2f}",
                f"MTF 50% (lp/mm): {self.mtf.relative_resolution(50):2.2f}",
                f"MTF 30% (lp/mm): {self.mtf.relative_resolution(30):2.2f}",
            ]
        if not as_list:
            text = "\n".join(text)
        return text

    def _generate_results_data(self) -> PlanarResult:
        data = PlanarResult(
            analysis_type=self.common_name,
            median_contrast=np.median([roi.contrast for roi in self.low_contrast_rois]),
            median_cnr=np.median(
                [roi.contrast_to_noise for roi in self.low_contrast_rois]
            ),
            num_contrast_rois_seen=sum(
                roi.passed_visibility for roi in self.low_contrast_rois
            ),
            phantom_center_x_y=(self.phantom_center.x, self.phantom_center.y),
            low_contrast_rois=[roi.as_dict() for roi in self.low_contrast_rois],
            percent_integral_uniformity=self.percent_integral_uniformity(),
            phantom_area=self.phantom_area,
        )

        if self.mtf is not None:
            data.mtf_lp_mm = [
                {p: self.mtf.relative_resolution(p)} for p in (80, 50, 30)
            ]
        return data

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
        canvas = pdf.PylinacCanvas(
            filename,
            page_title=f"{self.common_name} Phantom Analysis",
            metadata=metadata,
            logo=logo,
        )

        # write the text/numerical values
        text = self.results(as_list=True)
        canvas.add_text(text=text, location=(1.5, 25), font_size=14)
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 5.5), font_size=12)
            canvas.add_text(text=notes, location=(1, 5))

        # plot the image
        data = io.BytesIO()
        self.save_analyzed_image(
            data, image=True, low_contrast=False, high_contrast=False
        )
        canvas.add_image(data, location=(1, 3.5), dimensions=(19, 19))
        # plot the high contrast
        if self.high_contrast_rois:
            canvas.add_new_page()
            data = io.BytesIO()
            self.save_analyzed_image(
                data, image=False, low_contrast=False, high_contrast=True
            )
            canvas.add_image(data, location=(1, 7), dimensions=(19, 19))
        # plot the low contrast
        if self.low_contrast_rois:
            canvas.add_new_page()
            data = io.BytesIO()
            self.save_analyzed_image(
                data, image=False, low_contrast=True, high_contrast=False
            )
            canvas.add_image(data, location=(1, 7), dimensions=(19, 19))

        canvas.finish()
        if open_file:
            webbrowser.open(filename)

    @property
    def phantom_center(self) -> Point:
        return (
            Point(self._center_override)
            if self._center_override is not None
            else self._phantom_center_calc()
        )

    @property
    def phantom_radius(self) -> float:
        return (
            self._size_override
            if self._size_override is not None
            else self._phantom_radius_calc()
        )

    @property
    def phantom_angle(self) -> float:
        return (
            self._angle_override
            if self._angle_override is not None
            else self._phantom_angle_calc()
        )

    @property
    def phantom_area(self) -> float:
        """The area of the detected ROI in mm^2"""
        area_px = self._create_phantom_outline_object()[0].area
        return area_px / self.image.dpmm**2

    def _phantom_center_calc(self):
        return bbox_center(self.phantom_ski_region)

    def _phantom_angle_calc(self):
        pass

    def _phantom_radius_calc(self):
        return math.sqrt(self.phantom_ski_region.bbox_area)

    def _find_ssd(self):
        """If the SSD parameter is set to auto, search at SAD, then at -5cm SID"""
        if isinstance(self._ssd, str) and self._ssd.lower() == "auto":
            self._ssd = self.image.metadata.get("RadiationMachineSAD", 1000)
            try:
                # cached property; no error means it found it.
                self.phantom_ski_region
            except ValueError:
                # 5cm up from SID
                self._ssd = self.image.metadata.get("RTImageSID", 1500) - 50
                self.phantom_ski_region


class LightRadResult(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    field_size_x_mm: float  #:
    field_size_y_mm: float  #:
    field_epid_offset_x_mm: float  #:
    field_epid_offset_y_mm: float  #:
    field_bb_offset_x_mm: float  #:
    field_bb_offset_y_mm: float  #:


class StandardImagingFC2(ImagePhantomBase):
    common_name = "SI FC-2"
    _demo_filename = "fc2.dcm"
    # these positions are the offset in mm from the center of the image to the nominal position of the BBs
    bb_positions_10x10 = {
        "TL": [-40, -40],
        "BL": [-40, 40],
        "TR": [40, -40],
        "BR": [40, 40],
    }
    bb_positions_15x15 = {
        "TL": [-65, -65],
        "BL": [-65, 65],
        "TR": [65, -65],
        "BR": [65, 65],
    }
    bb_sampling_box_size_mm = 10
    field_strip_width_mm = 5
    bb_size_mm = 4
    bb_edge_threshold_mm: float
    bb_centers: dict[str, Point]

    @staticmethod
    def run_demo() -> None:
        """Run the Standard Imaging FC-2 phantom analysis demonstration."""
        fc2 = StandardImagingFC2.from_demo_image()
        fc2.analyze()
        fc2.plot_analyzed_image()

    def analyze(
        self, invert: bool = False, fwxm: int = 50, bb_edge_threshold_mm: float = 10
    ) -> None:
        """Analyze the FC-2 phantom to find the BBs and the open field and compare to each other as well as the EPID.

        Parameters
        ----------
        invert : bool
            Whether to force-invert the image from the auto-detected inversion.
        fwxm : int
            The FWXM value to use to detect the field. For flattened fields, the default of 50 should be fine.
            For FFF fields, consider using a lower value such as 25-30.
        bb_edge_threshold_mm : float
            The threshold in mm to use to determine if the BB is near the edge of the image. If the BB is within this threshold,
            a different algorithm is used to determine the BB position that is more robust to edge effects but
            can give uncertainty when in a flat region (i.e. away from the field edge).
        """
        self.bb_edge_threshold_mm = bb_edge_threshold_mm
        self._check_inversion()
        if invert:
            self.image.invert()
        (
            self.field_center,
            self.field_width_x,
            self.field_width_y,
        ) = self._find_field_info(fwxm=fwxm)
        self.bb_center = self._find_overall_bb_centroid(fwxm=fwxm)
        self.epid_center = self.image.center

    def results(self, as_list: bool = False) -> str | list[str]:
        """Return the results of the analysis."""
        text = [
            f"{self.common_name} results:",
            f"File: {self.image.truncated_path}",
            f"The detected inplane field size was {self.field_width_y:2.1f}mm",
            f"The detected crossplane field size was {self.field_width_x:2.1f}mm",
            f"The inplane field was {self.field_epid_offset_mm.y:2.1f}mm from the EPID CAX",
            f"The crossplane field was {self.field_epid_offset_mm.x:2.1f}mm from the EPID CAX",
            f"The inplane field was {self.field_bb_offset_mm.y:2.1f}mm from the BB inplane center",
            f"The crossplane field was {self.field_bb_offset_mm.x:2.1f}mm from the BB crossplane center",
        ]
        if as_list:
            return text
        else:
            text = "\n".join(text)
            return text

    @property
    def field_epid_offset_mm(self) -> Vector:
        """Field offset from CAX using vector difference"""
        return (
            self.epid_center.as_vector() - self.field_center.as_vector()
        ) / self.image.dpmm

    @property
    def field_bb_offset_mm(self) -> Vector:
        """Field offset from BB centroid using vector difference"""
        return (self.bb_center - self.field_center) / self.image.dpmm

    def _generate_results_data(self) -> LightRadResult:
        """Return the results as a dict or dataclass"""
        return LightRadResult(
            field_size_x_mm=self.field_width_x,
            field_size_y_mm=self.field_width_y,
            field_epid_offset_x_mm=self.field_epid_offset_mm.x,
            field_epid_offset_y_mm=self.field_epid_offset_mm.y,
            field_bb_offset_x_mm=self.field_bb_offset_mm.x,
            field_bb_offset_y_mm=self.field_bb_offset_mm.y,
        )

    def _check_inversion(self):
        """Perform a normal corner-check inversion. Since these are always 10x10 or 15x15 fields it seems unlikely the corners will be exposed."""
        self.image.check_inversion()

    def _find_field_info(self, fwxm: int) -> (Point, float, float):
        """Determine the center and field widths of the detected field by sampling a strip through the center of the image in inplane and crossplane"""
        sample_width = self.field_strip_width_mm / 2 * self.image.dpmm
        # sample the strip (nominally 5mm) centered about the image center. Average the strip to reduce noise.
        x_bounds = (
            int(self.image.center.x - sample_width),
            int(self.image.center.x + sample_width),
        )
        y_img = np.mean(self.image[:, x_bounds[0] : x_bounds[1]], 1)
        y_prof = FWXMProfilePhysical(
            values=y_img,
            dpmm=self.image.dpmm,
            normalization=Normalization.BEAM_CENTER,
            ground=True,
            fwxm_height=fwxm,
        )
        y = y_prof.center_idx
        field_width_y = y_prof.field_width_mm
        y_bounds = (
            int(self.image.center.y - sample_width),
            int(self.image.center.y + sample_width),
        )
        x_img = np.mean(self.image[y_bounds[0] : y_bounds[1], :], 0)
        x_prof = FWXMProfilePhysical(
            values=x_img,
            dpmm=self.image.dpmm,
            normalization=Normalization.BEAM_CENTER,
            ground=True,
            fwxm_height=fwxm,
        )
        x = x_prof.center_idx
        field_width_x = x_prof.field_width_mm
        return Point(x=x, y=y), field_width_x, field_width_y

    def _find_overall_bb_centroid(self, fwxm: int) -> Point:
        """Determine the geometric center of the 4 BBs"""
        self.bb_centers = bb_centers = self._detect_bb_centers(fwxm)
        central_x = np.mean([p.x for p in bb_centers.values()])
        central_y = np.mean([p.y for p in bb_centers.values()])
        return Point(x=central_x, y=central_y)

    def _detect_bb_centers(self, fwxm: int) -> dict:
        """Sample a 10x10mm square about each BB to detect it. Adjustable using self.bb_sampling_box_size_mm"""
        bb_positions = {}
        nominal_positions = self._determine_bb_set(fwxm=fwxm)
        dpmm = self.image.dpmm
        self.image.filter(size=3, kind="median")
        # sample the square, use skimage to find the ROI weighted centroid of the BBs
        for key, position in nominal_positions.items():
            near_edge = self._is_bb_near_edge(bb_position=position)
            if near_edge:
                # we apply a local histogram equalizer.
                # this is local contrast enhancer. It helps the BBs stand out from the background
                # and sharpen the gap between the BB and field edge.
                # we need to replace and reset the array since we're in the loop
                # This is a bit of BMF.
                original_array = np.copy(self.image.array)
                bb_radius_px = self.bb_size_mm / 2 * dpmm
                self.image.array = exposure.equalize_adapthist(
                    self.image.array, kernel_size=int(round(bb_radius_px * 2))
                )
                self.image.filter(size=3, kind="median")

            # now find the weighted centroid of the BB
            points = self.image.compute(
                SizedDiskLocator.from_center_physical(
                    expected_position_mm=position,
                    search_window_mm=(
                        self.bb_sampling_box_size_mm,
                        self.bb_sampling_box_size_mm,
                    ),
                    radius_mm=self.bb_size_mm / 2,
                    radius_tolerance_mm=self.bb_size_mm / 2,
                )
            )
            # if we applied the local histogram equalizer, revert the image back to normal for the next cycle
            if near_edge:
                self.image.array = original_array
            bb_positions[key] = points[0]
        return bb_positions

    def _determine_bb_set(self, fwxm: int) -> dict:
        """This finds the approximate field size to determine whether to check for the 10x10 BBs or the 15x15. Returns the BB positions"""
        if not np.allclose(self.field_width_x, self.field_width_y, atol=10):
            raise ValueError(
                f"The detected y and x field sizes were too different from one another. They should be within 1cm from each other. Detected field sizes: x={self.field_width_x:.2f}mm, y={self.field_width_y:.2f}mm"
            )
        if self.field_width_x > 140:
            return self.bb_positions_15x15
        else:
            return self.bb_positions_10x10

    def plot_analyzed_image(
        self, show: bool = True, **kwargs
    ) -> tuple[list[plt.Figure], list[str]]:
        """Plot the analyzed image.

        Parameters
        ----------
        show : bool
            Whether to actually show the image when called.
        """
        figs = []
        names = []
        fig, axes = plt.subplots(1)
        figs.append(fig)
        names.append("image")
        self.image.plot(ax=axes, show=False, metric_kwargs={"color": "g"}, **kwargs)
        axes.axis("off")
        axes.set_title(f"{self.common_name} Phantom Analysis")

        # plot the bb center as small lines
        axes.axhline(
            y=self.bb_center.y, color="g", xmin=0.25, xmax=0.75, label="BB Centroid"
        )
        axes.axvline(x=self.bb_center.x, color="g", ymin=0.25, ymax=0.75)

        # plot the epid center as image-sized lines
        axes.axhline(y=self.epid_center.y, color="b", label="EPID Center")
        axes.axvline(x=self.epid_center.x, color="b")

        # plot the field center as field-sized lines
        axes.axhline(
            y=self.field_center.y,
            xmin=0.15,
            xmax=0.85,
            color="red",
            label="Field Center",
        )
        axes.axvline(x=self.field_center.x, ymin=0.15, ymax=0.85, color="red")

        axes.legend()

        if show:
            plt.show()
        return figs, names

    def save_analyzed_image(
        self,
        filename: None | str | BinaryIO = None,
        to_streams: bool = False,
        **kwargs,
    ) -> dict[str, BinaryIO] | list[str] | None:
        """Save the analyzed image to disk or to stream. Kwargs are passed to plt.savefig()

        Parameters
        ----------
        filename : None, str, stream
            A string representing where to save the file to. If split_plots and to_streams are both true, leave as None as newly-created streams are returned.
        to_streams: bool
            This only matters if split_plots is True. If both of these are true, multiple streams will be created and returned as a dict.
        """
        if filename is None and to_streams is False:
            raise ValueError("Must pass in a filename unless saving to streams.")
        figs, names = self.plot_analyzed_image(show=False, **kwargs)
        if not to_streams:
            plt.savefig(filename, **kwargs)
        else:
            # append names to filename if it's file-like
            if not to_streams:
                filenames = []
                f, ext = osp.splitext(filename)
                for name in names:
                    filenames.append(f + "_" + name + ext)
            else:  # it's a stream buffer
                filenames = [io.BytesIO() for _ in names]
            for fig, name in zip(figs, filenames):
                fig.savefig(name, **kwargs)
            if to_streams:
                return {name: stream for name, stream in zip(names, filenames)}

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
        canvas = pdf.PylinacCanvas(
            filename,
            page_title=f"{self.common_name} Phantom Analysis",
            metadata=metadata,
            logo=logo,
        )

        # write the text/numerical values
        text = self.results(as_list=True)
        canvas.add_text(text=text, location=(1.5, 25), font_size=14)
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 5.5), font_size=12)
            canvas.add_text(text=notes, location=(1, 5))

        data = io.BytesIO()
        self.save_analyzed_image(data)
        canvas.add_image(data, location=(1, 3.5), dimensions=(19, 19))

        canvas.finish()
        if open_file:
            webbrowser.open(filename)

    def _is_bb_near_edge(self, bb_position: [float, float]) -> bool:
        """Check if the center of the nominal BB set is < threshold from the field edge."""
        half_width_x = self.field_width_x / 2
        half_width_y = self.field_width_y / 2
        threshold = self.bb_edge_threshold_mm
        # Check proximity to left or right edge
        is_near_horizontal_edge = abs(bb_position[0]) > half_width_x - threshold
        # Check proximity to top or bottom edge
        is_near_vertical_edge = abs(bb_position[1]) > half_width_y - threshold
        return is_near_horizontal_edge or is_near_vertical_edge


class IMTLRad(StandardImagingFC2):
    """The IMT light/rad phantom: https://www.imtqa.com/products/l-rad"""

    common_name = "IMT L-Rad"
    _demo_filename = "imtlrad.dcm"
    center_only_bb = {"Center": [0, 0]}
    bb_sampling_box_size_mm = 12
    field_strip_width_mm = 5
    bb_size_mm = 3

    def _determine_bb_set(self, fwxm: int) -> dict:
        return self.center_only_bb


class DoselabRLf(StandardImagingFC2):
    """The Doselab light/rad phantom"""

    common_name = "Doselab RLf"
    _demo_filename = "Doselab_RLf.dcm"
    # these positions are the offset in mm from the center of the image to the nominal position of the BBs
    bb_positions_10x10 = {
        "TL": [-17, -45],
        "BL": [-45, 17],
        "TR": [45, -17],
        "BR": [17, 45],
    }
    # 15x15 is not as robust as 10x10
    # bb_positions_15x15 = {
    #     "TL": [-45, -70],
    #     "BL": [-70, 45],
    #     "TR": [70, -45],
    #     "BR": [45, 70],
    # }

    def _determine_bb_set(self, fwxm: int) -> dict:
        return self.bb_positions_10x10

    @staticmethod
    def run_demo() -> None:
        """Run the Doselab RFl phantom analysis demonstration."""
        dl = DoselabRLf.from_demo_image()
        dl.analyze()
        dl.plot_analyzed_image()


class IsoAlign(StandardImagingFC2):
    """The PTW Iso-Align light/rad phantom"""

    common_name = "PTW Iso-Align"
    _demo_filename = "ptw_isoalign.dcm"
    # these positions are the offset in mm from the center of the image to the nominal position of the BBs
    bb_positions = {
        "Center": [0, 0],
        "Top": [0, -25],
        "Bottom": [0, 25],
        "Left": [-25, 0],
        "Right": [25, 0],
    }
    field_strip_width_mm = 10

    def _determine_bb_set(self, fwxm: int) -> dict:
        return self.bb_positions

    @staticmethod
    def run_demo() -> None:
        """Run the phantom analysis demonstration."""
        al = IsoAlign.from_demo_image()
        al.analyze()
        al.plot_analyzed_image()


class SNCFSQA(StandardImagingFC2):
    """SNC light/rad phantom. See the 'FSQA' phantom and specs: https://www.sunnuclear.com/products/suncheck-machine.

    Unlike other light/rad phantoms, this does not have at least a centered BB. The edge markers are in the penumbra
    and thus detecting them is difficult. We thus detect the one offset marker in the top right of the image.
    This is offset by 4cm in each direction. We can then assume that the phantom center is -4cm from this point,
    creating a 'virtual center' so we have an apples-to-apples comparison.
    """

    common_name = "SNC FSQA"
    _demo_filename = "FSQA_15x15.dcm"
    center_only_bb = {"TR": [40, -40]}
    # bb_sampling_box_size_mm = 8
    field_strip_width_mm = 5

    def _determine_bb_set(self, fwxm: int) -> dict:
        return self.center_only_bb

    def _find_overall_bb_centroid(self, fwxm: int) -> Point:
        """Determine the geometric center of the 4 BBs"""
        # detect the upper right BB
        self.bb_centers = self._detect_bb_centers(fwxm)
        # add another virtual bb at the center of the phantom, knowing it's offset by 4cm in each direction
        self.bb_centers["Virtual Center"] = self.bb_centers["TR"] - Point(
            40 * self.image.dpmm, -40 * self.image.dpmm
        )
        return self.bb_centers["Virtual Center"]


class LasVegas(ImagePhantomBase):
    _demo_filename = "lasvegas.dcm"
    common_name = "Las Vegas"
    phantom_bbox_size_mm2 = 20260
    detection_conditions = [is_centered, is_right_size]
    phantom_outline_object = {"Rectangle": {"width ratio": 0.62, "height ratio": 0.62}}
    low_contrast_background_roi_settings = {
        "roi 1": {"distance from center": 0.24, "angle": 0, "roi radius": 0.03},
        "roi 2": {"distance from center": 0.24, "angle": 90, "roi radius": 0.03},
        "roi 3": {"distance from center": 0.24, "angle": 180, "roi radius": 0.03},
        "roi 4": {"distance from center": 0.24, "angle": 270, "roi radius": 0.03},
    }
    low_contrast_roi_settings = {
        "roi 1": {"distance from center": 0.107, "angle": 0.5, "roi radius": 0.028},
        "roi 2": {"distance from center": 0.141, "angle": 39.5, "roi radius": 0.028},
        "roi 3": {"distance from center": 0.205, "angle": 58, "roi radius": 0.028},
        "roi 4": {"distance from center": 0.179, "angle": -76.5, "roi radius": 0.016},
        "roi 5": {"distance from center": 0.095, "angle": -63.5, "roi radius": 0.016},
        "roi 6": {"distance from center": 0.042, "angle": 0.5, "roi radius": 0.016},
        "roi 7": {"distance from center": 0.097, "angle": 65.5, "roi radius": 0.016},
        "roi 8": {"distance from center": 0.178, "angle": 76.5, "roi radius": 0.016},
        "roi 9": {"distance from center": 0.174, "angle": -97.5, "roi radius": 0.012},
        "roi 10": {"distance from center": 0.088, "angle": -105.5, "roi radius": 0.012},
        "roi 11": {"distance from center": 0.024, "angle": -183.5, "roi radius": 0.012},
        "roi 12": {"distance from center": 0.091, "angle": 105.5, "roi radius": 0.012},
        "roi 13": {"distance from center": 0.179, "angle": 97.5, "roi radius": 0.012},
        "roi 14": {"distance from center": 0.189, "angle": -113.5, "roi radius": 0.007},
        "roi 15": {"distance from center": 0.113, "angle": -131.5, "roi radius": 0.007},
        "roi 16": {
            "distance from center": 0.0745,
            "angle": -181.5,
            "roi radius": 0.007,
        },
        "roi 17": {"distance from center": 0.115, "angle": 130, "roi radius": 0.007},
        "roi 18": {"distance from center": 0.191, "angle": 113, "roi radius": 0.007},
        "roi 19": {
            "distance from center": 0.2085,
            "angle": -124.6,
            "roi radius": 0.003,
        },
        "roi 20": {"distance from center": 0.146, "angle": -144.3, "roi radius": 0.003},
    }

    @staticmethod
    def run_demo():
        """Run the Las Vegas phantom analysis demonstration."""
        lv = LasVegas.from_demo_image()
        lv.analyze()
        lv.plot_analyzed_image()

    def _preprocess(self):
        self._check_direction()

    def _check_inversion(self):
        """Check the inversion by using the histogram of the phantom region"""
        roi = self.phantom_ski_region
        phantom_array = self.image.array[
            roi.bbox[0] : roi.bbox[2], roi.bbox[1] : roi.bbox[3]
        ]
        phantom_sub_image = image.load(phantom_array)
        phantom_sub_image.crop(int(phantom_sub_image.shape[0] * 0.1))
        p5 = np.percentile(phantom_sub_image, 0.5)
        p50 = np.percentile(phantom_sub_image, 50)
        p95 = np.percentile(phantom_sub_image, 99.5)
        dist_to_5 = abs(p50 - p5)
        dist_to_95 = abs(p50 - p95)
        if dist_to_5 > dist_to_95:
            self.image.invert()

    def _check_direction(self) -> None:
        """Check that the phantom is facing the right direction and if not perform a left-right flip of the array."""
        circle = CollapsedCircleProfile(
            self.phantom_center,
            self.phantom_radius * 0.175,
            self.image,
            ccw=False,
            width_ratio=0.16,
            num_profiles=5,
        )
        roll_amount = np.where(circle.values == circle.values.min())[0][0]
        circle.roll(roll_amount)
        circle.filter(size=0.015, kind="median")
        valley_idxs, _ = circle.find_peaks(max_number=2)
        if valley_idxs[0] > valley_idxs[1]:
            self.image.array = np.fliplr(self.image.array)
            self._phantom_ski_region = None

    def _phantom_radius_calc(self) -> float:
        return math.sqrt(self.phantom_ski_region.bbox_area) * 1.626

    def _phantom_angle_calc(self) -> float:
        return 0.0

    def _plot_lowcontrast_graph(self, axes: plt.Axes):
        """Plot the low contrast ROIs to an axes, including visibility"""
        # plot contrast
        (line1,) = axes.plot(
            [roi.contrast for roi in self.low_contrast_rois],
            marker="o",
            color="m",
            label="Contrast",
        )
        axes.axhline(self._low_contrast_threshold, color="m")
        axes.grid(True)
        axes.set_title("Low-frequency Contrast")
        axes.set_xlabel("ROI #")
        axes.set_ylabel("Contrast")
        # plot CNR
        axes2 = axes.twinx()
        axes2.set_ylabel("CNR")
        (line2,) = axes2.plot(
            [roi.contrast_to_noise for roi in self.low_contrast_rois],
            marker="^",
            label="CNR",
        )
        # plot visibility; here's what different from the base method
        axes3 = axes.twinx()
        axes3.set_ylabel("Visibility")
        (line3,) = axes3.plot(
            [roi.visibility for roi in self.low_contrast_rois],
            marker="*",
            color="blue",
            label="Visibility",
        )
        axes3.axhline(self.visibility_threshold, color="blue")
        axes3.spines.right.set_position(("axes", 1.2))

        axes.legend(handles=[line1, line2, line3])

    def results(self, as_list: bool = False) -> str | list[str]:
        """Return the results of the analysis. Overridden because ROIs seen is based on visibility, not CNR.

        Parameters
        ----------
        as_list : bool
            Whether to return as a list of strings vs single string. Pretty much for internal usage.
        """
        text = [f"{self.common_name} results:", f"File: {self.image.truncated_path}"]
        text += [
            f"Median Contrast: {np.median([roi.contrast for roi in self.low_contrast_rois]):2.2f}",
            f"Median CNR: {np.median([roi.contrast_to_noise for roi in self.low_contrast_rois]):2.1f}",
            f'# Low contrast ROIs "seen": {sum(roi.passed_visibility for roi in self.low_contrast_rois):2.0f} of {len(self.low_contrast_rois)}',
        ]
        if not as_list:
            text = "\n".join(text)
        return text

    def _generate_results_data(self) -> PlanarResult:
        """Overridden because ROIs seen is based on visibility, not CNR"""
        return PlanarResult(
            analysis_type=self.common_name,
            median_contrast=np.median([roi.contrast for roi in self.low_contrast_rois]),
            median_cnr=np.median(
                [roi.contrast_to_noise for roi in self.low_contrast_rois]
            ),
            num_contrast_rois_seen=sum(
                roi.passed_visibility for roi in self.low_contrast_rois
            ),
            phantom_center_x_y=(self.phantom_center.x, self.phantom_center.y),
            low_contrast_rois=[r.as_dict() for r in self.low_contrast_rois],
            percent_integral_uniformity=self.percent_integral_uniformity(),
            phantom_area=self.phantom_area,
        )


class ElektaLasVegas(LasVegas):
    """Elekta's variant of the Las Vegas."""

    _demo_filename = "elekta_las_vegas.dcm"
    common_name = "Elekta Las Vegas"
    phantom_bbox_size_mm2 = 140 * 140
    phantom_outline_object = {"Rectangle": {"width ratio": 0.61, "height ratio": 0.61}}
    low_contrast_background_roi_settings = {
        "roi 1": {"distance from center": 0.24, "angle": 0, "roi radius": 0.03},
        "roi 2": {"distance from center": 0.24, "angle": 90, "roi radius": 0.03},
        "roi 3": {"distance from center": 0.24, "angle": 180, "roi radius": 0.03},
        "roi 4": {"distance from center": 0.24, "angle": 270, "roi radius": 0.03},
    }
    low_contrast_roi_settings = {
        "roi 1": {"distance from center": 0.161, "angle": 0.4, "roi radius": 0.024},
        "roi 2": {"distance from center": 0.181, "angle": 28.6, "roi radius": 0.024},
        "roi 3": {"distance from center": 0.238, "angle": 47.45, "roi radius": 0.024},
        "roi 4": {"distance from center": 0.183, "angle": -70.6, "roi radius": 0.015},
        "roi 5": {"distance from center": 0.107, "angle": -55.1, "roi radius": 0.015},
        "roi 6": {"distance from center": 0.061, "angle": 1, "roi radius": 0.015},
        "roi 7": {"distance from center": 0.107, "angle": 55.15, "roi radius": 0.015},
        "roi 8": {"distance from center": 0.185, "angle": 71.1, "roi radius": 0.015},
        "roi 9": {"distance from center": 0.175, "angle": -97.3, "roi radius": 0.011},
        "roi 10": {"distance from center": 0.09, "angle": -104.3, "roi radius": 0.011},
        "roi 11": {"distance from center": 0.022, "angle": -180, "roi radius": 0.011},
        "roi 12": {"distance from center": 0.088, "angle": 104.6, "roi radius": 0.011},
        "roi 13": {"distance from center": 0.1757, "angle": 97.26, "roi radius": 0.011},
        "roi 14": {
            "distance from center": 0.1945,
            "angle": -116.58,
            "roi radius": 0.006,
        },
        "roi 15": {
            "distance from center": 0.124,
            "angle": -135.11,
            "roi radius": 0.006,
        },
        "roi 16": {
            "distance from center": 0.0876,
            "angle": 179.85,
            "roi radius": 0.006,
        },
        "roi 17": {"distance from center": 0.1227, "angle": 135.4, "roi radius": 0.006},
        "roi 18": {
            "distance from center": 0.1947,
            "angle": 116.65,
            "roi radius": 0.006,
        },
        "roi 19": {
            "distance from center": 0.2258,
            "angle": -129.53,
            "roi radius": 0.003,
        },
        "roi 20": {
            "distance from center": 0.1699,
            "angle": -148.57,
            "roi radius": 0.003,
        },
        "roi 21": {
            "distance from center": 0.145,
            "angle": -179.82,
            "roi radius": 0.003,
        },
        "roi 22": {"distance from center": 0.1682, "angle": 149, "roi radius": 0.003},
    }

    @staticmethod
    def run_demo():
        """Run the Elekta Las Vegas phantom analysis demonstration."""
        lv = ElektaLasVegas.from_demo_image()
        lv.image.rot90(n=3)
        lv.analyze()
        lv.plot_analyzed_image()


class PTWEPIDQC(ImagePhantomBase):
    _demo_filename = "PTW-EPID-QC.dcm"
    common_name = "PTW EPID QC"
    phantom_bbox_size_mm2 = 250**2
    detection_conditions = [is_centered, is_right_size]
    detection_canny_settings = {"sigma": 4, "percentiles": (0.001, 0.01)}
    phantom_outline_object = {"Rectangle": {"width ratio": 8.55, "height ratio": 8.55}}
    high_contrast_roi_settings = {
        # angled rois
        "roi 1": {
            "distance from center": 1.5,
            "angle": -135,
            "roi radius": 0.35,
            "lp/mm": 0.15,
        },
        "roi 2": {
            "distance from center": 3.1,
            "angle": -109,
            "roi radius": 0.35,
            "lp/mm": 0.21,
        },
        "roi 3": {
            "distance from center": 3.4,
            "angle": -60,
            "roi radius": 0.3,
            "lp/mm": 0.27,
        },
        "roi 4": {
            "distance from center": 1.9,
            "angle": -60,
            "roi radius": 0.25,
            "lp/mm": 0.33,
        },
        # vertical rois
        "roi 5": {
            "distance from center": 3.68,
            "angle": -90,
            "roi radius": 0.18,
            "lp/mm": 0.5,
        },
        "roi 6": {
            "distance from center": 2.9,
            "angle": -90,
            "roi radius": 0.08,
            "lp/mm": 2,
        },
        "roi 7": {
            "distance from center": 2.2,
            "angle": -90,
            "roi radius": 0.04,
            "lp/mm": 3,
        },
    }
    low_contrast_roi_settings = {
        "roi 1": {"distance from center": 3.87, "angle": 31, "roi radius": 0.3},
        "roi 2": {"distance from center": 3.48, "angle": 17, "roi radius": 0.3},
        "roi 3": {"distance from center": 3.3, "angle": 0, "roi radius": 0.3},
        "roi 4": {"distance from center": 3.48, "angle": -17, "roi radius": 0.3},
        "roi 5": {"distance from center": 3.87, "angle": -31, "roi radius": 0.3},
        "roi 6": {"distance from center": 3.87, "angle": 180 - 31, "roi radius": 0.3},
        "roi 7": {"distance from center": 3.48, "angle": 180 - 17, "roi radius": 0.3},
        "roi 8": {"distance from center": 3.3, "angle": 180, "roi radius": 0.3},
        "roi 9": {"distance from center": 3.48, "angle": 180 + 17, "roi radius": 0.3},
    }
    low_contrast_background_roi_settings = {
        "roi 1": {"distance from center": 3.85, "angle": -148, "roi radius": 0.3},
    }

    @staticmethod
    def run_demo() -> None:
        """Run the Standard Imaging QC-3 phantom analysis demonstration."""
        ptw = PTWEPIDQC.from_demo_image()
        ptw.analyze()
        ptw.plot_analyzed_image()

    def _phantom_radius_calc(self) -> float:
        """The radius of the phantom in pixels; the value itself doesn't matter, it's just
        used for relative distances to ROIs.

        Returns
        -------
        radius : float
        """
        return math.sqrt(self.phantom_ski_region.bbox_area) * 0.116

    def _phantom_angle_calc(self) -> float:
        """The angle of the phantom. This assumes the user has placed the phantom with the high-contrast line pairs at the top
        and low contrast at the bottom at a fixed rotation of 0 degrees.

        Returns
        -------
        angle : float
            The angle in degrees.
        """
        return 0

    def _check_inversion(self):
        """The pixels inside the phantom should be mostly bright/high. If not, invert"""
        roi = self.phantom_ski_region
        phantom_array = self.image.array[
            roi.bbox[0] : roi.bbox[2], roi.bbox[1] : roi.bbox[3]
        ]
        p5, p50, p95 = np.percentile(phantom_array, [2, 50, 98])
        if abs(p50 - p5) < abs(p50 - p95):
            self.image.invert()


class IBAPrimusA(ImagePhantomBase):
    common_name = "IBA Primus A"
    _demo_filename = "iba_primus.dcm"
    phantom_bbox_size_mm2 = (
        15**2
    )  # with the Primus, we only search for the central crosshair
    detection_conditions = [is_centered, is_right_size, is_square]
    phantom_outline_object = {
        "Rectangle": {"width ratio": 10.75, "height ratio": 10.75}
    }
    high_contrast_roi_settings = {
        "roi 1": {
            "distance from center": 5.19,
            "angle": 86.65,
            "roi radius": 0.12,
            "lp/mm": 0.6,
        },
        "roi 2": {
            "distance from center": 4.92,
            "angle": 89.5,
            "roi radius": 0.1,
            "lp/mm": 0.7,
        },
        "roi 3": {
            "distance from center": 4.68,
            "angle": 92.3,
            "roi radius": 0.09,
            "lp/mm": 0.8,
        },
        "roi 4": {
            "distance from center": 4.45,
            "angle": 95.4,
            "roi radius": 0.08,
            "lp/mm": 0.9,
        },
        "roi 5": {
            "distance from center": 4.23,
            "angle": 99.5,
            "roi radius": 0.07,
            "lp/mm": 1,
        },
        "roi 6": {
            "distance from center": 4.07,
            "angle": 102.7,
            "roi radius": 0.06,
            "lp/mm": 1.2,
        },
        "roi 7": {
            "distance from center": 3.92,
            "angle": 105.73,
            "roi radius": 0.05,
            "lp/mm": 1.4,
        },
        "roi 8": {
            "distance from center": 3.82,
            "angle": 108.65,
            "roi radius": 0.04,
            "lp/mm": 1.6,
        },
        "roi 9": {
            "distance from center": 4.59,
            "angle": 74.4,
            "roi radius": 0.04,
            "lp/mm": 1.8,
        },
        "roi 10": {
            "distance from center": 4.4,
            "angle": 76.2,
            "roi radius": 0.035,
            "lp/mm": 2.0,
        },
        "roi 11": {
            "distance from center": 4.19,
            "angle": 77.77,
            "roi radius": 0.03,
            "lp/mm": 2.2,
        },
        "roi 12": {
            "distance from center": 4,
            "angle": 79.6,
            "roi radius": 0.03,
            "lp/mm": 2.5,
        },
        "roi 13": {
            "distance from center": 3.67,
            "angle": 83.1,
            "roi radius": 0.025,
            "lp/mm": 2.8,
        },
    }
    low_contrast_roi_settings = {
        "roi 1": {"distance from center": 3.95, "angle": 19, "roi radius": 0.15},
        "roi 2": {"distance from center": 3.95, "angle": 5, "roi radius": 0.15},
        "roi 3": {"distance from center": 3.95, "angle": -9, "roi radius": 0.15},
        "roi 4": {"distance from center": 3.95, "angle": -23, "roi radius": 0.15},
        "roi 5": {"distance from center": 3.95, "angle": -37, "roi radius": 0.15},
        "roi 6": {"distance from center": 3.95, "angle": -51, "roi radius": 0.15},
        "roi 7": {"distance from center": 3.95, "angle": -65, "roi radius": 0.15},
        "roi 8": {"distance from center": 3.95, "angle": -79, "roi radius": 0.15},
        "roi 9": {"distance from center": 3.95, "angle": -107, "roi radius": 0.15},
        "roi 10": {"distance from center": 3.95, "angle": -121, "roi radius": 0.15},
        "roi 11": {"distance from center": 3.95, "angle": -135, "roi radius": 0.15},
        "roi 12": {"distance from center": 3.95, "angle": -149, "roi radius": 0.15},
        "roi 13": {"distance from center": 3.95, "angle": -163, "roi radius": 0.15},
        "roi 14": {"distance from center": 3.95, "angle": -177, "roi radius": 0.15},
        "roi 15": {"distance from center": 3.95, "angle": -191, "roi radius": 0.15},
    }
    low_contrast_background_roi_settings = {
        "roi 1": {"distance from center": 3.95, "angle": -205, "roi radius": 0.15},
    }

    def _check_inversion(self):
        """The center region (crosshair) should be less intense than an area adjacent to it."""
        crosshair_disk = DiskROI(
            self.image.array,
            angle=0,
            roi_radius=self.phantom_radius / 2,
            dist_from_center=0,
            phantom_center=self.phantom_center,
        )
        adjacent_disk = DiskROI(
            self.image.array,
            angle=0,
            roi_radius=self.phantom_radius / 2,
            dist_from_center=self.phantom_radius,
            phantom_center=self.phantom_center,
        )
        if crosshair_disk.pixel_value < adjacent_disk.pixel_value:
            self.image.invert()

    def _wl_spread(self):
        """window/level spread based on low contrast ROI pixel values"""
        pixel_values = [roi.pixel_value for roi in self.low_contrast_rois]
        return abs(max(pixel_values) - min(pixel_values))

    def window_floor(self) -> float | None:
        return (
            min(roi.pixel_value for roi in self.low_contrast_rois) - self._wl_spread()
        )

    def window_ceiling(self) -> float | None:
        return (
            max(roi.pixel_value for roi in self.low_contrast_rois) + self._wl_spread()
        )

    @cached_property
    def phantom_angle(self) -> float:
        """Cache this; calculating the angle is expensive"""
        return super().phantom_angle

    def _phantom_angle_calc(self) -> float:
        """Fine-tune the angle by finding the two ends of the dynamic wedge steps and correcting"""
        prof = CollapsedCircleProfile(
            center=self.phantom_center,
            radius=self.phantom_radius * 4.37,
            image_array=self.image,
            start_angle=-np.pi / 2,
        )
        # get the points of max delta
        delta_array = np.argsort(np.diff(median_filter(prof.values, size=5)))
        # unfortunately, there may be several pixels of max gradient adjacent; we take the first
        # point and the next point that is not near the first
        first = delta_array[0]
        second = None
        one_degree = delta_array.size / 360
        for idx in delta_array:
            if first + one_degree < idx or idx < first - one_degree:
                second = idx
                break
        if not second:
            warnings.warn(
                "The phantom angle was not able to be fine-tuned; a default of 0 is being used instead. Ensure the image is not rotated."
            )
            return 0
        # now figure out the angle from the two deltas using the midpoint
        # perfect set up is when the midpoint is at 0.5. Use diff from that to get offset angle
        angle = (0.5 - ((second - first) / 2 + first) / prof.values.size) * 360
        near_cardinal = (-95 < angle < -85) or (85 < angle < 95) or (-5 < angle < 5)
        if near_cardinal:
            return angle
        else:
            # something is wrong; image likely rotated
            warnings.warn(
                "The phantom angle was not able to be fine-tuned; a default of 0 is being used instead. Ensure the image is not rotated."
            )
            return 0

    def _phantom_radius_calc(self):
        return math.sqrt(self.phantom_ski_region.area_bbox)

    @staticmethod
    def run_demo() -> None:
        """Run the Standard Imaging QC-3 phantom analysis demonstration."""
        primus = IBAPrimusA.from_demo_image()
        primus.analyze(ssd=1395)
        print(primus.results())
        primus.plot_analyzed_image()


class StandardImagingQC3(ImagePhantomBase):
    _demo_filename = "qc3.dcm"
    common_name = "SI QC-3"
    phantom_bbox_size_mm2 = 168**2
    detection_conditions = [is_centered, is_right_size]
    phantom_outline_object = {"Rectangle": {"width ratio": 7.5, "height ratio": 6}}
    high_contrast_roi_settings = {
        "roi 1": {
            "distance from center": 2.8,
            "angle": 0,
            "roi radius": 0.5,
            "lp/mm": 0.1,
        },
        "roi 2": {
            "distance from center": -2.8,
            "angle": 0,
            "roi radius": 0.5,
            "lp/mm": 0.2,
        },
        "roi 3": {
            "distance from center": 1.45,
            "angle": 0,
            "roi radius": 0.5,
            "lp/mm": 0.25,
        },
        "roi 4": {
            "distance from center": -1.45,
            "angle": 0,
            "roi radius": 0.5,
            "lp/mm": 0.45,
        },
        "roi 5": {
            "distance from center": 0,
            "angle": 0,
            "roi radius": 0.5,
            "lp/mm": 0.76,
        },
    }
    low_contrast_roi_settings = {
        "roi 1": {"distance from center": 2, "angle": -90, "roi radius": 0.5},
        "roi 2": {"distance from center": 2.4, "angle": 55, "roi radius": 0.5},
        "roi 3": {"distance from center": 2.4, "angle": -55, "roi radius": 0.5},
        "roi 4": {"distance from center": 2.4, "angle": 128, "roi radius": 0.5},
        "roi 5": {"distance from center": 2.4, "angle": -128, "roi radius": 0.5},
    }
    low_contrast_background_roi_settings = {
        "roi 1": {"distance from center": 2, "angle": 90, "roi radius": 0.5},
    }

    @classmethod
    def from_demo_image(cls):
        """Instantiate and load the demo image."""
        demo_file = retrieve_demo_file(name=cls._demo_filename)
        inst = cls(demo_file)
        inst.image.invert()
        return inst

    @staticmethod
    def run_demo() -> None:
        """Run the Standard Imaging QC-3 phantom analysis demonstration."""
        qc3 = StandardImagingQC3.from_demo_image()
        qc3.analyze()
        qc3.plot_analyzed_image()

    def _phantom_radius_calc(self) -> float:
        """The radius of the phantom in pixels; the value itself doesn't matter, it's just
        used for relative distances to ROIs.

        Returns
        -------
        radius : float
        """
        return math.sqrt(self.phantom_ski_region.bbox_area) * 0.0896

    @lru_cache()
    def _phantom_angle_calc(self) -> float:
        """The angle of the phantom. This assumes the user is using the stand that comes with the phantom,
        which angles the phantom at 45 degrees.

        Returns
        -------
        angle : float
            The angle in degrees.
        """
        angle = np.degrees(self.phantom_ski_region.orientation)
        if np.isclose(angle, 45, atol=5):
            return 45
        elif np.isclose(angle, -45, atol=5):
            return -45
        else:
            raise ValueError(
                "The phantom angle was not near +/-45 degrees. Please adjust the phantom."
            )


class StandardImagingQCkV(StandardImagingQC3):
    _demo_filename = "SI-QC-kV.dcm"
    common_name = "SI QC-kV"
    phantom_bbox_size_mm2 = 142**2
    detection_conditions = [is_centered, is_right_size]
    phantom_outline_object = {"Rectangle": {"width ratio": 7.8, "height ratio": 6.4}}
    high_contrast_roi_settings = {
        "roi 1": {
            "distance from center": 2.8,
            "angle": 0,
            "roi radius": 0.5,
            "lp/mm": 0.66,
        },
        "roi 2": {
            "distance from center": -2.8,
            "angle": 0,
            "roi radius": 0.5,
            "lp/mm": 0.98,
        },
        "roi 3": {
            "distance from center": 1.45,
            "angle": 0,
            "roi radius": 0.5,
            "lp/mm": 1.50,
        },
        "roi 4": {
            "distance from center": -1.45,
            "angle": 0,
            "roi radius": 0.5,
            "lp/mm": 2.00,
        },
        "roi 5": {
            "distance from center": 0,
            "angle": 0,
            "roi radius": 0.5,
            "lp/mm": 2.46,
        },
    }
    low_contrast_roi_settings = {
        "roi 1": {"distance from center": 2, "angle": -90, "roi radius": 0.5},
        "roi 2": {"distance from center": 2.4, "angle": 55, "roi radius": 0.5},
        "roi 3": {"distance from center": 2.4, "angle": -55, "roi radius": 0.5},
        "roi 4": {"distance from center": 2.4, "angle": 128, "roi radius": 0.5},
        "roi 5": {"distance from center": 2.4, "angle": -128, "roi radius": 0.5},
    }
    low_contrast_background_roi_settings = {
        "roi 1": {"distance from center": 2, "angle": 90, "roi radius": 0.5},
    }

    @staticmethod
    def run_demo() -> None:
        """Run the Standard Imaging QC-3 phantom analysis demonstration."""
        qc3 = StandardImagingQCkV.from_demo_image()
        qc3.analyze()
        qc3.plot_analyzed_image()

    def _phantom_radius_calc(self) -> float:
        """The radius of the phantom in pixels; the value itself doesn't matter, it's just
        used for relative distances to ROIs.

        Returns
        -------
        radius : float
        """
        return math.sqrt(self.phantom_ski_region.bbox_area) * 0.0989


class SNCkV(ImagePhantomBase):
    _demo_filename = "SNC-kV.dcm"
    common_name = "SNC kV-QA"
    phantom_bbox_size_mm2 = 134**2
    roi_match_condition = "closest"
    detection_conditions = [is_centered, is_right_size, is_square]
    phantom_outline_object = {"Rectangle": {"width ratio": 7.7, "height ratio": 5.6}}
    high_contrast_roi_settings = {
        "roi 1": {
            "distance from center": 1.8,
            "angle": 0,
            "roi radius": 0.7,
            "lp/mm": 0.6,
        },
        "roi 2": {
            "distance from center": -1.8,
            "angle": 90,
            "roi radius": 0.7,
            "lp/mm": 1.2,
        },
        "roi 3": {
            "distance from center": -1.8,
            "angle": 0,
            "roi radius": 0.7,
            "lp/mm": 1.8,
        },
        "roi 4": {
            "distance from center": 1.8,
            "angle": 90,
            "roi radius": 0.7,
            "lp/mm": 2.4,
        },
    }
    low_contrast_roi_settings = {
        "roi 1": {"distance from center": 2.6, "angle": -45, "roi radius": 0.6},
        "roi 2": {"distance from center": 2.6, "angle": -135, "roi radius": 0.6},
        "roi 3": {"distance from center": 2.6, "angle": 45, "roi radius": 0.6},
        "roi 4": {"distance from center": 2.6, "angle": 135, "roi radius": 0.6},
    }
    low_contrast_background_roi_settings = {
        "roi 1": {"distance from center": 0.5, "angle": 90, "roi radius": 0.25},
        "roi 2": {"distance from center": 0.5, "angle": -90, "roi radius": 0.25},
    }

    @staticmethod
    def run_demo() -> None:
        """Run the Sun Nuclear kV-QA phantom analysis demonstration."""
        snc = SNCkV.from_demo_image()
        snc.analyze()
        snc.plot_analyzed_image()

    def _phantom_radius_calc(self) -> float:
        """The radius of the phantom in pixels; the value itself doesn't matter, it's just
        used for relative distances to ROIs.

        Returns
        -------
        radius : float
        """
        return math.sqrt(self.phantom_ski_region.bbox_area) * 0.1071

    def _phantom_angle_calc(self) -> float:
        """The angle of the phantom. This assumes the user is using the stand that comes with the phantom,
        which angles the phantom at 45 degrees.

        Returns
        -------
        angle : float
            The angle in degrees.
        """
        angle = np.degrees(self.phantom_ski_region.orientation) + 180
        if np.isclose(angle, 135, atol=5):
            return angle
        else:
            raise ValueError(
                "The phantom angle was not near 135 degrees per manufacturer recommendations. Please adjust the phantom."
            )


class SNCMV(SNCkV):
    _demo_filename = "SNC-MV.dcm"
    common_name = "SNC MV-QA"
    phantom_bbox_size_mm2 = 118**2
    phantom_outline_object = {"Rectangle": {"width ratio": 7.5, "height ratio": 7.5}}
    high_contrast_roi_settings = {
        "roi 1": {
            "distance from center": -2.3,
            "angle": 0,
            "roi radius": 0.8,
            "lp/mm": 0.1,
        },
        "roi 2": {
            "distance from center": 2.3,
            "angle": 90,
            "roi radius": 0.8,
            "lp/mm": 0.2,
        },
        "roi 3": {
            "distance from center": 2.3,
            "angle": 0,
            "roi radius": 0.8,
            "lp/mm": 0.5,
        },
        "roi 4": {
            "distance from center": -2.3,
            "angle": 90,
            "roi radius": 0.8,
            "lp/mm": 1.0,
        },
    }
    low_contrast_roi_settings = {
        "roi 1": {"distance from center": 3.4, "angle": -45, "roi radius": 0.7},
        "roi 2": {"distance from center": 3.4, "angle": 45, "roi radius": 0.7},
        "roi 3": {"distance from center": 3.4, "angle": 135, "roi radius": 0.7},
        "roi 4": {"distance from center": 3.4, "angle": -135, "roi radius": 0.7},
    }
    low_contrast_background_roi_settings = {
        "roi 1": {"distance from center": 0.7, "angle": 0, "roi radius": 0.2},
        "roi 2": {"distance from center": -0.7, "angle": 0, "roi radius": 0.2},
    }

    @staticmethod
    def run_demo() -> None:
        """Run the Sun Nuclear MV-QA phantom analysis demonstration."""
        snc = SNCMV.from_demo_image()
        snc.analyze()
        snc.plot_analyzed_image()

    def _phantom_angle_calc(self) -> float:
        """The angle of the phantom. This assumes the user is using the stand that comes with the phantom,
        which angles the phantom at 45 degrees.

        Returns
        -------
        angle : float
            The angle in degrees.
        """
        return 45

    def _phantom_radius_calc(self) -> float:
        """The radius of the phantom in pixels; the value itself doesn't matter, it's just
        used for relative distances to ROIs.

        Returns
        -------
        radius : float
        """
        return math.sqrt(self.phantom_ski_region.bbox_area) * 0.095


class SNCMV12510(SNCMV):
    """The older SNC MV QA phantom w/ model number 1251000"""

    _demo_filename = "SNC_MV_12510.dcm"
    common_name = "SNC MV-QA (12510)"
    phantom_bbox_size_mm2 = 136**2
    phantom_outline_object = {"Rectangle": {"width ratio": 7.3, "height ratio": 6.2}}
    high_contrast_roi_settings = {
        "roi 1": {
            "distance from center": -1.7,
            "angle": 0,
            "roi radius": 0.7,
            "lp/mm": 0.1,
        },
        "roi 2": {
            "distance from center": 2.0,
            "angle": 80,
            "roi radius": 0.7,
            "lp/mm": 0.2,
        },
        "roi 3": {
            "distance from center": 2.4,
            "angle": 0,
            "roi radius": 0.7,
            "lp/mm": 0.5,
        },
        "roi 4": {
            "distance from center": -2.0,
            "angle": 100,
            "roi radius": 0.7,
            "lp/mm": 1.0,
        },
    }
    low_contrast_roi_settings = {
        "roi 1": {"distance from center": 3.1, "angle": -40, "roi radius": 0.7},
        "roi 2": {"distance from center": 3.1, "angle": 40, "roi radius": 0.7},
        "roi 3": {"distance from center": 2.5, "angle": 130, "roi radius": 0.7},
        "roi 4": {"distance from center": 2.5, "angle": -130, "roi radius": 0.7},
    }
    low_contrast_background_roi_settings = {
        "roi 1": {"distance from center": 1.0, "angle": 0, "roi radius": 0.2},
        "roi 2": {"distance from center": -0.2, "angle": 0, "roi radius": 0.2},
    }

    def _phantom_radius_calc(self) -> float:
        """The radius of the phantom in pixels; the value itself doesn't matter, it's just
        used for relative distances to ROIs.

        Returns
        -------
        radius : float
        """
        return math.sqrt(self.phantom_ski_region.bbox_area) * 0.105


class LeedsTOR(ImagePhantomBase):
    _demo_filename = "leeds.dcm"
    common_name = "Leeds"
    phantom_bbox_size_mm2 = 148**2
    _is_ccw = False
    phantom_outline_object = {"Circle": {"radius ratio": 0.97}}
    high_contrast_roi_settings = {
        "roi 1": {
            "distance from center": 0.2895,
            "angle": 54.62,
            "roi radius": 0.04,
            "lp/mm": 0.5,
        },
        "roi 2": {
            "distance from center": 0.187,
            "angle": 25.1,
            "roi radius": 0.04,
            "lp/mm": 0.56,
        },
        "roi 3": {
            "distance from center": 0.1848,
            "angle": 335.5,
            "roi radius": 0.04,
            "lp/mm": 0.63,
        },
        "roi 4": {
            "distance from center": 0.238,
            "angle": 80.06,
            "roi radius": 0.03,
            "lp/mm": 0.71,
        },
        "roi 5": {
            "distance from center": 0.0916,
            "angle": 62.96,
            "roi radius": 0.03,
            "lp/mm": 0.8,
        },
        "roi 6": {
            "distance from center": 0.093,
            "angle": -64,
            "roi radius": 0.02,
            "lp/mm": 0.9,
        },
        "roi 7": {
            "distance from center": 0.239,
            "angle": 101.98,
            "roi radius": 0.015,
            "lp/mm": 1.0,
        },
        "roi 8": {
            "distance from center": 0.0907,
            "angle": 122.62,
            "roi radius": 0.015,
            "lp/mm": 1.12,
        },
        "roi 9": {
            "distance from center": 0.09515,
            "angle": 239.07,
            "roi radius": 0.015,
            "lp/mm": 1.25,
        },
        "roi 10": {
            "distance from center": 0.2596,
            "angle": 115.8,
            "roi radius": 0.012,
            "lp/mm": 1.4,
        },
        "roi 11": {
            "distance from center": 0.138,
            "angle": 145,
            "roi radius": 0.012,
            "lp/mm": 1.6,
        },
        "roi 12": {
            "distance from center": 0.13967,
            "angle": 216.4,
            "roi radius": 0.010,
            "lp/mm": 1.8,
        },
    }
    low_contrast_background_roi_settings = {
        "roi 1": {"distance from center": 0.65, "angle": 30, "roi radius": 0.025},
        "roi 2": {"distance from center": 0.65, "angle": 120, "roi radius": 0.025},
        "roi 3": {"distance from center": 0.65, "angle": 210, "roi radius": 0.025},
        "roi 4": {"distance from center": 0.65, "angle": 300, "roi radius": 0.025},
    }
    low_contrast_roi_settings = {
        # set 1
        "roi 1": {"distance from center": 0.785, "angle": 30, "roi radius": 0.025},
        "roi 2": {"distance from center": 0.785, "angle": 45, "roi radius": 0.025},
        "roi 3": {"distance from center": 0.785, "angle": 60, "roi radius": 0.025},
        "roi 4": {"distance from center": 0.785, "angle": 75, "roi radius": 0.025},
        "roi 5": {"distance from center": 0.785, "angle": 90, "roi radius": 0.025},
        "roi 6": {"distance from center": 0.785, "angle": 105, "roi radius": 0.025},
        "roi 7": {"distance from center": 0.785, "angle": 120, "roi radius": 0.025},
        "roi 8": {"distance from center": 0.785, "angle": 135, "roi radius": 0.025},
        "roi 9": {"distance from center": 0.785, "angle": 150, "roi radius": 0.025},
        # set 2
        "roi 10": {"distance from center": 0.785, "angle": 210, "roi radius": 0.025},
        "roi 11": {"distance from center": 0.785, "angle": 225, "roi radius": 0.025},
        "roi 12": {"distance from center": 0.785, "angle": 240, "roi radius": 0.025},
        "roi 13": {"distance from center": 0.785, "angle": 255, "roi radius": 0.025},
        "roi 14": {"distance from center": 0.785, "angle": 270, "roi radius": 0.025},
        "roi 15": {"distance from center": 0.785, "angle": 285, "roi radius": 0.025},
        "roi 16": {"distance from center": 0.785, "angle": 300, "roi radius": 0.025},
        "roi 17": {"distance from center": 0.785, "angle": 315, "roi radius": 0.025},
        "roi 18": {"distance from center": 0.785, "angle": 330, "roi radius": 0.025},
    }

    @staticmethod
    def run_demo() -> None:
        """Run the Leeds TOR phantom analysis demonstration."""
        leeds = LeedsTOR.from_demo_image()
        leeds.analyze()
        leeds.plot_analyzed_image()

    @lru_cache()
    def _phantom_angle_calc(self) -> float:
        """Determine the angle of the phantom.

        This is done by searching for square-like boxes of the canny image. There are usually two: one lead and
        one copper. The box with the highest intensity (lead) is identified. The angle from the center of the lead
        square bounding box and the phantom center determines the phantom angle.

        Returns
        -------
        angle : float
            The angle in degrees
        """
        start_angle_deg = self._determine_start_angle_for_circle_profile()
        circle = self._circle_profile_for_phantom_angle(start_angle_deg, is_ccw=True)
        peak_idx, _ = circle.find_fwxm_peaks(threshold=0.6, max_number=1)

        shift_percent = peak_idx[0] / len(circle.values)
        shift_radians = shift_percent * 2 * np.pi
        shift_radians_corrected = 2 * np.pi - shift_radians

        angle = np.degrees(shift_radians_corrected) + start_angle_deg
        return angle

    def _phantom_radius_calc(self) -> float:
        """Determine the radius of the phantom.

        The radius is determined by finding the largest of the detected blobs of the canny image and taking
        its major axis length.

        Returns
        -------
        radius : float
            The radius of the phantom in pixels. The actual value is not important; it is used for scaling the
            distances to the low and high contrast ROIs.
        """
        return math.sqrt(self.phantom_ski_region.bbox_area) * 0.515

    def _determine_start_angle_for_circle_profile(self) -> float:
        """Determine an appropriate angle for starting the circular profile
        used to determine the phantom angle.

        In most cases we can just use 0 degs but for the case where the phantom
        is set up near 0 degs, the peak of the circular profile will be split
        between the left and right sides of the profile.  We can check for this
        case by looking at a few of the peak indexes and determining whether
        they are all on the left or right side of the profile or split left and
        right.  If they're split left and right, then we we need to use a
        different circular profile start angle  to get an accurate angle
        determination

        Returns
        -------
        start_angle_deg: float
            The start angle to be used for the circular profile used to determine the
            phantom rotation.
        """

        circle = self._circle_profile_for_phantom_angle(0)
        peak_idxs, _ = circle.find_fwxm_peaks(threshold=0.6, max_number=4)
        on_left_half = [x < len(circle.values) / 2 for x in peak_idxs]
        aligned_to_zero_deg = not (all(on_left_half) or not any(on_left_half))
        return 90 if aligned_to_zero_deg else 0

    def _preprocess(self) -> None:
        self._check_if_counter_clockwise()

    def _sample_high_contrast_rois(self) -> list[HighContrastDiskROI]:
        """Sample the high-contrast line pair regions. We overload to find
        the center of the high-res block which can be offset relative
        to the center depending on the model"""
        # find the high-res block ROI
        regions = self._get_canny_regions()
        high_res_block_size = self.phantom_bbox_size_px * 0.23
        sorted_regions = (
            Enumerable(regions)
            .where(
                lambda r: math.isclose(r.bbox_area, high_res_block_size, rel_tol=0.75)
            )
            .where(
                lambda r: bbox_center(r).distance_to(self.phantom_center)
                < 0.1 * self.phantom_radius
            )
            .order_by_descending(
                lambda r: bbox_center(r).distance_to(self.phantom_center)
            )
            .to_list()
        )
        if not sorted_regions:
            raise ValueError(
                "Could not find high-resolution block within the leeds phantom. Try rotating the image."
            )
        high_res_center = bbox_center(sorted_regions[0])
        self.high_res_center = high_res_center

        # do the same as the base method but centered on the high-res block
        hc_rois = []
        for stng in self.high_contrast_roi_settings.values():
            roi = HighContrastDiskROI(
                self.image,
                self.phantom_angle + stng["angle"],
                self.phantom_radius * stng["roi radius"],
                self.phantom_radius * stng["distance from center"],
                high_res_center,
                self._high_contrast_threshold,
            )
            hc_rois.append(roi)
        return hc_rois

    def _check_if_counter_clockwise(self) -> None:
        """Determine if the low-contrast bubbles go from high to low clockwise or counter-clockwise."""
        circle = self._circle_profile_for_phantom_angle(0)
        peak_idx, _ = circle.find_fwxm_peaks(threshold=0.6, max_number=1)
        circle.values = np.roll(circle.values, -peak_idx[0])
        _, first_set = circle.find_peaks(
            search_region=(0.05, 0.45), threshold=0, min_distance=0.025, max_number=9
        )
        _, second_set = circle.find_peaks(
            search_region=(0.55, 0.95), threshold=0, min_distance=0.025, max_number=9
        )
        self._is_ccw = max(first_set) > max(second_set)
        if not self._is_ccw:
            self.image.fliplr()
            del (
                self.phantom_ski_region
            )  # clear the property to calculate it again since we flipped it

    def _circle_profile_for_phantom_angle(
        self, start_angle_deg: float, is_ccw: bool = False
    ) -> CollapsedCircleProfile:
        """Create a circular profile centered at phantom origin

        Parameters
        ----------
        start_angle_deg: float

            Angle in degrees at which to start the profile
        Returns
        -------
        circle : CollapsedCircleProfile
            The circular profile centered on the phantom center and origin set to the given start angle.
        """
        circle = CollapsedCircleProfile(
            self.phantom_center,
            self.phantom_radius * 0.79,
            self.image.array,
            width_ratio=0.04,
            ccw=is_ccw,
            start_angle=np.deg2rad(start_angle_deg),
        )
        circle.ground()
        circle.filter(size=0.01)
        circle.invert()
        return circle

    def _check_inversion(self):
        """We recycle the circle profile used for angle detection to determine the correct inversion
        The profile is mostly even except the bright lead area. If the lead area is darker than the mean, it's inverted.
        """
        circle = self._circle_profile_for_phantom_angle(start_angle_deg=0)
        p2, p50, p98 = np.percentile(circle.values, [2, 50, 98])
        if abs(p50 - p98) < abs(p50 - p2):
            self.image.invert()


class LeedsTORBlue(LeedsTOR):
    """The Leeds TOR (Blue) is for analyzing older Leeds phantoms which have slightly offset ROIs compared to the newer, red-ring variant."""

    common_name = "Leeds (Blue)"
    high_contrast_roi_settings = {
        "roi 1": {
            "distance from center": 0.3,
            "angle": 54.8,
            "roi radius": 0.04,
            "lp/mm": 0.5,
        },
        "roi 2": {
            "distance from center": 0.187,
            "angle": 25.1,
            "roi radius": 0.04,
            "lp/mm": 0.56,
        },
        "roi 3": {
            "distance from center": 0.187,
            "angle": -27.5,
            "roi radius": 0.04,
            "lp/mm": 0.63,
        },
        "roi 4": {
            "distance from center": 0.252,
            "angle": 79.7,
            "roi radius": 0.03,
            "lp/mm": 0.71,
        },
        "roi 5": {
            "distance from center": 0.092,
            "angle": 63.4,
            "roi radius": 0.03,
            "lp/mm": 0.8,
        },
        "roi 6": {
            "distance from center": 0.094,
            "angle": -65,
            "roi radius": 0.02,
            "lp/mm": 0.9,
        },
        "roi 7": {
            "distance from center": 0.252,
            "angle": -260,
            "roi radius": 0.02,
            "lp/mm": 1.0,
        },
        "roi 8": {
            "distance from center": 0.094,
            "angle": -240,
            "roi radius": 0.018,
            "lp/mm": 1.12,
        },
        "roi 9": {
            "distance from center": 0.0958,
            "angle": -120,
            "roi radius": 0.018,
            "lp/mm": 1.25,
        },
        "roi 10": {
            "distance from center": 0.27,
            "angle": 115,
            "roi radius": 0.015,
            "lp/mm": 1.4,
        },
        "roi 11": {
            "distance from center": 0.13,
            "angle": 150,
            "roi radius": 0.011,
            "lp/mm": 1.6,
        },
        "roi 12": {
            "distance from center": 0.135,
            "angle": -150,
            "roi radius": 0.011,
            "lp/mm": 1.8,
        },
    }
    low_contrast_background_roi_settings = {
        "roi 1": {"distance from center": 0.6, "angle": 30, "roi radius": 0.025},
        "roi 2": {"distance from center": 0.6, "angle": 120, "roi radius": 0.025},
        "roi 3": {"distance from center": 0.6, "angle": 210, "roi radius": 0.025},
        "roi 4": {"distance from center": 0.6, "angle": 300, "roi radius": 0.025},
    }
    low_contrast_roi_settings = {
        # set 1
        "roi 1": {"distance from center": 0.83, "angle": 30, "roi radius": 0.025},
        "roi 2": {"distance from center": 0.83, "angle": 45, "roi radius": 0.025},
        "roi 3": {"distance from center": 0.83, "angle": 60, "roi radius": 0.025},
        "roi 4": {"distance from center": 0.83, "angle": 75, "roi radius": 0.025},
        "roi 5": {"distance from center": 0.83, "angle": 90, "roi radius": 0.025},
        "roi 6": {"distance from center": 0.83, "angle": 105, "roi radius": 0.025},
        "roi 7": {"distance from center": 0.83, "angle": 120, "roi radius": 0.025},
        "roi 8": {"distance from center": 0.83, "angle": 135, "roi radius": 0.025},
        "roi 9": {"distance from center": 0.83, "angle": 150, "roi radius": 0.025},
        # set 2
        "roi 10": {"distance from center": 0.83, "angle": 210, "roi radius": 0.025},
        "roi 11": {"distance from center": 0.83, "angle": 225, "roi radius": 0.025},
        "roi 12": {"distance from center": 0.83, "angle": 240, "roi radius": 0.025},
        "roi 13": {"distance from center": 0.83, "angle": 255, "roi radius": 0.025},
        "roi 14": {"distance from center": 0.83, "angle": 270, "roi radius": 0.025},
        "roi 15": {"distance from center": 0.83, "angle": 285, "roi radius": 0.025},
        "roi 16": {"distance from center": 0.83, "angle": 300, "roi radius": 0.025},
        "roi 17": {"distance from center": 0.83, "angle": 315, "roi radius": 0.025},
        "roi 18": {"distance from center": 0.83, "angle": 330, "roi radius": 0.025},
    }

    @classmethod
    def from_demo_image(cls):
        raise NotImplementedError("There is no demo file for this analysis")


class DoselabMC2kV(ImagePhantomBase):
    common_name = "Doselab MC2 kV"
    _demo_filename = "Doselab_kV.dcm"
    phantom_bbox_size_mm2 = 26300
    detection_conditions = [is_right_size]
    phantom_outline_object = {"Rectangle": {"width ratio": 0.55, "height ratio": 0.63}}
    low_contrast_background_roi_settings = {
        "roi 1": {"distance from center": 0.27, "angle": 48.5, "roi radius": 0.025},
    }
    low_contrast_roi_settings = {
        "roi 1": {"distance from center": 0.27, "angle": -48.5, "roi radius": 0.025},
        "roi 2": {"distance from center": 0.225, "angle": -65, "roi radius": 0.025},
        "roi 3": {"distance from center": 0.205, "angle": -88.5, "roi radius": 0.025},
        "roi 4": {"distance from center": 0.22, "angle": -110, "roi radius": 0.025},
        "roi 5": {"distance from center": 0.22, "angle": 110, "roi radius": 0.025},
        "roi 6": {"distance from center": 0.205, "angle": 88.5, "roi radius": 0.025},
        "roi 7": {"distance from center": 0.225, "angle": 65, "roi radius": 0.025},
    }
    high_contrast_roi_settings = {
        "roi 1": {
            "distance from center": 0.17,
            "angle": -20,
            "roi radius": 0.013,
            "lp/mm": 0.6,
        },
        "roi 2": {
            "distance from center": 0.16,
            "angle": -2,
            "roi radius": 0.007,
            "lp/mm": 1.2,
        },
        "roi 3": {
            "distance from center": 0.164,
            "angle": 12.8,
            "roi radius": 0.005,
            "lp/mm": 1.8,
        },
        "roi 4": {
            "distance from center": 0.175,
            "angle": 24.7,
            "roi radius": 0.0035,
            "lp/mm": 2.4,
        },
    }

    @staticmethod
    def run_demo() -> None:
        """Run the Doselab MC2 kV-area phantom analysis demonstration."""
        leeds = DoselabMC2kV.from_demo_image()
        leeds.analyze()
        leeds.plot_analyzed_image()

    def _phantom_radius_calc(self) -> float:
        return math.sqrt(self.phantom_ski_region.bbox_area) * 1.214

    def _phantom_angle_calc(self) -> float:
        roi = self.phantom_ski_region
        angle = np.degrees(roi.orientation) + 90
        if not np.isclose(angle, 45, atol=5):
            raise ValueError(
                "Angles not close enough to the ideal 45 degrees. Check phantom setup or override angle."
            )
        return angle


class DoselabMC2MV(DoselabMC2kV):
    common_name = "Doselab MC2 MV"
    _demo_filename = "Doselab_MV.dcm"
    low_contrast_background_roi_settings = {
        "roi 1": {"distance from center": 0.27, "angle": 48.5, "roi radius": 0.025},
    }
    low_contrast_roi_settings = {
        "roi 1": {"distance from center": 0.27, "angle": -48.5, "roi radius": 0.025},
        "roi 2": {"distance from center": 0.225, "angle": -65, "roi radius": 0.025},
        "roi 3": {"distance from center": 0.205, "angle": -88.5, "roi radius": 0.025},
        "roi 4": {"distance from center": 0.22, "angle": -110, "roi radius": 0.025},
        "roi 5": {"distance from center": 0.22, "angle": 110, "roi radius": 0.025},
        "roi 6": {"distance from center": 0.205, "angle": 88.5, "roi radius": 0.025},
        "roi 7": {"distance from center": 0.225, "angle": 65, "roi radius": 0.025},
    }
    high_contrast_roi_settings = {
        "roi 1": {
            "distance from center": 0.23,
            "angle": -135.3,
            "roi radius": 0.012,
            "lp/mm": 0.1,
        },
        "roi 2": {
            "distance from center": 0.173,
            "angle": 161,
            "roi radius": 0.012,
            "lp/mm": 0.2,
        },
        "roi 3": {
            "distance from center": 0.237,
            "angle": 133,
            "roi radius": 0.012,
            "lp/mm": 0.4,
        },
        "roi 4": {
            "distance from center": 0.298,
            "angle": 122.9,
            "roi radius": 0.01,
            "lp/mm": 0.8,
        },
    }

    @staticmethod
    def run_demo() -> None:
        """Run the Doselab MC2 MV-area phantom analysis demonstration."""
        leeds = DoselabMC2MV.from_demo_image()
        leeds.analyze()
        leeds.plot_analyzed_image()


def take_centermost_roi(rprops: list[RegionProperties], image_shape: tuple[int, int]):
    """Return the ROI that is closest to the center."""
    larger_rois = [
        rprop for rprop in rprops if rprop.area > 20 and rprop.eccentricity < 0.9
    ]  # drop stray pixel ROIs and line-like ROIs
    center_roi = sorted(
        larger_rois,
        key=lambda p: abs(p.centroid[0] - image_shape[0] / 2)
        + abs(p.centroid[1] - image_shape[1] / 2),
    )[0]
    return center_roi
