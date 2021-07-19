"""The picket fence module is meant for analyzing EPID images where a "picket fence" MLC pattern has been made.
Physicists regularly check MLC positioning through this test. This test can be done using film and one can
"eyeball" it, but this is the 21st century and we have numerous ways of quantifying such data. This module
attains to be one of them. It can load in an EPID dicom image (or superimpose multiple images) and determine the MLC peaks, error of each MLC
pair to the picket, and give a few visual indicators for passing/warning/failing.

Features:

* **Analyze any MLC type** - Both default MLCs and custom MLCs can be used.
* **Easy-to-read pass/warn/fail overlay** - Analysis gives you easy-to-read tools for determining the status of an MLC pair.
* **Any Source-to-Image distance** - Whatever your clinic uses as the SID for picket fence, pylinac can account for it.
* **Account for panel translation** - Have an off-CAX setup? No problem. Translate your EPID and pylinac knows.
* **Account for panel sag** - If your EPID sags at certain angles, just tell pylinac and the results will be shifted.
"""
import copy
import dataclasses
import enum
import io
import os.path as osp
import warnings
from dataclasses import dataclass
from itertools import cycle
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Union, Tuple, List, Optional, BinaryIO, Iterable, Sequence

import matplotlib.pyplot as plt
from cached_property import cached_property
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from py_linq import Enumerable

from .core import image, pdf
from .core.geometry import Line, Rectangle, Point
from .core.io import get_url, retrieve_demo_file
from .core.profile import MultiProfile, SingleProfile, Interpolation
from .core.utilities import open_path, ResultBase, convert_to_enum
from .log_analyzer import load_log
from .settings import get_dicom_cmap


LEFT_MLC_PREFIX = "A"
RIGHT_MLC_PREFIX = "B"


class Orientation(enum.Enum):
    """Possible orientations of the image"""
    UP_DOWN = 'Up-Down'  #:
    LEFT_RIGHT = 'Left-Right'  #:


class MLCArrangement:
    """Construct an MLC array"""

    def __init__(self, leaf_arrangement: List[Tuple[int, float]], offset: float = 0):
        """

        Parameters
        ----------
        leaf_arrangement: Description of the leaf arrangement. List of tuples containing the number of leaves and leaf width.
            E.g. (10, 5) is 10 leaves with 5mm widths.
        offset: The offset in mm of the leaves. Used for asymmetric arrangements. E.g. -2.5mm will shift the arrangement 2.5mm to the left.
        """
        self.centers = []
        self.widths = []
        rolling_edge = 0
        for leaf_num, width in leaf_arrangement:
            self.centers += np.arange(start=rolling_edge + width / 2, stop=leaf_num * width + rolling_edge + width / 2,
                                      step=width).tolist()
            rolling_edge = self.centers[-1] + width / 2
            self.widths += [width] * leaf_num
        self.centers = [c - np.mean(self.centers) + offset for c in self.centers]


class MLC(enum.Enum):
    """The pre-built MLC types"""
    MILLENNIUM = {'name': 'Millennium', 'arrangement': MLCArrangement([(10, 10), (40, 5), (10, 10)])}  #:
    HD_MILLENNIUM = {'name': 'HD Millennium', 'arrangement': MLCArrangement([(10, 5), (40, 2.5), (10, 5)])}  #:
    HALCYON_DISTAL = {'name': 'Halcyon distal', 'arrangement': MLCArrangement([(60, 5)], offset=2.5)}  #:
    HALCYON_PROXIMAL = {'name': 'Halcyon proximal', 'arrangement': MLCArrangement([(60, 5)])}  #:


@dataclass
class PFResult(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""
    tolerance_mm: float  #:
    action_tolerance_mm: float  #:
    percent_leaves_passing: float  #:
    number_of_pickets: int  #:
    absolute_median_error_mm: float  #:
    max_error_mm: float  #:
    mean_picket_spacing_mm: float  #:
    offsets_from_cax_mm: List[float]  #:
    passed: bool  #:
    failed_leaves: Union[List[str], List[int]]  #:
    mlc_skew: float  #:


class PFDicomImage(image.LinacDicomImage):
    """A subclass of a DICOM image that checks for noise and inversion when instantiated. Can also adjust for EPID sag."""

    def __init__(self, path: str, **kwargs):
        crop_mm = kwargs.pop("crop_mm", 3)
        super().__init__(path, **kwargs)
        # crop the images so that Elekta images don't fail. See #168
        crop_pixels = int(round(crop_mm * self.dpmm))
        self.crop(pixels=crop_pixels)
        self._check_for_noise()

    def _check_for_noise(self) -> None:
        """Check if the image has extreme noise (dead pixel, etc) by comparing
        min/max to 1/99 percentiles and smoothing if need be."""
        safety_stop = 5
        while self._has_noise() and safety_stop > 0:
            self.filter(size=3)
            safety_stop -= 1

    def _has_noise(self) -> bool:
        """Helper method to determine if there is spurious signal in the image."""
        min = self.array.min()
        max = self.array.max()
        near_min, near_max = np.percentile(self.array, [0.5, 99.5])
        max_is_extreme = max > near_max * 1.25
        min_is_extreme = (min < near_min * 0.75) and (abs(min - near_min) > 0.1 * (near_max - near_min))
        return max_is_extreme or min_is_extreme

    def adjust_for_sag(self, sag: int, orientation: Union[str, Orientation]) -> None:
        """Roll the image to adjust for EPID sag."""
        orient = convert_to_enum(orientation, Orientation)
        direction = 'y' if orient == Orientation.UP_DOWN else 'x'
        self.roll(direction, sag)


class PicketFence:
    """A class used for analyzing EPID images where radiation strips have been formed by the
    MLCs. The strips are assumed to be parallel to one another and normal to the image edge;
    i.e. a "left-right" or "up-down" orientation is assumed. Further work could follow up by accounting
    for any angle.
    """

    def __init__(self, filename: Union[str, Path, BinaryIO], filter: Optional[int] = None, log: Optional[str] = None,
                 use_filename: bool = False, mlc: Union[MLC, MLCArrangement, str] = MLC.MILLENNIUM, crop_mm: int = 3):
        """
        Parameters
        ----------
        filename
            Name of the file as a string or a file-like object.
        filter
            If None (default), no filtering will be done to the image.
            If an int, will perform median filtering over image of size ``filter``.
        log
            Path to a log file corresponding to the delivery. The expected fluence of the log file is
            used to construct the pickets. MLC peaks are then compared to an absolute reference instead of
            a fitted picket.
        use_filename
            If False (default), no action will be performed.
            If True, the filename will be searched for keywords that describe the gantry and/or collimator angle.
            For example, if set to True and the file name was "PF_gantry45.dcm" the gantry would be interpreted as being at 45 degrees.
        mlc
            The MLC model of the image. Must be an option from the enum :class:`~pylinac.picketfence.MLCs` or
            an :class:`~pylinac.picketfence.MLCArrangement`.
        crop_mm
            The number of mm to crop from all edges. Elekta is infamous for having columns of dead pixels on the side of their images.
            These need to be cleaned up first. For Varian images, this really shouldn't make a difference unless the pickets are
            very close to the edge. Generally speaking, they shouldn't be for the best accuracy.
        """
        leaf_analysis_width: float  #:
        mlc_meas: List  #:
        pickets: List  #:
        tolerance: float  #:
        action_tolerance: float  #:
        image: PFDicomImage  #:

        if filename is not None:
            self.image = PFDicomImage(filename, use_filenames=use_filename, crop_mm=crop_mm)
            if isinstance(filter, int):
                self.image.filter(size=filter)
            self.image.ground()
            self.image.normalize()
        if log is not None:
            self._load_log(log)
        else:
            self._log_fits = None
        self._is_analyzed = False
        self.mlc = self._get_mlc_arrangement(mlc)

    @staticmethod
    def _get_mlc_arrangement(value: Union[MLC, MLCArrangement, str]) -> MLCArrangement:
        if isinstance(value, MLC):
            return value.value['arrangement']
        if isinstance(value, MLCArrangement):
            return value
        if isinstance(value, str):
            return [member.value['arrangement'] for name, member in MLC.__members__.items() if member.value['name'] == value][0]

    @classmethod
    def from_url(cls, url: str, filter: int = None):
        """Instantiate from a URL."""
        filename = get_url(url, progress_bar=True)
        return cls(filename, filter=filter)

    @classmethod
    def from_demo_image(cls, filter: int = None):
        """Construct a PicketFence instance using the demo image."""
        demo_file = retrieve_demo_file(url='AS1200.dcm')
        return cls(demo_file, filter=filter)

    @classmethod
    def from_multiple_images(cls, path_list: Iterable[Union[str, Path]], dtype: np.dtype = np.uint16, **kwargs):
        """Load and superimpose multiple images and instantiate a Starshot object.

        Parameters
        ----------
        path_list : iterable
            An iterable of path locations to the files to be loaded/combined.
        kwargs
            Passed to :func:`~pylinac.core.image.load_multiples`.
        """
        obj = cls.from_demo_image()
        # save a combined image to a temporary dir, then load it back in as a PFDicomImage
        with TemporaryDirectory() as tmp:
            filename = osp.join(tmp, 'mydcm.dcm')
            image.load_multiples(path_list, dtype=dtype, **kwargs).save(filename)
            obj.image = PFDicomImage(filename)
        return obj

    @property
    def passed(self) -> bool:
        """Boolean specifying if all MLC positions were within tolerance."""
        return all(m.passed for m in self.mlc_meas)

    @property
    def percent_passing(self) -> float:
        """Return the percentage of MLC positions under tolerance."""
        num_meas = Enumerable(self.mlc_meas).select_many(lambda m: m.passed).count()
        num_pass = Enumerable(self.mlc_meas).select_many(lambda m: m.passed).count(lambda p: p == True)
        return float(100 * num_pass / num_meas)

    @property
    def max_error(self) -> float:
        """Return the maximum error found."""
        return float(np.max(self._flattened_errors()))

    @property
    def max_error_picket(self) -> int:
        """Return the picket number where the maximum error occurred."""
        return sorted([m for m in self.mlc_meas], key=lambda x: max(x.error))[0].picket_num

    @property
    def max_error_leaf(self) -> int:
        """Return the leaf that had the maximum error."""
        return sorted([m for m in self.mlc_meas], key=lambda x: max(x.error))[0].leaf_num

    def _flattened_errors(self) -> List[float]:
        return Enumerable(self.mlc_meas).select_many(lambda m: m.error).to_list()

    def failed_leaves(self) -> Union[List[int], List[str]]:
        """A list of the failed leaves. Either the leaf number or the bank+leaf number if using separate leaves."""
        if not self._is_analyzed:
            raise ValueError("It appears the PF image has not been analyzed yet. Use .analyze() first.")
        return Enumerable(self.mlc_meas)\
            .where(lambda m: not all(m.passed))\
            .select_many(lambda m: m.full_leaf_nums)\
            .distinct()\
            .to_list()

    @property
    def abs_median_error(self) -> float:
        """Return the median error found."""
        return float(np.median(np.abs(self._flattened_errors())))

    @property
    def num_pickets(self) -> int:
        """Return the number of pickets determined."""
        return len(self.pickets)

    @property
    def mean_picket_spacing(self) -> float:
        """The average distance between pickets in mm."""
        sorted_pickets = sorted(self.pickets, key=lambda x: x.dist2cax)
        return float(np.mean([abs(sorted_pickets[idx].dist2cax - sorted_pickets[idx + 1].dist2cax) for idx in
                        range(len(sorted_pickets) - 1)]))

    def _load_log(self, log: str) -> None:
        """Load a machine log that corresponds to the picket fence delivery.

        This log determines the location of the pickets. The MLC peaks are then compared to the expected log pickets,
        not a simple fit of the peaks."""
        # load the log fluence image
        mlog = load_log(log)
        fl = mlog.fluence.expected.calc_map(equal_aspect=True)
        fli = image.load(fl, dpi=254)  # 254 pix/in => 1 pix/0.1mm (default fluence calc)

        # equate them such that they're the same size & DPI
        fluence_img, img_array = image.equate_images(fli, self.image)
        self.image.array = img_array.array

        # get picket fits from the modified fluence image
        pf = PicketFence.from_demo_image()
        pf.image = fluence_img
        pf.analyze()
        self._log_fits = cycle([p.get_fit() for p in pf.pickets])

    @staticmethod
    def run_demo(tolerance: float = 0.5, action_tolerance: float = None) -> None:
        """Run the Picket Fence demo using the demo image. See analyze() for parameter info."""
        pf = PicketFence.from_demo_image()
        pf.analyze(tolerance, action_tolerance=action_tolerance)
        print(pf.results())
        pf.plot_analyzed_image(leaf_error_subplot=True)

    def analyze(self, tolerance: float = 0.5, action_tolerance: Optional[float] = None,
                num_pickets: Optional[int] = None, sag_adjustment: Union[float, int] = 0,
                orientation: Optional[Union[Orientation, str]] = None, invert: bool = False, leaf_analysis_width_ratio: float = 0.4,
                picket_spacing: Optional[float] = None, height_threshold: float = 0.5, edge_threshold: float = 1.5,
                peak_sort: str = 'peak_heights', required_prominence: float = 0.2, fwxm: int = 50,
                separate_leaves: bool = False, nominal_gap_mm: float = 1, dlg_mm: float = 1.5) -> None:
        """Analyze the picket fence image.

        Parameters
        ----------
        tolerance
            The tolerance of difference in mm between an MLC pair position and the
            picket fit line.
        action_tolerance
            If None (default), no action tolerance is set or compared to.
            If an int or float, the MLC pair measurement is also compared to this
            tolerance. Must be lower than tolerance. This value is usually meant
            to indicate that a physicist should take an "action" to reduce the error,
            but should not stop treatment.
        num_pickets
            The number of pickets in the image. A helper parameter to limit the total number of pickets,
            only needed if analysis is catching more pickets than there really are.
        sag_adjustment
            The amount of shift in mm to apply to the image to correct for EPID sag.
            For Up-Down picket images, positive moves the image down, negative up.
            For Left-Right picket images, positive moves the image left, negative right.
        orientation
            If None (default), the orientation is automatically determined. If for some reason the determined
            orientation is not correct, you can pass it directly using this parameter.
            If passed a string with 'u' (e.g. 'up-down', 'u-d', 'up') it will set the orientation of the pickets as
            going up-down. If passed a string with 'l' (e.g. 'left-right', 'lr', 'left') it will set it as going
            left-right.
        invert
            If False (default), the inversion of the image is automatically detected and used.
            If True, the image inversion is reversed from the automatic detection. This is useful when runtime errors
            are encountered.
        leaf_analysis_width_ratio
            The ratio of the leaf width to use as part of the evaluation. E.g. if the ratio is 0.5, the center half of
            the leaf will be used. This helps avoid tongue and groove influence.
        picket_spacing
            If None (default), the spacing between pickets is determined automatically.
            If given, it should be an int or float specifying the number of **PIXELS** apart the pickets are.
        height_threshold
            The threshold that the MLC peak needs to be above to be considered a picket (vs background).
            Lower if not all leaves are being caught. Note that for FFF beams this would very likely need to be lowered.
        edge_threshold
            The threshold of pixel value standard deviation within the analysis window of the MLC leaf to be considered a full leaf.
            This is how pylinac removes MLCs that are eclipsed by the jaw. This also is how to
            omit or catch leaves at the edge of the field. Raise to catch more edge leaves.
        peak_sort
            Either 'peak_heights' or 'prominences'. This is the method for determining the peaks. Usually not needed
            unless the wrong number of pickets have been detected.
            See the scipy.signal.find_peaks function for more information.
        required_prominence
            The required height of the picket (not individual MLCs) to be considered a peak.
            Pylinac takes a mean of the image axis perpendicular to the leaf motion to get an initial guess of the peak
            locations and also to determine picket spacing. Changing this can be useful for wide-gap tests where
            the shape of the beam horns can form two or more local maximums in the picket area. Increase if for wide-gap
            images that are catching too many pickets. Consider lowering for FFF beams if there are analysis issues.

            .. warning::

                We do not recommend performing FFF wide-gap PF tests. Make your FFF pickets narrow or measure with a flat beam instead.

        fwxm
            For each MLC kiss, the profile is a curve from low to high to low. The FWXM (0-100) is the height to use to measure
            to determine the center of the curve, which is the surrogate for MLC kiss position. I.e. for each MLC kiss,
            what height of the picket should you use to actually determine the center location? It is unusual to change this.
            If you have something in the way (we've seen crazy examples with a BB in the way) you may want to increase this.
        separate_leaves
            Whether to analyze leaves individually (each tip) or as a set (combined, center of the picket). False is
            the default for backwards compatibility.
        nominal_gap_mm
            The gap of the pickets in mm as set by the DICOM plan. Only used when separate leaves is True.
        dlg_mm
            The DLG value of the energy. Combined with the nominal gap, this calculates where the leaf tips should be.
            Only used when separate leaves is True.
        """
        if action_tolerance is not None and tolerance < action_tolerance:
            raise ValueError("Tolerance cannot be lower than the action tolerance")
        self.tolerance = tolerance
        self.action_tolerance = action_tolerance
        self.leaf_analysis_width = leaf_analysis_width_ratio
        self.separate_leaves = separate_leaves

        if invert:
            self.image.invert()

        self._orientation = orientation
        # adjust for sag
        if sag_adjustment != 0:
            sag_pixels = int(round(sag_adjustment * self.image.dpmm))
            self.image.adjust_for_sag(sag_pixels, self.orientation)

        if self.orientation == Orientation.UP_DOWN:
            leaf_prof = np.mean(self.image, 0)
        else:
            leaf_prof = np.mean(self.image, 1)
        leaf_prof = MultiProfile(leaf_prof)
        peak_idxs, peak_vals = leaf_prof.find_fwxm_peaks(min_distance=0.02, threshold=height_threshold, max_number=num_pickets,
                                                         peak_sort=peak_sort, required_prominence=required_prominence)
        # get picket spacing if not set by user
        if picket_spacing is None:
            picket_spacing = np.median(np.diff(np.sort(peak_idxs)))

        # loop through each leaf row and analyze each MLC kiss
        self.mlc_meas = []
        for leaf_num, center, width in self._leaves_in_view(leaf_analysis_width_ratio):
            for picket_num, (picket_idx, picket_peak_val) in enumerate(zip(peak_idxs, peak_vals)):
                window = self._get_mlc_window(leaf_center=center, leaf_width=width, approx_idx=picket_idx,
                                              spacing=picket_spacing)
                if self._is_mlc_peak_in_window(window, height_threshold, edge_threshold, picket_peak_val):
                    self.mlc_meas.append(MLCValue(picket_num=picket_num, approx_idx=picket_idx, leaf_width=width,
                                                  leaf_center=center, picket_spacing=picket_spacing,
                                                  orientation=self.orientation,
                                                  leaf_analysis_width_ratio=leaf_analysis_width_ratio,
                                                  tolerance=tolerance, action_tolerance=action_tolerance,
                                                  leaf_num=leaf_num, approx_peak_val=picket_peak_val,
                                                  image_window=window, image=self.image, fwxm=fwxm,
                                                  separate_leaves=separate_leaves, nominal_gap_mm=nominal_gap_mm,
                                                  dlg_mm=dlg_mm))
        if not self.mlc_meas:
            raise ValueError("No MLC measurements were found. This may be due to an incorrect inversion. Try setting invert=True. Or, you may have passed an incorrect orientation.")

        # drop any leaf rows that don't have the right amount of MLC kisses (i.e. near edge where one is dropped)
        median_num_leaves = Enumerable(self.mlc_meas) \
            .group_by(key=lambda m: m.leaf_num) \
            .median(lambda m: len(m))
        full_leaves = Enumerable(self.mlc_meas) \
            .group_by(key=lambda m: m.leaf_num) \
            .where(lambda m: len(m) == median_num_leaves) \
            .select_many(lambda m: m) \
            .select(lambda m: m.leaf_num) \
            .distinct() \
            .to_list()
        self.mlc_meas = [m for m in self.mlc_meas if m.leaf_num in full_leaves]
        if any([True for m in self.mlc_meas if m.leaf_num not in full_leaves]):
            warnings.warn("Some leaves were removed from analysis because they were not detected for all pickets. If some valid leaves are missing try adjusting height_threshold or edge_threshold")

        # retrospectively create the pickets and update the individual MLC measurements so error can be calculated
        self.pickets = []
        for picket_num, _ in enumerate(peak_idxs):
            self.pickets.append(Picket([m for m in self.mlc_meas if m.picket_num == picket_num],
                                       log_fits=self._log_fits, orientation=self.orientation,
                                       image=self.image, tolerance=tolerance, dlg=dlg_mm, nominal_gap=nominal_gap_mm,
                                       separate_leaves=separate_leaves))

        self._is_analyzed = True

    def _is_mlc_peak_in_window(self, window, height_threshold, edge_threshold, picket_peak_val) -> bool:
        """Whether the MLC peak is inside the given window. E.g. the jaw could be closed or at an edge."""
        if self.orientation == Orientation.UP_DOWN:
            std = np.std(window, axis=1)
        else:
            std = np.std(window, axis=0)
        is_above_height_threshold = np.max(window) > height_threshold * picket_peak_val
        is_not_at_edge = max(std) < edge_threshold * np.median(std)
        return is_above_height_threshold and is_not_at_edge

    def _get_mlc_window(self, leaf_center, leaf_width, approx_idx, spacing) -> np.ndarray:
        """A small 2D window of the image that contains the area around the picket."""
        leaf_width_px = leaf_width * self.image.dpmm
        leaf_center_px = leaf_center * self.image.dpmm + (
            self.image.shape[0] / 2 if self.orientation == Orientation.UP_DOWN else self.image.shape[1] / 2)
        if self.orientation == Orientation.UP_DOWN:
            # crop edges to image boundary if need be; if the pickets are too close to edge we could spill outside
            left_edge = max(int(approx_idx - spacing / 2), 0)
            right_edge = min(int(approx_idx + spacing / 2), self.image.shape[1])
            top_edge = max(int(leaf_center_px - leaf_width_px / 2), 0)
            bottom_edge = min(int(leaf_center_px + leaf_width_px / 2), self.image.shape[0])
            array = self.image[top_edge:bottom_edge, left_edge:right_edge]
        else:
            top_edge = max(int(approx_idx - spacing / 2), 0)
            bottom_edge = min(int(approx_idx + spacing / 2), self.image.shape[0])
            left_edge = max(int(leaf_center_px - leaf_width_px / 2), 0)
            right_edge = min(int(leaf_center_px + leaf_width_px / 2), self.image.shape[1])
            array = self.image[top_edge:bottom_edge, left_edge:right_edge]
        return array

    def _leaves_in_view(self, analysis_width) -> List[Tuple[int, int, int]]:
        """Crop the leaves if not all leaves are in view."""
        range = self.image.shape[0] / 2 if self.orientation == Orientation.UP_DOWN else self.image.shape[1] / 2
        # cut off the edge so that we're not halfway through a leaf.
        range -= max(self.mlc.widths[0] * analysis_width, self.mlc.widths[-1] * analysis_width) * self.image.dpmm
        leaves = [i for i, c in enumerate(self.mlc.centers) if abs(c) < (range / self.image.dpmm)]
        return [(leaf_num, center, width) for leaf_num, center, width in
                zip(leaves, self.mlc.centers[leaves[0]:leaves[-1] + 1], self.mlc.widths[leaves[0]:leaves[-1] + 1])]

    def plot_analyzed_image(self, guard_rails: bool = True, mlc_peaks: bool = True, overlay: bool = True,
                            leaf_error_subplot: bool = True, show: bool = True,
                            figure_size: Union[str, Tuple] = 'auto') -> None:
        """Plot the analyzed image.

        Parameters
        ----------
        guard_rails
            Do/don't plot the picket "guard rails" around the ideal picket
        mlc_peaks
            Do/don't plot the detected MLC peak positions.
        overlay
            Do/don't plot the alpha overlay of the leaf status.
        leaf_error_subplot
            If True, plots a linked leaf error subplot adjacent to the PF image plotting the average and standard
            deviation of leaf error.
        show
            Whether to display the plot. Set to false for saving to a figure, etc.
        figure_size
            Either 'auto' or a tuple. If auto, the figure size is set depending on the orientation. If a tuple, this is the
            figure size to use.
        """
        if not self._is_analyzed:
            raise RuntimeError("The image must be analyzed first. Use .analyze().")

        # plot the image
        if figure_size == 'auto':
            if self.orientation == Orientation.UP_DOWN:
                figure_size = (12, 8)
            else:
                figure_size = (9, 9)
        fig, ax = plt.subplots(figsize=figure_size)
        ax.imshow(self.image.array, cmap=get_dicom_cmap())

        if leaf_error_subplot:
            self._add_leaf_error_subplot(ax)

        if guard_rails:
            for picket in self.pickets:
                picket.add_guards_to_axes(ax.axes)
        if mlc_peaks:
            for mlc_meas in self.mlc_meas:
                mlc_meas.plot2axes(ax.axes, width=1.5)

        if overlay:
            for mlc_meas in self.mlc_meas:
                mlc_meas.plot_overlay2axes(ax.axes)

        # plot CAX
        ax.plot(self.image.center.x, self.image.center.y, 'r+', ms=12, markeredgewidth=3)

        # tighten up the plot view
        ax.set_xlim([0, self.image.shape[1]])
        ax.set_ylim([0, self.image.shape[0]])
        ax.axis('off')

        if show:
            plt.show()

    def _add_leaf_error_subplot(self, ax: plt.Axes) -> None:
        """Add a bar subplot showing the leaf error."""

        # make the new axis
        divider = make_axes_locatable(ax)
        if self.orientation == Orientation.UP_DOWN:
            axtop = divider.append_axes('right', 2, pad=1, sharey=ax)
        else:
            axtop = divider.append_axes('bottom', 2, pad=1, sharex=ax)

        # get leaf positions, errors, standard deviation, and leaf numbers
        if self.orientation == Orientation.UP_DOWN:
            pos = [position.marker_lines[0].center.y for position in self.pickets[0].mlc_meas]
        else:
            pos = [position.marker_lines[0].center.x for position in self.pickets[0].mlc_meas]

        # calculate the error and stdev values per MLC pair
        error_stdev = []
        error_vals = []
        for leaf_num in set([m.leaf_num for m in self.mlc_meas]):
            error_vals.append(np.mean([m.error for m in self.mlc_meas if m.leaf_num == leaf_num]))
            error_stdev.append(np.std([m.error for m in self.mlc_meas if m.leaf_num == leaf_num]))

        # plot the leaf errors as a bar plot
        if self.orientation == Orientation.UP_DOWN:
            axtop.barh(pos, error_vals, xerr=error_stdev, height=self.leaf_analysis_width * 2, alpha=0.4,
                       align='center')
            # plot the tolerance line(s)
            axtop.axvline(self.tolerance, color='r', linewidth=3)
            if self.action_tolerance is not None:
                axtop.axvline(self.action_tolerance, color='m', linewidth=3)
            # reset xlims to comfortably include the max error or tolerance value
            axtop.set_xlim([0, max([max(error_vals), self.tolerance]) + 0.1])
        else:
            axtop.bar(pos, error_vals, yerr=error_stdev, width=self.leaf_analysis_width * 2, alpha=0.4, align='center')
            # plot the tolerance line(s)
            axtop.axhline(self.tolerance, color='r', linewidth=3)
            if self.action_tolerance is not None:
                axtop.axhline(self.action_tolerance, color='m', linewidth=3)
            axtop.set_ylim([0, max([max(error_vals), self.tolerance]) + 0.1])

        axtop.grid(True)
        axtop.set_title("Average Error (mm)")

    def save_analyzed_image(self, filename: Union[str, io.BytesIO], guard_rails: bool = True, mlc_peaks: bool = True, overlay: bool = True,
                            leaf_error_subplot: bool = False, **kwargs) -> None:
        """Save the analyzed figure to a file. See :meth:`~pylinac.picketfence.PicketFence.plot_analyzed_image()` for
        further parameter info.
        """
        self.plot_analyzed_image(guard_rails, mlc_peaks, overlay, leaf_error_subplot=leaf_error_subplot, show=False)
        plt.savefig(filename, **kwargs)
        if isinstance(filename, str):
            print(f"Picket fence image saved to: {osp.abspath(filename)}")

    def results(self) -> str:
        """Return results of analysis. Use with print()."""
        pass_pct = self.percent_passing
        offsets = ' '.join('{:.1f}'.format(pk.dist2cax) for pk in self.pickets)
        string = f"Picket Fence Results: \n{pass_pct:2.1f}% " \
                 f"Passed\nMedian Error: {self.abs_median_error:2.3f}mm \n" \
                 f"Mean picket spacing: {self.mean_picket_spacing:2.1f}mm \n" \
                 f"Picket offsets from CAX (mm): {offsets}\n" \
                 f"Max Error: {self.max_error:2.3f}mm on Picket: {self.max_error_picket}, Leaf: {self.max_error_leaf}\n" \
                 f"Failing leaves: {self.failed_leaves()}"
        return string

    def results_data(self, as_dict=False) -> Union[PFResult, dict]:
        """Present the results data and metadata as a dataclass, dict, or tuple.
        The default return type is a dataclass."""
        data = PFResult(tolerance_mm=self.tolerance,
                        action_tolerance_mm=self.action_tolerance,
                        percent_leaves_passing=self.percent_passing,
                        number_of_pickets=self.num_pickets,
                        absolute_median_error_mm=self.abs_median_error,
                        max_error_mm=self.max_error,
                        mean_picket_spacing_mm=self.mean_picket_spacing,
                        offsets_from_cax_mm=[pk.dist2cax for pk in self.pickets],
                        passed=self.passed,
                        failed_leaves=self.failed_leaves(),
                        mlc_skew=self.mlc_skew())
        if as_dict:
            return dataclasses.asdict(data)
        return data

    def publish_pdf(self, filename: Union[str, io.BytesIO], notes: str = None, open_file: bool = False, metadata: dict = None) -> None:
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
        """
        plt.ioff()
        canvas = pdf.PylinacCanvas(filename, page_title="Picket Fence Analysis", metadata=metadata)
        data = io.BytesIO()
        self.save_analyzed_image(data, leaf_error_subplot=True)
        canvas.add_image(data, location=(3, 8), dimensions=(12, 12))
        text = [
            'Picket Fence results:',
            f'Magnification factor (SID/SAD): {self.image.metadata.RTImageSID / self.image.metadata.RadiationMachineSAD:2.2f}',
            f'Tolerance (mm): {self.tolerance}',
            f'Leaves passing (%): {self.percent_passing:2.1f}',
            f'Absolute median error (mm): {self.abs_median_error:2.3f}',
            f'Mean picket spacing (mm): {self.mean_picket_spacing:2.1f}',
            f'Maximum error (mm): {self.max_error:2.3f} on Picket {self.max_error_picket}, Leaf {self.max_error_leaf}',
            f"Failing leaves: {self.failed_leaves()}"
        ]
        text.append(f'Gantry Angle: {self.image.gantry_angle:2.2f}')
        text.append(f'Collimator Angle: {self.image.collimator_angle:2.2f}')
        canvas.add_text(text=text, location=(10, 25.5))
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 5.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 5))

        canvas.finish()

        if open_file:
            open_path(filename)

    def mlc_skew(self) -> float:
        """Apparent rotation in degrees of the MLC. This could be conflated with the EPID skew, so be careful when interpreting this value."""
        return Enumerable(self.pickets).avg(lambda p: np.rad2deg(np.arctan(p.fit.c[0])))

    def plot_histogram(self, bins: int = 10, show: bool = True) -> None:
        """Plot a histogram of the leaf errors"""
        if not self._is_analyzed:
            raise ValueError("It appears the PF image has not been analyzed yet. Use .analyze() first.")
        errors = Enumerable(self.mlc_meas).select_many(lambda m: m.error).to_list()
        fig, ax = plt.subplots()
        ax.set_title("Leaf error histogram")
        ax.set_ylabel("Counts")
        ax.set_xlabel("Error (mm)")
        ax.hist(errors, bins=bins)
        if show:
            plt.show()

    @cached_property
    def orientation(self) -> Orientation:
        """The orientation of the image, either Up-Down or Left-Right."""
        # if orientation was passed in, use it
        if self._orientation is not None:
            return convert_to_enum(self._orientation, Orientation)

        # replace any dead pixels with median value
        temp_image = self.image.array.copy()
        temp_image[temp_image < np.median(temp_image)] = np.median(temp_image)

        # find "range" of 80 to 90th percentiles
        row_sum = np.sum(temp_image, 0)
        col_sum = np.sum(temp_image, 1)
        row80, row90 = np.percentile(row_sum, [85, 95])
        col80, col90 = np.percentile(col_sum, [85, 95])
        row_range = row90 - row80
        col_range = col90 - col80

        # The true picket side will have a greater difference in
        # percentiles than will the non-picket size.
        if row_range < col_range:
            orientation = Orientation.LEFT_RIGHT
        else:
            orientation = Orientation.UP_DOWN
        return orientation


class MLCValue:

    def __init__(self, picket_num: int, approx_idx: int, leaf_width: float, leaf_center: float, picket_spacing: float, orientation: Orientation,
                 leaf_analysis_width_ratio: float, tolerance: float, action_tolerance: Optional[float], leaf_num: int,
                 approx_peak_val: float, image_window: np.ndarray, image: PFDicomImage, fwxm: int, separate_leaves: bool, nominal_gap_mm: float,
                 dlg_mm: float):
        """Representation of an MLC kiss or of each MLC about a kiss."""
        self._approximate_idx = approx_idx
        self.picket_num = picket_num
        self._approximate_peak_vale = approx_peak_val
        self.leaf_width_px = leaf_width * image.dpmm
        self._leaf_center = leaf_center
        self.leaf_center_px = leaf_center * image.dpmm + (
            image.shape[0] / 2 if orientation == Orientation.UP_DOWN else image.shape[1] / 2)
        self.leaf_num = leaf_num
        self._image_window = image_window
        self._image = image
        self._fwxm = fwxm
        self._analysis_ratio: float = leaf_analysis_width_ratio
        self._spacing: float = picket_spacing
        self._orientation = orientation
        self._tolerance: float = tolerance
        self._action_tolerance: float = action_tolerance
        self._separate_leaves = separate_leaves
        self._nominal_gap_mm = nominal_gap_mm
        self._dlg_mm = dlg_mm
        self.position = self.get_peak_positions()
        self._fit = None

    def __repr__(self) -> str:
        return f"Leaf: {self.leaf_num}, Picket: {self.picket_num}"

    @property
    def full_leaf_nums(self) -> Sequence[Union[str, int]]:
        """The fully-qualified leaf names. This will be the simple leaf number for traditional analysis or the bank+leaf num for separate leaves."""
        if not self._separate_leaves:
            return [self.leaf_num,]
        else:
            return [f"{LEFT_MLC_PREFIX}{self.leaf_num}", f"{RIGHT_MLC_PREFIX}{self.leaf_num}"]

    def plot2axes(self, axes: plt.Axes, width: Union[float, int] = 1) -> None:
        """Plot the measurement to the axes."""
        for idx, line in enumerate(self.marker_lines):
            line.plot2axes(axes, width, color=self.bg_color[idx])

    def get_peak_positions(self) -> Sequence[float]:
        if self._orientation == Orientation.UP_DOWN:
            pix_vals = np.median(self._image_window, axis=0)
        else:
            pix_vals = np.median(self._image_window, axis=1)
        interpolation_factor = 100
        prof = SingleProfile(pix_vals, interpolation=Interpolation.LINEAR, interpolation_factor=interpolation_factor)
        fwxm = prof.fwxm_data(self._fwxm)
        if self._separate_leaves:
            left = fwxm['left index (exact)']/interpolation_factor + max(self._approximate_idx - self._spacing / 2, 0)
            right = fwxm['right index (exact)'] / interpolation_factor + max(self._approximate_idx - self._spacing / 2, 0)
            return left, right
        else:
            return (fwxm['center index (exact)']/interpolation_factor + max(self._approximate_idx - self._spacing / 2, 0),)   # crop to left edge if need be

    @property
    def passed(self) -> Sequence[bool]:
        """Whether the MLC kiss or leaf was within tolerance."""
        return [abs(error) < self._tolerance for error in self.error]

    @property
    def passed_action(self) -> Optional[Sequence[bool]]:
        """Whether the MLC kiss or leaf was within the action tolerance."""
        return [abs(error) < self._action_tolerance for error in self.error] if self._action_tolerance is not None else [True, True]

    @property
    def bg_color(self) -> Sequence[str]:
        """The color of the measurement when the PF image is plotted, based on pass/fail status."""
        colors = []
        for idx, passed in enumerate(self.passed):
            if not passed:
                colors.append('r')
            elif self._action_tolerance is not None:
                colors.append('b' if self.passed_action[idx] else 'm')
            else:
                colors.append('b')
        return colors

    @property
    def error(self) -> Sequence[float]:
        """The error (difference) of the MLC measurement and the picket fit.
        If using individual leaf analysis, returns both errors otherwise return one."""
        errors = []
        for line, sign in zip(self.marker_lines, (-1, 1)):
            if self._orientation == Orientation.UP_DOWN:
                picket_pos = self._fit(line.center.y)
                mlc_pos = line.center.x
            else:
                picket_pos = self._fit(line.center.x)
                mlc_pos = line.center.y
            if self._separate_leaves:  # offset the picket position by the DLG and nominal gap
                mag_factor = self._image.sid / 1000
                picket_pos += sign * (self._dlg_mm + self._nominal_gap_mm) * mag_factor / 2 * self._image.dpmm
            errors.append((mlc_pos - picket_pos) / self._image.dpmm)
        return errors

    @property
    def max_abs_error(self) -> float:
        """The maximum absolute error"""
        return np.max(np.abs([self.error]))

    @property
    def marker_lines(self) -> List[Line]:
        """The line(s) representing the MLC measurement position. When using separated leaves
        there are two lines. Traditional analysis returns one."""
        upper_point = self.leaf_center_px - self.leaf_width_px / 2 * self._analysis_ratio
        lower_point = self.leaf_center_px + self.leaf_width_px / 2 * self._analysis_ratio

        lines = []
        for mlc_position in self.position:
            if self._orientation == Orientation.UP_DOWN:
                line = Line((mlc_position, upper_point), (mlc_position, lower_point))
            else:
                line = Line((upper_point, mlc_position), (lower_point, mlc_position))
            lines.append(line)
        return lines

    def plot_overlay2axes(self, axes) -> None:
        """Create a rectangle overlay with the width of the error. I.e. it stretches from the picket fit to the MLC position. Gives more visual size to the """
        # calculate height (based on leaf analysis ratio)
        upper_point = self.leaf_center_px - self.leaf_width_px / 2 * self._analysis_ratio
        lower_point = self.leaf_center_px + self.leaf_width_px / 2 * self._analysis_ratio
        height = abs(upper_point - lower_point) * 0.8

        for idx, line in enumerate(self.marker_lines):
            width = abs(self.error[idx]) * self._image.dpmm
            y = line.center.y
            x = self.position[idx] - (self.error[idx] * self._image.dpmm) / 2

            if self._orientation == Orientation.UP_DOWN:
                r = Rectangle(width, height, center=(x, y))
                # if any of the values are over tolerance, show another larger rectangle to draw the eye
                if not self.passed[idx] or not self.passed_action[idx]:
                    re = Rectangle(self._image_window.shape[1] * 0.2, height * 1.2, center=(x, y))
                    re.plot2axes(axes, edgecolor='none', fill=True, alpha=0.5, facecolor=self.bg_color[idx])
            else:
                r = Rectangle(height, width, center=(x, y))
                if not self.passed[idx] or not self.passed_action[idx]:
                    re = Rectangle(self._image_window.shape[1] * 0.2, height * 1.2, center=(x, y))
                    re.plot2axes(axes, edgecolor='none', fill=True, alpha=0.5, facecolor=self.bg_color[idx])
            r.plot2axes(axes, edgecolor='none', fill=True, alpha=1, facecolor=self.bg_color[idx])


class Picket:
    """Holds picket information in a Picket Fence test."""

    def __init__(self, mlc_measurements: List[MLCValue], log_fits, orientation: Orientation,
                 image: PFDicomImage, tolerance: float, separate_leaves: bool,
                 nominal_gap: float, dlg: float):
        self.mlc_meas: List[MLCValue] = mlc_measurements
        self.log_fits = log_fits
        self.tolerance = tolerance
        self.orientation = orientation
        self.image = image
        self._separate_leaves = separate_leaves
        self._nominal_gap = nominal_gap
        self._dlg = dlg
        self.fit = self.get_fit()
        for m in self.mlc_meas:
            m._fit = self.fit

    def get_fit(self) -> np.poly1d:
        """The fit of a polynomial to the MLC measurements."""
        if self.log_fits is not None:
            return next(self.log_fits)
        x = Enumerable(self.mlc_meas)\
            .select_many(lambda m: [l.point1.y for l in m.marker_lines])\
            .to_list()
        y = Enumerable(self.mlc_meas)\
            .select_many(lambda m: [l.point1.x for l in m.marker_lines])\
            .to_list()
        if self.orientation == Orientation.UP_DOWN:
            fit = np.polyfit(x, y, 1)
        else:
            fit = np.polyfit(y, x, 1)
        return np.poly1d(fit)

    @property
    def dist2cax(self) -> float:
        """The distance from the CAX to the picket, in mm."""
        # TODO: see about using line and built-in dist to point.
        center_fit = np.poly1d(self.fit)
        if self.orientation == Orientation.UP_DOWN:
            length = self.image.shape[0]
        else:
            length = self.image.shape[1]
        x_data = np.arange(length)
        y_data = center_fit(x_data)
        idx = int(round(len(x_data) / 2))
        if self.orientation == Orientation.UP_DOWN:
            axis = 'x'
            p1 = Point(y_data[idx], x_data[idx])
        else:
            axis = 'y'
            p1 = Point(x_data[idx], y_data[idx])
        return (getattr(self.image.center, axis) - getattr(p1, axis)) / self.image.dpmm

    @property
    def left_guard_separated(self) -> Iterable[np.poly1d]:
        """The line representing the left-sided guard rails.
        When not doing separate analysis, the left and right rails will overlap."""
        l_fit = np.copy(self.fit)
        l_fit[-1] += self.tolerance * self.image.dpmm
        if not self._separate_leaves:
            return [np.poly1d(l_fit), ]
        else:
            other_fit = copy.copy(l_fit)
            mag_factor = self.image.sid / 1000
            l_fit[-1] += (self._dlg + self._nominal_gap)*mag_factor/2 * self.image.dpmm
            other_fit[-1] -= (self._dlg + self._nominal_gap)*mag_factor/2 * self.image.dpmm
            return [np.poly1d(l_fit), np.poly1d(other_fit)]

    @property
    def right_guard_separated(self):
        """The line representing the right-sided guard rails."""
        r_fit = np.copy(self.fit)
        r_fit[-1] -= self.tolerance * self.image.dpmm
        if not self._separate_leaves:
            return [np.poly1d(r_fit), ]
        else:
            other_fit = copy.copy(r_fit)
            mag_factor = self.image.sid / 1000
            r_fit[-1] -= (self._dlg + self._nominal_gap)*mag_factor / 2 * self.image.dpmm
            other_fit[-1] += (self._dlg + self._nominal_gap)*mag_factor / 2 * self.image.dpmm
            return [np.poly1d(r_fit), np.poly1d(other_fit)]

    def add_guards_to_axes(self, axis: plt.Axes, color: str = 'g') -> None:
        """Plot guard rails to the axis."""
        if self.orientation == Orientation.UP_DOWN:
            length = self.image.shape[0]
        else:
            length = self.image.shape[1]
        x_data = np.arange(length)
        left_y_data = self.left_guard_separated
        right_y_data = self.right_guard_separated
        for left, right in zip(left_y_data, right_y_data):
            if self.orientation == Orientation.UP_DOWN:
                axis.plot(left(x_data), x_data, color=color)
                axis.plot(right(x_data), x_data, color=color)
            else:
                axis.plot(x_data, left(x_data), color=color)
                axis.plot(x_data, right(x_data), color=color)
