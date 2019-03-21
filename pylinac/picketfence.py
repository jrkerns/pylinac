"""The picket fence module is meant for analyzing EPID images where a "picket fence" MLC pattern has been made.
Physicists regularly check MLC positioning through this test. This test can be done using film and one can
"eyeball" it, but this is the 21st century and we have numerous ways of quantifying such data. This module
attains to be one of them. It can load in an EPID dicom image (or superimpose multiple images) and determine the MLC peaks, error of each MLC
pair to the picket, and give a few visual indicators for passing/warning/failing.

Features:

* **Analyze either HD or regular MLCs** - Just pass a flag and tell pylinac whether it's HD or not.
* **Easy-to-read pass/warn/fail overlay** - Analysis gives you easy-to-read tools for determining the status of an MLC pair.
* **Any Source-to-Image distance** - Whatever your clinic uses as the SID for picket fence, pylinac can account for it.
* **Account for panel translation** - Have an off-CAX setup? No problem. Translate your EPID and pylinac knows.
* **Account for panel sag** - If your EPID sags at certain angles, just tell pylinac and the results will be shifted.
"""
from collections import Sequence
from functools import lru_cache
import os.path as osp
import io
from itertools import cycle
from tempfile import TemporaryDirectory
from typing import Union, Tuple, List

import argue
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

from pylinac.core.typing import NumberLike
from pylinac.core.utilities import open_path
from .core import image
from .core.geometry import Line, Rectangle, Point
from .core.io import get_url, retrieve_demo_file
from .core import pdf
from .core.profile import MultiProfile, SingleProfile
from .log_analyzer import load_log
from .settings import get_dicom_cmap

# possible orientations of the pickets.
UP_DOWN = 'Up-Down'
LEFT_RIGHT = 'Left-Right'


class PFDicomImage(image.LinacDicomImage):
    """A subclass of a DICOM image that checks for noise and inversion when instantiated. Can also adjust for EPID sag."""

    def __init__(self, path: str, **kwargs):
        super().__init__(path, **kwargs)
        self._check_for_noise()
        self.check_inversion_by_histogram()

    def _check_for_noise(self):
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

    def adjust_for_sag(self, sag, orientation):
        """Roll the image to adjust for EPID sag."""
        direction = 'y' if orientation == UP_DOWN else 'x'
        self.roll(direction, sag)


class PicketFence:
    """A class used for analyzing EPID images where radiation strips have been formed by the
    MLCs. The strips are assumed to be parallel to one another and normal to the image edge;
    i.e. a "left-right" or "up-down" orientation is assumed. Further work could follow up by accounting
    for any angle.

    Attributes
    ----------
    pickets: :class:`~pylinac.picketfence.PicketHandler`
    image: :class:`~pylinac.core.image.DicomImage`

    Examples
    --------
    Run the demo::
        >>> PicketFence.run_demo()

    Load the demo image:
        >>> pf = PicketFence.from_demo_image()

    Load an image along with its machine log:
        >>> pf_w_log = PicketFence('my/pf.dcm', log='my/log.bin')

    Typical session:
        >>> img_path = r"C:/QA/June/PF.dcm"  # the EPID image
        >>> mypf = PicketFence(img_path)
        >>> mypf.analyze(tolerance=0.5, action_tolerance=0.3)
        >>> print(mypf.results())
        >>> mypf.plot_analyzed_image()
    """
    def __init__(self, filename: str, filter: int=None, log: str=None, use_filename: bool=False):
        """
        Parameters
        ----------
        filename : str, None
            Name of the file as a string. If None, image must be loaded later.
        filter : int, None
            If None (default), no filtering will be done to the image.
            If an int, will perform median filtering over image of size ``filter``.
        log : str
            Path to a log file corresponding to the delivery. The expected fluence of the log file is
            used to construct the pickets. MLC peaks are then compared to an absolute reference instead of
            a fitted picket.
        use_filename : bool
            If False (default), no action will be performed.
            If True, the filename will be searched for keywords that describe the gantry and/or collimator angle.
            For example, if set to True and the file name was "PF_gantry45.dcm" the gantry would be interpreted as being at 45 degrees.
        """
        if filename is not None:
            self.image = PFDicomImage(filename, use_filenames=use_filename)
            if isinstance(filter, int):
                self.image.filter(size=filter)
        if log is not None:
            self._load_log(log)
        else:
            self._log_fits = None
        self._is_analyzed = False

    @classmethod
    def from_url(cls, url: str, filter: int=None):
        """Instantiate from a URL."""
        filename = get_url(url, progress_bar=True)
        return cls(filename, filter=filter)

    @classmethod
    def from_demo_image(cls, filter: int=None):
        """Construct a PicketFence instance using the demo image."""
        demo_file = retrieve_demo_file(url='EPID-PF-LR.dcm')
        return cls(demo_file, filter=filter)

    @classmethod
    def from_multiple_images(cls, path_list: Sequence):
        """Load and superimpose multiple images and instantiate a Starshot object.

        Parameters
        ----------
        path_list : iterable
            An iterable of path locations to the files to be loaded/combined.
        """
        obj = cls.from_demo_image()
        # save a combined image to a temporary dir, then load it back in as a PFDicomImage
        with TemporaryDirectory() as tmp:
            filename = osp.join(tmp, 'mydcm.dcm')
            image.load_multiples(path_list, method='mean').save(filename)
            obj.image = PFDicomImage(filename)
        return obj

    @property
    def passed(self) -> bool:
        """Boolean specifying if all MLC positions were within tolerance."""
        return self.pickets.passed

    @property
    def percent_passing(self) -> float:
        """Return the percentage of MLC positions under tolerance."""
        num = 0
        num_pass = 0
        for picket in self.pickets:
            num += len(picket.error_array)
            num_pass += sum(picket.error_array < self.settings.tolerance)
        pct_pass = 100 * num_pass / num
        return pct_pass

    @property
    def max_error(self) -> float:
        """Return the maximum error found."""
        return max(picket.max_error for picket in self.pickets)

    @property
    def max_error_picket(self) -> int:
        """Return the picket number where the maximum error occurred."""
        return np.argmax([picket.max_error for picket in self.pickets])

    @property
    def max_error_leaf(self) -> int:
        """Return the leaf that had the maximum error."""
        picket = self.pickets[self.max_error_picket]
        return np.argmax(picket.error_array)

    @property
    @lru_cache()
    def abs_median_error(self) -> float:
        """Return the median error found."""
        return np.median(np.hstack([picket.error_array for picket in self.pickets]))

    @property
    def num_pickets(self) -> int:
        """Return the number of pickets determined."""
        return len(self.pickets)

    def _load_log(self, log: str):
        """Load a machine log that corresponds to the picket fence delivery.

        This log determines the location of the pickets. The MLC peaks are then compared to the expected log pickets,
        not a simple fit of the peaks."""
        # load the log fluence image
        mlog = load_log(log)
        fl = mlog.fluence.expected.calc_map(equal_aspect=True)
        fli = image.load(fl, dpi=254)  # 254 pix/in => 1 pix/0.1mm (default fluence calc)

        # equate them such that they're the same size & DPI
        fluence_img, self.image = image.equate_images(fli, self.image)

        # get picket fits from the modified fluence image
        pf = PicketFence.from_demo_image()
        pf.image = fluence_img
        pf.analyze()
        self._log_fits = cycle([p.fit for p in pf.pickets])

    @staticmethod
    def run_demo(tolerance: float=0.5, action_tolerance: float=None):
        """Run the Picket Fence demo using the demo image. See analyze() for parameter info."""
        pf = PicketFence.from_demo_image()
        pf.analyze(tolerance, action_tolerance=action_tolerance)
        print(pf.results())
        pf.plot_analyzed_image(leaf_error_subplot=True)

    def analyze(self, tolerance: float=0.5, action_tolerance: float=None, hdmlc: bool=False, num_pickets: int=None,
                sag_adjustment: Union[float, int]=0,
                orientation: str=None, invert: bool=False):
        """Analyze the picket fence image.

        Parameters
        ----------
        tolerance : int, float
            The tolerance of difference in mm between an MLC pair position and the
            picket fit line.
        action_tolerance : int, float, None
            If None (default), no action tolerance is set or compared to.
            If an int or float, the MLC pair measurement is also compared to this
            tolerance. Must be lower than tolerance. This value is usually meant
            to indicate that a physicist should take an "action" to reduce the error,
            but should not stop treatment.
        hdmlc : bool
            If False (default), a standard (5/10mm leaves) Millennium MLC model is assumed.
            If True, an HD (2.5/5mm leaves) Millennium is assumed.
        num_pickets : int, None

            .. versionadded:: 0.8

            The number of pickets in the image. A helper parameter to limit the total number of pickets,
            only needed if analysis is catching more pickets than there really are.
        sag_adjustment : float, int

            .. versionadded:: 0.8

            The amount of shift in mm to apply to the image to correct for EPID sag.
            For Up-Down picket images, positive moves the image down, negative up.
            For Left-Right picket images, positive moves the image left, negative right.
        orientation : None, str

            .. versionadded: 1.6

            If None (default), the orientation is automatically determined. If for some reason the determined
            orientation is not correct, you can pass it directly using this parameter.
            If passed a string with 'u' (e.g. 'up-down', 'u-d', 'up') it will set the orientation of the pickets as
            going up-down. If passed a string with 'l' (e.g. 'left-right', 'lr', 'left') it will set it as going
            left-right.
        invert : bool

            .. versionadded: 1.7

            If False (default), the inversion of the image is automatically detected and used.
            If True, the image inversion is reversed from the automatic detection. This is useful when runtime errors
            are encountered.
        """
        if action_tolerance is not None and tolerance < action_tolerance:
            raise ValueError("Tolerance cannot be lower than the action tolerance")

        # crop the images so that Elekta images don't fail. See #168
        if not self._is_analyzed:
            self.image.crop(pixels=2)

        if invert:
            self.image.invert()

        """Pre-analysis"""
        self._orientation = orientation
        self.settings = Settings(self.orientation, tolerance, action_tolerance, hdmlc, self.image, self._log_fits)
        # adjust for sag
        if sag_adjustment != 0:
            sag_pixels = int(round(sag_adjustment * self.settings.dpmm))
            self.image.adjust_for_sag(sag_pixels, self.orientation)

        """Analysis"""
        self.pickets = PicketManager(self.image, self.settings, num_pickets)
        self._is_analyzed = True

    def plot_analyzed_image(self, guard_rails: bool=True, mlc_peaks: bool=True, overlay: bool=True,
                            leaf_error_subplot: bool=True, show: bool=True):
        """Plot the analyzed image.

        Parameters
        ----------
        guard_rails : bool
            Do/don't plot the picket "guard rails" around the ideal picket
        mlc_peaks : bool
            Do/don't plot the MLC positions.
        overlay : bool
            Do/don't plot the alpha overlay of the leaf status.
        leaf_error_subplot : bool

            .. versionadded:: 1.0

            If True, plots a linked leaf error subplot adjacent to the PF image plotting the average and standard
            deviation of leaf error.
        """
        # plot the image
        fig, ax = plt.subplots(figsize=self.settings.figure_size)
        ax.imshow(self.image.array, cmap=get_dicom_cmap())

        # generate a leaf error subplot if desired
        if leaf_error_subplot:
            self._add_leaf_error_subplot(ax)

        # plot guard rails and mlc peaks as desired
        for p_num, picket in enumerate(self.pickets):
            if guard_rails:
                picket.add_guards_to_axes(ax.axes)
            if mlc_peaks:
                for idx, mlc_meas in enumerate(picket.mlc_meas):
                    mlc_meas.plot2axes(ax.axes, width=1.5)

        # plot the overlay if desired.
        if overlay:
            o = Overlay(self.image, self.settings, self.pickets)
            o.add_to_axes(ax)

        # plot CAX
        ax.plot(self.image.center.x, self.image.center.y, 'r+', ms=12, markeredgewidth=3)

        # tighten up the plot view
        ax.set_xlim([0, self.image.shape[1]])
        ax.set_ylim([0, self.image.shape[0]])
        ax.axis('off')

        if show:
            plt.show()

    def _add_leaf_error_subplot(self, ax: plt.Axes):
        """Add a bar subplot showing the leaf error."""
        tol_line_height = [self.settings.tolerance, self.settings.tolerance]
        tol_line_width = [0, max(self.image.shape)]

        # make the new axis
        divider = make_axes_locatable(ax)
        if self.settings.orientation == UP_DOWN:
            axtop = divider.append_axes('right', 2, pad=1, sharey=ax)
        else:
            axtop = divider.append_axes('bottom', 2, pad=1, sharex=ax)

        # get leaf positions, errors, standard deviation, and leaf numbers
        pos, vals, err, leaf_nums = self.pickets.error_hist()

        # plot the leaf errors as a bar plot
        if self.settings.orientation == UP_DOWN:
            axtop.barh(pos, vals, xerr=err, height=self.pickets[0].sample_width * 2, alpha=0.4, align='center')
            # plot the tolerance line(s)
            # TODO: replace .plot() calls with .axhline when mpld3 fixes funtionality
            axtop.plot(tol_line_height, tol_line_width, 'r-', linewidth=3)
            if self.settings.action_tolerance is not None:
                axtop.plot(tol_line_height, tol_line_width, 'y-', linewidth=3)

            # reset xlims to comfortably include the max error or tolerance value
            axtop.set_xlim([0, max(max(vals), self.settings.tolerance) + 0.1])
        else:
            axtop.bar(pos, vals, yerr=err, width=self.pickets[0].sample_width * 2, alpha=0.4, align='center')
            axtop.plot(tol_line_width, tol_line_height,
                       'r-', linewidth=3)
            if self.settings.action_tolerance is not None:
                axtop.plot(tol_line_width, tol_line_height, 'y-', linewidth=3)
            axtop.set_ylim([0, max(max(vals), self.settings.tolerance) + 0.1])

        # add formatting to axis
        axtop.grid(True)
        axtop.set_title("Average Error (mm)")

    def save_analyzed_image(self, filename: str, guard_rails: bool=True, mlc_peaks: bool=True, overlay: bool=True,
                            leaf_error_subplot: bool=False, **kwargs):
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
                 f"Mean picket spacing: {self.pickets.mean_spacing:2.1f}mm \n" \
                 f"Picket offsets from CAX (mm): {offsets}\n" \
                 f"Max Error: {self.max_error:2.3f}mm on Picket: {self.max_error_picket}, Leaf: {self.max_error_leaf}"
        return string

    def publish_pdf(self, filename: str, notes: str=None, open_file: bool=False, metadata: dict=None):
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
            f'Magnification factor (SID/SAD): {self.image.metadata.RTImageSID/self.image.metadata.RadiationMachineSAD:2.2f}',
            f'Tolerance (mm): {self.settings.tolerance}',
            f'Leaves passing (%): {self.percent_passing:2.1f}',
            f'Absolute median error (mm): {self.abs_median_error:2.3f}',
            f'Mean picket spacing (mm): {self.pickets.mean_spacing:2.1f}',
            f'Maximum error (mm): {self.max_error:2.3f} on Picket {self.max_error_picket}, Leaf {self.max_error_leaf}',
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

    @property
    @lru_cache(maxsize=1)
    def orientation(self) -> str:
        """The orientation of the image, either Up-Down or Left-Right."""
        # if orientation was passed in, use it
        if type(self._orientation) is str:
            if 'u' in self._orientation.lower():
                return UP_DOWN
            elif 'l' in self._orientation.lower():
                return LEFT_RIGHT

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
            orientation = LEFT_RIGHT
        else:
            orientation = UP_DOWN
        return orientation


class Overlay:
    """Class for handling the "overlay" feature of the plot."""

    def __init__(self, image, settings, pickets):
        self.image = image
        self.settings = settings
        self.pickets = pickets

    def add_to_axes(self, axes):
        """Add the overlay to the axes."""
        rect_width = self.pickets[0].sample_width*2
        for mlc_num, mlc in enumerate(sorted(self.pickets, key=lambda x: len(x.mlc_meas))[0].mlc_meas):
            # get pass/fail status of all measurements across pickets for that MLC
            if self.settings.action_tolerance is not None:
                if all(picket.mlc_passed_action(mlc_num) for picket in self.pickets):
                    color = 'b'
                elif all(picket.mlc_passed(mlc_num) for picket in self.pickets):
                    color = 'm'
                else:
                    color = 'r'
            elif all(picket.mlc_passed(mlc_num) for picket in self.pickets):
                color = 'b'
            else:
                color = 'r'

            # create a rectangle overlay
            if self.settings.orientation == UP_DOWN:
                r = Rectangle(self.image.shape[1], rect_width, center=(self.image.center.x, mlc.center.y))
            else:
                r = Rectangle(rect_width, self.image.shape[0], center=(mlc.center.x, self.image.center.y))
            r.plot2axes(axes.axes, edgecolor='none', fill=True, alpha=0.1, facecolor=color)


class Settings:
    """Simple class to hold various settings and info for PF analysis/plotting."""
    def __init__(self, orientation, tolerance, action_tolerance, hdmlc, image, log_fits):
        self.orientation = orientation
        self.tolerance = tolerance
        self.action_tolerance = action_tolerance
        self.hdmlc = hdmlc
        self.image = image
        self.dpmm = image.dpmm
        self.mmpd = 1/image.dpmm
        try:
            self.image_center = image.cax
        except AttributeError:
            self.image_center = image.center
        self.log_fits = log_fits

    @property
    def figure_size(self) -> Tuple[int, int]:
        """The size of the figure to draw; depends on the picket orientation."""
        if self.orientation == UP_DOWN:
            return (12, 8)
        else:
            return (9, 9)

    @property
    def small_leaf_width(self) -> int:
        """The width of a "small" leaf in pixels."""
        leaf_width_mm = 5
        leaf_width_pixels = leaf_width_mm * self.dpmm
        if self.hdmlc:
            leaf_width_pixels /= 2
        return leaf_width_pixels

    @property
    def large_leaf_width(self) -> int:
        """The width of a "large" leaf in pixels."""
        return self.small_leaf_width * 2

    @property
    def number_small_leaves(self) -> int:
        """The number of small leaves; depends on HDMLC status."""
        return 40 if not self.hdmlc else 32

    @property
    def number_large_leaves(self) -> int:
        """The number of large leaves; depends on HDMLC status."""
        return 20 if not self.hdmlc else 28

    @property
    @lru_cache()
    def leaf_centers(self) -> np.ndarray:
        """Return a set of leaf centers perpendicular to the leaf motion based on the position of the CAX."""
        # generate a set of leaf center points based on physical widths of large and small leaves
        first_shift = self.large_leaf_width * (self.number_large_leaves / 2 - 1) + self.large_leaf_width * 0.75
        second_shift = self.small_leaf_width * (self.number_small_leaves - 1) + self.large_leaf_width * 0.75

        large_leaf_section = np.arange(self.number_large_leaves / 2) * self.large_leaf_width
        small_leaf_section = (np.arange(self.number_small_leaves) * self.small_leaf_width) + first_shift
        large_leaf_section2 = (np.arange(
            self.number_large_leaves / 2) * self.large_leaf_width) + first_shift + second_shift
        leaf_centers = np.concatenate((large_leaf_section, small_leaf_section, large_leaf_section2))

        # now adjust them to align with the iso
        if self.orientation == UP_DOWN:
            leaf30_center = self.image_center.y - self.small_leaf_width / 2
            edge = self.image.shape[0]
        else:
            leaf30_center = self.image_center.x - self.small_leaf_width / 2
            edge = self.image.shape[1]
        adjustment = leaf30_center - leaf_centers[29]
        leaf_centers += adjustment

        # only include values that are reasonable as values might extend past image (e.g. with small SID)
        values_in_image = (leaf_centers > 0 + self.large_leaf_width / 2) & (
        leaf_centers < edge - self.large_leaf_width / 2)
        leaf_centers = leaf_centers[values_in_image]
        return np.round(leaf_centers).astype(int)


class PicketManager:
    """Finds and handles the pickets of the image."""
    def __init__(self, image, settings, num_pickets):
        self.pickets = []
        self.image = image
        self.settings = settings
        self.num_pickets = num_pickets
        self.find_pickets()

    def error_hist(self) -> Tuple[List, ...]:
        """Returns several lists of information about the MLC measurements. For use with plotting."""
        # for each MLC, get the average and standard deviation of the error across all the pickets
        error_means = []
        error_stds = []
        error_plot_positions = []
        mlc_leaves = []
        for mlc_num, mlc_meas in enumerate(sorted(self.pickets, key=lambda x: len(x.mlc_meas))[0].mlc_meas):
            errors = []
            for picket in self.pickets:
                errors.append(picket.mlc_meas[mlc_num].error)
            error_means.append(np.mean(errors))
            error_stds.append(np.std(errors))
            mlc_leaves.append(mlc_meas.leaf_pair)
            if self.settings.orientation == UP_DOWN:
                error_plot_positions.append(mlc_meas.center.y)
            else:
                error_plot_positions.append(mlc_meas.center.x)

        return error_plot_positions, error_means, error_stds, mlc_leaves

    def find_pickets(self):
        """Find the pickets of the image."""
        leaf_prof = self.image_mlc_inplane_mean_profile
        peak_idxs = leaf_prof.find_peaks(min_distance=0.02, threshold=0.5, max_number=self.num_pickets)
        peak_spacing = np.median(np.diff(np.sort(peak_idxs)))
        if np.isnan(peak_spacing):
            peak_spacing = 20

        for peak_idx in peak_idxs:
            self.pickets.append(Picket(self.image, self.settings, peak_idx, peak_spacing/2))

    @property
    def passed(self) -> bool:
        """Whether all the pickets passed tolerance."""
        return all(picket.passed for picket in self)

    def __getitem__(self, item):
        return self.pickets[item]

    def __len__(self):
        return len(self.pickets)

    @property
    def image_mlc_inplane_mean_profile(self) -> MultiProfile:
        """A profile of the image along the MLC travel direction."""
        if self.settings.orientation == UP_DOWN:
            leaf_prof = np.mean(self.image, 0)
        else:
            leaf_prof = np.mean(self.image, 1)
        return MultiProfile(leaf_prof)

    @property
    def mean_spacing(self) -> np.ndarray:
        """The average distance between pickets in mm."""
        sorted_pickets = sorted(self.pickets, key=lambda x: x.dist2cax)
        return np.mean([abs(sorted_pickets[idx].dist2cax - sorted_pickets[idx+1].dist2cax) for idx in range(len(sorted_pickets)-1)])


class Picket:
    """Holds picket information in a Picket Fence test.

    Attributes
    ----------
    mlc_meas : list
        Holds :class:`~pylinac.picketfence.MLCMeas` objects.
    """
    def __init__(self, image, settings, approximate_idx, spacing):
        self.mlc_meas = []
        self.image = image
        self.settings = settings
        self.approximate_idx = approximate_idx
        self.spacing = spacing
        self._get_mlc_positions()

    def _get_mlc_positions(self):
        """Calculate the positions of all the MLC pairs."""
        # for each MLC...
        for mlc_num, mlc_center in enumerate(self.settings.leaf_centers):
            # find the MLC peak
            mlc_position = self.find_mlc_peak(mlc_center)
            # add MLC measurement object
            if mlc_position is not None:
                self.add_mlc_meas(mlc_center, mlc_position)
        # now add the picket fit to the measurement so it can calculate error, etc.
        for idx, meas in enumerate(self.mlc_meas):
            meas.fit = self.fit

    def find_mlc_peak(self, mlc_center):
        """Determine the center of the picket."""
        mlc_rows = np.arange(mlc_center - self.sample_width, mlc_center + self.sample_width + 1)
        if self.settings.orientation == UP_DOWN:
            pix_vals = np.median(self.picket_array[mlc_rows, :], axis=0)
        else:
            pix_vals = np.median(self.picket_array[:, mlc_rows], axis=1)
        if max(pix_vals) > np.percentile(self.picket_array, 80):
            prof = SingleProfile(pix_vals)
            fw80mc = prof.fwxm_center(70, interpolate=True)
            return fw80mc + self.approximate_idx - self.spacing

    def add_mlc_meas(self, mlc_center, mlc_position):
        """Add an MLC measurement point."""
        upper_point = mlc_center - self.sample_width / 2
        lower_point = mlc_center + self.sample_width / 2

        if self.settings.orientation == UP_DOWN:
            meas = MLCMeas((mlc_position, upper_point), (mlc_position, lower_point), self.settings)
        else:
            meas = MLCMeas((upper_point, mlc_position), (lower_point, mlc_position), self.settings)
        self.mlc_meas.append(meas)

    @property
    def sample_width(self) -> float:
        """The width to sample the MLC leaf (~40% of the leaf width)."""
        return np.round(np.median(np.diff(self.settings.leaf_centers) * 2 / 5) / 2).astype(int)

    @property
    @lru_cache()
    def picket_array(self) -> np.ndarray:
        """A slice of the whole image that contains the area around the picket."""
        if self.settings.orientation == UP_DOWN:
            left_edge = int(self.approximate_idx - self.spacing)
            right_edge = int(self.approximate_idx + self.spacing)
            # see #167 & #174
            if left_edge < 0:
                self.spacing += left_edge
                left_edge = int(self.approximate_idx - self.spacing)
                right_edge = int(self.approximate_idx + self.spacing)
            array = self.image[:, left_edge:right_edge]
        else:
            top_edge = int(self.approximate_idx - self.spacing)
            bottom_edge = int(self.approximate_idx + self.spacing)
            # see #167 & #174
            if top_edge < 0:
                self.spacing += top_edge
                top_edge = int(self.approximate_idx - self.spacing)
                bottom_edge = int(self.approximate_idx + self.spacing)
            array = self.image[top_edge:bottom_edge, :]
        return array

    @property
    def abs_median_error(self) -> np.ndarray:
        """The absolute median error of the MLC measurements."""
        return np.median(np.abs(self.error_array))

    @property
    def max_error(self) -> float:
        """The max error of the MLC measurements."""
        return self.error_array.max()

    @property
    @lru_cache()
    def error_array(self) -> np.ndarray:
        """An array containing the error values of all the measurements."""
        return np.array([meas.error for meas in self.mlc_meas])

    @property
    def passed(self) -> bool:
        """Whether or not all the measurements passed."""
        return all(meas.passed for meas in self.mlc_meas)

    def mlc_passed(self, mlc) -> bool:
        """Return whether a specific MLC has passed tolerance."""
        return self.mlc_meas[mlc].passed

    def mlc_passed_action(self, mlc) -> bool:
        """Return whether a specific MLC has passed the action tolerance."""
        if self.settings.action_tolerance is not None:
            return self.mlc_meas[mlc].passed_action
        else:
            raise AttributeError("No action tolerance was specified")

    @property
    @lru_cache(maxsize=1)
    def fit(self):
        """The fit of a polynomial to the MLC measurements."""
        if self.settings.log_fits is not None:
            return next(self.settings.log_fits)
        x = np.array([mlc.point1.y for mlc in self.mlc_meas])
        y = np.array([mlc.point1.x for mlc in self.mlc_meas])
        if self.settings.orientation == UP_DOWN:
            fit = np.polyfit(x, y, 1)
        else:
            fit = np.polyfit(y, x, 1)
        return np.poly1d(fit)

    @property
    def dist2cax(self) -> float:
        """The distance from the CAX to the picket, in mm."""
        center_fit = np.poly1d(self.fit)
        if self.settings.orientation == UP_DOWN:
            length = self.image.shape[0]
        else:
            length = self.image.shape[1]
        x_data = np.arange(length)
        y_data = center_fit(x_data)
        idx = int(round(len(x_data) / 2))
        if self.settings.orientation == UP_DOWN:
            axis = 'x'
            p1 = Point(y_data[idx], x_data[idx])
        else:
            axis = 'y'
            p1 = Point(x_data[idx], y_data[idx])
        return (getattr(self.image.center, axis) - getattr(p1, axis)) * self.settings.mmpd

    @property
    def left_guard(self):
        """The line representing the left side guard rail."""
        l_fit = np.copy(self.fit)
        l_fit[-1] += self.settings.tolerance / self.settings.mmpd
        return np.poly1d(l_fit)

    @property
    def right_guard(self):
        """The line representing the right side guard rail."""
        r_fit = np.copy(self.fit)
        r_fit[-1] -= self.settings.tolerance / self.settings.mmpd
        return np.poly1d(r_fit)

    def add_guards_to_axes(self, axis: plt.Axes, color: str='g'):
        """Plot guard rails to the axis."""
        if self.settings.orientation == UP_DOWN:
            length = self.image.shape[0]
        else:
            length = self.image.shape[1]
        x_data = np.arange(length)
        left_y_data = self.left_guard(x_data)
        right_y_data = self.right_guard(x_data)
        if self.settings.orientation == UP_DOWN:
            axis.plot(left_y_data, x_data, color=color)
            axis.plot(right_y_data, x_data, color=color)
        else:
            axis.plot(x_data, left_y_data, color=color)
            axis.plot(x_data, right_y_data, color=color)


class MLCMeas(Line):
    """Represents an MLC measurement."""
    def __init__(self, point1, point2, settings):
        super().__init__(point1, point2)
        self.settings = settings
        self.fit = None

    def plot2axes(self, axes: plt.Axes, width: NumberLike=1):
        """Plot the measurement to the axes."""
        super().plot2axes(axes, width, color=self.bg_color)

    @property
    def bg_color(self) -> str:
        """The color of the measurement when the PF image is plotted, based on pass/fail status."""
        if not self.passed:
            return 'r'
        elif self.settings.action_tolerance is not None:
            if self.passed_action:
                return 'b'
            else:
                return 'm'
        else:
            return 'b'

    @property
    def error(self) -> float:
        """The error (difference) of the MLC measurement and the picket fit."""
        if self.settings.orientation == UP_DOWN:
            picket_pos = self.fit(self.center.y)
            mlc_pos = self.center.x
        else:
            picket_pos = self.fit(self.center.x)
            mlc_pos = self.center.y
        return abs(mlc_pos - picket_pos) * self.settings.mmpd

    @property
    def passed(self) -> bool:
        """Whether the MLC measurement was under tolerance."""
        return self.error < self.settings.tolerance

    @property
    def passed_action(self) -> bool:
        """Whether the MLC measurement was under the action level tolerance."""
        if self.settings.action_tolerance is not None:
            return self.error < self.settings.action_tolerance

    @property
    @lru_cache()
    def leaf_pair(self) -> Tuple[int, int]:
        """The leaf pair that formed the MLC measurement.

        Returns
        -------
        tuple : 2 elements which are the two leaf numbers
        """
        leaves = [0, 0]

        # get distance between MLC point and EPID center in *pixels*
        if self.settings.orientation == UP_DOWN:
            mlc_loc = self.center.y
            epid_center = self.settings.image_center.y
        else:
            mlc_loc = self.center.x
            epid_center = self.settings.image_center.x
        mlc_dist = mlc_loc - epid_center

        # determine leaf number based on if it's in/not in the "small leaf" region
        small_region_extent = self.settings.small_leaf_width * self.settings.number_small_leaves / 2

        # large leaf region
        if not small_region_extent > mlc_dist > -small_region_extent:
            if np.sign(mlc_dist) > 0:  # positive, meaning
                # offset MLC distance to the edge of the small leaf region
                mlc_dist -= small_region_extent
                # divide the MLC distance by the leaf width and convert to leaf number
                leaf = int(round((abs(mlc_dist) + self.settings.large_leaf_width / 2) / self.settings.large_leaf_width))
                starting_leaf = 14 if self.settings.hdmlc else 10 + 1
                leaves[0] = starting_leaf - leaf
            else:
                # offset MLC distance to the edge of the small leaf region
                mlc_dist += small_region_extent
                # divide the MLC distance by the leaf width and convert to leaf number
                leaf = int(round((abs(mlc_dist) + self.settings.large_leaf_width / 2) / self.settings.large_leaf_width))
                starting_leaf = 46 if self.settings.hdmlc else 50
                leaves[0] = starting_leaf + leaf

        # small leaf region
        else:
            # divide the MLC distance by the leaf width and convert to leaf number
            leaf = int(round((abs(mlc_dist) + self.settings.small_leaf_width / 2) / self.settings.small_leaf_width))
            if np.sign(mlc_dist) > 0:
                leaves[0] = 31 - leaf
            else:
                leaves[0] = 30 + leaf

        # set opposite leaf using an offset
        leaves[1] = 121 - leaves[0]

        return leaves
