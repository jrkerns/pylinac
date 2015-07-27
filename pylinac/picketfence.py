
"""The picket fence module is meant for analyzing EPID images where a "picket fence" MLC pattern has been made.
Physicists regularly check MLC positioning through this test. This test can be done using film and one can
"eyeball" it, but this is the 21st century and we have numerous ways of quantifying such data. This module
attains to be one of them. It will load in an EPID dicom image and determine the MLC peaks, error of each MLC
pair to the picket, and give a few visual indicators for passing/warning/failing.

Features:

* **Analyze either HD or regular MLCs** - Just pass a flag and tell pylinac whether it's HD or not.
* **Easy-to-read pass/warn/fail overlay** - Analysis gives you easy-to-read tools for determining the status of an MLC pair.
* **Any Source-to-Image distance** - Whatever your clinic uses as the SID for picket fence, pylinac can account for it.
* **Account for panel translation** - Have an off-CAX setup? No problem. Translate your EPID and pylinac knows.
* **Account for panel sag** - If your EPID sags at certain angles, just tell pylinac and the results will be shifted.
"""
import os.path as osp
from functools import lru_cache
from io import BytesIO

import numpy as np
import matplotlib.pyplot as plt

from pylinac.core.geometry import Line, Rectangle
from pylinac.core.io import get_filepath_UI
from pylinac.core.profile import Profile, SingleProfile
from pylinac.core.image import Image

orientations = {'UD': 'Up-Down', 'LR': 'Left-Right'}  # possible orientations of the pickets. UD is up-down, LR is left-right.


class PicketFence:
    """A class used for analyzing EPID images where radiation strips have been formed by the
    MLCs. The strips are assumed to be parallel to one another and normal to the image edge;
    i.e. a "left-right" or "up-down" orientation is assumed. Further work could follow up by accounting
    for any angle.

    Attributes
    ----------
    pickets: :class:`~pylinac.picketfence.PicketHandler`
    image: :class:`~pylinac.core.image.Image`

    Examples
    --------
    Run the demo::
        >>> PicketFence().run_demo()

    Typical session:
        >>> img_path = r"C:/QA/June/PF.dcm"  # the EPID image
        >>> mypf = PicketFence(img_path)
        >>> mypf.analyze(tolerance=0.5, action_tolerance=0.3)
        >>> print(mypf.return_results())
        >>> mypf.plot_analyzed_image()
    """
    def __init__(self, filename=None, filter=None):
        """
        Parameters
        ----------
        filename : str, None
            Name of the file as a string. If None, image must be loaded later.
        filter : int, None
            The filter size to apply to the image upon load.
        """
        if filename is not None:
            self.load_image(filename, filter)

    @classmethod
    def from_url(cls, url, filter=None):
        """Instantiate from a URL.

        .. versionadded:: 0.7.1
        """
        obj = cls()
        obj.load_url(url, filter=filter)
        return obj

    def load_url(self, url, filter=None):
        """Load from a URL.

        .. versionadded:: 0.7.1
        """
        try:
            import requests
        except ImportError:
            raise ImportError("Requests is not installed; cannot get the log from a URL")
        response = requests.get(url)
        if response.status_code != 200:
            raise ConnectionError("Could not connect to the URL")
        stream = BytesIO(response.content)
        self.load_image(stream, filter=filter)

    @property
    def passed(self):
        """Boolean specifying if all MLC positions were within tolerance."""
        return self.pickets.passed

    @property
    def percent_passing(self):
        """Return the percentage of MLC positions under tolerance."""
        num = 0
        num_pass = 0
        for picket in self.pickets:
            num += len(picket.error_array)
            num_pass += sum(picket.error_array < self.settings.tolerance)
        pct_pass = 100 * num_pass / num
        return pct_pass

    @property
    def max_error(self):
        """Return the maximum error found."""
        return max(picket.max_error for picket in self.pickets)

    @property
    def max_error_picket(self):
        """Return the picket number where the maximum error occured."""
        return np.argmax([picket.max_error for picket in self.pickets])

    @property
    def max_error_leaf(self):
        """Return the leaf that had the maximum error."""
        picket = self.pickets[self.max_error_picket]
        return np.argmax(picket.error_array)

    @property
    @lru_cache()
    def abs_median_error(self):
        """Return the median error found."""
        return np.median(np.hstack([picket.error_array for picket in self.pickets]))

    @property
    def num_pickets(self):
        """Return the number of pickets determined."""
        return len(self.pickets)

    @classmethod
    def from_demo_image(cls, filter=None):
        """Construct a PicketFence instance using the demo image.

        .. versionadded:: 0.6
        """
        obj = cls()
        obj.load_demo_image(filter=filter)
        return obj

    def load_demo_image(self, filter=None):
        """Load the demo image that is included with pylinac."""
        im_open_path = osp.join(osp.dirname(__file__), 'demo_files', 'picket_fence', 'EPID-PF-LR.dcm')
        self.load_image(im_open_path, filter=filter)

    def load_image(self, file_path, filter=None):
        """Load the image

        Parameters
        ----------
        file_path : str
            Path to the image file.
        filter : int, None
            If None (default), no filtering will be done to the image.
            If an int, will perform median filtering over image of size *filter*.
        """
        self.image = Image(file_path)
        if isinstance(filter, int):
            self.image.median_filter(size=filter)
        self._check_for_noise()
        self.image.check_inversion()

    @classmethod
    def from_image_UI(cls, filter=None):
        """Construct a PicketFence instance and load an image using a dialog box.

        .. versionadded:: 0.6
        """
        obj = cls()
        obj.load_image_UI(filter=filter)
        return obj

    def load_image_UI(self, filter=None):
        """Load the image using a UI dialog box."""
        path = get_filepath_UI()
        self.load_image(path, filter=filter)

    def _check_for_noise(self):
        """Check if the image has extreme noise (dead pixel, etc) by comparing
        min/max to 1/99 percentiles and smoothing if need be."""
        while self._has_noise():
            self.image.median_filter()

    def _has_noise(self):
        """Helper method to determine if there is spurious signal in the image."""
        min = self.image.array.min()
        max = self.image.array.max()
        near_min, near_max = np.percentile(self.image.array, [0.5, 99.5])
        max_is_extreme = max > near_max * 2
        min_is_extreme = (min < near_min) and (abs(near_min - min) > 0.2 * near_max)
        return max_is_extreme or min_is_extreme

    def _adjust_for_sag(self, sag):
        """Roll the image to adjust for EPID sag."""
        sag_pixels = int(round(sag * self.settings.dpmm))
        direction = 'y' if self.orientation == orientations['UD'] else 'x'
        self.image.roll(direction, sag_pixels)

    def run_demo(self, tolerance=0.5):
        """Run the Picket Fence demo using the demo image. See analyze() for parameter info."""
        self.load_demo_image()
        self.analyze(tolerance)
        print(self.return_results())
        self.plot_analyzed_image()

    def analyze(self, tolerance=0.5, action_tolerance=None, hdmlc=False, num_pickets=None, sag_adjustment=0):
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
            only needed if analysis is catching things that aren't pickets.
        sag_adjustment : float, int

            .. versionadded:: 0.8

            The amount of shift in mm to apply to the image to correct for EPID sag.
            For Up-Down picket images, positive moves the image down, negative up.
            For Left-Right picket images, positive moves the image left, negative right.
        """
        if action_tolerance is not None and tolerance < action_tolerance:
            raise ValueError("Tolerance cannot be lower than the action tolerance")

        """Pre-analysis"""
        self.settings = Settings(self.orientation, tolerance, action_tolerance, hdmlc, self.image.dpmm)
        self._adjust_for_sag(sag_adjustment)

        """Analysis"""
        self.pickets = PicketHandler(self.image, self.settings, num_pickets)

    def plot_analyzed_image(self, guard_rails=True, mlc_peaks=True, overlay=True, show=True):
        """Plot the analyzed image.

        Parameters
        ----------
        guard_rails : bool
            Do/don't plot the picket "guard rails".
        mlc_peaks : bool
            Do/don't plot the MLC positions.
        overlay : bool
            Do/don't plot the alpha overlay of the leaf status.
        """
        # plot the image
        plt.clf()
        ax = plt.imshow(self.image.array, cmap=plt.cm.Greys)

        # plot guard rails and mlc peaks as desired
        for p_num, picket in enumerate(self.pickets):
            if guard_rails:
                picket.add_guards_to_axes(ax.axes)
            if mlc_peaks:
                for idx, mlc_meas in enumerate(picket.mlc_meas):
                    mlc_meas.add_to_axes(ax.axes, width=1.5)
        # plot the overlay if desired.
        if overlay:
            o = Overlay(self.image, self.settings, self.pickets)
            o.add_to_axes(ax)

        plt.xlim([0, self.image.shape[1]])
        plt.ylim([0, self.image.shape[0]])
        plt.axis('off')

        if show:
            plt.show()

    def save_analyzed_image(self, filename, guard_rails=True, mlc_peaks=True, overlay=True, **kwargs):
        """Save the analyzed figure to a file."""
        self.plot_analyzed_image(guard_rails, mlc_peaks, overlay, show=False)
        plt.savefig(filename, **kwargs)

    def return_results(self):
        """Return results of analysis. Use with print()."""
        pass_pct = self.percent_passing
        string = "Picket Fence Results: \n{:2.1f}% " \
                 "Passed\nMedian Error: {:2.3f}mm \n" \
                 "Max Error: {:2.3f}mm on Picket: {}, Leaf: {}".format(pass_pct, self.abs_median_error, self.max_error,
                                                                                                   self.max_error_picket,
                                                                                                  self.max_error_leaf)
        return string

    @property
    def orientation(self):
        """The orientation of the image, either Up-Down or Left-Right."""
        # replace any dead pixels with median value
        temp_image = self.image.array.copy()
        temp_image[temp_image < np.median(temp_image)] = np.median(temp_image)

        # find "range" of 80 to 90th percentiles
        row_sum = np.sum(temp_image, 0)
        col_sum = np.sum(temp_image, 1)
        row80, row90 = np.percentile(row_sum, [80, 90])
        col80, col90 = np.percentile(col_sum, [80, 90])
        row_range = row90 - row80
        col_range = col90 - col80

        # The true picket side will have a greater difference in
        # percentiles than will the non-picket size.
        if row_range < col_range:
            orientation = orientations['LR']
        else:
            orientation = orientations['UD']
        return orientation


class Overlay:
    """Class for handling the "overlay" feature of the plot."""
    def __init__(self, image, settings, pickets):
        self.image = image
        self.settings = settings
        self.pickets = pickets

    def add_to_axes(self, axes):
        rect_width = self.pickets[0].sample_width*2
        for mlc_num, mlc in enumerate(self.pickets[0].mlc_meas):
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
            if self.settings.orientation == orientations['UD']:
                r = Rectangle(self.image.shape[1], rect_width, center=(self.image.center.x, mlc.center.y))
            else:
                r = Rectangle(rect_width, self.image.shape[0], center=(mlc.center.y, self.image.center.y))
            r.add_to_axes(axes.axes, edgecolor='none', fill=True, alpha=0.1, facecolor=color)


class Settings:
    """Simple class to hold various settings for PF."""
    def __init__(self, orientation, tolerance, action_tolerance, hdmlc, dpmm):
        self.orientation = orientation
        self.tolerance = tolerance
        self.action_tolerance = action_tolerance
        self.hdmlc = hdmlc
        self.dpmm = dpmm
        self.mmpd = 1/dpmm


class PicketHandler:
    """Finds and handles the pickets of the image."""
    def __init__(self, image_array, settings, num_pickets):
        self.pickets = []
        self.image_array = image_array
        self.settings = settings
        self.num_pickets = num_pickets
        self.find_pickets()

    def find_pickets(self):
        """Find the pickets of the image."""
        leaf_prof = self.image_mlc_inplane_mean_profile
        _, peak_idxs = leaf_prof.find_peaks(min_peak_distance=0.02, min_peak_height=0.5, max_num_peaks=self.num_pickets)
        peak_spacing = np.median(np.diff(peak_idxs))

        for peak_idx in peak_idxs:
            self.pickets.append(Picket(self.image_array, self.settings, peak_idx, peak_spacing/2))

    @property
    def passed(self):
        """Whether all the pickets passed tolerance."""
        return all(picket.passed for picket in self)

    def __getitem__(self, item):
        return self.pickets[item]

    def __len__(self):
        return len(self.pickets)

    @property
    def image_mlc_inplane_mean_profile(self):
        """A profile of the image along the MLC travel direction."""
        if self.settings.orientation == orientations['UD']:
            leaf_prof = np.mean(self.image_array, 0)
        else:
            leaf_prof = np.mean(self.image_array, 1)
        return Profile(leaf_prof)


class Picket:
    """Holds *Picket* information in a Picket Fence test."""
    def __init__(self, image, settings, approximate_idx, spacing):
        """
        Attributes
        ----------
        mlc_meas : list
            Holds :class:`~pylinac.picketfence.MLCMeas` objects.
        """
        self.mlc_meas = []
        self.image = image
        self.settings = settings
        self.approximate_idx = approximate_idx
        self.spacing = spacing
        self._get_mlc_positions()

    def _get_mlc_positions(self):
        """Calculate the positions of all the MLC pairs."""
        # for each MLC...
        for mlc_num, mlc_center in enumerate(self.leaf_centers):
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
        if self.settings.orientation == orientations['UD']:
            pix_vals = np.median(self.picket_array[mlc_rows, :], axis=0)
        else:
            pix_vals = np.median(self.picket_array[:, mlc_rows], axis=1)
        if max(pix_vals) > np.percentile(self.picket_array, 80):
            prof = SingleProfile(pix_vals)
            fw80mc = prof.get_FWXM_center(70, interpolate=True)
            return fw80mc + self.approximate_idx - self.spacing

    def add_mlc_meas(self, mlc_center, mlc_position):
        """Add an MLC measurement point."""
        upper_point = mlc_center - self.sample_width / 2
        lower_point = mlc_center + self.sample_width / 2

        if self.settings.orientation == orientations['UD']:
            meas = MLCMeas((mlc_position, upper_point), (mlc_position, lower_point), self.settings)
        else:
            meas = MLCMeas((upper_point, mlc_position), (lower_point, mlc_position), self.settings)
        self.mlc_meas.append(meas)

    @property
    def sample_width(self):
        """The width to sample the MLC leaf (~40% of the leaf width)."""
        return np.round(np.median(np.diff(self.leaf_centers) * 2 / 5) / 2).astype(int)

    @property
    @lru_cache()
    def picket_array(self):
        """A slice of the whole image that contains the area around the picket."""
        if self.settings.orientation == orientations['UD']:
            array = self.image.array[:, self.approximate_idx - self.spacing:self.approximate_idx + self.spacing]
        else:
            array = self.image.array[self.approximate_idx - self.spacing:self.approximate_idx + self.spacing, :]
        return array

    @property
    @lru_cache()
    def small_leaf_width(self):
        leaf_width_mm = 5
        leaf_width_pixels = leaf_width_mm * self.settings.dpmm
        if self.settings.hdmlc:
            leaf_width_pixels /= 2
        return leaf_width_pixels

    @property
    def large_leaf_width(self):
        return self.small_leaf_width * 2

    @property
    def number_small_leaves(self):
        return 40 if not self.settings.hdmlc else 32

    @property
    def number_large_leaves(self):
        return 20 if not self.settings.hdmlc else 28

    @property
    @lru_cache()
    def leaf_centers(self):
        """Return a set of leaf centers perpendicular to the leaf motion based on the position of the CAX."""
        # generate a set of leaf center points based on physical widths of large and small leaves
        first_shift = self.large_leaf_width * (self.number_large_leaves / 2 - 1) + self.large_leaf_width * 0.75
        second_shift = self.small_leaf_width * (self.number_small_leaves - 1) + self.large_leaf_width * 0.75

        large_leaf_section = np.arange(self.number_large_leaves / 2) * self.large_leaf_width
        small_leaf_section = (np.arange(self.number_small_leaves) * self.small_leaf_width) + first_shift
        large_leaf_section2 = (np.arange(self.number_large_leaves / 2) * self.large_leaf_width) + first_shift + second_shift
        leaf_centers = np.concatenate((large_leaf_section, small_leaf_section, large_leaf_section2))

        # now adjust them to align with the iso
        if self.settings.orientation == orientations['UD']:
            leaf30_center = self.image.cax.y - self.small_leaf_width / 2
            edge = self.image.shape[0]
        else:
            leaf30_center = self.image.cax.x - self.small_leaf_width / 2
            edge = self.image.shape[1]
        adjustment = leaf30_center - leaf_centers[29]
        leaf_centers += adjustment

        # only include values that are reasonable as values might extend past image (e.g. with small SID)
        values_in_image = (leaf_centers > 0 + self.large_leaf_width/2) & (leaf_centers < edge - self.large_leaf_width/2)
        leaf_centers = leaf_centers[values_in_image]
        return np.round(leaf_centers).astype(int)

    @property
    def abs_median_error(self):
        """The absolute median error of the MLC measurements."""
        return np.median(np.abs(self.error_array))

    @property
    def max_error(self):
        """The max error of the MLC measurements."""
        return self.error_array.max()

    @property
    @lru_cache()
    def error_array(self):
        """An array containing the error values of all the measurements."""
        return np.array([meas.error for meas in self.mlc_meas])

    @property
    def passed(self):
        """Whether or not all the measurements passed."""
        return all(meas.passed for meas in self.mlc_meas)

    def mlc_passed(self, mlc):
        """Return whether a specific MLC has passed tolerance."""
        return self.mlc_meas[mlc].passed

    def mlc_passed_action(self, mlc):
        """Return whether a specific MLC has passed the action tolerance."""
        if self.settings.action_tolerance is not None:
            return self.mlc_meas[mlc].passed_action
        else:
            raise AttributeError("No action tolerance was specified")

    @property
    def fit(self):
        """The fit of a polynomial to the MLC measurements."""
        x = np.array([mlc.point1.y for mlc in self.mlc_meas])
        y = np.array([mlc.point1.x for mlc in self.mlc_meas])
        if self.settings.orientation == orientations['UD']:
            fit = np.polyfit(x, y, 1)
        else:
            fit = np.polyfit(y, x, 1)
        return np.poly1d(fit)

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

    def add_guards_to_axes(self, axis, color='g'):
        """Plot guard rails to the axis."""
        if self.settings.orientation == orientations['UD']:
            length = self.image.shape[0]
        else:
            length = self.image.shape[1]
        x_data = np.arange(length)
        left_y_data = self.left_guard(x_data)
        right_y_data = self.right_guard(x_data)
        if self.settings.orientation == orientations['UD']:
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

    def add_to_axes(self, axes, width=1, color='w'):
        """Plot the measurement."""
        super().add_to_axes(axes, width, color=self.bg_color)

    @property
    def bg_color(self):
        """The color of the measurement when the PF image is plotted."""
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
    def error(self):
        """The error (difference) of the MLC measurement and the picket fit at that location."""
        if self.settings.orientation == orientations['UD']:
            picket_pos = self.fit(self.center.y)
            mlc_pos = self.center.x
        else:
            picket_pos = self.fit(self.center.x)
            mlc_pos = self.center.y
        return abs(mlc_pos - picket_pos) * self.settings.mmpd

    @property
    def passed(self):
        """Whether the MLC measurement was under tolerance."""
        return self.error < self.settings.tolerance

    @property
    def passed_action(self):
        """Whether the MLC measurement was under the action level tolerance."""
        if self.settings.action_tolerance is not None:
            return self.error < self.settings.action_tolerance


# -----------------------------------
# Picket Fence Demo
# -----------------------------------
if __name__ == '__main__':
    pf = PicketFence.from_demo_image()
    pf.analyze()
    print(pf.return_results())
    pf.plot_analyzed_image()
