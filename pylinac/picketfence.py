
"""The picketfence module is used for loading and analyzing EPID images of a "picket fence", a common MLC
pattern produced when performing linac QA."""
import os.path as osp

import numpy as np
import scipy.ndimage.filters as spfilt
from scipy import signal
import matplotlib.pyplot as plt

from pylinac.core.decorators import lazyproperty
from pylinac.core.geometry import Line, Rectangle
from pylinac.core.io import get_filepath_UI
from pylinac.core.profile import Profile
from pylinac.core.image import Image

orientations = {'UD': 'Up-Down', 'LR': 'Left-Right'}  # possible orientations of the pickets. UD is up-down, LR is left-right.


class PicketFence:
    """A class used for analyzing EPID images where radiation strips have been formed by the
    MLCs. The strips are assumed to be parallel to one another and normal to the image edge;
    i.e. a "left-right" or "up-down" orientation is assumed. Further work could follow up by accounting
    for any angle.

    Attributes
    ----------
    pickets: list
        Holds :class:`~pylinac.picketfence.Picket` objects.
    image: :class:`~pylinac.core.image.Image` object.

    Examples
    --------
    Run the demo::
        >>> PicketFence().run_demo()

    Typical session:
        >>> img_path = r"C:/QA/June/PF.dcm"  # the EPID image
        >>> mypf = PicketFence()
        >>> mypf.load_image(img_path)
        >>> mypf.analyze(tolerance=0.5, action_tolerance=0.3)
        >>> print(mypf.return_results())
        >>> mypf.plot_analyzed_image()
    """
    def __init__(self):
        self.pickets = []
        self._action_lvl = None

    def _clear_attrs(self):
        """Clear attributes; necessary when new image loaded or analysis done on same image."""
        self.pickets = []
        self._action_lvl = None

    @property
    def passed(self):
        """Boolean specifying if all MLC positions were within tolerance."""
        for picket in self.pickets:
            for meas in picket.mlc_meas:
                if not meas.passed:
                    return False
        return True

    @property
    def percent_passing(self):
        """Return the percentage of MLC positions under tolerance."""
        num = 0
        num_pass = 0
        for picket in self.pickets:
            num += len(picket._error_array)
            num_pass += sum(picket._error_array < picket._tolerance)
        pct_pass = 100 * num_pass / num
        return pct_pass

    @property
    def max_error(self):
        """Return the maximum error found."""
        max_error = 0
        for idx, picket in enumerate(self.pickets):
            if picket.max_error > max_error:
                max_error = picket.max_error
        return max_error

    @property
    def max_error_picket(self):
        """Return the picket number where the maximum error occured."""
        max_error = 0
        where_at = 0
        for idx, picket in enumerate(self.pickets):
            if picket.max_error > max_error:
                max_error = picket.max_error
                where_at = idx
        return where_at

    @property
    def max_error_leaf(self):
        """Return the leaf that had the maximum error."""
        picket = self.pickets[self.max_error_picket]
        return np.argmax(picket._error_array)

    @property
    def abs_median_error(self):
        """Return the median error found."""
        median_error = []
        for picket in self.pickets:
            median_error.append(picket.abs_median_error)
        return max(median_error)

    @property
    def _action_lvl_set(self):
        if self._action_lvl is not None:
            return True
        else:
            return False

    @property
    def num_pickets(self):
        """Return the number of pickets determined."""
        return len(self.pickets)

    def load_demo_image(self):
        """Load the demo image that is included with pylinac."""
        im_open_path = osp.join(osp.dirname(__file__), 'demo_files', 'picket_fence', 'EPID-PF-LR.dcm')
        self.load_image(im_open_path)

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
            self.image.array = spfilt.median_filter(self.image.array, size=filter)
        self._clear_attrs()

    def load_image_UI(self):
        """Load the image using a UI dialog box."""
        path = get_filepath_UI()
        self.load_image(path)

    def run_demo(self, tolerance=0.5):
        """Run the Picket Fence demo using the demo image. See analyze() for parameter info."""
        self.load_demo_image()
        self.analyze(tolerance)
        print(self.return_results())
        self.plot_analyzed_image()

    def analyze(self, tolerance=0.5, action_tolerance=None, hdmlc=False):
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
            to indicate an "action" is necessary on the part of the physicist to
            resolve the issue.
        hdmlc : bool
            If False (default), a standard (5/10mm leaves) Millennium MLC model is assumed.
            If True, an HD (2.5/5mm leaves) Millennium is assumed.
        """
        if action_tolerance is not None and tolerance < action_tolerance:
            raise ValueError("Tolerance cannot be lower than the action tolerance")

        """Pre-analysis"""
        self._clear_attrs()
        self._action_lvl = action_tolerance
        self._check_inversion()
        self._threshold()
        self._find_orientation()

        """Analysis"""
        self._construct_pickets(tolerance, action_tolerance)
        leaf_centers = self._find_leaf_centers(hdmlc)
        self.calc_mlc_positions(leaf_centers)
        self.calc_mlc_error()

    def _construct_pickets(self, tolerance, action_tolerance):
        """Construct the Picket instances."""
        if self.orientation == orientations['UD']:
            leaf_prof = np.median(self._analysis_array, 0)
        else:
            leaf_prof = np.median(self._analysis_array, 1)
        leaf_prof = Profile(leaf_prof)
        _, peak_idxs = leaf_prof.find_peaks(min_peak_distance=0.01, min_peak_height=0.5)
        for peak in range(len(peak_idxs)):
            self.pickets.append(Picket(self.image, tolerance, self.orientation, action_tolerance))

    def _find_leaf_centers(self, hdmlc):
        """Return the leaf centers perpendicular to the leaf motion."""
        # generate some settings
        sm_lf_wdth = 5 * self.image.dpmm
        bg_lf_wdth = sm_lf_wdth * 2
        if hdmlc:
            sm_lf_wdth /= 2
            bg_lf_wdth /= 2
        self._sm_lf_meas_wdth = slmw = int(round(sm_lf_wdth*3/4))
        self._bg_lf_meas_wdth = blmw = int(round(bg_lf_wdth*3/4))
        bl_ex = int(bg_lf_wdth/4)
        sm_ex = int(sm_lf_wdth/4)

        # generate leaf profile
        if self.orientation == orientations['UD']:
            leaf_prof = np.mean(self._analysis_array, 1)
            center = self.image.center.y
        else:
            leaf_prof = np.mean(self._analysis_array, 0)
            center = self.image.center.x
        leaf_prof = Profile(leaf_prof)

        # ground profile to reasonable level
        _, peak_idxs = leaf_prof.find_peaks(min_peak_distance=self._sm_lf_meas_wdth, exclude_lt_edge=sm_ex,
                                            exclude_rt_edge=sm_ex)
        min_val = leaf_prof.y_values[peak_idxs[0]:peak_idxs[-1]].min()
        leaf_prof.y_values[leaf_prof.y_values < min_val] = min_val

        # remove unevenness in signal
        leaf_prof.y_values = signal.detrend(leaf_prof.y_values, bp=[int(len(leaf_prof.y_values)/3), int(len(leaf_prof.y_values)*2/3)])
        _, peak_idxs = leaf_prof.find_peaks(min_peak_distance=self._sm_lf_meas_wdth, exclude_lt_edge=sm_ex, exclude_rt_edge=sm_ex)
        leaf_range = (peak_idxs[-1] - peak_idxs[0]) / self.image.dpmm  # mm
        sm_lf_range = 220  # mm

        # find leaf peaks
        if leaf_range > sm_lf_range:
            lt_biglittle_lf_bndry = int(round(center - 100 * self.image.dpmm))
            rt_biglittle_lf_bndry = int(round(center + 100 * self.image.dpmm))
            pp = leaf_prof.subdivide([lt_biglittle_lf_bndry, rt_biglittle_lf_bndry], slmw)
            if len(pp) != 3:
                raise ValueError("3 Profiles weren't found but should have been")
            # Left Big MLC region
            _, peak_idxs = pp[0].find_peaks(min_peak_distance=blmw, exclude_lt_edge=bl_ex)
            peak_diff = np.diff(peak_idxs).mean()
            lt_v_idx = np.array(peak_idxs[:-1]) + peak_diff/2

            # Middle, small MLC region
            _, peak_idxs = pp[1].find_peaks(min_peak_distance=slmw)
            peak_diff = np.diff(peak_idxs).mean()
            mid_v_idx = np.array(peak_idxs[:-1]) + peak_diff / 2

            # Right Big MLC region
            _, peak_idxs = pp[2].find_peaks(min_peak_distance=blmw,
                                            exclude_rt_edge=bl_ex)
            peak_diff = np.diff(peak_idxs).mean()
            rt_v_idx = np.array(peak_idxs[:-1]) + peak_diff / 2
            leaf_center_idxs = np.concatenate((lt_v_idx, mid_v_idx, rt_v_idx))
        else:
            _, peak_idxs = leaf_prof.find_peaks(min_peak_distance=slmw, exclude_lt_edge=sm_ex,
                                                exclude_rt_edge=sm_ex)
            _, peak_idxs = leaf_prof.find_FWXM_peaks(min_peak_distance=slmw, interpolate=True)
            peak_diff = np.diff(peak_idxs).mean()
            leaf_center_idxs = np.array(peak_idxs[:-1]) + peak_diff / 2
        return leaf_center_idxs

    def calc_mlc_positions(self, leaf_centers):
        """Calculate the positions of all the MLC pairs."""
        diff = np.diff(leaf_centers)
        sample_width = np.round(np.median(diff*2/5)/2).astype(int)

        for mlc_num, mlc_peak_loc in enumerate(np.round(leaf_centers).astype(int)):
            mlc_rows = np.arange(mlc_peak_loc-sample_width, mlc_peak_loc+sample_width+1)
            if self.orientation == orientations['UD']:
                pix_vals = np.median(self._analysis_array[mlc_rows, :], axis=0)
            else:
                pix_vals = np.median(self._analysis_array[:, mlc_rows], axis=1)
            prof = Profile(pix_vals)
            prof.find_FWXM_peaks(fwxm=80, min_peak_distance=0.01, min_peak_height=0.5, interpolate=True)
            for idx, peak in enumerate(prof.peaks):
                if self.orientation == orientations['UD']:
                    meas = MLC_Meas((peak.idx, mlc_rows[0]), (peak.idx, mlc_rows[-1]))
                else:
                    meas = MLC_Meas((mlc_rows[0], peak.idx), (mlc_rows[-1], peak.idx))
                self.pickets[idx].mlc_meas.append(meas)

    def calc_mlc_error(self):
        """Calculate the error of the MLC positions relative to the picket fit."""
        for picket in self.pickets:
            picket.fit_poly()
            picket.calc_mlc_errors()

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

                    if not mlc_meas.passed:
                        color = 'r'
                    elif self._action_lvl_set and not mlc_meas.passed_action:
                        color = 'm'
                    else:
                        color = 'b'
                    mlc_meas.add_to_axes(ax.axes, color=color, width=1.5)
        # plot the overlay if desired.
        if overlay:
            for mlc_num, mlc in enumerate(self.pickets[0].mlc_meas):

                below_tol = True
                if self._action_lvl_set:
                    below_action = True
                for picket in self.pickets:
                    if not picket.mlc_passed(mlc_num):
                        below_tol = False
                    if self._action_lvl_set and not picket.mlc_passed_action(mlc_num):
                        below_action = False
                if below_tol:
                    if self._action_lvl_set and not below_action:
                        color = 'm'
                    else:
                        color = 'g'
                else:
                    color = 'r'
                if self.orientation == orientations['UD']:
                    r = Rectangle(max(self.image.shape)*2, self._sm_lf_meas_wdth, (mlc.center.x, mlc.center.y))
                else:
                    r = Rectangle(self._sm_lf_meas_wdth, max(self.image.shape) * 2, (mlc.center.x, mlc.center.y))
                r.add_to_axes(ax.axes, edgecolor='none', fill=True, alpha=0.1, facecolor=color)

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
        """Print results of analysis."""
        pass_pct = self.percent_passing
        string = "Picket Fence Results: \n{:2.1f}% " \
                 "Passed\nMedian Error: {:2.3f}mm \n" \
                 "Max Error: {:2.3f}mm on Picket: {}, Leaf: {}".format(pass_pct, self.abs_median_error, self.max_error,
                                                                                                   self.max_error_picket,
                                                                                                  self.max_error_leaf)
        return string

    def _check_inversion(self):
        """Check the image for inversion (pickets are valleys, not peaks) by sampling the 4 image corners.
        If the average value of the four corners is above the average pixel value, then it is very likely inverted.
        """
        outer_edge = 10
        inner_edge = 30
        TL_corner = self.image.array[outer_edge:inner_edge, outer_edge:inner_edge]
        BL_corner = self.image.array[-inner_edge:-outer_edge, -inner_edge:-outer_edge]
        TR_corner = self.image.array[outer_edge:inner_edge, outer_edge:inner_edge]
        BR_corner = self.image.array[-inner_edge:-outer_edge, -inner_edge:-outer_edge]
        corner_avg = np.mean((TL_corner, BL_corner, TR_corner, BR_corner))
        if corner_avg > np.mean(self.image.array.flatten()):
            self.image.invert()

    def _threshold(self):
        """Threshold the image by subtracting the minimum value. Allows for more accurate image orientation determination.
        """
        col_prof = np.median(self.image.array, 0)
        col_prof = Profile(col_prof)
        row_prof = np.median(self.image.array, 1)
        row_prof = Profile(row_prof)
        _, r_peak_idx = row_prof.find_peaks(min_peak_distance=0.01, exclude_lt_edge=0.05, exclude_rt_edge=0.05)
        _, c_peak_idx = col_prof.find_peaks(min_peak_distance=0.01, exclude_lt_edge=0.05, exclude_rt_edge=0.05)
        min_val = self.image.array[r_peak_idx[0]:r_peak_idx[-1], c_peak_idx[0]:c_peak_idx[-1]].min()
        self._analysis_array = self.image.array.copy()
        self._analysis_array[self._analysis_array < min_val] = min_val
        self._analysis_array -= min_val

    # @property
    # def _analysis_array(self):
    #     return getattr(self, '_aa', self.image.array.copy())
    #
    # @_analysis_array.setter
    # def _analysis_array(self, array):
    #     if array.shape != self.image.shape:
    #         raise ValueError("Array size is not the same as the original image")
    #     self._aa = array

    def _find_orientation(self):
        """Determine the orientation of the radiation strips by examining percentiles of the sum of each axes of the image.
        A high standard deviation is a surrogate for the axis the pickets are along.
        """
        row_sum = np.sum(self._analysis_array, 0)
        col_sum = np.sum(self._analysis_array, 1)
        row80, row90 = np.percentile(row_sum, [80, 90])
        col80, col90 = np.percentile(col_sum, [80, 90])
        row_range = row90 - row80
        col_range = col90 - col80
        # The true picket side will have a greater difference in
        # percentiles than will the non-picket size.
        if row_range < col_range:
            self.orientation = orientations['LR']
        else:
            self.orientation = orientations['UD']


class Picket:
    """Holds *Picket* information in a Picket Fence test."""
    def __init__(self, image, tolerance, orientation, action_tolerance):
        """
        Attributes
        ----------
        mlc_meas : list
            Holds :class:`~pylinac.picketfence.MLC_Meas` objects.
        fit : numpy.poly1d
            The fit equation of the picket.
        """
        self.mlc_meas = []
        self._img = image
        self._tolerance = tolerance
        self._orientation = orientation
        self._action_tol = action_tolerance

    @property
    def mm_per_pixel(self):
        """The mm per pixel of the image."""
        return 1/self._img.dpmm

    @property
    def abs_median_error(self):
        """The absolute median error of the MLC measurements."""
        return np.median(np.abs(self._error_array))

    @property
    def max_error(self):
        """The max error of the MLC measurements."""
        return self._error_array.max()

    @lazyproperty
    def _error_array(self):
        err = []
        for meas in self.mlc_meas:
            err.append(meas.error)
        return np.array(err)

    def mlc_passed(self, mlc):
        """Return whether the MLC has passed tolerance."""
        if self.mlc_meas[mlc].passed:
            return True
        else:
            return False

    def mlc_passed_action(self, mlc):
        """Return whether the MLC has passed the action tolerance."""
        if self._action_tol is not None:
            if self.mlc_meas[mlc].passed_action:
                return True
            else:
                return False
        else:
            raise AttributeError("No action tolerance was specified")

    def fit_poly(self):
        """Fit a polynomial to the MLC measurements; also constructs the guard rails."""
        x = np.array([mlc.point1.y for mlc in self.mlc_meas])
        y = np.array([mlc.point1.x for mlc in self.mlc_meas])
        if self._orientation == orientations['UD']:
            fit = np.polyfit(x, y, 1)
        else:
            fit = np.polyfit(y, x, 1)
        self.fit = np.poly1d(fit)
        self._fit_guard_rails(fit)

    def _fit_guard_rails(self, fit):
        l_fit = np.copy(fit)
        l_fit[-1] += self._tolerance/self.mm_per_pixel
        self.left_guard = np.poly1d(l_fit)
        r_fit = np.copy(fit)
        r_fit[-1] -= self._tolerance/self.mm_per_pixel
        self.right_guard = np.poly1d(r_fit)

    def calc_mlc_errors(self):
        """Calculate the MLC error from the picket."""
        for idx, meas in enumerate(self.mlc_meas):
            if self._orientation == orientations['UD']:
                picket_pos = self.fit(meas.point1.y)
                mlc_pos = meas.point1.x
            else:
                picket_pos = self.fit(meas.point1.x)
                mlc_pos = meas.point1.y
            physical_error = (mlc_pos - picket_pos) * self.mm_per_pixel
            meas.error = physical_error
            if physical_error < self._tolerance:
                meas.passed = True
            if self._action_tol is not None:
                if physical_error < self._action_tol:
                    meas.passed_action = True

    def add_guards_to_axes(self, axis, color='g'):
        """Plot guard rails to the axis."""
        if self._orientation == orientations['UD']:
            length = self._img.shape[0]
        else:
            length = self._img.shape[1]
        x_data = np.arange(length)
        left_y_data = self.left_guard(x_data)
        right_y_data = self.right_guard(x_data)
        if self._orientation == orientations['UD']:
            axis.plot(left_y_data, x_data, color=color)
            axis.plot(right_y_data, x_data, color=color)
        else:
            axis.plot(x_data, left_y_data, color=color)
            axis.plot(x_data, right_y_data, color=color)


class MLC_Meas(Line):
    """Represents an MLC measurement."""
    def __init__(self, point1=None, point2=None, m=None, b=None):
        """
        Attributes
        ----------
        error : float
            The error of the MLC measurement and the picket fit.
        passed : bool
            Whether the MLC measurement was under tolerance.
        passed_action : bool
            Whether the MLC measurement was under the action level tolerance.
        """
        super().__init__(point1, point2, m, b)
        self.error = 0
        self.passed = False
        self.passed_action = False


# -----------------------------------
# Picket Fence Demo
# -----------------------------------
if __name__ == '__main__':
    # from scipy.ndimage.interpolation import rotate
    # import cProfile
    # cProfile.run('PicketFence().run_demo()', sort=1)
    # PicketFence().run_demo()
    pf = PicketFence()
    # pf.open_UI()
    pf.load_demo_image()
    # pf.image.rot90()
    # pf.image.array = rotate(pf.image.array, 0.5, reshape=False, mode='nearest')
    pf.analyze(tolerance=0.15, action_tolerance=0.03)
    print(pf.return_results())
    pf.plot_analyzed_image()