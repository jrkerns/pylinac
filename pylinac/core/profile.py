import copy

import numpy as np
from scipy import ndimage

from pylinac.core.common_functions import peak_detect
from pylinac.core.decorators import value_accept
from pylinac.core.geometry import Point, Circle


def _sanitize_input(ydata, xdata):
    """Ensure data input is valid and/or convert to valid types."""
    # init xdata if none passed in
    if xdata is None:
        xdata = list(range(len(ydata)))

    if len(ydata) != len(xdata):
        raise ValueError("X and Y data are not same length")

    # convert to numpy arrays
    ydata = np.array(ydata, dtype=float)
    xdata = np.array(xdata)
    return ydata, xdata


class Profile(object):
    """A class for analyzing 1-D profiles that contain signals. These can be radiation beam profiles, i.e. with 1 large,
        smooth signal, or multiple signals, e.g. star lines in a circle profile.

    """
    def __init__(self, y_values=None, x_values=None):
        if y_values is not None:
            y_values, x_values = _sanitize_input(y_values, x_values)
        self.y_values = y_values
        self.x_values = x_values

    def filter_profile(self, size=100):
        """Filter the profile with a median filter."""
        self.y_values = ndimage.median_filter(self.y_values, size=size)

    def ground_profile(self):
        """Adjust the profile so that the lowest value is 0."""
        self.y_values = self.y_values - np.min(self.y_values)

    def find_peaks(self, min_peak_height=0.3, min_peak_distance=10, max_num_peaks=None, exclude_edge_portion=0.0):
        """Find the peaks of the profile using a simple max value search."""
        peak_vals, peak_idxs = peak_detect(self.y_values, self.x_values, threshold=min_peak_height, min_peak_width=min_peak_distance,
                                           max_num_peaks=max_num_peaks,exclude_edge_peaks=exclude_edge_portion)
        self.peaks = [Point(idx=peak_idx) for peak_idx in peak_idxs]

    def find_valleys(self, min_peak_height=0.3, min_peak_distance=10, max_num_peaks=None):
        """Find the valleys (minimums) of the profile using a simple min value search."""
        valley_vals, valley_idxs = peak_detect(self.y_values, self.x_values, threshold=min_peak_height, min_peak_width=min_peak_distance,
                                           max_num_peaks=max_num_peaks, find_min_instead=True)
        self.valleys = [Point(idx=valley_idx) for valley_idx in valley_idxs]

    def find_FWHM_peaks(self, min_peak_height=0.3, min_peak_distance=10, max_num_peaks=None):
        """Find "peaks" using the center of the FWHM (rather than by max value).

        This search is a bit more complicated than a max value search. Peaks are first determined using the max-value technique.
        Then, those values are used as initial starting points for the FHWM calculation.
        """
        self.find_peaks(min_peak_height, min_peak_distance, max_num_peaks)

        subprofiles = self._subdivide_profiles()

        # update peak points with modified indices
        for peak, profile in zip(self.peaks, subprofiles):
            fwhmc = int(profile.get_FWXM_center(70))
            peak.idx = fwhmc

    def _subdivide_profiles(self):
        """Subdivide the profile data into smaller pieces that can be analyzed for a single peak.

        :returns: A list of BeamProfiles
        """
        # append the peak list to include the endpoints of the profile
        peaks = self.peaks.copy()
        peaks.insert(0, Point(idx=0))
        peaks.append(Point(idx=len(self.y_values)))

        # create a list of new profiles from segments of original profile data.
        # New profiles are segmented by initial peak locations.
        subprofiles = []
        for idx in range(len(peaks)-2):
            left_end = peaks[idx].idx
            peak_idx = peaks[idx+1].idx
            right_end = peaks[idx+2].idx

            ydata = self.y_values[left_end:right_end]
            xdata = np.arange(left_end, right_end)

            subprofile = SingleProfile(ydata, xdata, initial_peak=peak_idx)
            subprofiles.append(subprofile)

        return subprofiles


class CircleProfile(Profile, Circle):
    """A profile in the shape of a circle."""

    def __init__(self, center=None, radius=None):
        Circle.__init__(self, center, radius)
        Profile.__init__(self)
        self.x_locs = np.ndarray  # x-values of the circle profile's location
        self.y_locs = np.ndarray  # y-values of the circle profile's location

    def get_profile(self, image_array, interval=0.0005):
        """Extracts values of a circular profile atop an image matrix.

        Extraction starts on the left side (180 degrees on unit circle), going clockwise.
        """

        self._ensure_array_size(image_array, self.radius+self.center.x, self.radius+self.center.y)

        # create index and cos, sin points which will be the circle's rectilinear coordinates
        rads = np.arange(0, 2 * np.pi, interval)
        x = np.cos(rads) * self.radius + self.center.x
        y = np.sin(rads) * self.radius + self.center.y

        # Pull the values of the image along the y,x points defined above, creating a circular profile
        profile = ndimage.map_coordinates(image_array, [y, x], order=0)

        self.y_values = profile
        self.x_locs = x
        self.y_locs = y

    def map_peaks(self):
        """Map found peaks to the x,y locations on the image/array; i.e. adds x,y coords to the peak locations"""
        for peak in self.peaks:
            peak.x = self.x_locs[peak.idx]
            peak.y = self.y_locs[peak.idx]

    def _ensure_array_size(self, array, min_width, min_height):
        width = array.shape[0]
        height = array.shape[1]
        if width < min_width or height < min_height:
            raise ValueError("Array size not large enough to compute profile")

class SingleProfile(object):
    """A profile that has one large signal, e.g. a radiation beam profile.

    Numerous signal analysis methods are given, mostly based on FWHM calculations.
    """

    def __init__(self, y_values, x_values=None, normalize_sides=True, initial_peak=None):
        """
        :param normalize_sides: Flag specifying whether to ground each side of the profile separately. If False, the profile will be
            grounded by the profile global minimum.
        :param initial_peak: If the approximate peak of the profile is known it can be passed in. Not needed unless there is more than
            one major peak in the profile.
        """
        self.y_values, self.x_values = _sanitize_input(y_values, x_values)

        # get an initial peak to work with
        self.initial_peak = self._get_initial_peak(initial_peak)

        if normalize_sides:
            self.ydata_right = self.y_values - np.min(y_values[self.initial_peak:-1])
            self.ydata_left = self.y_values - np.min(y_values[0:self.initial_peak])
            self.ymax_left = self.ydata_left[self.initial_peak]
            self.ymax_right = self.ydata_right[self.initial_peak]
        else:
            self.ydata_left = self.y_values
            self.ydata_right = self.y_values
            self.ymax_left = np.max(self.ydata_left)
            self.ymax_right = np.max(self.ydata_right)

    def _get_initial_peak(self, initial_peak):
        """Determine an initial peak to use as a rough guideline."""
        # if not passed, get one by peak searching.
        if initial_peak is None:
            _, initial_peak_arr = peak_detect(self.y_values, self.x_values, max_num_peaks=1, exclude_edge_peaks=0.2)
            try:
                initial_peak = initial_peak_arr[0]
            except IndexError:
                raise ValueError("A reasonable initial peak was not found in the profile. Ensure peak is not at profile edge")
        # otherwise use the one passed.
        else:
            # ensure peak is within the x_data region and not near an edge
            peak_beyond_left_edge = initial_peak < self.x_values[int(0.2 * len(self.x_values))]
            peak_beyond_right_edge = initial_peak > self.x_values[int(0.8 * len(self.x_values))]
            if peak_beyond_left_edge or peak_beyond_right_edge:
                raise IndexError("Initial peak that was passed was not reasonably withn the profile x_data range")

            initial_peak = np.where(self.x_values == initial_peak)[0]
        return int(initial_peak)

    @value_accept(side=('left', 'right'))
    def get_X_penum_idx(self, side='left', penum_point=50):
        """
        the main method of ProfPenum

        returns: the index of the point and the value of the point

        side indicates which side the method is looking from: left or right

        penumpoint indicates what value percentage is being sought. E.g. if
        penumpoint=20 the method looks for the index & value of the point along
        the 1D array that closest matches the 20% of maximum point on the given
        side
        """
        found = False
        peak = copy.copy(self.initial_peak)
        if side == 'left':
            thresh = self.ymax_left * penum_point / 100
            while not found:
                if self.ydata_left[peak] < thresh:
                    found = True
                    peak += 1
                elif peak == 0:
                    return None
                peak -= 1
        elif side == 'right':
            thresh = self.ymax_right * penum_point / 100
            while not found:
                if self.ydata_right[peak] < thresh:
                    found = True
                    peak -= 1
                elif peak == 0:
                    return None
                peak += 1
        else:
            raise TypeError("side was not correctly specified; use 'left' or 'right'")

        return self.x_values[peak]

    def get_FWXM(self, X=50):
        """
        get the Full-Width X-Max, where X is the percentage height.
        E.g. X = 50 is 50% height, a.k.a half max, thus X=50 is the FWHM
        returns the width in number of elements of the FWXM
        """
        li = self.get_X_penum_idx('left', X)
        ri = self.get_X_penum_idx('right', X)
        fwxm = np.abs(ri - li)
        return fwxm

    def get_FWXM_center(self, X=50, rounded=False):
        """
        returns the center point (index) of FWXM
        """
        fwxm = self.get_FWXM(X)
        li = self.get_X_penum_idx('left', X)
        fwxmcen = np.abs(li + fwxm / 2)
        if rounded:
            fwxmcen = np.round(fwxmcen)
        return fwxmcen

    def get_penum_width(self, side='left', lower_penum=20, upper_penum=80):
        """
        return the actual penumbral width of the profile. This is the
        standard "penumbra width" that med. phys. talks about in
        radiation profiles. Standard is the 80/20 width, although 90/10
        is sometimes used.

        side options include: left, right, & both (average)
        """
        if lower_penum > upper_penum:
            raise ValueError("Upper penumbra value must be larger than the lower penumbra value")

        if side == 'left':
            li = self.get_X_penum_idx('left', lower_penum)
            ui = self.get_X_penum_idx('left', upper_penum)
            pen = np.abs(ui - li)
            return pen
        elif side == 'right':
            li = self.get_X_penum_idx('right', lower_penum)
            ui = self.get_X_penum_idx('right', upper_penum)
            pen = np.abs(ui - li)
            return pen
        elif side == 'both':
            li = self.get_X_penum_idx('left', lower_penum)
            ui = self.get_X_penum_idx('left', upper_penum)
            lpen = np.abs(ui - li)
            li = self.get_X_penum_idx('right', lower_penum)
            ui = self.get_X_penum_idx('right', upper_penum)
            rpen = np.abs(ui - li)
            pen = np.mean([lpen, rpen])
            return pen
        else:
            raise NameError("getpenumwidth input parameter not acceptable")

    def get_field_value(self, field_width_percent=80, value='mean'):
        """Get the value of the field in the profile, either mean, median, or max, within the given field width.

        :param field_width_percent: The percent width relative to the FWHM to sample.
        :type field_width_percent: int < 100
        :param value: Value type to extract. Either 'mean', 'median', or 'max' of field values.
        :type value: str
        """

        fwhmc = self.get_FWXM_center()
        fwhm = self.get_FWXM() * 0.8
        left = fwhmc - fwhm / 2
        right = fwhmc + fwhm / 2

        field_values = self.y_values[left:right]

        if value == 'mean':
            return field_values.mean()
        elif value == 'median':
            return np.median(field_values)
        elif value == 'max':
            return field_values.max()
        elif value == 'min':
            return field_values.min()
