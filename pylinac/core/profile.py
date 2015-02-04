
"""Module of objects that resemble or contain a profile, i.e. a 1 or 2-D f(x) representation."""

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


class Profile:
    """A class for analyzing 1-D profiles that contain signals. Most often this should contain
        profiles with multiple peaks. Methods are mostly for *finding & filtering* the signals, peaks, valleys, etc.

        1-D profiles with a single peak (e.g. radiation beam profiles) are better suited by the SingleProfile class.
    """
    def __init__(self, y_values=None, x_values=None):
        """
        Parameters
        ----------
        y_values : iterable
            Y-values of the profile.
        x_values : iterable, optional
            X-values of the profile. If None, will create a array from the range of y-values.
        """
        if y_values is not None:
            y_values, x_values = _sanitize_input(y_values, x_values)
        self.y_values = y_values
        self.x_values = x_values

    def filter(self, size=0.05):
        """Filter the profile with a median filter.

        Parameters
        ----------
        size : int, float
            Size of the median filter to apply.
            If a float, the size is the ratio of the length. Must be in the range 0-1.
            E.g. if size=0.1 for a 1000-element array, the filter will be 100 elements.
            If an int, the filter is the size passed.
        """
        if isinstance(size, float):
            if 0 < size < 1:
                size *= len(self.y_values)
            else:
                raise TypeError("Float was passed but was not between 0 and 1")

        self.y_values = ndimage.median_filter(self.y_values, size=size)

    def ground(self):
        """Ground the profile such that the lowest value is 0.

        .. note::
            This will also "ground" profiles that are negative or partially-negative.
            For such profiles, be careful that this is the behavior you desire.
        """
        self.y_values = self.y_values - np.min(self.y_values)

    def find_peaks(self, min_peak_height=0.3, min_peak_distance=0.05, max_num_peaks=None, exclude_lt_edge=0.0, exclude_rt_edge=0.0,
                   return_it=False):
        """Find the peaks (maximums) of the profile using a simple maximum value search.

        Parameters
        ----------
        return_it : bool
            If False (default), peaks are saved as an attr of the profile.
            If True, peak values and peak indices are returned.

        Returns
        -------
        peak_vals : numpy.array, numpy.array
            The peak values and the peak indices.

        See Also
        --------
        common_functions.peak_detect : Further parameter info
        """
        peak_vals, peak_idxs = peak_detect(self.y_values, self.x_values, min_peak_height, min_peak_distance,
                                           max_num_peaks, exclude_lt_edge, exclude_rt_edge)
        self.peaks = [Point(value=peak_val, idx=peak_idx) for peak_idx, peak_val in zip(peak_idxs, peak_vals)]
        if return_it:
            return self._return_extrema_val_as_list('peak'), self._return_extrema_idx_as_list('peak')

    def find_valleys(self, min_peak_height=0.3, min_peak_distance=10, max_num_peaks=None, exclude_lt_edge=0.0, exclude_rt_edge=0.0,
                     return_it=False):
        """Find the valleys (minimums) of the profile using a simple minimum value search.

        Parameters
        ----------
        return_it : bool
            If False (default), peaks are saved as an attr of the profile.
            If True, valley values and valley indices are returned.

        Returns
        -------
        numpy.array, numpy.array
            Two arrays are returned if the return_it flag is true: The valley values and the valley indices.

        See Also
        --------
        common_functions.peak_detect : Further parameter info
        """
        valley_vals, valley_idxs = peak_detect(self.y_values, self.x_values, min_peak_height, min_peak_distance,
                                               max_num_peaks, exclude_lt_edge, exclude_rt_edge, find_min_instead=True)
        self.valleys = [Point(value=valley_val, idx=valley_idx) for valley_idx, valley_val in zip(valley_idxs, valley_vals)]
        if return_it:
            return self._return_extrema_val_as_list('valley'), self._return_extrema_idx_as_list('valley')

    def find_FWXM_peaks(self, fwxm=70, min_peak_height=0.3, min_peak_distance=10, max_num_peaks=None, return_it=False):
        """Find peaks using the center of the FWHM (rather than by max value).

        This search is a bit more complicated than a max value search. Peaks are first determined using the max-value technique.
        Then, those values are used as initial starting points for the FHWM calculation.

        Parameters
        ----------
        fwxm : int, float
            The Full-Width-X-Maximum desired. E.g. 0.7 will return the FW70%M.
            Values must be between 0 and 100.

        See Also
        --------
        find_peaks : Further parameter info
        common_functions.peak_detect : Further parameter info

        """
        self.find_peaks(min_peak_height, min_peak_distance, max_num_peaks)

        if not self.peaks:
            raise AttributeError("No peaks were found; try lowering the minimum peak height or use a different region.")
        subprofiles = self._subdivide_profiles()

        # update peak points with modified indices
        for peak, profile in zip(self.peaks, subprofiles):
            fwhmc = int(profile.get_FWXM_center(fwxm))
            peak.idx = fwhmc

        if return_it:
            return self._return_extrema_val_as_list('peak'), self._return_extrema_idx_as_list('peak')

    def _subdivide_profiles(self):
        """Subdivide the profile data into smaller pieces that can be analyzed for a single peak.

        Returns
        -------
        list
            SingleProfiles
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

    @value_accept(type=('peak', 'valley'))
    def _return_extrema_idx_as_list(self, type):
        """Return the peak indices as a list."""
        if type == 'peak':
            return [peak.idx for peak in self.peaks]
        else:
            return [valley.idx for valley in self.valleys]

    @value_accept(type=('peak', 'valley'))
    def _return_extrema_val_as_list(self, type):
        """Return the peak values as a list."""
        if type == 'peak':
            return [peak.value for peak in self.peaks]
        else:
            return [valley.value for valley in self.valleys]


class CircleProfile(Profile, Circle):
    """A profile in the shape of a circle.

        CircleProfile inherits from Profile, meaning that it has
        x_values and y_values. Additionally, however, because of the
        circular nature, x- and y-coordinates must also be determined.

        .. warning::
            Be sure to keep the attributes straight when dealing with a
            CircleProfile. The y_values attr is the, e.g., pixel values of an
            image; the x_values are the indices of the profile (0, 1,...);
            the x_locs and y_locs are the 2D coordinates of the y_ and x_values.
    """

    def __init__(self, center=None, radius=None):
        Circle.__init__(self, center, radius)
        Profile.__init__(self)
        self.x_locs = np.ndarray  # x-values of the circle profile's location
        self.y_locs = np.ndarray  # y-values of the circle profile's location

    def get_profile(self, image_array, size=1000, start=0, ccw=True):
        """Extracts a profile of an image matrix along the circle.

        Parameters
        ----------
        image_array : numpy.ndarray
            2D Numpy array
        size : int
            Size in elements of the desired profile.
        start : int, float
            Starting position of the profile; 0 is right (0 on unit circle).

        .. warning:: Units should be in radians. Setting degrees will give unpredictable results.

        ccw : bool
            If True (default), the profile will proceed counter-clockwise (the direction on the unit circle).
            If False, will proceed clockwise.

        See Also
        --------
        numpy.ndimage.map_coordinates : Further algorithm details
        """
        self._ensure_array_size(image_array, self.radius+self.center.x, self.radius+self.center.y)

        # create index and cos, sin points which will be the circle's rectilinear coordinates
        interval = (2 * np.pi) / size
        rads = np.arange(0+start, (2*np.pi)+start-interval, interval)
        if ccw:
            rads = rads[::-1]
        x = np.cos(rads) * self.radius + self.center.x
        y = np.sin(rads) * self.radius + self.center.y

        # Pull the values of the image along the y,x points defined above, creating a circular profile
        profile = ndimage.map_coordinates(image_array, [y, x], order=0)

        self.y_values = profile
        self.x_locs = x
        self.y_locs = y

    def find_peaks(self, min_peak_height=0.3, min_peak_distance=10, max_num_peaks=None, exclude_lt_edge=0.0, exclude_rt_edge=0.0,
                   return_it=False):
        """Overloads Profile to also map peak locations to the image."""
        super().find_peaks(min_peak_height, min_peak_distance, max_num_peaks, exclude_lt_edge, exclude_rt_edge)
        self._map_peaks()
        if return_it:
            return self._return_extrema_val_as_list('peak'), self._return_extrema_idx_as_list('peak')

    def find_valleys(self, min_peak_height=0.3, min_peak_distance=10, max_num_peaks=None, exclude_lt_edge=0.0, exclude_rt_edge=0.0,
                     return_it=False):
        """Overload Profile to also map valley locations to the image."""
        super().find_valleys(min_peak_height, min_peak_distance, max_num_peaks, exclude_lt_edge, exclude_rt_edge)
        self._map_peaks()
        if return_it:
            return self._return_extrema_val_as_list('valley'), self._return_extrema_idx_as_list('valley')

    def find_FWXM_peaks(self, fwxm=70, min_peak_height=0.3, min_peak_distance=10, max_num_peaks=None, return_it=False):
        """Overloads Profile to also map the peak locations to the image."""
        super().find_FWXM_peaks(fwxm, min_peak_height, min_peak_distance, max_num_peaks)
        self._map_peaks()
        if return_it:
            return self._return_extrema_val_as_list('peak'), self._return_extrema_idx_as_list('peak')

    def _map_peaks(self):
        """Map found peaks to the x,y locations on the image/array; i.e. adds x,y coords to the peak locations"""
        for peak in self.peaks:
            peak.x = self.x_locs[peak.idx]
            peak.y = self.y_locs[peak.idx]

    #TODO: move this to utilities
    def _ensure_array_size(self, array, min_width, min_height):
        """Ensure the array size of inputs are greater than the minimums."""
        height = array.shape[0]
        width = array.shape[1]
        if width < min_width or height < min_height:
            raise ValueError("Array size not large enough to compute profile")

class SingleProfile:
    """A profile that has one large signal, e.g. a radiation beam profile.

    Signal analysis methods are given, mostly based on FWXM calculations.
    """

    def __init__(self, y_values, x_values=None, normalize_sides=True, initial_peak=None):
        """
        Parameters
        ----------
        normalize_sides : bool, optional
            If True (default), each side of the profile will be grounded independently.
            If False, the profile will be grounded by the profile global minimum.
        initial_peak : int, optional
            If the approximate peak of the profile is known it can be passed in. Not needed unless there is more than
            one major peak in the profile, e.g. a very high edge.

        See Also
        --------
        Profile : Further parameter info
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

    def _get_initial_peak(self, initial_peak, exclusion_region=0.1):
        """Determine an initial peak to use as a rough guideline.

        Parameters
        ----------
        exclusion_region : float
            The ratio of the profile to exclude from either side of the profile when searching
            for the initial peak. Must be between 0 and 1.
        """
        # if not passed, get one by peak searching.
        if initial_peak is None:
            _, initial_peak_arr = peak_detect(self.y_values, self.x_values, max_num_peaks=1, exclude_lt_edge=exclusion_region,
                                              exclude_rt_edge=exclusion_region)
            try:
                initial_peak = initial_peak_arr[0]
            except IndexError:
                raise ValueError("A reasonable initial peak was not found in the profile. Ensure peak is not at profile edge")
        # otherwise use the one passed.
        else:
            # ensure peak is within the x_data region and not near an edge
            peak_near_left_edge = initial_peak < self.x_values[int(exclusion_region * len(self.x_values))]
            peak_near_right_edge = initial_peak > self.x_values[int((1-exclusion_region) * len(self.x_values))]
            if peak_near_left_edge or peak_near_right_edge:
                raise IndexError("Initial peak that was passed was not reasonably within the profile x_data range")

            initial_peak = np.where(self.x_values == initial_peak)[0]
        return int(initial_peak)

    @value_accept(side=('left', 'right'))
    def get_X_penum_idx(self, side='left', penum_point=50):
        """Return the index of the given penumbra.

        Parameters
        ----------
        side : {'left', 'right'}
            Which side to look for the penumbra.
        penum_point : int
            The penumbra value to search for. E.g. if passed 20, the method finds
            the index of 0.2*max profile value.

        Returns
        -------
        int
            The index of the penumbra value
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

        return self.x_values[peak]

    def get_FWXM(self, x=50, round=False):
        """Return the width at X-Max, where X is the percentage height.

        Parameters
        ----------
        x : int
            The percent height of the profile. E.g. x = 50 is 50% height,
            i.e. FWHM.
        round : bool
            If False (default), the index is returned as-is, even if it's a float
            If True, converts the index to an int.

        Returns
        -------
        float
            The width in number of elements of the FWXM
        """
        li = self.get_X_penum_idx('left', x)
        ri = self.get_X_penum_idx('right', x)
        fwxm = np.abs(ri - li)
        if round:
            fwxm = int(fwxm)
        return fwxm

    def get_FWXM_center(self, x=50, round=False):
        """Return the center index of the FWXM.

        See Also
        --------
        get_FWXM : Further parameter info
        """
        fwxm = self.get_FWXM(x)
        li = self.get_X_penum_idx('left', x)
        fwxmcen = np.abs(li + fwxm / 2)
        if round:
            fwxmcen = np.round(fwxmcen)
        return fwxmcen

    @value_accept(side=('left', 'right', 'both'))
    def get_penum_width(self, side='left', lower=20, upper=80):
        """Return the penumbra width of the profile.

        This is the standard "penumbra width" calculation that medical physics talks about in
        radiation profiles. Standard is the 80/20 width, although 90/10
        is sometimes used.

        Parameters
        ----------
        lower : int
            The "lower" penumbra value used to calculate penumbra. Must be lower than upper.
        upper : int
            The "upper" penumbra value used to calculate penumbra.
        side : {'left', 'right', 'both'}
            Which side of the profile to determined penumbra.
            If 'both', the left and right sides are averaged.

        Raises
        ------
        ValueError
            If lower penumbra is larger than upper penumbra
        """
        if lower > upper:
            raise ValueError("Upper penumbra value must be larger than the lower penumbra value")

        if side == 'left':
            li = self.get_X_penum_idx('left', lower)
            ui = self.get_X_penum_idx('left', upper)
            pen = np.abs(ui - li)
        elif side == 'right':
            li = self.get_X_penum_idx('right', lower)
            ui = self.get_X_penum_idx('right', upper)
            pen = np.abs(ui - li)
        elif side == 'both':
            li = self.get_X_penum_idx('left', lower)
            ui = self.get_X_penum_idx('left', upper)
            lpen = np.abs(ui - li)
            li = self.get_X_penum_idx('right', lower)
            ui = self.get_X_penum_idx('right', upper)
            rpen = np.abs(ui - li)
            pen = np.mean([lpen, rpen])

        return pen

    @value_accept(field_width=(0, 1))
    def get_field_calculation(self, field_width=0.8, calculation='mean'):
        """Calculate the value of the field in the profile.

        This function is useful for determining field symmetry and flatness.

        Parameters
        ----------
        field_width : float
            The field width size relative to the FWHM.
            E.g. 0.8 will calculate in the rage of FWHM * 0.8. Must be between 0 and 1.
        calculation : {'mean', 'median', 'max', 'min}
            Calculation to perform on the field values.
        """

        fwhmc = self.get_FWXM_center()
        field_width = self.get_FWXM() * field_width
        left = int(fwhmc - field_width / 2)
        right = int(fwhmc + field_width / 2)

        field_values = self.y_values[left:right]

        if calculation == 'mean':
            return field_values.mean()
        elif calculation == 'median':
            return np.median(field_values)
        elif calculation == 'max':
            return field_values.max()
        elif calculation == 'min':
            return field_values.min()
