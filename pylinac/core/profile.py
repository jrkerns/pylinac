"""Module of objects that resemble or contain a profile, i.e. a 1 or 2-D f(x) representation."""
import copy
from functools import lru_cache

import numpy as np
from scipy import ndimage
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.patches import Circle as mpl_Circle

from pylinac.core.decorators import value_accept
from pylinac.core.geometry import Point, Circle

LEFT = 'left'
RIGHT = 'right'


class SingleProfile:
    """A profile that has one large signal, e.g. a radiation beam profile.

    Signal analysis methods are given, mostly based on FWXM calculations.
    """
    interpolation_factor = 100
    interpolation_type = 'linear'
    _values = np.ndarray

    def __init__(self, values, normalize_sides=True, initial_peak=None):
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
        self.values = values
        self._passed_initial_peak = initial_peak
        self._normalize_sides = normalize_sides

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError("Values must be a numpy array")
        self._values = value.astype(float)

    @property
    def _right_side_min(self):
        """The minimum on the right side of the peak."""
        return self.values[self._initial_peak_idx:].min()

    @property
    def _left_side_min(self):
        """The minimum on the left side of the peak."""
        return self.values[:self._initial_peak_idx].min()

    @property
    def _values_right(self):
        """The "right side" y data."""
        if self._normalize_sides:
            return self.values - self._right_side_min
        else:
            return self._grounded_values

    @property
    def _values_left(self):
        """The "left side" y data."""
        if self._normalize_sides:
            return self.values - self._left_side_min
        else:
            return self._grounded_values

    @property
    @lru_cache()
    def _grounded_values(self):
        """Ground the profile such that the lowest value is 0.
        """
        min_val = self.values.min()
        grounded_values = self.values - min_val
        return grounded_values

    @property
    @lru_cache()
    def _initial_peak_idx(self):
        x_idx = self._get_initial_peak(self._passed_initial_peak)
        return x_idx

    @_initial_peak_idx.setter
    def _initial_peak_idx(self, value):
        self._passed_initial_peak = value

    def _get_initial_peak(self, initial_peak):
        """Determine an initial peak to use as a rough guideline.

        Parameters
        ----------
        exclusion_region : float
            The ratio of the profile to exclude from either side of the profile when searching
            for the initial peak. Must be between 0 and 1.
        """
        # if not passed, get one by peak searching.
        lf_edge = 0.2
        rt_edge = 0.8
        if initial_peak is None:
            while True:
                _, initial_peak_arr = peak_detect(self.values, max_num_peaks=1, search_region=(lf_edge, rt_edge))
                try:
                    initial_peak = initial_peak_arr[0]
                    break
                except IndexError:
                    lf_edge -= 0.01
                    rt_edge -= 0.01
                    if lf_edge < 0:
                        raise ValueError("A reasonable initial peak was not found in the profile. Ensure peak is not at profile edge")
        # otherwise use the one passed.
        elif len(self.values) < initial_peak < 0:
            raise IndexError("Initial peak that was passed was not reasonably within the profile x_data range")

        return initial_peak

    # @value_accept(side=(LEFT, RIGHT), type=('value', 'index'))
    # @lru_cache()
    def _penumbra_point(self, side='left', penum_point=50, interpolate=False, type='index'):
        """Return the index of the given penumbra. Search starts at the peak and moves index-by-index
        outward until the penumbra value is hit.

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
        # get peak
        peak = copy.copy(self._initial_peak_idx)
        peak = peak*self.interpolation_factor if interpolate else peak

        # get y-data
        if side == LEFT:
            y_data = self._values_left_interp if interpolate else self._values_left
        else:
            y_data = self._values_right_interp if interpolate else self._values_right

        # flip data if on left side; index search moves to the right
        # if side == LEFT:
        #     y_data = y_data[::-1]
        #     offset = 1 if not interpolate else self.interpolation_factor
        #     peak = len(y_data) - offset - peak

        # get threshold
        max_point = y_data.max()
        threshold = max_point * (penum_point / 100)

        # find the index, moving 1 element at a time until the value is encountered
        found = False
        at_end = False
        try:
            while not found and not at_end:
                if y_data[peak] < threshold:
                    found = True
                    peak -= 1 if side == RIGHT else -1
                elif peak == 0:
                    at_end = True
                peak += 1 if side == RIGHT else -1
        except IndexError:
            raise IndexError("The point of interest was beyond the profile; i.e. the profile may be cut off on the side")

        # peak = peak/self.interpolation_factor if interpolate else peak
        if type == 'value':
            return self._values_interp[peak] if interpolate else self.values[peak]
        elif type == 'index':
            if interpolate:
                peak /= self.interpolation_factor
            return peak

    @property
    @lru_cache()
    def _values_left_interp(self):
        """Interpolated values of the "left side" data."""
        ydata_f = interp1d(self._indices, self._values_left, kind=self.interpolation_type)
        y_data = ydata_f(self._indices_interp)
        return y_data

    @property
    @lru_cache()
    def _values_right_interp(self):
        """Interpolated values of the "right side" data."""
        ydata_f = interp1d(self._indices, self._values_right, kind=self.interpolation_type)
        y_data = ydata_f(self._indices_interp)
        return y_data

    @property
    @lru_cache()
    def _values_interp(self):
        ydata_f = interp1d(self._indices, self.values, kind=self.interpolation_type)
        y_data = ydata_f(self._indices_interp)
        return y_data

    @property
    @lru_cache()
    def _indices_interp(self):
        """Interpolated values of the x data."""
        # return np.arange(start=0, stop=len(self.values)+1, step=1/self.interpolation_factor)
        return np.linspace(start=0, stop=len(self.values)-1, num=(len(self.values)-1) * self.interpolation_factor)

    @property
    @lru_cache()
    def _indices(self):
        # return np.arange(start=0, stop=len(self.values)+1)
        return np.linspace(start=0, stop=len(self.values)-1, num=len(self.values))

    def fwxm(self, x=50, interpolate=False):
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
        li = self._penumbra_point('left', x, interpolate)
        ri = self._penumbra_point('right', x, interpolate)
        fwxm = np.abs(ri - li)
        return fwxm

    def fwxm_center(self, x=50, interpolate=False, type='index'):
        """Return the center index of the FWXM.

        See Also
        --------
        get_FWXM : Further parameter info
        """
        fwxm = self.fwxm(x, interpolate=interpolate)
        li = self._penumbra_point('left', x, interpolate)
        fwxmcen = np.abs(li + fwxm / 2)
        if type == 'value':
            return self.values[fwxmcen] if not interpolate else self._values_interp[int(fwxmcen*self.interpolation_factor)]
        else:
            return fwxmcen

    @value_accept(side=('left', 'right', 'both'))
    def penumbra_width(self, side='left', lower=20, upper=80, interpolate=False):
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

        if side in (LEFT, RIGHT):
            li = self._penumbra_point(side, lower, interpolate)
            ui = self._penumbra_point(side, upper, interpolate)
            pen = np.abs(ui - li)
        elif side == 'both':
            li = self._penumbra_point('left', lower, interpolate)
            ui = self._penumbra_point('left', upper, interpolate)
            lpen = np.abs(ui - li)
            li = self._penumbra_point('right', lower, interpolate)
            ui = self._penumbra_point('right', upper, interpolate)
            rpen = np.abs(ui - li)
            pen = np.mean([lpen, rpen])

        return pen

    def field_values(self, field_width=0.8):
        """Return the values of the profile for the given field width."""
        left, right = self.field_edges(field_width)
        field_values = self.values[left:right]
        return field_values

    def field_edges(self, field_width=0.8):
        """Return the indices of the field width edges."""
        fwhmc = self.fwxm_center()
        field_width *= self.fwxm()
        left = fwhmc - field_width / 2
        right = fwhmc + field_width / 2
        return left, right

    @value_accept(field_width=(0, 1))
    def field_calculation(self, field_width=0.8, calculation='mean'):
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

        field_values = self.field_values(field_width)

        if calculation == 'mean':
            return field_values.mean()
        elif calculation == 'median':
            return np.median(field_values)
        elif calculation == 'max':
            return field_values.max()
        elif calculation == 'min':
            return field_values.min()
        elif calculation == 'area':
            cax = self.fwxm_center()
            lt_area = field_values[:cax+1]
            rt_area = field_values[cax:]
            return lt_area, rt_area

    def plot(self):
        """Plot the profile"""
        plt.plot(self.values)
        plt.show()

    def filter(self, size=0.05, type='median'):
        """Filter the profile with a median or gaussian filter.

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
                size *= len(self.values)
                size = max(size, 1)
            else:
                raise TypeError("Float was passed but was not between 0 and 1")

        if type == 'median':
            self.values = ndimage.median_filter(self.values, size=size)
        elif type == 'gaussian':
            self.values = ndimage.gaussian_filter(self.values, sigma=size)


class MultiProfile:
    """A class for analyzing 1-D profiles that contain multiple signals. Methods are mostly for *finding & filtering* the signals, peaks, valleys, etc.

        1-D profiles with a single peak (e.g. radiation beam profiles) are better suited by the SingleProfile class.
    """
    def __init__(self, values):
        """
        Parameters
        ----------
        y_values : iterable
            Y-values of the profile.
        x_values : iterable, optional
            X-values of the profile. If None, will create a array from the range of y-values.
        """
        self.values = values
        self.peaks = []
        self.valleys = []

    def filter(self, size=0.05, type='median'):
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
                size *= len(self.values)
            else:
                raise TypeError("Float was passed but was not between 0 and 1")

        if type == 'median':
            self.values = ndimage.median_filter(self.values, size=size)
        elif type == 'gaussian':
            self.values = ndimage.gaussian_filter(self.values, sigma=size)

    def ground(self):
        """Ground the profile such that the lowest value is 0.

        .. note::
            This will also "ground" profiles that are negative or partially-negative.
            For such profiles, be careful that this is the behavior you desire.
        """
        min_val = self.values.min()
        self.values = self.values - min_val
        return min_val

    def plot(self, show_peaks=True):
        """Plot the profile"""
        fig, ax = plt.subplots()
        ax.plot(self.values)
        # ax.show()
        if show_peaks:
            peaks_x = [peak.idx for peak in self.peaks]
            peaks_y = [peak.value for peak in self.peaks]
            ax.plot(peaks_x, peaks_y, 'go')

    # def plot_peaks(self):
    #     if hasattr(self, 'peaks'):
    #         fig, ax = plt.subplots()
    #         ax.plot(self.values)
    #         peaks_x = [peak.idx for peak in self.peaks]
    #         peaks_y = [peak.value for peak in self.peaks]
    #         ax.plot(peaks_x, peaks_y, 'g+')
            # ax.show()

    def find_peaks(self, min_peak_height=0.3, min_peak_distance=0.05, max_num_peaks=None, search_region=(0.0, 1.0), type='index'):
        """Find the peaks (maximums) of the profile using a simple maximum value search.

        Returns
        -------
        peak_vals : numpy.array, numpy.array
            The peak values and the peak indices.

        See Also
        --------
        common_functions.peak_detect : Further parameter info
        """
        peak_vals, peak_idxs = peak_detect(self.values, min_peak_height, min_peak_distance,
                                           max_num_peaks, search_region=search_region)
        self.peaks = [Point(value=peak_val, idx=peak_idx) for peak_idx, peak_val in zip(peak_idxs, peak_vals)]

        if type == 'index':
            return peak_idxs
        else:
            return peak_vals

    def find_valleys(self, min_peak_height=0.3, min_peak_distance=10, max_num_peaks=None, search_region=(0.0, 1.0), type='index'):
        """Find the valleys (minimums) of the profile using a simple minimum value search.

        Returns
        -------
        numpy.array, numpy.array
            Two arrays are returned: The valley values and the valley indices.

        See Also
        --------
        common_functions.peak_detect : Further parameter info
        """
        valley_vals, valley_idxs = peak_detect(self.values, min_peak_height, min_peak_distance,
                                               max_num_peaks, search_region=search_region, find_min_instead=True)
        self.valleys = [Point(value=valley_val, idx=valley_idx) for valley_idx, valley_val in zip(valley_idxs, valley_vals)]

        if type == 'index':
            return valley_idxs
        else:
            return valley_vals

    def find_fwxm_peaks(self, x=50, min_peak_height=0.3, min_peak_distance=10, max_num_peaks=None, interpolate=False,
                        search_region=(0.0, 1.0), type='index', interpolation_factor=100, interpolation_type='linear'):
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
        self.find_peaks(min_peak_height, min_peak_distance, max_num_peaks, search_region=search_region)
        if not self.peaks:
            return [], []

        # subdivide the profiles into SingleProfile's
        subprofiles = self.subdivide(interpolation_factor, interpolation_type)

        # update peak indices with the FWHM value instead of maximum value
        original_peaks = copy.deepcopy(self.peaks)
        for num, (peak, profile) in enumerate(zip(self.peaks, subprofiles)):
            shift = original_peaks[num - 1].idx if num > 0 else 0
            # shift = sum(len(profile.values) for profile in subprofiles[:num])
            fwhmc = profile.fwxm_center(x, interpolate=interpolate)
            peak.idx = fwhmc + shift

        if type == 'index':
            return [peak.idx for peak in self.peaks]
        else:
            return [peak.value for peak in self.peaks]

    def subdivide(self, interpolation_factor=100, interpolation_type='linear'):
        """Subdivide the profile data into smaller pieces that can be analyzed for a single peak.

        Returns
        -------
        list
            SingleProfiles
        """
        # append the peak list to include the endpoints of the profile
        peaks = self.peaks.copy()
        peaks.insert(0, Point(idx=0))
        peaks.append(Point(idx=len(self.values)))

        # create a list of single profiles from segments of original profile data.
        # New profiles are segmented by initial peak locations.
        subprofiles = []
        for idx in range(len(peaks)-2):
            left_end = peaks[idx].idx
            peak_idx = peaks[idx+1].idx - left_end
            right_end = peaks[idx+2].idx

            ydata = self.values[left_end:right_end]

            subprofile = SingleProfile(ydata, initial_peak=peak_idx)
            subprofile.interpolation_factor = interpolation_factor
            subprofile.interpolation_type = interpolation_type
            subprofiles.append(subprofile)

        return subprofiles


class CircleProfile(MultiProfile, Circle):
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
    def __init__(self, center, radius, image_array, size=1000, start=0, ccw=True):
        Circle.__init__(self, center, radius)
        self._ensure_array_size(image_array, self.radius + self.center.x, self.radius + self.center.y)
        self.image_array = image_array
        self._size = size
        self._start = start
        self.ccw = ccw
        self._x_locations = None
        self._y_locations = None
        MultiProfile.__init__(self, self._profile)

    @property
    @lru_cache()
    def _radians(self):
        interval = (2 * np.pi) / self._size
        rads = np.arange(0 + self._start, (2 * np.pi) + self._start - interval, interval)
        if self.ccw:
            rads = rads[::-1]
        return rads

    @property
    def x_locations(self):
        if self._x_locations is None:
            return np.cos(self._radians) * self.radius + self.center.x
        else:
            return self._x_locations

    @x_locations.setter
    def x_locations(self, array):
        self._x_locations = array

    @property
    def y_locations(self):
        if self._y_locations is None:
            return np.sin(self._radians) * self.radius + self.center.y
        else:
            return self._y_locations

    @y_locations.setter
    def y_locations(self, array):
        self._y_locations = array

    @property
    @lru_cache()
    def _profile(self):
        return ndimage.map_coordinates(self.image_array, [self.y_locations, self.x_locations], order=0)

    # def _get_profile(self, image_array, size=1000, start=0, ccw=True):
    #     """Extracts a profile of an image matrix along the circle.
    #
    #     Parameters
    #     ----------
    #     image_array : numpy.ndarray
    #         2D Numpy array
    #     size : int
    #         Size in elements of the desired profile.
    #     start : int, float
    #         Starting position of the profile; 0 is right (0 on unit circle).
    #
    #     .. warning:: Units should be in radians. Setting degrees will give unpredictable results.
    #
    #     ccw : bool
    #         If True (default), the profile will proceed counter-clockwise (the direction on the unit circle).
    #         If False, will proceed clockwise.
    #
    #     See Also
    #     --------
    #     numpy.ndimage.map_coordinates : Further algorithm details
    #     """

    def find_peaks(self, min_peak_height=0.3, min_peak_distance=10, max_num_peaks=None, search_region=(0.0, 1.0), type='index'):
        """Overloads Profile to also map peak locations to the image."""
        array = super().find_peaks(min_peak_height, min_peak_distance, max_num_peaks, search_region, type)
        self._map_peaks()
        return array

    def find_valleys(self, min_peak_height=0.3, min_peak_distance=10, max_num_peaks=None, search_region=(0.0, 1.0), type='index'):
        """Overload Profile to also map valley locations to the image."""
        array = super().find_valleys(min_peak_height, min_peak_distance, max_num_peaks, search_region, type)
        self._map_peaks()
        return array

    def find_fwxm_peaks(self, x=50, min_peak_height=0.3, min_peak_distance=10, max_num_peaks=None, interpolate=False, search_region=(0.0, 1.0), type='index', interpolation_type='linear'):
        """Overloads Profile to also map the peak locations to the image."""
        array = super().find_fwxm_peaks(x, min_peak_height, min_peak_distance, max_num_peaks, interpolate, search_region=search_region, type=type, interpolation_type=interpolation_type)
        self._map_peaks()
        return array

    def _map_peaks(self):
        """Map found peaks to the x,y locations on the image/array; i.e. adds x,y coordinates to the peak locations"""
        for peak in self.peaks:
            peak.x = self.x_locations[peak.idx]
            peak.y = self.y_locations[peak.idx]

    def roll_profile(self, amount):
        # Roll the profile and x and y coordinates
        self.values = np.roll(self.values, -amount)
        self.x_locations = np.roll(self.x_locations, -amount)
        self.y_locations = np.roll(self.y_locations, -amount)

    def add_to_axes(self, axes=None, edgecolor='black', fill=False, plot_peaks=True):
        """Add 2 circles to the axes: one at the maximum and minimum radius of the ROI."""
        if axes is None:
            fig, axes = plt.subplots()
            axes.imshow(self.image_array)
        axes.add_patch(
            mpl_Circle((self.center.x, self.center.y), edgecolor=edgecolor, radius=self.radius, fill=fill))
        if plot_peaks:
            x_locs = [peak.x for peak in self.peaks]
            y_locs = [peak.y for peak in self.peaks]
            axes.autoscale(enable=False)
            axes.scatter(x_locs, y_locs, s=40, marker='x', c=edgecolor)

    @staticmethod
    def _ensure_array_size(array, min_width, min_height):
            """Ensure the array size of inputs are greater than the minimums."""
            height = array.shape[0]
            width = array.shape[1]
            if width < min_width or height < min_height:
                raise ValueError("Array size not large enough to compute profile")


class CollapsedCircleProfile(CircleProfile):
    """A circular profile that collapses nearby pixels along the circle to create an composite profile.
    I.e. instead of simply profiling along the circle itself, it collapses the tangential pixels on the inner and
    outer sides for every point.
    """
    def __init__(self, center, radius, image_array, size=1000, start=0, ccw=True, width_ratio=0.1, num_profiles=20):
        Circle.__init__(self, center, radius)
        self._ensure_array_size(image_array, self.radius + self.center.x, self.radius + self.center.y)
        self.image_array = image_array
        self._size = size
        self._start = start
        self.ccw = ccw
        self.width_ratio = width_ratio
        self.num_profiles = num_profiles
        self._x_locations = None
        self._y_locations = None
        MultiProfile.__init__(self, self._profile)

    @property
    def _radii(self):
        return np.nditer(np.linspace(start=self.radius * (1 - self.width_ratio), stop=self.radius * (1 + self.width_ratio),
                  num=self.num_profiles))

    @property
    def _multi_x_locations(self):
        x = []
        cos = np.cos(self._radians)
        # extract profile for each circle radii
        for radius in self._radii:
            x.append(cos * radius + self.center.x)
        return x

    @property
    def _multi_y_locations(self):
        y = []
        sin = np.sin(self._radians)
        # extract profile for each circle radii
        for radius in self._radii:
            y.append(sin * radius + self.center.y)
        return y

    @property
    def _profile(self):
        profile = np.zeros(len(self._multi_x_locations[0]))
        for radius, x, y in zip(self._radii, self._multi_x_locations, self._multi_y_locations):
            profile += ndimage.map_coordinates(self.image_array, [y, x], order=0)
        profile /= self.num_profiles
        return profile

    def add_to_axes(self, axes=None, edgecolor='black', fill=False, plot_peaks=True):
        """Add 2 circles to the axes: one at the maximum and minimum radius of the ROI."""
        if axes is None:
            fig, axes = plt.subplots()
            axes.imshow(self.image_array)
        axes.add_patch(mpl_Circle((self.center.x, self.center.y), edgecolor=edgecolor, radius=self.radius*(1+self.width_ratio),
                                  fill=fill))
        axes.add_patch(mpl_Circle((self.center.x, self.center.y), edgecolor=edgecolor, radius=self.radius*(1-self.width_ratio),
                                  fill=fill))
        if plot_peaks:
            x_locs = [peak.x for peak in self.peaks]
            y_locs = [peak.y for peak in self.peaks]
            axes.autoscale(enable=False)
            axes.scatter(x_locs, y_locs, s=40, marker='x', c=edgecolor)


def peak_detect(values, threshold=0, min_peak_width=10, max_num_peaks=None, search_region=(0.0, 1.0), find_min_instead=False):
    """Find the peaks or valleys of a 1D signal.

    Uses the difference (np.diff) in signal to find peaks. Current limitations include:
        1) Only for use in 1-D data; 2D may be possible with the gradient function.
        2) Will not detect peaks at the very edge of array (i.e. 0 or -1 index)

    Parameters
    ----------
    y : array-like
        1D y-data of signal.
    x : array-like
        1D x-data of signal. If None, will create a uniform range the length of the y-data.
    threshold : int, float
        The value the peak must be above to be considered a peak. This removes "peaks"
        that are in a low-value region.
        If passed an int, the actual value is the threshold.
        E.g. when passed 15, any peak less with a value <15 is removed.
        If passed a float, it will threshold as a percent. Must be between 0 and 1.
        E.g. when passed 0.4, any peak <40% of the maximum value will be removed.
    min_peak_width : int, float
        If passed an int, parameter is the number of elements apart a peak must be from neighboring peaks.
        If passed a float, must be between 0 and 1 and represents the ratio of the profile to exclude.
        E.g. if passed 0.05 with a 1000-element profile, the minimum peak width will be 0.05*1000 = 50 elements.
    max_num_peaks : int
        Specify up to how many peaks will be returned. E.g. if 3 is passed in and 5 peaks are found, only the 3 largest
        peaks will be returned.
    find_min_instead : bool
        If False (default), peaks will be returned.
        If True, valleys will be returned.

    Returns
    -------
    max_vals : numpy.array
        The values of the peaks found.
    max_idxs : numpy.array
        The x-indices (locations) of the peaks.

    Raises
    ------
    ValueError
        If float not between 0 and 1 passed to threshold.
    """
    peak_vals = []  # a list to hold the y-values of the peaks. Will be converted to a numpy array
    peak_idxs = []  # ditto for x-values (index) of y data.

    if find_min_instead:
        values = -values

    """Limit search to search region"""
    left_end = search_region[0]
    if isinstance(left_end, float):
        left_index = int(left_end*len(values))
    elif isinstance(left_end, int):
        left_index = left_end

    right_end = search_region[1]
    if isinstance(right_end, float):
        right_index = int(right_end * len(values))
    elif isinstance(right_end, int):
        right_index = right_end

    values = values[left_index:right_index]

    y_diff = np.diff(values.astype(float))  # y and y_diff must be converted to signed type.

    if isinstance(threshold, float) and threshold < 1:
        data_range = values.max() - values.min()
        threshold = threshold * data_range + values.min()
    elif isinstance(threshold, float) and threshold >= 1:
        raise ValueError("When threshold is passed a float, value must be less than 1")

    """Find all potential peaks"""
    for idx in range(len(y_diff) - 1):
        # For each item of the diff array, check if:
        # 1) The y-value is above the threshold.
        # 2) The value of y_diff is positive (negative for valley search), it means the y-value changed upward.
        # 3) The next y_diff value is zero or negative (or positive for valley search); a positive-then-negative diff value means the value
        # is a peak of some kind. If the diff is zero it could be a flat peak, which still counts.

        # 1)
        if values[idx + 1] < threshold:
            continue

        y1_gradient = y_diff[idx] > 0
        y2_gradient = y_diff[idx + 1] <= 0

        # 2) & 3)
        if y1_gradient and y2_gradient:
            # If the next value isn't zero it's a single-pixel peak. Easy enough.
            if y_diff[idx + 1] != 0:
                peak_vals.append(values[idx + 1])
                peak_idxs.append(idx + 1 + left_index)
            # elif idx >= len(y_diff) - 1:
            #     pass
            # Else if the diff value is zero, it could be a flat peak, or it could keep going up; we don't know yet.
            else:
                # Continue on until we find the next nonzero diff value.
                try:
                    shift = 0
                    while y_diff[(idx + 1) + shift] == 0:
                        shift += 1
                        if (idx + 1 + shift) >= (len(y_diff) - 1):
                            break
                    # If the next diff is negative (or positive for min), we've found a peak. Also put the peak at the center of the flat
                    # region.
                    is_a_peak = y_diff[(idx + 1) + shift] < 0
                    if is_a_peak:
                        peak_vals.append(values[(idx + 1) + np.round(shift / 2)])
                        peak_idxs.append((idx + 1 + left_index) + np.round(shift / 2))
                except IndexError:
                    pass

    # convert to numpy arrays
    peak_vals = np.array(peak_vals)
    peak_idxs = np.array(peak_idxs)
    # try:
    #     peak_idxs += l_edge
    # except:
    #     pass

    """Enforce the min_peak_distance by removing smaller peaks."""
    # For each peak, determine if the next peak is within the min peak width range.
    if isinstance(min_peak_width, float):
        if 0 > min_peak_width >= 1:
            raise ValueError("When min_peak_width is passed a float, value must be between 0 and 1")
        else:
            min_peak_width = int(min_peak_width * len(values))

    index = 0
    while index < len(peak_idxs) - 1:

        # If the second peak is closer than min_peak_distance to the first peak, find the larger peak and remove the other one.
        if peak_idxs[index] > peak_idxs[index + 1] - min_peak_width:
            if peak_vals[index] > peak_vals[index + 1]:
                idx2del = index + 1
            else:
                idx2del = index
            peak_vals = np.delete(peak_vals, idx2del)
            peak_idxs = np.delete(peak_idxs, idx2del)
        else:
            index += 1

    """If Maximum Number passed, return only up to number given based on a sort of peak values."""
    if max_num_peaks is not None and len(peak_idxs) > max_num_peaks:
        sorted_peak_vals = peak_vals.argsort()  # sorts low to high
        peak_vals = peak_vals[sorted_peak_vals[-max_num_peaks:]]
        peak_idxs = peak_idxs[sorted_peak_vals[-max_num_peaks:]]

    # If we were looking for minimums, convert the values back to the original sign
    if find_min_instead:
        peak_vals = -peak_vals

    return peak_vals, peak_idxs


if __name__ == '__main__':
    import scipy.signal as sps
    xdata = np.linspace(0, 1.7 * np.pi, num=200)
    ydata = sps.sawtooth(xdata, width=0.5)
    s = SingleProfile(ydata, normalize_sides=True)
    s._penumbra_point('right', 50, interpolate=True, type='value')
    s.fwxm_center()
