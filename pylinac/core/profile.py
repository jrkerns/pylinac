"""Module of objects that resemble or contain a profile, i.e. a 1 or 2-D f(x) representation."""
from functools import lru_cache
from typing import Union, Tuple, Sequence, List, Optional

import argue
import numpy as np
from matplotlib.patches import Circle as mpl_Circle
import matplotlib.pyplot as plt
from scipy import ndimage, signal
from scipy.interpolate import interp1d

from .geometry import Point, Circle
from .typing import NumberLike


def stretch(array: np.ndarray, min: int=0, max: int=1, fill_dtype: Optional[np.dtype]=None) -> np.ndarray:
    """'Stretch' the profile to the fit a new min and max value and interpolate in between.
    From: http://www.labri.fr/perso/nrougier/teaching/numpy.100/  exercise #17

    Parameters
    ----------
    array: numpy.ndarray
        The numpy array to stretch.
    min : number
        The new minimum of the values.
    max : number
        The new maximum value.
    fill_dtype : numpy data type
        If None (default), the array will be stretched to the passed min and max.
        If a numpy data type (e.g. np.int16), the array will be stretched to fit the full range of values
        of that data type. If a value is given for this parameter, it overrides ``min`` and ``max``.
    """
    new_max = max
    new_min = min
    if fill_dtype is not None:
        try:
            di = np.iinfo(fill_dtype)
        except ValueError:
            di = np.finfo(fill_dtype)
        new_max = di.max
        new_min = di.min
    # perfectly normalize the array (0..1). ground, then div by range
    stretched_array = (array - array.min())/(array.max() - array.min())
    # stretch normalized array to new max/min
    stretched_array *= new_max
    # stretched_array += new_min
    if fill_dtype:
        stretched_array = stretched_array.astype(fill_dtype)
    return stretched_array


class ProfileMixin:
    """A mixin to provide various manipulations of 1D profile data."""
    values: np.ndarray

    def invert(self) -> None:
        """Invert (imcomplement) the profile."""
        orig_array = self.values
        self.values = -orig_array + orig_array.max() + orig_array.min()

    def normalize(self, norm_val: Union[str, NumberLike]='max') -> None:
        """Normalize the profile to the given value.

        Parameters
        ----------
        norm_val : str, number
            If a string, must be 'max', which normalizes the values to the maximum value.
            If a number, normalizes all values to that number.
        """
        if norm_val == 'max':
            val = self.values.max()
        else:
            val = norm_val
        self.values /= val

    def stretch(self, min: NumberLike=0, max: NumberLike=1) -> None:
        """'Stretch' the profile to the min and max parameter values.

        Parameters
        ----------
        min : number
            The new minimum of the values
        max : number
            The new maximum value.
        """
        self.values = stretch(self.values, min=min, max=max)

    def ground(self) -> float:
        """Ground the profile such that the lowest value is 0.

        Returns
        -------
        float
            The minimum value that was used as the grounding value.
        """
        min_val = self.values.min()
        self.values = self.values - min_val
        return min_val

    @argue.options(kind=('median', 'gaussian'))
    def filter(self, size: NumberLike=0.05, kind: str='median') -> None:
        """Filter the profile.

        Parameters
        ----------
        size : float, int
            Size of the median filter to apply.
            If a float, the size is the ratio of the length. Must be in the range 0-1.
            E.g. if size=0.1 for a 1000-element array, the filter will be 100 elements.
            If an int, the filter is the size passed.
        kind : {'median', 'gaussian'}
            The kind of filter to apply. If gaussian, `size` is the sigma value.
        """
        if isinstance(size, float):
            if 0 < size < 1:
                size = int(round(len(self.values)*size))
                size = max(size, 1)
            else:
                raise TypeError("Float was passed but was not between 0 and 1")

        if kind == 'median':
            self.values = ndimage.median_filter(self.values, size=size)
        elif kind == 'gaussian':
            self.values = ndimage.gaussian_filter(self.values, sigma=size)

    def __len__(self):
        return len(self.values)

    def __getitem__(self, items):
        return self.values[items]


class SingleProfile(ProfileMixin):
    """A profile that has one large signal, e.g. a radiation beam profile.
    Signal analysis methods are given, mostly based on FWXM calculations.
    Profiles with multiple peaks are better suited by the MultiProfile class.
    """
    interpolation_factor: int = 100
    interpolation_type: str = 'linear'
    _values: np.ndarray

    def __init__(self, values: np.ndarray):
        """
        Parameters
        ----------
        values : ndarray
            The profile numpy array. Must be 1D.
        """
        self.values = values

    @property
    def values(self) -> np.ndarray:
        """The profile array."""
        return self._values

    @values.setter
    def values(self, value):
        if not isinstance(value, np.ndarray):
            raise TypeError("Values must be a numpy array")
        self._values = value.astype(float)

    @property
    @lru_cache()
    def _values_interp(self) -> np.ndarray:
        """Interpolated values of the entire profile array."""
        ydata_f = interp1d(self._indices, self.values, kind=self.interpolation_type)
        y_data = ydata_f(self._indices_interp)
        return y_data

    @property
    def _indices_interp(self) -> np.ndarray:
        """Interpolated values of the profile index data."""
        return np.linspace(start=0, stop=len(self.values)-1, num=(len(self.values)-1) * self.interpolation_factor)

    @property
    def _indices(self) -> np.ndarray:
        """Values of the profile index data."""
        return np.linspace(start=0, stop=len(self.values)-1, num=len(self.values))

    @argue.bounds(x=(0, 100))
    def fwxm(self, x: int = 50) -> float:
        """Return the width at X-Max, where X is the percentage height.

        Parameters
        ----------
        x : int
            The percent height of the profile. E.g. x = 50 is 50% height,
            i.e. FWHM.

        Returns
        -------
        int, float
            The width in number of elements of the FWXM.
        """
        _, peak_props = find_peaks(self.values, fwxm_height=x/100, max_number=1)
        return peak_props['widths'][0]

    def fwxm_center(self, x: int=50, interpolate: bool=False) -> Tuple[NumberLike, NumberLike]:
        """Return the center of the FWXM.

        See Also
        --------
        fwxm() : Further parameter info
        """
        _, peak_props = find_peaks(self.values, fwxm_height=x/100, max_number=1)
        fwxm_center_idx = (peak_props['right_ips'][0] - peak_props['left_ips'][0])/2 + peak_props['left_ips'][0]
        if interpolate:
            return fwxm_center_idx, self._values_interp[int(fwxm_center_idx*self.interpolation_factor)]
        else:
            fwxm_center_idx = int(round(fwxm_center_idx))
            return fwxm_center_idx, self.values[fwxm_center_idx]

    @argue.bounds(lower=(0, 100), upper=(0, 100))
    def penumbra_width(self, lower: int=20, upper: int=80) -> Tuple[float, float]:
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

        Raises
        ------
        ValueError
            If lower penumbra is larger than upper penumbra
        """
        if lower > upper:
            raise ValueError("Upper penumbra value must be larger than the lower penumbra value")

        _, upper_peak_props = find_peaks(self.values, fwxm_height=upper/100, max_number=1)
        _, lower_peak_props = find_peaks(self.values, fwxm_height=lower/100, max_number=1)
        left_penum = np.abs(upper_peak_props['left_ips'][0] - lower_peak_props['left_ips'][0])
        right_penum = np.abs(upper_peak_props['right_ips'][0] - lower_peak_props['right_ips'][0])
        return left_penum, right_penum

    @argue.bounds(field_width=(0, 1))
    def field_values(self, field_width: float=0.8) -> np.ndarray:
        """Return a subarray of the values of the profile for the given field width.
        This is helpful for doing, e.g., flatness or symmetry calculations, where you
        want to calculate something over the field, not the whole profile.

        Parameters
        ----------
        field_width : float
            The field width of the profile, based on the fwhm. Must be between 0 and 1.

        Returns
        -------
        ndarray
        """
        left, right = self.field_edges(field_width)
        field_values = self.values[left:right]
        return field_values

    @argue.bounds(field_width=(0, 1))
    def field_edges(self, field_width: float=0.8, interpolate: bool=False) -> Tuple[NumberLike, NumberLike]:
        """Return the indices of the field width edges, based on the FWHM.

        See Also
        --------
        field_values() : Further parameter info.

        Returns
        -------
        left_index, right_index
        """
        fwhmc, _ = self.fwxm_center(interpolate=interpolate)
        field_width *= self.fwxm()
        if interpolate:
            left = fwhmc - (field_width / 2)
            right = fwhmc + (field_width / 2)
        else:
            left = int(round(fwhmc - field_width / 2))
            right = int(round(fwhmc + field_width / 2))
        return left, right

    @argue.options(calculation=('mean', 'median', 'max', 'min', 'area'))
    def field_calculation(self, field_width: float=0.8, calculation: str='mean') -> Union[float, Tuple[float, float]]:
        """Perform an operation on the field values of the profile.
        This function is useful for determining field symmetry and flatness.

        Parameters
        ----------
        calculation : {'mean', 'median', 'max', 'min', 'area'}
            Calculation to perform on the field values.

        Returns
        -------
        float

        See Also
        --------
        field_values() : Further parameter info.
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
            cax, _ = self.fwxm_center()
            lt_area = field_values[:cax+1]
            rt_area = field_values[cax:]
            return lt_area, rt_area

    def plot(self, x: int=50) -> None:
        """Plot the profile."""
        peak_idx, peak_props = find_peaks(self.values, fwxm_height=x/100, max_number=1)
        plt.plot(self.values)
        plt.plot(peak_idx, peak_props['peak_heights'][0], marker=7, color="green")
        plt.vlines(peak_idx, ymin=peak_props['peak_heights'][0]-peak_props['prominences'][0], ymax=peak_props['peak_heights'][0], color="red")
        plt.hlines(peak_props['width_heights'][0], xmin=peak_props['left_ips'], xmax=peak_props['right_ips'], color="red")
        plt.show()


class MultiProfile(ProfileMixin):
    """A class for analyzing 1-D profiles that contain multiple signals. Methods are mostly for *finding & filtering*
    the signals, peaks, valleys, etc. Profiles with a single peak (e.g. radiation beam profiles) are better suited by the SingleProfile class.

    Attributes
    ----------
    values : ndarray
        The array of values passed in on instantiation.
    peaks : list
        List of Points, containing value and index information.
    valleys : list
        Same as peaks, but for valleys.

    """
    values: Union[np.ndarray, Sequence]
    peaks: List
    valleys: List

    def __init__(self, values: Union[np.ndarray, Sequence]):
        """
        Parameters
        ----------
        values : iterable
            Array of profile values.
        """
        self.values = values
        self.peaks = []
        self.valleys = []

    def plot(self, ax: Optional[plt.Axes]=None) -> None:
        """Plot the profile.

        Parameters
        ----------
        ax: plt.Axes
            An axis to plot onto. Optional.
        """
        if ax is None:
            fig, ax = plt.subplots()
        ax.plot(self.values)
        peaks_x = [peak.idx for peak in self.peaks]
        peaks_y = [peak.value for peak in self.peaks]
        ax.plot(peaks_x, peaks_y, "gv")
        valley_x = [peak.idx for peak in self.valleys]
        valley_y = [peak.value for peak in self.valleys]
        ax.plot(valley_x, valley_y, "r^")

    def find_peaks(self, threshold: Union[float, int]=0.3, min_distance: Union[float, int]=0.05, max_number: int=None,
                   search_region: Tuple=(0.0, 1.0), peak_sort='prominences') -> Tuple[np.ndarray, np.ndarray]:
        """Find the peaks of the profile using a simple maximum value search. This also sets the `peaks` attribute.

        Parameters
        ----------
        threshold : int, float
            The value the peak must be above to be considered a peak. This removes "peaks"
            that are in a low-value region.
            If passed an int, the actual value is the threshold.
            E.g. when passed 15, any peak less with a value <15 is removed.
            If passed a float, it will threshold as a percent. Must be between 0 and 1.
            E.g. when passed 0.4, any peak <40% of the maximum value will be removed.
        min_distance : int, float
            If passed an int, parameter is the number of elements apart a peak must be from neighboring peaks.
            If passed a float, must be between 0 and 1 and represents the ratio of the profile to exclude.
            E.g. if passed 0.05 with a 1000-element profile, the minimum peak width will be 0.05*1000 = 50 elements.
        max_number : int, None
            Specify up to how many peaks will be returned. E.g. if 3 is passed in and 5 peaks are found, only the 3 largest
            peaks will be returned. If None, no limit will be applied.
        search_region : tuple of ints, floats, or both
            The region within the profile to search. The tuple specifies the (left, right) edges to search.
            This allows exclusion of edges from the search. If a value is an int, it is taken as is. If a float, must
            be between 0 and 1 and is the ratio of the profile length. The left value must be less than the right.

        Returns
        -------
        indices: ndarray, values, ndarray
            The indices and values of the peaks.
        """
        peak_idxs, peak_props = find_peaks(self.values, threshold=threshold, peak_separation=min_distance, max_number=max_number,
                                           search_region=search_region, peak_sort=peak_sort)
        self.peaks = [Point(value=peak_val, idx=peak_idx) for peak_idx, peak_val in zip(peak_idxs, peak_props['peak_heights'])]

        return peak_idxs, peak_props['peak_heights']

    def find_valleys(self, threshold: Union[float, int]=0.3, min_distance: Union[float, int]=0.05,
                     max_number: int=None, search_region: Tuple=(0.0, 1.0)) -> Tuple[np.ndarray, np.ndarray]:
        """Find the valleys (minimums) of the profile using a simple minimum value search.

        Returns
        -------
        indices: ndarray, values, ndarray
            The indices and values of the valleys.

        See Also
        --------
        :meth:`~pylinac.core.profile.MultiProfile.find_peaks` : Further parameter info.
        """
        valley_idxs, valley_props = find_peaks(-self.values, threshold=threshold, peak_separation=min_distance, max_number=max_number,
                                               search_region=search_region)
        self.valleys = [Point(value=self.values[valley_idx], idx=valley_idx) for valley_idx, valley_val in zip(valley_idxs, -valley_props['peak_heights'])]

        return valley_idxs, self.values[valley_idxs]

    @argue.bounds(x=(0, 100))
    def find_fwxm_peaks(self, x: int = 50, threshold: Union[float, int]=0.3, min_distance: Union[float, int]=0.05,
                        max_number: int=None, search_region: Tuple=(0.0, 1.0)) -> Tuple[np.ndarray, np.ndarray]:
        """Find peaks using the center of the FWXM (rather than by max value).

        Parameters
        ----------
        x : int, float
            The Full-Width-X-Maximum desired. E.g. 0.7 will return the FW70%M.
            Values must be between 0 and 100.

        See Also
        --------
        find_peaks : Further parameter info
        """
        _, peak_props = find_peaks(self.values, threshold=threshold, min_width=min_distance, max_number=max_number,
                                           search_region=search_region)
        fwxm_peak_idxs = []
        for lt, rt in zip(peak_props['left_ips'], peak_props['right_ips']):
            fwxm = int(round(lt + (rt - lt)/2))
            fwxm_peak_idxs.append(fwxm)

        fwxm_peak_vals = [self.values[fwxm] for fwxm in fwxm_peak_idxs]
        self.peaks = [Point(value=peak_val, idx=peak_idx) for peak_idx, peak_val in zip(fwxm_peak_idxs, fwxm_peak_vals)]

        return np.array(fwxm_peak_idxs), np.array(fwxm_peak_vals)


class CircleProfile(MultiProfile, Circle):
    """A profile in the shape of a circle.

    Attributes
    ----------
    image_array : ndarray
        The 2D image array.
    start_angle : int, float
        Starting position of the profile in radians; 0 is right (0 on unit circle).
    ccw : bool
        How the profile is/was taken; clockwise or counter-clockwise.
    """
    image_array: np.ndarray
    start_angle: Union[float, int]
    ccw: bool
    sampling_ratio: float
    _x_locations: Optional[np.ndarray]
    _y_locations: Optional[np.ndarray]

    def __init__(self, center: Point, radius: NumberLike, image_array: np.ndarray,
                 start_angle: Union[float, int]=0, ccw: bool=True, sampling_ratio: float=1.0):
        """
        Parameters
        ----------
        image_array : ndarray
            The 2D image array.
        start_angle : int, float
            Starting position of the profile in radians; 0 is right (0 on unit circle).
        ccw : bool
            If True (default), the profile will proceed counter-clockwise (the direction on the unit circle).
            If False, will proceed clockwise.
        sampling_ratio : float
            The ratio of pixel sampling to real pixels. E.g. if 1.0, the profile will have approximately
            the same number of elements as was encountered in the profile. A value of 2.0 will sample
            the profile at 2x the number of elements.

        See Also
        --------
        :class:`~pylinac.core.geometry.Circle` : Further parameter info.
        """
        Circle.__init__(self, center, radius)
        self._ensure_array_size(image_array, self.radius + self.center.x, self.radius + self.center.y)
        self.image_array = image_array
        self.start_angle = start_angle
        self.ccw = ccw
        self.sampling_ratio = sampling_ratio
        self._x_locations = None
        self._y_locations = None
        MultiProfile.__init__(self, self._profile)

    @property
    def size(self) -> float:
        """The elemental size of the profile."""
        return np.pi * self.radius * 2 * self.sampling_ratio

    @property
    def _radians(self) -> np.ndarray:
        interval = (2 * np.pi) / self.size
        rads = np.arange(0 + self.start_angle, (2 * np.pi) + self.start_angle - interval, interval)
        if self.ccw:
            rads = rads[::-1]
        return rads

    @property
    def x_locations(self) -> np.ndarray:
        """The x-locations of the profile values."""
        if self._x_locations is None:
            return np.cos(self._radians) * self.radius + self.center.x
        else:
            return self._x_locations

    @x_locations.setter
    def x_locations(self, array: np.ndarray):
        self._x_locations = array

    @property
    def y_locations(self) -> np.ndarray:
        """The x-locations of the profile values."""
        if self._y_locations is None:
            return np.sin(self._radians) * self.radius + self.center.y
        else:
            return self._y_locations

    @y_locations.setter
    def y_locations(self, array: np.ndarray):
        self._y_locations = array

    @property
    def _profile(self) -> np.ndarray:
        """The actual profile array; private attr that is passed to MultiProfile."""
        return ndimage.map_coordinates(self.image_array, [self.y_locations, self.x_locations], order=0)

    def find_peaks(self, threshold: Union[float, int]=0.3, min_distance: Union[float, int]=0.05,
                   max_number: int=None, search_region: Tuple[float, float]=(0.0, 1.0)) -> Tuple[np.ndarray, np.ndarray]:
        """Overloads Profile to also map peak locations to the image."""
        peak_idxs, peak_vals = super().find_peaks(threshold, min_distance, max_number, search_region)
        self._map_peaks()
        return peak_idxs, peak_vals

    def find_valleys(self, threshold: Union[float, int]=0.3, min_distance: Union[float, int]=0.05,
                     max_number: int=None, search_region: Tuple[float, float]=(0.0, 1.0)) -> Tuple[np.ndarray, np.ndarray]:
        """Overload Profile to also map valley locations to the image."""
        valley_idxs, valley_vals = super().find_valleys(threshold, min_distance, max_number, search_region)
        self._map_peaks()
        return valley_idxs, valley_vals

    @argue.bounds(x=(0, 100))
    def find_fwxm_peaks(self, x: int=50, threshold: Union[float, int]=0.3, min_distance: Union[float, int]=0.05,
                        max_number: int=None, search_region: Tuple[float, float]=(0.0, 1.0)) -> Tuple[np.ndarray, np.ndarray]:
        """Overloads Profile to also map the peak locations to the image."""
        peak_idxs, peak_vals = super().find_fwxm_peaks(x, threshold, min_distance, max_number,
                                                       search_region=search_region)
        self._map_peaks()
        return peak_idxs, peak_vals

    def _map_peaks(self) -> None:
        """Map found peaks to the x,y locations on the image/array; i.e. adds x,y coordinates to the peak locations"""
        for peak in self.peaks:
            peak.x = self.x_locations[int(peak.idx)]
            peak.y = self.y_locations[int(peak.idx)]

    def roll(self, amount: int) -> None:
        """Roll the profile and x and y coordinates."""
        self.values = np.roll(self.values, -amount)
        self.x_locations = np.roll(self.x_locations, -amount)
        self.y_locations = np.roll(self.y_locations, -amount)

    def plot2axes(self, axes: plt.Axes=None, edgecolor: str='black', fill: bool=False, plot_peaks: bool=True) -> None:
        """Plot the circle to an axes.

        Parameters
        ----------
        axes : matplotlib.Axes, None
            The axes to plot on. If None, will create a new figure of the image array.
        edgecolor : str
            Color of the Circle; must be a valid matplotlib color.
        fill : bool
            Whether to fill the circle. matplotlib keyword.
        plot_peaks : bool
            If True, plots the found peaks as well.
        """
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
    def _ensure_array_size(array: np.ndarray, min_width: int, min_height: int) -> None:
            """Ensure the array size of inputs are greater than the minimums."""
            height = array.shape[0]
            width = array.shape[1]
            if width < min_width or height < min_height:
                raise ValueError("Array size not large enough to compute profile")


class CollapsedCircleProfile(CircleProfile):
    """A circular profile that samples a thick band around the nominal circle, rather than just a 1-pixel-wide profile
    to give a mean value.
    """
    width_ratio: float
    num_profiles: int

    @argue.bounds(width_ratio=(0, 1))
    def __init__(self, center: Point, radius: NumberLike, image_array: np.ndarray, start_angle: int=0,
                 ccw: bool=True, sampling_ratio: float=1.0, width_ratio: float=0.1, num_profiles: int=20):
        """
        Parameters
        ----------
        width_ratio : float
            The "thickness" of the band to sample. The ratio is relative to the radius. E.g. if the radius is 20
            and the width_ratio is 0.2, the "thickness" will be 4 pixels.
        num_profiles : int
            The number of profiles to sample in the band. Profiles are distributed evenly within the band.

        See Also
        --------
        :class:`~pylinac.core.profile.CircleProfile` : Further parameter info.
        """
        self.width_ratio = width_ratio
        self.num_profiles = num_profiles
        super().__init__(center, radius, image_array, start_angle, ccw, sampling_ratio)

    @property
    def _radii(self) -> np.ndarray:
        return np.linspace(start=self.radius * (1 - self.width_ratio), stop=self.radius * (1 + self.width_ratio),
                  num=self.num_profiles)

    @property
    def size(self) -> float:
        return np.pi * max(self._radii) * 2 * self.sampling_ratio

    @property
    def _multi_x_locations(self) -> List:
        """List of x-locations of the sampling profiles"""
        x = []
        cos = np.cos(self._radians)
        # extract profile for each circle radii
        for radius in self._radii:
            x.append(cos * radius + self.center.x)
        return x

    @property
    def _multi_y_locations(self) -> List:
        """List of x-locations of the sampling profiles"""
        y = []
        sin = np.sin(self._radians)
        # extract profile for each circle radii
        for radius in self._radii:
            y.append(sin * radius + self.center.y)
        return y

    @property
    def _profile(self) -> np.ndarray:
        """The actual profile array; private attr that is passed to MultiProfile."""
        profile = np.zeros(len(self._multi_x_locations[0]))
        for radius, x, y in zip(self._radii, self._multi_x_locations, self._multi_y_locations):
            profile += ndimage.map_coordinates(self.image_array, [y, x], order=0)
        profile /= self.num_profiles
        return profile

    def plot2axes(self, axes: plt.Axes=None, edgecolor: str='black', fill: bool=False, plot_peaks: bool=True) -> None:
        """Add 2 circles to the axes: one at the maximum and minimum radius of the ROI.

        See Also
        --------
        :meth:`~pylinac.core.profile.CircleProfile.plot2axes` : Further parameter info.
        """
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
            axes.scatter(x_locs, y_locs, s=20, marker='x', c=edgecolor)


def find_peaks(values: np.ndarray, threshold: Union[float, int] = -np.inf, peak_separation: Union[float, int] = 0,
               max_number: int = None, fwxm_height: float = 0.5, min_width: int = 0,
               search_region: Tuple[float, float] = (0.0, 1.0), peak_sort='prominences') \
        -> Tuple[np.ndarray, dict]:
    """Find the peaks of a 1D signal. Heavily relies on the scipy implementation.

    Parameters
    ----------
    values : array-like
        Signal values to search for peaks within.
    threshold : int, float
        The value the peak must be above to be considered a peak. This removes "peaks"
        that are in a low-value region.
        If passed an int, the actual value is the threshold.
        E.g. when passed 15, any peak less with a value <15 is removed.
        If passed a float, it will threshold as a percent. Must be between 0 and 1.
        E.g. when passed 0.4, any peak <40% of the maximum value will be removed.
    peak_separation : int, float
        If passed an int, parameter is the number of elements apart a peak must be from neighboring peaks.
        If passed a float, must be between 0 and 1 and represents the ratio of the profile to exclude.
        E.g. if passed 0.05 with a 1000-element profile, the minimum peak width will be 0.05*1000 = 50 elements.
    max_number : int, None
        Specify up to how many peaks will be returned. E.g. if 3 is passed in and 5 peaks are found, only the 3 largest
        peaks will be returned.
    fwxm_height: float
        The relative height at which a FWXM calculation is performed. Although this function finds simple max values,
        the underlying function can provide fwxm information as well.
    min_width: int
        The minimum width of the peak.
    search_region: tuple
        The search region to use within the values.
        Using between 0 and 1 will convert to a ratio of the indices. E.g. to search the middle half of the passed values, use (0.25, 0.75).
        Using ints above 1 will use the indices directly. E.g. (33, 71) will search between those two indices.
    peak_sort
        Either 'peak_heights' or 'prominences'. This is the method for determining the peaks. Usually not needed
        unless the wrong number of pickets have been detected.
        See the scipy.signal.find_peaks function for more information.

    Returns
    -------
    peak_idxs : numpy.array
        The indices of the peaks found.
    peak_props : dict
        A dict containing contextual peak data.
    """
    peak_separation, shift_amount, threshold, trimmed_values = _parse_peak_args(peak_separation, search_region, threshold,
                                                                                values)

    peak_idxs, peak_props = signal.find_peaks(trimmed_values, rel_height=(1 - fwxm_height), width=min_width, height=threshold,
                                              distance=peak_separation)
    peak_idxs += shift_amount  # shift according to the search region left edge

    # get the "largest" peaks up to max number, and then re-sort to be left->right like it was originally
    largest_peak_idxs = sorted(list(np.argsort(peak_props[peak_sort]))[::-1][:max_number])

    # cut down prop arrays as need be
    for key, array_vals in peak_props.items():
        peak_props[key] = array_vals[largest_peak_idxs]
    return peak_idxs[largest_peak_idxs], peak_props


def _parse_peak_args(peak_separation: NumberLike, search_region: Tuple[float, float], threshold: NumberLike,
                     values: np.ndarray) -> Tuple[NumberLike, int, NumberLike, np.ndarray]:
    """Converts arguments as needed. E.g. converting a ratio to actual values"""
    # set threshold as % if between 0 and 1
    val_range = values.max() - values.min()
    if 0 <= threshold <= 1:
        threshold = values.min() + threshold * val_range
    # set separation as % if between 0 and 1
    if 0 <= peak_separation <= 1:
        peak_separation = max(int(peak_separation * len(values)), 1)
    # limit to search region
    if max(search_region) <= 1:
        shift_amount = int(search_region[0] * len(values))
        values = values[int(search_region[0] * len(values)):int(search_region[1] * len(values))]
    else:
        values = values[search_region[0]:search_region[1]]
        shift_amount = search_region[0]
    return peak_separation, shift_amount, threshold, values


