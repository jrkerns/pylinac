"""Module of objects that resemble or contain a profile, i.e. a 1 or 2-D f(x) representation."""
import enum
import warnings
from typing import Union, Tuple, Sequence, List, Optional

import argue
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle as mpl_Circle
from scipy import ndimage, signal
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import OptimizeWarning, minimize
from scipy.stats import linregress

from .geometry import Point, Circle
from .hill import Hill
from .typing import NumberLike

# for Hill fits of 2D device data the # of points can be small.
# This results in optimization warnings about the variance of the fit (the variance isn't of concern for us for that particular item)
warnings.simplefilter("ignore", OptimizeWarning)


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


class Interpolation(enum.Enum):
    """Interpolation Enum"""
    NONE = None
    LINEAR = 'Linear'
    SPLINE = 'Spline'


class Normalization(enum.Enum):
    """Normalization method Enum"""
    NONE = None
    GEOMETRIC_CENTER = 'Geometric center'
    BEAM_CENTER = 'Beam center'
    MAX = 'Max'


class Edge(enum.Enum):
    """Edge detection Enum"""
    FWHM = 'FWHM'
    INFLECTION_DERIVATIVE = 'Inflection Derivative'
    INFLECTION_HILL = 'Inflection Hill'


class SingleProfile(ProfileMixin):
    """A profile that has one large signal, e.g. a radiation beam profile.
    Signal analysis methods are given, mostly based on FWXM and on Hill function calculations.
    Profiles with multiple peaks are better suited by the MultiProfile class.
    """

    def __init__(self, values: np.ndarray, dpmm: float = None,
                 interpolation: Interpolation = Interpolation.LINEAR,
                 ground: bool = True,
                 interpolation_resolution_mm: float = 0.1,
                 interpolation_factor: float = 10,
                 normalization_method: Normalization = Normalization.BEAM_CENTER,
                 edge_detection_method: Edge = Edge.FWHM,
                 edge_smoothing_ratio: float = 0.003,
                 hill_window_ratio: float = 0.1):
        """
        Parameters
        ----------
        values
            The profile numpy array. Must be 1D.
        dpmm
            The dots (pixels) per mm. Pass to get info like beam width in distance units in addition to pixels
        interpolation
            Interpolation technique.
        ground
            Whether to ground the profile (set min value to 0). Helpful most of the time.
        interpolation_resolution_mm
            The resolution that the interpolation will scale to. *Only used if dpmm is passed and interpolation is set*.
            E.g. if the dpmm is 0.5 and the resolution is set to 0.1mm the data will be interpolated to have a new dpmm of 10 (1/0.1).
        interpolation_factor
            The factor to multiply the data by. **Only used if interpolation is used and dpmm is NOT passed**. E.g. 10
            will perfectly decimate the existing data according to the interpolation method passed.
        normalization_method
            How to pick the point to normalize the data to.
        edge_detection_method
            The method by which to detect the field edge. FWHM is reasonable most of the time except for FFF beams.
            Inflection-derivative will use the max gradient to determine the field edge. Note that this may not be the
            50% height. In fact, for FFF beams it shouldn't be. Inflection methods are better for FFF and other unusual
            beam shapes.
        edge_smoothing_ratio
            *Only applies to INFLECTION_DERIVATIVE and INFLECTION_HILL.*

            The ratio of the length of the values to use as the sigma for a Gaussian filter applied before searching for
            the inflection. E.g. 0.005 with a profile of 1000 points will result in a sigma of 5.
            This helps make the inflection point detection more robust to noise. Increase for noisy data.
        hill_window_ratio
            The ratio of the field size to use as the window to fit the Hill function. E.g. 0.2 will using a window
            centered about each edge with a width of 20% the size of the field width. *Only applies when the edge
            detection is ``INFLECTION_HILL`` *.
        """
        self._interp_method = interpolation
        self._interpolation_res = interpolation_resolution_mm
        self._interpolation_factor = interpolation_factor
        self._norm_method = normalization_method
        self._edge_method = edge_detection_method
        self._edge_smoothing_ratio = edge_smoothing_ratio
        self._hill_window_ratio = hill_window_ratio
        self.values = values  # set initial data so we can do things like find beam center
        self.dpmm = dpmm
        fitted_values, new_dpmm, x_indices = self._interpolate(values, dpmm, interpolation_resolution_mm,
                                                               interpolation_factor, interpolation)
        self.dpmm = new_dpmm  # update as needed
        self.values = fitted_values
        self.x_indices = x_indices
        if ground:
            fitted_values -= fitted_values.min()
        norm_values = self._normalize(fitted_values, normalization_method)
        self.values = norm_values  # update values

    @staticmethod
    def _interpolate(values, dpmm, interpolation_resolution, interpolation_factor, interp_method: Interpolation) -> (
            np.ndarray, float, float, float):
        """Fit the data to the passed interpolation method. Will also calculate the new values to correct the measurements such as dpmm"""
        x_indices = list(range(len(values)))
        if interp_method == Interpolation.NONE:
            return values, dpmm, x_indices  # do nothing
        elif interp_method == Interpolation.LINEAR:
            if dpmm is not None:
                samples = int(round(len(x_indices)/(dpmm*interpolation_resolution)))
                new_dpmm = 1/interpolation_resolution
            else:
                samples = int(round(len(x_indices)*interpolation_factor))
                new_dpmm = None
            f = interp1d(x_indices, values, kind='linear', bounds_error=False)
            new_x = np.linspace(0, len(x_indices)-1, num=samples)
            return f(new_x), new_dpmm, new_x
        elif interp_method == Interpolation.SPLINE:
            if dpmm is not None:
                samples = int(round(len(x_indices)/(dpmm*interpolation_resolution)))
                new_dpmm = 1 / interpolation_resolution
            else:
                samples = int(round(len(x_indices)*interpolation_factor))
                new_dpmm = None
            f = interp1d(x_indices, values, kind='cubic')
            new_x = np.linspace(0, len(x_indices)-1, num=samples)
            return f(new_x), new_dpmm, new_x

    def _normalize(self, values, method: Normalization) -> np.ndarray:
        """Normalize the data given a method."""
        if method == Normalization.NONE:
            return values
        elif method == Normalization.MAX:
            return values / values.max()
        elif method == Normalization.GEOMETRIC_CENTER:
            return values / self._geometric_center(values)['value (exact)']
        elif method == Normalization.BEAM_CENTER:
            return values / self.beam_center()['value (@rounded)']

    def _geometric_center(self, values) -> dict:
        """Returns the center index and value of the profile.

         If the profile has an even number of values the centre lies between the two centre indices and the centre
         value is the average of the two centre values else the centre index and value are returned."""
        plen = values.shape[0]
        # buffer overflow can cause the below addition to give strange results
        values = values.astype(np.float64)
        if plen % 2 == 0:  # plen is even and central detectors straddle CAX
            cax = (values[int(plen / 2)] + values[int(plen / 2) - 1]) / 2.0
        else:  # plen is odd and we have a central detector
            cax = values[int((plen - 1) / 2)]
        plen = (plen - 1)/2.0
        return {'index (exact)': plen, 'value (exact)': cax}

    def geometric_center(self) -> dict:
        """The geometric center (i.e. the device center)"""
        return self._geometric_center(self.values)

    def beam_center(self) -> dict:
        """The center of the detected beam. This can account for asymmetries in the beam position (e.g. offset jaws)"""
        if self._edge_method == Edge.FWHM:
            data = self.fwxm_data(x=50)
            return {'index (rounded)': data['center index (rounded)'],
                    'index (exact)': data['center index (exact)'],
                    'value (@rounded)': data['center value (@rounded)']}
        elif self._edge_method in (Edge.INFLECTION_DERIVATIVE, Edge.INFLECTION_HILL):
            infl = self.inflection_data()
            mid_point = infl['left index (exact)'] + (infl['right index (exact)'] - infl['left index (exact)']) / 2
            return {'index (rounded)': int(round(mid_point)),
                    'index (exact)': mid_point,
                    'value (@rounded)': self.values[int(round(mid_point))]}

    @argue.bounds(x=(0, 100))
    def fwxm_data(self, x: int = 50) -> dict:
        """Return the width at X-Max, where X is the percentage height.

        Parameters
        ----------
        x
            The percent height of the profile. E.g. x = 50 is 50% height,
            i.e. FWHM.
        """
        _, peak_props = find_peaks(self.values, fwxm_height=x/100, max_number=1)
        left_idx = peak_props['left_ips'][0]
        right_idx = peak_props['right_ips'][0]
        fwxm_center_idx = ((peak_props['right_ips'][0] - peak_props['left_ips'][0]) / 2 + peak_props['left_ips'][0])

        data = {'width (exact)': peak_props['widths'][0],
                'width (rounded)': int(round(right_idx)) - int(round(left_idx)),
                'center index (rounded)': int(round(fwxm_center_idx)),
                'center index (exact)': fwxm_center_idx,
                'center value (@rounded)': self.values[int(round(fwxm_center_idx))],
                'left index (exact)': left_idx,
                'left index (rounded)': int(round(left_idx)),
                'left value (@rounded)': self.values[int(round(left_idx))],
                'right index (exact)': right_idx,
                'right index (rounded)': int(round(right_idx)),
                'right value (@rounded)': self.values[int(round(right_idx))],
                'field values': self.values[int(round(left_idx)):
                                            int(round(right_idx))]}
        if self.dpmm:
            data['width (exact) mm'] = data['width (exact)'] / self.dpmm
            data['left distance (exact) mm'] = abs(data['center index (exact)'] - data['left index (exact)']) / self.dpmm
            data['right distance (exact) mm'] = abs(data['right index (exact)'] - data['center index (exact)']) / self.dpmm

        return data

    @argue.bounds(in_field_ratio=(0, 1.0), slope_exclusion_ratio=(0, 1.0))
    def field_data(self, in_field_ratio: float = 0.8, slope_exclusion_ratio=0.2) -> dict:
        """Return the width at X-Max, where X is the percentage height.

        Parameters
        ----------
        in_field_ratio
            In Field Ratio: 1.0 is the entire detected field; 0.8 would be the central 80%, etc.
        slope_exclusion_ratio
            Ratio of the field width to use as the cutoff between "top" calculation and "slope" calculation. Useful for FFF beams.
            This area is centrally located in the field. E.g. 0.2 will use the central 20% of the field to calculate
            the "top" value. To calculate the slope of each side, the field width between the edges of the in_field_ratio
            and the slope exclusion region are used.

            .. warning:: The "top" value is always calculated. For FFF beams this should be reasonable, but for flat beams
                         this value may end up being non-sensible.
        """
        if slope_exclusion_ratio >= in_field_ratio:
            raise ValueError("The exclusion region must be smaller than the field ratio")
        if self._edge_method == Edge.FWHM:
            data = self.fwxm_data(x=50)
            beam_center_idx = data['center index (exact)']
            full_width = data['width (exact)']

        elif self._edge_method in (Edge.INFLECTION_DERIVATIVE, Edge.INFLECTION_HILL):
            infl_data = self.inflection_data()
            beam_center_idx = self.beam_center()['index (exact)']
            full_width = infl_data['right index (exact)'] - infl_data['left index (exact)']
        beam_center_idx_r = int(round(beam_center_idx))

        cax_idx = self.geometric_center()['index (exact)']
        cax_idx_r = int(round(cax_idx))

        field_left_idx = beam_center_idx - in_field_ratio * full_width / 2
        field_left_idx_r = int(round(field_left_idx))
        field_right_idx = beam_center_idx + in_field_ratio * full_width / 2
        field_right_idx_r = int(round(field_right_idx))
        field_width = field_right_idx - field_left_idx

        # slope calcs
        inner_left_idx = beam_center_idx - slope_exclusion_ratio*field_width/2
        inner_left_idx_r = int(round(inner_left_idx))
        inner_right_idx = beam_center_idx + slope_exclusion_ratio*field_width/2
        inner_right_idx_r = int(round(inner_right_idx))
        left_fit = linregress(range(field_left_idx_r, inner_left_idx_r),
                              self.values[field_left_idx_r:inner_left_idx_r])
        right_fit = linregress(range(inner_right_idx_r, field_right_idx_r),
                               self.values[inner_right_idx_r:field_right_idx_r])

        # top calc
        fit_params = np.polyfit(range(inner_left_idx_r, inner_right_idx_r),
                                self.values[inner_left_idx_r:inner_right_idx_r], deg=2)
        width = abs(inner_right_idx_r - inner_left_idx_r)

        def poly_func(x):
            # return the negative since we're MINIMIZING and want the top value
            return -(fit_params[0] * (x ** 2) + fit_params[1] * x + fit_params[2])

        # minimize the polynomial function
        min_f = minimize(poly_func, x0=(inner_left_idx_r+width/2,), bounds=((inner_left_idx_r, inner_right_idx_r),))
        top_idx = min_f.x[0]
        top_val = -min_f.fun

        data = {'width (exact)': field_width,
                'beam center index (exact)': beam_center_idx,
                'beam center index (rounded)': beam_center_idx_r,
                'beam center value (@rounded)': self.values[int(round(beam_center_idx))],
                'cax index (exact)': cax_idx,
                'cax index (rounded)': cax_idx_r,
                'cax value (@rounded)': self.values[int(round(cax_idx))],
                'left index (exact)': field_left_idx,
                'left index (rounded)': field_left_idx_r,
                'left value (@rounded)': self.values[int(round(field_left_idx))],
                'left slope': left_fit.slope,
                'left intercept': left_fit.intercept,
                'right slope': right_fit.slope,
                'right intercept': right_fit.intercept,
                'left inner index (exact)': inner_left_idx,
                'left inner index (rounded)': inner_left_idx_r,
                'right inner index (exact)': inner_right_idx,
                'right inner index (rounded)': inner_right_idx_r,
                '"top" index (exact)': top_idx,
                '"top" index (rounded)': int(round(top_idx)),
                '"top" value (@exact)': top_val,
                'top params': fit_params,
                'right index (exact)': field_right_idx,
                'right index (rounded)': field_right_idx_r,
                'right value (@rounded)': self.values[int(round(field_right_idx))],
                'field values': self.values[int(round(field_left_idx)):
                                            int(round(field_right_idx))]}
        if self.dpmm:
            data['width (exact) mm'] = data['width (exact)'] / self.dpmm
            data['left slope (%/mm)'] = data['left slope'] * self.dpmm * 100
            data['right slope (%/mm)'] = data['right slope'] * self.dpmm * 100
            data['left distance->beam center (exact) mm'] = abs(data['beam center index (exact)'] - data['left index (exact)']) / self.dpmm
            data['right distance->beam center (exact) mm'] = abs(data['right index (exact)'] - data['beam center index (exact)']) / self.dpmm
            data['left distance->CAX (exact) mm'] = abs(data['cax index (exact)'] - data['left index (exact)']) / self.dpmm
            data['right distance->CAX (exact) mm'] = abs(data['cax index (exact)'] - data['right index (exact)']) / self.dpmm

            data['left distance->top (exact) mm'] = abs(data['"top" index (exact)'] - data['left index (exact)']) / self.dpmm
            data['right distance->top (exact) mm'] = abs(data['"top" index (exact)'] - data['right index (exact)']) / self.dpmm
            data['"top"->beam center (exact) mm'] = (data['"top" index (exact)'] - data['beam center index (exact)']) / self.dpmm
            data['"top"->CAX (exact) mm'] = abs(data['"top" index (exact)'] - data['cax index (exact)']) / self.dpmm
        return data

    def inflection_data(self) -> dict:
        """Calculate the profile inflection values using either the 2nd derivative or a fitted Hill function.

        .. note::
            This only applies if the edge detection method is `INFLECTION_...`.

        Parameters
        ----------

        """
        # get max/min of the gradient, which is basically the same as the 2nd deriv 0-crossing
        if self._edge_method == Edge.FWHM:
            raise ValueError("FWHM edge method does not have inflection points. Use a different edge detection method")
        d1 = np.gradient(gaussian_filter1d(self.values, sigma=self._edge_smoothing_ratio * len(self.values)))
        left_idx = np.argmax(d1)
        right_idx = np.argmin(d1)
        if self._edge_method == Edge.INFLECTION_DERIVATIVE:
            data = {'left index (rounded)': left_idx,
                    'left index (exact)': left_idx,
                    'right index (rounded)': right_idx,
                    'right index (exact)': right_idx,
                    'left value (@rounded)': self.values[int(round(left_idx))],
                    'left value (@exact)': self.values[int(round(left_idx))],
                    'right value (@rounded)': self.values[int(round(right_idx))],
                    'right value (@exact)': self.values[int(round(right_idx))]
                    }
            return data
        else:  # Hill
            # the 2nd deriv is a good approximation for the inflection point. Start there and fit Hill about it
            # penum_half_window = self.field_data()['width (exact)'] * self._hill_window_ratio / 2
            penum_half_window = int(round(self._hill_window_ratio * abs(right_idx - left_idx) / 2))

            # left side
            x_data = np.array([x for x in np.arange(left_idx - penum_half_window, left_idx + penum_half_window) if x >=0])
            y_data = self.values[x_data]
            # y_data = self.values[left_idx - penum_half_window: left_idx + penum_half_window]
            left_hill = Hill.fit(x_data, y_data)
            left_infl = left_hill.inflection_idx()

            # right side
            x_data = np.array([x for x in np.arange(right_idx - penum_half_window, right_idx + penum_half_window) if x < len(d1)])
            y_data = self.values[x_data]
            right_hill = Hill.fit(x_data, y_data)
            right_infl = right_hill.inflection_idx()

            data = {'left index (rounded)': left_infl['index (rounded)'],
                    'left index (exact)': left_infl['index (exact)'],
                    'right index (rounded)': right_infl['index (rounded)'],
                    'right index (exact)': right_infl['index (exact)'],
                    'left value (@exact)': left_hill.y(left_infl['index (exact)']),
                    'right value (@exact)': right_hill.y(right_infl['index (exact)']),
                    'left Hill params': left_hill.params, 'right Hill params': right_hill.params}
            return data

    def penumbra(self, lower: int = 20, upper: int = 80):
        """Calculate the penumbra of the field. Dependent on the edge detection method.

        Parameters
        ----------
        lower
            The lower % of the beam to use. If the edge method is FWHM, this is the typical % penumbra you're thinking.
            If the inflection method is used it will be the value/50 of the inflection point value. E.g. if the inflection
            point is perfectly at 50% with a ``lower`` of 20, then the penumbra value here will be 20% of the maximum.
            If the inflection point is at 30% of the max value (say for a FFF beam) then the lower penumbra will be ``lower/50``
            of the inflection point or ``0.3*lower/50``.
        upper
            Upper % of the beam to use. See lower for details.
        """
        if lower > upper:
            raise ValueError("Upper penumbra value must be larger than the lower penumbra value")
        if self._edge_method == Edge.FWHM:
            upper_data = self.fwxm_data(x=upper)
            lower_data = self.fwxm_data(x=lower)
            data = {f'left {lower}% index (exact)': lower_data['left index (exact)'],
                    f'left {lower}% value (@rounded)': lower_data['left value (@rounded)'],
                    f'left {upper}% index (exact)': upper_data['left index (exact)'],
                    f'left {upper}% value (@rounded)': upper_data['left value (@rounded)'],

                    f'right {lower}% index (exact)': lower_data['right index (exact)'],
                    f'right {lower}% value (@rounded)': lower_data['right value (@rounded)'],
                    f'right {upper}% index (exact)': upper_data['right index (exact)'],
                    f'right {upper}% value (@rounded)': upper_data['right value (@rounded)'],

                    'left values': self.values[lower_data['left index (rounded)']:upper_data['left index (rounded)']],
                    'right values': self.values[upper_data['right index (rounded)']:lower_data['right index (rounded)']],

                    f'left penumbra width (exact)': abs(upper_data['left index (exact)'] - lower_data['left index (exact)']),
                    f'right penumbra width (exact)': abs(upper_data['right index (exact)'] - lower_data['right index (exact)']),
                    }
            if self.dpmm:
                data['left penumbra width (exact) mm'] = data['left penumbra width (exact)'] / self.dpmm
                data['right penumbra width (exact) mm'] = data['right penumbra width (exact)'] / self.dpmm
            return data
        elif self._edge_method == Edge.INFLECTION_DERIVATIVE:
            infl_data = self.inflection_data()
            lower_left_value = infl_data['left value (@rounded)']*lower/50*100
            upper_left_value = infl_data['left value (@rounded)']*upper/50*100
            upper_left_data = self.fwxm_data(x=upper_left_value)
            lower_left_data = self.fwxm_data(x=lower_left_value)

            lower_right_value = infl_data['right value (@exact)']*lower/50*100
            upper_right_value = infl_data['right value (@exact)']*upper/50*100
            upper_right_data = self.fwxm_data(x=upper_right_value)
            lower_right_data = self.fwxm_data(x=lower_right_value)

            data = {f'left {lower}% index (exact)': lower_left_data['left index (exact)'],
                    f'left {upper}% index (exact)': upper_left_data['left index (exact)'],

                    f'right {lower}% index (exact)': lower_right_data['right index (exact)'],
                    f'right {upper}% index (exact)': upper_right_data['right index (exact)'],

                    'left values': self.values[lower_left_data['left index (rounded)']:upper_left_data['left index (rounded)']],
                    'right values': self.values[upper_right_data['right index (rounded)']:lower_right_data['right index (rounded)']],

                    f'left penumbra width (exact)': abs(upper_left_data['left index (exact)'] - lower_left_data['left index (exact)']),
                    f'right penumbra width (exact)': abs(upper_right_data['right index (exact)'] - lower_right_data['right index (exact)']),
                    }
            if self.dpmm:
                data['left penumbra width (exact) mm'] = data['left penumbra width (exact)'] / self.dpmm
                data['right penumbra width (exact) mm'] = data['right penumbra width (exact)'] / self.dpmm
            return data
        elif self._edge_method == Edge.INFLECTION_HILL:
            infl_data = self.inflection_data()
            left_hill = Hill.from_params(infl_data['left Hill params'])
            right_hill = Hill.from_params(infl_data['right Hill params'])

            lower_left_value = infl_data['left value (@exact)']*lower/50
            lower_left_index = left_hill.x(lower_left_value)
            upper_left_value = infl_data['left value (@exact)']*upper/50
            upper_left_index = left_hill.x(upper_left_value)

            lower_right_value = infl_data['right value (@exact)']*lower/50
            lower_right_index = right_hill.x(lower_right_value)
            upper_right_value = infl_data['right value (@exact)']*upper/50
            upper_right_index = right_hill.x(upper_right_value)

            data = {f'left {lower}% index (exact)': lower_left_index,
                    f'left {lower}% value (exact)': lower_left_value,
                    f'left {upper}% index (exact)': upper_left_index,
                    f'left {upper}% value (exact)': upper_left_value,

                    f'right {lower}% index (exact)': lower_right_index,
                    f'right {lower}% value (exact)': lower_right_value,
                    f'right {upper}% index (exact)': upper_right_index,
                    f'right {upper}% value (exact)': upper_right_value,

                    'left values': self.values[int(round(lower_left_index)):int(round(upper_left_index))],
                    'right values': self.values[int(round(upper_right_index)):int(round(lower_right_index))],

                    f'left penumbra width (exact)': abs(upper_left_index - lower_left_index),
                    f'right penumbra width (exact)': abs(upper_right_index - lower_right_index),

                    f'left gradient (exact)': left_hill.gradient_at(infl_data['left index (exact)']),
                    r'right gradient (exact)': right_hill.gradient_at(infl_data['right index (exact)']),
                    }
            if self.dpmm:
                data['left penumbra width (exact) mm'] = data['left penumbra width (exact)'] / self.dpmm
                data['left gradient (exact) %/mm'] = data['left gradient (exact)'] * self.dpmm * 100  # 100 to convert to %
                data['right penumbra width (exact) mm'] = data['right penumbra width (exact)'] / self.dpmm
                data['right gradient (exact) %/mm'] = data['right gradient (exact)'] * self.dpmm * 100
            return data

    @argue.options(calculation=('mean', 'median', 'max', 'min', 'area'))
    def field_calculation(self, in_field_ratio: float=0.8, calculation: str='mean') -> Union[float, Tuple[float, float]]:
        """Perform an operation on the field values of the profile.
        This function is useful for determining field symmetry and flatness.

        Parameters
        ----------
        in_field_ratio
            Ratio of the field width to use in the calculation.
        calculation : {'mean', 'median', 'max', 'min', 'area'}
            Calculation to perform on the field values.
        """
        field_values = self.field_data(in_field_ratio)

        if calculation == 'mean':
            return field_values['field values'].mean()
        elif calculation == 'median':
            return np.median(field_values['field values'])
        elif calculation == 'max':
            return field_values['field values'].max()
        elif calculation == 'min':
            return field_values['field values'].min()

    def plot(self) -> None:
        """Plot the profile."""
        plt.plot(self.values)
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


