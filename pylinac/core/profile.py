from __future__ import annotations  # noqa:I001

import enum
import math
import warnings
from abc import ABC, abstractmethod
from functools import cached_property
from typing import Any, Iterable, Literal, Sequence

import argue
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle as mpl_Circle
from scipy import ndimage, signal
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline, interp1d
from scipy.ndimage import gaussian_filter1d, zoom
from scipy.optimize import OptimizeWarning, minimize
from scipy.stats import linregress

# noqa:I001
from ..metrics.profile import FlatnessDifferenceMetric  # noqa:F401
from ..metrics.profile import FlatnessRatioMetric  # noqa:F401
from ..metrics.profile import PenumbraLeftMetric  # noqa:F401
from ..metrics.profile import PenumbraRightMetric  # noqa:F401
from ..metrics.profile import SymmetryAreaMetric  # noqa:F401
from ..metrics.profile import SymmetryPointDifferenceMetric  # noqa:F401
from ..metrics.profile import SymmetryPointDifferenceQuotientMetric  # noqa:F401
from ..metrics.profile import LEFT, RIGHT, ProfileMetric
from . import array_utils as utils
from . import validators
from .geometry import Circle, Point
from .hill import Hill
from .utilities import convert_to_enum

# for Hill fits of 2D device data the # of points can be small.
# This results in optimization warnings about the variance of the fit (the variance isn't of concern for us for that particular item)
warnings.simplefilter("ignore", OptimizeWarning)


def gamma_1d(
    reference: np.ndarray,
    evaluation: np.ndarray,
    dose_to_agreement: float = 1,
    distance_to_agreement: int = 1,
    gamma_cap_value: float = 2,
    global_dose: bool = True,
    dose_threshold: float = 5,
    fill_value: float = np.nan,
) -> np.ndarray:
    """Perform a 1D gamma of two 1D profiles/arrays. This does NOT check lengths or
    spatial consistency. It performs an element-by-element evaluation. It is the responsibility
    of the caller to ensure the reference and evaluation have comparable spatial resolution.

    The algorithm follows Table I of D. Low's 2004 paper: Evaluation of the gamma dose distribution comparison method: https://aapm.onlinelibrary.wiley.com/doi/epdf/10.1118/1.1598711

    Parameters
    ----------

    reference
        The reference profile.
    evaluation
        The evaluation profile.
    dose_to_agreement
        The dose to agreement in %. E.g. 1 is 1% of global reference max dose.
    distance_to_agreement
        The distance to agreement in **elements**. E.g. if the value is 4 this means 4 elements from the reference point under calculation.
        Must be >0
    gamma_cap_value
        The value to cap the gamma at. E.g. a gamma of 5.3 will get capped to 2. Useful for displaying data with a consistent range.
    global_dose
        Whether to evaluate the dose to agreement threshold based on the global max or the dose point under evaluation.
    dose_threshold
        The dose threshold as a number between 0 and 100 of the % of max dose under which a gamma is not calculated.
        This is not affected by the global/local dose normalization and the threshold value is evaluated against the global max dose, period.
    fill_value
        The value to give pixels that were not calculated because they were under the dose threshold. Default
        is NaN, but another option would be 0. If NaN, allows the user to calculate mean/median gamma over just the
        evaluated portion and not be skewed by 0's that should not be considered.
    """
    if reference.ndim != 1 or evaluation.ndim != 1:
        raise ValueError(
            f"Reference and evaluation arrays must be 1D. Got reference: {reference.ndim} and evaluation: {evaluation.ndim}"
        )
    threshold = reference.max() / 100 * dose_threshold
    # convert dose to agreement to % of global max; ignored later if local dose
    dose_ta = dose_to_agreement / 100 * reference.max()
    # pad eval array on both edges so our search does not go out of bounds
    eval_padded = np.pad(evaluation, distance_to_agreement, mode="edge")
    # iterate over each reference element, computing distance value and dose value
    gamma = []
    for r_idx, ref_point in enumerate(reference):
        # skip if below dose threshold
        if ref_point < threshold:
            gamma.append(fill_value)
            continue
        # we search at the same indices in eval_padded, but remember eval_padded has extra indices on each edge,
        # so this is actually searching from -DTA to +DTA because r_idx in eval_padded is off by distance_to_agreement.
        capital_gammas = []
        for e_idx, eval_point in enumerate(
            eval_padded[r_idx : r_idx + 2 * distance_to_agreement + 1]
        ):
            dist = abs(e_idx - distance_to_agreement)
            dose = eval_point - ref_point
            if not global_dose:
                dose_ta = dose_to_agreement / 100 * ref_point
            capital_gamma = math.sqrt(
                dist**2 / distance_to_agreement**2 + dose**2 / dose_ta**2
            )
            capital_gammas.append(capital_gamma)
        gamma.append(min(min(capital_gammas), gamma_cap_value))
    return np.asarray(gamma)


def stretch(
    array: np.ndarray, min: int = 0, max: int = 1, fill_dtype: np.dtype | None = None
) -> np.ndarray:
    """'Stretch' the profile to the fit a new min and max value and interpolate in between.
    From: http://www.labri.fr/perso/nrougier/teaching/numpy.100/  exercise #17

    .. deprecated:: 3.11

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
    warnings.warn(
        "Using stretch from the profile module is deprecated. Use 'stretch' from the pylinac.core.array_utils module",
        DeprecationWarning,
    )
    new_max = max
    if fill_dtype is not None:
        try:
            di = np.iinfo(fill_dtype)
        except ValueError:
            di = np.finfo(fill_dtype)
        new_max = di.max
    # perfectly normalize the array (0..1). ground, then div by range
    stretched_array = (array - array.min()) / (array.max() - array.min())
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
        """Invert the profile."""
        self.values = utils.invert(self.values)

    def bit_invert(self) -> None:
        """Invert the profile bit-wise."""
        self.values = utils.bit_invert(self.values)

    def normalize(self, norm_val: str | float | None = None) -> None:
        """Normalize the profile to the given value.

        Parameters
        ----------
        norm_val : number or 'max' or None
            If a number, normalize the array to that number. If None, normalizes to the maximum value.
        """
        # backwards compatibility
        if norm_val == "max":
            norm_val = None
        self.values = utils.normalize(self.values, value=norm_val)

    def stretch(self, min: float = 0, max: float = 1) -> None:
        """'Stretch' the profile to the min and max parameter values.

        Parameters
        ----------
        min : number
            The new minimum of the values
        max : number
            The new maximum value.
        """
        self.values = utils.stretch(self.values, min=min, max=max)

    def convert_to_dtype(self, dtype: type[np.dtype]) -> None:
        """Convert the profile datatype to another datatype while retaining the values relative to the datatype min/max"""
        self.values = utils.convert_to_dtype(self.values, dtype=dtype)

    def ground(self) -> float:
        """Ground the profile such that the lowest value is 0.

        Returns
        -------
        float
            The minimum value that was used as the grounding value.
        """
        min_val = self.values.min()
        self.values = utils.ground(self.values)
        return min_val

    def filter(self, size: float = 0.05, kind: str = "median") -> None:
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
        self.values = utils.filter(self.values, size=size, kind=kind)

    def __len__(self):
        return len(self.values)

    def __getitem__(self, items):
        return self.values[items]


class Interpolation(enum.Enum):
    """Interpolation Enum"""

    NONE = None  #:
    LINEAR = "Linear"  #:
    SPLINE = "Spline"  #:


class Normalization(enum.Enum):
    """Normalization method Enum"""

    NONE = None  #:
    GEOMETRIC_CENTER = "Geometric center"  #:
    BEAM_CENTER = "Beam center"  #:
    MAX = "Max"  #:


class Edge(enum.Enum):
    """Edge detection Enum"""

    FWHM = "FWHM"  #:
    INFLECTION_DERIVATIVE = "Inflection Derivative"  #:
    INFLECTION_HILL = "Inflection Hill"  #:


class Centering(enum.Enum):
    """See :ref:`centering`"""

    MANUAL = "Manual"  #:
    BEAM_CENTER = "Beam center"  #:
    GEOMETRIC_CENTER = "Geometric center"  #:


class ProfileBase(ProfileMixin, ABC):
    """The base class for multiple type of profiles. This class should not be instantiated directly.

    We use a base class to avoid having long stacked if statements for the different detection patterns.
    This is also more explicit and extensible."""

    metrics: list[ProfileMetric]
    metric_values: dict[str, float]

    def __init__(
        self,
        values: np.array,
        x_values: np.array | None = None,
        ground: bool = False,
        normalization: str | Normalization = Normalization.NONE,
        interpolation_order: int = 1,
    ):
        """A 1D profile that has one large signal, e.g. a radiation beam profile.
        Signal analysis methods are given, mostly based on FWXM and on Hill function calculations.
        Profiles with multiple peaks are better suited by the MultiProfile class.
        """
        validators.single_dimension(values)
        self.metrics = []
        self.metric_values = {}
        self._interp_order = interpolation_order
        if x_values is None:
            x_values = np.arange(len(values))
        x_diff = np.diff(x_values)
        if x_diff.max() > 0 > x_diff.min():
            raise ValueError("X values must be monotonically increasing or decreasing")
        # resort in case we are in descending order
        sort_idxs = np.argsort(x_values)
        x_values = x_values[sort_idxs]
        values = values[sort_idxs]
        self.values = values
        self.x_values = x_values
        if ground:
            self.values = utils.ground(values)
        if normalization == Normalization.MAX:
            self.normalize()
        elif normalization == Normalization.GEOMETRIC_CENTER:
            center_val = utils.geometric_center_value(self.values)
            self.normalize(center_val)
        elif normalization == Normalization.BEAM_CENTER:
            beam_center_val = self.y_at_x(self.center_idx)
            self.normalize(beam_center_val)

    def x_at_x(self, x: float) -> np.ndarray:
        """Deprecated alias for x_at_x_idx"""
        warnings.warn(
            "x_at_x is deprecated. Use x_at_x_idx instead", DeprecationWarning
        )
        return self.x_at_x_idx(x)

    def x_at_x_idx(self, x: float | np.ndarray) -> np.ndarray | float:
        """Return the physical x-value at the given index. When no x-values are provided, these are the same.
        However, physical dimensions can be different than the index."""
        # UnivariateSpline is the newer approach than interp1d
        # interp1d is legacy and may be removed in the future
        # we create an interpolation for x-values so we can
        # get the original x-value when the index is a float (very often the case)
        f = UnivariateSpline(
            x=np.arange(len(self.x_values)), y=self.x_values, k=self._interp_order, s=0
        )
        new_x = f(x)
        if new_x.size == 1:
            return float(new_x)
        return f(x)

    def x_idx_at_x(self, x: float) -> int:
        """Return the **index** of the x-value closest to the given x-value."""
        return int(np.argmin(np.abs(self.x_values - x)))

    def y_at_x(self, x: float | np.ndarray) -> np.ndarray | float:
        """Interpolated y-values. The x-value is the physical position, not the index. However, if no x-values were provided, these will be the same."""
        f = UnivariateSpline(x=self.x_values, y=self.values, k=self._interp_order, s=0)
        new_y = f(x)
        if new_y.size == 1:
            return float(new_y)
        return new_y

    def x_at_y(self, y: float | np.ndarray, side: str) -> np.ndarray | float:
        """Interpolated y-values. Can use floats as indices."""
        # I can't get UnivariateSpline to work here because it wants strictly increasing
        # data. So we use interp1d instead
        s = self.x_idx_at_x(self.center_idx)
        if side == LEFT:
            f = interp1d(x=self.values[:s], y=self.x_values[:s])
        elif side == RIGHT:
            f = interp1d(x=self.values[s:], y=self.x_values[s:])
        new_x = f(y)
        if new_x.size == 1:
            return float(new_x)
        return f(y)

    @abstractmethod
    def field_edge_idx(self, side: str) -> float:
        """The index of the field edge, given the side and edge detection method."""
        pass

    def field_indices(self, in_field_ratio: float) -> (int, int, int):
        """Return the indices of the left and right edge of the field, given the in-field ratio.
        Importantly, this will use the same rounding behavior as field_values.
        """
        xs = self.field_x_values(in_field_ratio)
        left = xs[0]
        right = xs[-1]
        width = max(right, left) - min(right, left)
        return left, right, width

    def field_x_values(self, in_field_ratio: float) -> np.ndarray:
        """Return the x-values of the field, given the in-field ratio.
        This is helpful when plotting the field to include the proper x-values as well.
        """
        left = self.field_edge_idx(side=LEFT)
        right = self.field_edge_idx(side=RIGHT)
        width = self.field_width_px
        f_left = left + (1 - in_field_ratio) / 2 * width
        f_right = right - (1 - in_field_ratio) / 2 * width
        # use floor/ceil to be inclusive of the edges
        lower_bound = math.floor(min((f_left, f_right)))
        upper_bound = math.ceil(max((f_left, f_right)))
        inner_range = np.nonzero(
            (self.x_values >= lower_bound) & (self.x_values <= upper_bound)
        )[0]
        return self.x_values[inner_range]

    @cached_property
    def center_idx(self) -> float:
        """The center index of the profile. Halfway between the field edges."""
        left = self.field_edge_idx(side=LEFT)
        right = self.field_edge_idx(side=RIGHT)
        return abs(right - left) / 2 + left

    @cached_property
    def field_width_px(self) -> float:
        """The field width of the profile in pixels"""
        left_idx = self.field_edge_idx(side=LEFT)
        right_idx = self.field_edge_idx(side=RIGHT)
        return max(right_idx, left_idx) - min(right_idx, left_idx)

    def field_values(
        self,
        in_field_ratio: float = 0.8,
    ) -> np.ndarray:
        """The array of values of the profile within the 'field' area. This is typically 80% of the detected
        field width."""
        x_values = self.field_x_values(in_field_ratio)
        return self.y_at_x(x_values)

    def as_resampled(
        self, interpolation_factor: float = 10, order: int = 3, **kwargs
    ) -> Any:
        """Resample the profile at a new resolution. Returns a new profile.

        Parameters
        ----------
        interpolation_factor : float
            The factor to zoom the profile by. E.g. 10 means the profile will be 10x larger.
        order : int
            The order of the spline interpolation. 1 is linear, 3 is cubic, etc.

        Warnings
        --------
        This method will respect the input datatype of the numpy array. If the input array is a float, the output array will be a float.
        This can cause issues for int arrays with a small range. E.g. if the range is only 10, interpolation
        will look more step-like than smooth. If this is the case, convert the array to a float before passing it to this method.
        The array is not automatically converted to float in this case to respect the original dtype. However,
        a warning will be produced.
        """
        arr_range = self.values.max() - self.values.min()
        if self.values.dtype != float and arr_range < 100:
            warnings.warn(
                f"Array range is small ({arr_range}) and is not a float. Interpolation may look step-like. "
                f"Consider converting the array to a float before passing it to this method.",
                UserWarning,
            )
        new_y = zoom(
            self.values,
            zoom=interpolation_factor,
            order=order,
            grid_mode=False,
            mode="nearest",
        )
        new_x = np.linspace(self.x_values.min(), self.x_values.max(), len(new_y))

        return type(self)(
            values=new_y,
            x_values=new_x,
            ground=False,
            normalization=Normalization.NONE,
            **kwargs,
        )

    def resample_to(
        self, target_profile: ProfileBase | PhysicalProfileMixin
    ) -> ProfileBase:
        """Resample a target profile to have the same sampling (x-values) rate as the source profile.
        This will return a new target profile with the same x-values as the source profile and with the values
        interpolated to match the source profile.

        For example, this can be used to resample an EPID profile to have the same sampling rate as an ion chamber
        profile or vice versa.

        Requirements
        ------------

        * The range of x-values for the target profile must be within the x-value range of the source profile.
          I.e. no extrapolation of the source profile. For example, if we have a profile of an IC Profile,
          that goes from -15cm to +15cm, we cannot resample onto an EPID profile that goes from -20cm to +20cm.
          To do so, go the other way: resample the EPID profile onto the IC profile.
        """
        # get the absolute x values; depends on the type of profile
        if isinstance(target_profile, PhysicalProfileMixin):
            target_x = target_profile.physical_x_values
        else:
            target_x = target_profile.x_values
        if isinstance(self, PhysicalProfileMixin):
            self_x = self.physical_x_values
        else:
            self_x = self.x_values
        f = InterpolatedUnivariateSpline(self_x, self.values, k=1, ext=2)
        try:
            target_y = f(target_x)
        except ValueError:
            raise ValueError(
                f"The target profile x-values are outside this profiles range. Extrapolation is not allowed. self x-values: {self_x.min()} to {self_x.max()}. target x-values: {target_x.min()} to {target_x.max()}. "
            )
        if isinstance(self, PhysicalProfileMixin):
            output_type = self.__class__.__bases__[-1]
        else:
            output_type = self.__class__
        return output_type(values=target_y, x_values=target_x)

    def plot(
        self,
        show: bool = True,
        axis: plt.Axes | None = None,
        show_field_edges: bool = True,
        show_grid: bool = True,
        show_center: bool = True,
    ) -> plt.Axes:
        """Plot the profile along with relevant overlays to point out features."""
        if axis is None:
            _, axis = plt.subplots()
        axis.plot(self.x_values, self.values, label="Data")
        if show_field_edges:
            axis.axvline(self.field_edge_idx(side=LEFT), ls="--", label="Field Edges")
            axis.axvline(self.field_edge_idx(side=RIGHT), ls="--")
        if show_center:
            axis.axvline(self.center_idx, ls=":", label="Center")
        for metric in self.metrics:
            metric.plot(axis)
        axis.grid(show_grid)
        axis.legend()
        if show:
            plt.show()
        return axis

    def compute(
        self, metrics: Iterable[ProfileMetric] | ProfileMetric
    ) -> Any | dict[str, Any]:
        """Compute metric(s) on the profile.

        Unlike other modules, calling ``compute`` is not strictly necessary.
        Only call it if there are metrics to calculate.

        Parameters
        ----------
        metrics: iterable of ProfileMetric | ProfileMetric
            List of metrics to calculate. If only one metric is desired, it can be passed directly.

        Returns
        -------
        dict | list
            A dictionary of metric names and values if multiple metrics were given.
            If only one metric was given, the value of that metric is returned.

        """
        self.metric_values = {}
        self.metrics = []
        values = {}
        if isinstance(metrics, ProfileMetric):
            metrics = [metrics]
        for metric in metrics:
            metric.inject_profile(self)
            self.metrics.append(metric)
            values[metric.name] = metric.calculate()
        self.metric_values |= values
        if len(values) == 1:
            return list(values.values())[0]
        else:
            return values


class FWXMProfile(ProfileBase):
    """A profile that has one large signal, e.g. a radiation beam profile and data derived from it is based on
    the Full-Width X-Maximum to find the edge indices"""

    def __init__(
        self,
        values: np.array,
        x_values: np.array | None = None,
        ground: bool = False,
        normalization: str | Normalization = Normalization.NONE,
        fwxm_height: float = 50,
    ):
        """A 1D profile that has one large signal, e.g. a radiation beam profile.
        Signal analysis methods are given, mostly based on FWXM and on Hill function calculations.
        Profiles with multiple peaks are better suited by the MultiProfile class.
        """
        self.fwxm_height = fwxm_height
        super().__init__(
            values=values,
            x_values=x_values,
            ground=ground,
            normalization=normalization,
        )

    def field_edge_idx(self, side: Literal["right", "left"]) -> float:
        """The edge index of the given side using the FWXM methodology"""
        _, peak_props = find_peaks(
            self.values, fwxm_height=self.fwxm_height / 100, max_number=1
        )
        if side == LEFT:
            idx = peak_props["left_ips"][0]
        elif side == RIGHT:
            idx = peak_props["right_ips"][0]
        return self.x_at_x_idx(idx)

    def as_resampled(
        self, interpolation_factor: float = 10, order: int = 3
    ) -> FWXMProfile:
        """Resample the profile at a new resolution. Returns a new profile.

        Parameters
        ----------
        interpolation_factor : float
            The factor to zoom the profile by. E.g. 10 means the profile will be 10x larger.
        order : int
            The order of the spline interpolation. 1 is linear, 3 is cubic, etc.
        """
        return super().as_resampled(
            interpolation_factor=interpolation_factor,
            order=order,
            fwxm_height=self.fwxm_height,
        )


class InflectionDerivativeProfile(ProfileBase):
    """A profile that has one large signal, e.g. a radiation beam profile and data derived from it is based on
    the Full-Width X-Maximum"""

    def __init__(
        self,
        values: np.array,
        x_values: np.array | None = None,
        ground: bool = False,
        normalization: str | Normalization = Normalization.NONE,
        edge_smoothing_ratio: float = 0.003,
    ):
        """A 1D profile that has one large signal, e.g. a radiation beam profile.
        Signal analysis methods are given, mostly based on FWXM and on Hill function calculations.
        Profiles with multiple peaks are better suited by the MultiProfile class.
        """
        self.edge_smoothing_ratio = edge_smoothing_ratio
        super().__init__(
            values=values,
            x_values=x_values,
            ground=ground,
            normalization=normalization,
        )

    def field_edge_idx(self, side: str) -> float:
        """The edge index of the given side using the second derivative methodology"""
        filtered_values = gaussian_filter1d(
            self.values, sigma=self.edge_smoothing_ratio * len(self.values)
        )
        diff = np.gradient(filtered_values)
        f_diff = interp1d(x=self.x_values, y=diff, kind="cubic")

        if side == LEFT:
            initial_guess = self.x_at_x_idx(np.argmax(diff))
            idx = minimize(lambda x: -f_diff(x), x0=initial_guess).x[0]
        else:
            initial_guess = self.x_at_x_idx(np.argmin(diff))
            idx = minimize(f_diff, x0=initial_guess).x[0]
        return idx

    def as_resampled(
        self, interpolation_factor: float = 10, order: int = 3
    ) -> InflectionDerivativeProfile:
        return super().as_resampled(
            interpolation_factor=interpolation_factor,
            order=order,
            edge_smoothing_ratio=self.edge_smoothing_ratio,
        )


class HillProfile(InflectionDerivativeProfile):
    """A profile that has one large signal, e.g. a radiation beam profile and data derived from it is based on
    the Full-Width X-Maximum"""

    def __init__(
        self,
        values: np.array,
        x_values: np.array | None = None,
        ground: bool = False,
        normalization: str = Normalization.NONE,
        edge_smoothing_ratio: float = 0.003,
        hill_window_ratio: float = 0.1,
    ):
        """A 1D profile that has one large signal, e.g. a radiation beam profile.
        Signal analysis methods are given, mostly based on FWXM and on Hill function calculations.
        Profiles with multiple peaks are better suited by the MultiProfile class.
        """
        self.hill_window_ratio = hill_window_ratio
        super().__init__(
            values=values,
            x_values=x_values,
            ground=ground,
            normalization=normalization,
            edge_smoothing_ratio=edge_smoothing_ratio,
        )

    def field_edge_idx(self, side: str) -> float:
        """The edge index of the given side using the FWXM methodology"""
        left_infl_idx = super().field_edge_idx(side=LEFT)
        right_infl_idx = super().field_edge_idx(side=RIGHT)
        window_size = (right_infl_idx - left_infl_idx) * self.hill_window_ratio
        if side == LEFT:
            left = int(round(left_infl_idx - window_size))
            right = int(round(left_infl_idx + window_size))
            x_data = self.x_values[left : right + 1]
            y_data = self.values[left : right + 1]
        else:
            left = int(round(right_infl_idx - window_size))
            right = int(round(right_infl_idx + window_size))
            x_data = self.x_values[left : right + 1]
            y_data = self.values[left : right + 1]
        hill_fit = Hill.fit(x_data=x_data, y_data=y_data)
        idx = hill_fit.inflection_idx()["index (exact)"]
        return self.x_at_x_idx(idx)

    def as_resampled(
        self, interpolation_factor: float = 10, order: int = 3
    ) -> HillProfile:
        return ProfileBase.as_resampled(
            self,
            interpolation_factor=interpolation_factor,
            order=order,
            edge_smoothing_ratio=self.edge_smoothing_ratio,
            hill_window_ratio=self.hill_window_ratio,
        )


class PhysicalProfileMixin:
    """A mixin when the profile has a physical component.
    This is pretty typical for EPID profiles, etc. The mixin
    adds a few methods that take physical distance into account."""

    x_values: np.ndarray
    values: np.ndarray
    field_width_px: float

    def __init__(
        self,
        dpmm: float,
    ):
        self.dpmm = dpmm

    @property
    def physical_x_values(self) -> np.array:
        """The x-values of the profile in absolute position, taking into account the dpmm."""
        half_pixel_offset = 0.5 / self.dpmm
        x_values = self.x_values / self.dpmm + half_pixel_offset
        return x_values

    @cached_property
    def field_width_mm(self) -> float:
        """The field width of the profile in mm"""
        return self.field_width_px / self.dpmm

    def as_simple_profile(self) -> ProfileBase:
        """Convert a physical profile into a simple profile where
        the x-values have been converted to the physical x-values.

        An example is converting an EPID profile into a simple profile
        where the x-values are in mm, not pixels.

        This can be useful when trying to compare physical profiles
        to simple profiles. E.g. an EPID vs an ion chamber profile acquisition.
        In the EPID's case, the x-values are indices w/ a dpmm component.
        The IC profile is usually already directly in absolute x-values.
        """
        base_profile_type = self.__class__.__bases__[-1]
        return base_profile_type(values=self.values, x_values=self.physical_x_values)

    def as_resampled(
        self,
        interpolation_resolution_mm: float = 0.1,
        order: int = 3,
        grid: bool = True,
    ) -> Any:
        """Resample the physical profile at a new resolution. Returns a new profile.

        Parameters
        ----------
        interpolation_resolution_mm : float
            The resolution to resample to in mm. E.g. 0.1 means the profile will be 0.1 mm resolution.
        order : int
            The order of the spline interpolation. 1 is linear, 3 is cubic, etc.
        grid : bool
            Whether to use grid mode when zooming. See parameter ``grid_mode`` in :func:`~scipy.ndimage.zoom` for more information.
            This should be true unless you are resampling an already-resampled physical array.

        Warnings
        --------
        This method will respect the input datatype of the numpy array. If the input array is a float, the output array will be a float.
        This can cause issues for int arrays with a small range. E.g. if the range is only 10, interpolation
        will look more step-like than smooth. If this is the case, convert the array to a float before passing it to this method.
        The array is not automatically converted to float in this case to respect the original dtype. However,
        a warning will be produced.
        """
        arr_range = self.values.max() - self.values.min()
        if self.values.dtype != float and arr_range < 100:
            warnings.warn(
                f"Array range is small ({arr_range}) and is not a float. Interpolation may look step-like. "
                f"Consider converting the array to a float before passing it to this method.",
                UserWarning,
            )

        factor = 1 / (self.dpmm * interpolation_resolution_mm)
        #  When dealing with physical arrays where each pixel/voxel is a physical distance, it is important to
        #  use grid mode when zooming. This is because each pixel/voxel has a physical size.
        #  See parameter ``grid_mode`` in :func:`~scipy.ndimage.zoom` for more information.
        new_y = zoom(
            self.values,
            zoom=factor,
            order=order,
            grid_mode=grid,
            mode="nearest",
        )
        # similarly, we assume that x-values are also physical
        # we thus have to offset the x-values by half a pixel/voxel
        # while accounting for the physical size of the pixel/voxel
        if grid:
            offset = 0.5 - 1 / (2 * factor)
            new_x = np.linspace(
                self.x_values.min() - offset, self.x_values.max() + offset, len(new_y)
            )
        else:
            new_x = np.linspace(self.x_values.min(), self.x_values.max(), len(new_y))

        return type(self)(
            values=new_y,
            x_values=new_x,
            ground=False,
            normalization=Normalization.NONE,
            dpmm=factor * self.dpmm,
        )


class FWXMProfilePhysical(PhysicalProfileMixin, FWXMProfile):
    def __init__(
        self,
        values: np.array,
        dpmm: float,
        x_values: np.array | None = None,
        ground: bool = False,
        normalization: str | Normalization = Normalization.NONE,
        fwxm_height: float = 50,
    ):
        FWXMProfile.__init__(
            self,
            values=values,
            x_values=x_values,
            ground=ground,
            normalization=normalization,
            fwxm_height=fwxm_height,
        )
        PhysicalProfileMixin.__init__(self, dpmm=dpmm)

    def as_resampled(
        self,
        interpolation_resolution_mm: float = 0.1,
        order: int = 3,
        grid: bool = True,
    ) -> FWXMProfilePhysical:
        return super().as_resampled(
            interpolation_resolution_mm=interpolation_resolution_mm,
            order=order,
            grid=grid,
        )


class InflectionDerivativeProfilePhysical(
    PhysicalProfileMixin, InflectionDerivativeProfile
):
    def __init__(
        self,
        values: np.array,
        dpmm: float,
        x_values: np.array | None = None,
        ground: bool = False,
        normalization: str | Normalization = Normalization.NONE,
        edge_smoothing_ratio: float = 0.003,
    ):
        InflectionDerivativeProfile.__init__(
            self,
            values=values,
            x_values=x_values,
            ground=ground,
            normalization=normalization,
            edge_smoothing_ratio=edge_smoothing_ratio,
        )
        PhysicalProfileMixin.__init__(self, dpmm=dpmm)

    def as_resampled(
        self,
        interpolation_resolution_mm: float = 0.1,
        order: int = 3,
        grid: bool = True,
    ) -> InflectionDerivativeProfilePhysical:
        return super().as_resampled(
            interpolation_resolution_mm=interpolation_resolution_mm,
            order=order,
            grid=grid,
        )


class HillProfilePhysical(PhysicalProfileMixin, HillProfile):
    def __init__(
        self,
        values: np.array,
        dpmm: float,
        x_values: np.array | None = None,
        ground: bool = False,
        normalization: str | Normalization = Normalization.NONE,
        edge_smoothing_ratio: float = 0.003,
        hill_window_ratio: float = 0.1,
    ):
        HillProfile.__init__(
            self,
            values=values,
            x_values=x_values,
            ground=ground,
            normalization=normalization,
            edge_smoothing_ratio=edge_smoothing_ratio,
            hill_window_ratio=hill_window_ratio,
        )
        PhysicalProfileMixin.__init__(self, dpmm=dpmm)

    def as_resampled(
        self,
        interpolation_resolution_mm: float = 0.1,
        order: int = 3,
        grid: bool = True,
    ) -> HillProfilePhysical:
        return super().as_resampled(
            interpolation_resolution_mm=interpolation_resolution_mm,
            order=order,
            grid=grid,
        )


class SingleProfile(ProfileMixin):
    """A profile that has one large signal, e.g. a radiation beam profile.
    Signal analysis methods are given, mostly based on FWXM and on Hill function calculations.
    Profiles with multiple peaks are better suited by the MultiProfile class.
    """

    def __init__(
        self,
        values: np.ndarray,
        dpmm: float = None,
        interpolation: Interpolation | str | None = Interpolation.LINEAR,
        ground: bool = True,
        interpolation_resolution_mm: float = 0.1,
        interpolation_factor: float = 10,
        normalization_method: Normalization | str = Normalization.BEAM_CENTER,
        edge_detection_method: Edge | str = Edge.FWHM,
        edge_smoothing_ratio: float = 0.003,
        hill_window_ratio: float = 0.1,
        x_values: np.ndarray | None = None,
    ):
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
            The resolution that the interpolation will scale to. **Only used if dpmm is passed and interpolation is set**.
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
            **Only applies to INFLECTION_DERIVATIVE and INFLECTION_HILL.**

            The ratio of the length of the values to use as the sigma for a Gaussian filter applied before searching for
            the inflection. E.g. 0.005 with a profile of 1000 points will result in a sigma of 5.
            This helps make the inflection point detection more robust to noise. Increase for noisy data.
        hill_window_ratio
            The ratio of the field size to use as the window to fit the Hill function. E.g. 0.2 will using a window
            centered about each edge with a width of 20% the size of the field width. **Only applies when the edge
            detection is INFLECTION_HILL**.
        x_values
            The x-values of the profile, if any. If None, will generate a simple range(len(values)).
        """
        self._interp_method = convert_to_enum(interpolation, Interpolation)
        self._interpolation_res = interpolation_resolution_mm
        self._interpolation_factor = interpolation_factor
        self._norm_method = convert_to_enum(normalization_method, Normalization)
        self._edge_method = convert_to_enum(edge_detection_method, Edge)
        self._edge_smoothing_ratio = edge_smoothing_ratio
        self._hill_window_ratio = hill_window_ratio
        self.values = (
            values  # set initial data so we can do things like find beam center
        )
        self.dpmm = dpmm
        fitted_values, new_dpmm, x_indices = self._interpolate(
            values,
            x_values,
            dpmm,
            interpolation_resolution_mm,
            interpolation_factor,
            self._interp_method,
        )
        self.values = fitted_values
        self.x_indices = x_indices
        self._x_interp1d = interp1d(list(range(len(x_indices))), x_indices)
        self._ground = ground
        if ground:
            fitted_values -= fitted_values.min()
        self._y_interp1d = interp1d(
            x_indices, fitted_values, bounds_error=False, fill_value="extrapolate"
        )
        norm_values = self._normalize(fitted_values, self._norm_method)
        self.values = norm_values  # update values
        self._y_interp1d = interp1d(
            x_indices, norm_values, bounds_error=False, fill_value="extrapolate"
        )

    def _x_interp_to_original(self, location: float | np.ndarray) -> float | np.ndarray:
        """Get the x-value of the (possibly) interpolated profile. The input value is in the original
        value range. E.g. a profile with x-range of 0-10 is interpolated to 10x. Asking for the location at 99 would scale back to 9.9.
        We need this function because peak finding is independent of the x-values. I.e. peaks are found and reported according
        to the (0, len(x_values)) range. If the x-values are interpolated we need to get back to the original x-value.
        """
        x = self._x_interp1d(location)
        if isinstance(location, (float, int)) or location.size == 1:
            return float(x)
        return x

    def _y_original_to_interp(self, location: float | np.ndarray) -> float | np.ndarray:
        """Get the interpolated y-value of the profile. This is a corollary to the _x_interp... function"""
        y = self._y_interp1d(location)
        if isinstance(location, (float, int)) or location.size == 1:
            return float(y)
        return y

    def resample(
        self, interpolation_factor: int = 10, interpolation_resolution_mm: float = 0.1
    ) -> SingleProfile:
        """Resample the profile at a new resolution. Returns a new profile"""
        # we have to set the dpmm to what it currently is after interpolating to resample.
        if self.dpmm:
            dpmm = 1 / self._interpolation_res
        else:
            dpmm = None
        return SingleProfile(
            values=self.values,
            x_values=self.x_indices,
            dpmm=dpmm,
            interpolation=self._interp_method,
            ground=self._ground,
            interpolation_resolution_mm=interpolation_resolution_mm,
            interpolation_factor=interpolation_factor,
            normalization_method=self._norm_method,
            edge_detection_method=self._edge_method,
            edge_smoothing_ratio=self._edge_smoothing_ratio,
            hill_window_ratio=self._hill_window_ratio,
        )

    @staticmethod
    def _interpolate(
        values,
        x_values,
        dpmm,
        interpolation_resolution,
        interpolation_factor,
        interp_method: Interpolation,
    ) -> (np.ndarray, float, float, float):
        """Fit the data to the passed interpolation method. Will also calculate the new values to correct the measurements such as dpmm"""
        if x_values is None:
            x_values = np.array(range(len(values)))
        if np.diff(x_values).min() < 0:
            raise ValueError("Profile values must be monotonically increasing")
        if interp_method == Interpolation.NONE:
            return values, dpmm, x_values  # do nothing
        else:
            if dpmm is not None:
                samples = int(round(len(x_values) / (dpmm * interpolation_resolution)))
                new_dpmm = 1 / interpolation_resolution
            else:
                samples = int(round(len(x_values) * interpolation_factor))
                new_dpmm = None
            # Warning: BMF ahead
            # the problem is that if we upsample, the left and right ends are not equally sampled.
            # E.g. upsampling a 3-pixel array (0, 1, 2) by 10 normally results in ~20 elements. You
            # interpolate between 0 and 1, and 1 and 2.
            # The first issue is that you do not have a simple X proportion of
            # elements (3 * 10 = 30 but we get 20). Additionally, if these are pixels they have a
            # finite, physical size and technically those values are at the center of the pixels.
            # Thus, you actually need to sample beyond the left and right edges. In the
            # above case you'd really need to sample from approximately -0.5 to 2.5 to get ~10 pixels
            # for each original pixel. We also need to offset the x-values to be back to 0 again from -0.5.
            # We solve this by offsetting the new x-values by a proportion of the sampling ratio.
            # A ratio of 1 (identical sampling) should not have any offset and return the same values
            # As the ratio goes up, we approach the limit of 0.5 pixels. This follows a proportional relationship
            # with the ratio.
            resampling_factor = samples / len(values)
            offset = 0.5 - 1 / (2 * resampling_factor)
            if interp_method == Interpolation.LINEAR:
                kind = "linear"
            elif interp_method == Interpolation.SPLINE:
                kind = "cubic"
            f = interp1d(
                x_values,
                values,
                kind=kind,
                bounds_error=False,
                fill_value="extrapolate",
            )
            new_x = np.linspace(
                x_values[0] - offset, x_values[-1] + offset, num=samples
            )
            new_y = f(new_x)
            return new_y, new_dpmm, new_x

    def _normalize(self, values, method: Normalization) -> np.ndarray:
        """Normalize the data given a method."""
        if method == Normalization.NONE:
            return values
        elif method == Normalization.MAX:
            return values / values.max()
        elif method == Normalization.GEOMETRIC_CENTER:
            return values / self._geometric_center(values)["value (exact)"]
        elif method == Normalization.BEAM_CENTER:
            return values / self.beam_center()["value (@rounded)"]

    def _geometric_center(self, values) -> dict:
        """Returns the center index and value of the profile.

        If the profile has an even number of values the centre lies between the two centre indices and the centre
        value is the average of the two centre values else the centre index and value are returned.
        """
        return {
            "index (exact)": self._x_interp_to_original(
                utils.geometric_center_idx(values)
            ),
            "value (exact)": utils.geometric_center_value(values),
        }

    def geometric_center(self) -> dict:
        """The geometric center (i.e. the device center)"""
        return self._geometric_center(self.values)

    def beam_center(self) -> dict:
        """The center of the detected beam. This can account for asymmetries in the beam position (e.g. offset jaws)"""
        if self._edge_method == Edge.FWHM:
            data = self.fwxm_data(x=50)
            return {
                "index (rounded)": data["center index (rounded)"],
                "index (exact)": data["center index (exact)"],
                "value (@rounded)": data["center value (@rounded)"],
            }
        elif self._edge_method in (Edge.INFLECTION_DERIVATIVE, Edge.INFLECTION_HILL):
            infl = self.inflection_data()
            mid_point = (
                infl["left index (exact)"]
                + (infl["right index (exact)"] - infl["left index (exact)"]) / 2
            )
            return {
                "index (rounded)": int(round(mid_point)),
                "index (exact)": mid_point,
                "value (@rounded)": self._y_original_to_interp(int(round(mid_point))),
            }

    @argue.bounds(x=(0, 100))
    def fwxm_data(self, x: int = 50) -> dict:
        """Return the width at X-Max, where X is the percentage height.

        Parameters
        ----------
        x
            The percent height of the profile. E.g. x = 50 is 50% height,
            i.e. FWHM.
        """
        _, peak_props = find_peaks(self.values, fwxm_height=x / 100, max_number=1)
        left_idx = float(self._x_interp_to_original(peak_props["left_ips"][0]))
        right_idx = float(self._x_interp_to_original(peak_props["right_ips"][0]))
        width = right_idx - left_idx
        fwxm_center_idx = (right_idx - left_idx) / 2 + left_idx

        data = {
            "width (exact)": width,
            "width (rounded)": int(round(width)),
            "center index (rounded)": int(round(fwxm_center_idx)),
            "center index (exact)": fwxm_center_idx,
            "center value (@rounded)": float(
                self._y_original_to_interp(int(round(fwxm_center_idx)))
            ),
            "left index (exact)": left_idx,
            "left index (rounded)": int(round(left_idx)),
            "left value (@rounded)": float(
                self._y_original_to_interp(int(round(left_idx)))
            ),
            "right index (exact)": right_idx,
            "right index (rounded)": int(round(right_idx)),
            "right value (@rounded)": float(
                self._y_original_to_interp(int(round(right_idx)))
            ),
            "field values": self._y_original_to_interp(
                self.x_indices[int(round(left_idx)) : int(round(right_idx))]
            ),
            "peak_props": peak_props,
        }
        if self.dpmm:
            data["width (exact) mm"] = data["width (exact)"] / self.dpmm
            data["left distance (exact) mm"] = (
                abs(data["center index (exact)"] - data["left index (exact)"])
                / self.dpmm
            )
            data["right distance (exact) mm"] = (
                abs(data["right index (exact)"] - data["center index (exact)"])
                / self.dpmm
            )

        return data

    @argue.bounds(in_field_ratio=(0, 1.0), slope_exclusion_ratio=(0, 1.0))
    def field_data(
        self, in_field_ratio: float = 0.8, slope_exclusion_ratio=0.2
    ) -> dict:
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
            raise ValueError(
                "The exclusion region must be smaller than the field ratio"
            )
        if self._edge_method == Edge.FWHM:
            data = self.fwxm_data(x=50)
            beam_center_idx = data["center index (exact)"]
            full_width = data["width (exact)"]

        elif self._edge_method in (Edge.INFLECTION_DERIVATIVE, Edge.INFLECTION_HILL):
            infl_data = self.inflection_data()
            beam_center_idx = self.beam_center()["index (exact)"]
            full_width = (
                infl_data["right index (exact)"] - infl_data["left index (exact)"]
            )
        beam_center_idx_r = int(round(beam_center_idx))

        cax_idx = self.geometric_center()["index (exact)"]
        cax_idx_r = int(round(cax_idx))

        field_left_idx = beam_center_idx - in_field_ratio * full_width / 2
        field_left_idx_r = int(round(field_left_idx))
        field_right_idx = beam_center_idx + in_field_ratio * full_width / 2
        field_right_idx_r = int(round(field_right_idx))
        field_width = field_right_idx - field_left_idx

        # slope calcs
        inner_left_idx = beam_center_idx - slope_exclusion_ratio * field_width / 2
        inner_left_idx_r = int(round(inner_left_idx))
        inner_right_idx = beam_center_idx + slope_exclusion_ratio * field_width / 2
        inner_right_idx_r = int(round(inner_right_idx))
        left_fit = linregress(
            range(field_left_idx_r, inner_left_idx_r + 1),
            self._y_original_to_interp(
                np.arange(field_left_idx_r, inner_left_idx_r + 1)
            ),
        )
        right_fit = linregress(
            range(inner_right_idx_r, field_right_idx_r + 1),
            self._y_original_to_interp(
                np.arange(inner_right_idx_r, field_right_idx_r + 1)
            ),
        )

        # top calc
        fit_params = np.polyfit(
            range(inner_left_idx_r, inner_right_idx_r + 1),
            self._y_original_to_interp(
                np.arange(inner_left_idx_r, inner_right_idx_r + 1)
            ),
            deg=2,
        )
        width = abs(inner_right_idx_r - inner_left_idx_r)

        def poly_func(x):
            # return the negative since we're MINIMIZING and want the top value
            return -(fit_params[0] * (x**2) + fit_params[1] * x + fit_params[2])

        # minimize the polynomial function
        min_f = minimize(
            poly_func,
            x0=(inner_left_idx_r + width / 2,),
            bounds=((inner_left_idx_r, inner_right_idx_r),),
        )
        top_idx = min_f.x[0]
        top_val = -min_f.fun

        data = {
            "width (exact)": field_width,
            "beam center index (exact)": beam_center_idx,
            "beam center index (rounded)": beam_center_idx_r,
            "beam center value (@rounded)": self._y_original_to_interp(
                round(beam_center_idx)
            ),
            "cax index (exact)": cax_idx,
            "cax index (rounded)": cax_idx_r,
            "cax value (@rounded)": self._y_original_to_interp(round(cax_idx)),
            "left index (exact)": field_left_idx,
            "left index (rounded)": field_left_idx_r,
            "left value (@rounded)": self._y_original_to_interp(round(field_left_idx)),
            "left slope": left_fit.slope,
            "left intercept": left_fit.intercept,
            "right slope": right_fit.slope,
            "right intercept": right_fit.intercept,
            "left inner index (exact)": inner_left_idx,
            "left inner index (rounded)": inner_left_idx_r,
            "right inner index (exact)": inner_right_idx,
            "right inner index (rounded)": inner_right_idx_r,
            '"top" index (exact)': top_idx,
            '"top" index (rounded)': int(round(top_idx)),
            '"top" value (@exact)': top_val,
            "top params": fit_params,
            "right index (exact)": field_right_idx,
            "right index (rounded)": field_right_idx_r,
            "right value (@rounded)": self._y_original_to_interp(
                round(field_right_idx)
            ),
            "field values": self._y_original_to_interp(
                np.arange(int(round(field_left_idx)), int(round(field_right_idx)) + 1)
            ),
        }
        if self.dpmm:
            data["width (exact) mm"] = data["width (exact)"] / self.dpmm
            data["left slope (%/mm)"] = data["left slope"] * self.dpmm * 100
            data["right slope (%/mm)"] = data["right slope"] * self.dpmm * 100
            data["left distance->beam center (exact) mm"] = (
                abs(data["beam center index (exact)"] - data["left index (exact)"])
                / self.dpmm
            )
            data["right distance->beam center (exact) mm"] = (
                abs(data["right index (exact)"] - data["beam center index (exact)"])
                / self.dpmm
            )
            data["left distance->CAX (exact) mm"] = (
                abs(data["cax index (exact)"] - data["left index (exact)"]) / self.dpmm
            )
            data["right distance->CAX (exact) mm"] = (
                abs(data["cax index (exact)"] - data["right index (exact)"]) / self.dpmm
            )

            data["left distance->top (exact) mm"] = (
                abs(data['"top" index (exact)'] - data["left index (exact)"])
                / self.dpmm
            )
            data["right distance->top (exact) mm"] = (
                abs(data['"top" index (exact)'] - data["right index (exact)"])
                / self.dpmm
            )
            data['"top"->beam center (exact) mm'] = (
                data['"top" index (exact)'] - data["beam center index (exact)"]
            ) / self.dpmm
            data['"top"->CAX (exact) mm'] = (
                abs(data['"top" index (exact)'] - data["cax index (exact)"]) / self.dpmm
            )
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
            raise ValueError(
                "FWHM edge method does not have inflection points. Use a different edge detection method"
            )
        d1 = np.gradient(
            gaussian_filter1d(
                self.values, sigma=self._edge_smoothing_ratio * len(self.values)
            )
        )
        (peak_idxs, _) = MultiProfile(d1).find_peaks(threshold=0.8)
        (valley_idxs, _) = MultiProfile(d1).find_valleys(threshold=0.8)
        left_idx = self._x_interp_to_original(peak_idxs[0])  # left-most index
        right_idx = self._x_interp_to_original(valley_idxs[-1])  # right-most index
        if self._edge_method == Edge.INFLECTION_DERIVATIVE:
            data = {
                "left index (rounded)": int(round(left_idx)),
                "left index (exact)": left_idx,
                "right index (rounded)": int(round(right_idx)),
                "right index (exact)": right_idx,
                "left value (@rounded)": self._y_original_to_interp(
                    int(round(left_idx))
                ),
                "left value (@exact)": self._y_original_to_interp(left_idx),
                "right value (@rounded)": self._y_original_to_interp(
                    int(round(right_idx))
                ),
                "right value (@exact)": self._y_original_to_interp(right_idx),
            }
            return data
        else:  # Hill
            # the 2nd deriv is a good approximation for the inflection point. Start there and fit Hill about it
            # penum_half_window = self.field_data()['width (exact)'] * self._hill_window_ratio / 2
            penum_half_window = int(
                round(self._hill_window_ratio * abs(right_idx - left_idx) / 2)
            )

            # left side
            x_data = np.array(
                [
                    x
                    for x in np.arange(
                        left_idx - penum_half_window, left_idx + penum_half_window
                    )
                    if x >= 0
                ]
            )
            y_data = self._y_original_to_interp(x_data)
            # y_data = self.values[left_idx - penum_half_window: left_idx + penum_half_window]
            left_hill = Hill.fit(x_data, y_data)
            left_infl = left_hill.inflection_idx()

            # right side
            x_data = np.array(
                [
                    x
                    for x in np.arange(
                        right_idx - penum_half_window, right_idx + penum_half_window
                    )
                    if x < len(d1)
                ]
            )
            y_data = self._y_original_to_interp(x_data)
            right_hill = Hill.fit(x_data, y_data)
            right_infl = right_hill.inflection_idx()

            data = {
                "left index (rounded)": left_infl["index (rounded)"],
                "left index (exact)": left_infl["index (exact)"],
                "right index (rounded)": right_infl["index (rounded)"],
                "right index (exact)": right_infl["index (exact)"],
                "left value (@exact)": left_hill.y(left_infl["index (exact)"]),
                "right value (@exact)": right_hill.y(right_infl["index (exact)"]),
                "left Hill params": left_hill.params,
                "right Hill params": right_hill.params,
            }
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
            raise ValueError(
                "Upper penumbra value must be larger than the lower penumbra value"
            )
        if self._edge_method == Edge.FWHM:
            upper_data = self.fwxm_data(x=upper)
            lower_data = self.fwxm_data(x=lower)
            data = {
                f"left {lower}% index (exact)": lower_data["left index (exact)"],
                f"left {lower}% value (@rounded)": lower_data["left value (@rounded)"],
                f"left {upper}% index (exact)": upper_data["left index (exact)"],
                f"left {upper}% value (@rounded)": upper_data["left value (@rounded)"],
                f"right {lower}% index (exact)": lower_data["right index (exact)"],
                f"right {lower}% value (@rounded)": lower_data[
                    "right value (@rounded)"
                ],
                f"right {upper}% index (exact)": upper_data["right index (exact)"],
                f"right {upper}% value (@rounded)": upper_data[
                    "right value (@rounded)"
                ],
                "left values": self.values[
                    lower_data["left index (rounded)"] : upper_data[
                        "left index (rounded)"
                    ]
                ],
                "right values": self.values[
                    upper_data["right index (rounded)"] : lower_data[
                        "right index (rounded)"
                    ]
                ],
                "left penumbra width (exact)": abs(
                    upper_data["left index (exact)"] - lower_data["left index (exact)"]
                ),
                "right penumbra width (exact)": abs(
                    upper_data["right index (exact)"]
                    - lower_data["right index (exact)"]
                ),
            }
            if self.dpmm:
                data["left penumbra width (exact) mm"] = (
                    data["left penumbra width (exact)"] / self.dpmm
                )
                data["right penumbra width (exact) mm"] = (
                    data["right penumbra width (exact)"] / self.dpmm
                )
            return data
        elif self._edge_method == Edge.INFLECTION_DERIVATIVE:
            infl_data = self.inflection_data()
            lower_left_percent = max(
                infl_data["left value (@exact)"] / self.values.max() * lower / 50 * 100,
                1,
            )
            upper_left_percent = min(
                infl_data["left value (@exact)"] / self.values.max() * upper / 50 * 100,
                99,
            )
            upper_left_data = self.fwxm_data(x=upper_left_percent)
            lower_left_data = self.fwxm_data(x=lower_left_percent)

            lower_right_value = max(
                infl_data["right value (@exact)"]
                / self.values.max()
                * lower
                / 50
                * 100,
                1,
            )
            upper_right_value = min(
                infl_data["right value (@exact)"]
                / self.values.max()
                * upper
                / 50
                * 100,
                99,
            )
            upper_right_data = self.fwxm_data(x=upper_right_value)
            lower_right_data = self.fwxm_data(x=lower_right_value)

            data = {
                f"left {lower}% index (exact)": lower_left_data["left index (exact)"],
                f"left {upper}% index (exact)": upper_left_data["left index (exact)"],
                f"right {lower}% index (exact)": lower_right_data[
                    "right index (exact)"
                ],
                f"right {upper}% index (exact)": upper_right_data[
                    "right index (exact)"
                ],
                "left values": self._y_original_to_interp(
                    np.arange(
                        lower_left_data["left index (rounded)"],
                        upper_left_data["left index (rounded)"],
                    )
                ),
                "right values": self._y_original_to_interp(
                    np.arange(
                        upper_right_data["right index (rounded)"],
                        lower_right_data["right index (rounded)"],
                    )
                ),
                "left penumbra width (exact)": abs(
                    upper_left_data["left index (exact)"]
                    - lower_left_data["left index (exact)"]
                ),
                "right penumbra width (exact)": abs(
                    upper_right_data["right index (exact)"]
                    - lower_right_data["right index (exact)"]
                ),
            }
            if self.dpmm:
                data["left penumbra width (exact) mm"] = (
                    data["left penumbra width (exact)"] / self.dpmm
                )
                data["right penumbra width (exact) mm"] = (
                    data["right penumbra width (exact)"] / self.dpmm
                )
            return data
        elif self._edge_method == Edge.INFLECTION_HILL:
            infl_data = self.inflection_data()
            left_hill = Hill.from_params(infl_data["left Hill params"])
            right_hill = Hill.from_params(infl_data["right Hill params"])

            lower_left_percent = infl_data["left value (@exact)"] * lower / 50
            lower_left_index = left_hill.x(lower_left_percent)
            upper_left_percent = infl_data["left value (@exact)"] * upper / 50
            upper_left_index = left_hill.x(upper_left_percent)

            lower_right_value = infl_data["right value (@exact)"] * lower / 50
            lower_right_index = right_hill.x(lower_right_value)
            upper_right_value = infl_data["right value (@exact)"] * upper / 50
            upper_right_index = right_hill.x(upper_right_value)

            data = {
                f"left {lower}% index (exact)": lower_left_index,
                f"left {lower}% value (exact)": lower_left_percent,
                f"left {upper}% index (exact)": upper_left_index,
                f"left {upper}% value (exact)": upper_left_percent,
                f"right {lower}% index (exact)": lower_right_index,
                f"right {lower}% value (exact)": lower_right_value,
                f"right {upper}% index (exact)": upper_right_index,
                f"right {upper}% value (exact)": upper_right_value,
                "left values": self.values[
                    int(round(lower_left_index)) : int(round(upper_left_index))
                ],
                "right values": self.values[
                    int(round(upper_right_index)) : int(round(lower_right_index))
                ],
                "left penumbra width (exact)": abs(upper_left_index - lower_left_index),
                "right penumbra width (exact)": abs(
                    upper_right_index - lower_right_index
                ),
                "left gradient (exact)": left_hill.gradient_at(
                    infl_data["left index (exact)"]
                ),
                r"right gradient (exact)": right_hill.gradient_at(
                    infl_data["right index (exact)"]
                ),
            }
            if self.dpmm:
                data["left penumbra width (exact) mm"] = (
                    data["left penumbra width (exact)"] / self.dpmm
                )
                data["left gradient (exact) %/mm"] = (
                    data["left gradient (exact)"] * self.dpmm * 100
                )  # 100 to convert to %
                data["right penumbra width (exact) mm"] = (
                    data["right penumbra width (exact)"] / self.dpmm
                )
                data["right gradient (exact) %/mm"] = (
                    data["right gradient (exact)"] * self.dpmm * 100
                )
            return data

    @argue.options(calculation=("mean", "median", "max", "min", "area"))
    def field_calculation(
        self,
        in_field_ratio: float = 0.8,
        calculation: str = "mean",
        slope_exclusion_ratio: float = 0.2,
    ) -> float | tuple[float, float]:
        """Perform an operation on the field values of the profile.
        This function is useful for determining field symmetry and flatness.

        Parameters
        ----------
        in_field_ratio
            Ratio of the field width to use in the calculation.
        calculation : {'mean', 'median', 'max', 'min', 'area'}
            Calculation to perform on the field values.
        """
        field_values = self.field_data(
            in_field_ratio, slope_exclusion_ratio=slope_exclusion_ratio
        )

        if calculation == "mean":
            return field_values["field values"].mean()
        elif calculation == "median":
            return float(np.median(field_values["field values"]))
        elif calculation == "max":
            return field_values["field values"].max()
        elif calculation == "min":
            return field_values["field values"].min()

    def gamma(
        self,
        evaluation_profile: SingleProfile,
        distance_to_agreement: float = 1,
        dose_to_agreement: float = 1,
        gamma_cap_value: float = 2,
        dose_threshold: float = 5,
        global_dose: bool = True,
        fill_value: float = np.nan,
    ) -> np.ndarray:
        """Calculate a 1D gamma. The passed profile is the evaluation profile. The instance calling this method is the reference profile.
        This profile must have the `dpmm` value given at instantiation so that physical spacing can be evaluated.
        The evaluation profile is resampled to be the same resolution as the reference profile.

        .. note::

            The difference between this method and the `gamma_1d` function is that 1) this is computed on Profile instances and 2)
            this validates the physical spacing of the profiles.

        Parameters
        ----------
        evaluation_profile
            The evaluation profile. This profile must have the `dpmm` value given at instantiation so that physical spacing can be evaluated.
        distance_to_agreement
            Distance in **mm** to search
        dose_to_agreement
            Dose in % of either global or local reference dose
        gamma_cap_value
            The value to cap the gamma at. E.g. a gamma of 5.3 will get capped to 2. Useful for displaying data with a consistent range.
        global_dose
            Whether to evaluate the dose to agreement threshold based on the global max or the dose point under evaluation.
        dose_threshold
            The dose threshold as a number between 0 and 100 of the % of max dose under which a gamma is not calculated.
            This is not affected by the global/local dose normalization and the threshold value is evaluated against the global max dose, period.
        fill_value
            The value to give pixels that were not calculated because they were under the dose threshold. Default
            is NaN, but another option would be 0. If NaN, allows the user to calculate mean/median gamma over just the
            evaluated portion and not be skewed by 0's that should not be considered.
        """
        if not self.dpmm or not evaluation_profile.dpmm:
            raise ValueError(
                "At least one profile does not have the dpmm attribute. Physical spacing cannot be determined. Set it before performing gamma analysis."
            )
        distance_to_agreement_px = int(round(distance_to_agreement * self.dpmm))
        # resample eval profile to be same resolution as reference
        resampled_evaluation = evaluation_profile.resample(
            interpolation_resolution_mm=self._interpolation_res
        )
        if len(resampled_evaluation.values) != len(self.values):
            warnings.warn(
                f"The number of elements in the reference and evaluation differ. Ref: {len(self.values)}, Eval: {len(resampled_evaluation.values)}"
            )
        # now that we've resampled, it's still possible that the x-values of the two profiles differ.
        # E.g. we may be at -0.475 and -0.37 for the first index depending on the amount of interpolation.
        # we thus need to evaluate the evaluation profile at the exact same x-indices as the reference.
        eval_at_ref_points = resampled_evaluation._y_original_to_interp(self.x_indices)
        return gamma_1d(
            reference=self.values,
            evaluation=eval_at_ref_points,
            dose_to_agreement=dose_to_agreement,
            distance_to_agreement=distance_to_agreement_px,
            gamma_cap_value=gamma_cap_value,
            global_dose=global_dose,
            dose_threshold=dose_threshold,
            fill_value=fill_value,
        )

    def plot(self, show: bool = True) -> None:
        """Plot the profile."""
        plt.plot(self.x_indices, self.values)
        if show:
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

    values: np.ndarray | Sequence
    peaks: list
    valleys: list

    def __init__(self, values: np.array | Sequence):
        """
        Parameters
        ----------
        values : iterable
            Array of profile values.
        """
        self.values = values
        self.peaks = []
        self.valleys = []

    def plot(self, ax: plt.Axes | None = None) -> None:
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

    def find_peaks(
        self,
        threshold: float | int = 0.3,
        min_distance: float | int = 0.05,
        max_number: int = None,
        search_region: tuple = (0.0, 1.0),
        peak_sort: str = "prominences",
    ) -> tuple[np.array, np.array]:
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
        peak_idxs, peak_props = find_peaks(
            self.values,
            threshold=threshold,
            peak_separation=min_distance,
            max_number=max_number,
            search_region=search_region,
            peak_sort=peak_sort,
        )
        self.peaks = [
            Point(value=peak_val, idx=peak_idx)
            for peak_idx, peak_val in zip(peak_idxs, peak_props["peak_heights"])
        ]

        return peak_idxs, peak_props["peak_heights"]

    def find_valleys(
        self,
        threshold: float | int = 0.3,
        min_distance: float | int = 0.05,
        max_number: int = None,
        search_region: tuple = (0.0, 1.0),
    ) -> tuple[np.array, np.array]:
        """Find the valleys (minimums) of the profile using a simple minimum value search.

        Returns
        -------
        indices: ndarray, values, ndarray
            The indices and values of the valleys.

        See Also
        --------
        :meth:`~pylinac.core.profile.MultiProfile.find_peaks` : Further parameter info.
        """
        valley_idxs, valley_props = find_peaks(
            -self.values,
            threshold=threshold,
            peak_separation=min_distance,
            max_number=max_number,
            search_region=search_region,
        )
        self.valleys = [
            Point(value=self.values[valley_idx], idx=valley_idx)
            for valley_idx, valley_val in zip(
                valley_idxs, -valley_props["peak_heights"]
            )
        ]

        return valley_idxs, self.values[valley_idxs]

    def find_fwxm_peaks(
        self,
        threshold: float | int = 0.3,
        min_distance: float | int = 0.05,
        max_number: int = None,
        search_region: tuple = (0.0, 1.0),
        peak_sort: str = "prominences",
        required_prominence=None,
    ) -> tuple[np.array, np.array]:
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
        _, peak_props = find_peaks(
            self.values,
            threshold=threshold,
            min_width=min_distance,
            max_number=max_number,
            search_region=search_region,
            peak_sort=peak_sort,
            required_prominence=required_prominence,
        )
        fwxm_peak_idxs = []
        for lt, rt in zip(peak_props["left_ips"], peak_props["right_ips"]):
            fwxm = int(round(lt + (rt - lt) / 2))
            fwxm_peak_idxs.append(fwxm)

        fwxm_peak_vals = [self.values[fwxm] for fwxm in fwxm_peak_idxs]
        self.peaks = [
            Point(value=peak_val, idx=peak_idx)
            for peak_idx, peak_val in zip(fwxm_peak_idxs, fwxm_peak_vals)
        ]

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

    image_array: np.array
    start_angle: float | int
    ccw: bool
    sampling_ratio: float
    _x_locations: np.array | None
    _y_locations: np.array | None

    def __init__(
        self,
        center: Point,
        radius: float,
        image_array: np.array,
        start_angle: float | int = 0,
        ccw: bool = True,
        sampling_ratio: float = 1.0,
    ):
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
        self._ensure_array_size(
            image_array, self.radius + self.center.x, self.radius + self.center.y
        )
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
    def _radians(self) -> np.array:
        interval = (2 * np.pi) / self.size
        rads = np.arange(
            0 + self.start_angle, (2 * np.pi) + self.start_angle - interval, interval
        )
        if self.ccw:
            rads = rads[::-1]
        return rads

    @property
    def x_locations(self) -> np.array:
        """The x-locations of the profile values."""
        if self._x_locations is None:
            return np.cos(self._radians) * self.radius + self.center.x
        else:
            return self._x_locations

    @x_locations.setter
    def x_locations(self, array: np.array):
        self._x_locations = array

    @property
    def y_locations(self) -> np.array:
        """The x-locations of the profile values."""
        if self._y_locations is None:
            return np.sin(self._radians) * self.radius + self.center.y
        else:
            return self._y_locations

    @y_locations.setter
    def y_locations(self, array: np.array):
        self._y_locations = array

    @property
    def _profile(self) -> np.array:
        """The actual profile array; private attr that is passed to MultiProfile."""
        return ndimage.map_coordinates(
            self.image_array, [self.y_locations, self.x_locations], order=0
        )

    def find_peaks(
        self,
        threshold: float | int = 0.3,
        min_distance: float | int = 0.05,
        max_number: int = None,
        search_region: tuple[float, float] = (0.0, 1.0),
    ) -> tuple[np.array, np.array]:
        """Overloads Profile to also map peak locations to the image."""
        peak_idxs, peak_vals = super().find_peaks(
            threshold, min_distance, max_number, search_region
        )
        self._map_peaks()
        return peak_idxs, peak_vals

    def find_valleys(
        self,
        threshold: float | int = 0.3,
        min_distance: float | int = 0.05,
        max_number: int = None,
        search_region: tuple[float, float] = (0.0, 1.0),
    ) -> tuple[np.array, np.array]:
        """Overload Profile to also map valley locations to the image."""
        valley_idxs, valley_vals = super().find_valleys(
            threshold, min_distance, max_number, search_region
        )
        self._map_peaks()
        return valley_idxs, valley_vals

    def find_fwxm_peaks(
        self,
        threshold: float | int = 0.3,
        min_distance: float | int = 0.05,
        max_number: int = None,
        search_region: tuple[float, float] = (0.0, 1.0),
    ) -> tuple[np.array, np.array]:
        """Overloads Profile to also map the peak locations to the image."""
        peak_idxs, peak_vals = super().find_fwxm_peaks(
            threshold, min_distance, max_number, search_region=search_region
        )
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

    def plot2axes(
        self,
        axes: plt.Axes = None,
        edgecolor: str = "black",
        fill: bool = False,
        plot_peaks: bool = True,
    ) -> None:
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
            mpl_Circle(
                (self.center.x, self.center.y),
                edgecolor=edgecolor,
                radius=self.radius,
                fill=fill,
            )
        )
        if plot_peaks:
            x_locs = [peak.x for peak in self.peaks]
            y_locs = [peak.y for peak in self.peaks]
            axes.autoscale(enable=False)
            axes.scatter(x_locs, y_locs, s=40, marker="x", c=edgecolor)

    @staticmethod
    def _ensure_array_size(
        array: np.array, min_width: float, min_height: float
    ) -> None:
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
    def __init__(
        self,
        center: Point,
        radius: float,
        image_array: np.ndarray,
        start_angle: int = 0,
        ccw: bool = True,
        sampling_ratio: float = 1.0,
        width_ratio: float = 0.1,
        num_profiles: int = 20,
    ):
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
    def _radii(self) -> np.array:
        return np.linspace(
            start=self.radius * (1 - self.width_ratio),
            stop=self.radius * (1 + self.width_ratio),
            num=self.num_profiles,
        )

    @property
    def size(self) -> float:
        return np.pi * max(self._radii) * 2 * self.sampling_ratio

    @property
    def _multi_x_locations(self) -> list:
        """List of x-locations of the sampling profiles"""
        x = []
        cos = np.cos(self._radians)
        # extract profile for each circle radii
        for radius in self._radii:
            x.append(cos * radius + self.center.x)
        return x

    @property
    def _multi_y_locations(self) -> list:
        """List of x-locations of the sampling profiles"""
        y = []
        sin = np.sin(self._radians)
        # extract profile for each circle radii
        for radius in self._radii:
            y.append(sin * radius + self.center.y)
        return y

    @property
    def _profile(self) -> np.array:
        """The actual profile array; private attr that is passed to MultiProfile."""
        profile = np.zeros(len(self._multi_x_locations[0]))
        for radius, x, y in zip(
            self._radii, self._multi_x_locations, self._multi_y_locations
        ):
            profile += ndimage.map_coordinates(self.image_array, [y, x], order=0)
        profile /= self.num_profiles
        return profile

    def plot2axes(
        self,
        axes: plt.Axes = None,
        edgecolor: str = "black",
        fill: bool = False,
        plot_peaks: bool = True,
    ) -> None:
        """Add 2 circles to the axes: one at the maximum and minimum radius of the ROI.

        See Also
        --------
        :meth:`~pylinac.core.profile.CircleProfile.plot2axes` : Further parameter info.
        """
        if axes is None:
            fig, axes = plt.subplots()
            axes.imshow(self.image_array)
        axes.add_patch(
            mpl_Circle(
                (self.center.x, self.center.y),
                edgecolor=edgecolor,
                radius=self.radius * (1 + self.width_ratio),
                fill=fill,
            )
        )
        axes.add_patch(
            mpl_Circle(
                (self.center.x, self.center.y),
                edgecolor=edgecolor,
                radius=self.radius * (1 - self.width_ratio),
                fill=fill,
            )
        )
        if plot_peaks:
            x_locs = [peak.x for peak in self.peaks]
            y_locs = [peak.y for peak in self.peaks]
            axes.autoscale(enable=False)
            axes.scatter(x_locs, y_locs, s=20, marker="x", c=edgecolor)


def find_peaks(
    values: np.array,
    threshold: float | int = -np.inf,
    peak_separation: float | int = 0,
    max_number: int | None = None,
    fwxm_height: float = 0.5,
    min_width: int = 0,
    search_region: tuple[float, float] = (0.0, 1.0),
    peak_sort: str = "prominences",
    required_prominence: float | np.array | None = None,
) -> tuple[np.array, dict]:
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
        The method of sorting peaks. See scipy find_peaks documentation.
    required_prominence
        The relative height of the peak compared to surrounding valleys to be considered a peak. See scipy's find_peaks documentation.

    Returns
    -------
    peak_idxs : numpy.array
        The indices of the peaks found.
    peak_props : dict
        A dict containing contextual peak data.
    """
    peak_separation, shift_amount, threshold, trimmed_values = _parse_peak_args(
        peak_separation, search_region, threshold, values
    )

    peak_idxs, peak_props = signal.find_peaks(
        trimmed_values,
        rel_height=(1 - fwxm_height),
        width=min_width,
        height=threshold,
        distance=peak_separation,
        prominence=required_prominence,
    )
    peak_idxs += shift_amount  # shift according to the search region left edge

    # get the "largest" peaks up to max number, and then re-sort to be left->right like it was originally
    largest_peak_idxs = sorted(
        list(np.argsort(peak_props[peak_sort]))[::-1][:max_number]
    )

    # cut down prop arrays as need be
    for key, array_vals in peak_props.items():
        peak_props[key] = array_vals[largest_peak_idxs]
    return peak_idxs[largest_peak_idxs], peak_props


def _parse_peak_args(
    peak_separation: float,
    search_region: tuple[float, float],
    threshold: float,
    values: np.array,
) -> tuple[float, int, float, np.array]:
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
        values = values[
            int(search_region[0] * len(values)) : int(search_region[1] * len(values))
        ]
    else:
        values = values[search_region[0] : search_region[1]]
        shift_amount = search_region[0]
    return peak_separation, shift_amount, threshold, values
