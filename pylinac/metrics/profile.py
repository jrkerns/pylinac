from __future__ import annotations

import math
import typing
from abc import ABC, abstractmethod
from functools import cached_property
from typing import Any, Literal

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from scipy.interpolate import UnivariateSpline
from scipy.optimize import minimize

if typing.TYPE_CHECKING:
    from ..core.profile import PhysicalProfileMixin, ProfileBase


RIGHT = "right"
LEFT = "left"


class ProfileMetric(ABC):
    """Abstract base class for profile metrics. A profile metric is a value that can be calculated from a profile
    and potentially has plot features associated with it.
    Examples include penumbra, flatness, and symmetry"""

    name: str
    unit: str = ""
    profile: ProfileBase | PhysicalProfileMixin

    def __init__(self, color: str | None = None, linestyle: str | None = None):
        self.color = color
        self.linestyle = linestyle

    def inject_profile(self, profile: ProfileBase) -> None:
        """Inject the profile into the metric class.
        We can't do this at instantiation because we don't have
        the profile yet. We also don't want to force the user
        to have to save it manually as they might forget.
        Finally, we want to have it around for any method we might use."""
        self.profile = profile

    def plot(self, axis: plt.Axes):
        """Plot the metric on the given axis."""
        pass

    @abstractmethod
    def calculate(self) -> Any:
        """Calculate the metric on the given profile."""
        pass


class FlatnessDifferenceMetric(ProfileMetric):
    """Flatness as defined by IAEA Rad Onc Handbook pg 196: https://www-pub.iaea.org/MTCD/Publications/PDF/Pub1196_web.pdf"""

    name = "Flatness (Difference)"
    unit = "%"

    def __init__(self, in_field_ratio: float = 0.8, color="g", linestyle="-."):
        self.in_field_ratio = in_field_ratio
        super().__init__(color=color, linestyle=linestyle)

    def calculate(self) -> float:
        """Calculate the flatness ratio of the profile."""
        return (
            100
            * (self.profile.field_values().max() - self.profile.field_values().min())
            / (self.profile.field_values().max() + self.profile.field_values().min())
        )

    def plot(self, axis: plt.Axes) -> None:
        """Plot the points of largest flattness difference as well as the search bounding box."""
        data = self.profile.field_values()
        left, _, width = self.profile.field_indices(in_field_ratio=self.in_field_ratio)
        # plot the search bounding box
        axis.add_patch(
            Rectangle(
                (left, np.min(data)),
                width,
                np.max(data) - np.min(data),
                fill=False,
                color=self.color,
                label=self.label + " Bounding box",
            )
        )
        # plot the max and min values
        axis.plot(
            [np.argmax(data) + left, np.argmin(data) + left],
            [np.max(data), np.min(data)],
            "o",
            color=self.color,
            label=self.name,
        )


class FlatnessRatioMetric(FlatnessDifferenceMetric):
    """Flatness as (apparently) defined by IEC."""

    name = "Flatness (Ratio)"

    def calculate(self) -> float:
        """Calculate the flatness ratio of the profile."""
        return (
            100 * self.profile.field_values().max() / self.profile.field_values().min()
        )


class SymmetryPointDifferenceMetric(ProfileMetric):
    """Symmetry using the point difference method."""

    unit = "%"
    name = "Point Difference Symmetry"

    def __init__(
        self,
        in_field_ratio: float = 0.8,
        color="magenta",
        linestyle="--",
        max_sym_range: float = 2,
        min_sym_range: float = -2,
    ):
        self.in_field_ratio = in_field_ratio
        self.max_sym = max_sym_range
        self.min_sym = min_sym_range
        super().__init__(color=color, linestyle=linestyle)

    @staticmethod
    def _calc_point(lt: float, rt: float, cax: float) -> float:
        return 100 * (lt - rt) / cax

    @cached_property
    def symmetry_values(self) -> list[float]:
        field_values = self.profile.field_values(in_field_ratio=self.in_field_ratio)
        cax_value = self.profile.y_at_x(self.profile.center_idx)
        return [
            self._calc_point(lt, rt, cax_value)
            for lt, rt in zip(field_values, field_values[::-1])
        ]

    def calculate(self) -> float:
        """Calculate the symmetry ratio of the profile."""
        max_sym_idx = np.argmax(np.abs(self.symmetry_values))
        return self.symmetry_values[max_sym_idx]

    def plot(self, axis: plt.Axes, markers: (str, str) = ("^", "v")) -> None:
        idx = np.argmax(self.symmetry_values)
        left_edge, right_edge, _ = self.profile.field_indices(
            in_field_ratio=self.in_field_ratio
        )
        # plot max sym value
        max_x = self.profile.x_at_x_idx(self.profile.x_idx_at_x(left_edge) + idx)
        axis.plot(
            max_x,
            self.profile.y_at_x(max_x),
            markers[0],
            color=self.color,
            label=self.name,
        )
        # plot min sym value
        min_x = self.profile.x_at_x_idx(self.profile.x_idx_at_x(right_edge) - idx)
        axis.plot(
            min_x,
            self.profile.y_at_x(min_x),
            markers[1],
            color=self.color,
        )

        # plot the symmetry on a secondary axis
        sec_ax = axis.twinx()
        sec_ax.set_ylabel(self.name)

        # plot the symmetry on the secondary axis
        # add some vertical padding and/or use the minimum/maximum symmetry values
        ylim_top = max((max(self.symmetry_values) + 0.5, self.max_sym + 0.5))
        ylim_bottom = min((min(self.symmetry_values) - 0.5, self.min_sym - 0.5))
        sec_ax.set_ylim(ylim_bottom, ylim_top)
        sec_ax.plot(
            self.profile.field_x_values(self.in_field_ratio),
            self.symmetry_values,
            color=self.color,
            linestyle=self.linestyle,
        )


class SymmetryPointDifferenceQuotientMetric(SymmetryPointDifferenceMetric):
    """Symmetry as defined by IEC."""

    name = "Point Difference Quotient Symmetry"

    def __init__(
        self,
        in_field_ratio: float = 0.8,
        color="magenta",
        linestyle="--",
        max_sym_range: float = 2,
        min_sym_range: float = 0,
    ):
        super().__init__(in_field_ratio, color, linestyle, max_sym_range, min_sym_range)

    @staticmethod
    def _calc_point(lt: float, rt: float, cax: float) -> float:
        """Calculate an individual point's symmetry."""
        return 100 * max((lt / rt), (rt / lt))

    def plot(self, axis: plt.Axes, markers: (str, str) = ("x", "x")) -> None:
        super().plot(axis, markers)


class SymmetryAreaMetric(ProfileMetric):
    """The symmetry using ratios of the areas of the left and right sides of the profile."""

    name = "Symmetry (Area)"

    def __init__(
        self,
        in_field_ratio: float = 0.8,
    ):
        self.in_field_ratio = in_field_ratio

    def calculate(self) -> float:
        """Calculate the symmetry ratio of the profile using the area of the left side vs the right side."""
        _, _, width = self.profile.field_indices(in_field_ratio=self.in_field_ratio)
        area_left = np.sum(
            self.profile.field_values(self.in_field_ratio)[: math.floor(width / 2) + 1]
        )
        area_right = np.sum(
            self.profile.field_values(self.in_field_ratio)[math.ceil(width / 2) :]
        )
        return 100 * (area_left - area_right) / (area_left + area_right)

    def plot(self, axis: plt.Axes):
        """Plot the symmetry by shading the left and right areas"""
        field_values = self.profile.field_values(self.in_field_ratio)
        x_values = self.profile.field_x_values(self.in_field_ratio)
        split = math.floor(len(field_values) / 2)
        left_data = field_values[: split + 1]
        left_x = x_values[: split + 1]
        right_data = field_values[split:]
        right_x = x_values[split:]
        axis.fill_between(left_x, left_data, alpha=0.2, label="Left Area")
        axis.fill_between(right_x, right_data, alpha=0.2, label="Right Area")


class PenumbraLeftMetric(ProfileMetric):
    unit = "%"
    name = "Left Penumbra"
    side = LEFT

    def __init__(self, lower: float = 20, upper: float = 80, color="pink", ls="-."):
        self.lower = lower
        self.upper = upper
        super().__init__(color=color, linestyle=ls)

    def calculate(self) -> float:
        """Calculate the left penumbra in mm.
        We first find the edge point and then return the
        distance from the lower penumbra value to upper penumbra value.
        The trick is that wherever the field edge is, is assumed to be 50%
        height. It's okay if it's not actually (like for FFF).
        """
        left_edge = self.profile.field_edge_idx(side=self.side)
        left_edge_value = self.profile.y_at_x(left_edge)
        lower_search_value = left_edge_value * 2 * self.lower / 100
        lower_index = self.profile.x_at_y(y=lower_search_value, side=self.side)
        upper_search_value = left_edge_value * 2 * self.upper / 100
        upper_index = self.profile.x_at_y(y=upper_search_value, side=self.side)
        self.lower_index = lower_index
        self.upper_index = upper_index
        return abs(upper_index - lower_index) / self.profile.dpmm

    def plot(self, axis: plt.Axes):
        axis.vlines(
            x=[self.lower_index, self.upper_index],
            ymin=self.profile.values.min(),
            ymax=self.profile.values.max(),
            color=self.color,
            linestyle=self.linestyle,
            label=self.name,
        )


class PenumbraRightMetric(PenumbraLeftMetric):
    side = RIGHT
    name = "Right Penumbra"


class TopDistanceMetric(ProfileMetric):
    """The distance from an FFF beam's "top" to the center of the field. Similar, although
    not 100% faithful to NCS-33. The NCS report uses the middle 5cm but we use a field ratio.
    In practice, this shouldn't make a difference."""

    name = "Top Distance"
    unit = "mm"

    def __init__(self, top_region_ratio: float = 0.2, color="orange"):
        self.top_region_ratio = top_region_ratio
        super().__init__(color=color)

    def calculate(self) -> float:
        """Calculate the distance from the top to the field center. Positive means the top is to the right,
        negative means the top is to the left."""
        values = self.profile.field_values(in_field_ratio=self.top_region_ratio)
        left, right, _ = self.profile.field_indices(
            in_field_ratio=self.top_region_ratio
        )
        fit_params = np.polyfit(
            range(left, right + 1),
            values,
            deg=2,
        )

        # minimize the polynomial function
        min_f = minimize(
            lambda x: -np.polyval(
                fit_params, x
            ),  # return the negative since we're MINIMIZING and want the top value
            method="Nelder-Mead",
            x0=self.profile.center_idx,
            bounds=((left, right),),
        )
        top_idx = min_f.x[0]
        self.top_idx = top_idx
        self.top_values = np.polyval(fit_params, range(left, right + 1))
        return (top_idx - self.profile.center_idx) / self.profile.dpmm

    def plot(self, axis: plt.Axes):
        """Plot the top point and the fitted curve."""
        axis.plot(
            self.top_idx,
            self.profile.y_at_x(self.top_idx),
            "o",
            color=self.color,
            label=self.name,
        )
        left, right, _ = self.profile.field_indices(
            in_field_ratio=self.top_region_ratio
        )
        axis.plot(
            range(left, right + 1),
            self.top_values,
            color=self.color,
            linestyle=self.linestyle,
            label=self.name + " Fit",
        )


class Dmax(ProfileMetric):
    """Find the Dmax of the profile. This is a special case of the PDD metric.

    Parameters
    ----------
    window_mm
        The width of the window to use for the fit. The window will be centered around the maximum value point, which
        is used as the initial guess for the fit.
    poly_order
        The order of the polynomial to use for the fit. See `UnivariateSpline <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html>`__ for more information.
        Generally, an order between 3 and 5 is recommended.
    color
        The color of the Dmax point.
    linestyle
        The linestyle of the fit line.
    """

    name = "Dmax"
    unit = "mm"
    fit_x: np.ndarray
    fit_y: np.ndarray
    point_x: float
    point_y: float
    window_mm: float

    def __init__(
        self,
        window_mm: float = 20,
        poly_order: int = 5,
        color: str | None = None,
        linestyle: str | None = "-.",
    ):
        super().__init__(color=color, linestyle=linestyle)
        self.window_mm = window_mm
        self.poly_order = poly_order

    def calculate(self) -> float:
        """Calculate the Dmax of the profile.

        We find the maximum value of the profile and then fit a polynomial to the profile in a window around the
        maximum value. The Dmax is the x-value of the polynomial's maximum value."""
        # find the approximate depth first via the maximum value.
        # we don't use profile.x_at_y because the max could be on the left or right side
        # of the center idx. We don't know; a hack is to just index the max y-value.
        dmax_idx = np.argmax(self.profile.values)
        appr_dmax_mm = self.profile.x_values[dmax_idx]
        f, fit_x = self._spline_fit(self.window_mm, appr_dmax_mm, self.poly_order)
        # now maximize the polynomial to find dmax
        fun = minimize(
            lambda x: -f(x), bounds=((fit_x.min(), fit_x.max()),), x0=fit_x.mean()
        )
        self.fit_x = fit_x
        self.fit_y = f(fit_x)
        self.point_x = fun.x[0]
        self.point_y = -fun.fun  # negative because we're minimizing the negative above
        return self.point_x

    def _spline_fit(
        self, window_mm: float, depth_mm: float, poly_order: int
    ) -> (UnivariateSpline, np.ndarray):
        """Fit a spline to the profile of a given window at the passed depth."""
        half_window = window_mm / 2
        start, end = max(depth_mm - half_window, 0), min(
            depth_mm + half_window, self.profile.x_values.max()
        )
        if abs(start - end) <= half_window or start > end:
            raise ValueError(
                f"The PDD/Dmax metric at {depth_mm} has a window that is at or past an edge and is too small to reliably fit the data. Make the window smaller or adjust the desired depth."
            )
        fit_x = np.arange(start, end + 1, 0.1)  # interpolate the fit to 0.1mm
        f = UnivariateSpline(fit_x, self.profile.y_at_x(fit_x), k=poly_order)
        return f, fit_x

    def plot(self, axis: plt.Axes):
        """Plot the PDD point and polynomial fit."""
        axis.plot(
            self.point_x,
            self.point_y,
            "D",
            color=self.color,
            label=f"{self.name} ({self.point_x:.2f}{self.unit})",
        )
        axis.plot(
            self.fit_x,
            self.fit_y,
            color=self.color,
            linestyle=self.linestyle,
        )


class PDD(Dmax):
    """The PDD at a given depth.

    This will fit a polynomial to the profile in a window around the depth of interest and
    calculate the y-value of the polynomial at the depth of interest. This is the
    un-normalized value. We then have to normalize to the Dmax.
    The original PDD is then set as PDD/Dmax to give a true percentage.

    Parameters
    ----------
    depth_mm
        The depth at which to calculate the PDD.
    window_mm
        The width of the window to use for the fit. The window will be centered around the depth of interest.
    poly_order
        The order of the polynomial to use for the fit. See `UnivariateSpline <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html>`__ for more information.
        Generally, an order between 1 and 2 is recommended.
    normalize_to
        The value to normalize the PDD to. Either "fit" or "max". If "fit", the Dmax is calculated using the
        default Dmax metric using the ``dmax_window_mm`` and ``dmax_poly_order`` parameters. If "max", the maximum value of the profile is used.
    dmax_window_mm
        The width of the window to use for the Dmax calculation. Only used if ``normalize_to`` is "fit".
    dmax_poly_order
        The order of the polynomial to use for the Dmax calculation. Only used if ``normalize_to`` is "fit".
    color
        The color of the PDD point.
    linestyle
        The linestyle of the fit line.
    """

    unit = "%"
    fit_x: np.ndarray
    fit_y: np.ndarray
    point_x: float
    point_y: float
    window_mm: float

    @property
    def name(self):
        return f"PDD@{self.depth_mm}mm"

    def __init__(
        self,
        depth_mm: float,
        window_mm: float = 10,
        poly_order: int = 2,
        normalize_to: Literal["fit", "max"] = "fit",
        dmax_window_mm: float = 20,
        dmax_poly_order: int = 5,
        color: str | None = None,
        linestyle: str | None = "-.",
    ):
        super().__init__(
            color=color, linestyle=linestyle, window_mm=window_mm, poly_order=poly_order
        )
        self.depth_mm = depth_mm
        self.window_mm = window_mm
        self.poly_order = poly_order
        self.normalize_to = normalize_to
        self.dmax_window = dmax_window_mm
        self.dmax_poly_order = dmax_poly_order

    def calculate(self) -> float:
        """Calculate the PDD of the profile.

        This fits a polynomial to the profile in a window around the depth of interest and
        returns the y-value of the polynomial at the depth of interest."""
        f, fit_x = self._spline_fit(self.window_mm, self.depth_mm, self.poly_order)
        self.fit_x = fit_x
        self.fit_y = f(fit_x)
        self.point_x = self.depth_mm
        self.point_y = f(self.depth_mm)
        # now we have to normalize to the dmax
        if self.normalize_to == "fit":
            dmax = Dmax(window_mm=self.dmax_window, poly_order=self.dmax_poly_order)
            dmax.inject_profile(self.profile)
            dmax.calculate()
            s = self.point_y / dmax.point_y
        elif self.normalize_to == "max":
            s = self.point_y / self.profile.values.max()
        else:
            raise ValueError(
                "The PDD normalization parameter must be either 'fit' or 'max'."
            )
        return s * 100

    def plot(self, axis: plt.Axes):
        """Plot the PDD point and polynomial fit."""
        axis.plot(
            self.point_x,
            self.point_y,
            "D",
            color=self.color,
            label=f"{self.name} ({self.calculate():.2f}{self.unit})",
        )
        axis.plot(
            self.fit_x,
            self.fit_y,
            color=self.color,
            linestyle=self.linestyle,
        )
