from __future__ import annotations

import math
import typing
from abc import ABC, abstractmethod
from functools import cached_property
from typing import Any

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
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
