from __future__ import annotations

import math
import warnings
from typing import Sequence

import argue
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

from .contrast import michelson
from .roi import HighContrastDiskROI


class MTF:
    """This class will calculate relative MTF"""

    def __init__(
        self,
        lp_spacings: Sequence[float],
        lp_maximums: Sequence[float],
        lp_minimums: Sequence[float],
    ):
        """
        Parameters
        ----------
        lp_spacings : sequence of floats
            These are the physical spacings per unit distance. E.g. 0.1 line pairs/mm.
        lp_maximums : sequence of floats
            These are the maximum values of the sample ROIs.
        lp_minimums : sequence of floats
            These are the minimum values of the sample ROIs.
        """
        self.spacings = lp_spacings
        self.maximums = lp_maximums
        self.minimums = lp_minimums
        if len(lp_spacings) != len(lp_maximums) != len(lp_minimums):
            raise ValueError(
                "The number of MTF spacings, maximums, and minimums must be equal."
            )
        if len(lp_spacings) < 2 or len(lp_maximums) < 2 or len(lp_minimums) < 2:
            raise ValueError(
                "The number of MTF spacings, maximums, and minimums must be greater than 1."
            )
        self.mtfs = {}
        self.norm_mtfs = {}
        for spacing, max, min in zip(lp_spacings, lp_maximums, lp_minimums):
            arr = np.array((max, min))
            self.mtfs[spacing] = michelson(arr)
        # sort according to spacings
        self.mtfs = {k: v for k, v in sorted(self.mtfs.items(), key=lambda x: x[0])}
        for key, value in self.mtfs.items():
            self.norm_mtfs[key] = (
                value / self.mtfs[lp_spacings[0]]
            )  # normalize to first region

        # check that the MTF drops monotonically by measuring the deltas between MTFs
        # if the delta is increasing it means the MTF rose on a subsequent value
        max_delta = np.max(np.diff(list(self.norm_mtfs.values())))
        if max_delta > 0:
            warnings.warn(
                "The MTF does not drop monotonically; be sure the ROIs are correctly aligned."
            )

    @argue.bounds(x=(0, 100))
    def relative_resolution(self, x: float = 50) -> float:
        """Return the line pair value at the given rMTF resolution value.

        Parameters
        ----------
        x : float
            The percentage of the rMTF to determine the line pair value. Must be between 0 and 100.
        """
        f = interp1d(
            list(self.norm_mtfs.values()),
            list(self.norm_mtfs.keys()),
            fill_value="extrapolate",
        )
        mtf = f(x / 100)
        if mtf > max(self.spacings):
            warnings.warn(
                f"MTF resolution wasn't calculated for {x}% that was asked for. The value returned is an extrapolation. Use a higher % MTF to get a non-interpolated value."
            )
        return float(mtf)

    @classmethod
    def from_high_contrast_diskset(
        cls, spacings: Sequence[float], diskset: Sequence[HighContrastDiskROI]
    ) -> MTF:
        """Construct the MTF using high contrast disks from the ROI module."""
        maximums = [roi.max for roi in diskset]
        minimums = [roi.min for roi in diskset]
        return cls(spacings, maximums, minimums)

    def plot(
        self,
        axis: plt.Axes | None = None,
        grid: bool = True,
        x_label: str = "Line pairs / mm",
        y_label: str = "Relative MTF",
        title: str = "RMTF",
        margins: float = 0.05,
        marker: str = "o",
        label: str = "rMTF",
    ) -> list[plt.Line2D]:
        """Plot the Relative MTF.

        Parameters
        ----------
        axis : None, matplotlib.Axes
            The axis to plot the MTF on. If None, will create a new figure.
        """
        if axis is None:
            fig, axis = plt.subplots()
        points = axis.plot(
            list(self.norm_mtfs.keys()),
            list(self.norm_mtfs.values()),
            marker=marker,
            label=label,
        )
        axis.margins(margins)
        axis.grid(grid)
        axis.set_xlabel(x_label)
        axis.set_ylabel(y_label)
        axis.set_title(title)

        # this whole below thing is to avoid zero division errors when plotting
        def invert(x: np.ndarray) -> np.ndarray:
            # 1/x with special treatment of x == 0
            n = np.copy(x).astype(float)
            near_zero = np.isclose(n, 0)
            n[near_zero] = np.inf
            n[~near_zero] = 1 / n[~near_zero]
            return n

        # in addition to setting the `functions`, we need to set the ticks manually.
        # if not, for some reason the ticks are bunched up on the left side and unreadable.
        # both are needed. 🤷‍♂️
        spacing_ax = axis.secondary_xaxis("top", functions=(invert, invert))
        x_ticks = axis.get_xticks()
        x2_ticks = 1 / np.clip(x_ticks, a_min=1e-2, a_max=None)
        spacing_ax.set_xticks(x2_ticks)
        spacing_ax.set_xlabel("Pair Distance (mm)")
        plt.tight_layout()
        return points


class PeakValleyMTF(MTF):
    pass


def moments_mtf(mean: float, std: float) -> float:
    """The moments-based MTF based on Hander et al 1997 Equation 8.

    See Also
    --------
    https://aapm.onlinelibrary.wiley.com/doi/epdf/10.1118/1.597928
    """
    return math.sqrt(2 * (std**2 - mean)) / mean


def moments_fwhm(width: float, mean: float, std: float) -> float:
    """The moments-based FWHM based on Hander et al 1997 Equation A8.

    Parameters
    ----------
    width : float
        The bar width in mm
    mean : float
        The mean of the ROI.
    std : float
        The standard deviation of the ROI.

    See Also
    --------
    https://aapm.onlinelibrary.wiley.com/doi/epdf/10.1118/1.597928
    """
    return 1.058 * width * math.sqrt(np.log(mean / (math.sqrt(2 * (std**2 - mean)))))


class MomentMTF:
    """A moments-based MTF. Based on the work of Hander et al 1997.

    Parameters
    ----------
    lpmms : sequence of floats
        The line pairs per mm.
    means : sequence of floats
        The means of the ROIs.
    stds : sequence of floats
        The standard deviations of the ROIs.

    See Also
    --------
    https://aapm.onlinelibrary.wiley.com/doi/epdf/10.1118/1.597928
    """

    mtfs: dict[float, float]
    fwhms: dict[float, float]

    def __init__(
        self, lpmms: Sequence[float], means: Sequence[float], stds: Sequence[float]
    ):
        self.mtfs = {}
        self.fwhms = {}
        for lpmm, mean, std in zip(lpmms, means, stds):
            bar_width = 1 / (2 * lpmm)  # lp is 2 bars
            self.mtfs[lpmm] = moments_mtf(mean, std)
            self.fwhms[lpmm] = moments_fwhm(bar_width, mean, std)

    @classmethod
    def from_high_contrast_diskset(
        cls, lpmms: Sequence[float], diskset: Sequence[HighContrastDiskROI]
    ) -> MomentMTF:
        """Construct the MTF using high contrast disks from the ROI module."""
        means = [roi.mean for roi in diskset]
        stds = [roi.std for roi in diskset]
        return cls(lpmms, means, stds)

    def plot(
        self,
        axis: plt.Axes | None = None,
        marker: str = "o",
    ) -> plt.Axes:
        """Plot the Relative MTF.

        Parameters
        ----------
        axis : None, matplotlib.Axes
            The axis to plot the MTF on. If None, will create a new figure.
        """
        if axis is None:
            fig, axis = plt.subplots()
        axis.plot(
            list(self.mtfs.keys()),
            list(self.mtfs.values()),
            marker=marker,
        )
        axis.grid(True)
        axis.set_xlabel("Line pairs / mm")
        axis.set_ylabel("MTF")
        axis.set_title("Moments-based MTF")
        spacing_ax = axis.secondary_xaxis(
            "top", functions=(lambda x: 1 / x, lambda x: 1 / x)
        )
        spacing_ax.set_xlabel("Pair Distance (mm)")
        plt.tight_layout()
        return axis

    def plot_fwhms(self, axis: plt.Axes | None = None, marker: str = "o") -> plt.Axes:
        if axis is None:
            fig, axis = plt.subplots()
        axis.plot(
            list(self.fwhms.keys()),
            list(self.fwhms.values()),
            marker=marker,
        )
        axis.grid(True)
        axis.set_xlabel("Line pairs / mm")
        axis.set_ylabel("FWHM")
        axis.set_title("Moments-based FWHM")
        spacing_ax = axis.secondary_xaxis(
            "top", functions=(lambda x: 1 / x, lambda x: 1 / x)
        )
        spacing_ax.set_xlabel("Pair Distance (mm)")
        plt.tight_layout()
        return axis
