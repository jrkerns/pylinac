from __future__ import annotations

import math
import warnings
from collections.abc import Sequence
from typing import Any, Callable, Literal

import argue
import numpy as np
import plotly.graph_objects as go
from matplotlib import pyplot as plt
from numpy import ndarray
from scipy.fft import fft, fftfreq
from scipy.interpolate import interp1d
from scipy.signal import windows

from .contrast import michelson
from .plotly_utils import add_title
from .roi import HighContrastDiskROI, RectangleROI


def _plot_invert(x: np.ndarray) -> np.ndarray:
    # avoid zero division errors when plotting
    # 1/x with special treatment of x == 0
    n = np.copy(x).astype(float)
    near_zero = np.isclose(n, 0)
    n[near_zero] = np.inf
    n[~near_zero] = 1 / n[~near_zero]
    return n


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
        cls,
        spacings: Sequence[float],
        diskset: Sequence[HighContrastDiskROI | RectangleROI],
    ) -> MTF:
        """Construct the MTF using high contrast disks from the ROI module."""
        maximums = [roi.max for roi in diskset]
        minimums = [roi.min for roi in diskset]
        return cls(spacings, maximums, minimums)

    def plotly(
        self,
        fig: go.Figure | None = None,
        x_label: str = "Line pairs / mm",
        y_label: str = "Relative MTF",
        title: str = "Relative MTF",
        name: str = "rMTF",
        **kwargs,
    ) -> go.Figure:
        """Plot the Relative MTF.

        Parameters
        ----------
        """
        fig = fig or go.Figure()
        fig.update_layout(
            showlegend=kwargs.pop("show_legend", True),
        )
        fig.add_scatter(
            x=list(self.norm_mtfs.keys()),
            y=list(self.norm_mtfs.values()),
            mode="markers+lines",
            name=name,
            **kwargs,
        )
        fig.update_layout(
            xaxis_title=x_label,
            yaxis_title=y_label,
        )
        add_title(fig, title)
        return fig

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

        # in addition to setting the `functions`, we need to set the ticks manually.
        # if not, for some reason the ticks are bunched up on the left side and unreadable.
        # both are needed. 🤷‍♂️
        spacing_ax = axis.secondary_xaxis("top", functions=(_plot_invert, _plot_invert))
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
        spacing_ax = axis.secondary_xaxis("top", functions=(_plot_invert, _plot_invert))
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
        spacing_ax = axis.secondary_xaxis("top", functions=(_plot_invert, _plot_invert))
        spacing_ax.set_xlabel("Pair Distance (mm)")
        plt.tight_layout()
        return axis


class EdgeSpreadFunctionMTF:
    """This class will calculate relative MTF from multiple edge spread function (ESF)
    The MTF is calculated for each ESF and the output is the average of all.

    Parameters
    ----------
    esf : list[ndarray]
        These are the edge spread functions (ESF). Each element of the list represents an ESF.
    sample_spacing : float | None
        This is the sample spacing in mm. If None, the frequency axis is cycles/pixel, otherwise it's converted to lp/mm. Default is None.
    padding_mode : Literal["none", "fixed", "auto"]
        This is the padding mode (adding zeros) to increase resolution. Default is "auto"

        * mode="none": array is unchanged
        * mode="fixed": pad to ``num_samples`` (must be larger than the largest array)
        * mode="auto": pad to the next power of two from the number of samples and ``num_samples``
    num_samples : int
        This is the size of the array after padding. Only applicable if padding_mode is "fixed" or "auto". Default is 1024
    windowing : Callable | None
        This is the function used to window the ESF. Default is Hann window.
    kwargs
        These are the parameters to be used when calling ``windowing``, ie windowing(kwargs)
    """

    def __init__(
        self,
        esf: list[ndarray],
        sample_spacing: float | None = None,
        padding_mode: Literal["none", "fixed", "auto"] = "auto",
        num_samples: int = 1024,
        windowing: Callable | None = windows.hann,
        **kwargs: Any,
    ):
        self.sample_spacing = sample_spacing

        # boxcar window is just a sequence of '1' so is the same as doing nothing
        windowing = windowing or windows.boxcar

        len_esf = np.unique([len(e) for e in esf])
        if padding_mode == "none":
            # validate that all arrays are the same size
            if len(len_esf) > 1:
                raise ValueError(
                    "If padding_mode='none', all ESF samples must have the same size"
                )
            num_samples = len_esf[0]
        elif padding_mode == "fixed":
            # validate that num_samples is larger the largest array
            if num_samples < max(len_esf):
                raise ValueError("num_samples must be larger than the largest array")
        elif padding_mode == "auto":
            # Select the next power of two (or num_samples if larger)
            next_power_of_two = max(2 ** np.ceil(np.log2(len_esf)))
            num_samples = int(max(next_power_of_two, num_samples))

        # frequency axis (lp/mm)
        pixel_spacing = 1 if sample_spacing is None else sample_spacing
        freq = fftfreq(num_samples, d=pixel_spacing)
        self.freq = freq[: num_samples // 2]

        # individual results
        results = [_compute_esf_mtf(e, num_samples, windowing, **kwargs) for e in esf]
        self._mtf, self._esf, self._lsf, self._lsf_windowed = (
            list(x) for x in zip(*results)
        )

        # overall mtf
        self.mtf = np.mean(np.array(self._mtf), axis=0)

    @argue.bounds(x=(0, 100))
    def relative_resolution(self, x: float = 50) -> float:
        """Return the line pair value at the given rMTF resolution value.

        Parameters
        ----------
        x : float
            The percentage of the rMTF to determine the line pair value. Must be between 0 and 100.
        """
        # invert x and mtf since interp requires xp to be increasing
        return float(np.interp(-x / 100, -self.mtf, self.freq))

    def plot(
        self,
        axis: plt.Axes | None = None,
        grid: bool = True,
        x_label: str | None = None,
        y_label: str = "Relative MTF",
        title: str = "RMTF",
        margins: float = 0.05,
        label: str = "rMTF",
    ) -> list[plt.Line2D]:
        if x_label is None:
            x_label = (
                "Cycles / sample" if self.sample_spacing is None else "Line pairs / mm"
            )

        if axis is None:
            fig, axis = plt.subplots()
        points = axis.plot(self.freq, self.mtf, label=label)
        axis.margins(margins)
        axis.grid(grid)
        axis.set_xlabel(x_label)
        axis.set_ylabel(y_label)
        axis.set_title(title)
        plt.tight_layout()
        return points

    def _plot_debug(self, plot_together: bool = False):
        n_esf = len(self._mtf)
        rows = 5
        cols = 1 if plot_together else n_esf

        fig, axis = plt.subplots(rows, cols)
        axis = (axis[np.newaxis]).transpose() if cols == 1 else axis
        for esf_idx in range(n_esf):
            col_idx = 0 if cols == 1 else esf_idx
            ax = axis[0, col_idx]
            ax.set_title("ESF")
            ax.plot(self._esf[esf_idx])

            ax = axis[1, col_idx]
            ax.set_title("LSF")
            ax.plot(self._lsf[esf_idx])

            ax = axis[2, col_idx]
            ax.set_title("LSF windowed")
            ax.plot(self._lsf_windowed[esf_idx])

            ax = axis[3, col_idx]
            ax.set_title("MTF")
            ax.plot(self.freq, self._mtf[esf_idx])

            ax = axis[4, col_idx]
            ax.set_title("MTF - diff")
            ax.plot(self.freq, self._mtf[esf_idx] - self.mtf)

        plt.tight_layout()
        plt.show()


def _compute_esf_mtf(
    esf: ndarray, num_samples: int, windowing: Callable, **kwargs
) -> tuple[ndarray, ndarray, ndarray, ndarray]:
    lsf = np.gradient(esf)
    lsf_windowed = lsf * windowing(len(esf), **kwargs)
    mtf = np.abs(fft(lsf_windowed, num_samples))
    mtf /= mtf[0]  # Normalize
    mtf = mtf[: num_samples // 2]
    return mtf, esf, lsf, lsf_windowed
