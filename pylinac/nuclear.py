from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from functools import cached_property
from pathlib import Path
from typing import Sequence, TypedDict

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from numpy.lib.stride_tricks import sliding_window_view
from scipy.optimize import curve_fit, minimize
from scipy.signal import convolve2d
from skimage.measure import block_reduce, label, regionprops
from skimage.morphology import (
    isotropic_erosion,
    remove_small_holes,
    remove_small_objects,
)
from skimage.segmentation import find_boundaries

from pylinac.core.contrast import michelson
from pylinac.core.geometry import Circle, Point, direction_to_coords
from pylinac.core.image import DicomImage, NMImageStack
from pylinac.core.mtf import MomentMTF
from pylinac.core.profile import find_peaks
from pylinac.core.roi import DiskROI, HighContrastDiskROI, RectangleROI
from pylinac.metrics.image import WeightedCentroid


@dataclass
class MaxCountRateResults:
    max_countrate: float  #:
    max_frame: int  #:
    frame_duration: float  #:
    sums: dict[int, float]  #:


class MaxCountRate:
    """Calculate the maximum countrate of a gamma camera.

    Reimplementation of the NMQC toolkit's MaxCountRate test (4.2)

    Parameters
    ----------
    path : str | Path
        The path to the DICOM file.

    See Also
    --------
    https://humanhealth.iaea.org/HHW/MedicalPhysics/NuclearMedicine/QualityAssurance/NMQC-Plugins/OperatorManual-2017-10-20.pdf
    """

    stack: NMImageStack
    frame_duration: float
    max_countrate: float
    sums: dict[int, float]

    def __init__(self, path: str | Path) -> None:
        self.stack = NMImageStack(path)

    def analyze(self, frame_duration: float = 1.0) -> float:
        """Analyze the DICOM file and return the maximum countrate.

        Parameters
        ----------
        frame_duration : float
            The duration of each frame in seconds.
        """
        self.frame_duration = frame_duration
        self.sums = {}
        for idx, img in enumerate(self.stack.frames):
            self.sums[idx] = img.array.sum() / frame_duration

    @property
    def max_countrate(self) -> float:
        """The maximum countrate in counts/second."""
        return max(self.sums.values())

    @property
    def max_frame(self) -> int:
        """The frame number that had the maximum countrate."""
        return max(self.sums, key=self.sums.get)

    @property
    def max_time(self) -> float:
        """The time of the maximum countrate."""
        return self.max_frame * self.frame_duration

    def plot(self, show: bool = True) -> None:
        """Plot the countrate over time."""
        fig, ax = plt.subplots()
        ax.plot(
            np.asarray(list(self.sums.keys())) * self.frame_duration, self.sums.values()
        )
        ax.grid(True)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Count Rate (cps)")
        # show the frames on a second x-axis
        ax2 = ax.twiny()
        ax2.set_xlabel("Frame")
        ax2.set_xlim(np.asarray(ax.get_xlim()) / self.frame_duration)
        plt.tight_layout()
        # plot the max point
        ax.plot(self.max_time, self.max_countrate, "ro")
        if show:
            plt.show()

    def results(self) -> str:
        """Return a string representation of the results."""
        return (
            f"Max countrate: {self.max_countrate:.0f} counts/second\n"
            f"Frame duration: {self.frame_duration:.2f} seconds\n"
            f"Max frame: {self.max_frame} out of {len(self.stack.frames)}\n"
        )

    def results_data(self, as_dict: bool = False) -> dict | MaxCountRateResults:
        """Return the results as a dict or dataclass."""
        d = MaxCountRateResults(
            max_countrate=self.max_countrate,
            frame_duration=self.frame_duration,
            max_frame=self.max_frame,
            sums=self.sums,
        )
        if as_dict:
            return asdict(d)
        else:
            return d


@dataclass
class PlanarUniformityResults:
    ufov_integral_uniformity: float  #:
    ufov_differential_uniformity: float  #:
    cfov_integral_uniformity: float  #:
    cfov_differential_uniformity: float  #:


@dataclass
class FOV:
    """Represents a field of view (FOV) of a gamma camera."""

    name: str
    fov: np.ndarray
    boundary_x: np.ndarray
    boundary_y: np.ndarray
    window_size: int

    @property
    def integral_uniformity(self) -> float:
        """The integral uniformity of the FOV."""
        # non-zero values only
        return integral_uniformity(self.fov[self.fov > 0])

    @cached_property
    def _differential_uniformities(
        self,
    ) -> (dict[tuple[int, int], float], dict[tuple[int, int], float]):
        """Helper to get the differential uniformity arrays of the FOV."""
        non_zero = np.where(self.fov > 0, self.fov, np.nan)
        y_view = sliding_window_view(non_zero, window_shape=self.window_size, axis=0)
        x_view = sliding_window_view(non_zero, window_shape=self.window_size, axis=1)
        y_diffs = {}
        for i in range(y_view.shape[0]):
            for j in range(y_view.shape[1]):
                # only add the differential uniformity if the integral uniformity is valid
                unif = integral_uniformity(y_view[i, j])
                if not np.isnan(unif):
                    y_diffs[(i, j)] = unif
        # same for x
        x_diffs = {}
        for i in range(x_view.shape[0]):
            for j in range(x_view.shape[1]):
                # only add the differential uniformity if the integral uniformity is valid
                unif = integral_uniformity(x_view[i, j])
                if not np.isnan(unif):
                    x_diffs[(i, j)] = unif
        return y_diffs, x_diffs

    @property
    def differential_uniformity(self) -> float:
        """The maximum differential uniformity of the FOV."""
        max_x = max(self._differential_uniformities[1].values())
        max_y = max(self._differential_uniformities[0].values())
        return max(max_x, max_y)

    @property
    def max_point(self) -> tuple[int, int]:
        """The point where the maximum value is located."""
        nan_array = np.where(self.fov == 0, np.nan, self.fov)
        local_point = np.unravel_index(np.nanargmax(nan_array), self.fov.shape)
        return int(local_point[0]), int(local_point[1])

    @property
    def min_point(self) -> tuple[int, int]:
        """The point where the minimum value is located."""
        nan_array = np.where(self.fov == 0, np.nan, self.fov)
        local_point = np.unravel_index(np.nanargmin(nan_array), self.fov.shape)
        return int(local_point[0]), int(local_point[1])

    def plot_to(self, axis: Axes, color: str) -> None:
        """Plot the FOV and it's related metrics to an axis"""
        axis.scatter(
            self.boundary_x,
            self.boundary_y,
            color=color,
            label=f"{self.name} Boundary",
            marker=".",
        )
        axis.scatter(
            self.max_point[1],
            self.max_point[0],
            color=color,
            marker="s",
            label=f"{self.name} Max",
        )
        axis.scatter(
            self.min_point[1],
            self.min_point[0],
            color=color,
            marker="x",
            label=f"{self.name} Min",
        )
        # plot the max differential uniformity 5x1 box
        max_x = max(self._differential_uniformities[1].values())
        max_y = max(self._differential_uniformities[0].values())
        if max_x > max_y:
            max_point = max(
                self._differential_uniformities[1],
                key=self._differential_uniformities[1].get,
            )
            width = self.window_size
            height = 1
        else:
            max_point = max(
                self._differential_uniformities[0],
                key=self._differential_uniformities[0].get,
            )
            width = 1
            height = self.window_size
        # offset by half pixel so we're centered about the target pixel
        rect = Rectangle(
            (max_point[1] - 0.5, max_point[0] - 0.5),
            width,
            height,
            linewidth=1,
            edgecolor=color,
            facecolor="none",
            label=f"{self.name} Max Diff. Window",
        )
        axis.add_patch(rect)
        axis.legend()


class PlanarUniformity:
    """Analyzes an image for its integral and differential uniformity."""

    stack: NMImageStack
    frame_results: dict[str, dict[str, FOV | np.ndarray]]

    def __init__(self, path: str | Path) -> None:
        self.stack = NMImageStack(path)
        self.path = Path(path)

    def analyze(
        self,
        ufov_ratio: float = 0.95,
        cfov_ratio: float = 0.75,
        window_size: int = 5,
        threshold: float = 0.75,
    ) -> None:
        """Analyze the field for the UFOV and CFOV's integral and differential uniformity.

        Parameters
        ----------
        ufov_ratio : float
            The ratio of the useful field of view (UFOV) to the detected size of the field of view.
        cfov_ratio : float
            The ratio of the central field of view (CFOV) to the UFOV. E.g. if this is 0.75 and the UFOV ratio is 0.95,
            then the CFOV will be 0.95*0.75=0.7125 of the total FOV.
        window_size : int
            The size of the window in pixels to use for the differential uniformity calculation.
        threshold : float
            The threshold to use for removing low-value pixels. This is a fraction of the mean value of the pixels
            that are > 0.
        """
        self.frame_results = {}
        for idx, frame in enumerate(self.stack.frames):
            # we need to block-reduce until the pixel size is within NEMA range
            cleaned_frame, _ = self.preprocess(frame, threshold=threshold)
            ufov_array, ufov_x, ufov_y = get_fov(cleaned_frame, ufov_ratio)
            ufov = FOV(
                name="UFOV",
                fov=ufov_array,
                boundary_x=ufov_x,
                boundary_y=ufov_y,
                window_size=window_size,
            )
            cfov_array, cfov_x, cfov_y = get_fov(cleaned_frame, cfov_ratio * ufov_ratio)
            cfov = FOV(
                name="CFOV",
                fov=cfov_array,
                boundary_x=cfov_x,
                boundary_y=cfov_y,
                window_size=window_size,
            )
            self.frame_results[str(idx + 1)] = {
                "ufov": ufov,
                "cfov": cfov,
                "binned_frame": cleaned_frame,
            }

    def results(self) -> str:
        """Return a string representation of the results."""
        s = []
        for key, result in self.frame_results.items():
            s.append(f"Frame {key}:\n")
            s.append(
                f"UFOV integral uniformity: {result['ufov'].integral_uniformity:.2f}%\n"
            )
            s.append(
                f"UFOV differential uniformity {result['ufov'].differential_uniformity:.2f}%\n"
            )
            s.append(
                f"CFOV integral uniformity: {result['cfov'].integral_uniformity:.2f}%\n"
            )
            s.append(
                f"CFOV differential uniformity {result['cfov'].differential_uniformity:.2f}%\n"
            )
            s.append("\n")
        return "".join(s)

    def results_data(self, as_dict: bool = False) -> dict:
        data = {}
        for key, result in self.frame_results.items():
            r = PlanarUniformityResults(
                ufov_integral_uniformity=result["ufov"].integral_uniformity,
                ufov_differential_uniformity=result["ufov"].differential_uniformity,
                cfov_integral_uniformity=result["cfov"].integral_uniformity,
                cfov_differential_uniformity=result["cfov"].differential_uniformity,
            )
            if as_dict:
                # make the key 1-based for ease of understanding by the user and match
                # dicom labeling
                data[f"Frame {key}"] = asdict(r)
            else:
                data[f"Frame {key}"] = r
        return data

    def plot(self, show: bool = True, cmap: str = "gray") -> (list[Figure], list[Axes]):
        """Plot each frame with the UFOV and CFOV boundaries, max points, and max differential uniformity window."""
        figs = []
        axes = []
        for key, result in self.frame_results.items():
            fig, axis = plt.subplots()
            nan_array = np.where(
                result["binned_frame"] == 0, np.nan, result["binned_frame"]
            )
            axis.imshow(
                result["binned_frame"],
                cmap=cmap,
                vmin=np.nanmin(nan_array),
                vmax=np.nanmax(nan_array),
            )
            result["ufov"].plot_to(axis, color="y")
            result["cfov"].plot_to(axis, color="r")
            axis.legend(loc="upper right")
            fig.suptitle(f"Frame {key}")
            figs.append(fig)
            axes.append(axis)
        if show:
            plt.show()
        return figs, axes

    @staticmethod
    def preprocess(frame: DicomImage, threshold: float) -> (np.ndarray, int):
        array = np.copy(frame.array)
        pixel_size = frame.metadata.PixelSpacing[0]
        bin_size = determine_binning(pixel_size)
        array = block_reduce(array, block_size=(bin_size, bin_size), func=np.sum)
        # filter image
        # kernel based on NEMA/IAEA
        # see https://www-pub.iaea.org/MTCD/Publications/PDF/Pub1394_web.pdf
        # pg 59/71
        kernel = np.array([[1, 2, 1], [2, 4, 2], [1, 2, 1]], dtype=float)
        kernel /= kernel.sum()
        array = convolve2d(array, kernel, mode="same")
        # clean the edges
        array[0, :] = 0
        array[-1, :] = 0
        array[:, 0] = 0
        array[:, -1] = 0
        # remove pixels that are <75% of mean of "meaningful" pixels
        # meaningful pixels are those > 10% of max
        # this helps remove the background
        threshold = array[array > np.max(array) * 0.10].mean() * threshold
        array[array < threshold] = 0
        # remove stray pixels
        binary_frame = array > 0
        remove_small_objects(binary_frame, min_size=2, out=binary_frame)
        remove_small_holes(binary_frame, area_threshold=2, out=binary_frame)
        array[binary_frame == 0] = 0
        return array, bin_size


def get_fov(array: np.ndarray, size: float) -> (np.ndarray, np.ndarray, np.ndarray):
    """Get the boundary of the FOV, either CFOV or UFOV.

    This will also return the boundaries of the FOV for plotting purposes.

    Parameters
    ----------
    array : np.ndarray
        The array to find the FOV of.
    size : float
        The size of the FOV as a ratio of the largest dimension (0-1).
    """
    binary_frame = array > 0
    labeled_frame, num_labels = label(binary_frame, connectivity=1, return_num=True)
    rois = regionprops(labeled_frame, intensity_image=array)
    largest_roi = max(rois, key=lambda x: x.area)
    # find the largest dimension, x or y
    longest_dim = max(largest_roi.image.shape)
    erosion = int(round((1 - size) * longest_dim))
    eroded_binary = isotropic_erosion(binary_frame, radius=erosion / 2)
    boundary = find_boundaries(eroded_binary, connectivity=1, mode="inner")
    boundary_y, boundary_x = np.nonzero(boundary)
    fov_array = np.where(eroded_binary, array, 0)
    # get the offset of the fov array from the original array
    return fov_array, boundary_x, boundary_y


def integral_uniformity(array: np.ndarray) -> float:
    """Calculate the integral uniformity of a FOV as defined by IAEA pg 7:
    https://humanhealth.iaea.org/HHW/MedicalPhysics/NuclearMedicine/QualityAssurance/NMQC-Plugins/OperatorManual-2017-10-20.pdf
    """
    return michelson(array) * 100


def determine_binning(pixel_size: float) -> int:
    """Determine the binning factor for a given pixel size.

    We want to bin the image until the pixel size is within the NEMA range of 4.48-8.32mm.
    """
    binning = 1
    while pixel_size < 4.48:
        pixel_size *= 2
        binning *= 2
    return binning


@dataclass
class CenterOfRotationResults:
    x_deviation_mm: float  #:
    y_deviation_mm: float  #:


class CenterOfRotation:
    """Analyze the center of rotation deviation of a gamma camera."""

    centroids: dict[float, Point]
    cor_x: dict[str, float | np.ndarray]
    cor_y: dict[str, float | np.ndarray]

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.stack = NMImageStack(path)

    def analyze(self) -> None:
        """Analyze the DICOM file to determine the deviation from the center of rotation."""
        # find the weighted centroids of all the images
        rot_dir = self.stack.metadata.RotationInformationSequence[0].RotationDirection
        rot_sign = -1 if rot_dir == "CW" else 1
        start_angle = self.stack.metadata.RotationInformationSequence[0].StartAngle
        step_size = self.stack.metadata.RotationInformationSequence[0].AngularStep
        centroids = {}
        for idx, frame in enumerate(self.stack.frames):
            centroid = frame.compute(WeightedCentroid())
            angle = start_angle + rot_sign * idx * step_size
            centroids[angle] = centroid

        x_values = np.radians(list(centroids.keys()))
        # the data is offset by half a pixel due to physical
        # pixel size. I.e. index 0 is really half a pixel from the edge
        half_pixel = self.stack.metadata.PixelSpacing[0] * 0.5
        # x-fit
        y_values = (
            np.asarray([p.x for p in centroids.values()])
            * self.stack.metadata.PixelSpacing[0]
            + half_pixel
        )
        params, _ = curve_fit(
            sinusoidal_fit, x_values, y_values, p0=[np.mean(y_values), 1, 1, 1]
        )
        fitted_y_data = sinusoidal_fit(
            x_values, params[0], params[1], params[2], params[3]
        )
        residuals = y_values - fitted_y_data
        self.cor_x = {
            "x_values": x_values,
            "y_values": y_values,
            "a": params[0],
            "b": params[1],
            "c": params[2],
            "phi": params[3],
            "fitted_y_values": fitted_y_data,
            "residuals": residuals,
        }
        # y-data
        y_values = (
            np.asarray([p.y for p in centroids.values()])
            * self.stack.metadata.PixelSpacing[0]
            + half_pixel
        )
        residuals = y_values - np.mean(y_values)
        self.cor_y = {
            "x_values": x_values,
            "residuals": residuals,
        }

    @property
    def x_cor_deviation_mm(self) -> float:
        """The deviation of the center of rotation from the center of the matrix in mm."""
        return np.max(np.abs(self.cor_x["residuals"]))

    @property
    def y_cor_deviation_mm(self) -> float:
        """The deviation of the center of rotation from the center of the matrix in mm."""
        return np.max(np.abs(self.cor_y["residuals"]))

    def results(self) -> str:
        """Return a string representation of the results."""
        return (
            f"Center of Rotation results for {self.path.name}\n"
            f"X-axis center of rotation deviation (mm): {self.x_cor_deviation_mm:.3f}\n"
            f"Y-axis center of rotation deviation (mm): {self.y_cor_deviation_mm:.3f}\n"
        )

    def results_data(
        self, as_dict: bool = False
    ) -> dict[str, float] | CenterOfRotationResults:
        """Return the results as a structure."""
        r = CenterOfRotationResults(
            x_deviation_mm=self.x_cor_deviation_mm,
            y_deviation_mm=self.y_cor_deviation_mm,
        )
        if as_dict:
            return asdict(r)
        else:
            return r

    def plot(self, show: bool = True) -> (list[Figure], list[Axes]):
        figs = []
        axes = []
        # plot the x-axis sine fit
        fig, ax = plt.subplots()
        ax.plot(self.cor_x["x_values"], self.cor_x["y_values"], "bo")
        ax.plot(
            self.cor_x["x_values"],
            self.cor_x["fitted_y_values"],
            "r-",
            label=f'{self.cor_x["a"]:2.2f}{self.cor_x["b"]:+2.3f}*sin({self.cor_x["c"]:2.2f}*\N{GREEK SMALL LETTER THETA}{self.cor_x["phi"]:+2.2f})',
        )
        ax.legend()
        ax.set_xlabel("Angle (radians)")
        ax.set_ylabel("Position (mm)")
        ax.grid(True)
        fig.suptitle("Sine fit (X-axis)")
        figs.append(fig)
        axes.append(ax)

        # plot the x and y-axis residuals
        for cor, axis in zip([self.cor_x, self.cor_y], ["X-axis", "Y-axis"]):
            fig, ax = plt.subplots()
            ax.plot(cor["x_values"], cor["residuals"], "bo")
            ax.set_xlabel("Angle (radians)")
            ax.set_ylabel("Residual Error (mm)")
            ax.grid(True)
            fig.suptitle(f"Residual error ({axis})")
            figs.append(fig)
            axes.append(ax)

        if show:
            plt.show()
        return figs, axes


def sinusoidal_fit(theta: float, a: float, b: float, c: float, phi: float) -> float:
    """Return a fitted sinusoidal function. See IAEA pg176, Method B (2)."""
    return a + b * np.sin(c * theta + phi)


def weighted_centroid_3d(arr: np.ndarray) -> tuple[float, float, float] | None:
    if np.sum(arr) == 0:
        return None  # Avoid division by zero

    # Get the indices of all elements in 3D
    z_indices, y_indices, x_indices = np.indices(arr.shape)

    # Calculate the sum of weights (total weight)
    total_weight = np.sum(arr)

    # Calculate the weighted sum of indices for each dimension
    x_weighted_sum = np.sum(x_indices * arr)
    y_weighted_sum = np.sum(y_indices * arr)
    z_weighted_sum = np.sum(z_indices * arr)

    # Calculate the centroid for each dimension
    centroid_x = x_weighted_sum / total_weight
    centroid_y = y_weighted_sum / total_weight
    centroid_z = z_weighted_sum / total_weight

    return centroid_x, centroid_y, centroid_z


@dataclass
class TomographicResolutionResults:
    x_fwhm: float  #:
    y_fwhm: float  #:
    z_fwhm: float  #:
    x_fwtm: float  #:
    y_fwtm: float  #:
    z_fwtm: float  #:


@dataclass
class TomographicResolutionAxisData:
    axis: str  #:
    profile_array: np.ndarray  #:
    pixel_size: float  #:

    def __post_init__(self):
        xs = np.arange(len(self.profile_array)) * self.pixel_size
        self.popt, _ = curve_fit(
            gaussian_fit,
            xs,
            self.profile_array,
            p0=[np.max(self.profile_array), np.mean(xs), self.pixel_size],
        )

    @property
    def fwhm(self) -> float:
        """The full width at half maximum of the x-profile."""
        return fwhm_from_gaussian(self.popt[2])

    @property
    def fwtm(self) -> float:
        """The full width at tenth maximum of the x-profile."""
        return fwtm_from_gaussian(self.popt[2])

    def plot(self) -> (Figure, Axes):
        fig, ax = plt.subplots()
        xs = np.arange(len(self.profile_array)) * self.pixel_size
        x_interpolated = (
            np.linspace(0, len(self.profile_array), num=len(self.profile_array) * 20)
            * self.pixel_size
        )
        fitted_y = gaussian_fit(x_interpolated, *self.popt)
        # plot raw
        ax.plot(xs, self.profile_array, "bo", label="Raw Data")
        # close limits to around the peak
        ax.set_xlim(
            (self.popt[1] - 10 * self.popt[2]), (self.popt[1] + 10 * self.popt[2])
        )
        # plot fitted gaussian
        ax.plot(x_interpolated, fitted_y, "r-", label="Gaussian Fit")
        ax.grid(True)
        ax.set_xlabel("Distance (mm)")
        ax.set_ylabel("Counts")
        fig.suptitle(f"{self.axis}-axis profile")
        return fig, ax


class TomographicResolution:
    """Analyze a tomographic resolution image for its x/y/z resolution. Based on IAEA test 4.3.4, pg 169

    Parameters
    ----------
    path : str | Path
        The path to the DICOM file.
    """

    x_axis: TomographicResolutionAxisData
    y_axis: TomographicResolutionAxisData
    z_axis: TomographicResolutionAxisData

    def __init__(self, path: str | Path) -> None:
        self.stack = NMImageStack(path)
        self.path = Path(path)

    def analyze(self) -> None:
        """Analyze the frames by finding the 3D weighted centroid and then taking a profile of the x/y/z axes
        and that position"""
        # find the weighted centroid of the stack in 3D
        array_3d = self.stack.as_3d_array()
        x, y, z = weighted_centroid_3d(array_3d)
        # the slice at the z-location will give the x/y resolution
        xy_frame = self.stack.frames[int(round(z))]
        p = xy_frame.compute(WeightedCentroid())
        # take a profile of the x/y resolution
        x_profile = xy_frame.array[int(round(p.y)), :]
        self.x_axis = TomographicResolutionAxisData(
            "X", x_profile, self.stack.metadata.PixelSpacing[0]
        )
        y_profile = xy_frame.array[:, int(round(p.x))]
        self.y_axis = TomographicResolutionAxisData(
            "Y", y_profile, self.stack.metadata.PixelSpacing[0]
        )
        # z axis. resolution is the spacing between slices.
        z_profile = array_3d[:, int(round(p.y)), int(round(p.x))]
        dpmm = abs(self.stack.metadata.SpacingBetweenSlices)
        self.z_axis = TomographicResolutionAxisData("Z", z_profile, dpmm)

    def results(self) -> str:
        """Return a string representation of the results."""
        return (
            f"Tomographic Resolution results for {self.path.name}\n"
            f"X-axis FWHM (mm): {self.x_axis.fwhm:.3f}\n"
            f"Y-axis FWHM (mm): {self.y_axis.fwhm:.3f}\n"
            f"Z-axis FWHM (mm): {self.z_axis.fwhm:.3f}\n"
            f"X-axis FWTM (mm): {self.x_axis.fwtm:.3f}\n"
            f"Y-axis FWTM (mm): {self.y_axis.fwtm:.3f}\n"
            f"Z-axis FWTM (mm): {self.z_axis.fwtm:.3f}\n"
        )

    def results_data(
        self, as_dict: bool = False
    ) -> dict[str, float] | TomographicResolutionResults:
        """Return the results as a structure."""
        r = TomographicResolutionResults(
            x_fwhm=self.x_axis.fwhm,
            y_fwhm=self.y_axis.fwhm,
            z_fwhm=self.z_axis.fwhm,
            x_fwtm=self.x_axis.fwtm,
            y_fwtm=self.y_axis.fwtm,
            z_fwtm=self.z_axis.fwtm,
        )
        if as_dict:
            return asdict(r)
        else:
            return r

    def plot(self) -> (list[Figure], list[Axes]):
        """Plot the x/y/z profiles and their gaussian fits."""
        figs = []
        axes = []
        # plot each axis: x, y, z
        for axis in (self.x_axis, self.y_axis, self.z_axis):
            fig, ax = axis.plot()
            figs.append(fig)
            axes.append(ax)
        return figs, axes


def fwhm_from_gaussian(std: float) -> float:
    """Return the FWHM of a gaussian given its standard deviation."""
    return 2 * math.sqrt(2 * math.log(2)) * std


def fwtm_from_gaussian(std: float) -> float:
    """Return the FWTM of a gaussian given its standard deviation."""
    return 2 * math.sqrt(2 * math.log(10)) * std


def gaussian_fit(
    x: np.ndarray, amplitude: float, mean: float, stddev: float
) -> np.ndarray:
    return amplitude * np.exp(-((x - mean) ** 2) / (2 * (stddev**2)))


def two_peak_gaussian_fit(
    x: np.ndarray,
    amplitude1: float,
    mean1: float,
    stddev1: float,
    amplitude2: float,
    mean2: float,
    stddev2: float,
) -> np.ndarray:
    return amplitude1 * np.exp(
        -((x - mean1) ** 2) / (2 * (stddev1**2))
    ) + amplitude2 * np.exp(-((x - mean2) ** 2) / (2 * (stddev2**2)))


class Nuclide:
    Tc99m = {
        "half_life_s": 6.0067 * 60 * 60
    }  # 6.0067 hours https://web.archive.org/web/20160804083942/http://www.nucleide.org/DDEP_WG/Nuclides/Tc-99m_tables.pdf
    Y90 = {"half_life_s": 64.1 * 60 * 60}  # 64.1 hours https://www.nndc.bnl.gov/nudat3/
    I131 = {
        "half_life_s": 8.019 * 24 * 60 * 60
    }  # 8.019 days https://www.nndc.bnl.gov/nudat3/
    Ga67 = {
        "half_life_s": 3.261 * 24 * 60 * 60
    }  # 3.261 days https://www.nndc.bnl.gov/nudat3/
    In111 = {
        "half_life_s": 2.804 * 24 * 60 * 60
    }  # 2.804 days https://www.nndc.bnl.gov/nudat3/
    Lu177 = {
        "half_life_s": 6.647 * 24 * 60 * 60
    }  # 6.647 days https://www.nndc.bnl.gov/nudat3/


@dataclass
class SimpleSensitivityResults:
    phantom_cps: float  #:
    background_cps: float  #:
    half_life_s: float  #:
    duration_s: float  #:
    decay_correction: float  #:
    sensitivity_mbq: float  #:
    sensitivity_uci: float  #:


class SimpleSensitivity:
    """The 'simple' sensitivity test as defined by IAEA 2.3.9. Equations come from the IAEA NMQC toolkit."""

    half_life_s: float
    activity_mbq: float

    def __init__(
        self, phantom_path: str | Path, background_path: str | Path | None = None
    ):
        self.phantom_path = Path(phantom_path)
        self.background_path = (
            Path(background_path) if background_path is not None else None
        )

    @property
    def phantom_cps(self) -> float:
        """The counts per second of the phantom."""
        phantom_img = DicomImage(self.phantom_path, raw_pixels=True)
        counts = phantom_img.array.sum()
        return counts / self.duration_s

    @property
    def duration_s(self) -> float:
        """The duration of the phantom image."""
        phantom_img = DicomImage(self.phantom_path, raw_pixels=True)
        return phantom_img.metadata.ActualFrameDuration / 1000

    @property
    def background_cps(self) -> float:
        """The counts per second of the background."""
        if self.background_path is None:
            return 0
        else:
            background_stack = NMImageStack(self.background_path)
            duration_s = background_stack.metadata.ActualFrameDuration / 1000
            # mean background
            avg_count = background_stack.as_3d_array().mean(axis=0).sum()
            return avg_count / duration_s

    def analyze(self, activity_mbq: float, nuclide: type[Nuclide]):
        self.half_life_s = nuclide["half_life_s"]
        self.activity_mbq = activity_mbq

    def results(self) -> str:
        """Return a string representation of the results."""
        return (
            f"Simple Sensitivity results for {self.phantom_path.name}\n"
            f"Phantom c/s: {self.phantom_cps:.0f}\n"
            f"Background c/p: {self.background_cps:.0f}\n"
            f"Half-life: {self.half_life_s:.0f}\n"
            f"Duration: {self.duration_s:.0f}\n"
            f"Decay Correction: {self.decay_correction:.3f}\n"
            f"Sensitivity (MBq): {self.sensitivity_mbq:.3f}\n"
            f"Sensitivity (uCi): {self.sensitivity_uci:.3f}\n"
        )

    def results_data(
        self, as_dict: bool = False
    ) -> dict[str, float] | SimpleSensitivityResults:
        """Return the results as a structure."""
        r = SimpleSensitivityResults(
            phantom_cps=self.phantom_cps,
            background_cps=self.background_cps,
            half_life_s=self.half_life_s,
            duration_s=self.duration_s,
            decay_correction=self.decay_correction,
            sensitivity_mbq=self.sensitivity_mbq,
            sensitivity_uci=self.sensitivity_uci,
        )
        if as_dict:
            return asdict(r)
        else:
            return r

    @property
    def decay_correction(self) -> float:
        """The decay correction factor."""
        x = np.log(2) * self.duration_s / self.half_life_s
        return 1 / x * (1 - np.exp(-x))

    @property
    def sensitivity_mbq(self) -> float:
        """The sensitivity in MBq."""
        return (
            self.phantom_cps / self.decay_correction - self.background_cps
        ) / self.activity_mbq

    @property
    def sensitivity_uci(self) -> float:
        """The sensitivity in uCi."""
        mqb_to_uci = 27.02702702702703  # 1 MBq = 27.02702702702703 uCi
        cpm = 60
        return self.sensitivity_mbq * cpm / mqb_to_uci


@dataclass
class DoubleGaussianProfile:
    axis: str  #:
    profile_array: np.ndarray  #:
    pixel_size: float  #:
    separation_mm: float  #:

    def __post_init__(self):
        xs = np.arange(len(self.profile_array)) * self.pixel_size
        peak_idxs, other = find_peaks(self.profile_array, max_number=2, threshold=0.1)
        self.popt, _ = curve_fit(
            two_peak_gaussian_fit,
            xs,
            self.profile_array,
            p0=[
                np.max(self.profile_array),
                peak_idxs[0],
                self.pixel_size,
                np.max(self.profile_array),
                peak_idxs[1],
                self.pixel_size,
            ],
        )

    @property
    def fwhm(self) -> float:
        """The full width at half maximum of the x-profile."""
        fwhm1 = fwhm_from_gaussian(self.popt[2])
        fwhm2 = fwhm_from_gaussian(self.popt[5])
        return (fwhm1 + fwhm2) / 2

    @property
    def fwtm(self) -> float:
        """The full width at tenth maximum of the x-profile."""
        fwtm1 = fwtm_from_gaussian(self.popt[2])
        fwtm2 = fwtm_from_gaussian(self.popt[5])
        return (fwtm1 + fwtm2) / 2

    @property
    def measured_pixel_size(self) -> float:
        """The measured pixel size of the image."""
        # the separation is already in physical units; we have to go back to pixels
        # to get the mm/px, which is what we're after
        separation_px = abs(self.popt[4] - self.popt[1]) / self.pixel_size
        return self.separation_mm / separation_px

    @property
    def pixel_size_difference(self) -> float:
        """The difference between the measured pixel size and the expected pixel size."""
        return (self.measured_pixel_size - self.pixel_size) / self.pixel_size * 100

    def plot(self) -> (Figure, Axes):
        fig, ax = plt.subplots()
        xs = np.arange(len(self.profile_array)) * self.pixel_size
        x_interpolated = (
            np.linspace(0, len(self.profile_array), num=len(self.profile_array) * 20)
            * self.pixel_size
        )
        fitted_y = two_peak_gaussian_fit(x_interpolated, *self.popt)
        # plot raw
        ax.plot(xs, self.profile_array, "bo", label="Raw Data")
        # plot fitted gaussian
        ax.plot(x_interpolated, fitted_y, "r-", label="Gaussian Fit")
        ax.grid(True)
        ax.legend()
        ax.set_xlabel("Distance (mm)")
        ax.set_ylabel("Counts")
        fig.suptitle(f"{self.axis}-axis profile")
        return fig, ax


@dataclass
class FourBarResolutionResults:
    x_fwhm: float
    y_fwhm: float
    x_fwtm: float
    y_fwtm: float
    x_measured_pixel_size: float
    y_measured_pixel_size: float
    x_pixel_size_difference: float
    y_pixel_size_difference: float


class FourBarResolution:
    """Spatial resolution in the X and Y direction as measured by a 'four-bar' phantom.

    Parameters
    ----------
    path : str | Path
        The path to the DICOM file.
    """

    y_prof: RectangleROI
    x_prof: RectangleROI
    y_axis: DoubleGaussianProfile
    x_axis: DoubleGaussianProfile

    def __init__(self, path: str | Path):
        self.stack = NMImageStack(path)
        self.path = Path(path)

    def analyze(self, separation_mm: float = 100, roi_width_mm: float = 10) -> None:
        """Take a vertical and horizontal profile about the center of the image.

        Fit two gaussians to find the two peaks in each direction.

        Parameters
        ----------
        separation_mm : float
            The distance between the two peaks in mm. The length of the ROI to sample
            will be 2x the separation **perpendicular to the sample lines**.
        roi_width_mm : float
            The width of the ROI. This is the width **in the direction** of the sample lines.
        """
        pixel_size = self.stack.metadata.PixelSpacing[0]
        width_px = roi_width_mm / self.stack.metadata.PixelSpacing[0]
        height_px = separation_mm * 2 / self.stack.metadata.PixelSpacing[0]
        center = Point(self.stack.metadata.Rows / 2, self.stack.metadata.Columns / 2)
        self.y_prof = RectangleROI(
            self.stack.frames[0],
            width=width_px,
            height=height_px,
            angle=0,
            dist_from_center=0,
            phantom_center=center,
        )
        v_array = self.y_prof.pixel_array.mean(axis=-1)
        self.y_axis = DoubleGaussianProfile(
            "Y/Vertical", v_array, pixel_size, separation_mm
        )
        self.x_prof = RectangleROI(
            self.stack.frames[0],
            width=height_px,
            height=width_px,
            angle=0,
            dist_from_center=0,
            phantom_center=center,
        )
        h_array = self.x_prof.pixel_array.mean(axis=0)
        self.x_axis = DoubleGaussianProfile(
            "X/Horizontal", h_array, pixel_size, separation_mm
        )

    def results(self) -> str:
        """Return a string representation of the results."""
        return (
            f"Four Bar Resolution results for {self.path.name}\n"
            f"X-axis FWHM (mm): {self.x_axis.fwhm:.3f}\n"
            f"X-axis FWTM (mm): {self.x_axis.fwtm:.3f}\n"
            f"X-axis Measured Pixel size (mm): {self.x_axis.measured_pixel_size:.3f}\n"
            f"X-axis Pixel size difference (%): {self.x_axis.pixel_size_difference:.2f}\n"
            f"Y-axis FWHM (mm): {self.y_axis.fwhm:.3f}\n"
            f"Y-axis FWTM (mm): {self.y_axis.fwtm:.3f}\n"
            f"Y-axis Measured Pixel size (mm): {self.y_axis.measured_pixel_size:.3f}\n"
            f"Y-axis Pixel size difference (%): {self.y_axis.pixel_size_difference:.2f}\n"
        )

    def results_data(self, as_dict: bool = False) -> dict | FourBarResolutionResults:
        """Return the results as a structure."""
        r = FourBarResolutionResults(
            x_fwhm=self.x_axis.fwhm,
            y_fwhm=self.y_axis.fwhm,
            x_fwtm=self.x_axis.fwtm,
            y_fwtm=self.y_axis.fwtm,
            x_measured_pixel_size=self.x_axis.measured_pixel_size,
            y_measured_pixel_size=self.y_axis.measured_pixel_size,
            x_pixel_size_difference=self.x_axis.pixel_size_difference,
            y_pixel_size_difference=self.y_axis.pixel_size_difference,
        )
        if as_dict:
            return asdict(r)
        else:
            return r

    def plot(self, show: bool = True) -> (list[Figure], list[Axes]):
        """Plot the image with the sample ROIs and the x/y profiles with gaussian fits"""
        figs = []
        axes = []
        fig, ax = plt.subplots()
        figs.append(fig)
        axes.append(ax)
        ax.imshow(self.stack.frames[0].array, cmap="gray")
        self.x_prof.plot2axes(ax, edgecolor="y")
        self.y_prof.plot2axes(ax, edgecolor="y")
        fig.suptitle(f"Four Bar Resolution for {self.path.name}")
        fig, ax = self.x_axis.plot()
        figs.append(fig)
        axes.append(ax)
        fig, ax = self.y_axis.plot()
        figs.append(fig)
        axes.append(ax)
        if show:
            plt.show()

        return figs, axes


@dataclass
class QuadrantResolutionResults:
    quadrants: dict[
        str, dict[str, float]
    ]  #:  quadrant idx: {'mtf': mtf, 'fwhm': fwhm, 'lpmm': lpmm}


class QuadrantResolution:
    """Analyze a 4-quadrant image of high-contrast line pairs to determine MTF and FWHM.

    Parameters
    ----------
    path : str | Path
        The path to the DICOM file.
    """

    rois: dict[str, HighContrastDiskROI]
    mtf: MomentMTF

    def __init__(self, path: str | Path) -> None:
        self.stack = NMImageStack(path)
        self.path = Path(path)

    def analyze(
        self,
        bar_widths: Sequence[float],
        roi_diameter_mm: float = 70,
        distance_from_center_mm: float = 130,
    ) -> None:
        """Analyze the image to determine the MTF resolution and FWHMs.

        Parameters
        ----------
        bar_widths : list
            The bar widths in mm. A line pair will be 2x this value.
        roi_diameter_mm : float
            The diameter of the ROI in mm.
        distance_from_center_mm : float
            The distance from the center of the image to the center of the ROI in mm.
        """
        if len(bar_widths) != 4:
            raise ValueError("Must have 4 bar widths")
        lpmm = 1 / (2 * np.asarray(bar_widths))
        self.rois = {}
        img_center = Point(
            self.stack.metadata.Rows / 2, self.stack.metadata.Columns / 2
        )
        angles = (45, -45, -135, 135)
        for angle, spacing in zip(angles, bar_widths):
            roi = HighContrastDiskROI(
                self.stack.frames[0],
                angle=angle,
                roi_radius=roi_diameter_mm,
                dist_from_center=distance_from_center_mm,
                phantom_center=img_center,
                contrast_threshold=0,
            )
            self.rois[spacing] = roi

        self.mtf = MomentMTF.from_high_contrast_diskset(lpmm, list(self.rois.values()))

    def results(self) -> str:
        """Return a string representation of the results."""
        s = f"Quadrant Resolution results for {self.path.name}\n"
        for quadrant, ((lpmm, mtf), fwhm) in enumerate(
            zip(self.mtf.mtfs.items(), self.mtf.fwhms.values())
        ):
            spacing = 1 / (lpmm * 2)
            s += f"Quadrant {quadrant+1}; Bar width: {spacing:.2f}mm; FWHM: {fwhm:.3f}mm; MTF: {mtf:.3f}\n"
        return s

    def results_data(self, as_dict: bool = False) -> dict | QuadrantResolutionResults:
        """Return the results as a structure."""
        r = QuadrantResolutionResults(
            quadrants={
                f"{idx+1}": {
                    "mtf": mtf,
                    "fwhm": fwhm,
                    "lpmm": lpmm,
                    "spacing": 1 / (lpmm * 2),
                }
                for idx, ((lpmm, mtf), fwhm) in enumerate(
                    zip(self.mtf.mtfs.items(), self.mtf.fwhms.values())
                )
            }
        )
        if as_dict:
            return asdict(r)
        else:
            return r

    def plot(self, show: bool = True) -> (list[Figure], list[Axes]):
        """Plot the image, the MTF, and the FWHMs."""
        figs = []
        axes = []
        fig, ax = plt.subplots()
        figs.append(fig)
        axes.append(ax)
        ax.imshow(self.stack.frames[0].array, cmap="gray")
        for idx, (spacing, roi) in enumerate(self.rois.items()):
            roi.plot2axes(ax, edgecolor="y", text=f"{idx+1}: {spacing:.2f}mm")
        fig.suptitle(f"Quadrant Resolution for {self.path.name}")
        # plot MTFs
        fig, ax = plt.subplots()
        figs.append(fig)
        axes.append(ax)
        self.mtf.plot(ax)
        # plot FWHMs
        fig, ax = plt.subplots()
        figs.append(fig)
        axes.append(ax)
        self.mtf.plot_fwhms(ax)
        if show:
            plt.show()

        return figs, axes


@dataclass
class TomographicUniformityResults:
    cfov_integral_uniformity: float
    cfov_differential_uniformity: float
    ufov_integral_uniformity: float
    ufov_differential_uniformity: float
    center_border_ratio: float
    first_frame: int
    last_frame: int


class TomographicUniformity(PlanarUniformity):
    """Evaluation of the tomographic uniformity of a SPECT image. Typically, a Jaszczak phantom or similar.
    This is similar to the Planar Uniformity test."""

    center_ratio: float
    first_frame: int
    last_frame: int
    threshold: float

    @property
    def frame_result(self) -> dict:
        """We always have a single result"""
        return self.frame_results[self.frame_key]

    @property
    def frame_key(self) -> str:
        """The key for the single frame result"""
        return f"{self.first_frame}:{self.last_frame}"

    def center_border_ratio(self, center_ratio: float, window_size: int) -> float:
        """The center-to-border ratio as defined by the NMQC toolkit.

        The center ROI is a 6cm diameter circle in the center of the phantom.
        The border ROI is the subtraction of the CFOV from the UFOV.
        """
        cleaned_frame, _ = self.preprocess(self.stack.frames[0], self.threshold)
        center_array, center_x, center_y = get_fov(cleaned_frame, size=center_ratio)
        center_fov = FOV(
            name="Center",
            fov=center_array,
            boundary_x=center_x,
            boundary_y=center_y,
            window_size=window_size,
        )
        self.frame_result["center_fov"] = center_fov
        # now get the UFOV - CFOV array
        mask = self.frame_result["cfov"].fov != 0
        ring = np.copy(self.frame_result["ufov"].fov)
        ring[mask] = np.nan
        ring[ring == 0] = np.nan
        center_array[center_array == 0] = np.nan
        return np.nanmean(center_array) / np.nanmean(ring)

    def analyze(
        self,
        first_frame: int = 0,
        last_frame: int = -1,
        ufov_ratio: float = 0.8,
        cfov_ratio: float = 0.75,
        center_ratio: float = 0.4,
        threshold: float = 0.75,
        window_size: int = 5,
    ) -> None:
        """Analyze the image to determine the uniformity.
        This will take a mean of pixel values for frames between the first and last stated frame.

        Parameters
        ----------
        first_frame : int
            The index of the first frame to analyze.
        last_frame : int
            The index of the last frame to analyze.
        ufov_ratio : float
            The ratio of the UFOV to the phantom.
        cfov_ratio : float
            The ratio of the central FOV to the UFOV.
        center_ratio : float
            The ratio of the center ROI to the phantom.
        threshold : float
            The threshold to use for the image.
        """
        # create a new single composite frame
        self.threshold = threshold
        array = self.stack.as_3d_array()
        if first_frame < 0:
            raise ValueError(
                "The first frame index is outside the array bounds. Increase the first frame index."
            )
        if last_frame < 0:
            last_frame += array.shape[0]
        if last_frame >= array.shape[0]:
            raise ValueError(
                "The last frame index is outside the array bounds. Decrease the last frame index."
            )
        if 0 < last_frame <= first_frame:
            raise ValueError(
                "The first frame index must be less than the last frame index."
            )
        # create the frankenstein frame
        new_array = array[first_frame:last_frame, :, :].mean(axis=0)
        new_frame = self.stack.frames[0]
        new_frame.array = new_array
        self.stack.frames = [
            new_frame,
        ]
        # we offset the indices by 1 since the first frame is 1, not 0
        self.first_frame = first_frame + 1
        self.last_frame = last_frame + 1
        super().analyze(
            ufov_ratio=ufov_ratio,
            threshold=threshold,
            cfov_ratio=cfov_ratio,
            window_size=window_size,
        )
        # rename from frame 1 to a dynamic name
        self.frame_results[self.frame_key] = self.frame_results.pop("1")
        self.center_ratio = self.center_border_ratio(
            center_ratio=center_ratio * ufov_ratio, window_size=window_size
        )

    def results_data(
        self, as_dict: bool = False
    ) -> dict | TomographicUniformityResults:
        """Return the results as a structure."""
        r = TomographicUniformityResults(
            cfov_integral_uniformity=self.frame_result["cfov"].integral_uniformity,
            cfov_differential_uniformity=self.frame_result[
                "cfov"
            ].differential_uniformity,
            ufov_integral_uniformity=self.frame_result["ufov"].integral_uniformity,
            ufov_differential_uniformity=self.frame_result[
                "ufov"
            ].differential_uniformity,
            center_border_ratio=self.center_ratio,
            first_frame=self.first_frame,
            last_frame=self.last_frame,
        )
        if as_dict:
            return asdict(r)
        else:
            return r

    def results(self) -> str:
        """Return a string representation of the results."""
        return (
            f"Tomographic Uniformity results for {self.path.name}\n"
            f"Frames: {self.first_frame}:{self.last_frame}\n"
            f"CFOV Integral Uniformity: {self.frame_result['cfov'].integral_uniformity:.3f}%\n"
            f"CFOV Differential Uniformity: {self.frame_result['cfov'].differential_uniformity:.3f}%\n"
            f"UFOV Integral Uniformity: {self.frame_result['ufov'].integral_uniformity:.3f}%\n"
            f"UFOV Differential Uniformity: {self.frame_result['ufov'].differential_uniformity:.3f}%\n"
            f"Center-to-Border ratio: {self.center_ratio:.3f}\n"
        )

    def plot(self, show: bool = True, cmap: str = "gray") -> (list[Figure], list[Axes]):
        """Plot the image with the sample ROIs"""
        figs, axes = super().plot(show=False, cmap=cmap)
        axis = axes[0]
        # plot the center ROI
        self.frame_result["center_fov"].plot_to(axis, color="b")
        if show:
            plt.show()


@dataclass
class TomographicROI:
    array3d: np.ndarray
    uniformity_baseline: float
    x: float
    y: float
    z: float
    radius: float
    number: str | int

    def __post_init__(self):
        self.sphere_array = (
            sample_sphere(
                self.array3d, col=self.x, row=self.y, zed=self.z, radius=self.radius
            ),
        )

    @property
    def mean_value(self) -> float:
        return float(np.nanmean(self.sphere_array))

    @property
    def min_value(self) -> float:
        return float(np.nanmin(self.sphere_array))

    @property
    def mean_contrast(self) -> float:
        return michelson(np.asarray([self.mean_value, self.uniformity_baseline])) * 100

    @property
    def max_contrast(self) -> float:
        return michelson(np.asarray([self.min_value, self.uniformity_baseline])) * 100

    def plot_to(self, axis: Axes):
        """Plot the ROI to the axis."""
        d = DiskROI(
            array=self.sphere_array,
            angle=0,
            roi_radius=self.radius,
            dist_from_center=0,
            phantom_center=Point(self.x, self.y),
        )
        d.plot2axes(axes=axis, edgecolor="r", text=str(self.number))


class TomgraphicSphere(TypedDict):
    x: float
    y: float
    z: float
    radius: float
    mean: float
    mean_contrast: float
    max_contrast: float


@dataclass
class TomographicContrastResults:
    uniformity_baseline: float
    spheres: dict[str, TomgraphicSphere]


class TomographicContrast:
    rois: dict[str, TomographicROI]

    def __init__(self, path: str | Path):
        self.stack = NMImageStack(path)
        self.path = Path(path)

    @cached_property
    def slice_data(self) -> dict[str, dict[str, float | Point]]:
        uniformities = {}
        array3d = self.stack.as_3d_array()
        global_max = array3d.max()
        for idx, frame in enumerate(self.stack.frames):
            # remove pixels that are <75% of mean of "meaningful" pixels
            # meaningful pixels are those > 10% of max
            # this helps remove the background
            arr = np.copy(frame.array)
            threshold = global_max * 0.10
            arr[arr < threshold] = 0
            binary_frame = arr > 0
            labeled_frame, num_labels = label(
                binary_frame, connectivity=1, return_num=True
            )
            if num_labels < 1:
                continue
            rois = regionprops(labeled_frame, intensity_image=arr)
            largest_roi = max(rois, key=lambda x: x.area)
            longest_dim = max(largest_roi.image.shape)
            erosion = int(round((1 - self.ufov_ratio) * longest_dim))
            eroded_binary = isotropic_erosion(binary_frame, radius=erosion / 2)
            fov_array = np.where(eroded_binary, arr, np.nan)
            u = michelson(fov_array)
            uniformities[str(idx + 1)] = {
                "fov diameter": longest_dim - erosion,
                "center": Point(x=largest_roi.centroid[1], y=largest_roi.centroid[0]),
                "area": np.count_nonzero(eroded_binary),
                "uniformity": u,
                "value": np.nanmean(fov_array),
            }
        # drop any frames where the area was significantly smaller than the median area
        # this mostly leads to dropping the edge frames, for good reason.
        median_area = np.median([v["area"] for v in uniformities.values()])
        std_area = np.std([v["area"] for v in uniformities.values()])
        return {
            k: v for k, v in uniformities.items() if v["area"] > median_area - std_area
        }

    @property
    def uniformity_frame(self) -> str:
        """The frame with the most uniformity."""
        return min(self.slice_data, key=lambda x: self.slice_data.get(x)["uniformity"])

    @property
    def uniformity_value(self) -> float:
        return self.slice_data[self.uniformity_frame]["value"]

    def analyze(
        self,
        sphere_diameters_mm: Sequence[float] = (38, 31.8, 25.4, 19.1, 15.9, 12.7),
        sphere_angles: Sequence[float] = (-10, -70, -130, -190, 110, 50),
        ufov_ratio: float = 0.8,
        search_window_px: int = 5,
        search_slices: int = 3,
    ) -> None:
        """Analyze the image to determine the contrast.

        Parameters
        ----------
        sphere_diameters_mm : list
            The diameters of the spheres in mm.
        sphere_angles : list
            The angles of the spheres in degrees.
        """
        self.ufov_ratio = ufov_ratio
        uniformities = self.slice_data
        if len(sphere_diameters_mm) != len(sphere_angles):
            raise ValueError(
                "The number of sphere diameters and angles must be the same."
            )

        # find the most non-uniform slice
        # this will usually be near the sphere slice
        # from this starting point, sample the spheres
        # optimize the location based on best contrast value
        max_uniformity_frame = max(
            uniformities, key=lambda x: uniformities[x]["uniformity"]
        )
        unif = uniformities[max_uniformity_frame]
        unif_z = int(max_uniformity_frame) - 1

        array3d = self.stack.as_3d_array()
        rois = {}
        for idx, (angle, diameter) in enumerate(
            zip(sphere_angles, sphere_diameters_mm)
        ):
            distance = math.sqrt(unif["area"] / math.pi) * 0.65
            radius = diameter / (2 * self.stack.metadata.PixelSpacing[0])
            col_x, row_y = direction_to_coords(
                unif["center"].x, unif["center"].y, distance, angle
            )
            # quicker but less accurate at the smallest ROIs
            res = minimize(
                contrast_f,
                x0=(col_x, row_y, unif_z),
                args=(array3d, radius, self.uniformity_value),
                method="Nelder-Mead",
                bounds=[
                    (col_x - search_window_px, col_x + search_window_px),
                    (row_y - search_window_px, row_y + search_window_px),
                    (unif_z - search_slices, unif_z + search_slices),
                ],
            )
            # alternatives left for future potential work
            # res = brute(contrast_f, ranges=[(col_x - search_window_px, col_x + search_window_px), (row_y - search_window_px, row_y + search_window_px), (unif_z - search_slices, unif_z + search_slices)], args=(array3d, radius, self.uniformity_value), Ns=search_window_px*2, full_output=False, finish=None)
            # res = differential_evolution(contrast_f, bounds=[(col_x - search_window_px, col_x + search_window_px), (row_y - search_window_px, row_y + search_window_px), (unif_z - search_slices, unif_z + search_slices)], args=(array3d, radius, self.uniformity_value), polish=False, x0=(col_x, row_y, unif_z), seed=1234)
            col, row, zed = res.x
            roi = TomographicROI(
                array3d=array3d,
                x=col,
                y=row,
                z=zed,
                radius=radius,
                uniformity_baseline=self.uniformity_value,
                number=idx + 1,
            )
            rois[str(idx + 1)] = roi
        self.rois = rois

    def results(self) -> str:
        """Return a string representation of the results."""
        s = f"Tomographic Contrast results for {self.path.name}\n"
        s += f"Uniformity baseline: {self.uniformity_value:.1f}\n"
        for idx, roi in self.rois.items():
            s += f"Sphere {idx}: X={roi.x:.2f},Y={roi.y:.2f},Z={roi.z:.2f} Mean: {roi.mean_value:.2f}; Mean Contrast: {roi.mean_contrast:.2f}; Max Contrast: {roi.max_contrast:.2f}\n"
        return s

    def results_data(self, as_dict: bool = False) -> dict | TomographicContrastResults:
        """Return the results as a structure."""
        r = TomographicContrastResults(
            uniformity_baseline=self.uniformity_value,
            spheres={
                idx: {
                    "x": roi.x,
                    "y": roi.y,
                    "z": roi.z,
                    "mean": roi.mean_value,
                    "mean_contrast": roi.mean_contrast,
                    "max_contrast": roi.max_contrast,
                }
                for idx, roi in self.rois.items()
            },
        )
        if as_dict:
            return asdict(r)
        else:
            return r

    def plot(self, show: bool = True) -> (list[Figure], list[Axes]):
        """Plot the uniformity frame, sphere ROI frame, and contrast vs sphere number."""
        # plot the ROIs
        roi_fig, roi_ax = plt.subplots()
        # show all ROIs on the most common sphere slice for simplicity
        median_slice = int(round(np.median([roi.z for roi in self.rois.values()])))
        roi_ax.imshow(self.stack.frames[median_slice].array, cmap="gray")
        for roi in self.rois.values():
            roi.plot_to(roi_ax)
        roi_ax.set_title(f"Sphere frame ({median_slice+1})")
        # plot the uniformity slice
        unif_fig, unif_ax = plt.subplots()
        unif_ax.imshow(
            self.stack.frames[int(self.uniformity_frame) - 1].array, cmap="gray"
        )
        un_data = self.slice_data[self.uniformity_frame]
        Circle(
            (un_data["center"].x, un_data["center"].y),
            radius=un_data["fov diameter"] / 2,
        ).plot2axes(unif_ax, edgecolor="b")
        unif_ax.set_title(f"Uniformity frame ({self.uniformity_frame})")
        # plot the contrast vs sphere number
        cont_fig, cont_ax = plt.subplots()
        cont_ax.plot(
            [int(idx) for idx in self.rois.keys()],
            [roi.mean_contrast for roi in self.rois.values()],
            color="b",
            marker="o",
            label="Mean Contrast",
        )
        cont_ax.plot(
            [int(idx) for idx in self.rois.keys()],
            [roi.max_contrast for roi in self.rois.values()],
            color="r",
            marker="o",
            label="Max Contrast",
        )
        cont_ax.set_xlabel("Sphere Number")
        cont_ax.set_ylabel("Contrast (Michelson * 100)")
        cont_ax.legend()
        cont_ax.grid(True)
        cont_ax.set_title("Contrast vs Sphere Number")
        if show:
            plt.show()
        return (roi_fig, unif_fig, cont_fig), (roi_ax, unif_ax, cont_ax)


def create_sphere_mask(
    array_shape: tuple[float, float, float],
    row: float,
    col: float,
    zed: float,
    radius: float,
) -> np.ndarray:
    """Create a mask of a sphere in an array."""
    z, y, x = np.ogrid[: array_shape[0], : array_shape[1], : array_shape[2]]
    mask = (x - col) ** 2 + (y - row) ** 2 + (z - zed) ** 2 <= radius**2
    return mask


def sample_sphere(
    array: np.ndarray, row: float, col: float, zed: float, radius: float
) -> np.ndarray:
    """Sample out a sphere from an array. Uses a mask."""
    sphere_mask = create_sphere_mask(
        array.shape, row=row, col=col, zed=zed, radius=radius
    )
    sphere_sample = np.full(array.shape, np.nan)
    sphere_sample[sphere_mask] = array[sphere_mask]
    return sphere_sample


def contrast_f(
    coords: np.ndarray, array: np.ndarray, radius: float, uniformity_baseline: float
) -> float:
    """The function to minimize for the sphere sampling."""
    col, row, zed = coords
    sample = sample_sphere(array, col=col, row=row, zed=zed, radius=radius)
    return -michelson(np.asarray([np.nanmean(sample), uniformity_baseline])) * 100
