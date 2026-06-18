from __future__ import annotations

import math
import warnings
from collections.abc import Sequence

import argue
import numpy as np
import plotly.graph_objects as go
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy.ndimage import sobel, gaussian_filter1d, gaussian_filter
from scipy.optimize import curve_fit
from skimage.feature import canny
from skimage.transform import hough_line, hough_line_peaks

from .contrast import michelson
from .plotly_utils import add_title
from .roi import HighContrastDiskROI


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
        diskset: Sequence[HighContrastDiskROI]
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


class EdgeMTF:
    """MTF calculation using the edge method according to IEC 62220-1-1:2015.
    
    This class implements the slanted edge method for calculating the Modulation 
    Transfer Function (MTF) as described in IEC 62220-1-1:2015. The method involves:
    
    1. Edge detection and angle determination (should be 3-5 degrees from vertical/horizontal)
    2. Edge Spread Function (ESF) extraction
    3. Differentiation to obtain Line Spread Function (LSF)
    4. Fourier Transform to calculate MTF
    
    The implementation provides robust edge detection, oversampling, and proper
    windowing to ensure accurate MTF calculation.
    
    Parameters
    ----------
    edge_data : np.ndarray
        2D array containing the edge image. The edge should be slanted at 
        approximately 3-5 degrees from vertical or horizontal.
    pixel_size : float
        Physical pixel size in mm. Used to convert spatial frequencies to mm^-1.
    edge_threshold : float, optional
        Threshold value for edge detection (0-1 range for normalized data).
        If None, will use automatic threshold (mean of image).
    edge_smoothing : float, optional
        Gaussian smoothing sigma for edge detection. Default is 1.0.
    """
    
    def __init__(
        self,
        edge_data: np.ndarray,
        pixel_size: float,
        edge_threshold: float | None = None,
        edge_smoothing: float = 1.0,
    ):
        """Initialize EdgeMTF with edge phantom image data."""
        if edge_data.ndim != 2:
            raise ValueError("Edge data must be a 2D array")
        if pixel_size <= 0:
            raise ValueError("Pixel size must be positive")
        
        self.edge_data = edge_data.astype(float)
        self.pixel_size = pixel_size
        self.edge_smoothing = edge_smoothing
        
        # Normalize edge data to 0-1 range for processing
        self.edge_data_norm = (self.edge_data - self.edge_data.min()) / (
            self.edge_data.max() - self.edge_data.min()
        )
        
        # Set threshold for edge detection
        if edge_threshold is None:
            self.edge_threshold = np.mean(self.edge_data_norm)
        else:
            if not 0 <= edge_threshold <= 1:
                raise ValueError("Edge threshold must be between 0 and 1")
            self.edge_threshold = edge_threshold
        
        # Calculate MTF
        self._calculate_mtf()
    
    def _find_edge_angle(self) -> tuple[float, bool]:
        """Find the angle of the edge using Hough Transform and determine orientation.
        
        Uses Hough line transform to detect the edge angle. This method is more robust
        than PCA because it:
        - Works with partial edges
        - Is insensitive to ROI boundary artifacts
        - Is independent of ROI aspect ratio
        - Handles noise better through voting mechanism
        
        Returns
        -------
        angle : float
            Angle of the edge in radians (with sign preserved for geometric calculations)
        is_vertical : bool
            True if edge is closer to vertical, False if closer to horizontal
        """
        # Apply Canny edge detection with adaptive parameters
        sigma = 2.0
        edges = canny(self.edge_data_norm, sigma=sigma, low_threshold=0.1, high_threshold=0.3)
        
        # If no edges found, try with lower threshold
        if not edges.any():
            edges = canny(self.edge_data_norm, sigma=sigma, low_threshold=0.05, high_threshold=0.2)
        
        # Fallback: use gradient-based edge detection
        if not edges.any():
            grad_y, grad_x = np.gradient(self.edge_data_norm)
            gradient_magnitude = np.sqrt(grad_x**2 + grad_y**2)
            threshold = np.percentile(gradient_magnitude, 90)
            edges = gradient_magnitude > threshold
        
        if not edges.any():
            raise ValueError(
                "Could not detect edge in image. Ensure edge has sufficient contrast."
            )
        
        # Hough line transform with high precision (0.1° angular resolution)
        angle_precision_deg = 0.1
        tested_angles = np.linspace(
            -np.pi / 2, np.pi / 2,
            num=int(180 / angle_precision_deg),
            endpoint=False
        )
        h, theta, d = hough_line(edges, theta=tested_angles)
        
        # Find the strongest line (peak in Hough space)
        hspace, angles, dists = hough_line_peaks(
            h, theta, d, num_peaks=1, threshold=0.3 * h.max()
        )
        
        if len(angles) == 0:
            # Fallback: find global maximum
            peak_idx = np.unravel_index(h.argmax(), h.shape)
            normal_angle_rad = theta[peak_idx[1]]
            hough_confidence = h[peak_idx] / h.max()
        else:
            normal_angle_rad = angles[0]
            hough_confidence = hspace[0] / h.max()
        
        # CRITICAL: Hough returns the angle of the NORMAL (perpendicular) to the edge
        # Convert from normal angle to edge angle by adding 90°
        edge_angle_rad = normal_angle_rad + np.pi / 2
        
        # Normalize to [-π/2, π/2] range
        if edge_angle_rad > np.pi / 2:
            edge_angle_rad -= np.pi
        elif edge_angle_rad < -np.pi / 2:
            edge_angle_rad += np.pi
        
        # Convert to degrees for diagnostics
        angle_deg = np.degrees(edge_angle_rad)
        angle_deg_abs = abs(angle_deg)
        
        # Determine orientation based on absolute angle (> 45° is vertical)
        is_vertical = angle_deg_abs > 45
        
        # Calculate edge strength: mean gradient magnitude at detected edge
        grad_y, grad_x = np.gradient(self.edge_data_norm)
        gradient_magnitude = gaussian_filter(
            np.sqrt(grad_x**2 + grad_y**2), sigma=1.0
        )
        edge_strength = np.mean(gradient_magnitude[edges])
        
        # Store diagnostic information
        self.edge_points_count = int(np.sum(edges))
        self.hough_confidence = float(hough_confidence)
        self.edge_strength = float(edge_strength)
        self.edge_angle_deg = float(angle_deg_abs)  # Absolute value for display
        self.angle_detection_method = "Hough Transform"
        
        # Check if angle is in acceptable range with orientation-specific thresholds
        if is_vertical:
            # Vertical edge: optimal range is 85-87°
            if angle_deg_abs < 85 or angle_deg_abs > 87:
                warnings.warn(
                    f"Vertical edge angle ({angle_deg_abs:.1f}°) is outside "
                    f"optimal range (85-87°). Results may be less accurate."
                )
        else:
            # Horizontal edge: optimal range is 3-5°
            if angle_deg_abs < 3 or angle_deg_abs > 5:
                warnings.warn(
                    f"Horizontal edge angle ({angle_deg_abs:.1f}°) is outside "
                    f"optimal range (3-5°). Results may be less accurate."
                )
        
        # Warn if Hough confidence is low
        if hough_confidence < 0.3:
            warnings.warn(
                f"Low Hough confidence ({hough_confidence:.2f}). "
                f"Edge may be poorly defined or noisy."
            )
        
        # Return angle with sign preserved for geometric calculations
        return edge_angle_rad, is_vertical
    
    def _extract_esf(self, angle: float, is_vertical: bool) -> tuple[np.ndarray, np.ndarray]:
        """Extract Edge Spread Function by projecting perpendicular to edge.
        
        Parameters
        ----------
        angle : float
            Angle of the edge in radians
        is_vertical : bool
            Whether the edge is primarily vertical
            
        Returns
        -------
        positions : np.ndarray
            Positions along the edge (in pixels, oversampled)
        esf : np.ndarray
            Edge Spread Function values
        """
        rows, cols = self.edge_data_norm.shape
        
        # Create coordinate grids
        y, x = np.mgrid[0:rows, 0:cols]
        
        # Calculate perpendicular distance to edge for each pixel
        # For a line at angle theta passing through center:
        # perpendicular distance = x*sin(theta) - y*cos(theta) + offset
        center_y, center_x = rows / 2, cols / 2
        perpendicular_dist = (x - center_x) * np.sin(angle) - (y - center_y) * np.cos(angle)
        
        # Flatten arrays
        distances = perpendicular_dist.flatten()
        intensities = self.edge_data_norm.flatten()
        
        # Create oversampled bins (4x oversampling as per IEC 62220-1-1:2015)
        oversampling_factor = 4
        min_dist, max_dist = distances.min(), distances.max()
        n_bins = int((max_dist - min_dist) * oversampling_factor)
        
        bin_edges = np.linspace(min_dist, max_dist, n_bins + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Bin the data to create ESF
        esf = np.zeros(n_bins)
        counts = np.zeros(n_bins)
        
        for i in range(n_bins):
            mask = (distances >= bin_edges[i]) & (distances < bin_edges[i + 1])
            if mask.any():
                esf[i] = np.mean(intensities[mask])
                counts[i] = np.sum(mask)
        
        # Remove bins with no data
        valid_bins = counts > 0
        if valid_bins.sum() < 10:
            raise ValueError(
                "Insufficient data for ESF calculation. Check edge image quality."
            )
        
        esf = esf[valid_bins]
        positions = bin_centers[valid_bins]
        
        # Sort by position
        sort_idx = np.argsort(positions)
        positions = positions[sort_idx]
        esf = esf[sort_idx]
        
        return positions, esf
    
    def _calculate_lsf(self, esf: np.ndarray) -> np.ndarray:
        """Calculate Line Spread Function by differentiating ESF.
        
        Parameters
        ----------
        esf : np.ndarray
            Edge Spread Function
            
        Returns
        -------
        lsf : np.ndarray
            Line Spread Function
        """
        # Differentiate ESF to get LSF
        # Use central differences for better accuracy
        lsf = np.gradient(esf)
        
        # Apply Hamming window to reduce ringing artifacts (IEC recommendation)
        window = np.hamming(len(lsf))
        lsf_windowed = lsf * window
        
        return lsf_windowed
    
    def _calculate_mtf(self):
        """Calculate MTF from edge data following IEC 62220-1-1:2015."""
        # Step 1: Find edge angle
        angle, is_vertical = self._find_edge_angle()
        self.edge_angle = angle
        self.is_vertical = is_vertical
        
        # Step 2: Extract ESF
        positions, esf = self._extract_esf(angle, is_vertical)
        self.esf_positions = positions
        self.esf = esf
        
        # Step 3: Calculate LSF
        lsf = self._calculate_lsf(esf)
        self.lsf = lsf
        
        # Step 4: Calculate MTF via FFT
        # Zero-pad to increase frequency resolution
        n_fft = 2 ** int(np.ceil(np.log2(len(lsf) * 2)))
        lsf_padded = np.zeros(n_fft)
        lsf_padded[:len(lsf)] = lsf
        
        # Take FFT
        mtf_complex = np.fft.fft(lsf_padded)
        mtf = np.abs(mtf_complex)
        
        # Normalize so MTF(0) = 1
        mtf = mtf / mtf[0] if mtf[0] != 0 else mtf
        
        # Calculate frequency axis
        # Pixel spacing in position array
        pixel_spacing = np.mean(np.diff(positions))
        # Nyquist frequency in cycles per mm
        nyquist_freq = 1 / (2 * self.pixel_size * pixel_spacing)
        # Frequency array
        freqs = np.fft.fftfreq(n_fft, d=pixel_spacing * self.pixel_size)
        
        # Take only positive frequencies up to Nyquist
        positive_freqs = freqs >= 0
        self.frequencies = freqs[positive_freqs]
        self.mtf_values = mtf[positive_freqs]
        
        # Limit to Nyquist frequency
        nyquist_mask = self.frequencies <= nyquist_freq
        self.frequencies = self.frequencies[nyquist_mask]
        self.mtf_values = self.mtf_values[nyquist_mask]
    
    @argue.bounds(percent=(0, 100))
    def spatial_resolution(self, percent: float = 50) -> float:
        """Return the spatial frequency at the given MTF percentage.
        
        Parameters
        ----------
        percent : float
            The MTF percentage (0-100) at which to find the spatial frequency.
            Common values are 50 (MTF50) and 10 (MTF10).
            
        Returns
        -------
        frequency : float
            The spatial frequency in cycles/mm at the given MTF percentage.
        """
        target_mtf = percent / 100
        
        # Find where MTF crosses the target value
        if target_mtf > self.mtf_values[0]:
            warnings.warn(
                f"Target MTF {percent}% is higher than MTF(0). Returning 0."
            )
            return 0.0
        
        if target_mtf < self.mtf_values[-1]:
            warnings.warn(
                f"MTF does not reach {percent}% within the measured frequency range. "
                f"Result is extrapolated."
            )
        
        # Interpolate to find frequency at target MTF
        # MTF decreases with frequency, so we need to reverse for interpolation
        f = interp1d(
            self.mtf_values[::-1],
            self.frequencies[::-1],
            kind='linear',
            fill_value='extrapolate'
        )
        
        frequency = f(target_mtf)
        return float(frequency)
    
    def plotly(
        self,
        fig: go.Figure | None = None,
        x_label: str = "Spatial Frequency (cycles/mm)",
        y_label: str = "MTF",
        title: str = "Edge-based MTF (IEC 62220-1-1:2015)",
        name: str = "MTF",
        **kwargs,
    ) -> go.Figure:
        """Plot the MTF using plotly.
        
        Parameters
        ----------
        fig : go.Figure, optional
            Existing figure to add trace to. If None, creates new figure.
        x_label : str
            Label for x-axis
        y_label : str
            Label for y-axis
        title : str
            Plot title
        name : str
            Name for the trace
        **kwargs
            Additional arguments passed to go.Scatter
            
        Returns
        -------
        fig : go.Figure
            The plotly figure
        """
        fig = fig or go.Figure()
        fig.update_layout(
            showlegend=kwargs.pop("show_legend", True),
        )
        fig.add_scatter(
            x=self.frequencies,
            y=self.mtf_values,
            mode="lines",
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
        x_label: str = "Spatial Frequency (cycles/mm)",
        y_label: str = "MTF",
        title: str = "Edge-based MTF (IEC 62220-1-1:2015)",
        label: str = "MTF",
    ) -> plt.Line2D:
        """Plot the MTF using matplotlib.
        
        Parameters
        ----------
        axis : plt.Axes, optional
            Matplotlib axis to plot on. If None, creates new figure.
        grid : bool
            Whether to show grid
        x_label : str
            Label for x-axis
        y_label : str
            Label for y-axis
        title : str
            Plot title
        label : str
            Label for the line
            
        Returns
        -------
        line : plt.Line2D
            The plotted line object
        """
        if axis is None:
            fig, axis = plt.subplots()
        
        line = axis.plot(self.frequencies, self.mtf_values, label=label)[0]
        axis.grid(grid)
        axis.set_xlabel(x_label)
        axis.set_ylabel(y_label)
        axis.set_title(title)
        axis.set_ylim([0, 1.05])
        axis.legend()
        plt.tight_layout()
        
        return line
    
    def plot_esf(
        self,
        axis: plt.Axes | None = None,
        grid: bool = True,
    ) -> plt.Line2D:
        """Plot the Edge Spread Function.
        
        Parameters
        ----------
        axis : plt.Axes, optional
            Matplotlib axis to plot on. If None, creates new figure.
        grid : bool
            Whether to show grid
            
        Returns
        -------
        line : plt.Line2D
            The plotted line object
        """
        if axis is None:
            fig, axis = plt.subplots()
        
        line = axis.plot(self.esf_positions, self.esf, label="ESF")[0]
        axis.grid(grid)
        axis.set_xlabel("Position (pixels)")
        axis.set_ylabel("Normalized Intensity")
        axis.set_title("Edge Spread Function")
        axis.legend()
        plt.tight_layout()
        
        return line
    
    def plot_lsf(
        self,
        axis: plt.Axes | None = None,
        grid: bool = True,
    ) -> plt.Line2D:
        """Plot the Line Spread Function.
        
        Parameters
        ----------
        axis : plt.Axes, optional
            Matplotlib axis to plot on. If None, creates new figure.
        grid : bool
            Whether to show grid
            
        Returns
        -------
        line : plt.Line2D
            The plotted line object
        """
        if axis is None:
            fig, axis = plt.subplots()
        
        line = axis.plot(self.lsf, label="LSF")[0]
        axis.grid(grid)
        axis.set_xlabel("Position (pixels)")
        axis.set_ylabel("Amplitude")
        axis.set_title("Line Spread Function")
        axis.legend()
        plt.tight_layout()

        return line


class BeadMTF:
    """PSF-based MTF derived from a sub-millimetre point-source bead.

    The 2D PSF patch (background-subtracted) is radially averaged to a 1D
    profile, Hanning-windowed, then Fourier-transformed to obtain the MTF.
    An optional sinc correction removes the blurring contribution of the
    finite bead diameter.
    """

    BEAD_DIAMETER_MM: float = 0.18

    def __init__(
        self,
        psf_patch: np.ndarray,
        pixel_size_mm: float,
        apply_bead_correction: bool = True,
    ):
        self.pixel_size = pixel_size_mm
        self.apply_bead_correction = apply_bead_correction
        self._calculate(psf_patch)

    def _calculate(self, psf: np.ndarray) -> None:
        psf_pos = np.clip(psf, 0.0, None)
        total = psf_pos.sum()
        if total <= 0:
            raise ValueError(
                "PSF patch has no positive values; bead localisation may have failed."
            )
        y_idx, x_idx = np.indices(psf.shape)
        cx = float((x_idx * psf_pos).sum() / total)
        cy = float((y_idx * psf_pos).sum() / total)

        r_px = np.sqrt((x_idx - cx) ** 2 + (y_idx - cy) ** 2)
        n_bins = int(np.floor(r_px.max()))

        radial_psf = np.zeros(n_bins)
        for i in range(n_bins):
            mask = (r_px >= i) & (r_px < i + 1)
            if mask.any():
                radial_psf[i] = psf[mask].mean()

        if radial_psf.max() > 0:
            radial_psf /= radial_psf.max()

        # Mirror the radial PSF to build a symmetric 1D signal so that the peak
        # sits at the centre of the Hanning window (weight = 1.0) rather than at
        # the edge (weight = 0.0), which would otherwise destroy the main signal.
        symmetric = np.concatenate([radial_psf[1:][::-1], radial_psf])
        window = np.hanning(len(symmetric))
        lsf_w = symmetric * window

        n = len(lsf_w)
        n_fft = 2 ** int(np.ceil(np.log2(max(n * 8, 128))))
        lsf_padded = np.zeros(n_fft)
        lsf_padded[:n] = lsf_w

        spectrum = np.fft.rfft(lsf_padded)
        mtf = np.abs(spectrum)
        freqs = np.fft.rfftfreq(n_fft, d=self.pixel_size)

        if mtf[0] > 0:
            mtf /= mtf[0]

        if self.apply_bead_correction:
            sinc_vals = np.sinc(self.BEAD_DIAMETER_MM * freqs)
            valid = np.abs(sinc_vals) > 0.05
            mtf[valid] /= sinc_vals[valid]
            mtf[~valid] = 0.0
            if mtf[0] > 0:
                mtf /= mtf[0]

        mtf = np.clip(mtf, 0.0, 1.0)

        if np.any(np.diff(mtf) > 0.01):
            warnings.warn(
                "Bead MTF is not monotonically decreasing; bead localisation "
                "or PSF extraction may be imprecise.",
                UserWarning,
            )

        self.frequencies = freqs
        self.mtf_values = mtf

    def relative_resolution(self, x: float) -> float:
        """Spatial frequency (lp/mm) where the MTF equals x% of its zero-frequency value."""
        target = x / 100.0
        f = interp1d(
            self.mtf_values[::-1],
            self.frequencies[::-1],
            kind="linear",
            bounds_error=False,
            fill_value=(self.frequencies[-1], self.frequencies[0]),
        )
        return float(f(target))

    def plotly(
        self,
        fig: go.Figure | None = None,
        x_label: str = "Spatial Frequency (lp/mm)",
        y_label: str = "MTF",
        title: str = "Bead PSF MTF",
        name: str = "MTF",
        **kwargs,
    ) -> go.Figure:
        fig = fig or go.Figure()
        fig.update_layout(showlegend=kwargs.pop("show_legend", True))
        fig.add_scatter(
            x=self.frequencies.tolist(),
            y=self.mtf_values.tolist(),
            mode="lines",
            name=name,
        )
        fig.update_layout(xaxis_title=x_label, yaxis_title=y_label)
        add_title(fig, title)
        return fig

    def plot(
        self,
        axis: plt.Axes | None = None,
        grid: bool = True,
        x_label: str = "Spatial Frequency (lp/mm)",
        y_label: str = "MTF",
        title: str = "Bead PSF MTF",
        label: str = "MTF",
    ) -> plt.Line2D:
        if axis is None:
            _, axis = plt.subplots()
        (line,) = axis.plot(self.frequencies, self.mtf_values, label=label)
        axis.grid(grid)
        axis.set_xlabel(x_label)
        axis.set_ylabel(y_label)
        axis.set_title(title)
        axis.set_ylim([0, 1.05])
        axis.legend()
        plt.tight_layout()
        return line

