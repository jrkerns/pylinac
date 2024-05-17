from __future__ import annotations

from typing import Iterable

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axis import Axis

from . import validators


def radial_average(arr: np.ndarray) -> np.ndarray:
    """Calculate the radial average of a 2D array about the center pixel"""
    # Determine the center of the array
    center = np.floor(np.array(arr.shape) / 2)

    # Calculate the Euclidean distance from the center for each element
    y, x = np.indices(arr.shape)
    r = np.sqrt((x - center[1]) ** 2 + (y - center[0]) ** 2)

    # Bin the elements based on their distance from the center
    r = r.astype(int)

    # Create an array to hold the sum of elements in each bin and their counts
    tbin = np.bincount(r.ravel(), arr.ravel())
    nr = np.bincount(r.ravel())

    # Avoid division by zero
    nonzero = nr != 0
    radial_mean = np.zeros(nr.shape)
    radial_mean[nonzero] = tbin[nonzero] / nr[nonzero]

    return radial_mean


def noise_power_spectrum_2d(
    pixel_size: float, rois: Iterable[np.ndarray]
) -> np.ndarray:
    """Calculate the noise power spectrum in 2D and 1D for a set of square ROIs.

    Based on ICRU 87, equation 11.1 and 11.2. Calculates the 2D FFT of each ROI, then averages the FFTs.
    FFTs are shifted so that the zero frequency is in the center and then radially averaged

    Notes
    -----
    The ROIs should be 2D and square (i.e. not a rectangle, circle, etc).
    The smallest dimension of the ROIs will be used. I.e. the ROIs can
    vary in size, but only the smallest dimension will be used. In practice,
    ROIs can sometimes vary depending on how they are extracted from an image (at least in pylinac).
    Sometimes the ROI may be 32 pixels wide and 30 pixels tall, for instance. In this case, the
    30x30 subset of the ROI will be used.


    Parameters
    ----------
    pixel_size : float
        The size of each pixel in mm.
    rois : list of ndarray
        A list of the ROIs to calculate the NPS over. Each ROI should be a 2D array.

    Returns
    -------
    nps2d : ndarray
        The 2D noise power spectrum.

    References
    ----------
    [ICRU 87] https://www.aapm.org/pubs/protected_files/ICRU/ICRU_Report_87_Radiation_Dose_and_Image-Quality_Assessment_in_Computed_Tomography_AAPM.pdf
    """
    length = min(min(roi.shape) for roi in rois)
    ffts = np.zeros((length, length, len(rois)))
    for idx, roi in enumerate(rois):
        rroi = roi[0:length, 0:length]
        b = np.abs(np.fft.fft2(rroi - np.mean(rroi))) ** 2
        s = np.fft.fftshift(b)
        ffts[:, :, idx] = s
    s = np.mean(ffts, axis=-1)
    nps2d = pixel_size**2 / length**2 * s
    return nps2d


def noise_power_spectrum_1d(spectrum_2d: np.ndarray) -> np.ndarray:
    """Calculate the 1D noise power spectrum from a 2D spectrum.

    Parameters
    ----------
    spectrum_2d : ndarray
        The 2D noise power spectrum.

    Returns
    -------
    ndarray
        The 1D noise power spectrum.
    """
    validators.double_dimension(spectrum_2d)
    return radial_average(spectrum_2d)


def average_power(nps1d: np.ndarray) -> float:
    """Average power given a 1D spectra.

    Parameters
    ----------
    nps1d : sequence of ndarray
        The 1D noise power spectra to average.

    Returns
    -------
    float
        The average noise power spectrum value
    """
    validators.single_dimension(nps1d)
    # Calculate the weighted average x position
    x_positions = np.linspace(0, 1, len(nps1d))
    return np.average(x_positions, weights=nps1d)


def max_frequency(nps1d: np.ndarray) -> float:
    """Determine the maximum frequency of the NPS given a 1D spectrum."""
    validators.single_dimension(nps1d)
    return np.argmax(nps1d) / len(nps1d)


def plot_nps1d(nps1d: np.ndarray, ax: Axis | None = None) -> Axis:
    """Plot the 1D noise power spectrum.

    This will plot the power spectrum as a function of frequency.
    We also need to remove frequencies above the Nyquist frequency.

    Parameters
    ----------
    nps1d : ndarray
        The 1D noise power spectrum.
    ax : matplotlib.axis.Axis, optional
        The axis to plot on. If None, a new figure will be created.

    Returns
    -------
    matplotlib.axis.Axis
        The axis the plot was created on.
    """
    validators.single_dimension(nps1d)
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(np.linspace(0, 1, len(nps1d)), nps1d)
    ax.set_title("1D Noise Power Spectrum")
    ax.set_xlabel("Frequency ($mm^{-1}$)")
    ax.set_ylabel("NPS / ($HU^2 mm^2$)")
    ax.grid(True)
    return ax
