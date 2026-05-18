from __future__ import annotations

from collections.abc import Iterable

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
    image_array: np.ndarray, pixel_size: float, big_roi_size: int, small_roi_size: int
) -> tuple[np.ndarray, float]:
    """Calculate the noise power spectrum in 2D and 1D for a set of square ROIs.

    Based on ICRU 87, equation 11.1 and 11.2, and on IEC 62220-1-1:2015 described method. For a given image,
    it selects a central big ROI to apply this function, then it creates ROIs of a given size l within the big ROI,
    moving it l/2 and calculates the 2D FFT of each ROI, then averages the FFTs.
    FFTs are shifted so that the zero frequency is in the center and then radially averaged

    Notes
    -----
    The ROIs should be 2D and square (i.e. not a rectangle, circle, etc).

    Parameters
    ----------
    image_array : ndarray
        Image to apply NPS to.
    pixel_size : float
        The size of each pixel in mm.
    big_roi_size : float
        Central large ROI size in pixels.
    small_roi_size: float
        Small ROIs size in pixels.

    Returns
    -------
    nps2d : ndarray
        The 2D noise power spectrum.
    mean_pv : float
        The mean pixel value of the large ROI. Used for NNPS calculations.

    References
    ----------
    [ICRU 87] https://www.aapm.org/pubs/protected_files/ICRU/ICRU_Report_87_Radiation_Dose_and_Image-Quality_Assessment_in_Computed_Tomography_AAPM.pdf
    """

    def extract_rois_from_image(image_array, big_roi_size, small_roi_size):
        if image_array.ndim != 2:
            raise ValueError("Input image must be a 2D array.")

        img_height, img_width = image_array.shape

        # 1. Define the "big central ROI". Calculate the top-left corner for the centered big ROI
        start_y = (img_height - big_roi_size) // 2
        start_x = (img_width - big_roi_size) // 2

        # Extract the big ROI from the original image
        big_roi = image_array[start_y : start_y + big_roi_size, start_x : start_x + big_roi_size]

        # 2. Divide the big ROI into smaller, overlapping ROIs
        rois_list = []

        for y in range(0, big_roi_size - small_roi_size + 1, int(small_roi_size / 2)):
            for x in range(0, big_roi_size - small_roi_size + 1, int(small_roi_size / 2)):
                small_roi = big_roi[y : y + small_roi_size, x : x + small_roi_size]
                rois_list.append(small_roi)

        return rois_list, big_roi

    rois, big_roi = extract_rois_from_image(image_array, big_roi_size, small_roi_size)

    ffts = np.zeros((small_roi_size, small_roi_size, len(rois)))
    for idx, roi in enumerate(rois):
        power_spectrum = np.abs(np.fft.fft2(roi - np.mean(roi))) ** 2
        shifted_spectrum = np.fft.fftshift(power_spectrum)
        ffts[:, :, idx] = shifted_spectrum
    s = np.mean(ffts, axis=-1)
    nps2d = pixel_size**2 / small_roi_size**2 * s
    mean_pv = np.mean(big_roi) # For NNPS calculation
    return nps2d, mean_pv


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
