from __future__ import annotations

import numpy as np

from ..core.utilities import OptionListMixin
from . import validators


class Contrast(OptionListMixin):
    """Contrast calculation technique. See :ref:`visibility`"""

    MICHELSON = "Michelson"  #:
    WEBER = "Weber"  #:
    RATIO = "Ratio"  #:
    RMS = "Root Mean Square"  #:
    DIFFERENCE = "Difference"  #:


def visibility(array: np.ndarray, radius: float, std: float, algorithm: str) -> float:
    """The visual perception of CNR. Uses the model from A Rose: https://www.osapublishing.org/josa/abstract.cfm?uri=josa-38-2-196.
    See also here: https://howradiologyworks.com/x-ray-cnr/.
    Finally, a review paper here: http://xrm.phys.northwestern.edu/research/pdf_papers/1999/burgess_josaa_1999.pdf
    Importantly, the Rose model is not applicable for high-contrast use cases.

    This uses the ``contrast`` function under the hood. Consult before using.

    Parameters
    ----------
    array
        The numpy array of the contrast ROI or a 2-element array containing the individual
        inputs. See ``contrast`` for more.
    radius
        The radius of the contrast ROI
    std
        Standard deviation of the array. This can sometimes be obtained from another ROI, so it
        is a separate parameter.
    algorithm
        The contrast method. See :class:`~pylinac.core.contrast.Contrast` for options.
    """
    c = contrast(array, algorithm)
    return c * np.sqrt(radius**2 * np.pi) / std


def contrast(array: np.ndarray, algorithm: str) -> float:
    """Generic contrast function. Different algorithms have different inputs, so caution is advised.
    When possible, the exact contrast function is preferred.

    For Michelson and RMS algorithms, the input array can be any ordinary numpy array.
    For Weber and Ratio algorithms, the array is assumed to be a 2-element array.

    Parameters
    ----------
    array
        The numpy array of the ROI or 2-element input array. This is used in combination with the
        method.
    algorithm
        The contrast method. See :class:`~pylinac.core.contrast.Contrast` for options.
    """
    algorithm = algorithm.lower()
    if algorithm == Contrast.MICHELSON.lower():
        return michelson(array)
    elif algorithm == Contrast.WEBER.lower():
        if array.size != 2:
            raise ValueError(
                "For Weber algorithm, the array must be exactly 2 elements. Consult the ``weber`` function for parameter details"
            )
        return weber(array[0], array[1])
    elif algorithm == Contrast.RMS.lower():
        return rms(array)
    elif algorithm == Contrast.RATIO.lower():
        if array.size != 2:
            raise ValueError(
                "For Ratio algorithm, the array must be exactly 2 elements. Consult the ``ratio`` function for parameter details"
            )
        return ratio(array[0], array[1])
    elif algorithm == Contrast.DIFFERENCE.lower():
        if array.size != 2:
            raise ValueError(
                "For Difference algorithm, the array must be exactly 2 elements. Consult the ``difference`` function for parameter details"
            )
        return difference(array[0], array[1])
    else:
        raise ValueError(
            f"Contrast input of {algorithm} did not match any valid options: {Contrast.__dict__.values()}"
        )


def rms(array: np.ndarray) -> float:
    """The root-mean-square contrast. Requires values be within 0 and 1."""
    if array.min() < 0 or array.max() > 1:
        raise ValueError(
            "RMS calculations require the input array to be normalized. I.e. only values between 0 and 1."
        )
    return np.sqrt(np.mean((array - array.mean()) ** 2))


def difference(feature: float, background: float) -> float:
    """The simple absolute difference between the feature ROI and background ROI.
    This can be useful if the default CNR formula is desired (since pylinac CNR is based
    on the contrast algorithm chosen.

    .. seealso::

        https://en.wikipedia.org/wiki/Contrast-to-noise_ratio
    """
    return abs(feature - background)


def michelson(array: np.ndarray) -> float:
    """The Michelson contrast. Used for sinusoidal patterns. Ranges from 0 to 1.

    .. seealso::

        https://en.wikipedia.org/wiki/Contrast_(vision)#Michelson_contrast
    """
    l_max, l_min = np.nanmax(array), np.nanmin(array)
    return (l_max - l_min) / (l_max + l_min)


def weber(feature: float, background: float) -> float:
    """The Weber contrast. Used for patterns with a small feature within a large background. Ranges from 0 to infinity.

    For backwards compatibility with previous versions, the absolute difference is used, making the range 0 to infinity vs -1 to infinity.

    .. seealso::

        https://en.wikipedia.org/wiki/Contrast_(vision)#Weber_contrast

    .. danger::

        The default definition does not use the absolute value. We only use it here for backwards compatibility.
    """
    return abs(feature - background) / background


def ratio(feature: float, reference: float) -> float:
    """The ratio of luminescence"""
    return feature / reference


def power_spectrum_1d(array: np.ndarray) -> np.ndarray:
    """Get the 1D noise power spectrum from a 2D numpy array. ChatGPT-made.

    Parameters
    ----------
    array: np.ndarray
        The 2D numpy array.

    Returns
    -------
    array: np.ndarray
        The 1D power spectrum as an array. This spectrum is radial averaged.

    References
    ----------
    https://bertvandenbroucke.netlify.app/2019/05/24/computing-a-power-spectrum-in-python/
    https://chat.openai.com/share/7d138eb7-bb20-428f-ad61-7c63f2ec435e
    """
    validators.double_dimension(array)
    # Apply Fourier Transform
    f_transform = np.fft.fft2(array)
    f_shift = np.fft.fftshift(f_transform)

    # Calculate 2D power spectrum
    power_spectrum_2d = np.abs(f_shift) ** 2

    # Convert to 1D power spectrum by radial averaging
    y, x = np.indices(power_spectrum_2d.shape)
    center = np.array([(x.max() - x.min()) / 2.0, (y.max() - y.min()) / 2.0])
    r = np.sqrt((x - center[0]) ** 2 + (y - center[1]) ** 2)
    r = r.astype(int)

    radial_spectrum = np.bincount(r.ravel(), power_spectrum_2d.ravel()) / np.bincount(
        r.ravel()
    )
    return radial_spectrum
