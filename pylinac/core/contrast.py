from __future__ import annotations

import numpy as np


class Contrast:
    """Contrast calculation technique. See :ref:`visibility`"""

    MICHELSON = "Michelson"  #:
    WEBER = "Weber"  #:
    RATIO = "Ratio"  #:
    RMS = "Root Mean Square"  #:


def visibility(array: np.array, radius: float, dqe: float, algorithm: str) -> float:
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
    dqe
        Detective Quantum Efficiency. A 1:1 corollary is not available, but generally the standard deviation
        of the ROI can be used as a surrogate.
    algorithm
        The contrast method. See :class:`~pylinac.core.contrast.Contrast` for options.
    """
    c = contrast(array, algorithm)
    return c * np.sqrt(radius**2 * np.pi) / dqe


def contrast(array: np.array, algorithm: str) -> float:
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
    if algorithm == Contrast.MICHELSON:
        return michelson(array)
    elif algorithm == Contrast.WEBER:
        if array.size != 2:
            raise ValueError(
                "For Weber algorithm, the array must be exactly 2 elements. Consult the ``weber`` function for parameter details"
            )
        return weber(array[0], array[1])
    elif algorithm == Contrast.RMS:
        return rms(array)
    elif algorithm == Contrast.RATIO:
        if array.size != 2:
            raise ValueError(
                "For Ratio algorithm, the array must be exactly 2 elements. Consult the ``ratio`` function for parameter details"
            )
        return ratio(array[0], array[1])
    else:
        raise ValueError(
            f"Contrast input of {algorithm} did not match any valid options: {Contrast.__dict__.values()}"
        )


def rms(array: np.array) -> float:
    """The root-mean-square contrast. Requires values be within 0 and 1."""
    if array.min() < 0 or array.max() > 1:
        raise ValueError(
            "RMS calculations require the input array to be normalized. I.e. only values between 0 and 1."
        )
    return np.sqrt(np.mean((array - array.mean()) ** 2))


def michelson(array: np.array) -> float:
    """The Michelson contrast. Used for sinusoidal patterns. Ranges from 0 to 1."""
    l_max, l_min = array.max(), array.min()
    return (l_max - l_min) / (l_max + l_min)


def weber(feature: float, background: float) -> float:
    """The Weber contrast. Used for patterns with a small feature within a large background. Ranges from -1 to infinity"""
    return (feature - background) / background


def ratio(feature: float, reference: float) -> float:
    """The ratio of luminescence"""
    return feature / reference
