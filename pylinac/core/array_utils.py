from __future__ import annotations

import numpy as np
from scipy import ndimage


def geometric_center_idx(array: np.ndarray) -> float:
    """Returns the center index and value of the profile.

    If the profile has an even number of array the centre lies between the two centre indices and the centre
    value is the average of the two centre array else the centre index and value are returned."""
    if array.ndim > 1:
        raise ValueError(
            f"Array was multidimensional. Must pass 1D array; found {array.ndim}"
        )
    return (array.shape[0] - 1) / 2.0


def geometric_center_value(array: np.ndarray) -> float:
    """Returns the center value of the profile.

    If the profile has an even number of elements the center lies between the two centre indices and the centre
    value is the average of the two center elements else the center index and value are returned."""
    if array.ndim > 1:
        raise ValueError(
            f"Array was multidimensional. Must pass 1D array; found {array.ndim}"
        )
    arr_len = array.shape[0]
    # buffer overflow can cause the below addition to give strange results
    if arr_len % 2 == 0:  # array is even and central detectors straddle CAX
        cax = (array[int(arr_len / 2)] + array[int(arr_len / 2) - 1]) / 2.0
    else:  # array is odd and we have a central detector
        cax = array[int((arr_len - 1) / 2)]
    return cax


def normalize(array: np.ndarray, value: float | None = None) -> np.ndarray:
    """Normalize an array to the passed value. If not value is passed, normalize to the maximum value"""
    if value is None:
        val = array.max()
    else:
        val = value
    array = array / val
    return array


def invert(array: np.ndarray) -> np.ndarray:
    """Invert the array. Makes the max the min and vic versa. Does NOT account for datatype"""
    return -array + array.max() + array.min()


def bit_invert(array: np.ndarray) -> np.ndarray:
    """Invert the array, ACCOUNTING for the datatype. I.e. 0 for an uint8 array goes to 255, whereas it goes to 65535 for unint16.
    I.e. this is a datatype-specific inversion."""
    try:
        return np.invert(array)
    except TypeError:
        raise ValueError(
            f"The datatype {array.dtype} could not be safely inverted. This usually means the array is a float-like datatype. Cast to an integer-like datatype first."
        )


def ground(array: np.ndarray, value: float = 0) -> np.ndarray:
    """Ground the profile. Note this will also work on profiles with negative values. I.e. this will always
    move the minimum value to 'value', regardless of whether the profile minimum was positive or negative

    Parameters
    ----------
    value
        The value to set the minimum value as.
    """
    return array - array.min() + value


def filter(array: np.ndarray, size: float = 0.05, kind: str = "median") -> np.ndarray:
    """Filter the profile.

    Parameters
    ----------
    array: np.ndarray
        The array to filter.
    size : float, int
        Size of the median filter to apply.
        If a float, the size is the ratio of the length. Must be in the range 0-1.
        E.g. if size=0.1 for a 1000-element array, the filter will be 100 elements.
        If an int, the filter is the size passed.
    kind : {'median', 'gaussian'}
        The kind of filter to apply. If gaussian, `size` is the sigma value.
    """
    if isinstance(size, float):
        if 0 < size < 1:
            size = int(round(len(array) * size))
            size = max(size, 1)
        else:
            raise ValueError("Float was passed but was not between 0 and 1")

    if kind == "median":
        filtered_array = ndimage.median_filter(array, size=size)
    elif kind == "gaussian":
        filtered_array = ndimage.gaussian_filter(array, sigma=size)
    else:
        raise ValueError(
            f"Filter type {kind} unsupported. Use one of 'median', 'gaussian'"
        )
    return filtered_array


def stretch(array: np.ndarray, min: int = 0, max: int = 1) -> np.ndarray:
    """'Stretch' the profile to the fit a new min and max value. This is a utility for grounding + normalizing.

    Parameters
    ----------
    array: numpy.ndarray
        The numpy array to stretch.
    min : number
        The new minimum of the array.
    max : number
        The new maximum value of the array.
    """
    if max <= min:
        raise ValueError(
            f"Max must be larger than min. Passed max of {max} was <= {min}"
        )
    dtype_info = get_dtype_info(array.dtype)
    if max > dtype_info.max:
        raise ValueError(
            f"Max of {max} was larger than the allowed datatype maximum of {dtype_info.max}"
        )
    if min < dtype_info.min:
        raise ValueError(
            f"Min of {min} was smaller than the allowed datatype minimum of {dtype_info.min}"
        )

    return ground(normalize(ground(array)) * (max - min), value=min)


def convert_to_dtype(array: np.ndarray, dtype: type[np.dtype]) -> np.ndarray:
    """Convert an array to another datatype, accounting for the array values.
    A normal numpy dtype conversion simply changes the datatype and leaves the values alone.
    This will convert an array and also convert the values to the same relative value of the new datatype.
    E.g. an element of value 100 on an uint8 array to be converted to an uint16 array will become ~25,690 (100/255 = 0.392 = x/65535, x = 25,690)
    """
    # original array info
    old_dtype_info = get_dtype_info(array.dtype)
    relative_values = array.astype(float) / old_dtype_info.max
    # new array info
    dtype_info = get_dtype_info(dtype)
    dtype_range = dtype_info.max - dtype_info.min
    return np.array(relative_values * dtype_range - dtype_info.max - 1, dtype=dtype)


def get_dtype_info(dtype: type[np.dtype]) -> np.iinfo | np.finfo:
    """Get the datatype of the array"""
    try:
        dtype_info = np.iinfo(dtype)
    except ValueError:
        dtype_info = np.finfo(dtype)
    return dtype_info
