
"""Common functions module."""

import numpy as np


def _datacheck_peakdetect(values, x_data):

    if x_data is None:
        x_data = list(range(len(values)))

    if len(values) != len(x_data):
        raise ValueError

    # needs to be a numpy array
    values = np.array(values)
    if np.ndim(values) != 1:
        raise IndexError("values passed was not a 1D array")
    x_data = np.array(x_data)
    return values, x_data

def peak_detect(y, x=None, threshold=0, min_peak_width=10, max_num_peaks=None, exclude_lt_edge=0.0,
                exclude_rt_edge=0.0, find_min_instead=False):
    """Find the peaks or valleys of a 1D signal.

    Uses the difference (np.diff) in signal to find peaks. Current limitations include:
        1) Only for use in 1-D data; 2D may be possible with the gradient function.
        2) Will not detect peaks at the very edge of array (i.e. 0 or -1 index)

    Parameters
    ----------
    y : array-like
        1D y-data of signal.
    x : array-like
        1D x-data of signal. If None, will create a uniform range the length of the y-data.
    threshold : int, float
        The value the peak must be above to be considered a peak. This removes "peaks"
        that are in a low-value region.
        If passed an int, the actual value is the threshold.
        E.g. when passed 15, any peak less with a value <15 is removed.
        If passed a float, it will threshold as a percent. Must be between 0 and 1.
        E.g. when passed 0.4, any peak <40% of the maximum value will be removed.
    min_peak_width : int, float
        If passed an int, parameter is the number of elements apart a peak must be from neighboring peaks.
        If passed a float, must be between 0 and 1 and represents the ratio of the profile to exclude.
        E.g. if passed 0.05 with a 1000-element profile, the minimum peak width will be 0.05*1000 = 50 elements.
    max_num_peaks : int
        Specify up to how many peaks will be returned. E.g. if 3 is passed in and 5 peaks are found, only the 3 largest
        peaks will be returned.
    find_min_instead : bool
        If False (default), peaks will be returned.
        If True, valleys will be returned.

    Returns
    -------
    max_vals : numpy.array
        The values of the peaks found.
    max_idxs : numpy.array
        The x-indices (locations) of the peaks.

    Raises
    ------
    ValueError
        If float not between 0 and 1 passed to threshold.
    """
    peak_vals = []  # a list to hold the y-values of the peaks. Will be converted to a numpy array
    peak_idxs = []  # ditto for x-values (index) of y data.

    y, x = _datacheck_peakdetect(y, x)  # check length and convert to numpy arrays

    if find_min_instead:
        y = -y

    """Exclude data if need be"""
    if exclude_rt_edge or exclude_lt_edge:
        if exclude_rt_edge < 1 and isinstance(exclude_rt_edge, float):
            r_edge = len(y) - int(len(y) * (1 - exclude_rt_edge))
        else:
            r_edge = exclude_rt_edge
        if exclude_lt_edge < 1 and isinstance(exclude_lt_edge, float):
            l_edge = int(len(y) * exclude_lt_edge)
        else:
            l_edge = exclude_lt_edge
        if exclude_rt_edge:
            y = y[:-r_edge]
            x = x[:-r_edge]
        if exclude_lt_edge:
            y = y[l_edge:]
            x = x[l_edge:]

    y_diff = np.diff(y.astype(float))  # y and y_diff must be converted to signed type.

    if isinstance(threshold, float):
        if threshold >= 1:
            raise ValueError("When threshold is passed a float, value must be less than 1")
        else:
            data_range = y.max() - y.min()
            threshold = threshold * data_range + y.min()

    """Find all potential peaks"""
    for idx in range(len(y_diff)-1):
        # For each item of the diff array, check if:
        # 1) The y-value is above the threshold.
        # 2) The value of y_diff is positive (negative for valley search), it means the y-value changed upward.
        # 3) The next y_diff value is zero or negative (or positive for valley search); a positive-then-negative diff value means the value
        # is a peak of some kind. If the diff is zero it could be a flat peak, which still counts.

        # 1)
        if y[idx + 1] < threshold:
            continue

        y1_gradient = y_diff[idx] > 0
        y2_gradient = y_diff[idx + 1] <= 0

        # 2) & 3)
        if y1_gradient and y2_gradient:
            # If the next value isn't zero it's a single-pixel peak. Easy enough.
            if y_diff[idx+1] != 0:
                peak_vals.append(y[idx + 1])
                peak_idxs.append(x[idx + 1])
            # elif idx >= len(y_diff) - 1:
            #     pass
            # Else if the diff value is zero, it could be a flat peak, or it could keep going up; we don't know yet.
            else:
                # Continue on until we find the next nonzero diff value.
                try:
                    shift = 0
                    while y_diff[(idx + 1) + shift] == 0:
                        shift += 1
                        if (idx + 1 + shift) >= (len(y_diff)-1):
                            break
                    # If the next diff is negative (or positive for min), we've found a peak. Also put the peak at the center of the flat
                    # region.
                    is_a_peak = y_diff[(idx + 1) + shift] < 0
                    if is_a_peak:
                        peak_vals.append(y[(idx + 1) + np.round(shift / 2)])
                        peak_idxs.append(x[(idx + 1) + np.round(shift / 2)])
                except IndexError:
                    pass

    # convert to numpy arrays
    peak_vals = np.array(peak_vals)
    peak_idxs = np.array(peak_idxs)
    # try:
    #     peak_idxs += l_edge
    # except:
    #     pass

    """Enforce the min_peak_distance by removing smaller peaks."""
    # For each peak, determine if the next peak is within the min peak width range.
    if isinstance(min_peak_width, float):
        if 0 > min_peak_width >= 1:
            raise ValueError("When min_peak_width is passed a float, value must be between 0 and 1")
        else:
            min_peak_width = int(min_peak_width*len(y))

    index = 0
    while index < len(peak_idxs) - 1:

        # If the second peak is closer than min_peak_distance to the first peak, find the larger peak and remove the other one.
        if peak_idxs[index] > peak_idxs[index+1] - min_peak_width:
            if peak_vals[index] > peak_vals[index+1]:
                idx2del = index + 1
            else:
                idx2del = index
            peak_vals = np.delete(peak_vals, idx2del)
            peak_idxs = np.delete(peak_idxs, idx2del)
        else:
            index += 1

    """If Maximum Number passed, return only up to number given based on a sort of peak values."""
    if max_num_peaks is not None and len(peak_idxs) > max_num_peaks:
        sorted_peak_vals = peak_vals.argsort()  # sorts low to high
        peak_vals = peak_vals[sorted_peak_vals[-max_num_peaks:]]
        peak_idxs = peak_idxs[sorted_peak_vals[-max_num_peaks:]]

    # If we were looking for minimums, convert the values back to the original sign
    if find_min_instead:
        peak_vals = -peak_vals

    return peak_vals, peak_idxs


