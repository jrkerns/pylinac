
import numpy as np

"""Common Functions used in Pylinac"""


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
    """Find the peaks or valleys of a 1-D signal. Uses the difference (np.diff) in signal to find peaks. Current limitations
    include:
    1) Only for use in 1-D data; 2-D may be possible with the gradient function. 2) Will not detect peaks at the very edge of array
    (i.e. 0 or -1 index)

    :param y: 1-D y-data of signal.
    :type y: numpy array
    :param x: 1-D x-data of signal. If left as None, will create a uniform range the length of the y-data.
    :type x: numpy array
    :param threshold: The value the peak must be above to be considered a peak.
        This removes "peaks" that are in a low-value region. If passed an int, the actual value is the threshold;
        if a float <1.0 is passed, it will threshold that percent. E.g. when passed 15, any peak less than 15 is
        removed. When passed 0.4, any peak less than 40% of the maximum value will be removed.
    :type threshold: int, float
    :param min_peak_width: The number of elements apart a peak must be from neighboring peaks.
    :type min_peak_width: int
    :param max_num_peaks: Specify up to how many peaks will be returned. E.g. if 3 is passed in and 5 peaks are found, only the 3 largest
        peaks will be returned.
    :param find_min_instead: If True, algorithm will find minimums of y instead of maximums.
    :type find_min_instead: bool
    :returns: two 1-D numpy arrays: max_vals, max_idxs; max_vals contains the y-values of the peaks, max_idxs contains the x-index of the
        peaks.
    """
    peak_vals = []  # a list to hold the y-values of the peaks. Will be converted to a numpy array
    peak_idxs = []  # ditto for x-values (index) of y data.

    y, x = _datacheck_peakdetect(y, x)  # check length and convert to numpy arrays

    if find_min_instead:
        y = -y

    """Exclude data if need be"""
    if exclude_rt_edge or exclude_lt_edge:
        if exclude_rt_edge < 1 and isinstance(exclude_rt_edge, float):
            r_edge = int(len(y) * (1 - exclude_rt_edge))
        else:
            r_edge = exclude_rt_edge
        if exclude_lt_edge < 1 and isinstance(exclude_lt_edge, float):
            l_edge = int(len(y) * exclude_lt_edge)
        else:
            l_edge = exclude_lt_edge
        if exclude_rt_edge:
            y = y[:r_edge]
            x = x[:r_edge]
        if exclude_lt_edge:
            y = y[l_edge:]
            x = x[l_edge:]

    y_diff = np.diff(y.astype(float))  # Had problems with uint input. y_diff *must* be converted to signed type.

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
            # Else if the diff value is zero, it could be a flat peak, or it could keep going up; we don't know yet.
            else:
                # Continue on until we find the next nonzero diff value.
                shift = 1
                while y_diff[(idx + 1) + shift] == 0:
                    shift += 1
                # If the next diff is negative (or positive for min), we've found a peak. Also put the peak at the center of the flat
                # region.
                is_a_peak = y_diff[(idx + 1) + shift] < 0
                if is_a_peak:
                    peak_vals.append(y[(idx + 1) + np.round(shift / 2)])
                    peak_idxs.append(x[(idx + 1) + np.round(shift / 2)])

    # convert to numpy arrays
    peak_vals = np.array(peak_vals)
    peak_idxs = np.array(peak_idxs)

    """Enforce the min_peak_distance by removing smaller peaks."""
    index = 0
    while index < len(peak_idxs) - 1:
        # For each peak, determine if the next peak is within the look_ahead range.

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


