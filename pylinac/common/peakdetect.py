from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()

from future.builtins import range
import numpy as np

def _datacheck_peakdetect(x_axis, y_axis):
    if x_axis is None:
        x_axis = list(range(len(y_axis)))

    if len(y_axis) != len(x_axis):
        raise ValueError

    #needs to be a numpy array
    y_axis = np.array(y_axis)
    x_axis = np.array(x_axis)
    return x_axis, y_axis

def peak_detect(y,x=None,threshold=0,min_peak_distance=10, delta=0, find_min_instead=False):
    """A function to find the peaks or valleys of a 1-D signal. Uses the difference (np.diff) in signal to find peaks. Current limitations
    include:
    1) Only for use in 1-D data; 2-D may be possible with the gradient function. 2) Will not detect peaks at the very edge of array
    (i.e. 0 or -1 index)

    :param y: 1-D numpy array; y-data of signal.
    :param x: 1-D numpy array; x-data of signal.
    :param threshold: scalar; the value the peak must be above to be considered a peak.
        This removes "peaks" that are in a low-value region.
    :param min_peak_distance: scalar; the number of elements apart a peak must be from neighboring peaks.
    :param delta: scalar; the value a peak candidate must be higher than its neighbors by to be considered a true peak.
    :param find_min_instead: boolean; If True, algorithm will find minimums of y instead of maximums
    :return: two 1-D numpy arrays: max_vals, max_idxs; max_vals contains the y-values of the peaks, max_idxs contains the x-index of the
        peaks.
    """
    peak_vals = []  # a list to hold the y-values of the peaks. Will be converted to a numpy array
    peak_idxs = []  # ditto for x-values (index) of y data.

    x, y = _datacheck_peakdetect(x,y)  # check length and convert to numpy arrays

    y_diff = np.diff(y.astype(float))  # Had problems with uint input. y_diff *must* be converted to signed type. Note: diff is faster than
    # ediff1d

    """Find all potential peaks"""
    for idx in range(len(y_diff)):
        # For each item of the diff array, check if:
        # 1) The index is not at the end,
        # 2) The value of y_diff is positive,
        # 3) The next y_diff value is zero or negative; a positive-then-negative diff value means the value is a peak of some kind. If
        #       the diff is zero it could be a flat peak, which still counts.
        # 4) That the y-value of the potential peak is above the threshold
        # Note: In the midst of noise this will catch a lot of peaks, and peakdetect is faster

        # 1)
        not_at_end = idx != len(y_diff) - 1
        # 2) & 3)
        if not find_min_instead and not_at_end:
            y1_gradient = y_diff[idx] > 0
            y2_gradient = y_diff[idx + 1] <= 0
        elif not_at_end:
            y1_gradient = y_diff[idx] < 0
            y2_gradient = y_diff[idx + 1] >= 0
        # 4)
        # above_threshold = y[idx + 1] >= threshold

        if not_at_end and y1_gradient and y2_gradient and y[idx + 1] >= threshold:
            # If the next value isn't zero it's a single-pixel peak. Easy enough.
            if y2_gradient != 0:
                peak_vals.append(y[idx+1])
                peak_idxs.append(idx+1)
            # Else if the diff value is zero, it could be a flat peak, or it could keep going up; we don't know yet.
            else:
                # Continue on until we find the next nonzero diff value.
                shift = 1
                while y_diff[(idx+1)+shift] == 0:
                    shift += 1
                # If the next diff is negative (or positive for min), we've found a peak. Also put the peak at the center of the flat
                # region.
                if not find_min_instead:
                    is_a_peak = y_diff[(idx + 1) + shift] < 0
                else:
                    is_a_peak = y_diff[(idx + 1) + shift] > 0

                if is_a_peak:
                    peak_vals.append(y[(idx + 1)+np.round(shift/2)])
                    peak_idxs.append((idx + 1)+np.round(shift/2))

    # convert to numpy arrays
    peak_vals = np.array(peak_vals)
    peak_idxs = np.array(peak_idxs)

    """Enforce the min_peak_distance by removing smaller peaks."""
    index = 1
    while index < len(peak_idxs) - 1:
        # For each peak, determine if the next peak is within the look_ahead range.

        # If the previous peak is closer than min_peak_distance, find the larger (or smaller for min) peak and remove other one.
        if peak_idxs[index-1] > peak_idxs[index] - min_peak_distance:
            if ((not find_min_instead and peak_vals[index-1] > peak_vals[index]) or
                (find_min_instead and peak_vals[index-1] < peak_vals[index])):
                idx2del = index
            else:
                idx2del = index-1
            peak_vals = np.delete(peak_vals, idx2del)
            peak_idxs = np.delete(peak_idxs, idx2del)
        # Else if the next peak is closer than min_peak_distance, do the same.
        elif peak_idxs[index+1] < peak_idxs[index] + min_peak_distance:
            if ((not find_min_instead and peak_vals[index + 1] > peak_vals[index]) or
                (find_min_instead and peak_vals[index + 1] < peak_vals[index])):
                idx2del = index
            else:
                idx2del = index + 1
            peak_vals = np.delete(peak_vals, idx2del)
            peak_idxs = np.delete(peak_idxs, idx2del)
        else:
            index += 1

    return peak_vals, peak_idxs

