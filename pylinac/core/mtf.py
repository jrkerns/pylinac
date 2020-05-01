import warnings

import numpy as np
from scipy.interpolate import interp1d

from .decorators import value_accept


class MTF:
    """This class will calculate relative MTF"""

    def __init__(self, lp_spacings, lp_maximums, lp_minimums):
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
        self.mtfs = {}
        self.norm_mtfs = {}
        for (spacing, max, min) in zip(lp_spacings, lp_maximums, lp_minimums):
            self.mtfs[spacing] = (max - min) / (max + min)
        # sort according to spacings
        self.mtfs = {k: v for k, v in sorted(self.mtfs.items(), key=lambda x: x[0])}
        for key, value in self.mtfs.items():
            self.norm_mtfs[key] = value / self.mtfs[lp_spacings[0]]  # normalize to first region

        # check that the MTF drops monotonically by measuring the deltas between MTFs
        # if the delta is increasing it means the MTF rose on a subsequent value
        max_delta = np.max(np.diff(list(self.norm_mtfs.values())))
        if max_delta > 0:
            warnings.warn("The MTF does not drop monotonically; be sure the ROIs are correctly aligned.")

    @value_accept(x=(0, 100))
    def relative_resolution(self, x=50):
        """Return the line pair value at the given rMTF resolution value.

        Parameters
        ----------
        x : float
            The percentage of the rMTF to determine the line pair value. Must be between 0 and 100.
        """
        f = interp1d(list(self.norm_mtfs.values()), list(self.norm_mtfs.keys()), fill_value='extrapolate')
        mtf = f(x / 100)
        if mtf > max(self.spacings):
            warnings.warn(f"MTF resolution wasn't calculated for {x}% that was asked for. The value returned is an extrapolation. Use a higher % MTF to get a non-interpolated value.")
        return float(mtf)

    @classmethod
    def from_high_contrast_diskset(cls, spacings, diskset):
        """Construct the MTF using high contrast disks from the ROI module."""
        maximums = [roi.max for roi in diskset]
        minimums = [roi.min for roi in diskset]
        return cls(spacings, maximums, minimums)

