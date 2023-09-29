from __future__ import annotations

from math import ceil, floor
from typing import Sequence

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import stats

from .core import image
from .core.array_utils import invert
from .core.profile import find_peaks
from .picketfence import MLC


class DLG:
    """Analyze a machine's dosimetric leaf gap by looking at profiles with various amounts of overlap. This is NOT the
    same procedure as the sweeping gaps as provided by Varian, although the determined value should be similar.
    """

    def __init__(self, path):
        self.image = image.LinacDicomImage(path)
        self.measured_dlg: float = -np.inf
        self.measured_dlg_per_leaf: list = []
        self.planned_dlg_per_leaf: list = []
        self._lin_fit = None

    def analyze(
        self,
        gaps: Sequence,
        mlc: MLC,
        y_field_size: float = 100,
        profile_width: int = 10,
    ):
        """Analyze an EPID image with varying MLC overlaps to determine the DLG.

        Parameters
        ----------
        gaps
            The gaps (i.e. overlap) of the leaves in mm.
            These should typically be in descending order and also be negative. E.g. (-1, ..., -2.2).

        mlc
            The MLC type/arrangement. This lets us know where the leaf centers are to take a profile along.

        y_field_size
            The field size along the y-dimension (perpendicular to the leaf travel). This will determined which leaves
            are associated with which gap.

        profile_width
            The width of the profile to take along the axes parallel to leaf motion. This should be a good bit wider
            than the gap values. The default is reasonable and it is unlikely it needs tweaking.
        """
        measured_dlg_per_leaf = []
        planned_dlg_per_leaf = []
        mlc = mlc.value["arrangement"]
        g = list(gaps)
        g.sort()
        profile_width_px = round(self.image.dpmm * profile_width)
        mid_width = self.image.shape[1] / 2
        mid_height = self.image.shape[0] / 2
        for idx, center in enumerate(mlc.centers):
            if -y_field_size / 2 < center < y_field_size / 2:
                # get the pixel window area
                center_px = center * self.image.dpmm
                width_px = mlc.widths[idx] / 4 * self.image.dpmm
                top = ceil(mid_height + center_px + width_px)
                bottom = floor(mid_height + center_px - width_px)
                # sample the window and take the average perpendicular to MLC motion
                window = self.image[
                    bottom:top,
                    int(mid_width - profile_width_px) : int(
                        mid_width + profile_width_px
                    ),
                ]
                width = self._determine_measured_gap(window.mean(axis=0))
                planned_dlg_per_leaf.append(
                    self._get_dlg_offset(y_field_size, center, g)
                )
                measured_dlg_per_leaf.append(width)
        # fit the data to a line and determine the DLG from the 0 intercept
        lin_fit = stats.linregress(planned_dlg_per_leaf, measured_dlg_per_leaf)
        dlg = lin_fit.intercept / lin_fit.slope
        self._lin_fit = lin_fit
        self.measured_dlg = dlg
        self.planned_dlg_per_leaf = planned_dlg_per_leaf
        self.measured_dlg_per_leaf = measured_dlg_per_leaf

    def plot_dlg(self, show: bool = True) -> None:
        """Plot the measured DLG values across the planned gaps"""
        if not self.measured_dlg_per_leaf:
            raise ValueError("Analyze the image before plotting with .analyze()")
        plt.plot(self.planned_dlg_per_leaf, self.measured_dlg_per_leaf, "gx")
        plt.plot(
            self.planned_dlg_per_leaf,
            self._lin_fit.intercept
            + self._lin_fit.slope * np.array(self.planned_dlg_per_leaf),
            "r",
            label="fitted line",
        )
        plt.title(f"Measured DLG: {self.measured_dlg:2.3f}mm")
        plt.grid()
        if show:
            plt.show()

    @staticmethod
    def _get_dlg_offset(field_size: float, leaf_center: float, dlgs: Sequence) -> float:
        """Return the planned leaf overlap for a given leaf"""
        roi_size = field_size / len(dlgs)
        y_bounds = [field_size / 2 - idx * roi_size for idx in range(len(dlgs) + 1)]
        for idx, gap in enumerate(dlgs):
            upper_bound = y_bounds[idx]
            lower_bound = y_bounds[idx + 1]
            if lower_bound < leaf_center < upper_bound:
                return gap

    @staticmethod
    def _determine_measured_gap(profile: np.ndarray) -> float:
        """Return the measured gap based on profile height"""
        mid_value = profile[int(len(profile) / 2)]
        if mid_value < profile.mean():
            profile = invert(profile)
        _, props = find_peaks(profile, max_number=1)
        if mid_value < profile.mean():
            return -props["prominences"][0]
        else:
            return props["prominences"][0]
