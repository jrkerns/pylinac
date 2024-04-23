from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from skimage.feature import canny
from skimage.transform import hough_line, hough_line_peaks

from pylinac.core.array_utils import stretch
from pylinac.core.image import load


class JawOrthogonality:
    """Determine the orthogonality of the jaws of a radiation field. Assumes a square field at a cardinal angle.

    Parameters
    ----------
    path : str
        The path to the image.
    """

    line_angles: dict[str, dict[str, float]]
    result: dict[str, float]

    def __init__(self, path: str | Path):
        self.image = load(path)

    def analyze(self):
        """Analyze the image for jaw orthogonality."""
        edge_image = stretch(self.image.array)
        edge_image = canny(edge_image)

        # Classic straight-line Hough transform
        # Sets a precision of 0.05 degree.
        tested_angles = np.linspace(-np.pi / 2, np.pi / 2, num=360 * 10, endpoint=False)
        h, theta, d = hough_line(edge_image, theta=tested_angles)

        hspace, angles, dists = hough_line_peaks(h, theta, d)
        sorted_angles_idx = np.argsort(np.abs(angles))
        sorted_angles = angles[sorted_angles_idx]
        sorted_dists = dists[sorted_angles_idx]
        # we now have the horizontal lines in the first two indices
        # and the vertical in the last two
        # but we don't know which one is top/bottom or left/right
        # so we need to sort them
        # we can do this by sorting the distances; the lower distance is will be the top/left
        # and the higher distance will be the bottom/right
        line_angles = {}
        if sorted_dists[0] < sorted_dists[1]:
            line_angles["left"] = {"angle": sorted_angles[0], "dist": sorted_dists[0]}
            line_angles["right"] = {"angle": sorted_angles[1], "dist": sorted_dists[1]}
        else:
            line_angles["left"] = {"angle": sorted_angles[1], "dist": sorted_dists[1]}
            line_angles["right"] = {"angle": sorted_angles[0], "dist": sorted_dists[0]}
        if sorted_dists[2] < sorted_dists[3]:
            line_angles["bottom"] = {"angle": sorted_angles[2], "dist": sorted_dists[2]}
            line_angles["top"] = {"angle": sorted_angles[3], "dist": sorted_dists[3]}
        else:
            line_angles["bottom"] = {"angle": sorted_angles[3], "dist": sorted_dists[3]}
            line_angles["top"] = {"angle": sorted_angles[2], "dist": sorted_dists[2]}

        top_left_angle = np.abs(
            np.rad2deg(line_angles["left"]["angle"] - line_angles["top"]["angle"])
        )
        top_right_angle = np.abs(
            np.rad2deg(line_angles["right"]["angle"] - line_angles["top"]["angle"])
        )
        bottom_left_angle = np.abs(
            np.rad2deg(line_angles["left"]["angle"] - line_angles["bottom"]["angle"])
        )
        bottom_right_angle = np.abs(
            np.rad2deg(line_angles["right"]["angle"] - line_angles["bottom"]["angle"])
        )

        result = {
            "top_left": top_left_angle,
            "top_right": top_right_angle,
            "bottom_left": bottom_left_angle,
            "bottom_right": bottom_right_angle,
        }
        self.line_angles = line_angles
        self.result = result

    def results(self) -> dict[str, float]:
        """Return a dict of the results. Keys are 'top_left', 'top_right', 'bottom_left', 'bottom_right'."""
        return self.result

    def plot_analyzed_image(self, show: bool = True):
        """Plot the image with the lines drawn. The lines are the detected jaw edges."""
        colors = ["r", "b", "c", "m"]
        fig, axes = plt.subplots()
        for idx, (key, data) in enumerate(self.line_angles.items()):
            (x0, y0) = data["dist"] * np.array(
                [np.cos(data["angle"]), np.sin(data["angle"])]
            )
            axes.axline(
                (x0, y0),
                slope=np.tan(data["angle"] + np.pi / 2),
                label=key,
                color=colors[idx],
            )
        axes.set_title("Jaw Orthogonality")
        axes.set_axis_off()
        axes.legend()
        self.image.plot(ax=axes, show=show)
