from __future__ import annotations

import numpy as np
from skimage.measure._regionprops import RegionProperties
from skimage.segmentation import find_boundaries

from ..core.geometry import Point


def deduplicate_points_and_boundaries(
    original_points: list[Point],
    new_points: list[Point],
    min_separation_px: float,
    original_boundaries: list[np.ndarray],
    new_boundaries: list[np.ndarray],
) -> (list[Point], list[np.ndarray]):
    """Deduplicate points that are too close together. The original points should be the
    starting point. The new point's x, y, and z values are compared to the existing points.
    If the new point is too close to the original point, it's dropped. If it's sufficiently
    far away, it is added. Will return a new combined list of points.

    We assume the original points are already deduplicated. When used in a loop starting from an empty list
    this is true."""
    combined_points = original_points
    combined_boundaries = original_boundaries
    for new_point, new_boundary in zip(new_points, new_boundaries):
        for original_point in original_points:
            if new_point.distance_to(original_point) < min_separation_px:
                break
        else:
            combined_points.append(new_point)
            combined_boundaries.append(new_boundary)
    return combined_points, combined_boundaries


def get_boundary(
    region: RegionProperties, top_offset: int, left_offset: int
) -> np.ndarray:
    """Find the boundary of the region as the absolute position in the image.
    This will calculate the outline of the region. Mostly used for plotting."""
    # we pad the boundaries to offset the ROI to the right and down to the absolute position
    # on the image.
    return np.pad(
        find_boundaries(
            # padding is needed as boundary edges aren't detected otherwise
            np.pad(
                region.image,
                pad_width=1,
                mode="constant",
                constant_values=0,
            ),
            connectivity=region.image.ndim,
            mode="inner",
            background=0,
        ),
        ((region.bbox[0] + top_offset - 1, 0), (region.bbox[1] + left_offset - 1, 0)),
        mode="constant",
        constant_values=0,
    )
