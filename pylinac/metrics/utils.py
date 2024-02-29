from __future__ import annotations

from typing import Callable

import numpy as np
from skimage import measure, segmentation
from skimage.measure._regionprops import RegionProperties
from skimage.segmentation import find_boundaries

from ..core.array_utils import stretch
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


def find_features(
    sample: np.ndarray,
    top_offset: int,
    left_offset: int,
    min_number: int,
    max_number: int,
    dpmm: float,
    detection_conditions: list[Callable],
    radius_mm: float,
    radius_tolerance_mm: float,
    min_separation_mm: float,
) -> (list[Point], list[np.ndarray], list[RegionProperties]):
    """Find up to max_number features in the image. This will look for features at every percentile range.
    Multiple features may be found at different threshold levels; they are combined and deduplicated.

    Parameters
    ----------
    sample : np.ndarray
        The image to search for BBs in. This can be an entire image or a cropped region.
    top_offset : int
        The top offset INDEX of the cropped region. This is used to calculate the absolute position of the BBs. This is
        used when the sample is a cropped region.
    left_offset : int
        The left offset INDEX of the cropped region. This is used to calculate the absolute position of the BBs.
        This is used when the sample is a cropped region.
    min_number : int
        The minimum number of BBs to search for.
    max_number : int
        The maximum number of BBs to search for. This is used to stop the search early if the number of BBs found.
    dpmm : float
        The dots per mm of the image. This is used to convert the radius and separation to pixels.
    detection_conditions : list[Callable]
        A list of callables that take a region and return True if the region is a BB.
    radius_mm : float
        The radius of the BB in mm. This is used to calculate the minimum separation between BBs.
    radius_tolerance_mm : float
        The tolerance of the BB radius in mm. This is used to calculate the minimum and maximum BB size.
    min_separation_mm : float
        The minimum separation between BBs in mm.

    Returns
    -------
    list[Point]
        The list of BBs found in the image.
    list[np.ndarray]
        The list of boundaries of the BBs found in the image.
    list[RegionProperties]
        The list of scikit-image regions found in the image.
    """
    total_features = []
    feature_boundaries = []
    feature_regions = {}
    # uses the same algo as original WL; this is better than a percentile method as the percentile method
    # can often be thrown off at the very ends of the distribution. It's more linear and faster to use the simple
    # spread of min/max.
    sample = stretch(sample, min=0, max=1)
    imin, imax = sample.min(), sample.max()
    spread = imax - imin
    step_size = spread / 50  # move in 1/50 increments; maximum of 50 passes per image
    cutoff = (
        imin + step_size
    )  # start at the min + 1 step; we know the min cutoff will be a blank, full image
    while cutoff <= imax and len(total_features) < max_number:
        try:
            binary_array = sample > cutoff
            labeled_arr = measure.label(binary_array, connectivity=1)
            cleared_labeled_arr = segmentation.clear_border(labeled_arr)
            regions = measure.regionprops(cleared_labeled_arr, intensity_image=sample)
            feature_regions = {i: r for i, r in enumerate(regions)}
            for condition in detection_conditions:
                to_pop = []
                for key, region in sorted(
                    feature_regions.items(),
                    key=lambda item: item[1].filled_area,
                    reverse=True,
                ):
                    if not condition(
                        region,
                        dpmm=dpmm,
                        bb_size=radius_mm,
                        tolerance=radius_tolerance_mm,
                        shape=binary_array.shape,
                    ):
                        to_pop.append(key)
                feature_regions = {
                    key: region
                    for key, region in feature_regions.items()
                    if key not in to_pop
                }
            if len(feature_regions) == 0:
                raise ValueError
            else:
                new_points = [
                    Point(region.weighted_centroid[1], region.weighted_centroid[0])
                    for region in feature_regions.values()
                ]
                new_boundaries = [
                    get_boundary(
                        detected_region,
                        top_offset=top_offset,
                        left_offset=left_offset,
                    )
                    for detected_region in feature_regions.values()
                ]
                total_features, feature_boundaries = deduplicate_points_and_boundaries(
                    original_points=total_features,
                    new_points=new_points,
                    min_separation_px=min_separation_mm * dpmm,
                    original_boundaries=feature_boundaries,
                    new_boundaries=new_boundaries,
                )
        except (IndexError, ValueError):
            pass
        finally:
            cutoff += step_size
    if len(total_features) < min_number:
        # didn't find the number we needed
        raise ValueError(
            f"Couldn't find the minimum number of disks in the image. Found {len(total_features)}; required: {min_number}"
        )
    # adjust points by the offset
    for feature in total_features:
        feature.x += left_offset
        feature.y += top_offset
    return total_features, feature_boundaries, list(feature_regions.values())
