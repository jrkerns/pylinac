from __future__ import annotations

import math

import numpy as np
from scipy.interpolate import interp1d
from skimage.draw import disk

from . import validators
from .array_utils import (
    is_monotonic,
    is_monotonically_decreasing,
)


def _construct_matrices(
    p: np.ndarray, vertices: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Construct the matrices V and vector P for use in the calculation of the
    projection weights. Ju et al, Equation 8.

    .. note::

        We return the transpose of V so that the dot product is correct. This is a difference from Low not made clear.

    Parameters:
    p : array_like
        The point p as an ndarray.
    vertices : list of array_like
        The vertices of the simplex as a list of ndarrays.

    Returns:
    tuple:
        V : ndarray
            Matrix V constructed from the simplex vertices.
        P : ndarray
            Vector P constructed from point p and the last vertex.
    """
    if vertices.shape[-2] != p.shape[-1]:
        raise ValueError("vertices.shape[-2] must be equal to p.shape[-1]")
    #vertices = np.moveaxis(vertices, source=-1, destination=-2)
    V = vertices[..., :-1, :] - vertices[..., -1, None, :, ]
    P = p - vertices[..., -1, :]
    V = np.moveaxis(V, source=-1, destination=-2)
    return V, P


def _calculate_weights(v: np.ndarray, p: np.ndarray) -> np.ndarray:
    """
    Calculate the weights for the projection of point p onto the simplex's
    support. Ju et al, Equation 7.

    Parameters:
    v : ndarray
        Matrix V from the vertices.
    p : ndarray
        Vector P from point p and the last vertex.

    Returns:
    ndarray
        The weights vector w.
    """
    shape = np.broadcast_shapes(v[..., 0].shape[:-1], p.shape[:-1])
    w_shape = (*shape, v.shape[-1]+1)
    w = np.empty(w_shape)
    VTV = np.einsum("...ij,...ik->...jk", v, v)
    VTV_inv = np.linalg.pinv(VTV)
    w[..., :-1] = np.einsum("...ij,...kj,...k", VTV_inv, v, p)
    w[..., -1] = 1 - np.sum(w[..., :-1], axis=-1)
    return w


def _compute_distance(p: np.ndarray, vertices: np.ndarray) -> np.ndarray:
    """Compute the distance from the point p to the projected point on the simplex's
    support using the weights. This will result in the distance to the support.
    Ju et al, Equation 6.

    .. note::

        Per Low et al, if any of the weights are negative, the projection is outside the simplex.
        In this case, we need to calculate the distance to the boundary.
        In the 1D case, this can be shortcut to be the minimum distance to the vertices.

    Parameters:
    p : array_like
        The point p.
    v : ndarray
        The vertices of the simplex.

    Returns:
    float
        The distance from p to the projected point.
    """
    V, P = _construct_matrices(p, vertices)
    weights = _calculate_weights(v=V, p=P)
    # if any weight is negative it means projection is outside the simplex
    # in the 1D case, this can be shortcut to be the minimum distance to the vertices
    # TBD: handle 2D+ cases; Low Equation 9 is the solution but unsure implementation
    neg_mask = np.any(weights < 0, axis=-1)

    vertices = np.moveaxis(vertices, source=-1, destination=-2)
    def pos_weighs_dist(vertices, p, weights):
        proj = np.sum(weights[..., None, :] * vertices, axis=-1)
        return np.linalg.norm(p - proj, axis=-1)

    def neg_weighs_dist(vertices, p):
        dists = np.sum((vertices - p[..., :, None]) ** 2, axis=-2)
        return np.sqrt(np.min(dists, axis=-1))

    dist = np.where(
        neg_mask, neg_weighs_dist(vertices, p), pos_weighs_dist(vertices, p, weights)
    )
    return dist


def gamma_geometric(
    reference: np.ndarray,
    evaluation: np.ndarray,
    reference_coordinates: np.ndarray | None = None,
    evaluation_coordinates: np.ndarray | None = None,
    dose_to_agreement: float = 1,
    distance_to_agreement: float = 1,
    gamma_cap_value: float = 2,
    dose_threshold: float = 5,
    fill_value: float = np.nan,
) -> np.ndarray:
    """Compute the Ju et al geometric gamma of two 1D profiles/arrays.
    2D support will come in the future and should be doable in-place.

    Parameters
    ----------
    reference
        The reference profile.
    evaluation
        The evaluation profile.
    reference_coordinates
        The x-values of the reference profile. If None, will be assumed to be evenly spaced from 0 to len(reference).
    evaluation_coordinates
        The x-values of the evaluation profile. If None, will be assumed to be evenly spaced from 0 to len(evaluation).
    dose_to_agreement
        The dose to agreement in %. E.g. 1 is 1% of global reference max dose.
    distance_to_agreement
        The distance to agreement in x-values (generally this should be mm).
    gamma_cap_value
        The value to cap the gamma at. E.g. a gamma of 5.3 will get capped to 2. Useful for displaying data with a consistent range.
    dose_threshold
        The dose threshold as a number between 0 and 100 of the % of max dose under which a gamma is not calculated.
    fill_value
        The value to give pixels that were not calculated because they were under the dose threshold. Default
        is NaN, but another option would be 0. If NaN, allows the user to calculate mean/median gamma over just the
        evaluated portion and not be skewed by 0's that should not be considered.

    Returns
    -------
    np.ndarray
        The gamma values. The dimensions will be the same as the evaluation profile.
    """
    if reference.ndim != 1 or evaluation.ndim != 1:
        raise ValueError(
            f"Reference and evaluation arrays must be 1D. Got reference: {reference.ndim} and evaluation: {evaluation.ndim}"
        )
    if distance_to_agreement <= 0:
        raise ValueError("Dose to agreement must be greater than 0")
    if dose_to_agreement <= 0:
        raise ValueError("Distance to agreement must be greater than 0")
    if reference_coordinates is None:
        reference_coordinates = np.arange(len(reference), dtype=float)
    if not is_monotonic(reference_coordinates):
        raise ValueError(
            "Reference x-values must be monotonically increasing or decreasing"
        )
    if len(reference) != len(reference_coordinates):
        raise ValueError(
            f"Reference and reference_x_values must be the same length. Got reference: {len(reference)} and reference_x_values: {len(reference_coordinates)}"
        )
    if evaluation_coordinates is None:
        evaluation_coordinates = np.arange(len(evaluation), dtype=float)
    if not is_monotonic(evaluation_coordinates):
        raise ValueError(
            "Evaluation x-values must be monotonically increasing or decreasing"
        )
    if len(evaluation) != len(evaluation_coordinates):
        raise ValueError(
            f"Evaluation and evaluation_x_values must be the same length. Got evaluation: {len(evaluation)} and evaluation_x_values: {len(evaluation_coordinates)}"
        )
    # normalize the dose threshold by the DTA
    threshold = float(dose_threshold) / float(dose_to_agreement)
    # convert dose to normalized distance of dose to agreement. I.e. D/delta(D) in Figure 1.
    normalized_reference = (
        reference.astype(float) * 100 / (reference.max() * dose_to_agreement)
    )
    normalized_evaluation = (
        evaluation.astype(float) * 100 / (reference.max() * dose_to_agreement)
    )
    # normalize the x-values; i.e. X/delta(d) in Figure 1.
    normalized_reference_x = reference_coordinates / distance_to_agreement
    normalized_evaluation_x = evaluation_coordinates / distance_to_agreement

    gamma = np.full(len(reference), fill_value)
    for idx, (ref_x, ref_point) in enumerate(
        zip(normalized_reference_x, normalized_reference)
    ):
        # skip if below dose threshold
        if ref_point < threshold:
            continue
        # We don't want to calculate gamma for all simplexs of the entire profile,
        # so we slice the vertices just beyond the edges of the DTA from the ref point.
        # For cases where the measurement spacing is 2x or larger than the DTA this
        # leads to the left and right indices being the same. We evaluate an extra
        # index for safety (+1/-1). This adds a small amount of computation but ensures we
        # are always evaluating a range of simplexes. Extra simplex vertices will
        # not change the gamma calculation.
        # This sub-sampling of the vertices is all for computational efficiency.
        left_diffs = np.abs(normalized_evaluation_x - (ref_x - distance_to_agreement))
        right_diffs = np.abs(normalized_evaluation_x - (ref_x + distance_to_agreement))
        if is_monotonically_decreasing(normalized_evaluation_x):
            # it can be the case that the x-values go from high to low vs low to high.
            # If the profiles are decreasing, we need to swap the left and right differences
            # We cannot sort/flip the array itself because then the gamma order
            # would be reversed from the profile order.
            left_diffs, right_diffs = right_diffs, left_diffs
        # we need to ensure we don't go out of bounds if evaluating at the edge, hence the max/min
        left_idx = max(np.argmin(left_diffs) - 1, 0)
        right_idx = min(np.argmin(right_diffs) + 1, len(normalized_evaluation) - 1)
        # the vertices are the (x, y) pairs of the reference profile
        vertices_x = normalized_evaluation_x[left_idx : right_idx + 1]
        vertices_y = normalized_evaluation[left_idx : right_idx + 1]
        vertices = np.stack([vertices_x, vertices_y], axis=1)
        # iterate in pairs of x, y
        vertices_window = np.lib.stride_tricks.sliding_window_view(vertices, 2, axis=0)
        vertices_window = np.moveaxis(vertices_window, source=-1, destination=-2)

        simplex_distances = _compute_distance(
            p=np.array([ref_x, ref_point]), vertices=vertices_window
        )
        gamma[idx] = min(min(simplex_distances), gamma_cap_value)
    return gamma


def gamma_2d(
    reference: np.ndarray,
    evaluation: np.ndarray,
    dose_to_agreement: float = 1,
    distance_to_agreement: int = 1,
    gamma_cap_value: float = 2,
    global_dose: bool = True,
    dose_threshold: float = 5,
    fill_value: float = np.nan,
) -> np.ndarray:
    """Compute a 2D gamma of two 2D numpy arrays. This does NOT do size or spatial resolution checking.
    It performs an element-by-element evaluation. It is the responsibility
    of the caller to ensure the reference and evaluation have comparable spatial resolution.

    The algorithm follows Table I of D. Low's 2004 paper: Evaluation of the gamma dose distribution comparison method: https://aapm.onlinelibrary.wiley.com/doi/epdf/10.1118/1.1598711

    This is similar to the gamma_1d function for profiles, except we must search a 2D grid around the reference point.

    Parameters
    ----------
    reference
        The reference 2D array.
    evaluation
        The evaluation 2D array.
    dose_to_agreement
        The dose to agreement criterion in %. E.g. 1 is 1% of global reference max dose.
    distance_to_agreement
        The distance to agreement criterion in **elements**. E.g. if the value is 4 this means 4 elements from the reference point under calculation.
        Must be >0
    gamma_cap_value
        The value to cap the gamma at. E.g. a gamma of 5.3 will get capped to 2. Useful for displaying data with a consistent range.
    global_dose
        Whether to evaluate the dose to agreement threshold based on the global max or the dose point under evaluation.
    dose_threshold
        The dose threshold as a number between 0 and 100 of the % of max dose under which a gamma is not calculated.
        This is not affected by the global/local dose normalization and the threshold value is evaluated against the global max dose, period.
    fill_value
        The value to give pixels that were not calculated because they were under the dose threshold. Default
        is NaN, but another option would be 0. If NaN, allows the user to calculate mean/median gamma over just the
        evaluated portion and not be skewed by 0's that should not be considered.
    """
    if reference.ndim != 2 or evaluation.ndim != 2:
        raise ValueError(
            f"Reference and evaluation arrays must be 2D. Got reference: {reference.ndim} and evaluation: {evaluation.ndim}"
        )

    # convert dose-to-agreement to % of global-max or % of local-value
    if global_dose:
        dose_ta = dose_to_agreement / 100 * reference.max()
    else:
        dose_ta = dose_to_agreement / 100 * reference

    # Work with normalized values to avoid normalization inside the loop
    eval_normalized = evaluation / dose_ta
    reference_normalized = reference / dose_ta
    threshold_normalized = dose_threshold / 100

    # pad eval array on both edges so our search does not go out of bounds
    eval_normalized = np.pad(eval_normalized, distance_to_agreement, mode="edge")

    # use scikit-image to compute the indices of a disk around the reference point
    # we can then compute gamma over the eval points at these indices
    # we use DTA+1 in disk because it looks like the results are exclusive of edges.
    disk_rr, disk_cc = disk((0, 0), distance_to_agreement + 1)

    # pre-calculate as much as possible not to repeat inside the loop
    # For each row/col these are the indexes which are covered by the disk
    row_r = np.array(range(reference.shape[0]))[np.newaxis].T + disk_rr
    col_r = np.array(range(reference.shape[1]))[np.newaxis].T + disk_cc
    # For the evaluation image the row/col indexes are offset by the padding distance
    row_e = row_r + distance_to_agreement
    col_e = col_r + distance_to_agreement
    # The spatial distance depends only on the disk so it can be precalculated
    dist_row = disk_rr / distance_to_agreement
    dist_col = disk_cc / distance_to_agreement
    dist_r_2 = dist_row**2 + dist_col**2

    gamma_cap_value_2 = gamma_cap_value**2
    gamma = np.full(reference.shape, float(gamma_cap_value))
    # iterate over each reference element, computing distance value and dose value
    for row_idx in range(reference.shape[0]):
        for col_idx in range(reference.shape[1]):
            ref_point = reference_normalized[row_idx, col_idx]

            # skip if below dose threshold
            if math.isnan(ref_point) or ref_point < threshold_normalized:
                gamma[row_idx, col_idx] = fill_value
                continue

            # roi from evaluation
            eval_roi = eval_normalized[row_e[row_idx, :], col_e[col_idx, :]]

            # Normalized dose difference between evaluated and reference dose points
            dist_dose = eval_roi - ref_point

            # capital gamma square (avoid sqrt and memory access if above cap)
            capital_gamma_2 = np.nanmin(dist_r_2 + dist_dose * dist_dose)
            if capital_gamma_2 >= gamma_cap_value_2:
                continue
            gamma[row_idx, col_idx] = np.sqrt(capital_gamma_2)

    return np.asarray(gamma)


def gamma_1d(
    reference: np.ndarray,
    evaluation: np.ndarray,
    reference_coordinates: np.ndarray | None = None,
    evaluation_coordinates: np.ndarray | None = None,
    dose_to_agreement: float = 1,
    distance_to_agreement: int = 1,
    gamma_cap_value: float = 2,
    global_dose: bool = True,
    dose_threshold: float = 5,
    resolution_factor: int = 3,
    fill_value: float = np.nan,
) -> (np.ndarray, np.ndarray, np.ndarray):
    """Perform a 1D gamma of two 1D profiles/arrays.

    The algorithm follows Table I of D. Low's 2004 paper: Evaluation of the gamma dose distribution comparison method: https://aapm.onlinelibrary.wiley.com/doi/epdf/10.1118/1.1598711

    Parameters
    ----------

    reference
        The reference profile.
    evaluation
        The evaluation profile.
    reference_coordinates
        The x-values of the reference profile. If None, will be assumed to be evenly spaced from 0 to len(reference).
    evaluation_coordinates
        The x-values of the evaluation profile. If None, will be assumed to be evenly spaced from 0 to len(evaluation).
    dose_to_agreement
        The dose to agreement in %. E.g. 1 is 1% of global reference max dose.
    distance_to_agreement
        The distance to agreement in x-values. If no x-values are passed, this is in elements.
        If x-values are passed, this is in the units of the x-values. E.g. if the x-values are in mm and the distance to agreement is 1, this is 1 mm.
        Must be >0.
    gamma_cap_value
        The value to cap the gamma at. E.g. a gamma of 5.3 will get capped to 2. Useful for displaying data with a consistent range.
    global_dose
        Whether to evaluate the dose to agreement threshold based on the global max or the dose point under evaluation.
    dose_threshold
        The dose threshold as a number between 0 and 100 of the % of max dose under which a gamma is not calculated.
        This is not affected by the global/local dose normalization and the threshold value is evaluated against the global max dose, period.
    resolution_factor
        The factor by which to resample the reference profile to be compared to the evaluation profile. Depends on the distance to agreement.
        If the distance to agreement is 1 mm and the resolution factor is 3, the reference profile will be resampled to 0.333 mm resolution.
        The rule of thumb is to use at least 3x the distance to agreement. Higher factors will increase computation time.
    fill_value
        The value to give pixels that were not calculated because they were under the dose threshold. Default
        is NaN, but another option would be 0. If NaN, allows the user to calculate mean/median gamma over just the
        evaluated portion and not be skewed by 0's that should not be considered.

    Returns
    -------
    np.ndarray
        The gamma values. The dimensions will be the same as the evaluation profile.
    np.ndarray
        The resampled reference profile. Useful for plotting.
    np.ndarray
        The x-values of the resampled reference profile. Useful for plotting.
    """
    if reference.ndim != 1 or evaluation.ndim != 1:
        raise ValueError(
            f"Reference and evaluation arrays must be 1D. Got reference: {reference.ndim} and evaluation: {evaluation.ndim}"
        )
    if reference_coordinates is None:
        reference_coordinates = np.arange(len(reference), dtype=float)
    if len(reference) != len(reference_coordinates):
        raise ValueError(
            f"Reference and reference_x_values must be the same length. Got reference: {len(reference)} and reference_x_values: {len(reference_coordinates)}"
        )
    if evaluation_coordinates is None:
        evaluation_coordinates = np.arange(len(evaluation), dtype=float)
    if len(evaluation) != len(evaluation_coordinates):
        raise ValueError(
            f"Evaluation and evaluation_x_values must be the same length. Got evaluation: {len(evaluation)} and evaluation_x_values: {len(evaluation_coordinates)}"
        )
    # we add some padding on the check because resampling SingleProfiles
    # can add ~1/2 pixel on each side to retain the same physical size
    # when upsampling.
    if min(evaluation_coordinates) - 1 > min(reference_coordinates) or max(
        evaluation_coordinates
    ) + 1 < max(reference_coordinates):
        raise ValueError(
            "The reference x-values must be within the range of the evaluation x-values"
        )
    if resolution_factor < 1 or not isinstance(resolution_factor, int):
        raise ValueError("Resolution factor must be an integer greater than 0")
    threshold = reference.max() / 100 * dose_threshold
    # convert dose to agreement to % of global max; ignored later if local dose
    dose_ta = dose_to_agreement / 100 * reference.max()

    # we need to interpolate the evaluation profile to the reference DTA search x-values
    # we allow extrapolation due to physical grid size at the edges.
    # I.e. at element 0 our DTA search will go to -1...+1.
    # this only comes into effect if the edge values of the x-values of the reference and evaluation are the same.
    ref_interp = interp1d(
        evaluation_coordinates, evaluation, kind="linear", fill_value="extrapolate"
    )

    eval_interp_array = []
    eval_x_vals = []
    gamma = []
    for ref_x, ref_point in zip(reference_coordinates, reference):
        # skip if below dose threshold
        if ref_point < threshold:
            gamma.append(fill_value)
            continue
        capital_gammas = []
        eval_xs = np.linspace(
            ref_x - distance_to_agreement,
            ref_x + distance_to_agreement,
            num=int(distance_to_agreement * resolution_factor * 2 + 1),
        )
        eval_x_vals.extend(eval_xs)
        eval_vals = ref_interp(eval_xs)
        eval_interp_array.extend(eval_vals)
        for eval_x, eval_point in zip(eval_xs, eval_vals):
            dist = abs(ref_x - eval_x)
            dose = float(ref_point) - float(
                eval_point
            )  # uints can cause overflow errors
            if not global_dose:
                dose_ta = dose_to_agreement / 100 * ref_point
            capital_gamma = math.sqrt(
                dist**2 / distance_to_agreement**2 + dose**2 / dose_ta**2
            )
            capital_gammas.append(capital_gamma)
        gamma.append(min(min(capital_gammas), gamma_cap_value))
    return np.asarray(gamma), np.asarray(eval_interp_array), np.asarray(eval_x_vals)
