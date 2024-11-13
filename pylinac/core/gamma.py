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
    p: np.ndarray, vertices: list[np.ndarray]
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
    validators.single_dimension(p)
    # check all vertices are the same length
    for vertex in vertices:
        validators.single_dimension(vertex)
        if len(vertex) != len(p):
            raise ValueError("All vertices must be the same length as the point")
    vertices = np.array(vertices)
    V = vertices[:-1] - vertices[-1]
    P = p - vertices[-1]
    return V.T, P


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
    VTV = np.dot(v.T, v)
    VTV_inv = np.linalg.pinv(VTV)
    w = np.dot(VTV_inv, np.dot(v.T, p))
    wk1 = 1 - np.sum(w)
    return np.append(w, wk1)


def _compute_distance(p: np.ndarray, vertices: list[np.ndarray]) -> float:
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
    if np.any(weights < 0):
        return min(math.dist(p, vertex) for vertex in vertices)
    proj = np.sum(weights[:, np.newaxis] * vertices, axis=0)
    return np.linalg.norm(p - proj)


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

    gamma = np.full(len(evaluation), fill_value)
    for idx, (eval_x, eval_point) in enumerate(
        zip(normalized_evaluation_x, normalized_evaluation)
    ):
        # skip if below dose threshold
        if eval_point < threshold:
            continue
        simplex_distances = []
        # We don't want to calculate gamma for all simplexs of the entire profile,
        # so we slice the vertices just beyond the edges of the DTA from the eval point.
        # For cases where the measurement spacing is 2x or larger than the DTA this
        # leads to the left and right indices being the same. We evaluate an extra
        # index for safety (+1/-1). This adds a small amount of computation but ensures we
        # are always evaluating a range of simplexes. Extra simplex vertices will
        # not change the gamma calculation.
        # This sub-sampling of the vertices is all for computational efficiency.
        left_diffs = np.abs(normalized_reference_x - (eval_x - distance_to_agreement))
        right_diffs = np.abs(normalized_reference_x - (eval_x + distance_to_agreement))
        if is_monotonically_decreasing(normalized_reference_x):
            # it can be the case that the x-values go from high to low vs low to high.
            # If the profiles are decreasing, we need to swap the left and right differences
            # We cannot sort/flip the array itself because then the gamma order
            # would be reversed from the profile order.
            left_diffs, right_diffs = right_diffs, left_diffs
        # we need to ensure we don't go out of bounds if evaluating at the edge, hence the max/min
        left_idx = max(np.argmin(left_diffs) - 1, 0)
        right_idx = min(np.argmin(right_diffs) + 1, len(normalized_reference) - 1)
        # the vertices are the (x, y) pairs of the reference profile
        vertices_x = normalized_reference_x[left_idx : right_idx + 1].tolist()
        vertices_y = normalized_reference[left_idx : right_idx + 1].tolist()
        vertices = list(np.array((x, y)) for x, y in zip(vertices_x, vertices_y))
        # iterate in pairs of x, y
        for v1, v2 in zip(vertices[:-1], vertices[1:]):
            distance = _compute_distance(
                p=np.array([eval_x, eval_point]), vertices=[v1, v2]
            )
            simplex_distances.append(distance)
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
        The dose to agreement in %. E.g. 1 is 1% of global reference max dose.
    distance_to_agreement
        The distance to agreement in **elements**. E.g. if the value is 4 this means 4 elements from the reference point under calculation.
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
    threshold = reference.max() / 100 * dose_threshold
    # convert dose to agreement to % of global max; ignored later if local dose
    dose_ta = dose_to_agreement / 100 * reference.max()
    # pad eval array on both edges so our search does not go out of bounds
    eval_padded = np.pad(evaluation, distance_to_agreement, mode="edge")
    # iterate over each reference element, computing distance value and dose value
    gamma = np.zeros(reference.shape)
    for row_idx, row in enumerate(reference):
        for col_idx, ref_point in enumerate(row):
            # skip if below dose threshold
            if ref_point < threshold:
                gamma[row_idx, col_idx] = fill_value
                continue
            # use scikit-image to compute the indices of a disk around the reference point
            # we can then compute gamma over the eval points at these indices
            # unlike the 1D computation, we have to search at an index offset by the distance to agreement
            # we use DTA+1 in disk because it looks like the results are exclusive of edges.
            # https://scikit-image.org/docs/stable/api/skimage.draw.html#disk
            rs, cs = disk(
                (row_idx + distance_to_agreement, col_idx + distance_to_agreement),
                distance_to_agreement + 1,
            )

            capital_gammas = []
            for r, c in zip(rs, cs):
                eval_point = eval_padded[r, c]
                # for the distance, we compare the ref row/col to the eval padded matrix
                # but remember the padded array is padded by DTA, so to compare distances, we
                # have to cancel the offset we used for dose purposes.
                dist = math.dist(
                    (row_idx, col_idx),
                    (r - distance_to_agreement, c - distance_to_agreement),
                )
                dose = float(eval_point) - float(ref_point)
                if not global_dose:
                    dose_ta = dose_to_agreement / 100 * ref_point
                capital_gamma = math.sqrt(
                    dist**2 / distance_to_agreement**2 + dose**2 / dose_ta**2
                )
                capital_gammas.append(capital_gamma)
            gamma[row_idx, col_idx] = min(np.nanmin(capital_gammas), gamma_cap_value)
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
    if min(reference_coordinates) - 1 > min(evaluation_coordinates) or max(
        reference_coordinates
    ) + 1 < max(evaluation_coordinates):
        raise ValueError(
            "The evaluation x-values must be within the range of the reference x-values"
        )
    if resolution_factor < 1 or not isinstance(resolution_factor, int):
        raise ValueError("Resolution factor must be an integer greater than 0")
    threshold = reference.max() / 100 * dose_threshold
    # convert dose to agreement to % of global max; ignored later if local dose
    dose_ta = dose_to_agreement / 100 * reference.max()

    # we need to interpolate the reference profile to the evaluation DTA search x-values
    # we allow extrapolation due to physical grid size at the edges.
    # I.e. at element 0 our DTA search will go to -1...+1.
    # this only comes into effect if the edge values of the x-values of the reference and evaluation are the same.
    ref_interp = interp1d(
        reference_coordinates, reference, kind="linear", fill_value="extrapolate"
    )

    ref_interp_array = []
    ref_x_vals = []
    gamma = []
    for eval_x, eval_point in zip(evaluation_coordinates, evaluation):
        # skip if below dose threshold
        if eval_point < threshold:
            gamma.append(fill_value)
            continue
        capital_gammas = []
        ref_xs = np.linspace(
            eval_x - distance_to_agreement,
            eval_x + distance_to_agreement,
            num=int(distance_to_agreement * resolution_factor * 2 + 1),
        )
        ref_x_vals.extend(ref_xs)
        ref_vals = ref_interp(ref_xs)
        ref_interp_array.extend(ref_vals)
        for ref_x, ref_point in zip(ref_xs, ref_vals):
            dist = abs(ref_x - eval_x)
            dose = float(eval_point) - float(
                ref_point
            )  # uints can cause overflow errors
            if not global_dose:
                dose_ta = dose_to_agreement / 100 * ref_point
            capital_gamma = math.sqrt(
                dist**2 / distance_to_agreement**2 + dose**2 / dose_ta**2
            )
            capital_gammas.append(capital_gamma)
        gamma.append(min(min(capital_gammas), gamma_cap_value))
    return np.asarray(gamma), np.asarray(ref_interp_array), np.asarray(ref_x_vals)
