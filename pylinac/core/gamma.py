from __future__ import annotations

import math

import numpy as np
from scipy.interpolate import interp1d
from skimage.draw import disk


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
    reference_x_values: np.ndarray | None = None,
    evaluation_x_values: np.ndarray | None = None,
    dose_to_agreement: float = 1,
    distance_to_agreement: int = 1,
    gamma_cap_value: float = 2,
    global_dose: bool = True,
    dose_threshold_percent: float = 5,
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
    reference_x_values
        The x-values of the reference profile. If None, will be assumed to be evenly spaced from 0 to len(reference).
    evaluation_x_values
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
    dose_threshold_percent
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
    if reference_x_values is None:
        reference_x_values = np.arange(len(reference), dtype=float)
    if len(reference) != len(reference_x_values):
        raise ValueError(
            f"Reference and reference_x_values must be the same length. Got reference: {len(reference)} and reference_x_values: {len(reference_x_values)}"
        )
    if evaluation_x_values is None:
        evaluation_x_values = np.arange(len(evaluation), dtype=float)
    if len(evaluation) != len(evaluation_x_values):
        raise ValueError(
            f"Evaluation and evaluation_x_values must be the same length. Got evaluation: {len(evaluation)} and evaluation_x_values: {len(evaluation_x_values)}"
        )
    if min(reference_x_values) > min(evaluation_x_values) or max(
        reference_x_values
    ) < max(evaluation_x_values):
        raise ValueError(
            "The evaluation x-values must be within the range of the reference x-values"
        )
    if resolution_factor < 1 or not isinstance(resolution_factor, int):
        raise ValueError("Resolution factor must be an integer greater than 0")
    threshold = reference.max() / 100 * dose_threshold_percent
    # convert dose to agreement to % of global max; ignored later if local dose
    dose_ta = dose_to_agreement / 100 * reference.max()
    # pad eval array on both edges so our search does not go out of bounds

    # resample the reference to be the desired resolution factor compared to the evaluation if not already
    ref_interp = interp1d(
        reference_x_values, reference, kind="linear", fill_value="extrapolate"
    )

    ref_interp_array = []
    ref_x_vals = []
    gamma = []
    for eval_x, eval_point in zip(evaluation_x_values, evaluation):
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
