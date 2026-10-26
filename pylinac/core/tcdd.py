from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit


def calculate_ctsm(rois: list[np.ndarray]) -> tuple[float, float, int]:
    """Calculate threshold contrast using the Statistical Method (Chao et al., 2000).

    Parameters
    ----------
    rois : list of ndarray
        ROIs of identical size sampled from a uniform phantom image.

    Returns
    -------
    c_t : float
        Threshold contrast at 99.9% confidence level (3.29 × σ).
    sigma_chi : float
        Standard deviation of the ROI means.
    n_rois : int
        Number of ROIs used.

    Notes
    -----
    C_T(d) = 3.29 × σ_χ(d), where σ_χ is the standard deviation of the mean
    pixel values across the n ROIs of size d. The factor 3.29 corresponds to
    the 99.9% confidence level of a normal distribution.

    References
    ----------
    Chao, E.H., et al. (2000). Radiology, 217(1), 162.
    Paruccini, N., et al. (2021). Physica Medica, 91, 105-112.
    """
    roi_means = np.array([np.mean(roi) for roi in rois])
    sigma_chi = float(np.std(roi_means, ddof=1))
    c_t = 3.29 * sigma_chi
    return c_t, sigma_chi, len(rois)


def contrast_threshold_model(d: np.ndarray, a: float, b: float, c: float) -> np.ndarray:
    """Contrast threshold model: C_T(d) = c/d² + b/d + a.

    Parameters
    ----------
    d : ndarray
        Detail diameter in mm.
    a, b, c : float
        Model parameters.

    Returns
    -------
    ndarray
        Threshold contrast values.
    """
    return c / (d**2) + b / d + a


def fit_contrast_threshold_curve(
    diameters: np.ndarray,
    c_t_values: np.ndarray,
) -> tuple[np.ndarray | None, dict | None]:
    """Fit the contrast threshold curve C_T(d) = c/d² + b/d + a.

    Parameters
    ----------
    diameters : ndarray
        Detail diameters in mm.
    c_t_values : ndarray
        Threshold contrast values (normalised, in %).

    Returns
    -------
    fitted_params : ndarray or None
        [a, b, c] parameters, or None if fitting failed.
    fit_quality : dict or None
        Keys: ``r_squared``, ``rmse``, ``d_smooth``, ``fitted_curve``.
        None if fitting failed.
    """
    try:
        popt, _ = curve_fit(
            contrast_threshold_model,
            diameters,
            c_t_values,
            p0=[0.01, 0.1, 1.0],
            maxfev=10000,
        )
        a, b, c = popt
        d_smooth = np.linspace(diameters.min(), diameters.max(), 100)
        fitted_curve = contrast_threshold_model(d_smooth, a, b, c)
        y_pred = contrast_threshold_model(diameters, a, b, c)
        ss_res = np.sum((c_t_values - y_pred) ** 2)
        ss_tot = np.sum((c_t_values - np.mean(c_t_values)) ** 2)
        r_squared = 1.0 - ss_res / ss_tot if ss_tot != 0 else float("nan")
        rmse = float(np.sqrt(np.mean((c_t_values - y_pred) ** 2)))
        fit_quality = {
            "r_squared": r_squared,
            "rmse": rmse,
            "d_smooth": d_smooth,
            "fitted_curve": fitted_curve,
        }
        return popt, fit_quality
    except Exception:
        return None, None
