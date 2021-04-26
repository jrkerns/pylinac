"""Perform non-linear regression using a Hill function."""

import numpy as np
import math
from scipy.optimize import curve_fit


class Hill:
    """A Hill function is for modeling a field falloff (i.e. penumbra). It's not a perfect model, but it fits to a
    function, which is not limited by resolution issues as may be experienced on low-res devices like ion chamber arrays."""

    @classmethod
    def fit(cls, x_data: np.ndarray, y_data: np.ndarray):
        """Fit x & y data to a Hill function."""
        fitted_parameters, _ = curve_fit(hill_func, x_data, y_data, p0=(min(y_data), max(y_data), np.median(x_data), 0))
        instance = cls()
        instance.params = fitted_parameters
        return instance

    def inflection_idx(self) -> dict:
        """Determine the x-value inflection point of the fitted Hill function."""
        idx = self.params[2] * math.pow((self.params[3] - 1) / (self.params[3] + 1), 1 / self.params[3])
        return {'index (exact)': idx, 'index (rounded)': int(round(idx))}

    @classmethod
    def from_params(cls, params):
        """Create a Hill function from pre-determined parameters. Useful to recreate a Hill function"""
        instance = cls()
        instance.params = params
        return instance

    def gradient_at(self, x: float) -> float:
        """Return the gradient of the Hill function at a given x-value"""
        cxd = math.pow(self.params[2]/x, self.params[3])
        return (self.params[1] - self.params[0])*self.params[3]*cxd/(math.pow(cxd + 1, 2)*x)

    def x(self, y: float) -> float:
        """Return the x-value given a y-value"""
        return self.params[2] * math.pow((y - self.params[0]) / (self.params[1] - y), 1 / self.params[3])

    def y(self, x: float) -> float:
        """Return the y-value given an x-value."""
        return self.params[0] + (self.params[1] - self.params[0]) / (1 + (self.params[2]/x) ** self.params[3])



def hill_func(x, a, b, c, d):  # Hill function
    """Calculates the Hill function at x.

        a : sigmoid low level
        b : sigmoid high level
        c : approximate inflection point
        d : slope of the sigmoid
    """
    return a + (b - a) / (1.0 + (c / x) ** d)


def inv_hill_func(y, fit_params):  # Inverse Hill function
    """Calculates the inverse Hill function at y.

        [0] : sigmoid low level
        [1] : sigmoid high level
        [2] : approximate inflection point
        [3] : slope of the sigmoid
    """
    if (y > min(fit_params[0], fit_params[1])) and (y < max(fit_params[0], fit_params[1])) and (fit_params[3] != 0):
        return fit_params[2]*math.pow((y - fit_params[0])/(fit_params[1] - y), 1/fit_params[3])
    else:
        return 0


def deriv_hill_func(x, fit_params) -> float:
    """calculates the tangent of the Hill function at X.

        [0] : sigmoid low level
        [1] : sigmoid high level
        [2] : approximate inflection point
        [3] : slope of the sigmoid
    """
    if x > 0:
        cxd = math.pow(fit_params[2]/x, fit_params[3])
        return (fit_params[1] - fit_params[0])*fit_params[3]*cxd/(math.pow(cxd + 1, 2)*x)
    else:
        return 0


def inflection(fit_params) -> float:
    """calculates the inflection point of the Hill function.

        [0] : sigmoid low level
        [1] : sigmoid high level
        [2] : approximate inflection point
        [3] : slope of the sigmoid
    """
    return fit_params[2]*math.pow((fit_params[3] - 1)/(fit_params[3] + 1), 1/fit_params[3])


def fit_to_hill(xData: np.ndarray, yData: np.ndarray):
    """Performs non-linear least squares regression on a Hill (sigmoid) function.

       Parameters
       ----------
       xData: X values of the function
       yData: Y values of the function

       Returns
       -------
       Fitted Parameters
            [0] : sigmoid low level
            [1] : sigmoid high level
            [2] : approximate inflection point
            [3] : slope of the sigmoid
       """
    fitted_parameters, _ = curve_fit(hill_func, xData, yData, p0=(min(yData), max(yData), np.median(xData), 0))
    return fitted_parameters
