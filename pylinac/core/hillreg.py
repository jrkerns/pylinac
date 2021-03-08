"""Perform non-linear regression using a Hill function."""

import numpy as np
import math
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution
import warnings


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


def infl_point_hill_func(fit_params) -> float:
    """calculates the inflection point of the Hill function.

        [0] : sigmoid low level
        [1] : sigmoid high level
        [2] : approximate inflection point
        [3] : slope of the sigmoid
    """
    return fit_params[2]*math.pow((fit_params[3] - 1)/(fit_params[3] + 1), 1/fit_params[3])


def hill_reg(xData: np.ndarray, yData: np.ndarray):
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

    def sumOfSquaredError(parameterTuple):
        """function for genetic algorithm to minimize (sum of squared error)"""
        warnings.filterwarnings("ignore")  # do not print warnings by genetic algorithm
        val = hill_func(xData, *parameterTuple)
        return np.sum((yData - val) ** 2.0)

    def generate_initial_parameters(xData, yData):
        # min and max used for bounds
        maxX = max(xData)
        minX = min(xData)
        maxY = max(yData)

        parameterBounds = []
        parameterBounds.append([0, maxY])  # search bounds for a
        parameterBounds.append([0, maxY])  # search bounds for b
        parameterBounds.append([minX, maxX])  # search bounds for c
        parameterBounds.append([-100, 100])  # search bounds for d

        # "seed" the numpy random number generator for repeatable results
        result = differential_evolution(sumOfSquaredError, parameterBounds, seed=3)
        return result.x

    # generate initial parameter values
    genetic_parameters = generate_initial_parameters(xData, yData)

    # curve fit the data
    fitted_parameters, _ = curve_fit(hill_func, xData, yData, genetic_parameters)
    return fitted_parameters
