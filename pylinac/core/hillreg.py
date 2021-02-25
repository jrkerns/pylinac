"""Perform non-linear regression using a Hill function"""
import numpy as np
import math
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution
import warnings


def hill_func(x, a, b, c, d):  # Hill function
    """ Calculates the Hill function at x
        a: sigmoid low level
        b: sigmoid high level
        c: approximate inflection point
        d: slope of the sigmoid
    """
    return a + (b - a) / (1.0 + (c / x) ** d)


def inv_hill_func(y, a, b, c, d):  # Inverse Hill function
    """ Calculates the inverse Hill function at y
        a: sigmoid low level
        b: sigmoid high level
        c: approximate inflection point
        d: slope of the sigmoid
    """
    if (y > min(a, b)) and (y < max(a, b)) and (d != 0):
        return c*math.pow((y - a)/(b - y), 1/d)
    else:
        return 0


def hill_reg(xData: np.ndarray, yData: np.ndarray):

    # function for genetic algorithm to minimize (sum of squared error)
    def sumOfSquaredError(parameterTuple):
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
    fitted_parameters, pcov = curve_fit(hill_func, xData, yData, genetic_parameters)
    return fitted_parameters
