

"""Objects that represent common geometric objects or patterns."""
import numpy as np

from pylinac.common.decorators import type_accept


class Line(object):
    """Class modeling a starshot line that is represented by two points."""

    def __init__(self, point1=None, point2=None, m=None, b=None, is_finite=False):
        """Create a line from *either* two distinct points, or an m*x+b definition.

        :param point1, point2: Points along the line
        :type point1, point2: Point
        """
        self.point1 = point1
        self.point2 = point2
        #TODO: add mx+b functionality

    def distance_to_point(self, point):
        """Calculate the distance from the line to a point.

        Equations are from here: http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html

        :type point: Point
        """

        # calculate from m*x+b definition
        if self.point2 is None:
            #TODO: work on this
            pass
        # calculate from 2 points definition
        else:
            lp1 = self.point1
            lp2 = self.point2
            numerator = np.abs((lp2.x - lp1.x)*(lp1.y - point.y) - (lp1.x - point.x)*(lp2.y - lp1.y))
            denominator = np.sqrt((lp2.x - lp1.x)**2 + (lp2.y - lp1.y)**2)
            return numerator/denominator

    def add_to_figure(self, fig, color='w'):
        """Plot the line to the passed figure."""
        fig.axes.plot((self.point1.x, self.point2.x), (self.point1.y, self.point2.y), color=color)



class Point(object):
    """A point with x and y coordinates.

    A namedtuple (Point = namedtuple('Point', ['x', 'y']) is probably more appropriate,
    but they aren't mutable, unlike an class attr, hence a class.
    """
    @type_accept(x=(float, int), y=(float, int))
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

