
"""Module for classes that represent common geometric objects or patterns."""

import numpy as np
from matplotlib.patches import Circle as mpl_Circle
from matplotlib.patches import Rectangle as mpl_Rectangle

from pylinac.core.decorators import type_accept, lazyproperty


class Point(object):
    """A point with x, y, and z coordinates.

    A namedtuple (Point = namedtuple('Point', ['x', 'y']) is probably more appropriate,
    but they aren't mutable, unlike an class attr, hence a class.
    """

    def __init__(self, x=0, y=0, z=0, as_int=False):
        if as_int:
            x = int(x)
            y = int(y)
            z = int(z)
        self.x = x
        self.y = y
        self.z = z


class Circle(object):
    """A circle with center Point and radius."""
    def __init__(self, center_point=None, radius=None):

        if center_point is not None or not isinstance(center_point, Point):
            raise TypeError("Circle center must be of type Point")

        self.center = center_point
        self.radius = radius

    def add_to_axes(self, axes, edgecolor='black', fill=False):
        """Plot the Circle on the axes."""
        axes.add_patch(mpl_Circle((self.center.x, self.center.y), edgecolor=edgecolor, radius=self.radius, fill=fill))


class Line(object):
    """Model a line that is represented by two points.

    Calculations of slope, etc are from here:
    http://en.wikipedia.org/wiki/Linear_equation
    and here:
    http://www.mathsisfun.com/algebra/line-equation-2points.html
    """
    def __init__(self, point1=None, point2=None, m=None, b=None, is_finite=False):
        """Create a line from *either* two distinct points, or an m*x+b definition.

        :param point1, point2: Points along the line
        :type point1, point2: Point
        """
        #TODO: incorporate is_finite
        # if created by passing two points...
        if isinstance(point1, Point) and isinstance(point2, Point):
            self.point1 = point1
            self.point2 = point2
        # otherwise by passing m and b...
        elif m is not None and b is not None:
            self.m = m
            self.b = b
        else:
            raise ValueError("Proper parameters not passed for proper Line instantiation")

    @lazyproperty
    def m(self):
        """Return the slope of the line.

        m = (y1 - y2)/(x1 - x2)

        From: http://www.purplemath.com/modules/slope.htm
        """
        return (self.point1.y - self.point2.y) / (self.point1.x - self.point2.x)

    @lazyproperty
    def b(self):
        """Return the intercept of the line.

        b = y - m*x
        """
        return self.point1.y - (self.m * self.point1.x)

    def y(self, x):
        """Return y-value along line given x."""
        return self.m * x + self.b

    def x(self, y):
        """Return x-value along line given y."""
        return (y - self.b)/self.m

    @type_accept(point=(Point, tuple))
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

    def add_to_axes(self, axes, color='w'):
        """Plot the line to the passed figure."""
        axes.plot((self.point1.x, self.point2.x), (self.point1.y, self.point2.y), color=color)


class Box(object):
    """A box object with width, height, center Point, top-left corner Point, and bottom-left corner Point."""
    @type_accept(center=(Point, None), tl_corner=(Point, None), bl_corner=(Point, None))
    def __init__(self, width, height, center=None, tl_corner=None, bl_corner=None, as_int=False):
        self._as_int = as_int
        self.width = width
        self.height = height
        if not any((center, tl_corner, bl_corner)):
            raise ValueError("Must at least specify one anchor point for the box.")
        elif center is not None:
            self.center = center
            self.tl_corner = Point(center.x - width/2, center.y + height/2)
            self.bl_corner = Point(center.x - width/2, center.y - height/2)
        elif tl_corner is not None:
            self.tl_corner = tl_corner
            self.center = Point(tl_corner.x + width/2, tl_corner.y - height/2)
            self.bl_corner = Point(tl_corner.x, tl_corner.y - height)
        elif bl_corner is not None:
            self.bl_corner = bl_corner
            self.center = Point(bl_corner.x + width / 2, bl_corner.y + height / 2)
            self.tl_corner = Point(bl_corner.x, bl_corner.y + height)

        if as_int:
            self.width = int(self.width)
            self.height = int(self.height)
            self.center = Point(int(self.center.x), int(self.center.y))
            self.bl_corner = Point(int(self.bl_corner.x), int(self.bl_corner.y))
            self.tl_corner = Point(int(self.tl_corner.x), int(self.tl_corner.y))

    def add_to_axes(self, axes, edgecolor='black', angle=0.0, fill=False):
        """Plot the Box to the axes."""
        axes.add_patch(mpl_Rectangle((self.center.x, self.center.y),
                                     width=self.width,
                                     height=self.height,
                                     angle=angle,
                                     edgecolor=edgecolor,
                                     fill=fill))

