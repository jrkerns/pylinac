
"""Module for classes that represent common geometric objects or patterns."""
import math

import numpy as np
from matplotlib.patches import Circle as mpl_Circle
from matplotlib.patches import Rectangle as mpl_Rectangle

from pylinac.core.decorators import type_accept, lazyproperty
from pylinac.core.utilities import isnumeric, is_iterable


class Point:
    """A point with x, y, and z coordinates.

    .. note:: A namedtuple (Point = namedtuple('Point', ['x', 'y']) is probably more appropriate,
              but they aren't mutable, unlike an class attr, hence a class.
    """

    def __init__(self, x=0, y=0, z=0, idx=None, value=None, as_int=False):
        """
        :param x: x-coordinate or iterable type containing all coordinates. If iterable, values are assumed to be in order: (x,y,z).
        :type x: int, float, iterable (list, tuple, etc)
        :param y: y-coordinate
        :type y: numeric, optional
        :param z: z-coordinate
        :type z: numeric, optional
        :param idx: Index of point. Useful for sequential coordinates; e.g. a point on a circle profile is sometimes easier to describe
            in terms of its index rather than x,y coords.
        :type idx: int, optional
        :param value: value at point location (e.g. pixel value of an image)
        :type value: numeric, optional
        :param as_int: flag specifying to convert coordinates to integers
        :type as_int: boolean
        """
        # Point object passed in
        if isinstance(x, Point):
            point = x
            x = point.x
            y = point.y
            z = point.z
            idx = point.idx
            value = point.value

        # if passed an iterable, separate out
        elif is_iterable(x):
            input_coords = x
            try:
                x = input_coords[0]
                y = input_coords[1]
                z = input_coords[2]
                idx = input_coords[3]
                value = input_coords[4]
            except IndexError:
                pass

        if as_int:
            x = int(x)
            y = int(y)
            z = int(z)

        self.x = x
        self.y = y
        self.z = z
        self.idx = idx
        self.value = value

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, val):
        if not isnumeric(val):
            if val is not None:
                raise TypeError("Point value was not a valid type. Must be numeric.")
        self._value = val

    def dist_to(self, point):
        """Calculate the distance to the given point."""
        return math.sqrt((self.x - point.x)**2 + (self.y - point.y)**2)


class Circle:
    """A circle with center Point and radius."""
    def __init__(self, center_point=None, radius=None):

        if center_point is not None:
            if is_iterable(center_point):
                # if iterable, convert to Point because Point will convert
                center_point = Point(center_point)
            elif not isinstance(center_point, Point):
                raise TypeError("Circle center must be of type Point or iterable")
        elif center_point is None:
            center_point = Point()

        self.center = center_point
        self.radius = radius

    def add_to_axes(self, axes, edgecolor='black', fill=False):
        """Plot the Circle on the axes."""
        axes.add_patch(mpl_Circle((self.center.x, self.center.y), edgecolor=edgecolor, radius=self.radius, fill=fill))


class Line:
    """Model a line that is represented by two points or an m*x+b representation.

    Calculations of slope, etc are from here:
    http://en.wikipedia.org/wiki/Linear_equation
    and here:
    http://www.mathsisfun.com/algebra/line-equation-2points.html
    """
    def __init__(self, point1=None, point2=None, m=None, b=None):
        """Create a line from *either* two distinct points, or an m*x+b definition.

        :param point1, point2: Points along the line
        :type point1, point2: Point
        """
        # if created by passing two points...
        if point1 is not None and point2 is not None:
            self.point1 = Point(point1)
            self.point2 = Point(point2)
        # otherwise by passing m and b...
        elif m is not None and b is not None:
            self.m = m
            self.b = b

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

    @property
    def length(self):
        """Return length of the line."""
        return self.point1.dist_to(self.point2)

    @type_accept(point=(Point, tuple))
    def distance_to(self, point):
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


class Rectangle(object):
    """A rectangle with width, height, center Point, top-left corner Point, and bottom-left corner Point."""
    @type_accept(center=(Point, None), tl_corner=(Point, None), bl_corner=(Point, None))
    def __init__(self, width, height, center=None, tl_corner=None, bl_corner=None, as_int=False):
        self._as_int = as_int
        self.width = width
        self.height = height
        if not any((center, tl_corner, bl_corner)):
            raise ValueError("Must specify at least one anchor point for the box.")
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
        """Plot the Rectangle to the axes."""
        axes.add_patch(mpl_Rectangle((self.bl_corner.x, self.bl_corner.y),
                                     width=self.width,
                                     height=self.height,
                                     angle=angle,
                                     edgecolor=edgecolor,
                                     fill=fill))


def sector_mask(shape, center, radius, angle_range=(0, 360)):
    """
    Return a boolean mask for a circular sector. The start/stop angles in
    `angle_range` should be given in clockwise order.

    Found here: https://stackoverflow.com/questions/18352973/mask-a-circular-sector-in-a-numpy-array/18354475#18354475
    """

    x, y = np.ogrid[:shape[0], :shape[1]]
    cy, cx = center.x, center.y
    # tmin, tmax = np.deg2rad(angle_range)
    tmin, tmax = angle_range

    # ensure stop angle > start angle
    if tmax < tmin:
        tmax += 2 * np.pi

    # convert cartesian --> polar coordinates
    r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy)
    theta = np.arctan2(x - cx, y - cy) - tmin

    # wrap angles between 0 and 2*pi
    theta %= (2 * np.pi)

    # circular mask
    circmask = r2 <= radius * radius

    # angular mask
    anglemask = theta <= (tmax - tmin)

    return circmask * anglemask