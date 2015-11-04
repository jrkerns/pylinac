
"""Module for classes that represent common geometric objects or patterns."""
from itertools import zip_longest
import math

import numpy as np
from matplotlib.patches import Circle as mpl_Circle
from matplotlib.patches import Rectangle as mpl_Rectangle

from pylinac.core.utilities import is_iterable


def tan(degrees):
    return math.tan(math.radians(degrees))


def cos(degrees):
    return math.cos(math.radians(degrees))


def sin(degrees):
    return math.sin(math.radians(degrees))


class Vector:
    """A vector with x, y, and z coordinates."""
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return "Vector(x={:.2f}, y={:.2f}, z={:.2f})".format(self.x, self.y, self.z)

    def as_scalar(self):
        """Return the scalar equivalent of the vector."""
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def distance_to(self, thing):
        """Calculate the distance to the given point.

        Parameters
        ----------
        thing : Circle, Point, 2 element iterable
            The other point to calculate distance to.
        """
        if isinstance(thing, Circle):
            return abs(np.sqrt((self.x - thing.center.x)**2 + (self.y - thing.center.y)**2) - thing.radius)
        else:
            p = Point(thing)
            return math.sqrt((self.x - p.x)**2 + (self.y - p.y)**2 + (self.z - p.z)**2)


def vector_is_close(vector1, vector2, delta=0.1):
    """Determine if two vectors are with delta of each other; this is a simple coordinate comparison check."""
    for attr in ('x', 'y', 'z'):
        if not getattr(vector2, attr) + delta >= getattr(vector1, attr) >= getattr(vector2, attr) - delta:
            return False
    return True


class Point:
    """A geometric point with x, y, and z coordinates/attributes."""
    _attr_list = ['x', 'y', 'z', 'idx', 'value']
    _coord_list = ['x', 'y', 'z']

    def __init__(self, x=0, y=0, z=0, idx=None, value=None, as_int=False):
        """
        Parameters
        ----------
        x : number-like, Point, iterable
            x-coordinate or iterable type containing all coordinates. If iterable, values are assumed to be in order: (x,y,z).
        y : number-like, optional
            y-coordinate
        idx : int, optional
            Index of point. Useful for sequential coordinates; e.g. a point on a circle profile is sometimes easier to describe
            in terms of its index rather than x,y coords.
        value : number-like, optional
            value at point location (e.g. pixel value of an image)
        as_int : boolean
            If True, coordinates are converted to integers.
        """
        if isinstance(x, Point):
            for attr in self._attr_list:
                item = getattr(x, attr, None)
                setattr(self, attr, item)
        elif is_iterable(x):
            for attr, item in zip_longest(self._attr_list, x, fillvalue=0):
                setattr(self, attr, item)
        else:
            self.x = x
            self.y = y
            self.z = z
            self.idx = idx
            self.value = value

        if as_int:
            self.x = int(round(self.x))
            self.y = int(round(self.y))
            self.z = int(round(self.z))

    def distance_to(self, thing):
        """Calculate the distance to the given point.

        Parameters
        ----------
        thing : Circle, Point, 2 element iterable
            The other thing to calculate distance to.
        """
        if isinstance(thing, Circle):
            return abs(np.sqrt((self.x - thing.center.x)**2 + (self.y - thing.center.y)**2) - thing.radius)
        p = Point(thing)
        return math.sqrt((self.x - p.x)**2 + (self.y - p.y)**2 + (self.z - p.z)**2)

    def as_array(self, only_coords=True):
        """Return the point as a numpy array."""
        if only_coords:
            return np.array([getattr(self, item) for item in self._coord_list])
        else:
            return np.array([getattr(self, item) for item in self._attr_list if (getattr(self, item) is not None)])

    def __repr__(self):
        return "Point(x={:3.2f}, y={:3.2f}, z={:3.2f})".format(self.x, self.y, self.z)

    def __eq__(self, other):
        # if all attrs equal, points considered equal
        return all(getattr(self, attr) == getattr(other, attr) for attr in self._attr_list)

    def __sub__(self, other):
        p = Point()
        for attr in self._attr_list:
            try:
                diff = getattr(self, attr) - getattr(other, attr)
            except TypeError:
                diff = None
            setattr(p, attr, diff)
        return p

    def __mul__(self, other):
        for attr in self._attr_list:
            try:
                self.__dict__[attr] *= other
            except TypeError:
                pass

    def __truediv__(self, other):
        for attr in self._attr_list:
            val = getattr(self, attr)
            try:
                setattr(self, attr, val/other)
            except TypeError:
                pass
        return self


class Circle:
    """A geometric circle with center Point, radius, and diameter."""

    def __init__(self, center_point=None, radius=None):
        """
        Parameters
        ----------
        center_point : Point, optional
            Center point of the wobble circle.
        radius : float, optional
            Radius of the wobble circle.
        """
        if center_point is None:
            center_point = Point()
        elif isinstance(center_point, Point) or is_iterable(center_point):
            center_point = Point(center_point)
        else:
            raise TypeError("Circle center must be of type Point or iterable")

        self.center = center_point
        self.radius = radius

    @property
    def diameter(self):
        """Get the diameter of the circle."""
        return self.radius*2

    def add_to_axes(self, axes, edgecolor='black', fill=False):
        """Plot the Circle on the axes.

        Parameters
        ----------
        axes : matplotlib.axes.Axes
            An MPL axes to plot to.
        edgecolor : str
            The color of the circle.
        fill : bool
            Whether to fill the circle with color or leave hollow.
        """
        axes.add_patch(mpl_Circle((self.center.x, self.center.y), edgecolor=edgecolor, radius=self.radius, fill=fill))


class Line:
    """A line that is represented by two points or by m*x+b.

    Notes
    -----
    Calculations of slope, etc are from here:
    http://en.wikipedia.org/wiki/Linear_equation
    and here:
    http://www.mathsisfun.com/algebra/line-equation-2points.html
    """
    def __init__(self, point1, point2):
        """
        Parameters
        ----------
        point1 : Point, optional
            One point of the line
        point2 : Point, optional
            Second point along the line.
        """
        self.point1 = Point(point1)
        self.point2 = Point(point2)

    @property
    def m(self):
        """Return the slope of the line.

        m = (y1 - y2)/(x1 - x2)

        From: http://www.purplemath.com/modules/slope.htm
        """
        return (self.point1.y - self.point2.y) / (self.point1.x - self.point2.x)

    @property
    def b(self):
        """Return the y-intercept of the line.

        b = y - m*x
        """
        return self.point1.y - (self.m * self.point1.x)

    def y(self, x):
        """Return y-value along line at position x."""
        return self.m * x + self.b

    def x(self, y):
        """Return x-value along line at position y."""
        return (y - self.b)/self.m

    @property
    def center(self):
        mid_x = np.abs((self.point2.x - self.point1.x)/2 + self.point1.x)
        mid_y = (self.point2.y - self.point1.y) / 2 + self.point1.y
        return Point(mid_x, mid_y)

    @property
    def length(self):
        """Return length of the line, if finite."""
        return self.point1.distance_to(self.point2)

    def distance_to(self, point):
        """Calculate the minimum distance from the line to a point.

        Equations are from here: http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html #14

        Parameters
        ----------
        point : Point, iterable
            The point to calculate distance to.
        """
        point = Point(point).as_array()
        lp1 = self.point1.as_array()
        lp2 = self.point2.as_array()
        # numerator = np.abs((lp2.x - lp1.x)*(lp1.y - point.y) - (lp1.x - point.x)*(lp2.y - lp1.y))
        # denominator = np.sqrt((lp2.x - lp1.x)**2 + (lp2.y - lp1.y)**2)
        numerator = sum(abs(x) for x in np.cross((lp2 - lp1), (lp1 - point)))
        denominator = sum(abs(x) for x in lp2 - lp1)
        return numerator/denominator

    def add_to_axes(self, axes, width=1, color='w'):
        """Plot the line to an axes.

        Parameters
        ----------
        axes : matplotlib.axes.Axes
            An MPL axes to plot to.
        color : str
            The color of the line.
        """
        axes.plot((self.point1.x, self.point2.x), (self.point1.y, self.point2.y), linewidth=width, color=color)


class Rectangle:
    """A rectangle with width, height, center Point, top-left corner Point, and bottom-left corner Point."""

    def __init__(self, width, height, center, as_int=False):
        """
        Parameters
        ----------
        width : number
            Width of the rectangle.
        height : number
            Height of the rectangle.
        center : Point, iterable, optional
            Center point of rectangle.
        as_int : bool
            If False (default), inputs are left as-is. If True, all inputs are converted to integers.
        """
        if as_int:
            self.width = int(np.round(width))
            self.height = int(np.round(height))
        else:
            self.width = width
            self.height = height
        self._as_int = as_int
        self.center = Point(center, as_int=as_int)

    @property
    def br_corner(self):
        return Point(self.center.x + self.width / 2, self.center.y - self.height / 2, as_int=self._as_int)

    @property
    def bl_corner(self):
        return Point(self.center.x - self.width / 2, self.center.y - self.height / 2, as_int=self._as_int)

    @property
    def tl_corner(self):
        return Point(self.center.x - self.width / 2, self.center.y + self.height / 2, as_int=self._as_int)

    @property
    def tr_corner(self):
        return Point(self.center.x + self.width / 2, self.center.y + self.height / 2, as_int=self._as_int)

    def add_to_axes(self, axes, edgecolor='black', angle=0.0, fill=False, alpha=1, facecolor='g'):
        """Plot the Rectangle to the axes.

        Parameters
        ----------
        axes : matplotlib.axes.Axes
            An MPL axes to plot to.
        edgecolor : str
            The color of the circle.
        angle : float
            Angle of the rectangle.
        fill : bool
            Whether to fill the rectangle with color or leave hollow.
        """
        axes.add_patch(mpl_Rectangle((self.bl_corner.x, self.bl_corner.y),
                                     width=self.width,
                                     height=self.height,
                                     angle=angle,
                                     edgecolor=edgecolor,
                                     alpha=alpha,
                                     facecolor=facecolor,
                                     fill=fill))


def sector_mask(shape, center, radius, angle_range=(0, 360)):
    """Return a circular arc-shaped boolean mask.

    Parameters
    ----------
    shape : tuple
        Shape of the image matrix. Usually easiest to pass something like array.shape
    center : Point, iterable
        The center point of the desired mask.
    radius : int, float
        Radius of the mask.
    angle_range : iterable
        The angle range of the mask. E.g. the default (0, 360) creates an entire circle.
        The start/stop angles should be given in clockwise order. 0 is right (0 on unit circle).

    References
    ----------
    https://stackoverflow.com/questions/18352973/mask-a-circular-sector-in-a-numpy-array/18354475#18354475
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