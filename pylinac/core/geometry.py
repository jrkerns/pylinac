
"""Module for classes that represent common geometric objects or patterns."""
from math import sqrt

import numpy as np
from matplotlib.patches import Circle as mpl_Circle
from matplotlib.patches import Rectangle as mpl_Rectangle

from pylinac.core.utilities import is_iterable, typed_property


class Point:
    """A geometric point with x, y, and z coordinates/attributes."""
    value = typed_property('value', (int, float, np.number, type(None)))

    def __init__(self, x=0, y=0, idx=None, value=None, as_int=False):
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
            for attr in ['x', 'y', 'idx', 'value']:
                item = getattr(x, attr)
                setattr(self, attr, item)
        elif is_iterable(x):
            for attr, item in zip(['x', 'y', 'idx', 'value'], x):
                setattr(self, attr, item)
        else:
            self.x = x
            self.y = y
            self.idx = idx
            self.value = value

        if as_int:
            self.x = int(round(self.x))
            self.y = int(round(self.y))

    def dist_to(self, point):
        """Calculate the distance to the given point.

        Parameters
        ----------
        point : Point, 2 element iterable
            The other point to calculate distance to.
        """
        p = Point(point)
        return sqrt((self.x - p.x)**2 + (self.y - p.y)**2)


class Scale:
    """A 'scale' object with x and y attrs. Used in conjunction with scaling images up or down."""
    def __init__(self, x, y):
        self.x = x
        self.y = y


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
    def __init__(self, point1=None, point2=None, m=None, b=None):
        """
        Parameters
        ----------
        point1 : Point, optional
            One point of the line
        point2 : Point, optional
            Second point along the line.
        m : int, float, optional
            slope of the line (rise/run)
        b : int, float, optional
            y-intercept of the line
        """
        # if created by passing two points...
        if point1 is not None and point2 is not None:
            self.point1 = Point(point1)
            self.point2 = Point(point2)
        # otherwise by passing m and b...
        # elif m is not None and b is not None:
        #     self.m = m
        #     self.b = b

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
    def is_finite(self):
        """Boolean property specifying if the line is finite."""
        if self.point1 is not None and self.point2 is not None:
            return True
        else:
            return False

    @property
    def length(self):
        """Return length of the line, if finite."""
        if self.is_finite:
            return self.point1.dist_to(self.point2)
        else:
            raise ValueError("Line is not finite")

    def distance_to(self, point):
        """Calculate the minimum distance from the line to a point.

        Equations are from here: http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html #14

        Parameters
        ----------
        point : Point, iterable
            The point to calculate distance to.
        """
        point = Point(point)
        # calculate from m*x+b definition
        if self.point2 is None:
            #TODO: work on this
            raise NotImplementedError
        # calculate from 2 points definition
        else:
            lp1 = self.point1
            lp2 = self.point2
            numerator = np.abs((lp2.x - lp1.x)*(lp1.y - point.y) - (lp1.x - point.x)*(lp2.y - lp1.y))
            denominator = np.sqrt((lp2.x - lp1.x)**2 + (lp2.y - lp1.y)**2)
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
    # @type_accept(center=Point, tl_corner=Point, bl_corner=Point)
    def __init__(self, width, height, center=None, tl_corner=None, bl_corner=None, as_int=False):
        """
        Parameters
        ----------
        width : number
            Width of the rectangle.
        height : number
            Height of the rectangle.
        center : Point, iterable, optional
            Center point of rectangle.
        tl_corner : Point, iterable, optional
            Top-Left corner of the rectangle.
        bl_corner : Point, iterable, optional
            Bottom-Left corner of the rectangle.
        as_int : bool
            If False (default), inputs are left as-is. If True, all inputs are converted to integers.
        """
        if as_int:
            self.width = int(np.round(width))
            self.height = int(np.round(height))
        else:
            self.width = width
            self.height = height

        if not any((center, tl_corner, bl_corner)):
            raise ValueError("Must specify at least one anchor point for the box.")
        elif center is not None:
            c = self.center = Point(center, as_int=as_int)
            self.tl_corner = Point(c.x - width/2, c.y + height/2, as_int=as_int)
            self.bl_corner = Point(c.x - width/2, c.y - height/2, as_int=as_int)
        elif tl_corner is not None:
            tl = self.tl_corner = Point(tl_corner, as_int=as_int)
            self.center = Point(tl.x + width/2, tl.y - height/2, as_int=as_int)
            self.bl_corner = Point(tl.x, tl.y - height, as_int=as_int)
        elif bl_corner is not None:
            bl = self.bl_corner = Point(bl_corner, as_int=as_int)
            self.center = Point(bl.x + width / 2, bl.y + height / 2, as_int=as_int)
            self.tl_corner = Point(bl.x, bl.y + height, as_int=as_int)

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