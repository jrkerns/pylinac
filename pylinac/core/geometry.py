"""Module for classes that represent common geometric objects or patterns."""
from __future__ import annotations

import math
from itertools import zip_longest
from typing import Iterable

import argue
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle as mpl_Circle
from matplotlib.patches import Rectangle as mpl_Rectangle
from mpl_toolkits.mplot3d.art3d import Line3D

from .utilities import is_iterable


def tan(degrees: float) -> float:
    """Calculate the tangent of the given degrees."""
    return math.tan(math.radians(degrees))


def atan(x: float, y: float) -> float:
    """Calculate the degrees of a given x/y from the origin"""
    return math.degrees(math.atan2(x, y))


def cos(degrees: float) -> float:
    """Calculate the cosine of the given degrees."""
    return math.cos(math.radians(degrees))


def sin(degrees: float) -> float:
    """Calculate the sine of the given degrees."""
    return math.sin(math.radians(degrees))


def direction_to_coords(
    start_x: float, start_y: float, distance: float, angle_degrees: float
) -> (float, float):
    """Calculate destination coordinates given a start position, distance, and angle.

    The 0-angle position is pointing to the right (i.e. unit circle)

    Parameters
    ----------
    start_x : float
        Starting x position
    start_y : float
        Starting y position
    distance : float
        Distance to travel
    angle_degrees : float
        Angle to travel in degrees.
    """
    # Convert angle from degrees to radians
    angle_radians = math.radians(angle_degrees)

    # Calculate destination coordinates
    end_x = start_x + distance * math.cos(angle_radians)
    end_y = start_y + distance * math.sin(angle_radians)
    return end_x, end_y


class Point:
    """A geometric point with x, y, and z coordinates/attributes."""

    z: float
    y: float
    x: float
    _attr_list: list[str] = ["x", "y", "z", "idx", "value"]
    _coord_list: list[str] = ["x", "y", "z"]

    def __init__(
        self,
        x: float | tuple | Point = 0,
        y: float = 0,
        z: float = 0,
        idx: int | None = None,
        value: float | None = None,
        as_int: bool = False,
    ):
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

    def distance_to(self, thing: Point | Circle) -> float:
        """Calculate the distance to the given point.

        Parameters
        ----------
        thing : Circle, Point, 2 element iterable
            The other thing to calculate distance to.
        """
        if isinstance(thing, Circle):
            return abs(
                np.sqrt((self.x - thing.center.x) ** 2 + (self.y - thing.center.y) ** 2)
                - thing.radius
            )
        p = Point(thing)
        return math.sqrt(
            (self.x - p.x) ** 2 + (self.y - p.y) ** 2 + (self.z - p.z) ** 2
        )

    def as_array(self, only_coords: bool = True) -> np.ndarray:
        """Return the point as a numpy array."""
        if only_coords:
            return np.array([getattr(self, item) for item in self._coord_list])
        else:
            return np.array(
                [
                    getattr(self, item)
                    for item in self._attr_list
                    if (getattr(self, item) is not None)
                ]
            )

    def as_vector(self) -> Vector:
        return Vector(x=self.x, y=self.y, z=self.z)

    def __repr__(self) -> str:
        return f"Point(x={self.x:3.2f}, y={self.y:3.2f}, z={self.z:3.2f})"

    def __eq__(self, other) -> bool:
        # if all attrs equal, points considered equal
        return all(
            getattr(self, attr) == getattr(other, attr) for attr in self._attr_list
        )

    def __add__(self, other) -> Vector:
        p = Vector()
        for attr in self._attr_list:
            try:
                sum = getattr(self, attr) + getattr(other, attr)
            except TypeError:
                sum = None
            setattr(p, attr, sum)
        return p

    def __sub__(self, other) -> Vector:
        p = Vector()
        for attr in self._attr_list:
            try:
                diff = getattr(self, attr) - getattr(other, attr)
            except TypeError:
                diff = None
            setattr(p, attr, diff)
        return p

    def __mul__(self, other: int | float) -> None:
        for attr in self._attr_list:
            try:
                self.__dict__[attr] *= other
            except TypeError:
                pass

    def __truediv__(self, other: int | float) -> Point:
        for attr in self._attr_list:
            val = getattr(self, attr)
            # sometimes not all attrs are defined (like index or value and only x,y,z). Skip dividing those.
            if val is not None:
                setattr(self, attr, val / other)
        return self


class Circle:
    """A geometric circle with center Point, radius, and diameter."""

    center: Point
    radius: float

    def __init__(self, center_point: Point | Iterable = (0, 0), radius: float = 0):
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
    def area(self) -> float:
        """The area of the circle."""
        return math.pi * self.radius**2

    @property
    def diameter(self) -> float:
        """Get the diameter of the circle."""
        return self.radius * 2

    def plot2axes(
        self,
        axes: plt.Axes,
        edgecolor: str = "black",
        fill: bool = False,
        text: str = "",
        fontsize: str = "medium",
        **kwargs,
    ) -> None:
        """Plot the Circle on the axes.

        Parameters
        ----------
        axes : matplotlib.axes.Axes
            An MPL axes to plot to.
        edgecolor : str
            The color of the circle.
        fill : bool
            Whether to fill the circle with color or leave hollow.
        text: str
            If provided, plots the given text at the center. Useful for identifying ROIs on a plotted image apart.
        fontsize: str
            The size of the text, if provided. See https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
            for options.
        """
        axes.add_patch(
            mpl_Circle(
                (self.center.x, self.center.y),
                edgecolor=edgecolor,
                radius=self.radius,
                fill=fill,
                **kwargs,
            )
        )
        if text:
            axes.text(
                x=self.center.x,
                y=self.center.y,
                s=text,
                fontsize=fontsize,
                color=edgecolor,
            )

    def as_dict(self) -> dict:
        """Convert to dict. Useful for dataclasses/Result"""
        return {
            "center_x": self.center.x,
            "center_y": self.center.y,
            "diameter": self.diameter,
        }


class Vector:
    """A vector with x, y, and z coordinates."""

    x: float
    y: float
    z: float

    def __init__(self, x: float = 0, y: float = 0, z: float = 0):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return f"Vector(x={self.x:.2f}, y={self.y:.2f}, z={self.z:.2f})"

    def as_scalar(self) -> float:
        """Return the scalar equivalent of the vector."""
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def distance_to(self, thing: Circle | Point) -> float:
        """Calculate the distance to the given point.

        Parameters
        ----------
        thing : Circle, Point, 2 element iterable
            The other point to calculate distance to.
        """
        if isinstance(thing, Circle):
            return abs(
                np.sqrt((self.x - thing.center.x) ** 2 + (self.y - thing.center.y) ** 2)
                - thing.radius
            )
        else:
            p = Point(thing)
            return math.sqrt(
                (self.x - p.x) ** 2 + (self.y - p.y) ** 2 + (self.z - p.z) ** 2
            )

    def __sub__(self, other: Vector) -> Vector:
        new_x = self.x - other.x
        new_y = self.y - other.y
        new_z = self.z - other.z
        return Vector(x=new_x, y=new_y, z=new_z)

    def __add__(self, other: Vector) -> Vector:
        new_x = self.x + other.x
        new_y = self.y + other.y
        new_z = self.z + other.z
        return Vector(x=new_x, y=new_y, z=new_z)

    def __truediv__(self, other: float) -> Vector:
        for attr in ("x", "y", "z"):
            val = getattr(self, attr)
            setattr(self, attr, val / other)
        return self


def vector_is_close(vector1: Vector, vector2: Vector, delta: float = 0.1) -> bool:
    """Determine if two vectors are with delta of each other; this is a simple coordinate comparison check."""
    for attr in ("x", "y", "z"):
        if np.isnan(getattr(vector1, attr)) and np.isnan(getattr(vector2, attr)):
            continue
        if (
            not getattr(vector2, attr) + delta
            >= getattr(vector1, attr)
            >= getattr(vector2, attr) - delta
        ):
            return False
    return True


class Line:
    """A line that is represented by two points or by m*x+b.

    Notes
    -----
    Calculations of slope, etc are from here:
    http://en.wikipedia.org/wiki/Linear_equation
    and here:
    http://www.mathsisfun.com/algebra/line-equation-2points.html
    """

    point1: Point
    point2: Point

    def __init__(
        self,
        point1: Point | tuple[float, float],
        point2: Point | tuple[float, float],
    ):
        """
        Parameters
        ----------
        point1 : Point
            One point of the line
        point2 : Point
            Second point along the line.
        """
        self.point1 = Point(point1)
        self.point2 = Point(point2)

    def __repr__(self) -> str:
        return (
            f"Line: p1:(x={self.point1.x:.1f}, y={self.point1.y:.1f}, z={self.point1.z:.1f}), "
            f"p2:(x={self.point2.x:.1f}, y={self.point2.y:.1f}, z={self.point2.z:.1f})"
        )

    @property
    def m(self) -> float:
        """Return the slope of the line.

        m = (y1 - y2)/(x1 - x2)

        From: http://www.purplemath.com/modules/slope.htm
        """
        return (self.point1.y - self.point2.y) / (self.point1.x - self.point2.x)

    @property
    def b(self) -> float:
        """Return the y-intercept of the line.

        b = y - m*x
        """
        return self.point1.y - (self.m * self.point1.x)

    def y(self, x) -> float:
        """Return y-value along line at position x."""
        return self.m * x + self.b

    def x(self, y) -> float:
        """Return x-value along line at position y."""
        return (y - self.b) / self.m

    @property
    def center(self) -> Point:
        """Return the center of the line as a Point."""
        mid_x = np.abs((self.point2.x - self.point1.x) / 2 + self.point1.x)
        mid_y = (self.point2.y - self.point1.y) / 2 + self.point1.y
        return Point(mid_x, mid_y)

    @property
    def length(self) -> float:
        """Return length of the line, if finite."""
        return self.point1.distance_to(self.point2)

    def distance_to(self, point: Point) -> float:
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
        numerator = np.sqrt(np.sum(np.power(np.cross((lp2 - lp1), (lp1 - point)), 2)))
        denominator = np.sqrt(np.sum(np.power(lp2 - lp1, 2)))
        return numerator / denominator

    def plot2axes(
        self, axes: plt.Axes, width: float = 1, color: str = "w", **kwargs
    ) -> Line3D:
        """Plot the line to an axes.

        Parameters
        ----------
        axes : matplotlib.axes.Axes
            An MPL axes to plot to.
        color : str
            The color of the line.
        """
        lines = axes.plot(
            (self.point1.x, self.point2.x),
            (self.point1.y, self.point2.y),
            (self.point1.z, self.point2.z),
            linewidth=width,
            color=color,
            **kwargs,
        )
        return lines[0]


class Rectangle:
    """A rectangle with width, height, center Point, top-left corner Point, and bottom-left corner Point."""

    width: int | float
    height: int | float
    _as_int: bool
    center: Point

    def __init__(
        self,
        width: float,
        height: float,
        center: Point | tuple,
        as_int: bool = False,
    ):
        """
        Parameters
        ----------
        width : number
            Width of the rectangle. Must be positive
        height : number
            Height of the rectangle. Must be positive.
        center : Point, iterable, optional
            Center point of rectangle.
        as_int : bool
            If False (default), inputs are left as-is. If True, all inputs are converted to integers.
        """
        argue.verify_bounds(width, argue.POSITIVE)
        argue.verify_bounds(height, argue.POSITIVE)
        if as_int:
            self.width = int(np.round(width))
            self.height = int(np.round(height))
        else:
            self.width = width
            self.height = height
        self._as_int = as_int
        self.center = Point(center, as_int=as_int)

    @property
    def area(self) -> float:
        """The area of the rectangle."""
        return self.width * self.height

    @property
    def br_corner(self) -> Point:
        """The location of the bottom right corner."""
        return Point(
            self.center.x + self.width / 2,
            self.center.y - self.height / 2,
            as_int=self._as_int,
        )

    @property
    def bl_corner(self) -> Point:
        """The location of the bottom left corner."""
        return Point(
            self.center.x - self.width / 2,
            self.center.y - self.height / 2,
            as_int=self._as_int,
        )

    @property
    def tl_corner(self) -> Point:
        """The location of the top left corner."""
        return Point(
            self.center.x - self.width / 2,
            self.center.y + self.height / 2,
            as_int=self._as_int,
        )

    @property
    def tr_corner(self) -> Point:
        """The location of the top right corner."""
        return Point(
            self.center.x + self.width / 2,
            self.center.y + self.height / 2,
            as_int=self._as_int,
        )

    def plot2axes(
        self,
        axes: plt.Axes,
        edgecolor: str = "black",
        angle: float = 0.0,
        fill: bool = False,
        alpha: float = 1,
        facecolor: str = "g",
        label=None,
        text: str = "",
        fontsize: str = "medium",
        text_rotation: float = 0,
        **kwargs,
    ):
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
        text: str
            If provided, plots the given text at the center. Useful for identifying ROIs on a plotted image apart.
        fontsize: str
            The size of the text, if provided. See https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
            for options.
        text_rotation: float
            The rotation of the text in degrees.
        """
        axes.add_patch(
            mpl_Rectangle(
                (self.bl_corner.x, self.bl_corner.y),
                width=self.width,
                height=self.height,
                angle=angle,
                edgecolor=edgecolor,
                alpha=alpha,
                facecolor=facecolor,
                fill=fill,
                label=label,
                **kwargs,
            )
        )
        if text:
            axes.text(
                x=self.center.x,
                y=self.center.y,
                s=text,
                fontsize=fontsize,
                color=edgecolor,
                rotation=text_rotation,
                horizontalalignment="center",
                verticalalignment="center",
            )
