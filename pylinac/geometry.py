

"""Objects that represent common geometric objects or patterns."""
import numpy as np

from decorators import type_accept


class Line(object):
    """Model a line that is represented by two points."""

    def __init__(self, point1=None, point2=None, m=None, b=None, is_finite=False):
        """Create a line from *either* two distinct points, or an m*x+b definition.

        :param point1, point2: Points along the line
        :type point1, point2: Point
        """
        self.point1 = point1
        self.point2 = point2
        #TODO: add m*x+b property functionality

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
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

class Box(object):

    @type_accept(center=(Point, None), dtype=(float, int))
    def __init__(self, width, height, center=None, tl_corner=None, bl_corner=None, dtype=float):
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

        # if int type specified, converted all properties to ints
        if dtype == int:
            self.width = int(self.width)
            self.height = int(self.height)
            self.center = Point(int(self.center.x), int(self.center.y))
            self.bl_corner = Point(int(self.bl_corner.x), int(self.bl_corner.y))
            self.tl_corner = Point(int(self.tl_corner.x), int(self.tl_corner.y))
