

"""Test the various geometric patterns in the pylinac.core.geometry module."""

import unittest

from pylinac.core.geometry import *


class Test_Point(unittest.TestCase):

    x = 5
    y = -14
    z = 3.4
    idx = 14
    value = 13.28
    coord_iter = (x, y, z)
    all_vals_iter = (x, y, z, idx, value)

    def test_inputs(self):
        # create an empty point
        p = Point()

        # assert properties are set
        p = Point(x=self.x, y=self.y, z=self.z, idx=self.idx, value=self.value)
        self.assertEqual(p.x, self.x)
        self.assertEqual(p.y, self.y)
        self.assertEqual(p.idx, self.idx)
        self.assertEqual(p.value, self.value)

    def test_point_input(self):
        inp = Point(self.all_vals_iter)
        p = Point(inp)
        self.assertEqual(p.x, self.x)
        self.assertEqual(p.y, self.y)
        self.assertEqual(p.idx, self.idx)
        self.assertEqual(p.value, self.value)


    def test_iterable_input(self):
        """Test when an interable (list, tuple, etc) is passed in it is parsed out."""
        p = Point(self.coord_iter)
        self.assertEqual(p.x, self.x)
        self.assertEqual(p.y, self.y)
        self.assertEqual(p.z, self.z)

        p = Point(self.all_vals_iter)
        self.assertEqual(p.idx, self.idx)
        self.assertEqual(p.value, self.value)

    def test_as_int(self):
        """Test that when as_int is passed, the coords are converted to ints."""
        p = Point(self.x, self.y, self.z, as_int=True)
        self.assertIsInstance(p.z, int)

    def test_nonnumeric_value(self):
        p = Point()
        self.assertRaises(TypeError, p.value, 'not numeric')


class Test_Circle(unittest.TestCase):

    radius = 5.6
    center_point_Point = Point()
    center_point_iter = [3, 4]

    def test_inputs(self):
        # create default
        c = Circle()
        self.assertIsInstance(c.center, Point)
        self.assertIsNone(c.radius)

        # input Point class
        c = Circle(self.center_point_Point, self.radius)
        self.assertEqual(c.center.x, self.center_point_Point.x)
        self.assertEqual(c.radius, self.radius)

        # input iterable
        c = Circle(self.center_point_iter)
        self.assertEqual(c.center.x, self.center_point_iter[0])
        self.assertEqual(c.center.y, self.center_point_iter[1])

class Test_Line(unittest.TestCase):

    point_1 = Point(1,1)
    point_2 = Point(2,3)
    # the slope (m) of the two points above
    m_from_points = 2
    # ditto for intercept (b)
    b_from_points = -1
    # the y value when x=4 (2*4-1)
    y_at_4 = 7
    # ditto for x when y=4
    x_at_4 = 2.5

    # m*x+b inputs
    m = 3
    b = -1
    y_at_1 = 2
    x_at_5 = 2

    def test_inputs(self):

        # create from two points and test properties
        l = Line(self.point_1, self.point_2)
        self.assertEqual(l.point1.x, self.point_1.x)
        self.assertEqual(l.point2.y, self.point_2.y)
        self.assertEqual(l.m, self.m_from_points)
        self.assertEqual(l.b, self.b_from_points)
        self.assertEqual(l.y(4), self.y_at_4)
        self.assertEqual(l.x(4), self.x_at_4)

        # create from m*x+b
        l = Line(m=self.m, b=self.b)
        self.assertEqual(l.m, self.m)
        self.assertEqual(l.y(1), self.y_at_1)
        self.assertEqual(l.x(5), self.x_at_5)

    def test_dist2point(self):

        point = Point(0,3)
        line = Line((0,0), (1,1))
        exp_dist = 2.12

        self.assertAlmostEqual(line.distance_to(point), exp_dist, delta=0.01)