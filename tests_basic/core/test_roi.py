from unittest import TestCase

import numpy as np

from pylinac.core.geometry import Point
from pylinac.core.roi import RectangleROI


class TestRectangle(TestCase):

    def test_array_shape(self):
        rect = RectangleROI(np.ones((500, 500)), width=20, height=50, angle=0,
                            dist_from_center=0, phantom_center=Point(250, 250))
        self.assertEqual(rect.pixel_array.shape[0], 50)  # rows
        self.assertEqual(rect.pixel_array.shape[1], 20)  # cols
