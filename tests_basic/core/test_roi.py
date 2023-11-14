from unittest import TestCase

import numpy as np

from pylinac.core.geometry import Point
from pylinac.core.roi import LowContrastDiskROI, RectangleROI


class TestRectangleROI(TestCase):
    def test_array_shape(self):
        rect = RectangleROI(
            np.ones((500, 500)),
            width=20,
            height=50,
            angle=0,
            dist_from_center=0,
            phantom_center=Point(250, 250),
        )
        self.assertEqual(rect.pixel_array.shape[0], 50)  # rows
        self.assertEqual(rect.pixel_array.shape[1], 20)  # cols

    def test_less_than_2x2_roi_not_allowed(self):
        with self.assertRaises(ValueError):
            RectangleROI(
                np.ones((500, 500)),
                width=0,
                height=50,
                angle=0,
                dist_from_center=0,
                phantom_center=Point(250, 250),
            )
        with self.assertRaises(ValueError):
            RectangleROI(
                np.ones((500, 500)),
                width=50,
                height=1,
                angle=0,
                dist_from_center=0,
                phantom_center=Point(250, 250),
            )


class TestRectangleStats(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        np.random.seed(12345)

        cls.random_array = np.random.rand(500, 500)

        cls.random_rect = RectangleROI(
            cls.random_array,
            width=500,
            height=500,
            angle=0,
            dist_from_center=0,
            phantom_center=Point(250, 250),
        )

    def test_mean(self):
        self.assertAlmostEqual(np.mean(self.random_array), self.random_rect.mean)

    def test_max(self):
        self.assertAlmostEqual(np.max(self.random_array), self.random_rect.max)

    def test_min(self):
        self.assertAlmostEqual(np.min(self.random_array), self.random_rect.min)

    def test_std(self):
        self.assertAlmostEqual(np.std(self.random_array), self.random_rect.std)


class TestEdgeCases(TestCase):
    def test_0_std_snr(self):
        # this will result in an infinite SNR; should be handled by numpy.
        # potentially need to address later to handle better; unsure
        homogeneous_array = np.ones((10, 10))
        disk = LowContrastDiskROI(
            array=homogeneous_array,
            angle=0,
            roi_radius=2,
            dist_from_center=0,
            phantom_center=Point(5, 5),
        )
        self.assertEqual(disk.signal_to_noise, float("inf"))

    def test_0_std_cnr(self):
        # this will result in an infinite SNR; should be handled by numpy.
        # potentially need to address later to handle better; unsure
        homogeneous_array = np.ones((10, 10))
        disk = LowContrastDiskROI(
            array=homogeneous_array,
            angle=0,
            roi_radius=2,
            dist_from_center=0,
            phantom_center=Point(5, 5),
            contrast_reference=0,
        )
        self.assertEqual(disk.contrast_to_noise, float("inf"))
