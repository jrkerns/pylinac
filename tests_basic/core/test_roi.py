from unittest import TestCase

import numpy as np

from pylinac.core.geometry import Point
from pylinac.core.roi import LowContrastDiskROI, RectangleROI


class TestDiskROI(TestCase):
    def test_array_size(self):
        disk = LowContrastDiskROI(
            np.ones((500, 500)),
            radius=20,
            center=Point(250, 250),
        )
        # it's a 1D array
        self.assertEqual(disk.pixel_values.shape[0], 1245)

    def test_from_phantom_center(self):
        disk = LowContrastDiskROI.from_phantom_center(
            np.ones((500, 500)),
            roi_radius=20,
            angle=0,
            dist_from_center=10,
            phantom_center=Point(250, 250),
        )
        self.assertEqual(disk.center.x, 260)
        self.assertEqual(disk.center.y, 250)


class TestRectangleROI(TestCase):
    def test_array_shape(self):
        rect = RectangleROI(
            np.ones((500, 500)),
            width=20,
            height=50,
            center=Point(250, 250),
        )
        self.assertEqual(rect.pixel_array.shape[0], 50)  # rows
        self.assertEqual(rect.pixel_array.shape[1], 20)  # cols
        self.assertEqual(rect.center.x, 250)
        self.assertEqual(rect.center.y, 250)

    def test_less_than_2x2_roi_not_allowed(self):
        with self.assertRaises(ValueError):
            RectangleROI(
                np.ones((500, 500)),
                width=0,
                height=50,
                center=Point(250, 250),
            )
        with self.assertRaises(ValueError):
            RectangleROI(
                np.ones((500, 500)),
                width=50,
                height=1,
                center=Point(250, 250),
            )

    def test_from_phantom_center(self):
        rect = RectangleROI.from_phantom_center(
            np.ones((500, 500)),
            width=20,
            height=50,
            angle=0,
            dist_from_center=10,
            phantom_center=Point(250, 250),
        )
        self.assertEqual(rect.pixel_array.shape[0], 50)  # rows
        self.assertEqual(rect.pixel_array.shape[1], 20)  # cols
        # center is shifted by 10 in x (angle=0)
        self.assertEqual(rect.center.x, 260)
        self.assertEqual(rect.center.y, 250)

    def test_rotation(self):
        # we add a spike near a corner. After we rotate, the spike should
        # no longer be in view of the rotated rectangle.
        array_with_spike = np.ones((100, 100))
        array_with_spike[33, 33] = 1000
        rect_rotated = RectangleROI(
            array_with_spike,
            width=40,
            height=40,
            center=Point(50, 50),
            rotation=45,
        )
        self.assertEqual(rect_rotated.max, 1)  # spike not in view
        # but an un-rotated rectangle should see the spike
        rect_unrotated = RectangleROI(
            array_with_spike,
            width=40,
            height=40,
            center=Point(50, 50),
        )
        self.assertEqual(rect_unrotated.max, 1000)  # spike in view

    def test_flat_pixels(self):
        """If there is no rotation, the statistics of the flat pixels should be the same as the pixel_array."""
        # see note in RectangleROI.pixel_array docstring; even widths/heights are rounded up.
        array = np.ones((100, 100))
        rect = RectangleROI(
            array,
            width=40,
            height=20,
            center=Point(50, 50),
        )
        array_2d = rect.pixel_array
        self.assertEqual(array_2d.size, 800)
        array_flat = rect.pixels_flat
        self.assertEqual(array_flat.size, 800)


class TestRectangleStats(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        np.random.seed(12345)  # noqa: NPY002

        cls.random_array = np.random.rand(500, 500)  # noqa: NPY002

        cls.random_rect = RectangleROI(
            cls.random_array,
            width=500,
            height=500,
            center=Point(250, 250),
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
        disk = LowContrastDiskROI.from_phantom_center(
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
        disk = LowContrastDiskROI.from_phantom_center(
            array=homogeneous_array,
            angle=0,
            roi_radius=2,
            dist_from_center=0,
            phantom_center=Point(5, 5),
            contrast_reference=0,
        )
        self.assertEqual(disk.contrast_to_noise, float("inf"))
