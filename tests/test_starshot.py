import os
import os.path as osp
import unittest

import matplotlib.pyplot as plt
import numpy as np

from pylinac.core.geometry import Point
from pylinac.starshot import Starshot
from tests.utils import save_file

TEST_DIR = osp.join(osp.dirname(__file__), 'test_files', 'Starshot')


class GeneralTests(unittest.TestCase):
    """Performs general tests (not image specific)."""

    def setUp(self):
        self.star = Starshot.from_demo_image()
        self.star.analyze()

    def test_analyze_without_images(self):
        star = Starshot()
        with self.assertRaises(AttributeError):
            star.analyze()

    def test_save_analyzed_image(self):
        """Test that saving an image does something."""
        save_file(self.star.save_analyzed_image)

    def test_save_analyzed_subimage(self):
        # save as normal file
        save_file(self.star.save_analyzed_subimage)
        # save into buffer
        save_file(self.star.save_analyzed_subimage, as_file_object='b')

    def test_from_url(self):
        url = 'https://s3.amazonaws.com/assuranceqa-staging/uploads/imgs/10X_collimator_dvTK5Jc.jpg'
        Starshot.from_url(url)  # shouldn't raise


class StarMixin:
    """Mixin for testing a starshot image."""
    star_file = ''
    wobble_diameter_mm = 0
    wobble_center = Point()
    num_rad_lines = 0
    recursive = True
    passes = True
    min_peak_height = 0.25
    test_all_radii = True
    fwxm = True
    wobble_tolerance = 0.2

    @classmethod
    def setUpClass(cls):
        cls.star = Starshot(cls.star_file)
        cls.star.analyze(recursive=cls.recursive, min_peak_height=cls.min_peak_height, fwhm=cls.fwxm)

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def test_passed(self):
        """Test that the demo image passed"""
        self.star.analyze(recursive=self.recursive, min_peak_height=self.min_peak_height)
        self.assertEqual(self.star.passed, self.passes, msg="Wobble was not within tolerance")

    def test_wobble_diameter(self):
        """Test than the wobble radius is similar to what it has been shown to be)."""
        self.assertAlmostEqual(self.star.wobble.diameter_mm, self.wobble_diameter_mm, delta=self.wobble_tolerance)

    def test_wobble_center(self):
        """Test that the center of the wobble circle is close to what it's shown to be."""
        # test y-coordinate
        y_coord = self.star.wobble.center.y
        self.assertAlmostEqual(y_coord, self.wobble_center.y, delta=3)
        # test x-coordinate
        x_coord = self.star.wobble.center.x
        self.assertAlmostEqual(x_coord, self.wobble_center.x, delta=3)

    def test_num_rad_lines(self):
        """Test than the number of radiation lines found is what is expected."""
        self.assertEqual(len(self.star.lines), self.num_rad_lines,
                         msg="The number of radiation lines found was not the number expected")

    def test_all_radii_give_same_wobble(self):
        """Test that the wobble stays roughly the same for all radii."""
        if self.test_all_radii:
            star = Starshot(self.star_file)
            for radius in np.linspace(0.9, 0.25, 8):
                star.analyze(radius=float(radius), min_peak_height=self.min_peak_height, recursive=self.recursive)
                self.assertAlmostEqual(star.wobble.diameter_mm, self.wobble_diameter_mm, delta=self.wobble_tolerance)


class Demo(StarMixin, unittest.TestCase):
    """Specific tests for the demo image"""
    star_file = osp.join(osp.dirname(osp.dirname(__file__)), 'pylinac', 'demo_files', 'starshot', 'starshot.tif')
    wobble_diameter_mm = 0.30
    wobble_center = Point(1270, 1437)
    num_rad_lines = 4

    def test_fails_with_tight_tol(self):
        star = Starshot.from_demo_image()
        star.analyze(tolerance=0.1)
        self.assertFalse(star.passed_constant)

    def test_bad_inputs_still_recovers(self):
        self.star.analyze(radius=0.3, min_peak_height=0.1)
        self.test_wobble_center()
        self.test_wobble_diameter()

    def test_demo_runs(self):
        """Test that the demo runs without error."""
        self.star.run_demo()

    def test_image_inverted(self):
        """Check that the demo image was actually inverted, as it needs to be."""
        star = Starshot.from_demo_image()
        top_left_corner_val_before = star.image.array[0,0]
        star._check_image_inversion()
        top_left_corner_val_after = star.image.array[0,0]
        self.assertNotEqual(top_left_corner_val_before, top_left_corner_val_after)

    def test_SID_can_be_overridden_for_nonEPID(self):
        self.star.analyze(SID=400)
        self.assertNotEqual(self.star.wobble.diameter, self.wobble_diameter_mm*2)

    def test_bad_start_point_recovers(self):
        """Test that even at a distance start point, the search algorithm recovers."""
        self.star.analyze(start_point=(1000, 1000))
        self.test_passed()
        self.test_wobble_center()
        self.test_wobble_diameter()


class Multiples(StarMixin, unittest.TestCase):
    """Test a starshot composed of multiple individual EPID images."""
    num_rad_lines = 9
    wobble_center = Point(254, 192)
    wobble_diameter_mm = 0.8
    test_all_radii = False

    @classmethod
    def setUpClass(cls):
        img_dir = osp.join(TEST_DIR, 'set')
        img_files = [osp.join(img_dir, filename) for filename in os.listdir(img_dir)]
        cls.star = Starshot.from_multiple_images(img_files)
        cls.star.analyze(radius=0.6)


class Starshot1(StarMixin, unittest.TestCase):
    star_file = osp.join(TEST_DIR, 'Starshot#1.tif')
    wobble_center = Point(508, 683)
    wobble_diameter_mm = 0.23
    num_rad_lines = 4


class Starshot1FWHM(Starshot1):
    fwhm = False

