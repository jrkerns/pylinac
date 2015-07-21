import os
import os.path as osp
import unittest

import numpy as np

from pylinac.starshot import Starshot
from pylinac.core.geometry import Point
from tests.utils import save_file


test_file_dir = osp.join(osp.dirname(__file__), 'test_files', 'Starshot')


class Star_general_tests(unittest.TestCase):
    """Performs general tests (not image specific)."""
    def setUp(self):
        self.star = Starshot.from_demo_image()
        self.star.analyze()

    def test_analyze_without_images(self):
        star = Starshot()
        self.assertRaises(AttributeError, star.analyze)

    def test_save_image(self):
        """Test that saving an image does something."""
        save_file('test.jpg', self.star.save_analyzed_image)


class StarMixin:
    # class attrs should be overridden by inheritors
    star_file = ''
    wobble_diameter_mm = 0
    wobble_center = Point()
    num_rad_lines = 0
    recursive = True
    passes = True

    @classmethod
    def setUpClass(cls):
        cls.star = Starshot(cls.star_file)
        cls.star.analyze(recursive=cls.recursive)

    def test_passed(self):
        """Test that the demo image passed"""
        self.star.analyze()
        self.assertEqual(self.star.passed, self.passes, msg="Wobble was not within tolerance")

    def test_wobble_diameter(self):
        """Test than the wobble radius is similar to what it has been shown to be)."""
        self.assertAlmostEqual(self.star.wobble.diameter_mm, self.wobble_diameter_mm, delta=0.25)

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
        for radius in np.linspace(0.9, 0.25, 8):
            self.star.analyze(radius=float(radius))
            self.assertAlmostEqual(self.star.wobble.diameter_mm, self.wobble_diameter_mm, delta=0.15)
            # print(self.star.wobble.diameter_mm)


class Demo(StarMixin, unittest.TestCase):
    """Specific tests for the demo image"""
    wobble_diameter_mm = 0.3
    wobble_center = Point(1270, 1436)
    num_rad_lines = 4

    @classmethod
    def setUpClass(cls):
        cls.star = Starshot.from_demo_image()
        cls.star.analyze()

    def test_fails_with_tight_tol(self):
        self.star.analyze(tolerance=0.1)
        self.assertFalse(self.star.passed)

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


class Collimator(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, '6XCollStar.tif')
    wobble_center = Point(1297, 1699)
    wobble_diameter_mm = 0.32
    num_rad_lines = 9

    def test_not_fwhm_passes(self):
        self.star.analyze(fwhm=False)
        self.test_passed()


class Collimator2(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, '10XCollStar.bmp')
    wobble_center = Point(1370, 1454)
    wobble_diameter_mm = 0.3
    num_rad_lines = 4


class Multiples(StarMixin, unittest.TestCase):
    """Test a starshot composed of multiple individual files."""
    num_rad_lines = 9
    wobble_center = Point(254, 192)
    wobble_diameter_mm = 0.8

    @classmethod
    def setUpClass(cls):
        img_dir = osp.join(test_file_dir, 'set')
        img_files = [osp.join(img_dir, filename) for filename in os.listdir(img_dir)]
        cls.star = Starshot.from_multiple_images(img_files)
        cls.star.analyze(radius=0.6)

    def test_radius_larger_than_passed(self):
        """Test the outcome of an analysis where the passed radius is outside the edge of the radiation lines."""
        # with recursive recovers
        self.star.analyze(radius=0.9)
        self.test_passed()
        self.test_wobble_center()

    def test_all_radii_give_same_wobble(self):
        pass


class Gantry(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'starshot_gantry.tif')
    wobble_center = Point(1302, 1513)
    wobble_diameter_mm = 0.95
    num_rad_lines = 9

    def test_bad_input_no_recursion_fails(self):
        """Test that without recursion, a bad setup fails."""
        with self.assertRaises(RuntimeError):
            self.star.analyze(radius=0.3, min_peak_height=0.95, recursive=False)

        # but will pass with recursion
        self.star.analyze()
        self.test_passed()


class Couch(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'couch.tif')
    wobble_center = Point(732, 944)
    wobble_diameter_mm = 2.1
    num_rad_lines = 6
    passes = False


@unittest.skip
class Couch2(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'AnnualGantryStarshot.tif')
    wobble_center = Point(732, 944)
    wobble_diameter_mm = 2.1
    num_rad_lines = 6
    passes = False