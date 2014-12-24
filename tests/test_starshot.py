import unittest
import os.path as osp

import numpy as np

from pylinac.core.utilities import go_up_dirlevel
from pylinac.starshot import Starshot, StarProfile, Wobble
from pylinac.core.geometry import Point


test_dir = osp.dirname(__file__)
test_file_dir = osp.join(test_dir, 'test_files', 'Starshot')

demo_file_dir = osp.join(go_up_dirlevel(1), 'pylinac', 'demo_files', 'starshot')
demo_filepath = osp.join(demo_file_dir, 'starshot_gantry.zip')



class Test_CircleProfile(unittest.TestCase):
    """Test the StarProfile class based on an analysis of demo image #1."""
    def setUp(self):
        self.cp = StarProfile()
        # self.star = Starshot().load_demo_image(1)

    def test_init(self):
        """Test that all the needed parameters are initialized (in the event code changes)."""
        needed_fields = ('profile', 'x', 'y', 'center', 'radius_pix', 'radius_perc')
        for field in needed_fields:
            self.assertTrue(hasattr(self.cp, field))
        self.assertIsInstance(self.cp.center, Point)


class Test_Wobble(unittest.TestCase):

    def setUp(self):
        self.wobble = Wobble()

    def test_init(self):
        needed_fields = ('center', 'radius', 'radius_pix')
        for field in needed_fields:
            self.assertTrue(hasattr(self.wobble, field))
        self.assertIsInstance(self.wobble.center, Point)


class Star_general_tests(unittest.TestCase):
    """Performs general tests (not demo specific)."""
    def setUp(self):
        self.star = Starshot()
        self.star.load_demo_image(1)

    def test_startpoint_autosets_if_unset(self):
        """Test that the mechanical isocenter will automatically set if not yet set."""
        # analyze image with the start point not yet set
        self.star.analyze()
        # the start point should now have been set
        self.assertNotEqual(self.star._algo_startpoint.x, 0, msg="The start point did not set automatically when analyzing")

    def test_extracted_file_deleted(self):
        """When using a demo image, the initial file is .zip. The algo extracts and should also delete
        the extracted file to save space.
        """
        extracted_file = demo_1_filepath.replace('.zip', '.tif')
        self.assertFalse(osp.exists(extracted_file))




class Star_test_demo1(unittest.TestCase):
    """Specific tests for the first demo image"""
    def setUp(self):
        self.star = Starshot()
        self.star.load_demo_image(number=1)

    def test_passed(self):
        """Test that the demo image passed"""
        self.star.analyze()
        self.assertTrue(self.star.passed, msg="Wobble was not within tolerance")

    def test_image_inverted(self):
        """Check that demo 1 was actually inverted."""
        top_left_corner_val_before = self.star.image.pixel_array[0,0]
        self.star._check_image_inversion()
        top_left_corner_val_after = self.star.image.pixel_array[0,0]
        self.assertNotEqual(top_left_corner_val_before, top_left_corner_val_after)

    def test_results(self):
        """Test that the wobble radius is similar to what it has been shown to be (0.495)."""
        self.star.analyze()
        self.assertAlmostEqual(self.star.wobble.radius, 0.49, delta=0.1)

        # Test that the center of the wobble circle is close to what it's shown to be (1511.5, 1302.2).
        # test y-coordinate
        y_coord = self.star.wobble.center.y
        self.assertAlmostEqual(y_coord, 1511.5, delta=1)
        # test x-coordinate
        x_coord = self.star.wobble.center.x
        self.assertAlmostEqual(x_coord, 1302.1, delta=1)

        """Test than the number of lines found is what is expected."""
        expected_lines = 9
        self.assertEqual(len(self.star.lines), expected_lines,
                         msg="The number of lines found was not the number expected")

    def test_still_pass_with_startpoint_off(self):
        """Test that the algo will still pass if the start point is set
            to a point somewhat different than actual center.
        """
        offset_point = Point(1250, 1600)
        self.star.set_start_point(offset_point, warn_if_far_away=False)
        self.star.analyze()

        self.assertTrue(self.star.passed)

        y_coord = self.star.wobble.center.y
        x_coord = self.star.wobble.center.x
        wobble_radius = self.star.wobble.radius_pix
        self.assertAlmostEqual(y_coord, 1511, delta=2)
        self.assertAlmostEqual(x_coord, 1302, delta=2)
        self.assertAlmostEqual(wobble_radius, 5, delta=0.5)


class Star_test_demo2(unittest.TestCase):
    """Tests specifically for demo image #2."""
    def setUp(self):
        self.star = Starshot()
        self.star.load_demo_image(2)

    def test_image_is_numpy(self):
        """The demo image should be numpy array when loaded."""
        self.assertIsInstance(self.star.image.pixel_array, np.ndarray)

    def test_passed(self):
        """Test that the demo image passed"""
        self.star.analyze()
        self.assertTrue(self.star.passed, msg="Wobble was not within tolerance")

    def test_wobble_radius(self):
        """Test than the wobble radius is similar to what it has been shown to be (0.495)."""
        self.star.analyze()
        self.assertAlmostEqual(self.star.wobble.radius, 0.17, delta=0.1)

    def test_wobble_center(self):
        """Test that the center of the wobble circle is close to what it's shown to be (1511.5, 1302.2)."""
        self.star.analyze()
        # test y-coordinate
        y_coord = self.star.wobble.center.y
        self.assertAlmostEqual(y_coord, 1698, delta=2)
        # test x-coordinate
        x_coord = self.star.wobble.center.x
        self.assertAlmostEqual(x_coord, 1296, delta=2)

    def test_num_rad_lines(self):
        """Test than the number of radiation lines found is what is expected."""
        expected_lines = 9
        self.star.analyze()
        self.assertEqual(len(self.star.lines), expected_lines,
                         msg="The number of lines found was not the number expected")


