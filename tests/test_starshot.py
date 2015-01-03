import unittest
import os.path as osp

from pylinac.core.utilities import go_up_dirlevel
from pylinac.starshot import Starshot, StarProfile, Wobble
from pylinac.core.geometry import Point


test_file_dir = osp.join(osp.dirname(__file__), 'test_files', 'Starshot')

demo_file_path = osp.join(go_up_dirlevel(1), 'pylinac', 'demo_files', 'starshot', 'starshot_gantry.zip')


class Test_CircleProfile(unittest.TestCase):
    """Test the StarProfile class based on an analysis of demo image #1."""
    def setUp(self):
        self.cp = StarProfile()

    # def test_init(self):
    #     """Test that all the needed parameters are initialized."""
    #     needed_fields = ('profile', 'x', 'y', 'center', 'radius_pix', 'radius_perc')
    #     for field in needed_fields:
    #         self.assertTrue(hasattr(self.cp, field))
    #     self.assertIsInstance(self.cp.center, Point)


class Test_Wobble(unittest.TestCase):

    def setUp(self):
        self.wobble = Wobble()

    # def test_init(self):
    #     needed_fields = ('center', 'radius', 'radius_mm')
    #     for field in needed_fields:
    #         self.assertTrue(hasattr(self.wobble, field))
    #     self.assertIsInstance(self.wobble.center, Point)


class Star_general_tests(unittest.TestCase):
    """Performs general tests (not demo specific)."""
    def setUp(self):
        self.star = Starshot()
        self.star.load_demo_image()

    def test_startpoint_autosets_if_unset(self):
        """Test that the mechanical isocenter will automatically set if not yet set."""
        # analyze image with the start point not yet set
        self.star.analyze()
        # the start point should now have been set
        self.assertNotEqual(self.star.start_point.x, 0, msg="The start point did not set automatically when analyzing")

    def test_demo_file_deleted(self):
        """When using a demo image, the initial file is .zip. The algo extracts and should also delete
        the extracted file to save space.
        """
        extracted_file = demo_file_path.replace('.zip', '.tif')
        self.assertFalse(osp.exists(extracted_file))


class Test_Star_Demo(unittest.TestCase):
    """Specific tests for the demo image"""
    def setUp(self):
        self.star = Starshot()
        self.star.load_demo_image()

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
        self.assertAlmostEqual(self.star.wobble.radius_mm, 0.49, delta=0.1)

        # Test that the center of the wobble circle is close to what it's shown to be (1511.5, 1302.2).
        # test y-coordinate
        y_coord = self.star.wobble.center.y
        self.assertAlmostEqual(y_coord, 1511.5, delta=1)
        # test x-coordinate
        x_coord = self.star.wobble.center.x
        self.assertAlmostEqual(x_coord, 1301.1, delta=1)

        """Test than the number of lines found is what is expected."""
        expected_lines = 9
        self.assertEqual(len(self.star.lines), expected_lines,
                         msg="The number of lines found was not the number expected")

    @unittest.expectedFailure
    def test_pass_with_startpoint_off(self):
        """Test that the algo will still pass if the start point is set
            to a point somewhat different than actual center.
        """
        offset_point = Point(1300, 1500)
        self.star.set_start_point(offset_point, warn_if_far_away=False)
        self.star.analyze()

        self.assertTrue(self.star.passed)

        y_coord = self.star.wobble.center.y
        self.assertAlmostEqual(y_coord, 1511, delta=3)
        x_coord = self.star.wobble.center.x
        self.assertAlmostEqual(x_coord, 1302, delta=3)
        wobble_radius = self.star.wobble.radius_pix
        self.assertAlmostEqual(wobble_radius, 5, delta=0.5)


class Test_Star_Test_Files(unittest.TestCase):
    """Tests specifically for the starshot test files."""

    _test_file_names = ('6XCollStar.tif', '10XCollStar.bmp')
    file_paths = [osp.join(test_file_dir, test_file_name) for test_file_name in _test_file_names]

    def setUp(self):
        self.stars = []
        for pth in self.file_paths:
            star = Starshot()
            star.load_image(pth)
            self.stars.append(star)

    def test_passed(self):
        """Test that the demo image passed"""
        for star in self.stars:
            star.analyze()
            self.assertTrue(star.passed, msg="Wobble was not within tolerance")

    def test_wobble_radius(self):
        """Test than the wobble radius is similar to what it has been shown to be (0.495)."""
        wobble_radii = (5.07, 0.76)
        for star, wobble_radius in zip(self.stars, wobble_radii):
            star.analyze()
            self.assertAlmostEqual(star.wobble.radius, wobble_radius, delta=0.1)

    def test_wobble_center(self):
        """Test that the center of the wobble circle is close to what it's shown to be (1511.5, 1302.2)."""
        wobble_centers = (Point(1297,1695.6), Point(1368.6,1453))
        for star, wob_cent in zip(self.stars, wobble_centers):
            star.analyze()
            # test y-coordinate
            y_coord = star.wobble.center.y
            self.assertAlmostEqual(y_coord, wob_cent.y, delta=2)
            # test x-coordinate
            x_coord = star.wobble.center.x
            self.assertAlmostEqual(x_coord, wob_cent.x, delta=2)

    def test_num_rad_lines(self):
        """Test than the number of radiation lines found is what is expected."""
        num_expected_lines = (9, 4)
        for star, num_exp_lines in zip(self.stars, num_expected_lines):
            star.analyze()
            self.assertEqual(len(star.lines), num_exp_lines,
                             msg="The number of radiation lines found was not the number expected")


