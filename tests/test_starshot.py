import unittest
import os.path as osp

from pylinac.starshot import Starshot
from pylinac.core.geometry import Point


test_file_dir = osp.join(osp.dirname(__file__), 'test_files', 'Starshot')


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

    def test_analyze_without_images(self):
        star = Starshot()
        self.assertRaises(AttributeError, star.analyze)


class Test_Star_Demo(unittest.TestCase):
    """Specific tests for the demo image"""
    demo_radius = 0.22
    demo_x = 1269
    demo_y = 1439
    demo_lines = 4

    def setUp(self):
        self.star = Starshot()
        self.star.load_demo_image()

    def test_demo_runs(self):
        """Test that the demo runs without error."""
        # TODO: come up with decorator that adds show parameter
        self.star.run_demo(show=False)

    def test_passed(self):
        """Test that the demo image passed"""
        self.star.analyze()
        self.assertTrue(self.star.passed, msg="Wobble was not within tolerance")

    def test_image_inverted(self):
        """Check that the demo image was actually inverted, as it needs to be."""
        top_left_corner_val_before = self.star.image.array[0,0]
        self.star._check_image_inversion()
        top_left_corner_val_after = self.star.image.array[0,0]
        self.assertNotEqual(top_left_corner_val_before, top_left_corner_val_after)

    def test_results(self):
        """Test that the wobble diameter is similar to what it has been shown to be (0.495)."""
        self.star.analyze()
        self.assertAlmostEqual(self.star.wobble.radius_mm, self.demo_radius, delta=0.1)

        # Test that the center of the wobble circle is close to what it's shown to be (1511.5, 1302.2).
        # test y-coordinate
        y_coord = self.star.wobble.center.y
        self.assertAlmostEqual(y_coord, self.demo_y, delta=1)
        # test x-coordinate
        x_coord = self.star.wobble.center.x
        self.assertAlmostEqual(x_coord, self.demo_x, delta=1)

        """Test than the number of lines found is what is expected."""
        self.assertEqual(len(self.star.lines), self.demo_lines,
                         msg="The number of lines found was not the number expected")

    def test_pass_with_startpoint_off(self):
        """Test that the algo will still pass if the start point is set
            to a point somewhat different than actual center.
        """
        offset_point = Point(1400, 1300)
        self.star.set_start_point(offset_point, warn_if_far_away=False)
        self.star.analyze()

        self.assertTrue(self.star.passed)

        y_coord = self.star.wobble.center.y
        self.assertAlmostEqual(y_coord, self.demo_y, delta=3)
        x_coord = self.star.wobble.center.x
        self.assertAlmostEqual(x_coord, self.demo_x, delta=3)
        wobble_radius = self.star.wobble.radius_mm
        self.assertAlmostEqual(wobble_radius, self.demo_radius, delta=0.5)


class Test_Star_Test_Files(unittest.TestCase):
    """Tests specifically for the starshot test files."""
    # TODO: convert this to base/sublcasses like CBCT test suite

    _test_file_names = ('6XCollStar.tif', '10XCollStar.bmp', 'starshot_gantry.tif')
    file_paths = [osp.join(test_file_dir, test_file_name) for test_file_name in _test_file_names]

    # @classmethod
    # def setUpClass(cls):
    #     cls.stars = []
    #     for pth in cls.file_paths:
    #         star = Starshot()
    #         star.load_image(pth)
    #         cls.stars.append(star)

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
        wobble_radii = (2.3, 0.9, 5.2)
        for star, wobble_radius in zip(self.stars, wobble_radii):
            star.analyze()
            self.assertAlmostEqual(star.wobble.radius, wobble_radius, delta=0.1)

    def test_wobble_center(self):
        """Test that the center of the wobble circle is close to what it's shown to be (1511.5, 1302.2)."""
        wobble_centers = (Point(1297,1698), Point(1368.6,1453), Point(1302, 1512))
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
        num_expected_lines = (9, 4, 9)
        for star, num_exp_lines in zip(self.stars, num_expected_lines):
            star.analyze()
            self.assertEqual(len(star.lines), num_exp_lines,
                             msg="The number of radiation lines found was not the number expected")


