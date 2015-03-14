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

    def test_analyze_without_images(self):
        star = Starshot()
        self.assertRaises(AttributeError, star.analyze)


class Star_Test:
    star_file = ''
    wobble_diameter_mm = 0
    wobble_center = Point()
    num_rad_lines = 0

    def test_passed(self):
        """Test that the demo image passed"""
        self.assertTrue(self.star.passed, msg="Wobble was not within tolerance")

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


class Test_Star_Demo(unittest.TestCase, Star_Test):
    """Specific tests for the demo image"""
    wobble_diameter_mm = 0.35
    wobble_center = Point(1270, 1438)
    num_rad_lines = 4

    def setUp(self):
        self.star = Starshot()
        self.star.load_demo_image()

    def test_passed(self):
        self.star.analyze()
        super().test_passed()

    def test_failed_with_tight_tol(self):
        self.star._tolerance = 0.1
        self.star.analyze()
        self.assertFalse(self.star.passed)

    def test_wobble_center(self):
        self.star.analyze()
        super().test_wobble_center()

    def test_num_rad_lines(self):
        self.star.analyze()
        super().test_num_rad_lines()

    def test_wobble_diameter(self):
        self.star.analyze()
        super().test_wobble_diameter()

    def test_bad_inputs(self):
        self.star.analyze(radius=0.3, min_peak_height=0.1)
        self.test_wobble_center()
        self.test_wobble_diameter()

    def test_demo_runs(self):
        """Test that the demo runs without error."""
        # TODO: come up with decorator that adds show parameter
        self.star.run_demo()

    def test_image_inverted(self):
        """Check that the demo image was actually inverted, as it needs to be."""
        top_left_corner_val_before = self.star.image.array[0,0]
        self.star._check_image_inversion()
        top_left_corner_val_after = self.star.image.array[0,0]
        self.assertNotEqual(top_left_corner_val_before, top_left_corner_val_after)

    def test_SID_can_be_overridden_for_nonEPID(self):
        self.star.analyze(SID=400)
        self.assertNotEqual(self.star.wobble.diameter, self.wobble_diameter_mm*2)


class Test_Coll(unittest.TestCase, Star_Test):
    star_file = osp.join(test_file_dir, '6XCollStar.tif')
    wobble_center = Point(1297, 1699)
    wobble_diameter_mm = 0.37
    num_rad_lines = 9

    @classmethod
    def setUpClass(cls):
        cls.star = Starshot()
        cls.star.load_image(cls.star_file)
        cls.star.analyze(recursive=True)

    def test_not_fwhm_passes(self):
        self.star.analyze(fwhm=False)
        self.test_passed()


class Test_Coll2(unittest.TestCase, Star_Test):
    star_file = osp.join(test_file_dir, '10XCollStar.bmp')
    wobble_center = Point(1370, 1454)
    wobble_diameter_mm = 0.41
    num_rad_lines = 4

    @classmethod
    def setUpClass(cls):
        cls.star = Starshot()
        cls.star.load_image(cls.star_file)
        cls.star.analyze(recursive=True)


class Test_Gantry(unittest.TestCase, Star_Test):
    star_file = osp.join(test_file_dir, 'starshot_gantry.tif')
    wobble_center = Point(1302, 1513)
    wobble_diameter_mm = 1.15
    num_rad_lines = 9

    @classmethod
    def setUpClass(cls):
        cls.star = Starshot()
        cls.star.load_image(cls.star_file)
        cls.star.analyze(recursive=True)

    # def test_passed(self):
    #     # this image does *not* pass
    #     self.assertFalse(self.star.passed)

    # def test_print_fail(self):
    #     self.star.analyze(0.3, recursive=False)
    #     self.assertTrue('FAIL' in self.star.return_results())

    def test_bad_input(self):
        self.star.analyze(radius=0.3, min_peak_height=0.1)
        self.test_wobble_diameter()