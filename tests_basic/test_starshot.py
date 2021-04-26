import os
import os.path as osp
import tempfile
from unittest import TestCase

import matplotlib.pyplot as plt
import numpy as np

from pylinac.core.geometry import Point
from pylinac import Starshot
from tests_basic.utils import save_file, LoadingTestBase, LocationMixin

plt.close('all')
TEST_DIR = osp.join(osp.dirname(__file__), 'test_files', 'Starshot')


class TestStarshotLoading(LoadingTestBase, TestCase):
    klass = Starshot
    url = 'starshot.tif'
    kwargs = {'dpi': 30, 'sid': 1000}

    def test_no_dpi(self):
        # raise error when DPI isn't in image or given
        with self.assertRaises(ValueError):
            Starshot.from_url(self.full_url)
        # but is fine when DPI is given
        Starshot.from_url(self.full_url, **self.kwargs)


class TestPlottingSaving(TestCase):

    def setUp(self):
        self.star = Starshot.from_demo_image()
        self.star.analyze()

    def test_save_analyzed_image(self):
        """Test that saving an image does something."""
        save_file(self.star.save_analyzed_image)

    def test_save_analyzed_subimage(self):
        # save as normal file
        save_file(self.star.save_analyzed_subimage)
        # save into buffer
        save_file(self.star.save_analyzed_subimage, as_file_object='b')


class StarMixin(LocationMixin):
    """Mixin for testing a starshot image."""
    dir_location = TEST_DIR
    is_dir = False  # whether the starshot is a single file (False) or directory of images to combine (True)
    wobble_diameter_mm = 0
    wobble_center = Point()
    num_rad_lines = 0
    recursive = True
    passes = True
    min_peak_height = 0.25
    radius = 0.85
    test_all_radii = True
    fwxm = True
    wobble_tolerance = 0.1
    kwargs = {'sid': 1000}
    verbose = False

    @classmethod
    def setUpClass(cls):
        cls.star = cls.construct_star()
        cls.star.analyze(recursive=cls.recursive, min_peak_height=cls.min_peak_height, fwhm=cls.fwxm, radius=cls.radius)

    @classmethod
    def construct_star(cls):
        filename = cls.get_filename()
        if cls.is_dir:
            files = [osp.join(filename, file) for file in os.listdir(filename)]
            star = Starshot.from_multiple_images(files, **cls.kwargs)
        else:
            star = Starshot(filename, **cls.kwargs)
        return star

    def test_passed(self):
        """Test that the demo image passed"""
        self.star.analyze(recursive=self.recursive, min_peak_height=self.min_peak_height, fwhm=self.fwxm, radius=self.radius)
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
            star = self.construct_star()
            radii = []
            for radius in np.linspace(0.9, 0.25, 8):
                star.analyze(radius=float(radius), min_peak_height=self.min_peak_height, recursive=self.recursive, fwhm=self.fwxm)
                self.assertAlmostEqual(star.wobble.diameter_mm, self.wobble_diameter_mm, delta=self.wobble_tolerance)
                radii.append(star.wobble.diameter_mm)
            if self.verbose:
                print(f"Radii mean: {np.mean(radii):2.2f}, range: {np.max(radii) - np.min(radii):2.2f}")


class Demo(StarMixin, TestCase):
    """Specific tests_basic for the demo image"""
    wobble_diameter_mm = 0.30
    wobble_center = Point(1270, 1437)
    num_rad_lines = 4
    wobble_tolerance = 0.15
    # outside program: 0.24-0.26mm

    @classmethod
    def construct_star(cls):
        return Starshot.from_demo_image()


class Multiples(StarMixin, TestCase):
    """Test a starshot composed of multiple individual EPID images."""
    num_rad_lines = 9
    wobble_center = Point(254, 192)
    wobble_diameter_mm = 0.7
    wobble_tolerance = 0.2
    file_path = ['set']
    is_dir = True

    def test_loading_from_zip(self):
        img_zip = osp.join(TEST_DIR, 'set.zip')
        star = Starshot.from_zip(img_zip)
        # shouldn't raise
        star.analyze()


class Starshot1(StarMixin, TestCase):
    file_path = ['Starshot#1.tif']
    wobble_center = Point(508, 683)
    wobble_diameter_mm = 0.23
    num_rad_lines = 4
    # outside 0.20-0.27mm


class StarshotPerfect30Deg(StarMixin, TestCase):
    file_path = ['Starshot 30 deg perfect.dcm']
    wobble_center = Point(639.5, 639.5)
    wobble_diameter_mm = 0.0
    num_rad_lines = 6


class Starshot1FWHM(Starshot1):
    fwhm = False


class CRStarshot(StarMixin, TestCase):
    file_path = ['CR-Starshot.dcm']
    wobble_center = Point(1030.5, 1253.6)
    wobble_diameter_mm = 0.3
    num_rad_lines = 6
    # outside 0.25-0.26mm


class GeneralTests(Demo, TestCase):

    def test_demo_runs(self):
        """Test that the demo runs without error."""
        self.star.run_demo()

    def test_fails_with_tight_tol(self):
        star = Starshot.from_demo_image()
        star.analyze(tolerance=0.1)
        self.assertFalse(star.passed)

    def test_bad_inputs_still_recovers(self):
        self.star.analyze(radius=0.3, min_peak_height=0.1)
        self.test_wobble_center()
        self.test_wobble_diameter()

    def test_image_inverted(self):
        """Check that the demo image was actually inverted, as it needs to be."""
        star = Starshot.from_demo_image()
        top_left_corner_val_before = star.image.array[0, 0]
        star.image.check_inversion_by_histogram(percentiles=[4, 50, 96])
        top_left_corner_val_after = star.image.array[0, 0]
        self.assertNotEqual(top_left_corner_val_before, top_left_corner_val_after)

    def test_bad_start_point_recovers(self):
        """Test that even at a distance start point, the search algorithm recovers."""
        self.star.analyze(start_point=(1000, 1000))
        self.test_passed()
        self.test_wobble_center()
        self.test_wobble_diameter()

    def test_publish_pdf(self):
        with tempfile.TemporaryFile() as t:
            self.star.publish_pdf(t, notes='stuff', metadata={"Unit": 'TB1'})

