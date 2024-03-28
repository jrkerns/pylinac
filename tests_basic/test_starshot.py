import io
import json
import os
import os.path as osp
import tempfile
from unittest import TestCase

import matplotlib.pyplot as plt
import numpy as np

from pylinac import Starshot
from pylinac.core.geometry import Point
from pylinac.starshot import StarshotResults
from tests_basic.utils import (
    CloudFileMixin,
    FromURLTesterMixin,
    get_file_from_cloud_test_repo,
    get_folder_from_cloud_test_repo,
    save_file,
)

plt.close("all")
TEST_DIR = "Starshot"


class TestStarshotLoading(TestCase, FromURLTesterMixin):
    klass = Starshot
    url = "starshot.tif"
    kwargs = url_kwargs = {"dpi": 30, "sid": 1000}

    def test_load_from_file_object(self):
        with open(
            get_file_from_cloud_test_repo([TEST_DIR, "Starshot-30-deg-perfect.dcm"]),
            "rb",
        ) as f:
            star = Starshot(f)
            star.analyze()
        self.assertIsInstance(star, Starshot)

    def test_load_from_stream(self):
        with open(
            get_file_from_cloud_test_repo([TEST_DIR, "Starshot-30-deg-perfect.dcm"]),
            "rb",
        ) as f:
            s = io.BytesIO(f.read())
            star = Starshot(s)
            star.analyze()
        self.assertIsInstance(star, Starshot)

    def test_no_dpi(self):
        # raise error when DPI isn't in image or given
        with self.assertRaises(ValueError):
            Starshot.from_url(self.full_url)
        # but is fine when DPI is given
        Starshot.from_url(self.full_url, **self.kwargs)


class TestPlottingSaving(TestCase):
    @classmethod
    def tearDownClass(cls):
        plt.close("all")

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
        save_file(self.star.save_analyzed_subimage, as_file_object="b")


class StarMixin(CloudFileMixin):
    """Mixin for testing a starshot image."""

    # dir_location = TEST_DIR
    dir_path = ["Starshot"]
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
    wobble_tolerance = 0.2
    kwargs = {"sid": 1000}
    verbose = False

    @classmethod
    def setUpClass(cls):
        cls.star = cls.construct_star()
        cls.star.analyze(
            recursive=cls.recursive,
            min_peak_height=cls.min_peak_height,
            fwhm=cls.fwxm,
            radius=cls.radius,
        )

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
        self.star.analyze(
            recursive=self.recursive,
            min_peak_height=self.min_peak_height,
            fwhm=self.fwxm,
            radius=self.radius,
        )
        self.assertEqual(
            self.star.passed, self.passes, msg="Wobble was not within tolerance"
        )

    def test_wobble_diameter(self):
        """Test than the wobble radius is similar to what it has been shown to be)."""
        self.assertAlmostEqual(
            self.star.wobble.diameter_mm,
            self.wobble_diameter_mm,
            delta=self.wobble_tolerance,
        )

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
        self.assertEqual(
            len(self.star.lines),
            self.num_rad_lines,
            msg="The number of radiation lines found was not the number expected",
        )

    def test_all_radii_give_same_wobble(self):
        """Test that the wobble stays roughly the same for all radii."""
        if self.test_all_radii:
            star = self.construct_star()
            radii = []
            for radius in np.linspace(0.9, 0.25, 8):
                star.analyze(
                    radius=float(radius),
                    min_peak_height=self.min_peak_height,
                    recursive=self.recursive,
                    fwhm=self.fwxm,
                )
                self.assertAlmostEqual(
                    star.wobble.diameter_mm,
                    self.wobble_diameter_mm,
                    delta=self.wobble_tolerance,
                )
                radii.append(star.wobble.diameter_mm)
            if self.verbose:
                print(
                    f"Radii mean: {np.mean(radii):2.2f}, range: {np.max(radii) - np.min(radii):2.2f}"
                )


class Demo(StarMixin, TestCase):
    """Specific tests_basic for the demo image"""

    wobble_diameter_mm = 0.30
    wobble_center = Point(1270, 1437)
    num_rad_lines = 4
    wobble_tolerance = 0.15
    # independently verified: 0.24-0.26mm
    delete_file = False

    @classmethod
    def construct_star(cls):
        return Starshot.from_demo_image()


class Multiples(StarMixin, TestCase):
    """Test a starshot composed of multiple individual EPID images."""

    num_rad_lines = 9
    wobble_center = Point(254, 192)
    wobble_diameter_mm = 0.7
    wobble_tolerance = 0.2
    dir_path = ["Starshot", "set"]
    is_dir = True
    delete_file = False

    @classmethod
    def get_filename(cls):
        """Return the canonical path to the file."""
        return get_folder_from_cloud_test_repo([*cls.dir_path])

    def test_loading_from_zip(self):
        img_zip = get_file_from_cloud_test_repo(["Starshot", "multiples.zip"])
        star = Starshot.from_zip(img_zip)
        # shouldn't raise
        star.analyze()


class Starshot1(StarMixin, TestCase):
    file_name = "Starshot-1.tif"
    wobble_center = Point(508, 683)
    wobble_diameter_mm = 0.23
    num_rad_lines = 4
    # outside 0.20-0.27mm


class StarshotPerfect30Deg(StarMixin, TestCase):
    file_name = "Starshot-30-deg-perfect.dcm"
    wobble_center = Point(639.5, 639.5)
    wobble_diameter_mm = 0.0
    num_rad_lines = 6


class Starshot1FWHM(Starshot1):
    fwhm = False


class CRStarshot(StarMixin, TestCase):
    file_name = "CR-Starshot.dcm"
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
            self.star.publish_pdf(t, notes="stuff", metadata={"Unit": "TB1"})

    def test_results(self):
        data = self.star.results()
        self.assertIsInstance(data, str)

        data = self.star.results(as_list=True)
        self.assertIsInstance(data, list)
        self.assertIsInstance(data[0], str)

    def test_results_data(self):
        data = self.star.results_data()
        self.assertIsInstance(data, StarshotResults)
        self.assertEqual(data.circle_radius_mm, self.star.wobble.radius_mm)

        data_dict = self.star.results_data(as_dict=True)
        self.assertIsInstance(data_dict, dict)

        data_json = self.star.results_data(as_json=True)
        self.assertIsInstance(data_json, str)
        # shouldn't raise
        json.loads(data_json)

    def test_set_figure_size(self):
        self.star.plot_analyzed_image(figsize=(7, 11))
        fig = plt.gcf()
        self.assertEqual(fig.bbox_inches.height, 11)
        self.assertEqual(fig.bbox_inches.width, 7)

    def test_set_figure_size_subimage(self):
        self.star.plot_analyzed_subimage(figsize=(7, 11))
        fig = plt.gcf()
        self.assertEqual(fig.bbox_inches.height, 11)
        self.assertEqual(fig.bbox_inches.width, 7)


class Starshot2(StarMixin, TestCase):
    file_name = "Starshot#2.tif"
    wobble_center = Point(566, 590)
    wobble_diameter_mm = 0.2
    num_rad_lines = 4
    # outside: 0.18-0.19


class Starshot3(StarMixin, TestCase):
    file_name = "Starshot#3.tif"
    wobble_center = Point(466, 595)
    wobble_diameter_mm = 0.32
    num_rad_lines = 6
    # outside 0.33


class Starshot4(StarMixin, TestCase):
    file_name = "Starshot#4.tif"
    wobble_center = Point(446, 565)
    wobble_diameter_mm = 0.38
    num_rad_lines = 6
    # outside 0.39


class Starshot5(StarMixin, TestCase):
    file_name = "Starshot#5.tif"
    wobble_center = Point(557, 580)
    wobble_diameter_mm = 0.15
    num_rad_lines = 4
    wobble_tolerance = 0.2
    test_all_radii = False
    # outside: 0.14


class Starshot6(StarMixin, TestCase):
    file_name = "Starshot#6.tif"
    wobble_center = Point(528, 607)
    wobble_diameter_mm = 0.3
    num_rad_lines = 7


class Starshot7(StarMixin, TestCase):
    file_name = "Starshot#7.tif"
    wobble_center = Point(469, 646)
    wobble_diameter_mm = 0.2
    num_rad_lines = 4
    wobble_tolerance = 0.2


class Starshot8(StarMixin, TestCase):
    file_name = "Starshot#8.tiff"
    wobble_center = Point(686, 669)
    wobble_diameter_mm = 0.35
    num_rad_lines = 5


class Starshot9(StarMixin, TestCase):
    file_name = "Starshot#9.tiff"
    wobble_center = Point(714, 611)
    wobble_diameter_mm = 0.3
    num_rad_lines = 5


class Starshot10(StarMixin, TestCase):
    file_name = "Starshot#10.tiff"
    wobble_center = Point(725, 802)
    wobble_diameter_mm = 0.65
    num_rad_lines = 5


class Starshot11(StarMixin, TestCase):
    file_name = "Starshot#11.tiff"
    wobble_center = Point(760, 650)
    wobble_diameter_mm = 0.6
    num_rad_lines = 4


class Starshot12(StarMixin, TestCase):
    file_name = "Starshot#12.tiff"
    wobble_center = Point(315, 292)
    wobble_diameter_mm = 0.88
    num_rad_lines = 4


class Starshot13(StarMixin, TestCase):
    file_name = "Starshot#13.tiff"
    wobble_center = Point(376, 303)
    wobble_diameter_mm = 0.2
    num_rad_lines = 4


class Starshot14(StarMixin, TestCase):
    file_name = "Starshot#14.tiff"
    wobble_center = Point(334, 282)
    wobble_diameter_mm = 0.55
    num_rad_lines = 4


class Starshot15(StarMixin, TestCase):
    file_name = "Starshot#15.tiff"
    wobble_center = Point(346, 309)
    wobble_diameter_mm = 0.6
    num_rad_lines = 4


class Starshot16(StarMixin, TestCase):
    file_name = "Starshot#16.tiff"
    wobble_center = Point(1444, 1452)
    wobble_diameter_mm = 0.6
    num_rad_lines = 6


class Starshot17(StarMixin, TestCase):
    file_name = "Starshot#17.tiff"
    wobble_center = Point(1475, 1361)
    wobble_diameter_mm = 0.44
    num_rad_lines = 6


class Starshot18(StarMixin, TestCase):
    file_name = "Starshot#18.tiff"
    wobble_center = Point(1516, 1214)
    wobble_diameter_mm = 0.6
    num_rad_lines = 6


class Starshot19(StarMixin, TestCase):
    file_name = "Starshot#19.tiff"
    wobble_center = Point(1475, 1276)
    wobble_diameter_mm = 0.6
    num_rad_lines = 6


class Starshot20(StarMixin, TestCase):
    file_name = "Starshot#20.tiff"
    wobble_center = Point(347, 328)
    wobble_diameter_mm = 0.75
    num_rad_lines = 4


class Starshot21(StarMixin, TestCase):
    file_name = "Starshot#21.tiff"
    wobble_center = Point(354, 294)
    wobble_diameter_mm = 1.1
    wobble_tolerance = 0.2
    num_rad_lines = 4
    passes = False


class Starshot22(StarMixin, TestCase):
    file_name = "Starshot#22.tiff"
    wobble_center = Point(1305, 1513)
    wobble_diameter_mm = 0.9
    num_rad_lines = 9
    # outside 0.93mm


class Starshot23(StarMixin, TestCase):
    file_name = "Starshot#23.tiff"
    wobble_center = Point(1297, 1699)
    wobble_diameter_mm = 0.38
    num_rad_lines = 9


class Starshot24(StarMixin, TestCase):
    file_name = "Starshot#24.tiff"
    wobble_center = Point(1370, 1454)
    wobble_diameter_mm = 0.3
    num_rad_lines = 4


class Starshot25(StarMixin, TestCase):
    file_name = "Starshot#25.tiff"
    wobble_center = Point(286, 279)
    wobble_diameter_mm = 0.3
    num_rad_lines = 4


class Starshot26(StarMixin, TestCase):
    file_name = "Starshot#26.tiff"
    wobble_center = Point(1511, 1452)
    wobble_diameter_mm = 0.55
    num_rad_lines = 4
    wobble_tolerance = 0.15


class Starshot27(StarMixin, TestCase):
    file_name = "Starshot#27.tiff"
    wobble_center = Point(1105, 1306)
    wobble_diameter_mm = 0.4
    num_rad_lines = 6


class ChicagoSet(StarMixin, TestCase):
    file_name = "Chicago"
    wobble_center = Point(638, 639.3)
    wobble_diameter_mm = 0.7
    num_rad_lines = 5
    is_dir = True

    @classmethod
    def get_filename(cls):
        """Return the canonical path to the file."""
        return get_folder_from_cloud_test_repo([*cls.dir_path, cls.file_name])
