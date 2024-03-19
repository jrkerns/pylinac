from __future__ import annotations

from typing import Iterable
from unittest import TestCase

from pylinac import WinstonLutzMultiTargetMultiField
from pylinac.winston_lutz import BBArrangement
from tests_basic.utils import CloudFileMixin

TEST_DIR = "Winston-Lutz"


class TestWLMultiImage(TestCase):
    def test_demo_images(self):
        wl = WinstonLutzMultiTargetMultiField.from_demo_images()
        # shouldn't raise
        wl.analyze(BBArrangement.DEMO)

    def test_demo(self):
        # shouldn't raise
        WinstonLutzMultiTargetMultiField.run_demo()

    def test_publish_pdf(self):
        wl = WinstonLutzMultiTargetMultiField.from_demo_images()
        wl.analyze(BBArrangement.DEMO)
        wl.publish_pdf("output.pdf")

    def test_save_images(self):
        wl = WinstonLutzMultiTargetMultiField.from_demo_images()
        wl.analyze(BBArrangement.DEMO)
        wl.save_images()

    def test_save_images_to_stream(self):
        wl = WinstonLutzMultiTargetMultiField.from_demo_images()
        wl.analyze(BBArrangement.DEMO)
        wl.save_images_to_stream()

    def test_no_axis_plot(self):
        wl = WinstonLutzMultiTargetMultiField.from_demo_images()
        wl.analyze(BBArrangement.DEMO)
        with self.assertRaises(NotImplementedError):
            wl.plot_axis_images()

    def test_no_summary_plot(self):
        wl = WinstonLutzMultiTargetMultiField.from_demo_images()
        wl.analyze(BBArrangement.DEMO)
        with self.assertRaises(NotImplementedError):
            wl.plot_summary()


class WinstonLutzMultiTargetMultFieldMixin(CloudFileMixin):
    dir_path = ["Winston-Lutz"]
    num_images = 0
    zip = True
    bb_size = 5
    print_results = False
    arrangement: Iterable[dict]
    wl: WinstonLutzMultiTargetMultiField
    max_2d_distance: float
    mean_2d_distance: float
    median_2d_distance: float

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        if cls.zip:
            cls.wl = WinstonLutzMultiTargetMultiField.from_zip(filename)
        else:
            cls.wl = WinstonLutzMultiTargetMultiField(filename)
        cls.wl.analyze(cls.arrangement)
        if cls.print_results:
            print(cls.wl.results())

    def test_number_of_images(self):
        self.assertEqual(self.num_images, len(self.wl.images))

    def test_bb_max_distance(self):
        self.assertAlmostEqual(
            self.wl.max_bb_deviation_2d, self.max_2d_distance, delta=0.15
        )

    def test_bb_median_distance(self):
        self.assertAlmostEqual(
            self.wl.median_bb_deviation_2d,
            self.median_2d_distance,
            delta=0.1,
        )

    def test_bb_mean_distance(self):
        self.assertAlmostEqual(
            self.wl.mean_bb_deviation_2d, self.mean_2d_distance, delta=0.1
        )


class SNCMultiMet(WinstonLutzMultiTargetMultFieldMixin, TestCase):
    dir_path = ["Winston-Lutz", "multi_target_multi_field"]
    file_name = "SNC_MM_KB.zip"
    num_images = 13
    arrangement = BBArrangement.SNC_MULTIMET
    max_2d_distance = 0.78
    median_2d_distance = 0.56
    mean_2d_distance = 0.58


class MPCSubset(WinstonLutzMultiTargetMultFieldMixin, TestCase):
    dir_path = ["MPC"]
    file_name = "6xsubset.zip"
    num_images = 3
    arrangement = BBArrangement.ISOCAL
    max_2d_distance = 0.78
    median_2d_distance = 0.56
    mean_2d_distance = 0.58
