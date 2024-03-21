from __future__ import annotations

from typing import Iterable
from unittest import TestCase, skip

from pylinac import WinstonLutzMultiTargetMultiField
from pylinac.winston_lutz import BBArrangement, WinstonLutzMultiTargetSingleField
from tests_basic.utils import CloudFileMixin

TEST_DIR = "Winston-Lutz"


class TestWLMultiImage(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.wl = WinstonLutzMultiTargetMultiField.from_demo_images()
        cls.wl.analyze(BBArrangement.DEMO)

    def test_demo(self):
        # shouldn't raise
        WinstonLutzMultiTargetMultiField.run_demo()

    def test_publish_pdf(self):
        self.wl.publish_pdf("output.pdf")

    def test_save_images(self):
        self.wl.save_images()

    def test_save_images_to_stream(self):
        self.wl.save_images_to_stream()

    def test_no_axis_plot(self):
        with self.assertRaises(NotImplementedError):
            self.wl.plot_axis_images()

    def test_no_summary_plot(self):
        with self.assertRaises(NotImplementedError):
            self.wl.plot_summary()

    def test_no_plot_location(self):
        with self.assertRaises(NotImplementedError):
            self.wl.plot_location()

    def test_results(self):
        results = self.wl.results()
        self.assertIn("Multi-Target Multi-Field", results)
        self.assertIn("Max 2D distance of any BB->Field: 0.00 mm", results)

    def test_results_data(self):
        results = self.wl.results_data()
        self.assertEqual(results.max_2d_field_to_bb_mm, 0.0)
        self.assertEqual(results.bb_maxes["Iso"], 0.0)
        self.assertEqual(results.num_total_images, 4)

    def test_no_gantry_iso_size(self):
        with self.assertRaises(NotImplementedError):
            self.wl.gantry_iso_size

    def test_no_collimator_iso_size(self):
        with self.assertRaises(NotImplementedError):
            self.wl.collimator_iso_size

    def test_no_couch_iso_size(self):
        with self.assertRaises(NotImplementedError):
            self.wl.couch_iso_size

    def test_no_gantry_coll_iso_size(self):
        with self.assertRaises(NotImplementedError):
            self.wl.gantry_coll_iso_size


class WinstonLutzMultiTargetMultFieldMixin(CloudFileMixin):
    dir_path = ["Winston-Lutz"]
    num_images = 0
    zip = True
    bb_size = 5
    print_results = False
    arrangement: Iterable[dict]
    wl: WinstonLutzMultiTargetMultiField
    loader = WinstonLutzMultiTargetMultiField
    max_2d_distance: float
    mean_2d_distance: float
    median_2d_distance: float
    bb_maxes: dict[str, float] = {}

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        if cls.zip:
            cls.wl = cls.loader.from_zip(filename)
        else:
            cls.wl = cls.loader(filename)
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

    def test_bb_maxes(self):
        results = self.wl.results_data()
        for key, value in self.bb_maxes.items():
            self.assertAlmostEqual(results.bb_maxes[key], value, delta=0.05)


class SNCMultiMet(WinstonLutzMultiTargetMultFieldMixin, TestCase):
    dir_path = ["Winston-Lutz", "multi_target_multi_field"]
    file_name = "SNC_MM_KB.zip"
    num_images = 13
    arrangement = BBArrangement.SNC_MULTIMET
    max_2d_distance = 0.78
    median_2d_distance = 0.25
    mean_2d_distance = 0.27
    bb_maxes = {"Iso": 0.42, "1": 0.63}


class WinstonLutzMultiTargetSingleFieldMixin(WinstonLutzMultiTargetMultFieldMixin):
    loader = WinstonLutzMultiTargetSingleField
    arrangement = BBArrangement.ISOCAL
    is_open_field: bool = False
    wl: WinstonLutzMultiTargetSingleField

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        if cls.zip:
            cls.wl = cls.loader.from_zip(filename)
        else:
            cls.wl = cls.loader(filename)
        cls.wl.analyze(cls.arrangement, is_open_field=cls.is_open_field)
        if cls.print_results:
            print(cls.wl.results())


@skip("MPC/Single-Field not yet supported")
class MPCSubset(WinstonLutzMultiTargetSingleFieldMixin, TestCase):
    dir_path = ["MPC"]
    file_name = "6xsubset.zip"
    num_images = 3
    arrangement = BBArrangement.ISOCAL
    max_2d_distance = 0.78
    median_2d_distance = 0.56
    mean_2d_distance = 0.58
    is_open_field = True


# Test if no BBs found on an image
# Test if matches cannot be made
