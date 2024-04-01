from __future__ import annotations

import shutil
import tempfile
from typing import Sequence
from unittest import TestCase

import matplotlib.pyplot as plt

from pylinac import WinstonLutzMultiTargetMultiField
from pylinac.core.image_generator import (
    AS1200Image,
    GaussianFilterLayer,
    PerfectFieldLayer,
    generate_winstonlutz_multi_bb_multi_field,
)
from pylinac.winston_lutz import BBArrangement, BBConfig
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
        fig, ax = self.wl.plot_location()
        self.assertIsInstance(fig, plt.Figure)
        self.assertIsInstance(ax, plt.Axes)

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
    arrangement: tuple[BBConfig, ...]
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


class SyntheticMultiMetMixin(WinstonLutzMultiTargetMultFieldMixin):
    tmp_path: str = ""
    zip = False
    arrangement: BBArrangement
    field_offsets: Sequence[Sequence[float, float]]
    bb_offsets: Sequence[Sequence[float, float]]
    images_axes: Sequence[Sequence[int, int, int]] = (
        (0, 0, 0),
        (90, 0, 0),
        (180, 0, 0),
        (270, 0, 0),
        # (0, 0, 45),
        # (0, 0, 90),
        # (0, 0, 270),
        # (0, 0, 315),
    )
    field_size = (20.0, 20.0)
    bb_size = 5.0
    max_2d_distance = 0
    median_2d_distance = 0
    mean_2d_distance = 0

    @classmethod
    def tearDownClass(cls):
        # clean up the folder we created;
        # in BB space can be at a premium.
        shutil.rmtree(cls.tmp_path, ignore_errors=True)

    @classmethod
    def get_filename(cls) -> str:
        """We generate the files and return a local temp path.
        This may get called multiple times so we do a poor-man's caching"""
        if not cls.tmp_path:
            cls.tmp_path = tempfile.mkdtemp()
            generate_winstonlutz_multi_bb_multi_field(
                simulator=AS1200Image(1000),
                field_layer=PerfectFieldLayer,
                dir_out=cls.tmp_path,
                field_offsets=cls.field_offsets,
                bb_offsets=cls.bb_offsets,
                field_size_mm=cls.field_size,
                bb_size_mm=cls.bb_size,
                final_layers=[GaussianFilterLayer(sigma_mm=1)],
                image_axes=cls.images_axes,
            )
        return cls.tmp_path

    @property
    def num_images(self):
        """Shortcut the num of images check since we are creating them. No need to check."""
        return len(self.images_axes)


class SyntheticPerfect1BB(SyntheticMultiMetMixin, TestCase):
    arrangement = (
        BBConfig(
            name="Iso",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=0,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
    )
    bb_size = 5
    field_size = (20, 20)
    field_offsets = [(0, 0, 0)]
    bb_offsets = [(0, 0, 0)]
    max_2d_distance = 0
    median_2d_distance = 0
    mean_2d_distance = 0


class Synthetic1BBOffsetIn(SyntheticMultiMetMixin, TestCase):
    arrangement = (
        BBConfig(
            name="Iso",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=0,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
    )
    bb_size = 5
    field_size = (20, 20)
    field_offsets = [(0, 0, 0)]
    bb_offsets = [(0, 0, 1)]  # 1mm offset
    max_2d_distance = 1
    median_2d_distance = 1
    mean_2d_distance = 1


class Synthetic1BBOffsetLeft(SyntheticMultiMetMixin, TestCase):
    arrangement = (
        BBConfig(
            name="Iso",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=0,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
    )
    bb_size = 5
    field_size = (20, 20)
    field_offsets = [(0, 0, 0)]
    bb_offsets = [(1, 0, 0)]  # 1mm offset
    max_2d_distance = 1
    median_2d_distance = 0.5
    mean_2d_distance = 0.5


class Synthetic1BBOffsetUp(SyntheticMultiMetMixin, TestCase):
    arrangement = (
        BBConfig(
            name="Iso",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=0,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
    )
    bb_size = 5
    field_size = (20, 20)
    field_offsets = [(0, 0, 0)]
    bb_offsets = [(0, 1, 0)]  # 1mm offset
    max_2d_distance = 1
    median_2d_distance = 0.5
    mean_2d_distance = 0.5


class Synthetic2BBPerfect(SyntheticMultiMetMixin, TestCase):
    arrangement = (
        BBConfig(
            name="Iso",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=0,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
        BBConfig(
            name="1",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=-30,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
    )
    bb_size = 5
    field_size = (20, 20)
    field_offsets = [(0, 0, 0), (0, 0, -30)]
    bb_offsets = [(0, 0, 0), (0, 0, -30)]
    max_2d_distance = 0
    median_2d_distance = 0
    mean_2d_distance = 0


class Synthetic2BB1OffLeft(SyntheticMultiMetMixin, TestCase):
    arrangement = (
        BBConfig(
            name="Iso",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=0,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
        BBConfig(
            name="Left",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=-30,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
    )
    bb_size = 5
    field_size = (20, 20)
    field_offsets = [(0, 0, 0), (0, 0, -30)]
    bb_offsets = [(0, 0, 0), (1, 0, -30)]
    max_2d_distance = 1
    median_2d_distance = 0
    mean_2d_distance = 0.25


# class WinstonLutzMultiTargetSingleFieldMixin(WinstonLutzMultiTargetMultFieldMixin):
#     loader = WinstonLutzMultiTargetSingleField
#     arrangement = BBArrangement.ISOCAL
#     is_open_field: bool = False
#     wl: WinstonLutzMultiTargetSingleField
#
#     @classmethod
#     def setUpClass(cls):
#         filename = cls.get_filename()
#         if cls.zip:
#             cls.wl = cls.loader.from_zip(filename)
#         else:
#             cls.wl = cls.loader(filename)
#         cls.wl.analyze(cls.arrangement, is_open_field=cls.is_open_field)
#         if cls.print_results:
#             print(cls.wl.results())


# @skip("MPC/Single-Field not yet supported")
# class MPCSubset(WinstonLutzMultiTargetSingleFieldMixin, TestCase):
#     dir_path = ["MPC"]
#     file_name = "6xsubset.zip"
#     num_images = 3
#     arrangement = BBArrangement.ISOCAL
#     max_2d_distance = 0.78
#     median_2d_distance = 0.56
#     mean_2d_distance = 0.58
#     is_open_field = True


# Test if no BBs found on an image
# Test if matches cannot be made
