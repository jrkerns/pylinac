from __future__ import annotations

import shutil
import tempfile
from typing import Sequence
from unittest import TestCase

import matplotlib.pyplot as plt

from pylinac import WinstonLutzMultiTargetMultiField
from pylinac.core.geometry import Vector, cos, sin
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
    bb_shift_vector = Vector()
    bb_roll: float = 0
    bb_pitch: float = 0
    bb_yaw: float = 0
    bb_max_2d_yaw: float = 0

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

    def test_bb_cartesian_shift_vector(self):
        translation, y, p, r = self.wl.bb_shift_vector
        self.assertAlmostEqual(translation.x, self.bb_shift_vector.x, delta=0.05)
        self.assertAlmostEqual(translation.y, self.bb_shift_vector.y, delta=0.05)
        self.assertAlmostEqual(translation.z, self.bb_shift_vector.z, delta=0.05)

    def test_bb_shift_instructions(self):
        inst = self.wl.bb_shift_instructions()
        self.assertIsInstance(inst, str)
        if self.bb_shift_vector.x > 0:
            self.assertIn("RIGHT", inst)
        elif self.bb_shift_vector.x < 0:
            self.assertIn("LEFT", inst)
        if self.bb_shift_vector.y > 0:
            self.assertIn("IN", inst)
        elif self.bb_shift_vector.y < 0:
            self.assertIn("OUT", inst)
        if self.bb_shift_vector.z > 0:
            self.assertIn("UP", inst)
        elif self.bb_shift_vector.z < 0:
            self.assertIn("DOWN", inst)

    def test_roll(self):
        _, _, _, roll = self.wl.bb_shift_vector
        self.assertAlmostEqual(roll, self.bb_roll, delta=0.1)

    def test_yaw(self):
        _, yaw, _, _ = self.wl.bb_shift_vector
        self.assertAlmostEqual(yaw, self.bb_yaw, delta=0.1)

    def test_pitch(self):
        _, _, pitch, _ = self.wl.bb_shift_vector
        self.assertAlmostEqual(pitch, self.bb_pitch, delta=0.1)

    def test_2d_couch_yaw_error(self):
        d = self.wl._couch_rotation_error()
        self.assertAlmostEqual(
            max(v["yaw error"] for v in d.values()), self.bb_max_2d_yaw, delta=0.1
        )


class SNCMultiMet(WinstonLutzMultiTargetMultFieldMixin, TestCase):
    dir_path = ["Winston-Lutz", "multi_target_multi_field"]
    file_name = "SNC_MM_KB.zip"
    num_images = 19
    arrangement = BBArrangement.SNC_MULTIMET
    max_2d_distance = 0.78
    median_2d_distance = 0.25
    mean_2d_distance = 0.27
    bb_maxes = {"Iso": 0.58, "1": 0.68}
    bb_shift_vector = Vector(0.2, 0.0, 0.05)
    bb_pitch = 0.1
    bb_yaw = 0.1


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
        (0, 0, 45),
        (0, 0, 90),
        (0, 0, 270),
        (0, 0, 315),
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
    bb_offsets = [(0, 0, 1)]  # 1mm in
    max_2d_distance = 1
    median_2d_distance = 1
    mean_2d_distance = 1
    bb_shift_vector = Vector(0, -1, 0)  # vector is opposite of offset


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
    median_2d_distance = 1
    mean_2d_distance = 0.75
    bb_shift_vector = Vector(1, 0, 0)


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
    median_2d_distance = 0.0
    mean_2d_distance = 0.25
    bb_shift_vector = Vector(0, 0, -1)


class Synthetic3BBPerfect(SyntheticMultiMetMixin, TestCase):
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
            name="Out",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=-30,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
        BBConfig(
            name="Up/In",
            offset_left_mm=0,
            offset_up_mm=40,
            offset_in_mm=30,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
    )
    bb_size = 5
    field_size = (20, 20)
    field_offsets = [(0, 0, 0), (0, 0, -30), (0, 40, 30)]
    bb_offsets = [(0, 0, 0), (0, 0, -30), (0, 40, 30)]
    max_2d_distance = 0
    median_2d_distance = 0
    mean_2d_distance = 0


class Synthetic2BBYaw(SyntheticMultiMetMixin, TestCase):
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
        BBConfig(
            name="Up",
            offset_left_mm=0,
            offset_up_mm=40,
            offset_in_mm=0,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
    )
    bb_size = 5
    field_size = (20, 20)
    field_offsets = [(0, 0, 0), (0, 0, -30), (0, 40, 0)]  # left, up in
    bb_offsets = [(0, 0, 0), (1, 0, -30), (0, 40, 0)]
    max_2d_distance = 1
    median_2d_distance = 0
    mean_2d_distance = 0.25
    bb_yaw = 1.9  # one bb is off in the x, creating a yaw rotation
    bb_max_2d_yaw = 91.8  # not realistic due to choice of BB placement symmetry (roll is 180 degrees in one case), but it's the max; note the 91.8-90=1.8


class Synthetic2BBRoll(SyntheticMultiMetMixin, TestCase):
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
            name="Up",
            offset_left_mm=0,
            offset_up_mm=30,
            offset_in_mm=-30,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
        BBConfig(
            name="In",
            offset_left_mm=0,
            offset_up_mm=0,
            offset_in_mm=40,
            bb_size_mm=5,
            rad_size_mm=20,
        ),
    )
    bb_size = 5
    field_size = (20, 20)
    field_offsets = [(0, 0, 0), (0, 0, 40), (0, 30, -30)]  # left, up in
    bb_offsets = [(0, 0, 0), (0, 0, 40), (30 * sin(5), 30 * cos(5), -30)]
    max_2d_distance = 2.85
    median_2d_distance = 0
    mean_2d_distance = 0.7
    bb_roll = 5.2
    bb_max_2d_yaw = -2.1


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
