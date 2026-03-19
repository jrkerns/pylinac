import io
import os
from pathlib import Path
from unittest import TestCase

from matplotlib import pyplot as plt

from pylinac import GEHeliosCTDaily
from pylinac.core.io import TemporaryZipDirectory
from pylinac.helios import GEHeliosResult
from tests_basic.core.test_utilities import QuaacTestBase
from tests_basic.utils import (
    CloudFileMixin,
    FromZipTesterMixin,
    InitTesterMixin,
    PlotlyTestMixin,
    get_file_from_cloud_test_repo,
    get_folder_from_cloud_repo,
    save_file,
)

TEST_DIR = ["Helios"]

get_folder_from_cloud_repo(TEST_DIR)


class TestLoading(TestCase, FromZipTesterMixin, InitTesterMixin):
    klass = GEHeliosCTDaily
    zip = [*TEST_DIR, "GEHeliosCTDaily1.zip"]
    init_file = [*TEST_DIR, "GEHeliosCTDaily1"]
    is_folder = True

    def test_load_from_list_of_paths(self):
        with TemporaryZipDirectory(get_file_from_cloud_test_repo(self.zip)) as zfolder:
            paths = [Path(zfolder) / f for f in os.listdir(zfolder)]
            GEHeliosCTDaily(paths)

    def test_load_from_list_of_streams(self):
        with TemporaryZipDirectory(get_file_from_cloud_test_repo(self.zip)) as zfolder:
            paths = [Path(zfolder, f) for f in os.listdir(zfolder)]
            paths = [io.BytesIO(open(p, "rb").read()) for p in paths]
            GEHeliosCTDaily(paths)


class TestGeneral(TestCase):
    def setUp(self):
        self.path = get_file_from_cloud_test_repo([*TEST_DIR, "GEHeliosCTDaily1.zip"])
        self.ct = GEHeliosCTDaily.from_zip(self.path)

    def test_results_data(self):
        self.ct.analyze()
        data = self.ct.results_data()
        self.assertIsInstance(data, GEHeliosResult)
        self.assertEqual(data.num_images, self.ct.num_images)

    def test_results_warnings(self):
        self.ct.analyze()
        data = self.ct.results_data()
        self.assertEqual(len(data.warnings), 0)


class TestPlottingSaving(TestCase):
    @classmethod
    def setUpClass(cls):
        path = get_file_from_cloud_test_repo([*TEST_DIR, "GEHeliosCTDaily1.zip"])
        cls.ct = GEHeliosCTDaily.from_zip(path)
        cls.ct.analyze()

    @classmethod
    def tearDownClass(cls):
        plt.close("all")
        del cls.ct

    def test_plot_images(self):
        """Test that saving an image does something."""
        save_file(self.ct.plot_images)

    def test_save_images(self):
        """Test that saving an image does something."""
        save_file(self.ct.save_images, to_single_file=False)

    def test_subimages_errors(self):
        """We don't use subimages here. easier to pass as a list of figs"""
        with self.assertRaises(NotImplementedError):
            self.ct.plot_analyzed_subimage("sr")
        with self.assertRaises(NotImplementedError):
            self.ct.save_analyzed_subimage("sr")

    def test_set_figure_size(self):
        self.ct.plot_analyzed_image(figsize=(8, 13))
        fig = plt.gcf()
        self.assertEqual(fig.bbox_inches.height, 13)
        self.assertEqual(fig.bbox_inches.width, 8)


class TestQuaac(QuaacTestBase, CloudFileMixin, TestCase):
    dir_path = TEST_DIR
    file_name = "GEHeliosCTDaily1.zip"

    def quaac_instance(self):
        filename = self.get_filename()
        ct = GEHeliosCTDaily.from_zip(filename)
        ct.analyze()
        return ct


class HeliosMixin(CloudFileMixin):
    """Mixin that loads a dataset, analyzes it, and validates expected values."""

    dir_path = TEST_DIR
    origin_slice: int
    phantom_roll: float = 0
    contrast_difference: float
    noise_stdev: float
    uniformity_difference: float
    high_contrast_mtf_50: float
    low_contrast_mean: float
    x_adjustment: float = 0
    y_adjustment: float = 0
    angle_adjustment: float = 0
    roi_size_factor: float = 1
    scaling_factor: float = 1

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        cls.ct = GEHeliosCTDaily.from_zip(filename)
        cls.ct.analyze(
            x_adjustment=cls.x_adjustment,
            y_adjustment=cls.y_adjustment,
            angle_adjustment=cls.angle_adjustment,
            roi_size_factor=cls.roi_size_factor,
            scaling_factor=cls.scaling_factor,
        )

    @classmethod
    def tearDownClass(cls):
        plt.close("all")
        del cls.ct
        super().tearDownClass()

    def test_roll(self):
        self.assertAlmostEqual(self.ct.catphan_roll, self.phantom_roll, delta=0.3)

    def test_origin_slice(self):
        self.assertEqual(self.ct.origin_slice, self.origin_slice)

    def test_contrast_difference(self):
        self.assertAlmostEqual(
            self.ct.contrast_scale_module.contrast_difference,
            self.contrast_difference,
            delta=2,
        )

    def test_noise_stdev(self):
        self.assertAlmostEqual(
            self.ct.noise_uniformity_module.noise_center_std,
            self.noise_stdev,
            delta=0.3,
        )

    def test_uniformity_difference(self):
        self.assertAlmostEqual(
            self.ct.noise_uniformity_module.uniformity_difference,
            self.uniformity_difference,
            delta=3,
        )

    def test_high_contrast_mtf(self):
        self.assertAlmostEqual(
            self.ct.high_contrast_module.mtf.relative_resolution(50),
            self.high_contrast_mtf_50,
            delta=0.1,
        )

    def test_low_contrast_mean(self):
        self.assertAlmostEqual(
            self.ct.low_contrast_module.mean,
            self.low_contrast_mean,
            delta=0.1,
        )


class Helios_1(HeliosMixin, PlotlyTestMixin, TestCase):
    file_name = "GEHeliosCTDaily1.zip"
    origin_slice = 5
    contrast_difference = 124.66
    noise_stdev = 4.40
    uniformity_difference = -0.01
    high_contrast_mtf_50 = 0.54
    low_contrast_mean = 0.83
    num_figs = 6

    def setUp(self) -> None:
        self.instance = self.ct


class Helios_2(HeliosMixin, TestCase):
    file_name = "GEHeliosCTDaily2.zip"
    origin_slice = 5
    contrast_difference = 126.15
    noise_stdev = 4.06
    uniformity_difference = -0.18
    high_contrast_mtf_50 = 0.55
    low_contrast_mean = 0.40
