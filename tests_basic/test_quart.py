import io
import json
import os
from pathlib import Path
from unittest import TestCase

from matplotlib import pyplot as plt
from scipy import ndimage

from pylinac import QuartDVT
from pylinac.core.geometry import Point
from pylinac.core.io import TemporaryZipDirectory
from pylinac.quart import ACRYLIC, POLY, QuartDVTResult
from tests_basic.utils import (
    CloudFileMixin,
    FromZipTesterMixin,
    InitTesterMixin,
    get_file_from_cloud_test_repo,
    save_file,
)

TEST_DIR = ["CBCT", "Quart"]


class TestQuartDVT(TestCase, FromZipTesterMixin, InitTesterMixin):
    klass = QuartDVT
    zip = [*TEST_DIR, "Head_Quart.zip"]
    init_file = [*TEST_DIR, "Head"]
    is_folder = True

    def test_load_from_list_of_paths(self):
        # shouldn't raise
        with TemporaryZipDirectory(get_file_from_cloud_test_repo(self.zip)) as zfolder:
            paths = [Path(zfolder) / f for f in os.listdir(zfolder)]
            QuartDVT(paths)

    def test_load_from_list_of_streams(self):
        # shouldn't raise
        with TemporaryZipDirectory(get_file_from_cloud_test_repo(self.zip)) as zfolder:
            paths = [Path(zfolder, f) for f in os.listdir(zfolder)]
            paths = [io.BytesIO(open(p, "rb").read()) for p in paths]
            QuartDVT(paths)


class TestQuartDVTGeneral(TestCase):
    def setUp(self):
        self.path = get_file_from_cloud_test_repo([*TEST_DIR, "Head_Quart.zip"])
        self.quart = QuartDVT.from_zip(self.path)

    def test_phan_center(self):
        """Test locations of the phantom center."""
        known_phan_center = Point(256, 256)
        self.quart.analyze()
        self.assertAlmostEqual(
            self.quart.hu_module.phan_center.x, known_phan_center.x, delta=0.7
        )
        self.assertAlmostEqual(
            self.quart.hu_module.phan_center.y, known_phan_center.y, delta=0.7
        )

    def test_results_data(self):
        self.quart.analyze()
        data = self.quart.results_data()
        self.assertIsInstance(data, QuartDVTResult)
        self.assertEqual(data.num_images, self.quart.num_images)

        # check the additional modules got added
        self.assertIsInstance(data.hu_module.rois, dict)
        self.assertIsInstance(data.geometric_module.mean_high_contrast_distance, float)

        data_dict = self.quart.results_data(as_dict=True)
        self.assertIsInstance(data_dict, dict)

        data_str = self.quart.results_data(as_json=True)
        self.assertIsInstance(data_str, str)
        # shouldn't raise
        json.loads(data_str)

    def test_lazy_same_as_default(self):
        """Test that the results are the same from a lazy load vs default"""
        lazy_quart = QuartDVT.from_zip(self.path, memory_efficient_mode=True)
        self.quart.analyze()
        lazy_quart.analyze()
        self.assertEqual(self.quart.results(), lazy_quart.results())
        # results data should be the same except the time of evaluation (differs by ms)
        eager_results = self.quart.results_data(as_dict=True)
        eager_results.pop("date_of_analysis")
        lazy_results = lazy_quart.results_data(as_dict=True)
        lazy_results.pop("date_of_analysis")
        self.assertEqual(eager_results, lazy_results)


class TestPlottingSaving(TestCase):
    @classmethod
    def setUpClass(cls):
        path = get_file_from_cloud_test_repo([*TEST_DIR, "Head_Quart.zip"])
        cls.quart = QuartDVT.from_zip(path)
        cls.quart.analyze()

    @classmethod
    def tearDownClass(cls):
        plt.close("all")

    def test_plot_images(self):
        """Test that saving an image does something."""
        save_file(self.quart.plot_images)

    def test_pdf(self):
        # shouldn't raise
        self.quart.publish_pdf(io.BytesIO())

    def test_save_images(self):
        """Test that saving an image does something."""
        save_file(self.quart.save_images, to_single_file=False)

    def test_save_as_stream(self):
        stream_dict = self.quart.save_images(to_stream=True)
        self.assertIsInstance(stream_dict, dict)
        self.assertIsInstance(stream_dict["HU linearity"], io.BytesIO)

    def test_subimages_errors(self):
        """We don't use subimages here. easier to pass as a list of figs"""
        with self.assertRaises(NotImplementedError):
            self.quart.plot_analyzed_subimage("sr")
        with self.assertRaises(NotImplementedError):
            self.quart.save_analyzed_subimage("sr")

    def test_set_figure_size(self):
        self.quart.plot_analyzed_image(figsize=(8, 13))
        fig = plt.gcf()
        self.assertEqual(fig.bbox_inches.height, 13)
        self.assertEqual(fig.bbox_inches.width, 8)


class QuartDVTMixin(CloudFileMixin):
    dir_path = ["CBCT", "Quart"]
    origin_slice: int
    phantom_roll: float = 0
    snr: float
    cnr: float
    slice_thickness: float
    high_contrast_distance: float = 0.0
    horiz_dist = float
    vert_dist = float
    hu_values: dict
    unif_values = {
        "Center": ACRYLIC,
        "Left": ACRYLIC,
        "Right": ACRYLIC,
        "Top": ACRYLIC,
        "Bottom": ACRYLIC,
    }

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        cls.quart = QuartDVT.from_zip(filename, memory_efficient_mode=True)
        cls.quart.analyze()

    def test_roll(self):
        self.assertAlmostEqual(self.quart.catphan_roll, self.phantom_roll, delta=0.3)

    def test_snr(self):
        self.assertAlmostEqual(self.quart.hu_module.signal_to_noise, self.snr, delta=1)

    def test_cnr(self):
        self.assertAlmostEqual(
            self.quart.hu_module.contrast_to_noise, self.cnr, delta=0.2
        )

    def test_distances(self):
        self.assertAlmostEqual(
            self.horiz_dist,
            self.quart.geometry_module.distances()["horizontal mm"],
            delta=0.3,
        )
        self.assertAlmostEqual(
            self.vert_dist,
            self.quart.geometry_module.distances()["vertical mm"],
            delta=0.3,
        )

    def test_mean_high_contrast_distance(self):
        self.assertAlmostEqual(
            self.high_contrast_distance,
            self.quart.geometry_module.mean_high_contrast_resolution(),
            delta=0.05,
        )

    def test_HU_values(self):
        """Test HU values."""
        for key, roi in self.quart.hu_module.rois.items():
            exp_val = self.hu_values[key]
            meas_val = roi.pixel_value
            self.assertAlmostEqual(exp_val, meas_val, delta=5)

    def test_uniformity_values(self):
        """Test Uniformity HU values."""
        for key, exp_val in self.unif_values.items():
            meas_val = self.quart.uniformity_module.rois[key].pixel_value
            self.assertAlmostEqual(exp_val, meas_val, delta=5)


class TestQuartHead(QuartDVTMixin, TestCase):
    file_name = "Head_Quart.zip"
    phantom_roll = 0.2
    slice_thickness = 1.9
    high_contrast_distance = 0.9
    snr = 50
    cnr = 6.45
    horiz_dist = 159.3
    vert_dist = 159.6
    hu_values = {"Poly": POLY, "Acrylic": 126, "Air": -999, "Teflon": 981}
    unif_values = {"Center": 114, "Left": 114, "Right": 136, "Top": 125, "Bottom": 127}


class TestQuartHeadOffset(QuartDVTMixin, TestCase):
    """Shift the phantom over by several pixels to ensure no row/col algorithm issues

    Unfortunately, I can't move them that far because the FOV is very tight
    """

    file_name = "Head_Quart.zip"
    phantom_roll = 0.2
    slice_thickness = 1.9
    horiz_dist = 159.3
    vert_dist = 159.6
    high_contrast_distance = 0.9
    snr = 50
    cnr = 6.45
    hu_values = {"Poly": POLY, "Acrylic": 126, "Air": -999, "Teflon": 981}
    unif_values = {"Center": 114, "Left": 114, "Right": 136, "Top": 125, "Bottom": 127}

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        cls.quart = QuartDVT.from_zip(filename)
        for img in cls.quart.dicom_stack:
            img.roll(direction="x", amount=40)
        cls.quart.localize()
        cls.quart.analyze()


class TestQuartHeadRotated(QuartDVTMixin, TestCase):
    """Rotate the phantom over by several pixels to ensure no row/col algorithm issues

    Unfortunately, I can't move them that far because the FOV is very tight
    """

    file_name = "Head_Quart.zip"
    phantom_roll = -2.8
    slice_thickness = 1.9
    horiz_dist = 159.3
    vert_dist = 159.6
    high_contrast_distance = 0.80
    snr = 50
    cnr = 6.45
    hu_values = {"Poly": POLY, "Acrylic": 126, "Air": -999, "Teflon": 981}
    unif_values = {"Center": 114, "Left": 114, "Right": 136, "Top": 125, "Bottom": 127}

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        cls.quart = QuartDVT.from_zip(filename)
        for img in cls.quart.dicom_stack:
            img.array = ndimage.rotate(img.array, angle=3, mode="nearest")
        cls.quart.localize()
        cls.quart.analyze()


class TestQuartPelvis(QuartDVTMixin, TestCase):
    file_name = "Pelvis_Quart.zip"
    phantom_roll = 0.2
    slice_thickness = 1.9
    snr = 173.8
    cnr = 28.0
    horiz_dist = 159.3
    vert_dist = 159.6
    high_contrast_distance = 0.82
    hu_values = {"Poly": -29, "Acrylic": 140, "Air": -1000, "Teflon": 989}
    unif_values = {"Center": 120, "Left": 132, "Right": 142, "Top": 136, "Bottom": 137}


class TestHypersightQuart(QuartDVTMixin, TestCase):
    file_name = "Hypersight Quart (w water).zip"
    phantom_roll = 0.0
    slice_thickness = 1.9
    snr = 400
    cnr = 56.0
    horiz_dist = 159.7
    vert_dist = 159.6
    high_contrast_distance = 0.6
    hu_values = {"Poly": -47, "Acrylic": 106, "Air": -1000, "Teflon": 963, "Water": 0}
    unif_values = {"Center": 98, "Left": 96, "Right": 92, "Top": 92, "Bottom": 93}


class TestTableQuart(QuartDVTMixin, TestCase):
    """This dataset has a foam separator that still catches sometimes. See RAM-3226"""

    file_name = "Jan-Quart.zip"
    phantom_roll = 0.0
    slice_thickness = 1.87
    snr = 37
    cnr = 6.2
    horiz_dist = 159.3
    vert_dist = 159.7
    high_contrast_distance = 0.6
    hu_values = {"Poly": -63, "Acrylic": 111, "Air": -1000, "Teflon": 970, "Water": 0}
    unif_values = {"Center": 112, "Left": 104, "Right": 104, "Top": 97, "Bottom": 124}
