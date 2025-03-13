import copy
import json
import pytest
import math
from unittest import TestCase

from matplotlib import pyplot as plt
from skimage.transform import rotate

from pylinac.cheese import CIRS062M, TomoCheese, TomoCheeseResult
from tests_basic.core.test_utilities import QuaacTestBase, ResultsDataBase
from tests_basic.utils import (
    CloudFileMixin,
    FromDemoImageTesterMixin,
    FromURLTesterMixin,
    FromZipTesterMixin,
    InitTesterMixin,
    PlotlyTestMixin,
    save_file,
)

TEST_DIR = "Tomo"


@pytest.mark.proprietary
class TestInstantiation(
    TestCase,
    InitTesterMixin,
    FromDemoImageTesterMixin,
    FromURLTesterMixin,
    FromZipTesterMixin,
):
    klass = TomoCheese
    init_file = [TEST_DIR, "TomoTherapy Cheese Phantom"]
    demo_load_method = "from_demo_images"
    url = "TomoCheese.zip"
    zip = [TEST_DIR, "TomoCheese.zip"]
    is_folder = True


class TestResults(TestCase, ResultsDataBase):
    model = TomoCheese

    @classmethod
    def setUpClass(cls) -> None:
        cls.cheese = TomoCheese.from_demo_images()
        cls.cheese.analyze()

    @classmethod
    def tearDownClass(cls) -> None:
        plt.close("all")
        del cls.cheese

    def test_results_as_str(self):
        assert isinstance(self.cheese.results(), str)

    def test_results_as_list(self):
        r = self.cheese.results(as_list=True)
        assert isinstance(r, list)

    def test_results_data(self):
        r = self.cheese.results_data()
        assert isinstance(r, TomoCheeseResult)
        assert self.cheese.module.rois["5"].pixel_value == r.roi_5["median"]
        r = self.cheese.results_data(as_dict=True)
        assert isinstance(r, dict)
        assert self.cheese.module.rois["9"].std == r["roi_9"]["std"]

        data_json = self.cheese.results_data(as_json=True)
        assert isinstance(data_json, str)
        # shouldn't raise
        json.loads(data_json)

    def test_results_warnings(self):
        self.cheese.analyze()
        data = self.cheese.results_data()
        self.assertEqual(len(data.warnings), 0)


class TestGeneral(TestCase):
    def test_demo(self):
        TomoCheese.run_demo()

    def test_passing_origin_slice(self):
        cheese = TomoCheese.from_demo_images()
        cheese.analyze(origin_slice=26)  # automatic is 24
        assert cheese.origin_slice == 26


class TestAnalysis(TestCase):
    def test_cropping_before_analysis(self):
        cheese = TomoCheese.from_demo_images()
        for img in cheese.dicom_stack:
            img.crop(pixels=20, edges=("bottom",))
        # shouldn't raise
        cheese.analyze()

    def test_rolling_before_analysis(self):
        """Rolling (shifting) the phantom by a nominal amount shouldn't affect analysis"""
        cheese = TomoCheese.from_demo_images()
        cheese.analyze()
        original_roi_1 = copy.copy(cheese.module.rois["1"].pixel_value)
        for img in cheese.dicom_stack:
            img.roll(direction="x", amount=20)
        cheese.analyze()
        new_roi_1 = cheese.module.rois["1"].pixel_value
        assert math.isclose(original_roi_1, new_roi_1, abs_tol=3)

    def test_rotating_phantom(self):
        """Check that a roll is corrected"""
        cheese = TomoCheese.from_demo_images()
        cheese.analyze()
        assert math.isclose(cheese.catphan_roll, -0.25, abs_tol=0.05)
        for img in cheese.dicom_stack:
            img.array = rotate(img.array, angle=3, mode="edge")
        cheese.analyze()
        assert math.isclose(cheese.catphan_roll, -3.25, abs_tol=0.05)

    def test_too_much_rotation_resets(self):
        # too much roll will reset to 0 however
        cheese = TomoCheese.from_demo_images()
        for img in cheese.dicom_stack:
            img.array = rotate(img.array, angle=13, mode="edge")
        cheese.analyze()
        assert cheese.catphan_roll == 0

    def test_roi_config(self):
        cheese = TomoCheese.from_demo_images()
        config = {"3": {"density": 4.12}}
        cheese.analyze(roi_config=config)
        self.assertEqual(cheese.roi_config, config)


class TestPlottingSaving(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cheese = TomoCheese.from_demo_images()
        cls.cheese.analyze()

    @classmethod
    def tearDownClass(cls):
        plt.close("all")
        del cls.cheese

    def test_save_pdf(self):
        # shouldn't raise
        save_file(self.cheese.publish_pdf, "temp")

    def test_set_figure_size(self):
        self.cheese.plot_analyzed_image(figsize=(8, 13))
        fig = plt.gcf()
        self.assertEqual(fig.bbox_inches.height, 13)
        self.assertEqual(fig.bbox_inches.width, 8)

    def test_save_image(self):
        save_file(self.cheese.save_analyzed_image)

    def test_save_subimage_fails(self):
        """There is no sub-images for the tomo"""
        with self.assertRaises(NotImplementedError):
            self.cheese.save_analyzed_subimage()

    def test_plotting_without_density_fails(self):
        cheese = TomoCheese.from_demo_images()
        cheese.analyze()  # no roi config
        with self.assertRaises(ValueError):
            cheese.plot_density_curve()

    def test_plotting_density(self):
        cheese = TomoCheese.from_demo_images()
        cheese.analyze(roi_config={"1": {"density": 1.2}})
        cheese.plot_density_curve()


class TestTomoQuaac(QuaacTestBase, TestCase):
    def quaac_instance(self):
        t = TomoCheese.from_demo_images()
        t.analyze()
        return t


@pytest.mark.proprietary
class TestCIRS062Quaac(QuaacTestBase, CloudFileMixin, TestCase):
    dir_path = ["Tomo", "CIRS062M"]
    file_name = "CIRS062M - Erogluer.zip"

    def quaac_instance(self):
        t = CIRS062M.from_zip(self.get_filename())
        t.analyze()
        return t


class CheeseMixin(CloudFileMixin):
    model = TomoCheese
    origin_slice = 0
    dir_path = ["Tomo"]
    hu_values = {}
    expected_roll = -0.2
    x_adjustment = 0
    y_adjustment = 0
    roi_size_factor = 1
    scaling_factor = 1
    angle_adjustment = 0

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        cls.cheese = cls.model.from_zip(filename)
        cls.cheese.analyze(
            x_adjustment=cls.x_adjustment,
            y_adjustment=cls.y_adjustment,
            roi_size_factor=cls.roi_size_factor,
            scaling_factor=cls.scaling_factor,
            angle_adjustment=cls.angle_adjustment,
        )

    @classmethod
    def tearDownClass(cls):
        plt.close("all")
        del cls.cheese
        super().tearDownClass()

    def test_slice_locations(self):
        """Test the locations of the slices of interest."""
        self.assertAlmostEqual(self.cheese.origin_slice, self.origin_slice, delta=1)

    def test_phantom_roll(self):
        """Test the roll of the phantom."""
        self.assertAlmostEqual(self.cheese.catphan_roll, self.expected_roll, delta=0.3)

    def test_HU_values(self):
        """Test HU values."""
        for roi_num, expected_hu in self.hu_values.items():
            meas_hu = self.cheese.module.rois[roi_num].pixel_value
            self.assertAlmostEqual(expected_hu, meas_hu, delta=5)

    def test_pdf(self):
        save_file(self.cheese.publish_pdf, "temp")


class TestTomoCheeseDemo(CheeseMixin, PlotlyTestMixin, TestCase):
    origin_slice = 24
    expected_roll = -0.23
    hu_values = {
        "1": 16,
        "2": 20,
        "3": 23,
        "4": 2,
        "5": 16,
        "6": -669,
        "7": 15,
        "8": 25,
        "9": 653,
        "10": 25,
        "11": 24,
        "12": 102,
        "13": 7,
        "14": -930,
        "15": 23,
        "16": 14,
        "17": -516,
        "18": 447,
        "19": 269,
        "20": 14,
    }
    num_figs = 1
    fig_data = {
        0: {
            "title": "Tomo Cheese",
            "num_traces": 21,
        }
    }

    def setUp(self) -> None:
        self.instance = self.cheese

    @classmethod
    def setUpClass(cls):
        cls.cheese = TomoCheese.from_demo_images()
        cls.cheese.analyze()


@pytest.mark.proprietary
class TestHighHURodTomo(CheeseMixin, TestCase):
    dir_path = ["Tomo"]
    file_name = "High HU rod tomo.zip"
    origin_slice = 38
    expected_roll = -0.23
    hu_values = {
        "1": -501,
        "2": -24,
        "11": -616,
        "16": 437,
    }


@pytest.mark.proprietary
class TestCIRS062MErogluer(CheeseMixin, TestCase):
    model = CIRS062M
    origin_slice = 32
    dir_path = ["Tomo", "CIRS062M"]
    file_name = "CIRS062M - Erogluer.zip"
    expected_roll = 0
    hu_values = {
        "1": 5,
        "2": -55,
        "3": -811,
        "4": 66,
        "5": 73,
        "6": -20,
        "7": 1337,
        "8": 67,
        "9": 245,
        "10": -806,
        "11": 911,
        "12": -490,
        "13": 882,
        "14": -16,
        "15": 67,
        "16": -494,
        "17": 274,
    }


@pytest.mark.proprietary
class TestCIRS062MButsonCareDose(CheeseMixin, TestCase):
    model = CIRS062M
    origin_slice = 38
    dir_path = ["Tomo", "CIRS062M"]
    file_name = "Butson - caredose.zip"
    expected_roll = 0
    hu_values = {
        "1": -963,
        "2": -480,
        "3": 57,
        "4": 46,
        "5": 847,
        "6": -22,
        "7": 207,
        "8": -780,
        "9": -60,
        "10": 62,
        "11": 1291,
        "12": -20,
        "13": -60,
        "14": -476,
        "15": 44,
        "16": 214,
        "17": 856,
    }


@pytest.mark.proprietary
class TestCIRS062MButson160mAsB(CheeseMixin, TestCase):
    # same as caredose but different mAs, thus same HU values
    model = CIRS062M
    origin_slice = 38
    dir_path = ["Tomo", "CIRS062M"]
    file_name = "Butson - 160mAs.zip"
    expected_roll = 0
    hu_values = {
        "1": -963,
        "2": -480,
        "3": 57,
        "4": 46,
        "5": 847,
        "6": -22,
        "7": 207,
        "8": -780,
        "9": -60,
        "10": 62,
        "11": 1291,
        "12": -20,
        "13": -60,
        "14": -476,
        "15": 44,
        "16": 214,
        "17": 856,
    }


@pytest.mark.proprietary
class TestCIRS062MButsonLungOuter(CheeseMixin, TestCase):
    model = CIRS062M
    origin_slice = 24
    dir_path = ["Tomo", "CIRS062M"]
    file_name = "Butson - lung outer.zip"
    expected_roll = 0
    hu_values = {
        "1": -963,
        "2": -480,
        "3": 57,
        "4": 46,
        "5": 847,
        "6": -22,
        "7": 207,
        "8": -60,
        "9": -60,
        "10": 62,
        "11": 1291,
        "12": -20,
        "13": -780,
        "14": -476,
        "15": 44,
        "16": 214,
        "17": 856,
    }


@pytest.mark.proprietary
class TestCIRS062MButsonEmptyButoOuter(CheeseMixin, TestCase):
    model = CIRS062M
    origin_slice = 30
    dir_path = ["Tomo", "CIRS062M"]
    file_name = "Butson - empty outer.zip"
    expected_roll = 0
    hu_values = {
        "1": -60,
        "2": -480,
        "3": 57,
        "4": 46,
        "5": 847,
        "6": -22,
        "7": 207,
        "8": -780,
        "9": -60,
        "10": 62,
        "11": 1291,
        "12": -20,
        "13": -974,
        "14": -476,
        "15": 44,
        "16": 214,
        "17": 856,
    }


@pytest.mark.proprietary
class TestCIRS062MShifted(CheeseMixin, TestCase):
    model = CIRS062M
    origin_slice = 30
    dir_path = ["Tomo", "CIRS062M"]
    file_name = "CIRS bad, shifted.zip"
    expected_roll = -0.55
    hu_values = {
        "1": -3,
        "8": 45,
        "17": -489,
    }


@pytest.mark.proprietary
class TestCIRS062MApplyShifts(CheeseMixin, TestCase):
    # we apply some crazy shifts that will cause the ROIs to be off but also border
    # the ROI edges. This should make testing constancy robust since any small change
    # of the adjustments will cause a large difference in the results.
    model = CIRS062M
    origin_slice = 30
    dir_path = ["Tomo", "CIRS062M"]
    file_name = "CIRS bad, shifted.zip"
    expected_roll = 4.44
    hu_values = {
        "2": -741,  # air roi
        "3": -505,  # air roi
        "17": -345,
    }
    x_adjustment = 1
    y_adjustment = -2
    angle_adjustment = 5
    roi_size_factor = 1.2
    scaling_factor = 1.03
