import io
import json
import os
import os.path as osp
import tempfile
from unittest import TestCase, skip

import matplotlib.pyplot as plt
import numpy as np

from pylinac import CatPhan503, CatPhan504, CatPhan600, CatPhan604
from pylinac.core.geometry import Point
from pylinac.core.io import TemporaryZipDirectory
from pylinac.ct import CTP404CP503, CTP404CP504, CTP528CP503, CTP528CP504, CatphanResult
from tests_basic.utils import (
    CloudFileMixin,
    FromDemoImageTesterMixin,
    FromURLTesterMixin,
    FromZipTesterMixin,
    InitTesterMixin,
    get_file_from_cloud_test_repo,
    save_file,
)

TEST_DIR = "CBCT"


class TestInstantiation(
    TestCase,
    InitTesterMixin,
    FromDemoImageTesterMixin,
    FromURLTesterMixin,
    FromZipTesterMixin,
):
    klass = CatPhan504
    init_file = [TEST_DIR, "Pelvis"]
    demo_load_method = "from_demo_images"
    url = "CatPhan504.zip"
    zip = [TEST_DIR, "CBCT_4.zip"]
    is_folder = True

    def test_load_from_list_of_paths(self):
        # shouldn't raise
        with TemporaryZipDirectory(get_file_from_cloud_test_repo(self.zip)) as zfolder:
            paths = [osp.join(zfolder, f) for f in os.listdir(zfolder)]
            CatPhan504(paths)

    def test_load_from_list_of_streams(self):
        # shouldn't raise
        with TemporaryZipDirectory(get_file_from_cloud_test_repo(self.zip)) as zfolder:
            paths = [osp.join(zfolder, f) for f in os.listdir(zfolder)]
            paths = [io.BytesIO(open(p, "rb").read()) for p in paths]
            CatPhan504(paths)


class TestGeneral(TestCase):
    """Test general things when using cbct module."""

    def setUp(self):
        self.cbct = CatPhan504.from_demo_images()

    def test_load_from_stream(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, "CBCT_4.zip"])
        ref_cbct = CatPhan504.from_zip(path)
        ref_cbct.analyze()
        with open(path, "rb") as f:
            s = io.BytesIO(f.read())
            cbct = CatPhan504.from_zip(s)
            cbct.analyze()
        self.assertIsInstance(cbct, CatPhan504)
        self.assertEqual(cbct.origin_slice, ref_cbct.origin_slice)

    def test_load_from_file_object(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, "CBCT_4.zip"])
        ref_cbct = CatPhan504.from_zip(path)
        ref_cbct.analyze()
        with open(path, "rb") as f:
            cbct = CatPhan504.from_zip(f)
            cbct.analyze()
        self.assertIsInstance(cbct, CatPhan504)
        self.assertEqual(cbct.origin_slice, ref_cbct.origin_slice)

    def test_scan_extent_not_long_enough_far_side(self):
        """Test that if a scan doesn't include all the modules it raises an error"""
        path = get_file_from_cloud_test_repo([TEST_DIR, "CBCT_4.zip"])
        ref_cbct = CatPhan504.from_zip(path)
        # artificially chop the dicom stack
        ref_cbct.dicom_stack.images = ref_cbct.dicom_stack.images[
            :-25
        ]  # chop off the back
        ref_cbct.dicom_stack.metadatas = ref_cbct.dicom_stack.metadatas[
            :-25
        ]  # chop off the back
        with self.assertRaises(ValueError):
            ref_cbct.localize()

    def test_scan_extent_not_long_enough_near_side(self):
        """Test that if a scan doesn't include all the modules it raises an error"""
        path = get_file_from_cloud_test_repo([TEST_DIR, "CBCT_4.zip"])
        ref_cbct = CatPhan504.from_zip(path)
        # artificially chop the dicom stack
        ref_cbct.dicom_stack.images = ref_cbct.dicom_stack.images[
            20:
        ]  # chop off the front
        ref_cbct.dicom_stack.metadatas = ref_cbct.dicom_stack.metadatas[
            20:
        ]  # chop off the front
        with self.assertRaises(ValueError):
            ref_cbct.localize()

    def test_crop_before_analysis(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, "CBCT_4.zip"])
        cbct = CatPhan504.from_zip(path)
        for img in cbct.dicom_stack:
            img.crop(pixels=20, edges=("bottom",))
        # shouldn't raise
        cbct.analyze()

    def test_demo(self):
        """Run the demo to make sure it works."""
        self.cbct.run_demo()

    def test_helpers(self):
        """Test the various helper methods."""
        self.cbct.analyze()
        self.cbct._results()

    def test_phan_center(self):
        """Test locations of the phantom center."""
        known_phan_center = Point(257, 255)
        self.cbct.analyze()
        self.assertAlmostEqual(
            self.cbct.ctp404.phan_center.x, known_phan_center.x, delta=0.7
        )
        self.assertAlmostEqual(
            self.cbct.ctp404.phan_center.y, known_phan_center.y, delta=0.7
        )

    def test_results(self):
        self.cbct.analyze()
        data = self.cbct.results()
        self.assertIsInstance(data, str)

        data_list = self.cbct.results(as_list=True)
        self.assertIsInstance(data_list, list)
        self.assertIsInstance(data_list[0], list)
        self.assertIsInstance(data_list[0][0], str)
        self.assertEqual(len(data_list), 4)

    def test_results_data(self):
        self.cbct.analyze()
        data = self.cbct.results_data()
        self.assertIsInstance(data, CatphanResult)
        self.assertEqual(data.num_images, self.cbct.num_images)

        # check the additional modules got added
        self.assertAlmostEqual(data.ctp528.start_angle_radians, np.pi, delta=0.02)
        self.assertEqual(
            data.ctp486.nps_max_freq, self.cbct.ctp486.max_noise_power_frequency
        )
        self.assertEqual(data.ctp486.nps_avg_power, self.cbct.ctp486.avg_noise_power)

        for p in range(10, 91, 10):
            self.assertEqual(
                data.ctp528.mtf_lp_mm[p], self.cbct.ctp528.mtf.relative_resolution(p)
            )

        data_dict = self.cbct.results_data(as_dict=True)
        self.assertIsInstance(data_dict, dict)

        data_json = self.cbct.results_data(as_json=True)
        self.assertIsInstance(data_json, str)
        # shouldn't raise
        json.loads(data_json)

    def test_contrast_str(self):
        # shouldn't raise
        self.cbct.analyze(contrast_method="Michelson")
        self.cbct.results_data()

    def test_same_results_for_lazy_load(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, "CBCT_4.zip"])
        ct = CatPhan504.from_zip(path, memory_efficient_mode=False)
        ct.analyze()
        lazy_ct = CatPhan504.from_zip(path, memory_efficient_mode=True)
        lazy_ct.analyze()
        self.assertEqual(ct.results(), lazy_ct.results())


class TestCustomPhantom(TestCase):
    def test_removing_module(self):
        class CatPhan504Modified(CatPhan504):
            modules = {
                CTP404CP504: {"offset": 0},
            }

        ct = CatPhan504Modified.from_demo_images()  # shouldn't raise
        ct.analyze()

    def test_omitting_404_module(self):
        class CatPhan504Modified(CatPhan504):
            modules = {
                CTP528CP504: {"offset": -30},
            }

        ct = CatPhan504Modified.from_demo_images()  # shouldn't raise
        with self.assertRaises(ValueError):
            ct.analyze()

    def test_modifying_module(self):
        class CTP528Custom(CTP528CP503):
            start_angle = 0.1

        class CatPhan503Modified(CatPhan503):
            modules = {
                CTP404CP503: {"offset": 0},
                CTP528Custom: {"offset": -30},
            }

        ct = CatPhan503Modified.from_demo_images()  # shouldn't raise
        ct.analyze()
        ct.plot_analyzed_image()


class TestDemos(TestCase):
    def test_504(self):
        CatPhan504.run_demo()

    def test_503(self):
        CatPhan503.run_demo()

    def test_600(self):
        CatPhan600.run_demo()

    def test_604(self):
        CatPhan604.run_demo()


class TestCustomHUValues(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cbct = CatPhan504.from_demo_images()

    def test_override(self):
        self.cbct.analyze(expected_hu_values={"Air": -1000, "Poly": 321})
        self.assertEqual(self.cbct.ctp404.rois["Air"].nominal_val, -1000)
        self.assertEqual(self.cbct.ctp404.rois["Poly"].nominal_val, 321)

    def test_extra_keys_harmless(self):
        # shouldn't raise
        self.cbct.analyze(expected_hu_values={"stuffcicles": 1111})

    def test_nothing_passed(self):
        # shouldn't raise; normal; backwards-compatible
        self.cbct.analyze(expected_hu_values=None)
        self.assertEqual(self.cbct.ctp404.rois["Air"].nominal_val, -1000)


class TestPlottingSaving(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cbct = CatPhan504.from_demo_images()
        cls.cbct.analyze()

    @classmethod
    def tearDownClass(cls):
        plt.close("all")

    def test_save_image(self):
        """Test that saving an image does something."""
        for method in ["save_analyzed_image", "save_analyzed_subimage"]:
            methodcall = getattr(self.cbct, method)
            save_file(methodcall)

    def test_plot_images(self):
        """Test the various plotting functions."""
        self.cbct.plot_analyzed_image()

    def test_plot_subimages(self):
        for item in ["hu", "un", "mtf", "sp", "prof", "lin", "lc", "side"]:
            self.cbct.plot_analyzed_subimage(item)

        self.cbct.plot_analyzed_subimage("lin", delta=False)

        with self.assertRaises(ValueError):
            self.cbct.plot_analyzed_subimage("sr")

    def test_save_subimages(self):
        with tempfile.TemporaryDirectory() as tdir:
            file = os.path.join(tdir, "dummy.png")
            for item in ["hu", "un", "mtf", "sp", "prof", "lin", "lc"]:
                self.cbct.save_analyzed_subimage(file, item)

            self.cbct.save_analyzed_subimage(file, "lin", delta=False)

            with self.assertRaises(ValueError):
                self.cbct.save_analyzed_subimage(file, "sr")

    def test_set_figure_size(self):
        self.cbct.plot_analyzed_image(figsize=(8, 13))
        fig = plt.gcf()
        self.assertEqual(fig.bbox_inches.height, 13)
        self.assertEqual(fig.bbox_inches.width, 8)


class CatPhanMixin(CloudFileMixin):
    """A mixin to use for testing Varian CBCT scans; does not inherit from TestCase as it would be run
    otherwise."""

    catphan = CatPhan504
    check_uid = True
    origin_slice = 0
    dir_path = ["CBCT"]
    hu_tolerance = 40
    scaling_tolerance = 1
    zip = True
    expected_roll = 0
    hu_values = {}
    unif_values = {}
    mtf_values = {}
    avg_line_length = 50
    slice_thickness = 2
    thickness_slice_straddle = "auto"
    lowcon_visible = 0
    memory_efficient = True
    print_debug = False
    avg_noise_power = None
    max_noise_power_frequency = 0

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        if cls.zip:
            cls.cbct = cls.catphan.from_zip(
                filename, memory_efficient_mode=cls.memory_efficient
            )
        else:
            cls.cbct = cls.catphan(filename, memory_efficient_mode=cls.memory_efficient)
        cls.cbct.analyze(
            cls.hu_tolerance,
            cls.scaling_tolerance,
            thickness_slice_straddle=cls.thickness_slice_straddle,
        )
        if cls.print_debug:
            print(cls.cbct._results())

    @classmethod
    def tearDownClass(cls):
        plt.close("all")
        super().tearDownClass()

    def test_slice_thickness(self):
        """Test the slice thickness."""
        self.assertAlmostEqual(
            self.cbct.ctp404.meas_slice_thickness, self.slice_thickness, delta=0.3
        )

    def test_lowcontrast_bubbles(self):
        """Test the number of low contrast bubbles visible."""
        if not isinstance(self.cbct, CatPhan503):
            self.assertAlmostEqual(
                self.cbct.ctp515.rois_visible, self.lowcon_visible, delta=1
            )

    def test_slice_locations(self):
        """Test the locations of the slices of interest."""
        self.assertAlmostEqual(self.cbct.origin_slice, self.origin_slice, delta=1)

    def test_phantom_roll(self):
        """Test the roll of the phantom."""
        self.assertAlmostEqual(self.cbct.catphan_roll, self.expected_roll, delta=0.3)

    def test_HU_values(self):
        """Test HU values."""
        for key, roi in self.cbct.ctp404.rois.items():
            exp_val = self.hu_values[key]
            meas_val = roi.pixel_value
            self.assertAlmostEqual(exp_val, meas_val, delta=5)

    def test_uniformity_values(self):
        """Test Uniformity HU values."""
        for key, exp_val in self.unif_values.items():
            meas_val = self.cbct.ctp486.rois[key].pixel_value
            self.assertAlmostEqual(exp_val, meas_val, delta=5)

    def test_geometry_line_length(self):
        """Test the geometry distances."""
        self.assertAlmostEqual(
            self.avg_line_length, self.cbct.ctp404.avg_line_length, delta=0.1
        )

    def test_avg_noise_power(self):
        if self.avg_noise_power:
            self.assertAlmostEqual(
                self.avg_noise_power, self.cbct.ctp486.avg_noise_power, delta=0.1
            )

    def test_max_noise_frequency(self):
        if self.max_noise_power_frequency:
            self.assertEqual(
                self.cbct.ctp486.max_noise_power_frequency,
                self.max_noise_power_frequency,
            )

    def test_MTF_values(self):
        """Test MTF values."""
        for key, exp_mtf in self.mtf_values.items():
            meas_mtf = self.cbct.ctp528.mtf.relative_resolution(key)
            self.assertAlmostEqual(exp_mtf, meas_mtf, delta=0.1)

    def test_pdf(self):
        save_file(self.cbct.publish_pdf, "temp")


class CatPhanDemo(CatPhanMixin, TestCase):
    """Test the CBCT demo (Varian high quality head protocol)."""

    expected_roll = -0.3
    origin_slice = 32
    hu_values = {
        "Poly": -45,
        "Acrylic": 117,
        "Delrin": 341,
        "Air": -998,
        "Teflon": 997,
        "PMP": -200,
        "LDPE": -103,
    }
    unif_values = {"Center": 17, "Left": 10, "Right": 0, "Top": 6, "Bottom": 6}
    mtf_values = {50: 0.56}
    avg_line_length = 49.92
    lowcon_visible = 3  # changed w/ visibility refactor in v3.0
    slice_thickness = 2.5
    delete_file = False

    @classmethod
    def setUpClass(cls):
        cls.cbct = CatPhan504.from_demo_images()
        cls.cbct.analyze()


class CatPhan4(CatPhanMixin, TestCase):
    """A Varian CBCT dataset"""

    file_name = "CBCT_4.zip"
    expected_roll = -2.57
    origin_slice = 31
    hu_values = {
        "Poly": -33,
        "Acrylic": 119,
        "Delrin": 335,
        "Air": -979,
        "Teflon": 970,
        "PMP": -185,
        "LDPE": -94,
    }
    unif_values = {"Center": 17, "Left": 10, "Right": 22, "Top": 13, "Bottom": 18}
    mtf_values = {80: 0.24, 90: 0.2, 60: 0.31, 70: 0.23, 95: 0.15}
    lowcon_visible = 1
    slice_thickness = 2.4


class Elekta2(CatPhanMixin, TestCase):
    """An Elekta CBCT dataset"""

    catphan = CatPhan503
    file_name = "Elekta_2.zip"
    origin_slice = 162
    hu_values = {
        "Poly": -319,
        "Acrylic": -224,
        "Delrin": -91,
        "Air": -863,
        "Teflon": 253,
        "PMP": -399,
        "LDPE": -350,
    }
    unif_values = {
        "Center": -285,
        "Left": -279,
        "Right": -278,
        "Top": -279,
        "Bottom": -279,
    }
    mtf_values = {80: 0.28, 90: 0.22, 60: 0.37, 70: 0.31, 95: 0.18}
    slice_thickness = 1

    def test_plot_low_contrast_is_none(self):
        """There is no low contrast section to plot so ensure it doesn't"""
        self.assertIsNone(self.cbct.plot_analyzed_subimage("lc"))

    def test_save_low_contrast_is_none(self):
        """There is no low contrast section to plot so ensure it doesn't"""
        with io.BytesIO() as bio:
            self.assertIsNone(
                self.cbct.save_analyzed_subimage(filename=bio, subimage="lc")
            )


class CatPhan600_2(CatPhanMixin, TestCase):
    """An Elekta CBCT dataset"""

    catphan = CatPhan600
    file_name = "zzCAT201602.zip"
    expected_roll = -0.64
    origin_slice = 34
    hu_values = {
        "Poly": -29,
        "Acrylic": 123,
        "Delrin": 336,
        "Air": -932,
        "Teflon": 897,
        "PMP": -164,
        "LDPE": -80,
    }
    hu_passed = False
    unif_values = {"Center": 14, "Left": 15, "Right": 15, "Top": 16, "Bottom": 13}
    mtf_values = {50: 0.4}
    avg_line_length = 50.02
    slice_thickness = 4.5
    lowcon_visible = 2  # changed w/ visibility refactor in v3.0


class CatPhan600WaterVial(CatPhanMixin, TestCase):
    """Scan with the water vial in place"""

    catphan = CatPhan600
    dir_path = ["CBCT", "CatPhan_600"]
    file_name = "Catphan600_water_vial.zip"
    expected_roll = -0.25
    origin_slice = 96
    hu_values = {
        "Poly": -55,
        "Acrylic": 115,
        "Delrin": 339,
        "Air": -1000,
        "Teflon": 960,
        "PMP": -200,
        "LDPE": -123,
        "Vial": 17,
    }
    hu_passed = False
    unif_values = {"Center": -20, "Left": -10, "Right": -13, "Top": -3, "Bottom": -12}
    mtf_values = {50: 0.49}
    avg_line_length = 50.07
    slice_thickness = 2.1
    lowcon_visible = 1

    def test_vial_roi(self):
        self.assertIn("Vial", self.cbct.ctp404.rois)


class CatPhan600WaterVial2(CatPhanMixin, TestCase):
    """Scan with the water vial in place. This broke the MTF algo previously."""

    catphan = CatPhan600
    dir_path = ["CBCT", "CatPhan_600"]
    file_name = "Catphan_600_mtf.zip"
    expected_roll = -0.25
    origin_slice = 104
    hu_values = {
        "Poly": -15,
        "Acrylic": 166,
        "Delrin": 400,
        "Air": -1000,
        "Teflon": 1108,
        "PMP": -186,
        "LDPE": -95,
        "Vial": 63,
    }
    hu_passed = False
    unif_values = {"Center": 14, "Left": 4, "Right": 8, "Top": 24, "Bottom": 10}
    mtf_values = {50: 0.58}
    avg_line_length = 50.07
    slice_thickness = 2.1
    lowcon_visible = 1

    def test_vial_roi(self):
        self.assertIn("Vial", self.cbct.ctp404.rois)


class CatPhan604Test(CatPhanMixin, TestCase):
    catphan = CatPhan604
    file_name = "CBCTCatPhan604.zip"
    origin_slice = 47
    hu_values = {
        "Poly": -47,
        "Acrylic": 105,
        "Delrin": 338,
        "Air": -981,
        "Teflon": 942,
        "PMP": -194,
        "LDPE": -105,
        "50% Bone": 771,
        "20% Bone": 263,
    }
    unif_values = {"Center": -3, "Left": 0, "Right": 0, "Top": 0, "Bottom": 0}
    mtf_values = {50: 0.43}
    lowcon_visible = 1  # changed w/ visibility refactor in v3.0
    avg_noise_power = 0.252


class CatPhan504Mixin(CatPhanMixin):
    catphan = CatPhan504
    dir_path = [TEST_DIR, "CatPhan_504"]


class CatPhan503Mixin(CatPhanMixin):
    catphan = CatPhan503
    dir_path = [TEST_DIR, "CatPhan_503"]


class CatPhan600Mixin(CatPhanMixin):
    catphan = CatPhan600
    dir_path = [TEST_DIR, "CatPhan_600"]


class CatPhan604Mixin(CatPhanMixin):
    catphan = CatPhan604
    dir_path = [TEST_DIR, "CatPhan_604"]


class VarianPelvis(CatPhan504Mixin, TestCase):
    """Test the Varian Pelvis protocol CBCT."""

    file_name = "Pelvis.zip"
    expected_roll = -0.26
    origin_slice = 32
    hu_values = {
        "Poly": -36,
        "Acrylic": 114,
        "Delrin": 345,
        "Air": -998,
        "Teflon": 992,
        "PMP": -188,
        "LDPE": -95,
    }
    unif_values = {"Center": 17, "Left": 5, "Right": 4, "Top": 4, "Bottom": 4}
    mtf_values = {30: 0.42, 50: 0.34, 80: 0.24}
    avg_line_length = 49.8
    lowcon_visible = 3
    slice_thickness = 2.5


class VarianPelvisSpotlight(CatPhan504Mixin, TestCase):
    """Test the Varian Pelvis Spotlight protocol CBCT."""

    file_name = "Pelvis spotlight.zip"
    expected_roll = -0.26
    origin_slice = 32
    hu_values = {
        "Poly": -43,
        "Acrylic": 118,
        "Delrin": 341,
        "Air": -998,
        "Teflon": 967,
        "PMP": -198,
        "LDPE": -100,
    }
    unif_values = {"Center": 19, "Left": 3, "Right": -1, "Top": -1, "Bottom": 0}
    mtf_values = {50: 0.55}
    avg_line_length = 49.94
    lowcon_visible = 5
    slice_thickness = 2.4


class VarianLowDoseThorax(CatPhan504Mixin, TestCase):
    """Test the Varian Low-Dose Thorax protocol CBCT."""

    file_name = "Low dose thorax.zip"
    expected_roll = -0.29
    origin_slice = 32
    hu_values = {
        "Poly": -42,
        "Acrylic": 119,
        "Delrin": 341,
        "Air": -998,
        "Teflon": 992,
        "PMP": -191,
        "LDPE": -94,
    }
    unif_values = {"Center": 16, "Left": 7, "Right": -1, "Top": 3, "Bottom": 2}
    mtf_values = {50: 0.3}
    avg_line_length = 49.7
    lowcon_visible = 0
    slice_thickness = 2.4


class VarianStandardHead(CatPhan504Mixin, TestCase):
    """Test the Varian Standard Head protocol CBCT."""

    file_name = "Standard head.zip"
    expected_roll = -0.19
    origin_slice = 32
    hu_values = {
        "Poly": -40,
        "Acrylic": 127,
        "Delrin": 350,
        "Air": -997,
        "Teflon": 997,
        "PMP": -191,
        "LDPE": -101,
    }
    unif_values = {"Center": 17, "Left": 15, "Right": 4, "Top": 9, "Bottom": 9}
    mtf_values = {50: 0.54}
    avg_line_length = 49.94
    lowcon_visible = 1
    slice_thickness = 2.4


class VarianLowDoseHead(CatPhan504Mixin, TestCase):
    """Test the Varian Low-Dose Head protocol CBCT."""

    file_name = "Low dose head.zip"
    expected_roll = -0.4
    origin_slice = 32
    hu_values = {
        "Poly": -41,
        "Acrylic": 121,
        "Delrin": 350,
        "Air": -997,
        "Teflon": 1003,
        "PMP": -197,
        "LDPE": -103,
    }
    unif_values = {"Center": 13, "Left": 11, "Right": 3, "Top": 7, "Bottom": 6}
    mtf_values = {50: 0.55}
    lowcon_visible = 1
    avg_line_length = 49.93
    slice_thickness = 1.3


@skip  # SR slice fails; TODO: make slices modular, so if 1 fails, the others still run
class GEMonthlyCT(CatPhan504Mixin, TestCase):
    """Test a monthly CT scan from GE."""

    file_name = "GE_CT.zip"
    hu_tolerance = 90
    origin_slice = 143
    hu_values = {
        "Poly": -32,
        "Acrylic": 119,
        "Delrin": 333,
        "Air": -944,
        "Teflon": 909,
        "PMP": -173,
        "LDPE": -87,
    }
    unif_values = {"Center": 11, "Left": 11, "Right": 11, "Top": 11, "Bottom": 11}
    mtf_values = {80: 0.2}
    lowcon_visible = 4


class ToshibaMonthlyCT(CatPhan504Mixin, TestCase):
    """Test a monthly CT scan from Toshiba."""

    file_name = "Toshiba.zip"
    origin_slice = 36
    hu_values = {
        "Poly": -32,
        "Acrylic": 106,
        "Delrin": 467,
        "Air": -992,
        "Teflon": 1207,
        "PMP": -165,
        "LDPE": -85,
    }
    unif_values = {"Center": 8, "Left": 7, "Right": 7, "Top": 7, "Bottom": 6}
    mtf_values = {50: 0.23}
    lowcon_visible = 4
    slice_thickness = 2.8


class CBCT1(CatPhan504Mixin, TestCase):
    """A Varian CBCT dataset"""

    file_name = "CBCT_1.zip"
    expected_roll = -0.53
    origin_slice = 32
    hu_values = {
        "Poly": -35,
        "Acrylic": 130,
        "Delrin": 347,
        "Air": -996,
        "Teflon": 1004,
        "PMP": -186,
        "LDPE": -94,
    }
    unif_values = {"Center": 13, "Left": 17, "Right": 5, "Top": 10, "Bottom": 9}
    mtf_values = {50: 0.48}
    avg_line_length = 49.9
    lowcon_visible = 1
    slice_thickness = 1.77

    def test_plot_low_contrast(self):
        """There is no low contrast section to plot so ensure it doesn't"""
        self.assertIsInstance(self.cbct.plot_analyzed_subimage("lc"), plt.Figure)

    def test_save_low_contrast(self):
        """There is no low contrast section to plot so ensure it doesn't"""
        with io.BytesIO() as bio:
            self.assertIsInstance(
                self.cbct.save_analyzed_subimage(filename=bio, subimage="lc"),
                plt.Figure,
            )


class CBCT2(CatPhan504Mixin, TestCase):
    """A Varian CBCT dataset"""

    file_name = "CBCT_2.zip"
    expected_roll = -0.3
    origin_slice = 34
    hu_values = {
        "Poly": -16,
        "Acrylic": 135,
        "Delrin": 367,
        "Air": -965,
        "Teflon": 1017,
        "PMP": -163,
        "LDPE": -71,
    }
    unif_values = {"Center": 47, "Left": 35, "Right": 37, "Top": 36, "Bottom": 37}
    mtf_values = {50: 0.34}
    lowcon_visible = 2
    avg_line_length = 49.9
    slice_thickness = 2.4
    avg_noise_power = 0.395


class CBCT3(CatPhan504Mixin, TestCase):
    """A Varian CBCT dataset"""

    file_name = "CBCT_3.zip"
    expected_roll = -2.66
    origin_slice = 36
    hu_values = {
        "Poly": -44,
        "Acrylic": 113,
        "Delrin": 325,
        "Air": -982,
        "Teflon": 952,
        "PMP": -194,
        "LDPE": -103,
    }
    unif_values = {"Center": -2, "Left": -3, "Right": 6, "Top": -2, "Bottom": 6}
    mtf_values = {50: 0.33}
    avg_line_length = 49.9
    lowcon_visible = 1
    slice_thickness = 2.4


# CBCT4 is in the regular test_cbct.py file


class CBCT5(CatPhan504Mixin, TestCase):
    """A Varian CBCT dataset"""

    file_name = "CBCT_5.zip"
    origin_slice = 34
    hu_values = {
        "Poly": -56,
        "Acrylic": 101,
        "Delrin": 328,
        "Air": -999,
        "Teflon": 977,
        "PMP": -201,
        "LDPE": -110,
    }
    unif_values = {"Center": 19, "Left": -10, "Right": -5, "Top": -7, "Bottom": -8}
    mtf_values = {50: 0.33}
    avg_line_length = 49.55
    lowcon_visible = 5  # changed w/ visibility refactor in v3.0
    slice_thickness = 2.35


class CBCT6(CatPhan504Mixin, TestCase):
    """A Varian CBCT dataset"""

    file_name = "CBCT_6.zip"
    expected_roll = 0.2
    origin_slice = 38
    hu_values = {
        "Poly": -44,
        "Acrylic": 107,
        "Delrin": 327,
        "Air": -994,
        "Teflon": 972,
        "PMP": -192,
        "LDPE": -100,
    }
    unif_values = {"Center": -3, "Left": -3, "Right": -13, "Top": -7, "Bottom": -6}
    mtf_values = {50: 0.54}
    avg_line_length = 49.94
    lowcon_visible = 4  # changed w/ visibility refactor in v3.0
    slice_thickness = 2.35


class CBCT7(CatPhan504Mixin, TestCase):
    """A Varian CBCT dataset"""

    file_name = "CBCT_7.zip"
    expected_roll = 0.5
    origin_slice = 36
    hu_values = {
        "Poly": -48,
        "Acrylic": 108,
        "Delrin": 331,
        "Air": -999,
        "Teflon": 984,
        "PMP": -198,
        "LDPE": -107,
    }
    unif_values = {"Center": 12, "Left": -7, "Right": -7, "Top": -8, "Bottom": -7}
    mtf_values = {50: 0.33}
    avg_line_length = 49.6
    lowcon_visible = 4  # changed w/ visibility refactor in v3.0
    slice_thickness = 2.3


class CBCT8(CatPhan504Mixin, TestCase):
    """A Varian CBCT dataset"""

    file_name = "CBCT_8.zip"
    use_classifier = False
    expected_roll = 0.55
    origin_slice = 40
    hu_values = {
        "Poly": -37,
        "Acrylic": 114,
        "Delrin": 334,
        "Air": -994,
        "Teflon": 982,
        "PMP": -186,
        "LDPE": -97,
    }
    unif_values = {"Center": -4, "Left": 2, "Right": -5, "Top": 0, "Bottom": -1}
    mtf_values = {50: 0.53}
    avg_line_length = 49.95
    lowcon_visible = 3
    slice_thickness = 2.3


class CBCT9(CatPhan504Mixin, TestCase):
    """A Varian CBCT dataset"""

    file_name = "CBCT_9.zip"
    expected_roll = 0.4
    origin_slice = 35
    hu_values = {
        "Poly": -53,
        "Acrylic": 107,
        "Delrin": 330,
        "Air": -999,
        "Teflon": 980,
        "PMP": -199,
        "LDPE": -107,
    }
    unif_values = {"Center": 10, "Left": -8, "Right": -7, "Top": -6, "Bottom": -5}
    mtf_values = {50: 0.34}
    avg_line_length = 49.6
    lowcon_visible = 3
    slice_thickness = 2.4


class CBCT10(CatPhan504Mixin, TestCase):
    """A Varian CBCT dataset"""

    file_name = "CBCT_10.zip"
    expected_roll = 0.4
    origin_slice = 38
    hu_values = {
        "Poly": -37,
        "Acrylic": 109,
        "Delrin": 334,
        "Air": -992,
        "Teflon": 985,
        "PMP": -186,
        "LDPE": -93,
    }
    unif_values = {"Center": -1, "Left": 4, "Right": -5, "Top": 0, "Bottom": -1}
    mtf_values = {50: 0.53}
    lowcon_visible = 4
    slice_thickness = 2.3


class CBCT11(CatPhan504Mixin, TestCase):
    """A Varian CBCT dataset"""

    file_name = "CBCT_11.zip"
    expected_roll = 0.5
    origin_slice = 38
    hu_values = {
        "Poly": -37,
        "Acrylic": 108,
        "Delrin": 332,
        "Air": -994,
        "Teflon": 982,
        "PMP": -187,
        "LDPE": -95,
    }
    unif_values = {"Center": -3, "Left": 2, "Right": -7, "Top": -2, "Bottom": -3}
    mtf_values = {50: 0.53}
    avg_line_length = 49.94
    lowcon_visible = 4
    slice_thickness = 2.35
    avg_noise_power = 0.226


class CBCT12(CatPhan504Mixin, TestCase):
    """A Varian CBCT dataset"""

    file_name = "CBCT_12.zip"
    origin_slice = 36
    hu_values = {
        "Poly": -55,
        "Acrylic": 112,
        "Delrin": 329,
        "Air": -999,
        "Teflon": 982,
        "PMP": -201,
        "LDPE": -107,
    }
    unif_values = {"Center": 3, "Left": -5, "Right": -11, "Top": -9, "Bottom": -8}
    mtf_values = {50: 0.23}
    avg_line_length = 49.59
    lowcon_visible = 2
    slice_thickness = 2.35


class CBCT13(CatPhan504Mixin, TestCase):
    """A Varian CBCT dataset"""

    file_name = "CBCT_13.zip"
    expected_roll = 0.2
    origin_slice = 36
    hu_values = {
        "Poly": -51,
        "Acrylic": 106,
        "Delrin": 329,
        "Air": -999,
        "Teflon": 976,
        "PMP": -198,
        "LDPE": -107,
    }
    unif_values = {"Center": 3, "Left": -9, "Right": -9, "Top": -10, "Bottom": -8}
    mtf_values = {50: 0.3}
    avg_line_length = 49.66
    lowcon_visible = 5  # changed w/ visibility refactor in v3.0
    slice_thickness = 2.5


class CBCT14(CatPhan504Mixin, TestCase):
    """A Varian CBCT dataset"""

    file_name = "CBCT_14.zip"
    expected_roll = 0.84
    origin_slice = 32
    hu_values = {
        "Poly": -41,
        "Acrylic": 125,
        "Delrin": 334,
        "Air": -995,
        "Teflon": 986,
        "PMP": -184,
        "LDPE": -89,
    }
    unif_values = {"Center": 18, "Left": 13, "Right": 15, "Top": 14, "Bottom": 14}
    mtf_values = {50: 0.3}
    avg_line_length = 49.4
    lowcon_visible = 1
    slice_thickness = 2.4


class CBCT15(CatPhan504Mixin, TestCase):
    """A Varian CBCT dataset."""

    file_name = "CBCT_15.zip"
    origin_slice = 61
    hu_values = {
        "Poly": -32,
        "Acrylic": 121,
        "Delrin": 353,
        "Air": -998,
        "Teflon": 945,
        "PMP": -186,
        "LDPE": -93,
    }
    unif_values = {"Center": -2, "Left": 6, "Right": 5, "Top": 3, "Bottom": 11}
    mtf_values = {50: 0.35}
    lowcon_visible = 6


class CBCT16(CatPhan504Mixin, TestCase):
    """A Varian CBCT dataset"""

    file_name = "CBCT_16.zip"
    expected_roll = 0.2
    origin_slice = 32
    hu_values = {
        "Poly": -37,
        "Acrylic": 128,
        "Delrin": 342,
        "Air": -995,
        "Teflon": 1000,
        "PMP": -181,
        "LDPE": -87,
    }
    unif_values = {"Center": 17, "Left": 20, "Right": 18, "Top": 19, "Bottom": 19}
    mtf_values = {50: 0.31}
    avg_line_length = 49.7
    lowcon_visible = 1
    slice_thickness = 2.4


class CBCT17(CatPhan504Mixin, TestCase):
    """A Varian CBCT dataset"""

    file_name = "CBCT_17.zip"
    expected_roll = 0.45
    origin_slice = 34
    hu_values = {
        "Poly": -43,
        "Acrylic": 114,
        "Delrin": 336,
        "Air": -995,
        "Teflon": 989,
        "PMP": -192,
        "LDPE": -101,
    }
    unif_values = {"Center": 11, "Left": 0, "Right": -7, "Top": -6, "Bottom": -2}
    mtf_values = {50: 0.57}
    avg_line_length = 49.9
    lowcon_visible = 1


class Catphan504Ring(CatPhan504Mixin, TestCase):
    """A Varian CBCT dataset with ring artifact"""

    file_name = "ringartifact.zip"
    expected_roll = 0.1
    origin_slice = 32
    hu_values = {
        "Poly": -52,
        "Acrylic": 103,
        "Delrin": 321,
        "Air": -995,
        "Teflon": 950,
        "PMP": -201,
        "LDPE": -105,
    }
    unif_values = {"Center": -2, "Left": -4, "Right": -16, "Top": -12, "Bottom": -8}
    mtf_values = {50: 0.44}
    avg_line_length = 50.1
    lowcon_visible = 1
    slice_thickness = 2.5


class Katy1(CatPhan504Mixin, TestCase):
    """CBCT with very high HU values."""

    file_name = "Katy-iX-Monday, March 10, 2014 1-05-47 PM (super high HU).zip"
    origin_slice = 44
    hu_values = {
        "Poly": 584,
        "Acrylic": 797,
        "Delrin": 1046,
        "Air": -404,
        "Teflon": 1720,
        "PMP": 424,
        "LDPE": 522,
    }
    unif_values = {"Center": 612, "Left": 628, "Right": 645, "Top": 631, "Bottom": 642}
    mtf_values = {50: 0.51}
    lowcon_visible = 0  # changed w/ visibility refactor in v3.0
    slice_thickness = 2.4


class CTWithCloseCouch(CatPhan503Mixin, TestCase):
    """A CT where the couch is super close"""

    file_name = "CT with close couch.zip"
    expected_roll = 0.43
    origin_slice = 133
    hu_values = {
        "Poly": -48,
        "Acrylic": 87,
        "Delrin": 272,
        "Air": -809,
        "Teflon": 789,
        "PMP": -145,
        "LDPE": -87,
    }
    unif_values = {"Center": -4, "Left": -12.5, "Right": -15, "Top": 3, "Bottom": -23}
    mtf_values = {50: 0.44}
    avg_line_length = 49.8
    slice_thickness = 1.42


class AGElekta1(CatPhan503Mixin, TestCase):
    """An Elekta CBCT dataset"""

    file_name = "AG-Elekta-F0-M20-mA25.zip"
    expected_roll = -1.2
    origin_slice = 189
    hu_values = {
        "Poly": 341,
        "Acrylic": 447,
        "Delrin": 590,
        "Air": -339,
        "Teflon": 1039,
        "PMP": 241,
        "LDPE": 308,
    }
    unif_values = {"Center": 397, "Left": 355, "Right": 346, "Top": 348, "Bottom": 350}
    mtf_values = {50: 0.19}
    avg_line_length = 50.2
    slice_thickness = 1.15


class AGElekta2(CatPhan503Mixin, TestCase):
    """An Elekta CBCT dataset"""

    file_name = "AG-Elekta-F1-S10-mA10.zip"
    expected_roll = -0.44
    origin_slice = 189
    hu_values = {
        "Poly": 746,
        "Acrylic": 834,
        "Delrin": 969,
        "Air": 181,
        "Teflon": 1324,
        "PMP": 670,
        "LDPE": 722,
    }
    unif_values = {"Center": 707, "Left": 758, "Right": 748, "Top": 750, "Bottom": 758}
    mtf_values = {50: 0.22}
    avg_line_length = 50.2
    slice_thickness = 1


class Elekta4(CatPhan503Mixin, TestCase):
    """An Elekta CBCT dataset"""

    file_name = "Elekta_4.zip"
    origin_slice = 129
    hu_values = {
        "Poly": -319,
        "Acrylic": -224,
        "Delrin": -91,
        "Air": -863,
        "Teflon": 255,
        "PMP": -401,
        "LDPE": -352,
    }
    unif_values = {
        "Center": -273,
        "Left": -267,
        "Right": -266,
        "Top": -267,
        "Bottom": -267,
    }
    mtf_values = {50: 0.35}
    slice_thickness = 1


class Elekta7(CatPhan503Mixin, TestCase):
    """An Elekta CBCT dataset"""

    file_name = "Elekta_7.zip"
    origin_slice = 159
    hu_values = {
        "Poly": -128,
        "Acrylic": -16,
        "Delrin": 141,
        "Air": -782,
        "Teflon": 541,
        "PMP": -226,
        "LDPE": -164,
    }
    unif_values = {"Center": -81, "Left": -73, "Right": -72, "Top": -73, "Bottom": -73}
    mtf_values = {50: 0.4}
    slice_thickness = 1


class Elekta8(CatPhan503Mixin, TestCase):
    """An Elekta CBCT dataset"""

    file_name = "Elekta_8.zip"
    origin_slice = 161
    hu_values = {
        "Poly": -324,
        "Acrylic": -229,
        "Delrin": -97,
        "Air": -868,
        "Teflon": 245,
        "PMP": -406,
        "LDPE": -356,
    }
    unif_values = {
        "Center": -293,
        "Left": -286,
        "Right": -285,
        "Top": -286,
        "Bottom": -286,
    }
    mtf_values = {50: 0.35}
    slice_thickness = 1


class UNC100kV(CatPhan503Mixin, TestCase):
    file_name = "UNC-100kV_CBCT_Feb2016.zip"
    origin_slice = 131
    hu_values = {
        "Poly": -112,
        "Acrylic": 35,
        "Delrin": 245,
        "Air": -973,
        "Teflon": 856,
        "PMP": -239,
        "LDPE": -168,
    }
    unif_values = {"Center": -45, "Left": -72, "Right": -73, "Top": -61, "Bottom": -71}
    mtf_values = {50: 0.55}
    slice_thickness = 1.4


class UNC120kV(CatPhan503Mixin, TestCase):
    file_name = "UNC-120kV_CBCT_Feb2016.zip"
    origin_slice = 131
    hu_values = {
        "Poly": -274,
        "Acrylic": -150,
        "Delrin": 22,
        "Air": -996,
        "Teflon": 486,
        "PMP": -380,
        "LDPE": -315,
    }
    unif_values = {
        "Center": -223,
        "Left": -236,
        "Right": -238,
        "Top": -232,
        "Bottom": -239,
    }
    mtf_values = {50: 0.56}
    slice_thickness = 1.5


class CatPhan600_1(CatPhan600Mixin, TestCase):
    file_name = "zzCAT201601.zip"
    expected_roll = -1.1
    origin_slice = 158
    hu_values = {
        "Poly": -31,
        "Acrylic": 124,
        "Delrin": 344,
        "Air": -958,
        "Teflon": 921,
        "PMP": -173,
        "LDPE": -87,
    }
    unif_values = {"Center": 14, "Left": 14, "Right": 13, "Top": 14, "Bottom": 12}
    mtf_values = {50: 0.22}
    avg_line_length = 50.1
    slice_thickness = 1.25
    lowcon_visible = 1  # after roll adjustment 1 more pixel is included, pulling down contrast detection

    def test_no_vial_roi(self):
        self.assertNotIn("Vial", self.cbct.ctp404.rois)


class TOHPhilips2mm(CatPhan504Mixin, TestCase):
    file_name = "H TOH - Philips120kV2mm.zip"
    expected_roll = -0.2
    origin_slice = 61
    hu_values = {
        "Poly": -30,
        "Acrylic": 130,
        "Delrin": 292,
        "Air": -1000,
        "Teflon": 874,
        "PMP": -178,
        "LDPE": -87,
    }
    unif_values = {"Center": 22, "Left": 20, "Right": 20, "Top": 20, "Bottom": 20}
    mtf_values = {50: 0.35}
    avg_line_length = 49.9
    lowcon_visible = 6


class DBCatPhan503Roll(CatPhan503Mixin, TestCase):
    file_name = "DenisBrojan-Catphan503_MFOV_resolution_problem.zip"
    expected_roll = -0.12
    origin_slice = 79
    hu_values = {
        "Poly": -95,
        "Acrylic": 12,
        "Delrin": 166,
        "Air": -740,
        "Teflon": 564,
        "PMP": -189,
        "LDPE": -130,
    }
    unif_values = {"Center": -75, "Left": -60, "Right": -57, "Top": -56, "Bottom": -56}
    mtf_values = {50: 0.25}
    avg_line_length = 48.9
    lowcon_visible = 6


@skip("Is not a 504; figure out what it is")
class Philips1MM(CatPhan504Mixin, TestCase):
    use_classifier = False
    file_name = "Philips 1mm.zip"
    expected_roll = -0.49
    origin_slice = 284
    hu_values = {
        "Poly": -32,
        "Acrylic": 125,
        "Delrin": 343,
        "Air": -967,
        "Teflon": 914,
        "PMP": -175,
        "LDPE": -85,
    }
    unif_values = {"Center": 95, "Left": 96, "Right": 96, "Top": 99, "Bottom": 94}
    mtf_values = {50: 0.48}
    avg_line_length = 49.8
    lowcon_visible = 0


@skip("Is not a 504; figure out what it is")
class Siemens5MM(CatPhan504Mixin, TestCase):
    use_classifier = False
    check_uid = False
    file_name = "Siemens 5mm.zip"
    expected_roll = 1.26
    origin_slice = 41
    hu_values = {
        "Poly": -40,
        "Acrylic": 120,
        "Delrin": 347,
        "Air": -1003,
        "Teflon": 953,
        "PMP": -187,
        "LDPE": -95,
    }
    unif_values = {"Center": 47, "Left": 50, "Right": 49, "Top": 53, "Bottom": 49}
    mtf_values = {50: 0.65}
    lowcon_visible = 0


@skip("Investigate later")
class CatPhan604Sen(CatPhan604Mixin, TestCase):
    expected_roll = -0.5
    file_name = "SiemensSenCatPhan604.zip"
    origin_slice = 191
    hu_values = {
        "Poly": -44,
        "Acrylic": 122,
        "Delrin": 392,
        "Air": -1000,
        "Teflon": 1061,
        "PMP": -220,
        "LDPE": -113,
        "50% Bone": 700,
        "20% Bone": 230,
    }
    unif_values = {"Center": 10, "Left": 10, "Right": 10, "Top": 10, "Bottom": 10}
    mtf_values = {50: 0.26}
    slice_thickness = 0.25
    lowcon_visible = 2


class CatPhan604Som(CatPhan604Mixin, TestCase):
    file_name = "SiemensSomCatPhan604.zip"
    origin_slice = 179
    hu_values = {
        "Poly": -34,
        "Acrylic": 120,
        "Delrin": 343,
        "Air": -1000,
        "Teflon": 928,
        "PMP": -183,
        "LDPE": -92,
        "50% Bone": 633,
        "20% Bone": 216,
    }
    unif_values = {"Center": 10, "Left": 10, "Right": 10, "Top": 10, "Bottom": 10}
    mtf_values = {50: 0.35}
    lowcon_visible = 3  # changed w/ visibility refactor in v3.0
    slice_thickness = 1


class CatPhan604wJig(CatPhan604Mixin, TestCase):
    file_name = "Catphan604-with-jig.zip"
    origin_slice = 43
    avg_line_length = 49.9
    hu_values = {
        "Poly": -19,
        "Acrylic": 130,
        "Delrin": 371,
        "Air": -999,
        "Teflon": 975,
        "PMP": -179,
        "LDPE": -91,
        "50% Bone": 717,
        "20% Bone": 246,
    }
    unif_values = {"Center": 6, "Left": 4, "Right": 13, "Top": 10, "Bottom": 7}
    mtf_values = {50: 0.28}
    lowcon_visible = 1
    avg_noise_power = 0.287


class CatPhan604wJig2(CatPhan604Mixin, TestCase):
    file_name = "Catphan604wJig2.zip"
    expected_roll = -0.61
    origin_slice = 49
    slice_thickness = 1.95
    avg_line_length = 49.9
    hu_values = {
        "Poly": -46,
        "Acrylic": 120,
        "Delrin": 361,
        "Air": -999,
        "Teflon": 965,
        "PMP": -192,
        "LDPE": -107,
        "50% Bone": 691,
        "20% Bone": 233,
    }
    unif_values = {"Center": 10, "Left": -2, "Right": 9, "Top": 7, "Bottom": 1}
    mtf_values = {50: 0.28}
    lowcon_visible = 2


class CatPhan503SliceOverlap(CatPhan503Mixin, TestCase):
    """This dataset is 6mm slices with a 2mm overlap. The slice thickness
    default algorithm gives incorrect values unless the straddle is set explicitly"""

    file_name = "catphan_slice_overlap.zip"
    expected_roll = 0.1
    origin_slice = 115
    hu_values = {
        "Poly": -48,
        "Acrylic": 123,
        "Delrin": 327,
        "Air": -1010,
        "Teflon": 901,
        "PMP": -190,
        "LDPE": -105,
    }
    thickness_slice_straddle = 0
    slice_thickness = 5.9
    unif_values = {"Center": 8, "Left": 7, "Right": 8, "Top": 8, "Bottom": 8}
    mtf_values = {50: 0.40}


class CatPhan604NegativeSliceOverlap(CatPhan604Mixin, TestCase):
    """Has a negative value for slice overlap. Has to do with stack order, but is irrelevant for our needs:
    https://dicom.innolitics.com/ciods/nm-image/nm-reconstruction/00180088"""

    file_name = "negative_spacing.zip"
    expected_roll = 0.1
    origin_slice = 65
    hu_values = {
        "Poly": -20,
        "Acrylic": 128,
        "Delrin": 351,
        "Air": -967,
        "Teflon": 910,
        "PMP": -165,
        "LDPE": -75,
        "50% Bone": 604,
        "20% Bone": 211,
    }
    thickness_slice_straddle = 0
    slice_thickness = 2.15
    unif_values = {"Center": 24, "Left": 23, "Right": 22, "Top": 23, "Bottom": 21}
    mtf_values = {50: 0.40}
    lowcon_visible = 6
