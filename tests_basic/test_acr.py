import io
import os
import pytest
from pathlib import Path
from unittest import TestCase

from matplotlib import pyplot as plt
from scipy import ndimage

from pylinac import ACRMRILarge
from pylinac.acr import ACRCT, ACRCTResult, ACRMRIResult
from pylinac.core.geometry import Point
from pylinac.core.io import TemporaryZipDirectory
from tests_basic.core.test_utilities import QuaacTestBase, ResultsDataBase
from tests_basic.utils import (
    CloudFileMixin,
    FromZipTesterMixin,
    InitTesterMixin,
    PlotlyTestMixin,
    get_file_from_cloud_test_repo,
    save_file,
)

pytestmark = pytest.mark.proprietary

TEST_DIR_CT = ["ACR", "CT"]
TEST_DIR_MR = ["ACR", "MRI"]


class TestACRCT(TestCase, FromZipTesterMixin, InitTesterMixin):
    klass = ACRCT
    zip = [*TEST_DIR_CT, "Philips.zip"]
    init_file = [*TEST_DIR_CT, "Philips_CT"]
    is_folder = True

    def test_load_from_list_of_paths(self):
        # shouldn't raise
        with TemporaryZipDirectory(get_file_from_cloud_test_repo(self.zip)) as zfolder:
            paths = [Path(zfolder) / f for f in os.listdir(zfolder)]
            ACRCT(paths)

    def test_load_from_list_of_streams(self):
        # shouldn't raise
        with TemporaryZipDirectory(get_file_from_cloud_test_repo(self.zip)) as zfolder:
            paths = [Path(zfolder, f) for f in os.listdir(zfolder)]
            paths = [io.BytesIO(open(p, "rb").read()) for p in paths]
            ACRCT(paths)


class TestCTGeneral(TestCase):
    def setUp(self):
        self.path = get_file_from_cloud_test_repo([*TEST_DIR_CT, "Philips.zip"])
        self.ct = ACRCT.from_zip(self.path)

    def test_phan_center(self):
        """Test locations of the phantom center."""
        known_phan_center = Point(258, 253)
        self.ct.analyze()
        self.assertAlmostEqual(
            self.ct.ct_calibration_module.phan_center.x, known_phan_center.x, delta=0.7
        )
        self.assertAlmostEqual(
            self.ct.ct_calibration_module.phan_center.y, known_phan_center.y, delta=0.7
        )

    def test_results_data(self):
        self.ct.analyze()
        data = self.ct.results_data()
        self.assertIsInstance(data, ACRCTResult)
        self.assertEqual(data.num_images, self.ct.num_images)

        # check the additional modules got added
        self.assertIsInstance(data.ct_module.rois, dict)

    def test_results_warnings(self):
        self.ct.analyze()
        data = self.ct.results_data()
        self.assertEqual(len(data.warnings), 0)

    def test_lazy_is_same_as_default(self):
        self.ct.analyze()
        lazy_ct = ACRCT.from_zip(self.path, memory_efficient_mode=True)
        lazy_ct.analyze()
        self.assertEqual(self.ct.results(), lazy_ct.results())

    def test_passing_origin_slice(self):
        ct = ACRCT.from_zip(self.path)
        ct.analyze(origin_slice=3)  # automatic is 2
        self.assertEqual(ct.origin_slice, 3)


class TestPlottingSaving(TestCase):
    @classmethod
    def setUpClass(cls):
        path = get_file_from_cloud_test_repo([*TEST_DIR_CT, "Philips.zip"])
        cls.ct = ACRCT.from_zip(path)
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


class TestACRCTQuaac(QuaacTestBase, CloudFileMixin, TestCase):
    dir_path = ["ACR", "CT"]
    file_name = "Philips.zip"

    def quaac_instance(self):
        filename = self.get_filename()
        ct = ACRCT.from_zip(filename)
        ct.analyze()
        return ct


class ACRCTMixin(CloudFileMixin):
    dir_path = ["ACR", "CT"]
    origin_slice: int
    phantom_roll: float = 0
    mtf_50: float
    slice_thickness: float
    hu_values: dict
    unif_values: dict
    x_adjustment: float = 0
    y_adjustment: float = 0
    angle_adjustment: float = 0
    roi_size_factor: float = 1
    scaling_factor: float = 1

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        cls.ct = ACRCT.from_zip(filename)
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

    def test_mtf(self):
        self.assertAlmostEqual(
            self.ct.spatial_resolution_module.mtf.relative_resolution(50),
            self.mtf_50,
            delta=0.1,
        )

    def test_HU_values(self):
        """Test HU values."""
        for key, roi in self.ct.ct_calibration_module.rois.items():
            exp_val = self.hu_values[key]
            meas_val = roi.pixel_value
            self.assertAlmostEqual(exp_val, meas_val, delta=5)

    def test_uniformity_values(self):
        """Test Uniformity HU values."""
        for key, exp_val in self.unif_values.items():
            meas_val = self.ct.uniformity_module.rois[key].pixel_value
            self.assertAlmostEqual(exp_val, meas_val, delta=5)


class ACRPhilips(ACRCTMixin, PlotlyTestMixin, TestCase):
    file_name = "Philips.zip"
    mtf_50 = 0.54
    phantom_roll = -0.3
    hu_values = {"Poly": -87, "Acrylic": 126, "Bone": 904, "Air": -987, "Water": 4}
    unif_values = {"Center": 1, "Left": 1, "Right": 1, "Top": 1, "Bottom": 1}
    num_figs = 6
    fig_data = {
        0: {
            "title": "HU Linearity",
            "num_traces": 6,
        },
        1: {
            "title": "HU Uniformity",
            "num_traces": 6,
        },
        2: {
            "title": "Spatial Resolution",
            "num_traces": 9,
        },
        3: {
            "title": "Low Contrast",
            "num_traces": 3,
        },
        4: {
            "title": "Relative MTF",
            "num_traces": 1,
            "x_label": "Line pairs / mm",
            "y_label": "Relative MTF",
        },
        5: {
            "title": "Side View",
            "num_traces": 5,
        },
    }

    def setUp(self) -> None:
        self.instance = self.ct


class ACRCTApplyROIOffset(ACRPhilips):
    # crazy shifts to create a unique set of results that's sensitive to changes
    x_adjustment = 1
    y_adjustment = -2
    angle_adjustment = 5
    roi_size_factor = 1.5
    scaling_factor = 1.02
    mtf_50 = 46
    phantom_roll = 4.75
    hu_values = {"Poly": -79, "Acrylic": 119, "Bone": 891, "Air": -971, "Water": 4}


class ACRPhilipsOffset(ACRCTMixin, TestCase):
    """Shift the phantom over by several pixels to ensure no row/col algorithm issues

    Unfortunately, I can't move them that far because the FOV is very tight
    """

    file_name = "Philips.zip"
    mtf_50 = 0.54
    phantom_roll = -0.3
    hu_values = {"Poly": -87, "Acrylic": 126, "Bone": 904, "Air": -987, "Water": 4}
    unif_values = {"Center": 1, "Left": 1, "Right": 1, "Top": 1, "Bottom": 1}

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        cls.ct = ACRCT.from_zip(filename)
        for img in cls.ct.dicom_stack:
            img.roll(direction="x", amount=5)
        cls.ct.analyze()


class ACRPhilipsRotated(ACRCTMixin, TestCase):
    """Rotate the phantom over by several pixels to ensure no row/col algorithm issues

    Unfortunately, I can't move them that far because the FOV is very tight
    """

    file_name = "Philips.zip"
    mtf_50 = 0.54
    phantom_roll = -3.3
    hu_values = {"Poly": -87, "Acrylic": 126, "Bone": 904, "Air": -987, "Water": 4}
    unif_values = {"Center": 1, "Left": 1, "Right": 1, "Top": 1, "Bottom": 1}

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        cls.ct = ACRCT.from_zip(filename)
        for img in cls.ct.dicom_stack:
            img.array = ndimage.rotate(img.array, angle=3, mode="nearest")
        cls.ct.analyze()


class TestACRMRI(TestCase, FromZipTesterMixin, InitTesterMixin):
    klass = ACRMRILarge
    zip = [*TEST_DIR_MR, "GE 3T.zip"]
    init_file = [*TEST_DIR_MR, "GE - 3T"]
    is_folder = True

    def test_load_from_list_of_paths(self):
        # shouldn't raise
        with TemporaryZipDirectory(get_file_from_cloud_test_repo(self.zip)) as zfolder:
            paths = [Path(zfolder) / f for f in os.listdir(zfolder)]
            ACRMRILarge(paths)

    def test_load_from_list_of_streams(self):
        # shouldn't raise
        with TemporaryZipDirectory(get_file_from_cloud_test_repo(self.zip)) as zfolder:
            paths = [Path(zfolder, f) for f in os.listdir(zfolder)]
            paths = [io.BytesIO(open(p, "rb").read()) for p in paths]
            ACRMRILarge(paths)

    def test_geometric_distortion_profile_lengths(self):
        path = get_file_from_cloud_test_repo([*TEST_DIR_MR, "GE 3T.zip"])
        mri = ACRMRILarge.from_zip(path)
        mri.analyze()
        data = mri.results_data()
        self.assertAlmostEqual(
            data.geometric_distortion_module.profiles["negative diagonal"][
                "width (mm)"
            ],
            190.63,
            delta=0.05,
        )
        self.assertAlmostEqual(
            data.geometric_distortion_module.profiles["positive diagonal"][
                "width (mm)"
            ],
            190.59,
            delta=0.05,
        )
        self.assertAlmostEqual(
            data.geometric_distortion_module.profiles["vertical"]["width (mm)"],
            190.44,
            delta=0.05,
        )
        self.assertAlmostEqual(
            data.geometric_distortion_module.profiles["horizontal"]["width (mm)"],
            190.44,
            delta=0.05,
        )


class TestACRMRIResultData(TestCase, ResultsDataBase):
    def construct_analyzed_instance(self):
        path = get_file_from_cloud_test_repo([*TEST_DIR_MR, "GE 3T.zip"])
        mri = ACRMRILarge.from_zip(path)
        mri.analyze()
        return mri

    def test_row_mtf_keys(self):
        phantom = self.construct_analyzed_instance()
        data = phantom.results_data()
        for key, value in {10: 1.13, 80: 0.564}.items():
            self.assertAlmostEqual(data.slice1.row_mtf_lp_mm[key], value, delta=0.01)

    def test_col_mtf_keys(self):
        phantom = self.construct_analyzed_instance()
        data = phantom.results_data()
        for key, value in {10: 1.11, 80: 0.826}.items():
            self.assertAlmostEqual(data.slice1.col_mtf_lp_mm[key], value, delta=0.01)


class TestMRGeneral(TestCase):
    def setUp(self):
        path = get_file_from_cloud_test_repo([*TEST_DIR_MR, "GE 3T.zip"])
        self.mri = ACRMRILarge.from_zip(path)

    def test_phan_center(self):
        """Test locations of the phantom center."""
        known_phan_center = Point(129, 127)
        self.mri.analyze()
        self.assertAlmostEqual(
            self.mri.slice1.phan_center.x, known_phan_center.x, delta=0.7
        )
        self.assertAlmostEqual(
            self.mri.slice1.phan_center.y, known_phan_center.y, delta=0.7
        )

    def test_results_data(self):
        self.mri.analyze()
        data = self.mri.results_data()
        self.assertIsInstance(data, ACRMRIResult)
        self.assertEqual(data.num_images, self.mri.num_images)

        # check the additional modules got added
        self.assertIsInstance(data.slice11.rois, dict)

    def test_results_warnings(self):
        self.mri.analyze()
        data = self.mri.results_data()
        self.assertEqual(len(data.warnings), 0)

    def test_echo_number(self):
        """Test analyzing a specific echo number works"""
        path = get_file_from_cloud_test_repo([*TEST_DIR_MR, "AXIAL_DUAL_ECHO.zip"])
        mri1 = ACRMRILarge.from_zip(path)
        mri1.analyze(echo_number=1)
        echo_number_1 = mri1.dicom_stack[0].metadata.EchoNumbers
        path = get_file_from_cloud_test_repo([*TEST_DIR_MR, "AXIAL_DUAL_ECHO.zip"])
        mri2 = ACRMRILarge.from_zip(path)
        mri2.analyze(echo_number=2)
        echo_number_2 = mri2.dicom_stack[0].metadata.EchoNumbers
        self.assertNotEqual(echo_number_1, echo_number_2)

    def test_echo_number_invalid(self):
        path = get_file_from_cloud_test_repo([*TEST_DIR_MR, "AXIAL_DUAL_ECHO.zip"])
        mri = ACRMRILarge.from_zip(path)
        with self.assertRaises(ValueError):
            mri.analyze(echo_number=3)  # only 2 echoes

    def test_echo_number_defaults_to_first(self):
        path = get_file_from_cloud_test_repo([*TEST_DIR_MR, "AXIAL_DUAL_ECHO.zip"])
        mri = ACRMRILarge.from_zip(path)
        self.assertEqual(len({s.metadata.EchoNumbers for s in mri.dicom_stack}), 2)
        mri.analyze(echo_number=None)
        self.assertEqual(mri.dicom_stack[0].metadata.EchoNumbers, "1")

    def test_config_extent_rounds(self):
        """Test that the extent check rounds the config extent to the nearest slice"""
        path = get_file_from_cloud_test_repo([*TEST_DIR_MR, "Config rounding.zip"])
        mri = ACRMRILarge.from_zip(path)
        self.assertTrue(mri._ensure_physical_scan_extent())


class TestMRPlottingSaving(TestCase):
    @classmethod
    def setUpClass(cls):
        path = get_file_from_cloud_test_repo([*TEST_DIR_MR, "GE 3T.zip"])
        cls.mri = ACRMRILarge.from_zip(path)
        cls.mri.analyze()

    @classmethod
    def tearDownClass(cls):
        plt.close("all")
        del cls.mri

    def test_plot_images(self):
        """Test that saving an image does something."""
        save_file(self.mri.plot_images)

    def test_save_images(self):
        """Test that saving an image does something."""
        save_file(self.mri.save_images, to_single_file=False)

    def test_subimages_errors(self):
        """We don't use subimages here. easier to pass as a list of figs"""
        with self.assertRaises(NotImplementedError):
            self.mri.plot_analyzed_subimage("sr")
        with self.assertRaises(NotImplementedError):
            self.mri.save_analyzed_subimage("sr")

    def test_set_figure_size(self):
        self.mri.plot_analyzed_image(figsize=(8, 13))
        fig = plt.gcf()
        self.assertEqual(fig.bbox_inches.height, 13)
        self.assertEqual(fig.bbox_inches.width, 8)


class TestACRMRIQuaac(QuaacTestBase, CloudFileMixin, TestCase):
    dir_path = ["ACR", "MRI"]
    file_name = "T1-Single.zip"

    def quaac_instance(self):
        filename = self.get_filename()
        ct = ACRMRILarge.from_zip(filename)
        ct.analyze()
        return ct


class ACRMRMixin(CloudFileMixin):
    dir_path = ["ACR", "MRI"]
    phantom_roll: float = 0
    row_mtf_50: float
    col_mtf_50: float
    slice_thickness: float
    slice1_shift: float
    slice11_shift: float
    psg: float
    results: list[str] = []
    x_adjustment: float = 0
    y_adjustment: float = 0
    angle_adjustment: float = 0
    roi_size_factor: float = 1
    scaling_factor: float = 1

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        cls.mri = ACRMRILarge.from_zip(filename, memory_efficient_mode=True)
        cls.mri.analyze(
            x_adjustment=cls.x_adjustment,
            y_adjustment=cls.y_adjustment,
            angle_adjustment=cls.angle_adjustment,
            roi_size_factor=cls.roi_size_factor,
            scaling_factor=cls.scaling_factor,
        )

    @classmethod
    def tearDownClass(cls):
        plt.close("all")
        del cls.mri
        super().tearDownClass()

    def test_roll(self):
        self.assertAlmostEqual(self.mri.catphan_roll, self.phantom_roll, delta=0.3)

    def test_mtf(self):
        self.assertAlmostEqual(
            self.mri.slice1.row_mtf.relative_resolution(50), self.row_mtf_50, delta=0.1
        )
        self.assertAlmostEqual(
            self.mri.slice1.col_mtf.relative_resolution(50), self.col_mtf_50, delta=0.1
        )

    def test_slice_thickness(self):
        print(self.mri.slice1.measured_slice_thickness_mm)
        print(self.mri.slice1.thickness_rois["Top"].wire_fwhm)
        print(self.mri.slice1.thickness_rois["Bottom"].wire_fwhm)
        print(self.mri.slice1.mm_per_pixel)
        self.assertAlmostEqual(
            self.mri.slice1.measured_slice_thickness_mm, self.slice_thickness, delta=0.5
        )

    def test_slice1_shift(self):
        self.assertAlmostEqual(
            self.mri.slice1.slice_shift_mm, self.slice1_shift, delta=0.2
        )

    def test_slice11_shift(self):
        self.assertAlmostEqual(
            self.mri.slice11.slice_shift_mm, self.slice11_shift, delta=0.2
        )

    def test_psg(self):
        self.assertAlmostEqual(self.mri.uniformity_module.psg, self.psg, delta=0.3)

    def test_results(self):
        results = self.mri.results()
        for result in self.results:
            self.assertIn(result, results)


class ACRT1Single(ACRMRMixin, PlotlyTestMixin, TestCase):
    file_name = "T1-Single.zip"
    row_mtf_50 = 0.96
    col_mtf_50 = 0.96
    phantom_roll = -0.5
    slice_thickness = 5
    slice1_shift = -1
    slice11_shift = 0
    psg = 0.3
    num_figs = 6
    fig_data = {
        0: {"title": "Slice 1 (Thickness, Offset, Resolution)", "num_traces": 13},
        1: {"title": "Geometric Distortion", "num_traces": 5},
        5: {
            "title": "Relative MTF",
            "num_traces": 2,
            "x_label": "Line pairs / mm",
            "y_label": "Relative MTF",
        },
    }

    def setUp(self) -> None:
        self.instance = self.mri


class ACRDualEcho(ACRMRMixin, TestCase):
    file_name = "AXIAL_DUAL_ECHO.zip"
    row_mtf_50 = 0.96
    col_mtf_50 = 0.96
    phantom_roll = 0
    slice_thickness = 5
    slice1_shift = -1
    slice11_shift = 0
    psg = 0.3


class ACRDualEcho2(ACRMRMixin, TestCase):
    file_name = "AXIAL_DUAL_ECHO.zip"
    row_mtf_50 = 0.96
    col_mtf_50 = 0.96
    phantom_roll = 0
    slice_thickness = 4.4
    slice1_shift = -1
    slice11_shift = 0
    psg = 0.3

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        cls.mri = ACRMRILarge.from_zip(filename)
        cls.mri.analyze(echo_number=2)


class ACRUVMSliceLocation(ACRMRMixin, TestCase):
    """Test that the image patient position is the primary measurement of slice location."""

    file_name = "UVM_slice_location_mismatch.zip"
    row_mtf_50 = 0.96
    col_mtf_50 = 0.96
    phantom_roll = 0
    slice_thickness = 4.5
    slice1_shift = 0.5
    slice11_shift = 2
    psg = 0.3


class ACRGE3T(ACRMRMixin, TestCase):
    file_name = "GE 3T.zip"
    row_mtf_50 = 0.96
    col_mtf_50 = 0.96
    phantom_roll = 0
    slice_thickness = 5
    slice1_shift = 0
    slice11_shift = 1.5
    psg = 0.3
    results = ["5.50mm"]


class ACRGEROIOffsetsApplied(ACRGE3T):
    # crazy shifts to create a unique set of results that's sensitive to changes
    x_adjustment = 1
    y_adjustment = -2
    angle_adjustment = 5
    roi_size_factor = 1.2
    scaling_factor = 1.02
    slice1_shift = 3.9
    slice11_shift = 5.86
    col_mtf_50 = 1.06
    phantom_roll = 2.83
    results = ["5.04mm"]


class ACRGE3TOffset(ACRGE3T):
    """Shift the phantom over by several pixels to ensure no row/col algorithm issues

    Unfortunately, I can't move them that far because the FOV is very tight
    """

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        cls.mri = ACRMRILarge.from_zip(filename)
        for img in cls.mri.dicom_stack:
            img.roll(direction="x", amount=5)
        cls.mri.analyze()


class ACRGE3TRotated(ACRGE3T):
    """Rotate the phantom over by a bit. Sadly, this does mess up the algorithm slightly as
    many ROIs are rectangles and cannot be truly rotated by the determined roll.
    Adding this test so for constancy but also so that in the future if the
    ROI analysis is improved this test can be fixed.
    """

    phantom_roll = -0.4
    results = ["4.80mm"]  # induced rotation does change this a bit. See above.

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        cls.mri = ACRMRILarge.from_zip(filename)
        for img in cls.mri.dicom_stack:
            img.array = ndimage.rotate(img.array, angle=0.5, mode="nearest")
        cls.mri.analyze()


class ACRBadProfile(ACRMRMixin, TestCase):
    """This dataset was causing long_profile errors due to the profile not
    being a true profile w/ low values at the edges.

    There are also other issues going on (different geometry module).
    The only point here is to test the long profile error."""

    file_name = "Config rounding.zip"
    row_mtf_50 = 0.96
    col_mtf_50 = 0.96
    phantom_roll = 0
    slice_thickness = 4.2
    slice1_shift = 1.5
    slice11_shift = 1.95
    psg = 0.3


class ACRMRIUnfilled(ACRMRMixin, TestCase):
    # a phantom with significant air bubbles and also a rotation
    # results aren't great, but are provided for benchmarking
    file_name = "Sola.zip"
    row_mtf_50 = 1
    col_mtf_50 = 0.1
    phantom_roll = 0.2
    slice_thickness = 3
    slice1_shift = -0.5
    slice11_shift = 1.5
    psg = 0.3

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        cls.mri = ACRMRILarge.from_zip(filename)
        for img in cls.mri.dicom_stack:
            img.array = ndimage.rotate(img.array, angle=1, mode="nearest")
        cls.mri.analyze()
