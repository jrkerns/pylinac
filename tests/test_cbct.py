import unittest
import os
import os.path as osp
import time

import numpy as np
from dicom.errors import InvalidDicomError

from pylinac.cbct import CBCT
from pylinac.core.geometry import Point


varian_test_file_dir = osp.join(osp.dirname(__file__), 'test_files', 'CBCT', 'Varian')
other_test_file_dir = osp.join(osp.dirname(__file__), 'test_files', 'CBCT')


class GeneralTests(unittest.TestCase):
    """Test general things when using cbct module."""

    def setUp(self):
        self.cbct = CBCT()

    def test_demo(self):
        """Run the demo to make sure it works."""
        self.cbct.run_demo()

    def test_loading(self):
        """Test various loading schemes."""
        # load demo images
        CBCT.from_demo_images()

        # load from folder directly
        folder = osp.join(varian_test_file_dir, 'Pelvis')
        CBCT(folder)

        # load from zip file
        zfile = osp.join(varian_test_file_dir, 'Pelvis.zip')
        CBCT.from_zip_file(zfile)

    def test_images_not_loaded(self):
        """Raise error if trying to analyze when images aren't loaded yet."""
        self.assertRaises(AttributeError, self.cbct.analyze)

    def test_bad_files(self):
        """Test bad file inputs."""
        not_a_folder = "notafolder"
        self.assertRaises(NotADirectoryError, self.cbct.load_folder, not_a_folder)

        not_a_zip = "notazip.zip"
        self.assertRaises(FileExistsError, self.cbct.load_zip_file, not_a_zip)

        not_image_folder = osp.join(osp.dirname(__file__), 'core')
        self.assertRaises(FileNotFoundError, self.cbct.load_folder, not_image_folder)

        no_CT_images_zip = osp.join(varian_test_file_dir, 'dummy.zip')
        self.assertRaises(FileNotFoundError, self.cbct.load_zip_file, no_CT_images_zip)

    def test_images_not_from_same_study(self):
        """Loading images from different studies should raise and error."""
        mixed_zip = osp.join(varian_test_file_dir, 'mixed_studies.zip')
        with self.assertRaises(InvalidDicomError):
            self.cbct.load_zip_file(mixed_zip)

    def test_phan_center(self):
        """Test locations of the phantom center."""
        self.cbct.load_demo_images()

        known_phan_center = Point(257, 255)
        self.cbct._construct_HU()
        self.assertAlmostEqual(self.cbct.HU.phan_center.x, known_phan_center.x, delta=0.7)
        self.assertAlmostEqual(self.cbct.HU.phan_center.y, known_phan_center.y, delta=0.7)

        # test a shifted image set
        shifted_phan_center = Point(287, 255)
        self.cbct.settings.images = np.roll(self.cbct.settings.images, 30, axis=1)
        self.cbct._construct_HU()
        self.assertAlmostEqual(self.cbct.HU.phan_center.x, shifted_phan_center.x, delta=0.7)
        self.assertAlmostEqual(self.cbct.HU.phan_center.y, shifted_phan_center.y, delta=0.7)

    @unittest.skip
    def test_finding_HU_slice(self):
        """Test the robustness of the algorithm to find the HU linearity slice."""
        self.cbct.load_demo_images()

        self.assertEqual(self.cbct.settings.hu_slice_num, 32)

        # roll the phantom data by 4 slices
        np.roll(self.cbct.settings.images, 4, axis=2)

    def test_save_image(self):
        """Test that saving an image does something."""
        filename = 'saved_img.jpg'

        self.cbct.load_demo_images()
        self.cbct.analyze(hu_tolerance=10, scaling_tolerance=0.01)
        for method in ['save_analyzed_image', 'save_analyzed_subimage']:
            methodcall = getattr(self.cbct, method)
            methodcall(filename)
            time.sleep(0.1)  # sleep just to let OS work
            self.assertTrue(osp.isfile(filename), "Save file did not successfully save the image")

            # cleanup
            os.remove(filename)
            self.assertFalse(osp.isfile(filename), "Save file test did not clean up saved image")

    def test_plot_images(self):
        """Test the various plotting functions."""
        self.cbct.load_demo_images()
        self.cbct.analyze()

        self.cbct.plot_analyzed_image()
        for item in ['hu', 'unif', 'mtf', 'sr']:
            self.cbct.plot_analyzed_subimage(item)

        with self.assertRaises(ValueError):
            self.cbct.plot_analyzed_subimage('srunif')


class CBCTMixin:
    """A mixin to use for Varian CBCT scans; does not inherit from TestCase as it would be run
        otherwise."""
    hu_tolerance = 40
    scaling_tolerance = 1
    zip = True
    read_all = False

    @classmethod
    def setUpClass(cls):
        if cls.zip:
            cls.cbct = CBCT.from_zip_file(cls.location, read_all=cls.read_all)
        else:
            cls.cbct = CBCT(cls.location, read_all=cls.read_all)

    def test_all_passed(self):
        """Test the pass flags for all tests."""
        self.cbct.analyze(self.hu_tolerance, self.scaling_tolerance)
        self.assertTrue(self.cbct.HU.overall_passed)
        self.assertTrue(self.cbct.UN.overall_passed)
        self.assertTrue(self.cbct.GEO.overall_passed)

    def test_slice_locations(self):
        """Test the locations of the slices of interest."""
        for attr, slice_name in zip(('hu_slice_num', 'un_slice_num', 'sr_slice_num', 'lc_slice_num'), ('HU', 'UN', 'SR', 'LC')):
            self.assertEqual(getattr(self.cbct.settings, attr), self.slice_locations[slice_name])

    def test_phantom_roll(self):
        """Test the roll of the demo phantom."""
        self.assertAlmostEqual(self.cbct.settings.phantom_roll, self.expected_roll, delta=0.001)

    def test_HU_values(self):
        """Test HU values."""
        self.cbct.analyze()
        for key, roi in self.cbct.HU.ROIs.items():
            exp_val = self.hu_values[key]
            meas_val = roi.pixel_value
            self.assertAlmostEqual(exp_val, meas_val, delta=5)

    def test_uniformity_values(self):
        """Test Uniformity HU values."""
        self.cbct.analyze()
        for key, roi in self.cbct.UN.ROIs.items():
            exp_val = self.unif_values[key]
            meas_val = roi.pixel_value
            self.assertAlmostEqual(exp_val, meas_val, delta=5)

    def test_geometry_line_lengths(self):
        """Test the geometry distances."""
        self.cbct.analyze()
        for key, line in self.cbct.GEO.lines.items():
            exp_dist = self.line_lengths[key]
            meas_dist = line.length_mm(self.cbct.settings.mm_per_pixel)
            self.assertAlmostEqual(exp_dist, meas_dist, delta=0.08)

    def test_MTF_values(self):
        """Test MTF values."""
        self.cbct.analyze()
        for key, exp_mtf in self.mtf_values.items():
            meas_mtf = self.cbct.SR.get_MTF(key)
            self.assertAlmostEqual(exp_mtf, meas_mtf, delta=0.05)


class CBCTDemo(CBCTMixin, unittest.TestCase):
    """Test the CBCT demo (Varian high quality head protocol)."""
    expected_roll = 0.004
    slice_locations = {'HU': 32, 'UN': 3, 'SR': 44, 'LC': 20}
    hu_values = {'Poly': -45, 'Acrylic': 117, 'Delrin': 341, 'Air': -998, 'Teflon': 997, 'PMP': -200, 'LDPE': -103}
    unif_values = {'Center': 17, 'Left': 10, 'Right': 0, 'Top': 6, 'Bottom': 6}
    mtf_values = {60: 1.19, 70: 1.19, 80: 1.09, 90: 0.75, 95: 0.58}
    line_lengths = {'Right-Vert': 49.9, 'Left-Vert': 50.05, 'Top-Horiz': 49.95, 'Bottom-Horiz': 49.96}

    @classmethod
    def setUpClass(cls):
        cls.cbct = CBCT.from_demo_images()


class VarianPelvis(CBCTMixin, unittest.TestCase):
    """Test the Varian Pelvis protocol CBCT."""
    location = osp.join(varian_test_file_dir, 'Pelvis.zip')
    expected_roll = 0.004
    slice_locations = {'HU': 32, 'UN': 3, 'SR': 44, 'LC': 20}
    hu_values = {'Poly': -36, 'Acrylic': 114, 'Delrin': 342, 'Air': -993, 'Teflon': 992, 'PMP': -188, 'LDPE': -95}
    unif_values = {'Center': 23, 'Left': 5, 'Right': 4, 'Top': 4, 'Bottom': 4}
    mtf_values = {60: 0.94, 70: 0.84, 80: 0.74, 90: 0.64, 95: 0.57}
    line_lengths = {'Right-Vert': 49.76, 'Left-Vert': 49.6, 'Top-Horiz': 50.02, 'Bottom-Horiz': 49.8}


class VarianLowDoseThorax(CBCTMixin, unittest.TestCase):
    """Test the Varian Low-Dose Thorax protocol CBCT."""
    location = osp.join(varian_test_file_dir, 'Low dose thorax.zip')
    expected_roll = 0.005
    slice_locations = {'HU': 32, 'UN': 3, 'SR': 44, 'LC': 20}
    hu_values = {'Poly': -46, 'Acrylic': 119, 'Delrin': 341, 'Air': -998, 'Teflon': 992, 'PMP': -193, 'LDPE': -97}
    unif_values = {'Center': 23, 'Left': 7, 'Right': -1, 'Top': 3, 'Bottom': 2}
    mtf_values = {60: 0.79, 70: 0.67, 80: 0.56, 90: 0.46, 95: 0.41}
    line_lengths = {'Right-Vert': 49.7, 'Left-Vert': 49.7, 'Top-Horiz': 49.7, 'Bottom-Horiz': 50.0}


class VarianPelvisSpotlight(CBCTMixin, unittest.TestCase):
    """Test the Varian Pelvis Spotligh protocol CBCT."""
    location = osp.join(varian_test_file_dir, 'Pelvis spotlight.zip')
    expected_roll = 0.004
    slice_locations = {'HU': 32, 'UN': 3, 'SR': 44, 'LC': 20}
    hu_values = {'Poly': -43, 'Acrylic': 118, 'Delrin': 341, 'Air': -998, 'Teflon': 967, 'PMP': -198, 'LDPE': -100}
    unif_values = {'Center': 19, 'Left': 3, 'Right': -1, 'Top': -1, 'Bottom': 0}
    mtf_values = {60: 1.19, 70: 1.19, 80: 1.02, 90: 0.71, 95: 0.58}
    line_lengths = {'Right-Vert': 50.08, 'Left-Vert': 49.91, 'Top-Horiz': 49.97, 'Bottom-Horiz': 49.92}


class VarianStandardHead(CBCTMixin, unittest.TestCase):
    """Test the Varian Standard Head protocol CBCT."""
    location = osp.join(varian_test_file_dir, 'Standard head.zip')
    expected_roll = 0.004
    slice_locations = {'HU': 31, 'UN': 2, 'SR': 43, 'LC': 19}
    hu_values = {'Poly': -43, 'Acrylic': 124, 'Delrin': 345, 'Air': -991, 'Teflon': 997, 'PMP': -199, 'LDPE': -101}
    unif_values = {'Center': 17, 'Left': 15, 'Right': 4, 'Top': 9, 'Bottom': 9}
    mtf_values = {60: 1.19, 70: 1.19, 80: 1.00, 90: 0.77, 95: 0.54}
    line_lengths = {'Right-Vert': 49.84, 'Left-Vert': 50.24, 'Top-Horiz': 49.91, 'Bottom-Horiz': 49.93}


class VarianLowDoseHead(CBCTMixin, unittest.TestCase):
    """Test the Varian Low-Dose Head protocol CBCT."""
    location = osp.join(varian_test_file_dir, 'Low dose head.zip')
    expected_roll = 0.007
    slice_locations = {'HU': 32, 'UN': 3, 'SR': 44, 'LC': 20}
    hu_values = {'Poly': -41, 'Acrylic': 123, 'Delrin': 350, 'Air': -990, 'Teflon': 998, 'PMP': -200, 'LDPE': -103}
    unif_values = {'Center': 16, 'Left': 11, 'Right': 3, 'Top': 7, 'Bottom': 6}
    mtf_values = {60: 1.19, 70: 1.08, 80: 0.89, 90: 0.66, 95: 0.53}
    line_lengths = {'Right-Vert': 49.63, 'Left-Vert': 50.25, 'Top-Horiz': 50.22, 'Bottom-Horiz': 49.90}


class GEMonthlyCT(CBCTMixin, unittest.TestCase):
    """Test a monthly CT scan from GE."""
    location = osp.join(other_test_file_dir, 'GE_CT.zip')
    expected_roll = 0
    hu_tolerance = 60
    slice_locations = {'HU': 143, 'UN': 85, 'SR': 167, 'LC': 119}
    hu_values = {'Poly': -32, 'Acrylic': 119, 'Delrin': 333, 'Air': -946, 'Teflon': 909, 'PMP': -173, 'LDPE': -87}
    unif_values = {'Center': 12, 'Left': 11, 'Right': 11, 'Top': 11, 'Bottom': 12}
    mtf_values = {60: 0.69, 70: 0.57, 80: 0.48, 90: 0.40, 95: 0.30}
    line_lengths = {'Right-Vert': 50.08, 'Left-Vert': 50.12, 'Top-Horiz': 49.8, 'Bottom-Horiz': 50.05}


class ToshibaMonthlyCT(CBCTMixin, unittest.TestCase):
    """Test a monthly CT scan from Toshiba."""
    location = osp.join(other_test_file_dir, 'Toshiba.zip')
    expected_roll = 0.002
    hu_tolerance = 240
    slice_locations = {'HU': 36, 'UN': 12, 'SR': 46, 'LC': 26}
    hu_values = {'Poly': -32, 'Acrylic': 106, 'Delrin': 467, 'Air': -994, 'Teflon': 1214, 'PMP': -165, 'LDPE': -85}
    unif_values = {'Center': 8, 'Left': 7, 'Right': 7, 'Top': 7, 'Bottom': 6}
    mtf_values = {60: 0.58, 70: 0.46, 80: 0.36, 90: 0.28, 95: 0.24}
    line_lengths = {'Right-Vert': 50.1, 'Left-Vert': 50.08, 'Top-Horiz': 50.12, 'Bottom-Horiz': 49.9}


class CBCT1(CBCTMixin, unittest.TestCase):
    """A Varian CBCT dataset"""
    location = osp.join(varian_test_file_dir, 'CBCT_1.zip')
    read_all = True
    expected_roll = 0.01
    slice_locations = {'HU': 31, 'UN': 2, 'SR': 43, 'LC': 19}
    hu_values = {'Poly': -39, 'Acrylic': 130, 'Delrin': 347, 'Air': -986, 'Teflon': 1002, 'PMP': -189, 'LDPE': -90}
    unif_values = {'Center': 13, 'Left': 17, 'Right': 5, 'Top': 13, 'Bottom': 11}
    mtf_values = {60: 1.19, 70: 1.19, 80: 0.91, 90: 0.41, 95: 0.30}
    line_lengths = {'Right-Vert': 49.5, 'Left-Vert': 49.93, 'Top-Horiz': 50.16, 'Bottom-Horiz': 50}

    def test_read_in_failure(self):
        """Test that importing the dataset in without the read_all parameter fails since files are poorly named."""
        with self.assertRaises(FileNotFoundError):
            CBCT.from_zip_file(self.location)


class CBCT2(CBCTMixin, unittest.TestCase):
    """A Varian CBCT dataset"""
    location = osp.join(varian_test_file_dir, 'CBCT_2.zip')
    read_all = True
    expected_roll = 0.004
    hu_tolerance = 50
    slice_locations = {'HU': 34, 'UN': 5, 'SR': 46, 'LC': 22}
    hu_values = {'Poly': -16, 'Acrylic': 135, 'Delrin': 367, 'Air': -967, 'Teflon': 1017, 'PMP': -163, 'LDPE': -71}
    unif_values = {'Center': 47, 'Left': 35, 'Right': 41, 'Top': 39, 'Bottom': 37}
    mtf_values = {60: 0.96, 70: 0.89, 80: 0.81, 90: 0.61, 95: 0.50}
    line_lengths = {'Right-Vert': 50.10, 'Left-Vert': 50.11, 'Top-Horiz': 49.93, 'Bottom-Horiz': 50.13}


class CBCT3(CBCTMixin, unittest.TestCase):
    """A Varian CBCT dataset"""
    location = osp.join(varian_test_file_dir, 'CBCT_3.zip')
    expected_roll = 0.046
    hu_tolerance = 50
    slice_locations = {'HU': 35, 'UN': 6, 'SR': 47, 'LC': 23}
    hu_values = {'Poly': -44, 'Acrylic': 110, 'Delrin': 325, 'Air': -979, 'Teflon': 949, 'PMP': -194, 'LDPE': -107}
    unif_values = {'Center': 2, 'Left': -1, 'Right': 11, 'Top': 9, 'Bottom': 2}
    mtf_values = {60: 1.19, 70: 1.12, 80: 0.97, 90: 0.82, 95: 0.46}
    line_lengths = {'Right-Vert': 49.93, 'Left-Vert': 49.57, 'Top-Horiz': 49.91, 'Bottom-Horiz': 50.18}


class CBCT4(CBCTMixin, unittest.TestCase):
    """A Varian CBCT dataset"""
    location = osp.join(varian_test_file_dir, 'CBCT_4.zip')
    expected_roll = 0.045
    slice_locations = {'HU': 31, 'UN': 2, 'SR': 43, 'LC': 19}
    hu_values = {'Poly': -33, 'Acrylic': 119, 'Delrin': 335, 'Air': -979, 'Teflon': 970, 'PMP': -185, 'LDPE': -94}
    unif_values = {'Center': 21, 'Left': 13, 'Right': 30, 'Top': 23, 'Bottom': 20}
    mtf_values = {60: 1.19, 70: 1.19, 80: 1.03, 90: 0.40, 95: 0.30}
    line_lengths = {'Right-Vert': 50.05, 'Left-Vert': 50.08, 'Top-Horiz': 50.2, 'Bottom-Horiz': 50.15}


class CBCT5(CBCTMixin, unittest.TestCase):
    """A Varian CBCT dataset"""
    location = osp.join(varian_test_file_dir, 'CBCT_5.zip')
    expected_roll = 0.0
    slice_locations = {'HU': 35, 'UN': 6, 'SR': 47, 'LC': 23}
    hu_values = {'Poly': -54, 'Acrylic': 101, 'Delrin': 328, 'Air': -999, 'Teflon': 975, 'PMP': -203, 'LDPE': -110}
    unif_values = {'Center': 19, 'Left': -8, 'Right': -5, 'Top': -7, 'Bottom': -6}
    mtf_values = {60: 0.86, 70: 0.75, 80: 0.61, 90: 0.50, 95: 0.45}
    line_lengths = {'Right-Vert': 49.99, 'Left-Vert': 50.10, 'Top-Horiz': 49.63, 'Bottom-Horiz': 49.48}


class CBCT6(CBCTMixin, unittest.TestCase):
    """A Varian CBCT dataset"""
    location = osp.join(varian_test_file_dir, 'CBCT_6.zip')
    expected_roll = -0.003
    slice_locations = {'HU': 38, 'UN': 9, 'SR': 50, 'LC': 26}
    hu_values = {'Poly': -42, 'Acrylic': 107, 'Delrin': 327, 'Air': -994, 'Teflon': 972, 'PMP': -192, 'LDPE': -100}
    unif_values = {'Center': -5, 'Left': 0, 'Right': -13, 'Top': -7, 'Bottom': -6}
    mtf_values = {60: 1.19, 70: 1.19, 80: 0.87, 90: 0.43, 95: 0.31}
    line_lengths = {'Right-Vert': 49.90, 'Left-Vert': 49.87, 'Top-Horiz': 49.03, 'Bottom-Horiz': 49.44}


class CBCT7(CBCTMixin, unittest.TestCase):
    """A Varian CBCT dataset"""
    location = osp.join(varian_test_file_dir, 'CBCT_7.zip')
    expected_roll = -0.008
    slice_locations = {'HU': 36, 'UN': 7, 'SR': 48, 'LC': 24}
    hu_values = {'Poly': -50, 'Acrylic': 108, 'Delrin': 334, 'Air': -999, 'Teflon': 982, 'PMP': -200, 'LDPE': -107}
    unif_values = {'Center': 14, 'Left': -5, 'Right': -5, 'Top': -5, 'Bottom': -5}
    mtf_values = {60: 0.85, 70: 0.76, 80: 0.67, 90: 0.56, 95: 0.48}
    line_lengths = {'Right-Vert': 49.92, 'Left-Vert': 49.04, 'Top-Horiz': 49.81, 'Bottom-Horiz': 50.17}


class CBCT12(CBCTMixin, unittest.TestCase):
    """A Varian CBCT dataset"""
    location = osp.join(varian_test_file_dir, 'CBCT_12.zip')
    expected_roll = -0.0015
    slice_locations = {'HU': 36, 'UN': 7, 'SR': 48, 'LC': 24}
    hu_values = {'Poly': -55, 'Acrylic': 112, 'Delrin': 335, 'Air': -999, 'Teflon': 982, 'PMP': -201, 'LDPE': -107}
    unif_values = {'Center': 5, 'Left': -5, 'Right': -9, 'Top': -7, 'Bottom': -6}
    mtf_values = {60: 0.67, 70: 0.56, 80: 0.48, 90: 0.39, 95: 0.30}
    line_lengths = {'Right-Vert': 49.67, 'Left-Vert': 49.57, 'Top-Horiz': 49.38, 'Bottom-Horiz': 49.73}


class CBCT13(CBCTMixin, unittest.TestCase):
    """A Varian CBCT dataset"""
    location = osp.join(varian_test_file_dir, 'CBCT_13.zip')
    expected_roll = -0.003
    slice_locations = {'HU': 35, 'UN': 6, 'SR': 47, 'LC': 23}
    hu_values = {'Poly': -53, 'Acrylic': 106, 'Delrin': 329, 'Air': -999, 'Teflon': 976, 'PMP': -200, 'LDPE': -107}
    unif_values = {'Center': 3, 'Left': -7, 'Right': -6, 'Top': -8, 'Bottom': -6}
    mtf_values = {60: 0.93, 70: 0.84, 80: 0.73, 90: 0.58, 95: 0.49}
    line_lengths = {'Right-Vert': 49.91, 'Left-Vert': 49.71, 'Top-Horiz': 49.96, 'Bottom-Horiz': 49.96}


class CBCT14(CBCTMixin, unittest.TestCase):
    """A Varian CBCT dataset"""
    location = osp.join(varian_test_file_dir, 'CBCT_14.zip')
    expected_roll = -0.014
    slice_locations = {'HU': 32, 'UN': 3, 'SR': 44, 'LC': 20}
    hu_values = {'Poly': -41, 'Acrylic': 125, 'Delrin': 334, 'Air': -995, 'Teflon': 986, 'PMP': -184, 'LDPE': -89}
    unif_values = {'Center': 18, 'Left': 13, 'Right': 15, 'Top': 14, 'Bottom': 14}
    mtf_values = {60: 0.78, 70: 0.61, 80: 0.54, 90: 0.46, 95: 0.43}
    line_lengths = {'Right-Vert': 49.92, 'Left-Vert': 49.32, 'Top-Horiz': 49.62, 'Bottom-Horiz': 50.07}


class CBCT15(CBCTMixin, unittest.TestCase):
    """A Varian CBCT dataset."""
    location = osp.join(varian_test_file_dir, 'CBCT_15.zip')
    expected_roll = 0
    slice_locations = {'HU': 60, 'UN': 23, 'SR': 75, 'LC': 45}
    hu_values = {'Poly': -32, 'Acrylic': 121, 'Delrin': 353, 'Air': -995, 'Teflon': 945, 'PMP': -186, 'LDPE': -93}
    unif_values = {'Center': -2, 'Left': 6, 'Right': 5, 'Top': 11, 'Bottom': 3}
    mtf_values = {60: 0.86, 70: 0.75, 80: 0.62, 90: 0.50, 95: 0.44}
    line_lengths = {'Right-Vert': 50.13, 'Left-Vert': 50.05, 'Top-Horiz': 49.95, 'Bottom-Horiz': 50}


class CBCT16(CBCTMixin, unittest.TestCase):
    """A Varian CBCT dataset"""
    location = osp.join(varian_test_file_dir, 'CBCT_16.zip')
    expected_roll = -0.003
    slice_locations = {'HU': 31, 'UN': 2, 'SR': 43, 'LC': 19}
    hu_values = {'Poly': -37, 'Acrylic': 128, 'Delrin': 342, 'Air': -995, 'Teflon': 1000, 'PMP': -181, 'LDPE': -87}
    unif_values = {'Center': 17, 'Left': 22, 'Right': 23, 'Top': 22, 'Bottom': 22}
    mtf_values = {60: 0.92, 70: 0.79, 80: 0.69, 90: 0.56, 95: 0.45}
    line_lengths = {'Right-Vert': 49.73, 'Left-Vert': 49.80, 'Top-Horiz': 49.65, 'Bottom-Horiz': 50}


class CBCT17(CBCTMixin, unittest.TestCase):
    """A Varian CBCT dataset"""
    location = osp.join(varian_test_file_dir, 'CBCT_17.zip')
    expected_roll = -0.008
    slice_locations = {'HU': 35, 'UN': 6, 'SR': 47, 'LC': 23}
    hu_values = {'Poly': -46, 'Acrylic': 117, 'Delrin': 344, 'Air': -989, 'Teflon': 989, 'PMP': -197, 'LDPE': -101}
    unif_values = {'Center': 5, 'Left': 0, 'Right': -7, 'Top': -6, 'Bottom': -2}
    mtf_values = {60: 1.19, 70: 1.19, 80: 1.19, 90: 0.74, 95: 0.64}
    line_lengths = {'Right-Vert': 49.65, 'Left-Vert': 50.12, 'Top-Horiz': 49.95, 'Bottom-Horiz': 50}