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


class Test_General(unittest.TestCase):
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
    """A base class to use for Varian CBCT scans; does not inherit from TestCase as it would be run
        otherwise."""

    @classmethod
    def setUpClass(cls):
        if cls.zip:
            cls.cbct = CBCT.from_zip_file(cls.location)
        else:
            cls.cbct = CBCT.from_folder(cls.location)

    def test_all_passed(self, hu_tolerance=40, scaling_tolerance=1):
        """Test the pass flags for all tests."""
        self.cbct.analyze(hu_tolerance, scaling_tolerance)
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


class Test_CBCT_demo(CBCTMixin, unittest.TestCase):
    """Test the CBCT demo (Varian high quality head protocol)."""
    expected_roll = 0.004
    slice_locations = {'HU': 32, 'UN': 2, 'SR': 44, 'LC': 20}
    hu_values = {'Poly': -45, 'Acrylic': 117, 'Delrin': 341, 'Air': -998, 'Teflon': 997, 'PMP': -200, 'LDPE': -103}
    unif_values = {'Center': 17, 'Left': 10, 'Right': 0, 'Top': 6, 'Bottom': 6}
    mtf_values = {60: 1.19, 70: 1.19, 80: 1.09, 90: 0.75, 95: 0.58}
    line_lengths = {'Right-Vert': 49.9, 'Left-Vert': 50.05, 'Top-Horiz': 49.95, 'Bottom-Horiz': 49.96}

    @classmethod
    def setUpClass(cls):
        cls.cbct = CBCT.from_demo_images()


class Test_Varian_Pelvis(CBCTMixin, unittest.TestCase):
    """Test the Varian Pelvis protocol CBCT."""
    location = osp.join(varian_test_file_dir, 'Pelvis.zip')
    zip = True
    expected_roll = 0.004
    slice_locations = {'HU': 32, 'UN': 2, 'SR': 44, 'LC': 20}
    hu_values = {'Poly': -36, 'Acrylic': 114, 'Delrin': 342, 'Air': -993, 'Teflon': 992, 'PMP': -188, 'LDPE': -95}
    unif_values = {'Center': 23, 'Left': 5, 'Right': 4, 'Top': 4, 'Bottom': 4}
    mtf_values = {60: 0.94, 70: 0.84, 80: 0.74, 90: 0.64, 95: 0.57}
    line_lengths = {'Right-Vert': 49.76, 'Left-Vert': 49.6, 'Top-Horiz': 50.02, 'Bottom-Horiz': 49.8}


class Test_Varian_Low_Dose_Thorax(CBCTMixin, unittest.TestCase):
    """Test the Varian Low-Dose Thorax protocol CBCT."""
    location = osp.join(varian_test_file_dir, 'Low dose thorax.zip')
    zip = True
    expected_roll = 0.005
    slice_locations = {'HU': 32, 'UN': 2, 'SR': 44, 'LC': 20}
    hu_values = {'Poly': -46, 'Acrylic': 119, 'Delrin': 341, 'Air': -998, 'Teflon': 992, 'PMP': -193, 'LDPE': -97}
    unif_values = {'Center': 23, 'Left': 7, 'Right': -1, 'Top': 3, 'Bottom': 2}
    mtf_values = {60: 0.79, 70: 0.67, 80: 0.56, 90: 0.46, 95: 0.41}
    line_lengths = {'Right-Vert': 49.7, 'Left-Vert': 49.7, 'Top-Horiz': 49.7, 'Bottom-Horiz': 50.0}


class TestGEMonthlyCT(CBCTMixin, unittest.TestCase):
    """Test a monthly CT scan from GE."""
    location = osp.join(other_test_file_dir, 'GE.zip')
    zip = True
    expected_roll = 0
    slice_locations = {'HU': 143, 'UN': 83, 'SR': 167, 'LC': 119}
    hu_values = {'Poly': -32, 'Acrylic': 119, 'Delrin': 333, 'Air': -946, 'Teflon': 909, 'PMP': -173, 'LDPE': -87}
    unif_values = {'Center': 12, 'Left': 11, 'Right': 11, 'Top': 11, 'Bottom': 12}
    mtf_values = {60: 0.69, 70: 0.57, 80: 0.48, 90: 0.40, 95: 0.30}
    line_lengths = {'Right-Vert': 50.08, 'Left-Vert': 50.12, 'Top-Horiz': 49.8, 'Bottom-Horiz': 50.05}

    def test_all_passed(self, hu_tolerance=60):
        super().test_all_passed(hu_tolerance)


class TestToshibaMonthlyCT(CBCTMixin, unittest.TestCase):
    """Test a monthly CT scan from Toshiba."""
    location = osp.join(other_test_file_dir, 'Toshiba.zip')
    zip = True
    expected_roll = 0.002
    slice_locations = {'HU': 36, 'UN': 11, 'SR': 46, 'LC': 26}
    hu_values = {'Poly': -32, 'Acrylic': 106, 'Delrin': 467, 'Air': -994, 'Teflon': 1214, 'PMP': -165, 'LDPE': -85}
    unif_values = {'Center': 8, 'Left': 7, 'Right': 7, 'Top': 7, 'Bottom': 6}
    mtf_values = {60: 0.58, 70: 0.46, 80: 0.36, 90: 0.28, 95: 0.24}
    line_lengths = {'Right-Vert': 50.1, 'Left-Vert': 50.08, 'Top-Horiz': 50.12, 'Bottom-Horiz': 49.9}

    def test_all_passed(self, hu_tolerance=240):
        super().test_all_passed(hu_tolerance)