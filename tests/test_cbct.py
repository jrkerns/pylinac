import unittest
import os
import os.path as osp
import time

import numpy as np

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

    def test_helpers(self):
        """Test the various helper methods."""
        self.cbct.load_demo_images()
        self.cbct.analyze()
        self.cbct._return_results()

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
        with self.assertRaises(ValueError):
            self.cbct.load_zip_file(mixed_zip)

    def test_phan_center(self):
        """Test locations of the phantom center."""
        self.cbct.load_demo_images()

        known_phan_center = Point(257, 255)
        self.cbct.analyze()
        self.assertAlmostEqual(self.cbct.hu.phan_center.x, known_phan_center.x, delta=0.7)
        self.assertAlmostEqual(self.cbct.hu.phan_center.y, known_phan_center.y, delta=0.7)

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
        for item in ['hu', 'un', 'mtf', 'sp', 'prof', 'lin', 'lc']:
            self.cbct.plot_analyzed_subimage(item)

        self.cbct.plot_analyzed_subimage('lin', delta=False)

        with self.assertRaises(ValueError):
            self.cbct.plot_analyzed_subimage('sr')


class CBCTMixin:
    """A mixin to use for Varian CBCT scans; does not inherit from TestCase as it would be run
        otherwise."""
    hu_tolerance = 40
    scaling_tolerance = 1
    zip = True
    expected_roll = 0
    slice_locations = {}
    hu_values = {}
    unif_values = {}
    mtf_values = {}
    avg_line_length = 50
    thickness_passed = True
    lowcon_visible = 0

    @classmethod
    def setUpClass(cls):
        if cls.zip:
            cls.cbct = CBCT.from_zip_file(cls.location)
        else:
            cls.cbct = CBCT(cls.location)
        cls.cbct.analyze(cls.hu_tolerance, cls.scaling_tolerance)

    def test_slice_thickness(self):
        """Test the slice thickness."""
        self.assertEqual(self.cbct.thickness.passed, self.thickness_passed)

    def test_lowcontrast_bubbles(self):
        """Test the number of low contrast bubbles visible."""
        self.assertEqual(self.cbct.lowcontrast.rois_visible, self.lowcon_visible)

    def test_all_passed(self):
        """Test the pass flags for all tests."""
        self.assertTrue(self.cbct.hu.overall_passed)
        self.assertTrue(self.cbct.uniformity.overall_passed)
        self.assertTrue(self.cbct.geometry.overall_passed)

    def test_slice_locations(self):
        """Test the locations of the slices of interest."""
        for attr, slice_name in zip(('hu_slice_num', 'un_slice_num', 'sr_slice_num', 'lc_slice_num'), ('HU', 'UN', 'SR', 'LC')):
            self.assertAlmostEqual(getattr(self.cbct.settings, attr), self.slice_locations[slice_name], delta=1)

    def test_phantom_roll(self):
        """Test the roll of the phantom."""
        self.assertAlmostEqual(self.cbct.settings.phantom_roll, self.expected_roll, delta=0.1)

    def test_HU_values(self):
        """Test HU values."""
        for key, roi in self.cbct.hu.rois.items():
            exp_val = self.hu_values[key]
            meas_val = roi.pixel_value
            self.assertAlmostEqual(exp_val, meas_val, delta=5)

    def test_uniformity_values(self):
        """Test Uniformity HU values."""
        for key, roi in self.cbct.uniformity.rois.items():
            exp_val = self.unif_values[key]
            meas_val = roi.pixel_value
            self.assertAlmostEqual(exp_val, meas_val, delta=5)

    def test_geometry_line_length(self):
        """Test the geometry distances."""
        self.assertAlmostEqual(self.avg_line_length, self.cbct.geometry.avg_line_length, delta=0.05)

    def test_MTF_values(self):
        """Test MTF values."""
        for key, exp_mtf in self.mtf_values.items():
            meas_mtf = self.cbct.spatialres.mtf(key)
            self.assertAlmostEqual(exp_mtf, meas_mtf, delta=0.1)


class CBCTDemo(CBCTMixin, unittest.TestCase):
    """Test the CBCT demo (Varian high quality head protocol)."""
    expected_roll = 0.3
    slice_locations = {'HU': 32, 'UN': 3, 'SR': 44, 'LC': 20}
    hu_values = {'Poly': -45, 'Acrylic': 117, 'Delrin': 341, 'Air': -998, 'Teflon': 997, 'PMP': -200, 'LDPE': -103}
    unif_values = {'Center': 17, 'Left': 10, 'Right': 0, 'Top': 6, 'Bottom': 6}
    mtf_values = {80: 0.76, 90: 0.61, 60: 0.99, 70: 0.88, 95: 0.45}
    avg_line_length = 49.92
    lowcon_visible = 4

    @classmethod
    def setUpClass(cls):
        cls.cbct = CBCT.from_demo_images()
        cls.cbct.analyze()


class CBCT4(CBCTMixin, unittest.TestCase):
    """A Varian CBCT dataset"""
    location = osp.join(varian_test_file_dir, 'CBCT_4.zip')
    expected_roll = 2.57
    slice_locations = {'HU': 31, 'UN': 2, 'SR': 43, 'LC': 19}
    hu_values = {'Poly': -33, 'Acrylic': 119, 'Delrin': 335, 'Air': -979, 'Teflon': 970, 'PMP': -185, 'LDPE': -94}
    unif_values = {'Center': 21, 'Left': 13, 'Right': 30, 'Top': 23, 'Bottom': 20}
    mtf_values = {80: 0.47, 90: 0.39, 60: 0.63, 70: 0.55, 95: 0.3}
    avg_line_length = 49.94
    thickness_passed = False
