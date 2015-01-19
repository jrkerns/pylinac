import unittest
import os.path as osp

import numpy as np

from pylinac.cbct import CBCT
from pylinac.core.geometry import Point


varian_test_file_dir = osp.join(osp.dirname(__file__), 'test_files', 'CBCT', 'Varian')


class Test_General(unittest.TestCase):
    """Test general things when using cbct module."""

    def setUp(self):
        self.cbct = CBCT()

    def test_demo(self):
        """Run the demo to make sure it works."""
        self.cbct.run_demo(show=False)

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

    def test_too_few_images_loaded(self):
        """Test when too few images are given to be loaded."""
        pass

    def test_phan_center(self):
        """Test locations of the phantom center."""
        self.cbct.load_demo_images()

        known_phan_center = Point(257, 255)
        self.cbct.construct_HU()
        self.assertAlmostEqual(self.cbct.HU.phan_center.x, known_phan_center.x, delta=0.7)
        self.assertAlmostEqual(self.cbct.HU.phan_center.y, known_phan_center.y, delta=0.7)

        # test a shifted image set
        shifted_phan_center = Point(287, 255)
        self.cbct.algo_data.images = np.roll(self.cbct.algo_data.images, 30, axis=1)
        self.cbct.construct_HU()
        self.assertAlmostEqual(self.cbct.HU.phan_center.x, shifted_phan_center.x, delta=0.7)
        self.assertAlmostEqual(self.cbct.HU.phan_center.y, shifted_phan_center.y, delta=0.7)


class Varian_CBCT:
    """A base class to use for Varian CBCT scans; does not inherit from TestCase as it would be run
        otherwise."""

    def test_all_passed(self):
        """Test the pass flags for all tests."""
        self.cbct.analyze()
        self.assertTrue(self.cbct.HU.overall_passed)
        self.assertTrue(self.cbct.UN.overall_passed)
        self.assertTrue(self.cbct.GEO.overall_passed)

    def test_slice_locations(self):
        """Test the locations of the slices of interest for a Varian scan."""
        self.assertEqual(self.cbct.algo_data.HU_slice_num, 32)
        self.assertEqual(self.cbct.algo_data.UN_slice_num, 9)
        self.assertEqual(self.cbct.algo_data.SR_slice_num, 44)
        self.assertEqual(self.cbct.algo_data.LC_slice_num, 23)

    def test_phantom_roll(self, expected_roll):
        """Test the roll of the demo phantom.

        Parameters
        ----------
        expected_roll : float
            The expected roll of the phantom in **degrees**.
        """
        self.assertAlmostEqual(self.cbct.algo_data.phantom_roll, expected_roll, delta=0.001)

    def test_HU_values(self, known_HU_dict):
        """Test HU values."""
        self.cbct.analyze()
        for key, roi in self.cbct.HU.ROIs.items():
            exp_val = known_HU_dict[key]
            meas_val = roi.pixel_value
            self.assertAlmostEqual(exp_val, meas_val, delta=3)

    def test_uniformity_values(self, known_HU_dict):
        """Test Uniformity HU values."""
        self.cbct.analyze()
        for key, roi in self.cbct.UN.ROIs.items():
            exp_val = known_HU_dict[key]
            meas_val = roi.pixel_value
            self.assertAlmostEqual(exp_val, meas_val, delta=3)

    def test_geometry_line_lengths(self, known_distances_dict):
        """Test the geometry distances."""
        self.cbct.analyze()
        for key, line in self.cbct.GEO.lines.items():
            exp_dist = known_distances_dict[key]
            meas_dist = line.length_mm(self.cbct.algo_data.mm_per_pixel)
            self.assertAlmostEqual(exp_dist, meas_dist, delta=0.05)

    def test_MTF_values(self, known_MTF_dict):
        """Test MTF values."""
        self.cbct.analyze()
        for key, exp_mtf in known_MTF_dict.items():
            meas_mtf = self.cbct.SR.get_MTF(key)
            self.assertAlmostEqual(exp_mtf, meas_mtf, delta=0.05)


class Test_CBCT_demo(unittest.TestCase, Varian_CBCT):
    """Test the CBCT demo (Varian high quality head protocol)."""
    @classmethod
    def setUpClass(cls):
        cls.cbct = CBCT()
        cls.cbct.load_demo_images(cleanup=False)

    @classmethod
    def tearDownClass(cls):
        """Remove the demo directory upon test class finish."""
        cls.cbct.load_demo_images(cleanup=True)

    def test_phantom_roll(self):
        exp_roll = 0.004
        super().test_phantom_roll(exp_roll)

    def test_HU_values(self):
        known_HU_dict = {'Poly': -45, 'Acrylic': 117, 'Delrin': 341, 'Air': -998, 'Teflon': 997, 'PMP': -200, 'LDPE': -103}
        super().test_HU_values(known_HU_dict)

    def test_uniformity_values(self):
        known_unifHUs = {'Center': 21, 'Left': 10, 'Right': 0, 'Top': 3, 'Bottom': 6}
        super().test_uniformity_values(known_unifHUs)

    def test_MTF_values(self):
        known_MTFs = {60: 1.19, 70: 1.19, 80: 1.09, 90: 0.75, 95: 0.58}
        super().test_MTF_values(known_MTFs)

    def test_geometry_line_lengths(self):
        known_dists = {'Right-Vert': 49.9, 'Left-Vert': 50.05, 'Top-Horiz': 49.95, 'Bottom-Horiz': 49.96}
        super().test_geometry_line_lengths(known_dists)

    def test_demo_cleanup(self):
        """Test that the demo images are/aren't cleaned up."""
        demo_dir = osp.join(osp.dirname(osp.dirname(__file__)), 'pylinac', 'demo_files', 'cbct')
        unzipped_demo_dir = osp.join(demo_dir, 'High quality head')

        # files should be left unzipped given flag
        self.cbct.load_demo_images(cleanup=False)
        self.assertTrue(osp.isdir(unzipped_demo_dir), 'Demo files not left unpacked')

        # files should automatically cleanup
        self.cbct.load_demo_images()
        self.assertFalse(osp.isdir(unzipped_demo_dir), 'Demo files not cleaned up after loading')


class Test_Varian_Pelvis(unittest.TestCase, Varian_CBCT):
    """Test the Varian Pelvis protocol CBCT."""
    pelvis_zip = osp.join(varian_test_file_dir, 'Pelvis.zip')
    pelvis_dir = osp.join(varian_test_file_dir, 'Pelvis')

    @classmethod
    def setUpClass(cls):
        cls.cbct = CBCT()
        cls.cbct.load_zip_file(cls.pelvis_zip)

    def test_phantom_roll(self):
        exp_roll = 0.004
        super().test_phantom_roll(exp_roll)

    def test_HU_values(self):
        known_HU_dict = {'Poly': -36, 'Acrylic': 114, 'Delrin': 342, 'Air': -993, 'Teflon': 992, 'PMP': -188,
                         'LDPE': -95}
        super().test_HU_values(known_HU_dict)

    def test_uniformity_values(self):
        known_unifHUs = {'Center': 14, 'Left': 5, 'Right': 4, 'Top': 4, 'Bottom': 4}
        super().test_uniformity_values(known_unifHUs)

    def test_MTF_values(self):
        known_MTFs = {60: 0.94, 70: 0.84, 80: 0.74, 90: 0.64, 95: 0.57}
        super().test_MTF_values(known_MTFs)

    def test_geometry_line_lengths(self):
        known_dists = {'Right-Vert': 49.76, 'Left-Vert': 49.6, 'Top-Horiz': 50.02, 'Bottom-Horiz': 49.8}
        super().test_geometry_line_lengths(known_dists)


class Test_Varian_Low_Dose_Thorax(unittest.TestCase, Varian_CBCT):
    """Test the Varian Low-Dose Thorax protocol CBCT."""
    thorax_zip = osp.join(varian_test_file_dir, 'Low dose thorax.zip')
    thorax_dir, _ = osp.splitext(thorax_zip)

    @classmethod
    def setUpClass(cls):
        cls.cbct = CBCT()
        cls.cbct.load_zip_file(cls.thorax_zip)

    def test_phantom_roll(self):
        exp_roll = 0.005
        super().test_phantom_roll(exp_roll)

    def test_HU_values(self):
        known_HU_dict = {'Poly': -46, 'Acrylic': 119, 'Delrin': 341, 'Air': -998, 'Teflon': 992, 'PMP': -193,
                         'LDPE': -97}
        super().test_HU_values(known_HU_dict)

    def test_uniformity_values(self):
        known_unifHUs = {'Center': 14, 'Left': 7, 'Right': -1, 'Top': 3, 'Bottom': 2}
        super().test_uniformity_values(known_unifHUs)

    def test_MTF_values(self):
        known_MTFs = {60: 0.79, 70: 0.67, 80: 0.56, 90: 0.46, 95: 0.41}
        super().test_MTF_values(known_MTFs)

    def test_geometry_line_lengths(self):
        known_dists = {'Right-Vert': 49.7, 'Left-Vert': 49.7, 'Top-Horiz': 49.7, 'Bottom-Horiz': 50.0}
        super().test_geometry_line_lengths(known_dists)
