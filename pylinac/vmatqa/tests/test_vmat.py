from __future__ import division, absolute_import, print_function

import unittest

import numpy as np

from pylinac.vmatqa.vmat import VMAT


class VMAT_general_tests(unittest.TestCase):
    """Generic tests for VMAT class."""
    def setUp(self):
        self.vmat = VMAT()

    def test_analyze_without_both_images_loaded(self):
        """Raise an error if both images aren't loaded when analyzing."""
        self.assertRaises(AttributeError, self.vmat.analyze, 'drmlc')

    def test_throw_tolerance_error(self):
        """Raise an error if tolerance value passed to analyze is bad."""
        self.vmat.load_demo_image()
        self.assertRaises(ValueError, self.vmat.analyze, 'drgs', 0)  # below allowance
        self.assertRaises(ValueError, self.vmat.analyze, 'drgs', 9)  # above allowance
        self.assertRaises(TypeError, self.vmat.analyze, 'drgs', '1')  # string instead of number

    def test_throw_test_error(self):
        """Raise an error if the test string passed is bad."""
        self.vmat.load_demo_image()
        # self.assertRaises(ValueError, self.vmat.analyze, 'dgsr', 3)  # User forgets to input test type
        self.assertRaises(ValueError, self.vmat.analyze, 'drsg', 3)  # misspelled test name

class VMAT_demo_test_drmlc(unittest.TestCase):
    """Tests of the DRMLC demo images."""
    def setUp(self):
        self.vmat = VMAT()
        self.vmat.load_demo_image('drmlc')

    def test_images_are_numpy(self):
        """Test that the demo images are numpy arrays."""
        self.assertIsInstance(self.vmat.image_open, np.ndarray)
        self.assertIsInstance(self.vmat.image_dmlc, np.ndarray)

    def test_drmlc_passed(self):
        """Test that the demo images for MLC Speed test pass analysis."""
        self.vmat.analyze('drmlc', 3)
        self.assertTrue(self.vmat._passed_analysis)


class VMAT_demo_test_drgs(unittest.TestCase):
    """Tests of the DRGS demo images."""
    def setUp(self):
        self.vmat = VMAT()
        self.vmat.load_demo_image('drgs')

    def test_images_are_numpy(self):
        """Test that the demo images are numpy arrays."""
        self.assertIsInstance(self.vmat.image_open, np.ndarray)
        self.assertIsInstance(self.vmat.image_dmlc, np.ndarray)

    def test_passed(self):
        """Test that the demo images for Dose-Rate/Gantry Speed test pass analysis."""
        self.vmat.analyze('drgs')
        self.assertTrue(self.vmat._passed_analysis)

    def test_overall_results(self):
        """Test that Overall Results are close to what they should be."""
        self.vmat.analyze('drgs')
        self.assertAlmostEqual(self.vmat._dev_max, 0.9, delta=0.1)
        self.assertAlmostEqual(self.vmat._dev_mean, 0.6, delta=0.1)
        self.assertAlmostEqual(self.vmat._dev_min, -1.6, delta=0.1)

    def test_segment_results(self):
        """Test that the segment results are close to what they should be."""
        self.vmat.analyze('drgs')
        seg_max = (-1.1, 0.06, 0.66, 0.9, 0.8, 0.3, -0.6)
        seg_min = (-1.6, 0.01, 0.4, 0.55, 0.5, 0.2, -1.0)
        seg_mean = (1.4, 0.05, 0.54, 0.7, 0.65, 0.26, 0.83)
        for seg in range(self.vmat._num_segments):
            self.assertAlmostEqual(self.vmat._seg_dev_max[seg], seg_max[seg], delta=0.1)
            self.assertAlmostEqual(self.vmat._seg_dev_min[seg], seg_min[seg], delta=0.1)
            self.assertAlmostEqual(self.vmat._seg_dev_mean[seg], seg_mean[seg], delta=0.1)
