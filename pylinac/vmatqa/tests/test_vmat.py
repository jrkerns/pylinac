from __future__ import division, absolute_import, print_function

import unittest

import numpy as np

from pylinac.vmatqa.vmat import VMAT


class VMAT_general_tests(unittest.TestCase):
    """Generic tests for VMAT class."""
    def setUp(self):
        self.vmat = VMAT()

    def test_analyze_without_test_set(self):
        """raise an error if the test is not set when performing analysis"""
        self.vmat.load_demo_image()
        # trying to analyze without setting test should raise error
        self.assertRaises(AttributeError, self.vmat.analyze)

    def test_analyze_without_both_images_loaded(self):
        """Raise an error if both images aren't loaded when analyzing."""
        self.assertRaises(AttributeError, self.vmat.analyze)

class VMAT_demo_test_mlcs(unittest.TestCase):
    """Tests of the MLCS demo images."""
    def setUp(self):
        self.vmat = VMAT()

    def test_images_are_numpy(self):
        """Test that the demo images are numpy arrays."""
        self.vmat.load_demo_image('mlcs')
        self.assertIsInstance(self.vmat.image_open, np.ndarray)
        self.assertIsInstance(self.vmat.image_mlc, np.ndarray)

    def test_mlcs_passed(self):
        """Test that the demo images for MLC Speed test pass analysis."""
        self.vmat.load_demo_image('mlcs')
        self.vmat.analyze('mlcs', 3)
        self.assertTrue(self.vmat._passed_analysis)

    def test_segment_stdev(self):
        """Test that individual segment ratio standard deviations are close to what they should be."""
        self.vmat.load_demo_image('mlcs')
        self.vmat.analyze('mlcs')
        for seg_std in self.vmat._segment_stds:
            self.assertAlmostEqual(seg_std, 0.006, delta=0.003)

class VMAT_demo_test_drgs(unittest.TestCase):
    """Tests of the DRGS demo images."""
    def setUp(self):
        self.vmat = VMAT()

    def test_images_are_numpy(self):
        """Test that the demo images are numpy arrays."""
        self.vmat.load_demo_image('drgs')
        self.assertIsInstance(self.vmat.image_open, np.ndarray)
        self.assertIsInstance(self.vmat.image_mlc, np.ndarray)

    def test_passed(self):
        """Test that the demo images for Dose-Rate/Gantry Speed test pass analysis."""
        self.vmat.load_demo_image('drgs')
        self.vmat.analyze('drgs')
        self.assertTrue(self.vmat._passed_analysis)

    def test_sample_values(self):
        """Test that the individual sample ratios are close to what they should be."""
        self.vmat.load_demo_image('drgs')
        self.vmat.analyze('drgs')
        for sample in self.vmat._samples.ravel():
            self.assertAlmostEqual(sample, 1, delta=0.03)

    def test_segment_means(self):
        """Test the segment mean values."""

    def test_segment_stdev(self):
        """Test the segment standard deviations are close to what they should be."""
        self.vmat.load_demo_image('drgs')
        self.vmat.analyze('drgs')
        for seg_std in self.vmat._segment_stds:
            self.assertLessEqual(seg_std, 0.01)