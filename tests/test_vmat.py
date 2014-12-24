from __future__ import division, absolute_import, print_function

import unittest
import os.path as osp

from pylinac.vmat import VMAT


test_files_dir = osp.join(osp.dirname(__file__), 'test_files')
drgs_105_open = osp.join(test_files_dir, 'DRGSopen-105-example.dcm')
drgs_105_dmlc = osp.join(test_files_dir, 'DRGSmlc-105-example.dcm')
drmlc_105_open = osp.join(test_files_dir, 'DRMLCopen-105-example.dcm')
drmlc_105_dmlc = osp.join(test_files_dir, 'DRMLCmlc-105-example.dcm')


class Test_general(unittest.TestCase):
    """Generic tests for VMAT class."""
    def setUp(self):
        self.vmat = VMAT()

    def test_analyze_without_both_images_loaded(self):
        """Raise an error if both images aren't loaded when analyzing."""
        self.assertRaises(AttributeError, self.vmat.analyze, 'drmlc')

    def test_throw_tolerance_error(self):
        """Raise an error if tolerance value passed to analyze is bad."""
        self.vmat.load_demo_image()
        # below allowance
        self.assertRaises(ValueError, self.vmat.analyze, 'drgs', 0)
        # above allowance
        self.assertRaises(ValueError, self.vmat.analyze, 'drgs', 9)
        # string instead of number
        self.assertRaises(TypeError, self.vmat.analyze, 'drgs', '1')
        # valid tolerance; shouldn't throw error
        self.vmat.analyze('drgs', 4)

    def test_throw_test_type_error(self):
        """Raise an error if the test string passed is bad."""
        self.vmat.load_demo_image()
        self.assertRaises(ValueError, self.vmat.analyze, 'drsg', 3)


class Test_DRGS_demo(unittest.TestCase):
    """Tests of the result values of the DRMLC demo images."""
    def setUp(self):
        self.vmat = VMAT()
        self.vmat.load_demo_image('drgs')

    def test_overall_passed(self):
        """Test that the overall pass flag is true for default settings"""
        self.vmat.analyze('drgs')
        self.assertTrue(self.vmat.passed)

    def test_img_inversion(self):
        """Check that the demo images indeed get inverted."""
        top_corner_before = self.vmat.image_open.pixel_array[:20, :20].mean()
        self.vmat._check_img_inversion()
        top_corner_after = self.vmat.image_open.pixel_array[:20, :20].mean()
        self.assertNotEqual(top_corner_before, top_corner_after)

    def test_scaling(self):
        """Test scaling values are 1; i.e. 384x512@150cmSID."""
        SID_scale, img_scaling = self.vmat._calc_im_scaling_factors()
        self.assertEqual(SID_scale, 1)
        self.assertEqual(img_scaling.x, 1)
        self.assertEqual(img_scaling.y, 1)

    def test_segments(self):
        """Test various values of the Segments."""
        self.vmat.analyze('drgs')

        # test center locations
        # ...segment 0
        self.assertEqual(self.vmat.segments[0].center.x, 159)
        self.assertEqual(self.vmat.segments[0].center.y, 192)
        # ...segment 4
        self.assertEqual(self.vmat.segments[4].center.x, 311)
        self.assertEqual(self.vmat.segments[4].center.y, 192)

        # test segment properties
        self.assertAlmostEqual(self.vmat.segments[1].abs_mean_dev, 0.0315, delta=0.003)
        self.assertAlmostEqual(self.vmat.segments[2].max_dev, -0.283, delta=0.003)
        self.assertAlmostEqual(self.vmat.segments[5].mean_ratio, 1.003, delta=0.0005)
        self.assertAlmostEqual(self.vmat.segments[6].min_dev, 0.420, delta=0.003)
        self.assertAlmostEqual(self.vmat.segments[0].deviations.max(), 1.234, delta=0.002)

    def test_samples(self):
        """Test various property values of Samples."""
        self.vmat.analyze('drgs')

        # correct number (38 for Millenium MLC; 76 for HDMLC)
        self.assertEqual(len(self.vmat.num_samples), 38)

        # test some center locations
        self.assertEqual(self.vmat.segments[0].samples[0].center.x, 159)
        self.assertEqual(self.vmat.segments[0].samples[0].center.y, 15)

        self.assertEqual(self.vmat.segments[2].samples[8].center.x, 235)
        self.assertEqual(self.vmat.segments[2].samples[8].center.y, 91)

        self.assertEqual(self.vmat.segments[6].samples[37].center.x, 390)
        self.assertEqual(self.vmat.segments[6].samples[37].center.y, 368)

        # test some sample ratios
        self.assertAlmostEqual(self.vmat.segments[6].samples[29].ratio, 1.009, delta=0.0005)
        self.assertAlmostEqual(self.vmat.segments[1].samples[3].ratio, 1.007, delta=0.0005)
        self.assertAlmostEqual(self.vmat.segments[4].samples[10].ratio, 0.9985, delta=0.0005)
        self.assertAlmostEqual(self.vmat.segments[5].samples[32].ratio, 1.0016, delta=0.0005)


class Test_DRMLC_demo(unittest.TestCase):
    """Tests of the result values of the DRMLC demo images."""

    def setUp(self):
        self.vmat = VMAT()
        self.vmat.load_demo_image('drmlc')

    def test_overall_passed(self):
        """Test that the overall pass flag is true for default settings"""
        self.vmat.analyze('drmlc')
        self.assertTrue(self.vmat.passed)

    def test_img_inversion(self):
        """Check that the demo images indeed get inverted."""
        top_corner_before = self.vmat.image_open.pixel_array[:20, :20].mean()
        self.vmat._check_img_inversion()
        top_corner_after = self.vmat.image_open.pixel_array[:20, :20].mean()
        self.assertNotEqual(top_corner_before, top_corner_after)

    def test_scaling(self):
        """Test scaling values are 1; i.e. 384x512@150cmSID."""
        SID_scale, img_scaling = self.vmat._calc_im_scaling_factors()
        self.assertEqual(SID_scale, 1)
        self.assertEqual(img_scaling.x, 1)
        self.assertEqual(img_scaling.y, 1)

    def test_segments(self):
        """Test various values of the Segments."""
        self.vmat.analyze('drmlc')

        # test center locations
        # ...segment 0
        self.assertEqual(self.vmat.segments[0].center.x, 170)
        self.assertEqual(self.vmat.segments[0].center.y, 192)
        # ...segment 2
        self.assertEqual(self.vmat.segments[2].center.x, 288)
        self.assertEqual(self.vmat.segments[2].center.y, 192)

        # test segment properties
        self.assertAlmostEqual(self.vmat.segments[0].abs_mean_dev, 0.406, delta=0.003)
        self.assertAlmostEqual(self.vmat.segments[1].max_dev, -0.273, delta=0.003)
        self.assertAlmostEqual(self.vmat.segments[3].mean_ratio, 1.008, delta=0.0005)
        self.assertAlmostEqual(self.vmat.segments[2].min_dev, -0.500, delta=0.003)
        self.assertAlmostEqual(self.vmat.segments[0].deviations.max(), 0.497, delta=0.002)

    def test_samples(self):
        """Test various property values of Samples."""
        self.vmat.analyze('drmlc')

        # correct number (38 for Millennium MLC; 76 for HDMLC)
        self.assertEqual(self.vmat.num_samples, 38)

        # test some center locations
        self.assertEqual(self.vmat.segments[0].samples[0].center.x, 170)
        self.assertEqual(self.vmat.segments[0].samples[0].center.y, 15)

        self.assertEqual(self.vmat.segments[2].samples[8].center.x, 288)
        self.assertEqual(self.vmat.segments[2].samples[8].center.y, 91)

        self.assertEqual(self.vmat.segments[3].samples[37].center.x, 345)
        self.assertEqual(self.vmat.segments[3].samples[37].center.y, 368)

        # test some sample ratios
        self.assertAlmostEqual(self.vmat.segments[0].samples[29].ratio, 1.006, delta=0.0005)
        self.assertAlmostEqual(self.vmat.segments[1].samples[3].ratio, 1.003, delta=0.0005)
        self.assertAlmostEqual(self.vmat.segments[2].samples[10].ratio, 0.999, delta=0.0005)
        self.assertAlmostEqual(self.vmat.segments[3].samples[32].ratio, 1.006, delta=0.0005)


class Test_DRMLC_105(unittest.TestCase):
    """Tests of the result values of DRMLC images at 105cm SID."""

    def setUp(self):
        self.vmat = VMAT()
        self.vmat.load_image(drmlc_105_open, 'open')
        self.vmat.load_image(drmlc_105_dmlc, 'dmlc')

    def test_overall_passed(self):
        """Test that the overall pass flag is true for default settings"""
        self.vmat.analyze('drmlc')
        self.assertTrue(self.vmat.passed)

    def test_img_inversion(self):
        """Check that the images indeed get inverted (applicable to most EPID images)."""
        top_corner_before = self.vmat.image_open.pixel_array[:20, :20].mean()
        self.vmat._check_img_inversion()
        top_corner_after = self.vmat.image_open.pixel_array[:20, :20].mean()
        self.assertNotEqual(top_corner_before, top_corner_after)

    def test_scaling(self):
        """Test scaling values."""
        SID_scale, img_scaling = self.vmat._calc_im_scaling_factors()
        # SID scale is 0.7; i.e. 150 (reference SID) / 105 (measured SID).
        self.assertEqual(SID_scale, 0.7)
        # image scales are 2.0; 768x1024 (measurement on AS1000) vs 384x512 (reference on AS500)
        self.assertEqual(img_scaling.x, 2)
        self.assertEqual(img_scaling.y, 2)

    def test_segments(self):
        """Test various values of the Segments."""
        self.vmat.analyze('drmlc')

        # test center locations
        # ...segment 0
        self.assertEqual(self.vmat.segments[0].center.x, 426)
        self.assertEqual(self.vmat.segments[0].center.y, 384)
        # ...segment 2
        self.assertEqual(self.vmat.segments[2].center.x, 544)
        self.assertEqual(self.vmat.segments[2].center.y, 384)

        # test segment properties
        self.assertAlmostEqual(self.vmat.segments[0].abs_mean_dev, 0.098, delta=0.002)
        self.assertAlmostEqual(self.vmat.segments[1].max_dev, 0.204, delta=0.003)
        self.assertAlmostEqual(self.vmat.segments[3].mean_ratio, 1.001, delta=0.001)
        self.assertAlmostEqual(self.vmat.segments[2].min_dev, -0.283, delta=0.002)
        self.assertAlmostEqual(self.vmat.segments[0].deviations.max(), 0.257, delta=0.002)

    def test_samples(self):
        """Test various property values of Samples."""
        self.vmat.analyze('drmlc')

        # correct number (38 for Millennium MLC; 76 for HDMLC)
        self.assertEqual(self.vmat.num_samples, 38)

        # test some center locations
        self.assertEqual(self.vmat.segments[0].samples[0].center.x, 426)
        self.assertEqual(self.vmat.segments[0].samples[0].center.y, 136)

        self.assertEqual(self.vmat.segments[2].samples[8].center.x, 544)
        self.assertEqual(self.vmat.segments[2].samples[8].center.y, 243)

        self.assertEqual(self.vmat.segments[3].samples[37].center.x, 601)
        self.assertEqual(self.vmat.segments[3].samples[37].center.y, 631)

        # test some sample ratios
        self.assertAlmostEqual(self.vmat.segments[0].samples[29].ratio, 0.995, delta=0.0005)
        self.assertAlmostEqual(self.vmat.segments[1].samples[3].ratio, 1.000, delta=0.0005)
        self.assertAlmostEqual(self.vmat.segments[2].samples[10].ratio, 0.996, delta=0.0005)
        self.assertAlmostEqual(self.vmat.segments[3].samples[32].ratio, 0.995, delta=0.0005)


class Test_DRGS_105(unittest.TestCase):
    """Tests of the result values of DRMLC images at 105cm SID."""

    def setUp(self):
        self.vmat = VMAT()
        self.vmat.load_image(drgs_105_open, 'open')
        self.vmat.load_image(drgs_105_dmlc, 'dmlc')

    def test_overall_passed(self):
        """Test that the overall pass flag is true for default settings"""
        self.vmat.analyze('drgs')
        self.assertTrue(self.vmat.passed)

    def test_img_inversion(self):
        """Check that the images indeed get inverted (applicable to most EPID images)."""
        top_corner_before = self.vmat.image_open.pixel_array[:20, :20].mean()
        self.vmat._check_img_inversion()
        top_corner_after = self.vmat.image_open.pixel_array[:20, :20].mean()
        self.assertNotEqual(top_corner_before, top_corner_after)

    def test_scaling(self):
        """Test scaling values."""
        SID_scale, img_scaling = self.vmat._calc_im_scaling_factors()
        # SID scale is 0.7; i.e. 150 (reference SID) / 105 (measured SID).
        self.assertEqual(SID_scale, 0.7)
        # image scales are 2.0; 768x1024 (measurement on AS1000) vs 384x512 (reference on AS500)
        self.assertEqual(img_scaling.x, 2)
        self.assertEqual(img_scaling.y, 2)

    def test_segments(self):
        """Test various values of the Segments."""
        self.vmat.analyze('drgs')

        # test center locations
        self.assertEqual(self.vmat.segments[0].center.x, 415)
        self.assertEqual(self.vmat.segments[0].center.y, 384)

        self.assertEqual(self.vmat.segments[2].center.x, 491)
        self.assertEqual(self.vmat.segments[2].center.y, 384)

        # test segment properties
        self.assertAlmostEqual(self.vmat.segments[0].abs_mean_dev, 0.290, delta=0.002)
        self.assertAlmostEqual(self.vmat.segments[1].max_dev, 0.389, delta=0.003)
        self.assertAlmostEqual(self.vmat.segments[5].mean_ratio, 0.995, delta=0.001)
        self.assertAlmostEqual(self.vmat.segments[2].min_dev, -0.018, delta=0.002)
        self.assertAlmostEqual(self.vmat.segments[0].deviations.max(), 0.584, delta=0.002)

    def test_samples(self):
        """Test various property values of Samples."""
        self.vmat.analyze('drgs')

        # correct number (38 for Millennium MLC; 76 for HDMLC)
        self.assertEqual(self.vmat.num_samples, 38)

        # test some center locations
        self.assertEqual(self.vmat.segments[0].samples[0].center.x, 415)
        self.assertEqual(self.vmat.segments[0].samples[0].center.y, 136)

        self.assertEqual(self.vmat.segments[2].samples[8].center.x, 491)
        self.assertEqual(self.vmat.segments[2].samples[8].center.y, 243)

        self.assertEqual(self.vmat.segments[6].samples[37].center.x, 646)
        self.assertEqual(self.vmat.segments[6].samples[37].center.y, 631)

        # test some sample ratios
        self.assertAlmostEqual(self.vmat.segments[0].samples[29].ratio, 0.996, delta=0.0005)
        self.assertAlmostEqual(self.vmat.segments[1].samples[6].ratio, 0.998, delta=0.0005)
        self.assertAlmostEqual(self.vmat.segments[5].samples[23].ratio, 0.999, delta=0.0005)
        self.assertAlmostEqual(self.vmat.segments[6].samples[37].ratio, 0.994, delta=0.0005)