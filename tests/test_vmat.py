from functools import partial
import unittest
import os.path as osp

from pylinac.vmat import VMAT
from pylinac.core.geometry import Point


_vmat_test_files_dir = osp.join(osp.dirname(__file__), 'test_files', 'VMAT')

within_1 = partial(unittest.TestCase().assertAlmostEqual, delta=1)
within_01 = partial(unittest.TestCase().assertAlmostEqual, delta=0.1)

class Test_general(unittest.TestCase):
    """Generic tests for VMAT class."""
    def setUp(self):
        self.vmat = VMAT()

    def test_demo_image_loads(self):
        """Test the that demo images load properly."""
        # shouldn't raise
        self.vmat.load_demo_image('drgs')
        self.vmat.load_demo_image('mlcs')

    def test_analyze_without_both_images_loaded(self):
        """Raise an error if both images aren't loaded when analyzing."""
        self.assertRaises(AttributeError, self.vmat.analyze, 'mlcs')
        self.vmat.load_demo_image('drgs')
        self.vmat.analyze('drgs')  # shouldn't raise

    def test_img_loaded_tags(self):
        """Test the 'is_loaded' type tags."""
        self.assertFalse(self.vmat.open_img_is_loaded)
        self.assertFalse(self.vmat.dmlc_img_is_loaded)
        self.vmat.load_demo_image()
        self.assertTrue(self.vmat.open_img_is_loaded)
        self.assertTrue(self.vmat.dmlc_img_is_loaded)

    def test_number_of_segments(self):
        """Test that the right amount of segments are constructed based on the given test."""
        self.vmat.load_demo_image()
        self.vmat.analyze('drgs')
        self.assertEqual(len(self.vmat.segments), 7)
        self.vmat.analyze('mlcs')
        self.assertEqual(len(self.vmat.segments), 4)


class VMATMixin:

    def test_overall_passed(self, test):
        self.vmat.analyze(test)
        self.assertTrue(self.vmat.passed)

    def test_img_inversion(self):
        """Check that the demo images indeed get inverted."""
        top_corner_before = self.vmat.image_open.array[:20, :20].mean()
        self.vmat._check_img_inversion()
        top_corner_after = self.vmat.image_open.array[:20, :20].mean()
        self.assertNotEqual(top_corner_before, top_corner_after)

    def test_segment_positions(self, segment_dict):
        for key, value in segment_dict.items():
            within_1(self.vmat.segments[key].center.x, value.x)
            within_1(self.vmat.segments[key].center.y, value.y)

    def test_segment_values(self, segment_dict):
        for key, value in segment_dict.items():
            within_01(self.vmat.segments[key].r_dev, value['r_dev'])
            within_01(self.vmat.segments[key].r_corr, value['r_corr'])


class Test_DRGS_demo(VMATMixin, unittest.TestCase):
    """Tests of the result values of the DRGS demo images."""
    def setUp(self):
        self.vmat = VMAT.from_demo_images('drgs')
        self.vmat.settings.x_offset = 20

    def test_demo(self):
        """Run the demo; no errors should arise."""
        self.vmat.run_demo_drgs()

    def test_overall_passed(self):
        """Test that the overall pass flag is true for default settings"""
        super().test_overall_passed('drgs')

    def test_segment_positions(self):
        """Test various values of the Segments."""
        self.vmat.analyze('drgs')
        segment_dict = {0: Point(161, 192), 4: Point(314, 192)}
        super().test_segment_positions(segment_dict)

    def test_segment_values(self):
        self.vmat.analyze('drgs')
        segment_dict = {
            0: {'r_dev': 0.965, 'r_corr': 101.85},
            4: {'r_dev': -0.459, 'r_corr': 100.42},
        }
        super().test_segment_values(segment_dict)

class Test_MLCS_demo(VMATMixin, unittest.TestCase):
    """Tests of the result values of the DRMLC demo images."""

    def setUp(self):
        self.vmat = VMAT.from_demo_images('mlcs')

    def test_demo(self):
        self.vmat.run_demo_mlcs()

    def test_overall_passed(self):
        """Test that the overall pass flag is true for default settings"""
        super().test_overall_passed('mlcs')

    def test_segment_positions(self):
        """Test various values of the Segments."""
        self.vmat.analyze('mlcs')
        segment_dict = {0: Point(170, 192), 2: Point(285, 192)}
        super().test_segment_positions(segment_dict)

    def test_segment_values(self):
        self.vmat.analyze('mlcs')
        segment_dict = {
            0: {'r_dev': 0.437, 'r_corr': 100.89},
            2: {'r_dev': -0.405, 'r_corr': 100.04},
        }
        super().test_segment_values(segment_dict)


class Test_MLCS_105(VMATMixin, unittest.TestCase):
    """Tests of the result values of MLCS images at 105cm SID."""

    drmlc_105_open = osp.join(_vmat_test_files_dir, 'DRMLCopen-105-example.dcm')
    drmlc_105_dmlc = osp.join(_vmat_test_files_dir, 'DRMLCmlc-105-example.dcm')

    def setUp(self):
        self.vmat = VMAT.from_images((self.drmlc_105_dmlc, self.drmlc_105_open))

    def test_overall_passed(self):
        """Test that the overall pass flag is true for default settings"""
        super().test_overall_passed('mlcs')

    def test_segment_positions(self):
        """Test various values of the Segments."""
        self.vmat.analyze('mlcs')
        segment_dict = {0: Point(391, 384), 2: Point(552, 384)}
        super().test_segment_positions(segment_dict)

    def test_segment_values(self):
        self.vmat.analyze('mlcs')
        segment_dict = {
            0: {'r_dev': -0.040, 'r_corr': 100.83},
            2: {'r_dev': -0.021, 'r_corr': 100.85},
        }
        super().test_segment_values(segment_dict)

class Test_DRGS_105(VMATMixin, unittest.TestCase):
    """Tests of the result values of DRMLC images at 105cm SID."""
    drgs_105_open = osp.join(_vmat_test_files_dir, 'DRGSopen-105-example.dcm')
    drgs_105_dmlc = osp.join(_vmat_test_files_dir, 'DRGSmlc-105-example.dcm')

    def setUp(self):
        self.vmat = VMAT((self.drgs_105_dmlc, self.drgs_105_open))
        self.vmat.settings.x_offset = 20

    def test_overall_passed(self):
        """Test that the overall pass flag is true for default settings"""
        super().test_overall_passed('drgs')

    def test_segment_positions(self):
        """Test various values of the Segments."""
        self.vmat.analyze('drgs')
        segment_dict = {0: Point(371, 384), 2: Point(478, 384)}
        super().test_segment_positions(segment_dict)

    def test_segment_values(self):
        self.vmat.analyze('drgs')
        segment_dict = {
            0: {'r_dev': 0.780, 'r_corr': 102.43},
            4: {'r_dev': -0.282, 'r_corr': 101.357},
        }
        super().test_segment_values(segment_dict)