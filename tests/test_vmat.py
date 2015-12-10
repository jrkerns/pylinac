from functools import partial
import os.path as osp
from unittest import TestCase

import matplotlib.pyplot as plt

from pylinac.core.geometry import Point
from pylinac.vmat import VMAT, DMLC, OPEN, PROFILE, DRGS, DRMLC
from tests.utils import save_file

DEMO_DIR = osp.join(osp.dirname(osp.dirname(__file__)), 'pylinac', 'demo_files', 'vmat')
TEST_DIR = osp.join(osp.dirname(__file__), 'test_files', 'VMAT')

within_1 = partial(TestCase().assertAlmostEqual, delta=1)
within_01 = partial(TestCase().assertAlmostEqual, delta=0.1)


class TestGeneral(TestCase):
    """Generic tests for VMAT class."""
    def setUp(self):
        self.vmat = VMAT()

    def test_analyze_without_both_images_loaded(self):
        """Raise an error if both images aren't loaded when analyzing."""
        self.assertRaises(AttributeError, self.vmat.analyze, DRMLC)
        self.vmat.load_demo_image('drgs')
        self.vmat.analyze('drgs')  # shouldn't raise

    def test_img_loaded_tags(self):
        """Test the 'is_loaded' type tags."""
        self.assertFalse(self.vmat.open_img_is_loaded)
        self.assertFalse(self.vmat.dmlc_img_is_loaded)
        self.vmat.load_demo_image()
        self.assertTrue(self.vmat.open_img_is_loaded)
        self.assertTrue(self.vmat.dmlc_img_is_loaded)

    def test_img_inversion(self):
        """Check that the demo images indeed get inverted."""
        self.vmat.load_demo_image()
        top_corner_before = self.vmat.image_open.array[:20, :20].mean()
        self.vmat._check_img_inversion()
        top_corner_after = self.vmat.image_open.array[:20, :20].mean()
        self.assertNotEqual(top_corner_before, top_corner_after)

    def test_analyze_without_test_type(self):
        dmlc = osp.join(TEST_DIR, 'no_test_type_dmlc.dcm')
        opn = osp.join(TEST_DIR, 'no_test_type_open.dcm')
        self.vmat.load_images((dmlc, opn))

        with self.assertRaises(ValueError):
            self.vmat.analyze()

        # but will run when test type is passed
        self.vmat.analyze('drmlc')

    def test_failure_with_tight_tolerance(self):
        self.vmat.load_demo_image()
        self.vmat.analyze(tolerance=0.1)
        self.vmat.return_results()


class TestLoading(TestCase):
    """Tests of the various loading schemas."""

    def setUp(self):
        self.vmat = VMAT()

    def test_load_image(self):
        good_name = osp.join(TEST_DIR, 'no_test_type_dmlc.dcm')
        bad_name = osp.join(TEST_DIR, 'no_test_or_image_type_1.dcm')

        # image type can be determined from a good name
        self.vmat.load_image(good_name)

        # but not a bad one
        with self.assertRaises(ValueError):
            self.vmat.load_image(bad_name)

    def test_from_urls(self):
        urls = ['https://s3.amazonaws.com/assuranceqa-staging/uploads/imgs/DRGS_dmlc.dcm',
                'https://s3.amazonaws.com/assuranceqa-staging/uploads/imgs/DRGS_open.dcm']
        VMAT.from_urls(urls)  # shouldn't raise

    def test_passing_3_images(self):
        """Test passing the wrong number of images."""
        with self.assertRaises(ValueError):
            self.vmat.load_images(('', '', ''))

    def test_demo_image_loads(self):
        """Test the that demo images load properly."""
        # shouldn't raise
        self.vmat.load_demo_image(DRGS)
        self.vmat.load_demo_image(DRMLC)

    def test_from_zip(self):
        path = osp.join(TEST_DIR, 'DRMLC.zip')
        v = VMAT.from_zip(path)
        v.analyze()

    def test_loading_with_bad_names(self):
        one = osp.join(TEST_DIR, 'no_test_or_image_type_1.dcm')
        two = osp.join(TEST_DIR, 'no_test_or_image_type_2.dcm')
        with self.assertRaises(ValueError):
            self.vmat.load_images((one, two))

        # but will work when everything is specified
        self.vmat.load_image(one, OPEN)
        self.vmat.load_image(two, DMLC)
        self.vmat.analyze(DRMLC)


class TestPlottingSaving(TestCase):
    """Test the plotting and plot saving methods."""

    @classmethod
    def setUpClass(cls):
        cls.vmat = VMAT.from_demo_images()
        cls.vmat.analyze()

    def test_plot_analyzed_image(self):
        self.vmat.plot_analyzed_image()

    def test_save_analyzed_image(self):
        # save as normal file
        save_file(self.vmat.save_analyzed_image)
        # save from buffer
        save_file(self.vmat.save_analyzed_image, as_file_object='b')

    def test_plot_subimage(self):
        for subimage in (DMLC, OPEN, PROFILE):
            self.vmat.plot_analyzed_image(subimage)

    def test_save_subimage(self):
        for subimage in (DMLC, OPEN, PROFILE):
            save_file(self.vmat.save_analyzed_subimage, subimage)
            save_file(self.vmat.save_analyzed_subimage, subimage, as_file_object='bytes')
        # also do interactive
        save_file(self.vmat.save_analyzed_subimage, PROFILE, interactive=True)
        save_file(self.vmat.save_analyzed_subimage, PROFILE, interactive=True, as_file_object='str')



class VMATMixin:
    filepaths = ('open', 'dmlc')
    test_type = ''
    x_offset = 0
    segment_positions = {1: Point(100, 200)}
    segment_values = {
        0: {'r_dev': 0, 'r_corr': 100},
        4: {'r_dev': 0, 'r_corr': 100},
    }
    avg_abs_r_deviation = 0
    avg_r_deviation = 0
    max_r_deviation = 0
    passes = True

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def setUp(self):
        self.vmat = VMAT(self.filepaths)
        self.vmat.analyze(self.test_type, x_offset=self.x_offset)

    def test_overall_passed(self):
        self.vmat.analyze(self.test_type, x_offset=self.x_offset)
        self.assertEqual(self.vmat.passed, self.passes)

    def test_fail_with_tight_tolerance(self):
        self.vmat.analyze(self.test_type, tolerance=0.01, x_offset=self.x_offset)
        self.assertFalse(self.vmat.passed)

    def test_segment_positions(self):
        for key, value in self.segment_positions.items():
            within_1(self.vmat.segments[key].center.x, value.x)
            within_1(self.vmat.segments[key].center.y, value.y)

    def test_segment_values(self):
        for key, value in self.segment_values.items():
            within_01(self.vmat.segments[key].r_dev, value['r_dev'])
            within_01(self.vmat.segments[key].r_corr, value['r_corr'])

    def test_deviations(self):
        self.assertAlmostEqual(self.vmat.avg_abs_r_deviation, self.avg_abs_r_deviation, delta=0.05)
        self.assertAlmostEqual(self.vmat.avg_r_deviation, self.avg_r_deviation, delta=0.02)
        self.assertAlmostEqual(self.vmat.max_r_deviation, self.max_r_deviation, delta=0.1)


class TestDRGSDemo(VMATMixin, TestCase):
    """Tests of the result values of the DRGS demo images."""
    filepaths = (osp.join(DEMO_DIR, 'DRGS_dmlc.dcm'),
                 osp.join(DEMO_DIR, 'DRGS_open.dcm'))
    test_type = DRGS
    segment_positions = {0: Point(161, 192), 4: Point(314, 192)}
    segment_values = {
        0: {'r_dev': 0.965, 'r_corr': 101.85},
        4: {'r_dev': -0.459, 'r_corr': 100.42},
    }
    avg_abs_r_deviation = 0.46
    avg_r_deviation = 0
    max_r_deviation = 0.96
    x_offset = 20

    def test_demo(self):
        """Run the demo; no errors should arise."""
        self.vmat.run_demo_drgs()


class TestDRMLCDemo(VMATMixin, TestCase):
    """Tests of the result values of the DRMLC demo images."""
    filepaths = (osp.join(DEMO_DIR, 'DRMLC_dmlc.dcm'),
                 osp.join(DEMO_DIR, 'DRMLC_open.dcm'))
    test_type = DRMLC
    segment_positions = {0: Point(170, 192), 2: Point(285, 192)}
    segment_values = {
        0: {'r_dev': 0.437, 'r_corr': 100.89},
        2: {'r_dev': -0.405, 'r_corr': 100.04},
    }
    avg_abs_r_deviation = 0.38
    avg_r_deviation = 0
    max_r_deviation = 0.44

    def test_demo(self):
        self.vmat.run_demo_drmlc()


class TestDRMLC105(VMATMixin, TestCase):
    """Tests of the result values of MLCS images at 105cm SID."""
    filepaths = (osp.join(TEST_DIR, 'DRMLCopen-105-example.dcm'),
                 osp.join(TEST_DIR, 'DRMLCdmlc-105-example.dcm'))
    test_type = DRMLC
    segment_positions = {0: Point(391, 384), 2: Point(552, 384)}
    segment_values = {
        0: {'r_dev': -0.040, 'r_corr': 100.83},
        2: {'r_dev': -0.021, 'r_corr': 100.85},
    }
    avg_abs_r_deviation = 0.03
    avg_r_deviation = 0
    max_r_deviation = 0.04


class TestDRGS105(VMATMixin, TestCase):
    """Tests of the result values of DRMLC images at 105cm SID."""
    filepaths = (osp.join(TEST_DIR, 'DRGSopen-105-example.dcm'),
                 osp.join(TEST_DIR, 'DRGSdmlc-105-example.dcm'))
    test_type = DRGS
    x_offset = 20
    segment_positions = {0: Point(371, 384), 2: Point(478, 384)}
    segment_values = {
        0: {'r_dev': 0.780, 'r_corr': 102.43},
        4: {'r_dev': -0.282, 'r_corr': 101.357},
    }
    avg_abs_r_deviation = 0.34
    avg_r_deviation = 0
    max_r_deviation = 0.78


class TestDRMLC2(VMATMixin, TestCase):
    """Tests of the result values of MLCS images at 105cm SID."""
    filepaths = (osp.join(TEST_DIR, 'DRMLC#2_open.dcm'),
                 osp.join(TEST_DIR, 'DRMLC#2_dmlc.dcm'))
    test_type = DRMLC
    segment_positions = {0: Point(199, 192), 2: Point(275, 192)}
    segment_values = {
        0: {'r_dev': 0.40, 'r_corr': 101.06},
        2: {'r_dev': -0.49, 'r_corr': 100.16},
    }
    avg_abs_r_deviation = 0.4
    avg_r_deviation = 0
    max_r_deviation = -0.49


class TestDRGS2(VMATMixin, TestCase):
    """Tests of the result values of DRMLC images at 105cm SID."""
    filepaths = (osp.join(TEST_DIR, 'DRGS#2_open.dcm'),
                 osp.join(TEST_DIR, 'DRGS#2_dmlc.dcm'))
    test_type = DRGS
    x_offset = 12
    segment_positions = {0: Point(191, 192), 2: Point(242, 192)}
    segment_values = {
        0: {'r_dev': 1.3, 'r_corr': 103.0},
        4: {'r_dev': -0.8, 'r_corr': 100.86},
    }
    avg_abs_r_deviation = 0.7
    avg_r_deviation = 0
    max_r_deviation = 1.3