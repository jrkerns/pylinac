from functools import partial
from os import path as osp
from typing import Union, List
from unittest import TestCase

from pylinac.core.geometry import Point
from pylinac.core.io import retrieve_demo_file
from pylinac import DRGS, DRMLC

from tests_basic.utils import save_file, has_www_connection

TEST_DIR = osp.join(osp.dirname(__file__), 'test_files', 'VMAT')

within_5 = partial(TestCase().assertAlmostEqual, delta=5)
within_1 = partial(TestCase().assertAlmostEqual, delta=1)


class TestLoadingBase:
    demo_name = ''
    klass = object

    def test_demo_is_reachable(self):
        if has_www_connection():
            file = retrieve_demo_file(url=self.demo_name)
            self.assertTrue(osp.isfile(file))

    def test_can_load_demo(self):
        instance = self.klass.from_demo_images()
        self.assertIsInstance(instance, self.klass)

    def test_normal_instantiation(self):
        one = osp.join(TEST_DIR, 'no_test_or_image_type_1.dcm')
        two = osp.join(TEST_DIR, 'no_test_or_image_type_2.dcm')
        instance = self.klass(image_paths=(one, two))
        self.assertIsInstance(instance, self.klass)

    def test_from_url(self):
        url = f'https://s3.amazonaws.com/pylinac/{self.demo_name}'
        instance = self.klass.from_url(url=url)
        self.assertIsInstance(instance, self.klass)

    def test_passing_3_images_fails(self):
        """Test passing the wrong number of images."""
        with self.assertRaises(ValueError):
            self.klass(image_paths=('', '', ''))

    def test_print_results(self):
        instance = self.klass.from_demo_images()
        instance.analyze()
        self.assertIsInstance(instance.results(), str)

    def test_plot_analyzed_image(self):
        instance = self.klass.from_demo_images()
        instance.analyze()
        instance.plot_analyzed_image()  # shouldn't raise

    def test_publish_pdf(self):
        instance = self.klass.from_demo_images()
        instance.analyze()
        save_file(instance.publish_pdf)


class TestDRGSLoading(TestLoadingBase, TestCase):
    demo_name = 'drgs.zip'
    klass = DRGS


class TestDRMLCLoading(TestLoadingBase, TestCase):
    demo_name = 'drmlc.zip'
    klass = DRMLC


class VMATMixin:
    klass = object
    filepaths = Union[str, List]
    is_zip = False
    segment_positions = {1: Point(100, 200)}
    segment_values = {
        0: {'r_dev': 0, 'r_corr': 100},
        4: {'r_dev': 0, 'r_corr': 100},
    }
    avg_abs_r_deviation = 0
    avg_r_deviation = 0
    max_r_deviation = 0
    passes = True
    print_debug = False

    @classmethod
    def absolute_path(cls):
        if cls.is_zip:
            path = osp.join(TEST_DIR, cls.filepaths)
        else:
            path = [osp.join(TEST_DIR, path) for path in cls.filepaths]
        return path

    def setUp(self):
        if self.is_zip:
            self.vmat = self.klass.from_zip(self.absolute_path())
        else:
            self.vmat = self.klass(self.absolute_path())
        self.vmat.analyze()
        if self.print_debug:
            print(self.vmat.results())
            print(f"Segment 0: rdev {self.vmat.segments[0].r_dev:2.3f}, rcorr {self.vmat.segments[0].r_corr:2.3f}")
            if self.klass == DRGS:
                print(f"Segment 4: rdev {self.vmat.segments[4].r_dev:2.3f}, rcorr {self.vmat.segments[4].r_corr:2.3f}")
            else:
                print(f"Segment 2: rdev {self.vmat.segments[2].r_dev:2.3f}, rcorr {self.vmat.segments[2].r_corr:2.3f}")
            print("Max dev", self.vmat.max_r_deviation)

    def test_overall_passed(self):
        self.assertEqual(self.vmat.passed, self.passes)

    def test_fail_with_tight_tolerance(self):
        self.vmat.analyze(tolerance=0.01)
        self.assertFalse(self.vmat.passed)

    def test_segment_positions(self):
        for key, value in self.segment_positions.items():
            within_5(self.vmat.segments[key].center.x, value.x)
            within_5(self.vmat.segments[key].center.y, value.y)

    def test_segment_values(self):
        for key, value in self.segment_values.items():
            within_1(self.vmat.segments[key].r_dev, value['r_dev'])
            within_1(self.vmat.segments[key].r_corr, value['r_corr'])

    def test_deviations(self):
        self.assertAlmostEqual(self.vmat.avg_abs_r_deviation, self.avg_abs_r_deviation, delta=0.05)
        self.assertAlmostEqual(self.vmat.avg_r_deviation, self.avg_r_deviation, delta=0.02)
        self.assertAlmostEqual(self.vmat.max_r_deviation, self.max_r_deviation, delta=0.1)

    def test_different_segment_size_is_nearly_the_same(self):
        segment_width_mm = 10
        segment_height_mm = 50
        self.vmat.analyze(segment_size_mm=(segment_width_mm, segment_height_mm))
        self.assertTrue(self.vmat.segments[0]._nominal_width_mm, segment_width_mm)
        self.assertTrue(self.vmat.segments[0]._nominal_height_mm, segment_height_mm)


class TestDRGSDemo(VMATMixin, TestCase):
    """Tests of the result values of the DRGS demo images."""
    segment_positions = {0: Point(161, 192), 4: Point(314, 192)}
    segment_values = {
        0: {'r_dev': 0.965, 'r_corr': 6.2},
        4: {'r_dev': -0.459, 'r_corr': 6},
    }
    avg_abs_r_deviation = 0.66
    max_r_deviation = 1.8
    passes = False

    def setUp(self):
        self.vmat = DRGS.from_demo_images()
        self.vmat.analyze()

    def test_demo(self):
        """Run the demo; no errors should arise."""
        self.vmat.run_demo()


class TestDRMLCDemo(VMATMixin, TestCase):
    """Tests of the result values of the DRMLC demo images."""
    segment_positions = {0: Point(170, 192), 2: Point(285, 192)}
    segment_values = {
        0: {'r_dev': -0.7, 'r_corr': 5.7},
        2: {'r_dev': -0.405, 'r_corr': 5.8},
    }
    avg_abs_r_deviation = 0.44
    max_r_deviation = 0.89

    def setUp(self):
        self.vmat = DRMLC.from_demo_images()
        self.vmat.analyze()

    def test_demo(self):
        self.vmat.run_demo()


class TestDRMLC105(VMATMixin, TestCase):
    """Tests of the result values of MLCS images at 105cm SID."""
    klass = DRMLC
    filepaths = ('DRMLCopen-105-example.dcm', 'DRMLCdmlc-105-example.dcm')
    segment_positions = {0: Point(391, 384), 2: Point(552, 384)}
    segment_values = {
        0: {'r_dev': -2.1, 'r_corr': 13.6},
        2: {'r_dev': 0.22, 'r_corr': 14},
    }
    avg_abs_r_deviation = 1.06
    max_r_deviation = 2.11
    passes = False


class TestDRGS105(VMATMixin, TestCase):
    """Tests of the result values of DRMLC images at 105cm SID."""
    filepaths = ('DRGSopen-105-example.dcm', 'DRGSdmlc-105-example.dcm')
    klass = DRGS
    segment_positions = {0: Point(371, 384), 2: Point(478, 384)}
    segment_values = {
        0: {'r_dev': 1.385, 'r_corr': 15.12},
        4: {'r_dev': -0.8, 'r_corr': 14.8},
    }
    avg_abs_r_deviation = 0.68
    max_r_deviation = 1.38


class TestDRMLC2(VMATMixin, TestCase):
    """Tests of the result values of MLCS images at 105cm SID."""
    filepaths = ('DRMLC#2_open.dcm', 'DRMLC#2_dmlc.dcm')
    klass = DRMLC
    segment_positions = {0: Point(199, 192), 2: Point(275, 192)}
    segment_values = {
        0: {'r_dev': 0.77, 'r_corr': 6.1},
        2: {'r_dev': -1.1, 'r_corr': 6},
    }
    avg_abs_r_deviation = 1.4
    max_r_deviation = 1.98
    passes = False


class TestDRGS2(VMATMixin, TestCase):
    """Tests of the result values of DRMLC images at 105cm SID."""
    filepaths = ('DRGS#2_open.dcm', 'DRGS#2_dmlc.dcm')
    klass = DRGS
    segment_positions = {0: Point(191, 192), 2: Point(242, 192)}
    segment_values = {
        0: {'r_dev': 1.5, 'r_corr': 6.4},
        4: {'r_dev': -0.7, 'r_corr': 6.3},
    }
    avg_abs_r_deviation = 0.7
    max_r_deviation = 1.5
