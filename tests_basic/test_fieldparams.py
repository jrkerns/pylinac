"""Tests for the flatsym module of pylinac."""
from unittest import TestCase
import os.path as osp
from functools import partial
import numpy as np

from pylinac.core.exceptions import NotAnalyzed
from pylinac.core.io import retrieve_demo_file
from pylinac.core.profile import SingleProfile
from pylinac.fieldparams import FieldParams
import pylinac.fieldparams as fp

from tests_basic.utils import has_www_connection, LocationMixin, save_file

TEST_DIR = osp.join(osp.dirname(__file__), 'test_files', 'Flatness & Symmetry')


def create_instance():
    fs = FieldParams.from_demo_image()
    fs.analyze('varian')
    return fs


class CalcParamTests(TestCase):

    """Values verified from BeamScheme file 2-Sep-2011-A.txt X profile"""

    profile = SingleProfile(np.array([
        2.18,2.68,3.27,4.36,5.83,9.12,15.44,63.73,94.96,97.26,98.38,98.78,98.86,99,98.89,98.98,98.8,98.95,98.9,98.52,
        98.05,97.31,96.26,95.38,94.59,94.53,94.47,94.46,94.49,94.57,94.7,95.18,95.51,96.51,97.32,97.82,97.95,97.99,
        97.98,98.2,98.33,98.31,98.33,98.1,97.7,95.9,92.2,36.68,12.18,8.02,4.92,3.97,3.01]))
    profile.dpmm = 0.2
    delta = 0.2

    def test_left_edge_50(self):
        fp.interpolate = True
        fp.norm = 'cax'
        self.assertAlmostEqual(fp.left_edge_50(self.profile), 96.71, delta=self.delta)

    def test_right_edge_50(self):
        fp.interpolate = True
        fp.norm = 'cax'
        self.assertAlmostEqual(fp.right_edge_50(self.profile), 103.96, delta=self.delta)

    def test_left_edge_infl(self):
        fp.interpolate = True
        fp.norm = 'cax'
        self.assertAlmostEqual(fp.left_edge_infl(self.profile), 96.2, delta=self.delta)

    def test_right_edge_infl(self):
        fp.interpolate = True
        fp.norm = 'cax'
        self.assertAlmostEqual(fp.right_edge_infl(self.profile), 103.9, delta=self.delta)

    def test_field_size_edge_50(self):
        fp.interpolate = True
        fp.norm = 'cax'
        self.assertAlmostEqual(fp.field_size_edge_50(self.profile), 200.75, delta=self.delta)

    def test_field_size_edge_infl(self):
        fp.interpolate = True
        fp.norm = 'cax'
        fp.pen_width = 0
        self.assertAlmostEqual(fp.field_size_edge_infl(self.profile), 200.29, delta=self.delta)

    def test_field_center_edge_50(self):
        fp.interpolate = True
        fp.norm = 'cax'
        self.assertAlmostEqual(fp.field_center_edge_50(self.profile), 3.67, delta=self.delta)

    def test_field_center_edge_infl(self):
        fp.interpolate = True
        fp.norm = 'cax'
        fp.pen_width = 0
        self.assertAlmostEqual(fp.field_center_edge_infl(self.profile), 3.67, delta=self.delta)

    def test_penumbra_left_80_20(self):
        fp.interpolate = True
        fp.norm = 'cax'
        self.assertAlmostEqual(fp.penumbra_left_80_20(self.profile), 6.54, delta=self.delta)

    def test_penumbra_right_80_20(self):
        fp.interpolate = True
        fp.norm = 'cax'
        self.assertAlmostEqual(fp.penumbra_right_80_20(self.profile), 7.13, delta=self.delta)

    def test_penumbra_left_infl(self):
        fp.interpolate = True
        fp.norm = 'cax'
        fp.pen_width = 0
        self.assertAlmostEqual(fp.penumbra_left_infl(self.profile), 6.37, delta=self.delta)

    def test_penumbra_right_infl(self):
        fp.interpolate = True
        fp.norm = 'cax'
        fp.pen_width = 0
        self.assertAlmostEqual(fp.penumbra_right_infl(self.profile), 5.37, delta=self.delta)

    def test_penumbra_slope_left_infl(self):
        fp.interpolate = True
        fp.norm = 'cax'
        fp.pen_width = 0
        self.assertAlmostEqual(fp.penumbra_slope_left_infl(self.profile), 12.48, delta=self.delta)

    def test_penumbra_slope_right_infl(self):
        fp.interpolate = True
        fp.norm = 'cax'
        fp.pen_width = 0
        self.assertAlmostEqual(fp.penumbra_slope_right_infl(self.profile), -15.12, delta=self.delta)

    def test_flatness_dose_difference(self):
        self.assertAlmostEqual(fp.flatness_dose_difference(self.profile, 0.8), 2.35, delta=self.delta)

    def test_flatness_dose_ratio(self):
        self.assertAlmostEqual(fp.flatness_dose_ratio(self.profile, 0.8), 104.80, delta=self.delta)

    def test_symmetry_point_difference(self):
        self.assertAlmostEqual(fp.symmetry_point_difference(self.profile, 0.8), 1.91, delta=self.delta)

    def test_symmetry_pdq_iec(self):
        self.assertAlmostEqual(fp.symmetry_pdq_iec(self.profile, 0.8), 101.88, delta=self.delta)

    def test_symmetry_area(self):
        self.assertAlmostEqual(fp.symmetry_area(self.profile, 1.0), 0.44, delta=self.delta)

    def test_deviation_max(self):
        fp.norm = 'cax'
        self.assertAlmostEqual(fp.deviation_max(self.profile, 0.8), 104.80, delta=self.delta)

    def test_deviation_diff(self):
        fp.norm = 'cax'
        self.assertAlmostEqual(fp.deviation_diff(self.profile, 0.8), 4.80, delta=self.delta)


class FieldParamTests(TestCase):

    def test_demo_is_reachable(self):
        if has_www_connection():
            file = retrieve_demo_file(url='flatsym_demo.dcm')
            self.assertTrue(osp.isfile(file))

    def test_demo_loads_properly(self):
        """Loading the demo shouldn't raise an error."""
        FieldParams.from_demo_image()  # shouldn't raise

    def test_demo_runs(self):
        FieldParams.run_demo()

    def test_profile_limits(self):
        """Extreme profile limits should not raise an error."""
        fs = FieldParams.from_demo_image()
        fs.analyze(protocol="varian", vert_position=0.5, horiz_position=0.5, vert_width=0, horiz_width=0)
        fs.analyze(protocol="varian", vert_position=0.0, horiz_position=0.0, vert_width=1, horiz_width=1)
        fs.analyze(protocol="varian", vert_position=1.0, horiz_position=1.0, vert_width=1, horiz_width=1)

    def test_analyze_sets_analyzed_flag(self):
        fs = create_instance()
        self.assertTrue(fs._is_analyzed)

    # profile.find peak does not raise error if peak is inverted.
    # def test_analyze_fails_when_incorrectly_inverted(self):
    #    fs = create_instance()
    #    with self.assertRaises(ValueError):
    #        fs.analyze('varian', invert=True)
    #    fs = create_instance()
    #    with self.assertRaises(ValueError):
    #        fs.analyze('elekta', invert=True)

    def test_protocols(self):
        fs = FieldParams.from_demo_image()
        analyze = partial(fs.analyze, protocol='varian')
        for method in ('all', 'default', 'varian', 'elekta', 'siemens', 'vom80', 'iec9076', 'din', 'afssaps-jorf','fff'):
            analyze(protocol=method)  # shouldn't raise

    def test_results(self):
        fs = create_instance()
        self.assertIsInstance(fs.results(), str)

    def test_results_fails_if_not_analyzed(self):
        fs = FieldParams.from_demo_image()
        with self.assertRaises(NotAnalyzed):
            fs.results()

    def test_plot_works(self):
        fs = create_instance()
        fs.plot_analyzed_image()
        fs.plot_analyzed_image(show=True)

    def test_plot_fails_if_not_analysed(self):
        fs = FieldParams.from_demo_image()
        with self.assertRaises(NotAnalyzed):
            fs.plot_analyzed_image()

    def test_pdf_gets_generated(self):
        fs = create_instance()
        save_file(fs.publish_pdf)

    def test_pdf_fails_if_not_analyzed(self):
        fs = FieldParams.from_demo_image()
        with self.assertRaises(NotAnalyzed):
            fs.publish_pdf('dummy.pdf')


class FieldParamsBase(LocationMixin):
    dir_location = TEST_DIR
    sym_tolerance = 0.05
    flat_tolerance = 0.05
    apply_smoothing = None
    vert_symmetry = 0
    vert_flatness = 0
    horiz_symmetry = 0
    horiz_flatness = 0
    protocol = 'varian'
    vert_position = 0.5
    horiz_position = 0.5
    vert_width = 0
    horiz_width = 0
    print_results = False

    @classmethod
    def setUpClass(cls):
        fp.norm = 'max grounded'
        fp.interpolate = False
        cls.fs = FieldParams(cls.get_filename(), filter=cls.apply_smoothing)
        cls.fs.analyze(protocol=cls.protocol, vert_position=cls.vert_position, horiz_position=cls.horiz_position,
                       vert_width=cls.vert_width, horiz_width=cls.horiz_width)
        if cls.print_results:
            print(cls.fs.results())

    def test_vert_symmetry(self):
        self.assertAlmostEqual(self.fs.parameters['vertical']['symmetry: {:.2f} %'], self.vert_symmetry, delta=self.sym_tolerance)

    def test_horiz_symmetry(self):
        self.assertAlmostEqual(self.fs.parameters['horizontal']['symmetry: {:.2f} %'], self.horiz_symmetry, delta=self.sym_tolerance)

    def test_vert_flatness(self):
        self.assertAlmostEqual(self.fs.parameters['vertical']['flatness: {:.2f} %'], self.vert_flatness, delta=self.flat_tolerance)

    def test_horiz_flatness(self):
        self.assertAlmostEqual(self.fs.parameters['horizontal']['flatness: {:.2f} %'], self.horiz_flatness, delta=self.flat_tolerance)


class FieldParamsDemo(FieldParamsBase, TestCase):
    fp.norm = 'max grounded'
    fp.interpolate = False
    vert_flatness = 1.93
    vert_symmetry = 2.46
    horiz_flatness = 1.86
    horiz_symmetry = 2.99

    @classmethod
    def get_filename(cls):
        return retrieve_demo_file(url='flatsym_demo.dcm')


class FlatSym6X(FieldParamsBase, TestCase):
    file_path = ['6x auto bulb 2.dcm']
    apply_smoothing = 5
    # independently verified
    vert_flatness = 1.5
    vert_symmetry = 0.4
    horiz_flatness = 1.4
    horiz_symmetry = 0.44


class FlatSym18X(FieldParamsBase, TestCase):
    file_path = ['18x auto bulb2.dcm']
    # independently verified
    apply_smoothing = 5
    vert_flatness = 1.4
    vert_symmetry = 0.44
    horiz_flatness = 1.5
    horiz_symmetry = 0.5


class FieldParamsWideDemo(FieldParamsBase, TestCase):
    vert_width = 0.025
    horiz_width = 0.025
    vert_flatness = 1.84
    vert_symmetry = 2.47
    horiz_flatness = 1.84
    horiz_symmetry = 2.96

    @classmethod
    def get_filename(cls):
        return retrieve_demo_file(url='flatsym_demo.dcm')
