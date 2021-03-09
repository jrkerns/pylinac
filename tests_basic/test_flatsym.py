"""Tests for the flatsym module of pylinac."""
from unittest import TestCase
import os.path as osp
from functools import partial

from pylinac.core.exceptions import NotAnalyzed
from pylinac.core.io import retrieve_demo_file
from pylinac.flatsym import FlatSym

from tests_basic.utils import has_www_connection, LocationMixin, save_file

TEST_DIR = osp.join(osp.dirname(__file__), 'test_files', 'Flatness & Symmetry')


def create_instance():
    fs = FlatSym.from_demo_image()
    fs.analyze('varian', 'varian')
    return fs


class FlatSymTests(TestCase):

    def test_demo_is_reachable(self):
        if has_www_connection():
            file = retrieve_demo_file(url='flatsym_demo.dcm')
            self.assertTrue(osp.isfile(file))

    def test_demo_loads_properly(self):
        """Loading the demo shouldn't raise an error"""
        FlatSym.from_demo_image()  # shouldn't raise

    def test_demo_runs(self):
        FlatSym.run_demo()

    def test_profile_limits(self):
        """Extreme profile limits should not raise an error"""
        fs = FlatSym.from_demo_image()
        fs.analyze(flatness_method="varian", symmetry_method="varian", vert_position=0.5, horiz_position=0.5, vert_width=0, horiz_width=0)
        fs.analyze(flatness_method="varian", symmetry_method="varian", vert_position=0.0, horiz_position=0.0, vert_width=1, horiz_width=1)
        fs.analyze(flatness_method="varian", symmetry_method="varian", vert_position=1.0, horiz_position=1.0, vert_width=1, horiz_width=1)

    def test_analyze_sets_analyzed_flag(self):
        fs = create_instance()
        self.assertTrue(fs._is_analyzed)

    def test_flatness_methods(self):
        fs = FlatSym.from_demo_image()
        analyze = partial(fs.analyze, symmetry_method='varian')
        for method in ('varian', 'elekta', 'siemens', 'vom80', 'iec'):
            analyze(flatness_method=method)  # shouldn't raise

    def test_symmetry_methods(self):
        fs = FlatSym.from_demo_image()
        analyze = partial(fs.analyze, flatness_method='varian')
        for method in ('point difference', 'elekta', 'pdq iec'):
            analyze(symmetry_method=method)  # shouldn't raise

    def test_results(self):
        fs = create_instance()
        self.assertIsInstance(fs.results(), str)

    def test_results_fails_if_not_analyzed(self):
        fs = FlatSym.from_demo_image()
        with self.assertRaises(NotAnalyzed):
            fs.results()

    def test_plot_works(self):
        fs = create_instance()
        fs.plot_analyzed_image()
        fs.plot_analyzed_image(show=True)

    def test_plot_fails_if_not_analysed(self):
        fs = FlatSym.from_demo_image()
        with self.assertRaises(NotAnalyzed):
            fs.plot_analyzed_image()

    def test_pdf_gets_generated(self):
        fs = create_instance()
        save_file(fs.publish_pdf)

    def test_pdf_fails_if_not_analyzed(self):
        fs = FlatSym.from_demo_image()
        with self.assertRaises(NotAnalyzed):
            fs.publish_pdf('dummy.pdf')


class FlatSymBase(LocationMixin):
    dir_location = TEST_DIR
    sym_tolerance = 0.05
    flat_tolerance = 0.05
    apply_smoothing = None
    vert_symmetry = 0
    vert_flatness = 0
    horiz_symmetry = 0
    horiz_flatness = 0
    symmetry_method = 'varian'
    flatness_method = 'varian'
    vert_position = 0.5
    horiz_position = 0.5
    vert_width = 0
    horiz_width = 0
    print_results = False

    @classmethod
    def setUpClass(cls):
        cls.fs = FlatSym(cls.get_filename(), filter=cls.apply_smoothing)
        cls.fs.analyze(flatness_method=cls.flatness_method, symmetry_method=cls.symmetry_method,
                       vert_position=cls.vert_position, horiz_position=cls.horiz_position,
                       vert_width=cls.vert_width, horiz_width=cls.horiz_width)
        if cls.print_results:
            print(cls.fs.results())

    def test_vert_symmetry(self):
        self.assertAlmostEqual(self.fs.symmetry['vertical']['value'], self.vert_symmetry, delta=self.sym_tolerance)

    def test_horiz_symmetry(self):
        self.assertAlmostEqual(self.fs.symmetry['horizontal']['value'], self.horiz_symmetry, delta=self.sym_tolerance)

    def test_vert_flatness(self):
        self.assertAlmostEqual(self.fs.flatness['vertical']['value'], self.vert_flatness, delta=self.flat_tolerance)

    def test_horiz_flatness(self):
        self.assertAlmostEqual(self.fs.flatness['horizontal']['value'], self.horiz_flatness, delta=self.flat_tolerance)


class FlatSymDemo(FlatSymBase, TestCase):
    vert_flatness = 1.93
    vert_symmetry = 2.46
    horiz_flatness = 1.86
    horiz_symmetry = 2.99

    @classmethod
    def get_filename(cls):
        return retrieve_demo_file(url='flatsym_demo.dcm')


class FlatSym6X(FlatSymBase, TestCase):
    file_path = ['6x auto bulb 2.dcm']
    apply_smoothing = 5
    # independently verified
    vert_flatness = 1.5
    vert_symmetry = 0.4
    horiz_flatness = 1.4
    horiz_symmetry = 0.5


class FlatSym18X(FlatSymBase, TestCase):
    file_path = ['18x auto bulb2.dcm']
    # independently verified
    apply_smoothing = 5
    vert_flatness = 1.4
    vert_symmetry = 0.5
    horiz_flatness = 1.5
    horiz_symmetry = 0.5


class FlatSymWideDemo(FlatSymBase, TestCase):
    vert_width = 0.025
    horiz_width = 0.025
    vert_flatness = 1.84
    vert_symmetry = 2.47
    horiz_flatness = 1.84
    horiz_symmetry = 2.96

    @classmethod
    def get_filename(cls):
        return retrieve_demo_file(url='flatsym_demo.dcm')

