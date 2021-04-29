"""Tests for the flatsym module of pylinac."""
import enum
import os.path as osp
from unittest import TestCase

from pylinac.core.exceptions import NotAnalyzed
from pylinac.core.io import retrieve_demo_file
from pylinac.core.profile import Edge, Normalization, Interpolation
from pylinac.field_analysis import FieldAnalysis, Protocol, DeviceFieldAnalysis, Device, Centering, \
    symmetry_point_difference, flatness_dose_difference, plot_flatness, plot_symmetry_point_difference
from tests_basic.utils import has_www_connection, LocationMixin, save_file

TEST_DIR = osp.join(osp.dirname(__file__), 'test_files', 'Flatness & Symmetry')


def create_instance():
    fs = FieldAnalysis.from_demo_image()
    fs.analyze(protocol=Protocol.VARIAN)
    return fs


class FieldAnalysisTests(TestCase):

    def test_demo_is_reachable(self):
        if has_www_connection():
            file = retrieve_demo_file(url='flatsym_demo.dcm')
            self.assertTrue(osp.isfile(file))

    def test_demo_loads_properly(self):
        """Loading the demo shouldn't raise an error"""
        FieldAnalysis.from_demo_image()  # shouldn't raise

    def test_demo_runs(self):
        FieldAnalysis.run_demo()  # shouldn't raise

    def test_profile_limits(self):
        """Extreme profile limits should not raise an error"""
        fs = FieldAnalysis.from_demo_image()
        fs.analyze(protocol=Protocol.VARIAN, centering=Centering.MANUAL, vert_position=0.5, horiz_position=0.5, vert_width=0, horiz_width=0)
        fs.analyze(protocol=Protocol.VARIAN, centering=Centering.MANUAL, vert_position=0.0, horiz_position=0.0, vert_width=1, horiz_width=1)
        fs.analyze(protocol=Protocol.VARIAN, centering=Centering.MANUAL, vert_position=1.0, horiz_position=1.0, vert_width=1, horiz_width=1)

    def test_analyze_sets_analyzed_flag(self):
        fs = create_instance()
        self.assertTrue(fs._is_analyzed)

    def test_results(self):
        fs = create_instance()
        self.assertIsInstance(fs.results(), str)

    def test_results_fails_if_not_analyzed(self):
        fs = FieldAnalysis.from_demo_image()
        with self.assertRaises(NotAnalyzed):
            fs.results()

    def test_plot_works(self):
        fs = create_instance()
        fs.plot_analyzed_image()
        fs.plot_analyzed_image(show=True)

    def test_plot_fails_if_not_analysed(self):
        fs = FieldAnalysis.from_demo_image()
        with self.assertRaises(NotAnalyzed):
            fs.plot_analyzed_image()

    def test_pdf_gets_generated(self):
        fs = create_instance()
        save_file(fs.publish_pdf)

    def test_pdf_fails_if_not_analyzed(self):
        fs = FieldAnalysis.from_demo_image()
        with self.assertRaises(NotAnalyzed):
            fs.publish_pdf('dummy.pdf')


class FieldAnalysisBase(LocationMixin):
    dir_location = TEST_DIR
    sym_tolerance = 0.05
    flat_tolerance = 0.05
    apply_smoothing = None
    vert_symmetry = 0
    vert_flatness = 0
    horiz_symmetry = 0
    horiz_flatness = 0
    vert_field_size = 100
    horiz_field_size = 100
    cax_to_top = 50
    cax_to_bottom = 50
    cax_to_left = 50
    cax_to_right = 50
    penum_top = 3.0
    penum_bottom = 3.0
    penum_right = 3.0
    penum_left = 3.0
    invert = False
    is_FFF = False
    penumbra = (20, 80)
    protocol = Protocol.VARIAN
    interpolation_method = Interpolation.LINEAR
    edge_detection_method = Edge.INFLECTION_DERIVATIVE
    centering = Centering.BEAM_CENTER
    normalization_method = Normalization.GEOMETRIC_CENTER
    in_field_ratio = 0.8
    slope_exclusion_ratio = 0.2
    vert_position = 0.5
    horiz_position = 0.5
    vert_width = 0
    horiz_width = 0
    top_slope = 0
    bottom_slope = 0
    left_slope = 0
    right_slope = 0
    top_to_beam_center_vert = 0
    top_to_beam_center_horiz = 0
    print_results = False

    @classmethod
    def setUpClass(cls):
        cls.fs = FieldAnalysis(cls.get_filename(), filter=cls.apply_smoothing)
        cls.fs.analyze(protocol=cls.protocol, centering=cls.centering,
                       vert_position=cls.vert_position, horiz_position=cls.horiz_position,
                       vert_width=cls.vert_width, horiz_width=cls.horiz_width, in_field_ratio=cls.in_field_ratio,
                       slope_exclusion_ratio=cls.slope_exclusion_ratio, invert=cls.invert, is_FFF=cls.is_FFF,
                       penumbra=cls.penumbra, interpolation=cls.interpolation_method, edge_detection_method=cls.edge_detection_method,
                       )
        if cls.print_results:
            print(cls.fs.results())

    def test_top_slope(self):
        if self.is_FFF:
            self.assertAlmostEqual(self.fs.results_data()['top slope'], self.top_slope, delta=0.1)

    def test_bottom_slope(self):
        if self.is_FFF:
            self.assertAlmostEqual(self.fs.results_data()['bottom slope'], self.bottom_slope, delta=0.1)

    def test_left_slope(self):
        if self.is_FFF:
            self.assertAlmostEqual(self.fs.results_data()['left slope'], self.left_slope, delta=0.1)

    def test_right_slope(self):
        if self.is_FFF:
            self.assertAlmostEqual(self.fs.results_data()['right slope'], self.right_slope, delta=0.1)

    def test_cax_to_top(self):
        self.assertAlmostEqual(self.fs.results_data()['CAX->Top (mm)'], self.cax_to_top, delta=0.3)

    def test_cax_to_bottom(self):
        self.assertAlmostEqual(self.fs.results_data()['CAX->Bottom (mm)'], self.cax_to_bottom, delta=0.3)

    def test_cax_to_left(self):
        self.assertAlmostEqual(self.fs.results_data()['CAX->Left (mm)'], self.cax_to_left, delta=0.3)

    def test_cax_to_right(self):
        self.assertAlmostEqual(self.fs.results_data()['CAX->Right (mm)'], self.cax_to_right, delta=0.3)

    def test_penumbra_bottom(self):
        self.assertAlmostEqual(self.fs.results_data()['bottom penumbra (mm)'], self.penum_bottom, delta=0.3)

    def test_penumbra_top(self):
        self.assertAlmostEqual(self.fs.results_data()['top penumbra (mm)'], self.penum_top, delta=0.3)

    def test_penumbra_left(self):
        self.assertAlmostEqual(self.fs.results_data()['left penumbra (mm)'], self.penum_left, delta=0.3)

    def test_penumbra_right(self):
        self.assertAlmostEqual(self.fs.results_data()['right penumbra (mm)'], self.penum_right, delta=0.3)

    def test_horiz_field_size(self):
        self.assertAlmostEqual(self.fs.results_data()['field size horizontal (mm)'], self.horiz_field_size, delta=1)

    def test_vert_field_size(self):
        self.assertAlmostEqual(self.fs.results_data()['field size vertical (mm)'], self.vert_field_size, delta=1)

    def test_vert_symmetry(self):
        self.assertAlmostEqual(self.fs.results_data()['symmetry vertical'], self.vert_symmetry, delta=self.sym_tolerance)

    def test_horiz_symmetry(self):
        self.assertAlmostEqual(self.fs.results_data()['symmetry horizontal'], self.horiz_symmetry, delta=self.sym_tolerance)

    def test_vert_flatness(self):
        self.assertAlmostEqual(self.fs.results_data()['flatness vertical'], self.vert_flatness, delta=self.flat_tolerance)

    def test_horiz_flatness(self):
        self.assertAlmostEqual(self.fs.results_data()['flatness horizontal'], self.horiz_flatness, delta=self.flat_tolerance)


class DeviceAnalysisBase(FieldAnalysisBase):
    edge_detection_method = Edge.INFLECTION_HILL

    @classmethod
    def setUpClass(cls):
        cls.fs = DeviceFieldAnalysis(cls.get_filename(), device=Device.PROFILER)
        cls.fs.analyze(protocol=cls.protocol,
                       in_field_ratio=cls.in_field_ratio,
                       slope_exclusion_ratio=cls.slope_exclusion_ratio, is_FFF=cls.is_FFF,
                       penumbra=cls.penumbra, interpolation=cls.interpolation_method, edge_detection_method=cls.edge_detection_method,
                       )
        if cls.print_results:
            print(cls.fs.results())


class Profiler6FFF(DeviceAnalysisBase, TestCase):
    file_path = ['6fff.prm']
    is_FFF = True
    vert_flatness = 14.35
    vert_symmetry = -0.9
    horiz_flatness = 14.6
    horiz_symmetry = 0.5
    vert_position = 0.6
    vert_field_size = 255.8
    horiz_field_size = 246.0
    cax_to_top = 128.3
    cax_to_left = 122.7
    cax_to_right = 123.3
    cax_to_bottom = 127.4
    penum_top = 6.1
    penum_bottom = 6.9
    penum_right = 5.8
    penum_left = 6.4
    top_slope = 0.3
    bottom_slope = -0.3
    left_slope = 0.3
    right_slope = -0.3


class FlatSymDemo(FieldAnalysisBase, TestCase):
    # independently verified
    vert_flatness = 1.7
    vert_symmetry = -2.6
    horiz_flatness = 1.85
    horiz_symmetry = -2.99
    vert_position = 0.6
    vert_field_size = 200.2
    horiz_field_size = 141.0
    cax_to_top = 99.8
    cax_to_left = 60.4
    cax_to_right = 80.5
    cax_to_bottom = 100.5
    penum_top = 3.9
    penum_bottom = 2.8
    penum_right = 3.0
    penum_left = 2.7

    @classmethod
    def get_filename(cls):
        return retrieve_demo_file(url='flatsym_demo.dcm')


class FlatSymWideDemo(FlatSymDemo, TestCase):
    vert_width = 0.025
    horiz_width = 0.025
    vert_flatness = 1.6


class FlatSym6X(FieldAnalysisBase, TestCase):
    file_path = ['6x auto bulb 2.dcm']
    # independently verified
    horiz_width = 0.01
    vert_width = 0.01
    vert_flatness = 1.45
    vert_symmetry = -0.4
    horiz_flatness = 1.4
    horiz_symmetry = -0.4
    vert_field_size = 99.5
    horiz_field_size = 99.5
    cax_to_top = 49.9
    cax_to_left = 49.2
    cax_to_right = 50.3
    cax_to_bottom = 49.6
    penum_top = 3.0
    penum_bottom = 3.0
    penum_right = 2.5
    penum_left = 2.8


class FlatSym18X(FieldAnalysisBase, TestCase):
    file_path = ['18x auto bulb2.dcm']
    # independently verified
    horiz_width = 0.01
    vert_width = 0.01
    vert_flatness = 1.4
    vert_symmetry = 0.4
    horiz_flatness = 1.5
    horiz_symmetry = -0.5
    vert_field_size = 99.4
    horiz_field_size = 99.5
    cax_to_top = 49.9
    cax_to_left = 49.2
    cax_to_right = 50.3
    cax_to_bottom = 49.5
    penum_top = 3.7
    penum_bottom = 3.4
    penum_right = 3.0
    penum_left = 3.0


class TestCustomProtocol(TestCase):

    def test_custom_protocol(self):

        class MyProtocol(enum.Enum):
            Awesomeness = {
                'symmetry': {'calc': symmetry_point_difference, 'unit': '%', 'plot': plot_symmetry_point_difference},
                'flatness': {'calc': flatness_dose_difference, 'unit': '%', 'plot': plot_flatness},
            }

        fa = FieldAnalysis.from_demo_image()
        fa.analyze(protocol=MyProtocol.Awesomeness)  # shouldn't raise
        fa.results()

