"""Tests for the flatsym module of pylinac."""
import enum
import io
import os.path as osp
from unittest import TestCase

from pylinac.core.exceptions import NotAnalyzed
from pylinac.core.io import retrieve_demo_file
from pylinac.core.profile import Edge, Normalization, Interpolation
from pylinac.field_analysis import FieldAnalysis, Protocol, DeviceFieldAnalysis, Device, Centering, \
    symmetry_point_difference, flatness_dose_difference, plot_flatness, plot_symmetry_point_difference, FieldResult
from tests_basic.utils import has_www_connection, CloudFileMixin, save_file, get_file_from_cloud_test_repo

TEST_DIR = 'flatness_symmetry'


def create_instance():
    fs = FieldAnalysis.from_demo_image()
    fs.analyze()
    return fs


class FieldAnalysisTests(TestCase):

    def test_load_from_file_object(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, '6x-auto-bulb-2.dcm'])
        ref_fa = FieldAnalysis(path)
        ref_fa.analyze()
        with open(path, 'rb') as f:
            fa = FieldAnalysis(f)
            fa.analyze()
        self.assertIsInstance(fa, FieldAnalysis)
        self.assertEqual(fa.image.shape, ref_fa.image.shape)

    def test_load_from_stream(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, '6x-auto-bulb-2.dcm'])
        ref_fa = FieldAnalysis(path)
        ref_fa.analyze()
        with open(path, 'rb') as f:
            s = io.BytesIO(f.read())
            fa = FieldAnalysis(s)
            fa.analyze()
        self.assertIsInstance(fa, FieldAnalysis)
        self.assertEqual(fa.image.shape, ref_fa.image.shape)

    def test_demo_is_reachable(self):
        if has_www_connection():
            file = retrieve_demo_file(url='flatsym_demo.dcm', force=True)
            self.assertTrue(osp.isfile(file))

    def test_demo_loads_properly(self):
        """Loading the demo shouldn't raise an error"""
        FieldAnalysis.from_demo_image()  # shouldn't raise

    def test_demo_runs(self):
        FieldAnalysis.run_demo()  # shouldn't raise

    def test_profile_limits(self):
        """Extreme profile limits should not raise an error"""
        fs = FieldAnalysis.from_demo_image()
        fs.analyze()
        fs.analyze()
        fs.analyze()

    def test_analyze_sets_analyzed_flag(self):
        fs = create_instance()
        self.assertTrue(fs._is_analyzed)

    def test_results(self):
        fs = create_instance()
        self.assertIsInstance(fs.results(), str)

    def test_results_data(self):
        fs = create_instance()
        fs.analyze()
        data = fs.results_data()
        self.assertIsInstance(data, FieldResult)
        self.assertEqual(data.field_size_vertical_mm, fs._results['field_size_vertical_mm'])
        self.assertEqual(data.protocol_results['flatness_vertical'], fs._extra_results['flatness_vertical'])

        data_dict = fs.results_data(as_dict=True)
        self.assertIsInstance(data_dict, dict)
        self.assertEqual(data_dict['protocol_results']['flatness_vertical'], fs._extra_results['flatness_vertical'])

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

    def test_save_analyzed_image(self):
        fa = create_instance()
        save_file(fa.plot_analyzed_image)


class FieldAnalysisBase(CloudFileMixin):
    dir_path = ['flatness_symmetry']
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
            self.assertAlmostEqual(self.fs.results_data().top_slope_percent_mm, self.top_slope, delta=0.1)

    def test_bottom_slope(self):
        if self.is_FFF:
            self.assertAlmostEqual(self.fs.results_data().bottom_slope_percent_mm, self.bottom_slope, delta=0.1)

    def test_left_slope(self):
        if self.is_FFF:
            self.assertAlmostEqual(self.fs.results_data().left_slope_percent_mm, self.left_slope, delta=0.1)

    def test_right_slope(self):
        if self.is_FFF:
            self.assertAlmostEqual(self.fs.results_data().right_slope_percent_mm, self.right_slope, delta=0.1)

    def test_cax_to_top(self):
        self.assertAlmostEqual(self.fs.results_data().cax_to_top_mm, self.cax_to_top, delta=0.3)

    def test_cax_to_bottom(self):
        self.assertAlmostEqual(self.fs.results_data().cax_to_bottom_mm, self.cax_to_bottom, delta=0.3)

    def test_cax_to_left(self):
        self.assertAlmostEqual(self.fs.results_data().cax_to_left_mm, self.cax_to_left, delta=0.3)

    def test_cax_to_right(self):
        self.assertAlmostEqual(self.fs.results_data().cax_to_right_mm, self.cax_to_right, delta=0.3)

    def test_penumbra_bottom(self):
        self.assertAlmostEqual(self.fs.results_data().bottom_penumbra_mm, self.penum_bottom, delta=0.3)

    def test_penumbra_top(self):
        self.assertAlmostEqual(self.fs.results_data().top_penumbra_mm, self.penum_top, delta=0.3)

    def test_penumbra_left(self):
        self.assertAlmostEqual(self.fs.results_data().left_penumbra_mm, self.penum_left, delta=0.3)

    def test_penumbra_right(self):
        self.assertAlmostEqual(self.fs.results_data().right_penumbra_mm, self.penum_right, delta=0.3)

    def test_horiz_field_size(self):
        self.assertAlmostEqual(self.fs.results_data().field_size_horizontal_mm, self.horiz_field_size, delta=1)

    def test_vert_field_size(self):
        self.assertAlmostEqual(self.fs.results_data().field_size_vertical_mm, self.vert_field_size, delta=1)

    def test_vert_symmetry(self):
        self.assertAlmostEqual(self.fs.results_data().protocol_results['symmetry_vertical'], self.vert_symmetry, delta=self.sym_tolerance)

    def test_horiz_symmetry(self):
        self.assertAlmostEqual(self.fs.results_data().protocol_results['symmetry_horizontal'], self.horiz_symmetry, delta=self.sym_tolerance)

    def test_vert_flatness(self):
        self.assertAlmostEqual(self.fs.results_data().protocol_results['flatness_vertical'], self.vert_flatness, delta=self.flat_tolerance)

    def test_horiz_flatness(self):
        self.assertAlmostEqual(self.fs.results_data().protocol_results['flatness_horizontal'], self.horiz_flatness, delta=self.flat_tolerance)


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
    file_name = '6fff.prm'
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
    file_name = '6x-auto-bulb-2.dcm'
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
    file_name = '18x-auto-bulb2.dcm'
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
        fa.analyze()  # shouldn't raise
        fa.results()

