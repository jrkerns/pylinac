from unittest import TestCase
from pathlib import Path

import matplotlib.pyplot as plt

from pylinac import CatPhan503, CatPhan504, CatPhan600, CatPhan604
from pylinac.ct import CTP404CP504, CTP404CP503, CTP515, CTP528CP503, CTP528CP504
from pylinac.core.geometry import Point
from tests_basic.utils import save_file, LoadingTestBase, LocationMixin

TEST_DIR = Path('./test_files/CBCT')


class CBCTLoading(LoadingTestBase, TestCase):
    klass = CatPhan504
    constructor_input = Path(TEST_DIR, 'Pelvis')
    demo_load_method = 'from_demo_images'
    url = 'CatPhan504.zip'
    zip = Path(TEST_DIR, 'CBCT_4.zip')


class GeneralTests(TestCase):
    """Test general things when using cbct module."""

    def setUp(self):
        self.cbct = CatPhan504.from_demo_images()

    def test_demo(self):
        """Run the demo to make sure it works."""
        self.cbct.run_demo()

    def test_helpers(self):
        """Test the various helper methods."""
        self.cbct.analyze()
        self.cbct._results()

    def test_phan_center(self):
        """Test locations of the phantom center."""
        known_phan_center = Point(257, 255)
        self.cbct.analyze()
        self.assertAlmostEqual(self.cbct.ctp404.phan_center.x, known_phan_center.x, delta=0.7)
        self.assertAlmostEqual(self.cbct.ctp404.phan_center.y, known_phan_center.y, delta=0.7)


class CustomPhantom(TestCase):

    def test_removing_module(self):
        class CatPhan504Modified(CatPhan504):
            modules = {
                CTP404CP504: {'offset': 0},
            }
        ct = CatPhan504Modified.from_demo_images()  # shouldn't raise
        ct.analyze()

    def test_omitting_404_module(self):
        class CatPhan504Modified(CatPhan504):
            modules = {
                CTP528CP504: {'offset': -30},
            }
        ct = CatPhan504Modified.from_demo_images()  # shouldn't raise
        with self.assertRaises(ValueError):
            ct.analyze()

    def test_modifying_module(self):
        class CTP528Custom(CTP528CP503):
            start_angle = 0.1

        class CatPhan503Modified(CatPhan503):
            modules = {
                CTP404CP503: {'offset': 0},
                CTP528Custom: {'offset': -30},
            }
        ct = CatPhan503Modified.from_demo_images()  # shouldn't raise
        ct.analyze()
        ct.plot_analyzed_image()


class DemoTests(TestCase):

    def test_504(self):
        CatPhan504.run_demo()

    def test_503(self):
        CatPhan503.run_demo()

    def test_600(self):
        CatPhan600.run_demo()

    def test_604(self):
        CatPhan604.run_demo()


class PlottingSaving(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cbct = CatPhan504.from_demo_images()
        cls.cbct.analyze()

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def test_save_image(self):
        """Test that saving an image does something."""
        for method in ['save_analyzed_image', 'save_analyzed_subimage']:
            methodcall = getattr(self.cbct, method)
            save_file(methodcall)

    def test_plot_images(self):
        """Test the various plotting functions."""
        self.cbct.plot_analyzed_image()
        for item in ['hu', 'un', 'mtf', 'sp', 'prof', 'lin', 'lc']:
            self.cbct.plot_analyzed_subimage(item)

        self.cbct.plot_analyzed_subimage('lin', delta=False)

        with self.assertRaises(ValueError):
            self.cbct.plot_analyzed_subimage('sr')


class CatPhanMixin(LocationMixin):
    """A mixin to use for testing Varian CBCT scans; does not inherit from TestCase as it would be run
        otherwise."""
    catphan = CatPhan504
    check_uid = True
    origin_slice = 0
    file_path = []
    dir_location = TEST_DIR
    hu_tolerance = 40
    scaling_tolerance = 1
    zip = True
    expected_roll = 0
    hu_values = {}
    unif_values = {}
    mtf_values = {}
    avg_line_length = 50
    slice_thickness = 2
    lowcon_visible = 0
    print_debug = False

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        if cls.zip:
            cls.cbct = cls.catphan.from_zip(filename)
        else:
            cls.cbct = cls.catphan(filename)
        cls.cbct.analyze(cls.hu_tolerance, cls.scaling_tolerance)
        if cls.print_debug:
            print(cls.cbct._results())

    @classmethod
    def tearDownClass(cls):
        # somewhere there is a memory leak if ``cbct`` isn't deleted.
        delattr(cls, 'cbct')

    def test_slice_thickness(self):
        """Test the slice thickness."""
        self.assertAlmostEqual(self.cbct.ctp404.meas_slice_thickness, self.slice_thickness, delta=0.3)

    def test_lowcontrast_bubbles(self):
        """Test the number of low contrast bubbles visible."""
        if not isinstance(self.cbct, CatPhan503):
            self.assertAlmostEqual(self.cbct.ctp515.rois_visible, self.lowcon_visible, delta=1)

    def test_slice_locations(self):
        """Test the locations of the slices of interest."""
        self.assertAlmostEqual(self.cbct.origin_slice, self.origin_slice, delta=1)

    def test_phantom_roll(self):
        """Test the roll of the phantom."""
        self.assertAlmostEqual(self.cbct.catphan_roll, self.expected_roll, delta=0.3)

    def test_HU_values(self):
        """Test HU values."""
        for key, roi in self.cbct.ctp404.rois.items():
            exp_val = self.hu_values[key]
            meas_val = roi.pixel_value
            self.assertAlmostEqual(exp_val, meas_val, delta=5)

    def test_uniformity_values(self):
        """Test Uniformity HU values."""
        for key, exp_val in self.unif_values.items():
            meas_val = self.cbct.ctp486.rois[key].pixel_value
            self.assertAlmostEqual(exp_val, meas_val, delta=5)

    def test_geometry_line_length(self):
        """Test the geometry distances."""
        self.assertAlmostEqual(self.avg_line_length, self.cbct.ctp404.avg_line_length, delta=0.1)

    def test_MTF_values(self):
        """Test MTF values."""
        for key, exp_mtf in self.mtf_values.items():
            meas_mtf = self.cbct.ctp528.mtf.relative_resolution(key)
            self.assertAlmostEqual(exp_mtf, meas_mtf, delta=0.1)

    def test_pdf(self):
        save_file(self.cbct.publish_pdf, 'temp')


class CatPhanDemo(CatPhanMixin, TestCase):
    """Test the CBCT demo (Varian high quality head protocol)."""
    expected_roll = -0.3
    origin_slice = 32
    hu_values = {'Poly': -45, 'Acrylic': 117, 'Delrin': 341, 'Air': -998, 'Teflon': 997, 'PMP': -200, 'LDPE': -103}
    unif_values = {'Center': 17, 'Left': 10, 'Right': 0, 'Top': 6, 'Bottom': 6}
    mtf_values = {50: 0.56}
    avg_line_length = 49.92
    lowcon_visible = 3
    slice_thickness = 2.5

    @classmethod
    def setUpClass(cls):
        cls.cbct = CatPhan504.from_demo_images()
        cls.cbct.analyze()


class CatPhan4(CatPhanMixin, TestCase):
    """A Varian CBCT dataset"""
    file_path = ['CBCT_4.zip']
    expected_roll = -2.57
    origin_slice = 31
    hu_values = {'Poly': -33, 'Acrylic': 119, 'Delrin': 335, 'Air': -979, 'Teflon': 970, 'PMP': -185, 'LDPE': -94}
    unif_values = {'Center': 17, 'Left': 10, 'Right': 22, 'Top': 13, 'Bottom': 18}
    mtf_values = {80: 0.24, 90: 0.2, 60: 0.31, 70: 0.23, 95: 0.15}
    lowcon_visible = 1
    slice_thickness = 2.4


class Elekta2(CatPhanMixin, TestCase):
    """An Elekta CBCT dataset"""
    catphan = CatPhan503
    file_path = ['Elekta_2.zip']
    origin_slice = 162
    hu_values = {'Poly': -319, 'Acrylic': -224, 'Delrin': -91, 'Air': -863, 'Teflon': 253, 'PMP': -399, 'LDPE': -350}
    unif_values = {'Center': -285, 'Left': -279, 'Right': -278, 'Top': -279, 'Bottom': -279}
    mtf_values = {80: 0.28, 90: 0.22, 60: 0.37, 70: 0.31, 95: 0.18}
    slice_thickness = 1


class CatPhan600_2(CatPhanMixin, TestCase):
    """An Elekta CBCT dataset"""
    catphan = CatPhan600
    file_path = ['zzCAT201602.zip']
    expected_roll = -0.64
    origin_slice = 34
    hu_values = {'Poly': -29, 'Acrylic': 123, 'Delrin': 336, 'Air': -932, 'Teflon': 897, 'PMP': -164, 'LDPE': -80}
    hu_passed = False
    unif_values = {'Center': 14, 'Left': 15, 'Right': 15, 'Top': 16, 'Bottom': 13}
    mtf_values = {50: 0.4}
    avg_line_length = 50.02
    slice_thickness = 4.5


class CatPhan604Test(CatPhanMixin, TestCase):
    catphan = CatPhan604
    file_path = ['CBCTCatPhan604.zip']
    origin_slice = 45
    hu_values = {'Poly': -47, 'Acrylic': 105, 'Delrin': 338, 'Air': -981, 'Teflon': 942, 'PMP': -194, 'LDPE': -105, '50% Bone': 771, '20% Bone': 263}
    unif_values = {'Center': -3, 'Left': 0, 'Right': 0, 'Top': 0, 'Bottom': 0}
    mtf_values = {50: 0.43}
