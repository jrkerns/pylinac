import os.path as osp
import unittest
from unittest import TestCase

import matplotlib.pyplot as plt

from pylinac import CBCT
from pylinac.core.geometry import Point
from tests.utils import save_file, LoadingTestBase, LocationMixin

VARIAN_DIR = osp.join(osp.dirname(__file__), 'test_files', 'CBCT', 'Varian')
ELEKTA_DIR = osp.join(osp.dirname(__file__), 'test_files', 'CBCT', 'Elekta')
plt.close('all')


class CBCTLoading(LoadingTestBase, TestCase):
    klass = CBCT
    constructor_input = osp.join(VARIAN_DIR, 'Pelvis')
    demo_method = 'from_demo_images'
    url = 'CBCT_4.zip'
    zip = osp.join(VARIAN_DIR, 'CBCT_4.zip')


class GeneralTests(TestCase):
    """Test general things when using cbct module."""

    def setUp(self):
        self.cbct = CBCT.from_demo_images()

    def test_demo(self):
        """Run the demo to make sure it works."""
        self.cbct.run_demo()

    def test_helpers(self):
        """Test the various helper methods."""
        self.cbct.analyze()
        self.cbct._return_results()

    def test_phan_center(self):
        """Test locations of the phantom center."""
        known_phan_center = Point(257, 255)
        self.cbct.analyze()
        self.assertAlmostEqual(self.cbct.hu.phan_center.x, known_phan_center.x, delta=0.7)
        self.assertAlmostEqual(self.cbct.hu.phan_center.y, known_phan_center.y, delta=0.7)


class PlottingSaving(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cbct = CBCT.from_demo_images()
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


class CBCTMixin(LocationMixin):
    """A mixin to use for testing Varian CBCT scans; does not inherit from TestCase as it would be run
        otherwise."""
    file_path = []
    dir_location = VARIAN_DIR
    hu_tolerance = 40
    scaling_tolerance = 1
    zip = True
    expected_roll = 0
    slice_locations = {}
    hu_values = {}
    hu_passed = True
    unif_values = {}
    unif_passed = True
    mtf_values = {}
    avg_line_length = 50
    length_passed = True
    thickness_passed = True
    lowcon_visible = 0

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        if cls.zip:
            cls.cbct = CBCT.from_zip(filename)
        else:
            cls.cbct = CBCT(filename)
        cls.cbct.analyze(cls.hu_tolerance, cls.scaling_tolerance)

    @classmethod
    def tearDownClass(cls):
        # somewhere there is a memory leak if ``cbct`` isn't deleted.
        delattr(cls, 'cbct')

    def test_slice_thickness(self):
        """Test the slice thickness."""
        self.assertEqual(self.cbct.thickness.passed, self.thickness_passed)

    def test_lowcontrast_bubbles(self):
        """Test the number of low contrast bubbles visible."""
        self.assertAlmostEqual(self.cbct.lowcontrast.rois_visible, self.lowcon_visible, delta=1)

    def test_all_passed(self):
        """Test the pass flags for all tests."""
        self.assertEqual(self.cbct.hu.overall_passed, self.hu_passed)
        self.assertEqual(self.cbct.uniformity.overall_passed, self.unif_passed)
        self.assertEqual(self.cbct.geometry.overall_passed, self.length_passed)

    def test_slice_locations(self):
        """Test the locations of the slices of interest."""
        for attr, slice_name in zip(('hu_slice_num', 'un_slice_num', 'sr_slice_num', 'lc_slice_num'), ('HU', 'UN', 'SR', 'LC')):
            self.assertAlmostEqual(getattr(self.cbct.settings, attr), self.slice_locations[slice_name], delta=1)

    def test_phantom_roll(self):
        """Test the roll of the phantom."""
        self.assertAlmostEqual(self.cbct.settings.phantom_roll, self.expected_roll, delta=0.3)

    def test_HU_values(self):
        """Test HU values."""
        for key, roi in self.cbct.hu.rois.items():
            exp_val = self.hu_values[key]
            meas_val = roi.pixel_value
            self.assertAlmostEqual(exp_val, meas_val, delta=5)

    def test_uniformity_values(self):
        """Test Uniformity HU values."""
        for key, roi in self.cbct.uniformity.rois.items():
            exp_val = self.unif_values[key]
            meas_val = roi.pixel_value
            self.assertAlmostEqual(exp_val, meas_val, delta=5)

    def test_geometry_line_length(self):
        """Test the geometry distances."""
        self.assertAlmostEqual(self.avg_line_length, self.cbct.geometry.avg_line_length, delta=0.1)

    def test_MTF_values(self):
        """Test MTF values."""
        for key, exp_mtf in self.mtf_values.items():
            meas_mtf = self.cbct.spatialres.mtf(key)
            self.assertAlmostEqual(exp_mtf, meas_mtf, delta=0.1)


class CBCTDemo(CBCTMixin, TestCase):
    """Test the CBCT demo (Varian high quality head protocol)."""
    expected_roll = -0.3
    slice_locations = {'HU': 32, 'UN': 6, 'SR': 44, 'LC': 20}
    hu_values = {'Poly': -45, 'Acrylic': 117, 'Delrin': 341, 'Air': -998, 'Teflon': 997, 'PMP': -200, 'LDPE': -103}
    unif_values = {'Center': 17, 'Left': 10, 'Right': 0, 'Top': 6, 'Bottom': 6}
    mtf_values = {80: 0.76, 90: 0.61, 60: 0.99, 70: 0.88, 95: 0.45}
    avg_line_length = 49.92
    lowcon_visible = 3

    @classmethod
    def setUpClass(cls):
        cls.cbct = CBCT.from_demo_images()
        cls.cbct.analyze()


class CBCT4(CBCTMixin, TestCase):
    """A Varian CBCT dataset"""
    file_path = ['CBCT_4.zip']
    expected_roll = -2.57
    slice_locations = {'HU': 31, 'UN': 6, 'SR': 43, 'LC': 19}
    hu_values = {'Poly': -33, 'Acrylic': 119, 'Delrin': 335, 'Air': -979, 'Teflon': 970, 'PMP': -185, 'LDPE': -94}
    unif_values = {'Center': 17, 'Left': 10, 'Right': 22, 'Top': 18, 'Bottom': 13}
    mtf_values = {80: 0.47, 90: 0.39, 60: 0.63, 70: 0.55, 95: 0.3}
    avg_line_length = 49.94
    # thickness_passed = False
    lowcon_visible = 3


class Elekta2(CBCTMixin, TestCase):
    """An Elekta CBCT dataset"""
    file_path = ['Elekta_2.zip']
    dir_location = ELEKTA_DIR
    slice_locations = {'HU': 162, 'UN': 52, 'SR': 132, 'LC': 132}
    hu_values = {'Poly': -319, 'Acrylic': -224, 'Delrin': -91, 'Air': -863, 'Teflon': 253, 'PMP': -399, 'LDPE': -350}
    hu_passed = False
    unif_values = {'Center': -285, 'Left': -279, 'Right': -278, 'Top': -279, 'Bottom': -279}
    unif_passed = False
    mtf_values = {80: 0.53, 90: 0.44, 60: 0.74, 70: 0.63, 95: 0.36}
    lowcon_visible = 2
