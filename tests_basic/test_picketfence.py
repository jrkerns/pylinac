import os.path as osp
import tempfile
from unittest import TestCase, skip

import matplotlib.pyplot as plt

from pylinac.picketfence import PicketFence, UP_DOWN, LEFT_RIGHT
from tests_basic.utils import save_file, LoadingTestBase, LocationMixin

TEST_DIR = osp.join(osp.dirname(__file__), 'test_files', 'Picket Fence')


class TestLoading(LoadingTestBase, TestCase):
    klass = PicketFence
    constructor_input = osp.join(TEST_DIR, 'AS500_PF.dcm')
    url = 'EPID-PF-LR.dcm'

    def test_filter_on_load(self):
        PicketFence(self.constructor_input, filter=3)  # shouldn't raise

    def test_load_with_log(self):
        log_file = osp.join(TEST_DIR, 'PF_log.bin')
        pf_file = osp.join(TEST_DIR, 'PF.dcm')
        pf = PicketFence(pf_file, log=log_file)
        pf.analyze(hdmlc=True)


class GeneralTests(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.pf = PicketFence.from_demo_image()
        cls.pf.analyze()

    def test_bad_tolerance_values(self):
        self.assertRaises(ValueError, self.pf.analyze, 0.2, 0.3)

    def test_demo(self):
        PicketFence.run_demo()

    def test_publish_pdf(self):
        with tempfile.TemporaryFile() as t:
            self.pf.publish_pdf(t, notes='stuff', metadata={'Unit': 'TB1'})


class TestPlottingSaving(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.pf = PicketFence.from_demo_image()
        cls.pf.analyze()
        cls.pf_updown = PicketFence.from_demo_image()
        cls.pf_updown.image.rot90()
        cls.pf_updown.analyze()

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def test_plotting(self):
        self.pf.plot_analyzed_image()
        self.pf_updown.plot_analyzed_image()

    def test_saving_image(self):
        save_file(self.pf.save_analyzed_image)
        save_file(self.pf_updown.save_analyzed_image)


class PFTestMixin(LocationMixin):
    """Base Mixin for testing a picketfence image."""
    dir_location = TEST_DIR
    picket_orientation = UP_DOWN
    hdmlc = False
    num_pickets = 10
    pass_num_pickets = False
    percent_passing = 100
    max_error = 0
    abs_median_error = 0
    sag_adjustment = 0
    passes = True
    log = None
    mean_picket_spacing = 15

    @classmethod
    def get_logfile(cls):
        """Return the canonical path to the log file."""
        if cls.log is not None:
            return osp.join(cls.dir_location, *cls.log)

    @classmethod
    def setUpClass(cls):
        cls.pf = PicketFence(cls.get_filename(), log=cls.get_logfile())
        if cls.pass_num_pickets:
            cls.pf.analyze(hdmlc=cls.hdmlc, sag_adjustment=cls.sag_adjustment, num_pickets=cls.num_pickets)
        else:
            cls.pf.analyze(hdmlc=cls.hdmlc, sag_adjustment=cls.sag_adjustment)

    def test_passed(self):
        self.assertEqual(self.pf.passed, self.passes)

    def test_picket_orientation(self):
        self.assertEqual(self.pf.orientation, self.picket_orientation)

    def test_num_pickets(self):
        self.assertEqual(self.pf.num_pickets, self.num_pickets)

    def test_percent_passing(self):
        self.assertAlmostEqual(self.pf.percent_passing, self.percent_passing, delta=1)

    def test_max_error(self):
        self.assertAlmostEqual(self.pf.max_error, self.max_error, delta=0.1)

    def test_abs_median_error(self):
        self.assertAlmostEqual(self.pf.abs_median_error, self.abs_median_error, delta=0.05)

    def test_picket_spacing(self):
        self.assertAlmostEqual(self.pf.pickets.mean_spacing, self.mean_picket_spacing, delta=0.5)


class PFDemo(PFTestMixin, TestCase):
    """Tests specifically for the EPID demo image."""
    picket_orientation = LEFT_RIGHT
    max_error = 0.217
    abs_median_error = 0.06

    @classmethod
    def setUpClass(cls):
        cls.pf = PicketFence.from_demo_image()
        cls.pf.analyze(hdmlc=cls.hdmlc, sag_adjustment=cls.sag_adjustment)

    def test_demo_lower_tolerance(self):
        pf = PicketFence.from_demo_image()
        pf.analyze(0.15, action_tolerance=0.05)
        pf.plot_analyzed_image()
        self.assertAlmostEqual(pf.percent_passing, 95, delta=1)


class AS1200(PFTestMixin, TestCase):
    """Tests for the AS1200 image."""
    file_path = ['AS1200.dcm']
    max_error = 0.08
    abs_median_error = 0.02


class ClinacWeirdBackground(PFTestMixin, TestCase):
    file_path = ['Clinac weird background.dcm']
    max_error = 0.12
    abs_median_error = 0.02
    num_pickets = 5
    mean_picket_spacing = 50


class ElektaCloseEdges(PFTestMixin, TestCase):
    file_path = ['PF, Elekta, pickets near edges.dcm']
    max_error = 0.23
    abs_median_error = 0.07
    num_pickets = 9
    mean_picket_spacing = 30


class ElektaCloseEdgesRot90(PFTestMixin, TestCase):
    file_path = ['PF, Elekta, pickets near edges.dcm']
    max_error = 0.23
    abs_median_error = 0.07
    num_pickets = 9
    mean_picket_spacing = 30
    picket_orientation = LEFT_RIGHT

    @classmethod
    def setUpClass(cls):
        cls.pf = PicketFence(cls.get_filename(), log=cls.get_logfile())
        cls.pf.image.rot90()
        cls.pf.analyze(hdmlc=cls.hdmlc, sag_adjustment=cls.sag_adjustment)


class MultipleImagesPF(PFTestMixin, TestCase):
    """Test of a multiple image picket fence; e.g. EPID images."""
    max_error = 0.112
    abs_median_error = 0.019
    picket_orientation = LEFT_RIGHT
    num_pickets = 5
    mean_picket_spacing = 30

    @classmethod
    def setUpClass(cls):
        path1 = osp.join(TEST_DIR, 'combo-jaw.dcm')
        path2 = osp.join(TEST_DIR, 'combo-mlc.dcm')
        cls.pf = PicketFence.from_multiple_images([path1, path2], stretch_each=True)
        cls.pf.analyze(hdmlc=cls.hdmlc, sag_adjustment=cls.sag_adjustment, orientation='left')
