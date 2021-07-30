import io
import os
import os.path as osp
import tempfile
from unittest import TestCase

import matplotlib.pyplot as plt

from pylinac.picketfence import PicketFence, Orientation, PFResult, MLCArrangement
from tests_basic.utils import save_file, LoadingTestBase, CloudFileMixin, get_file_from_cloud_test_repo

TEST_DIR = 'picket_fence'


class TestLoading(LoadingTestBase, TestCase):
    klass = PicketFence
    constructor_input = ['picket_fence', 'AS500_PF.dcm']
    url = 'EPID-PF-LR.dcm'

    def test_filter_on_load(self):
        PicketFence(self.get_constructor_input(), filter=3)  # shouldn't raise

    def test_load_with_log(self):
        log_file = get_file_from_cloud_test_repo([TEST_DIR, 'PF_log.bin'])
        pf_file = get_file_from_cloud_test_repo([TEST_DIR, 'PF.dcm'])
        pf = PicketFence(pf_file, log=log_file)
        pf.analyze()

    def test_load_from_file_object(self):
        pf_file = get_file_from_cloud_test_repo([TEST_DIR, 'PF.dcm'])
        ref_pf = PicketFence(pf_file)
        ref_pf.analyze()
        with open(pf_file, 'rb') as f:
            pf = PicketFence(f)
            pf.analyze()
        self.assertIsInstance(pf, PicketFence)
        self.assertEqual(pf.percent_passing, ref_pf.percent_passing)

    def test_load_from_stream(self):
        pf_file = get_file_from_cloud_test_repo([TEST_DIR, 'PF.dcm'])
        ref_pf = PicketFence(pf_file)
        ref_pf.analyze()
        with open(pf_file, 'rb') as f:
            s = io.BytesIO(f.read())
            pf = PicketFence(s)
            pf.analyze()
        self.assertIsInstance(pf, PicketFence)
        self.assertEqual(pf.percent_passing, ref_pf.percent_passing)


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
        with tempfile.NamedTemporaryFile(delete=False) as t:
            self.pf.publish_pdf(t.name, notes='stuff', metadata={'Unit': 'TB1'})
        os.remove(t.name)

    def test_results_data(self):
        data = self.pf.results_data()
        self.assertIsInstance(data, PFResult)
        self.assertEqual(data.max_error_mm, self.pf.max_error)

        data_dict = self.pf.results_data(as_dict=True)
        self.assertIsInstance(data_dict, dict)
        self.assertIn('pylinac_version', data_dict)

    def test_no_measurements_suggests_inversion(self):
        file_loc = get_file_from_cloud_test_repo([TEST_DIR, 'noisy-FFF-wide-gap-pf.dcm'])
        pf = PicketFence(file_loc)
        with self.assertRaises(ValueError):
            pf.analyze(invert=False)

    def test_orientation_passing_as(self):
        # below shouldn't raise
        # as enum
        pf = PicketFence.from_demo_image()
        pf.analyze(orientation=Orientation.LEFT_RIGHT)

        # as str
        pf = PicketFence.from_demo_image()
        pf.analyze(orientation="Left-Right")

    def test_custom_MLC_arrangement(self):
        mlc_setup = MLCArrangement(leaf_arrangement=[(10, 10), (40, 5), (10, 10)])

        # pass it in to the mlc parameter
        path = get_file_from_cloud_test_repo([TEST_DIR, 'AS500_PF.dcm'])
        pf = PicketFence(path, mlc=mlc_setup)

        # shouldn't raise
        pf.analyze()
        pf.results()
        pf.results_data()

    def test_mlc_string(self):
        mlc_setup = 'Millennium'

        # pass it in to the mlc parameter
        path = get_file_from_cloud_test_repo([TEST_DIR, 'AS500_PF.dcm'])
        pf = PicketFence(path, mlc=mlc_setup)

        # shouldn't raise
        pf.analyze()
        pf.results()
        pf.results_data()


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


class PFTestMixin(CloudFileMixin):
    """Base Mixin for testing a picketfence image."""
    dir_path = ['picket_fence']
    picket_orientation = Orientation.UP_DOWN
    mlc = 'Millennium'
    num_pickets = 10
    pass_num_pickets = False
    percent_passing = 100
    max_error = 0
    abs_median_error = 0
    sag_adjustment = 0
    invert = False
    passes = True
    log = None
    mean_picket_spacing = 15

    @classmethod
    def get_logfile(cls):
        """Return the canonical path to the log file."""
        if cls.log is not None:
            return osp.join(*cls.dir_path, *cls.log)

    @classmethod
    def setUpClass(cls):
        cls.pf = PicketFence(cls.get_filename(), log=cls.get_logfile())
        if cls.pass_num_pickets:
            cls.pf.analyze(sag_adjustment=cls.sag_adjustment, num_pickets=cls.num_pickets, invert=cls.invert)
        else:
            cls.pf.analyze(sag_adjustment=cls.sag_adjustment, invert=cls.invert)

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
        self.assertAlmostEqual(self.pf.mean_picket_spacing, self.mean_picket_spacing, delta=0.5)


class PFDemo(PFTestMixin, TestCase):
    """Tests specifically for the EPID demo image."""
    picket_orientation = Orientation.LEFT_RIGHT
    max_error = 0.217
    abs_median_error = 0.06

    @classmethod
    def setUpClass(cls):
        cls.pf = PicketFence.from_demo_image()
        cls.pf.analyze(sag_adjustment=cls.sag_adjustment)

    @classmethod
    def tearDownClass(cls):
        pass  # override delete behavior

    def test_demo_lower_tolerance(self):
        pf = PicketFence.from_demo_image()
        pf.analyze(0.15, action_tolerance=0.05)
        pf.plot_analyzed_image()
        self.assertAlmostEqual(pf.percent_passing, 94, delta=1)


class WideGapSimulation(PFTestMixin, TestCase):
    file_name = 'noisy-wide-gap-pf.dcm'
    max_error = 0.11
    invert = True
    abs_median_error = 0.06
    num_pickets = 7
    mean_picket_spacing = 30


class FFFWideGapSimulation(PFTestMixin, TestCase):
    file_name = 'noisy-FFF-wide-gap-pf.dcm'
    max_error = 0.17
    invert = True
    abs_median_error = 0.06
    num_pickets = 7
    mean_picket_spacing = 30


class AS1200(PFTestMixin, TestCase):
    """Tests for the AS1200 image."""
    file_name = 'AS1200.dcm'
    max_error = 0.08
    abs_median_error = 0.02


class ClinacWeirdBackground(PFTestMixin, TestCase):
    file_name = 'Clinac-weird-background.dcm'
    max_error = 0.12
    abs_median_error = 0.02
    num_pickets = 5
    mean_picket_spacing = 50


class ElektaCloseEdges(PFTestMixin, TestCase):
    file_name = 'PF,-Elekta,-pickets-near-edges.dcm'
    max_error = 0.23
    abs_median_error = 0.07
    num_pickets = 9
    mean_picket_spacing = 30


class ElektaCloseEdgesRot90(PFTestMixin, TestCase):
    file_name = 'PF,-Elekta,-pickets-near-edges.dcm'
    max_error = 0.23
    abs_median_error = 0.07
    num_pickets = 9
    mean_picket_spacing = 30
    picket_orientation = Orientation.LEFT_RIGHT

    @classmethod
    def setUpClass(cls):
        cls.pf = PicketFence(cls.get_filename(), log=cls.get_logfile())
        cls.pf.image.rot90()
        cls.pf.analyze(sag_adjustment=cls.sag_adjustment)


class MultipleImagesPF(PFTestMixin, TestCase):
    """Test of a multiple image picket fence; e.g. EPID images."""
    max_error = 0.112
    abs_median_error = 0.019
    picket_orientation = Orientation.LEFT_RIGHT
    num_pickets = 5
    mean_picket_spacing = 30

    @classmethod
    def setUpClass(cls):
        path1 = get_file_from_cloud_test_repo([TEST_DIR, 'combo-jaw.dcm'])
        path2 = get_file_from_cloud_test_repo([TEST_DIR, 'combo-mlc.dcm'])
        cls.pf = PicketFence.from_multiple_images([path1, path2], stretch_each=True)
        cls.pf.analyze(sag_adjustment=cls.sag_adjustment, orientation=Orientation.LEFT_RIGHT)
