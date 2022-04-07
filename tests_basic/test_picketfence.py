import io
import os
import os.path as osp
import tempfile
from unittest import TestCase, skip

import matplotlib.pyplot as plt

from pylinac.core import image
from pylinac.picketfence import PicketFence, Orientation, PFResult, MLCArrangement
from tests_basic.utils import save_file, CloudFileMixin, get_file_from_cloud_test_repo, InitTesterMixin, \
    FromURLTesterMixin, FromDemoImageTesterMixin

TEST_DIR = 'picket_fence'


class TestInstantiation(TestCase, InitTesterMixin, FromURLTesterMixin, FromDemoImageTesterMixin):
    klass = PicketFence
    init_file = ['picket_fence', 'AS500_PF.dcm']
    url = 'EPID-PF-LR.dcm'

    def test_filter_on_load(self):
        PicketFence(self.full_init_file, filter=3)  # shouldn't raise

    def test_load_with_log(self):
        log_file = get_file_from_cloud_test_repo([TEST_DIR, 'PF_log.bin'])
        pf_file = get_file_from_cloud_test_repo([TEST_DIR, 'PF.dcm'])
        pf = PicketFence(pf_file, log=log_file)
        pf.analyze()

    def test_load_from_file_object(self):
        pf_file = get_file_from_cloud_test_repo([TEST_DIR, 'PF.dcm'])
        ref_pf = PicketFence(pf_file)
        ref_pf.analyze()
        with open(pf_file, "rb") as f:
            pf = PicketFence(f)
            pf.analyze()
        self.assertIsInstance(pf, PicketFence)
        self.assertEqual(pf.percent_passing, ref_pf.percent_passing)

    def test_load_from_stream(self):
        pf_file = get_file_from_cloud_test_repo([TEST_DIR, 'PF.dcm'])
        ref_pf = PicketFence(pf_file)
        ref_pf.analyze()
        with open(pf_file, "rb") as f:
            s = io.BytesIO(f.read())
            pf = PicketFence(s)
            pf.analyze()
        self.assertIsInstance(pf, PicketFence)
        self.assertEqual(pf.percent_passing, ref_pf.percent_passing)

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

    def test_image_kwargs(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, 'AS500_PF.dcm'])

        # do normal analysis
        phan = PicketFence(path)
        phan.analyze()
        offset = phan.results_data().offsets_from_cax_mm[0]

        # pass kwarg; use same dpi as image; CAX offset should be the same, would be different with different DPI
        img = image.load(path)
        phan = PicketFence(path, image_kwargs={'dpi': img.dpi})
        phan.analyze()
        offset_manual_dpi = phan.results_data().offsets_from_cax_mm[0]

        self.assertEqual(offset, offset_manual_dpi)


class TestAnalyze(TestCase):

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
            self.pf.publish_pdf(t.name, notes="stuff", metadata={"Unit": "TB1"})
        os.remove(t.name)

    def test_results_data(self):
        data = self.pf.results_data()
        self.assertIsInstance(data, PFResult)
        self.assertEqual(data.max_error_mm, self.pf.max_error)

        data_dict = self.pf.results_data(as_dict=True)
        self.assertIsInstance(data_dict, dict)
        self.assertIn("pylinac_version", data_dict)

    def test_no_measurements_suggests_inversion(self):
        file_loc = get_file_from_cloud_test_repo([TEST_DIR, 'noisy-FFF-wide-gap-pf.dcm'])
        pf = PicketFence(file_loc)
        pf.image.invert()
        with self.assertRaises(ValueError):
            pf.analyze(invert=False)

    def test_orientation_passing_as(self):
        # below shouldn't raise
        # as enum
        pf = PicketFence.from_demo_image()
        pf.analyze(orientation=Orientation.UP_DOWN)

        # as str
        pf = PicketFence.from_demo_image()
        pf.analyze(orientation="Up-Down")

    def test_histogram(self):
        pf = PicketFence.from_demo_image()
        pf.analyze()
        pf.plot_histogram()

        pf2 = PicketFence.from_demo_image()
        # can't plot before analyzing
        with self.assertRaises(ValueError):
            pf2.plot_histogram()

    def test_failed_leaves_before_analyzed(self):
        pf = PicketFence.from_demo_image()
        with self.assertRaises(ValueError):
            pf.failed_leaves()

    def test_failed_leaves_traditional(self):
        pf = PicketFence.from_demo_image()
        pf.analyze(separate_leaves=False, tolerance=0.05)
        self.assertEqual(
            set(pf.failed_leaves()),
            {12, 14, 16, 17, 18, 19, 26, 29, 31, 32, 33, 35, 39, 42, 43, 47},
        )

    def test_failed_leaves_separate(self):
        pf = PicketFence.from_demo_image()
        pf.analyze(separate_leaves=True, tolerance=0.15, nominal_gap_mm=3)
        self.assertEqual(
            set(pf.failed_leaves()),
            {
                "A12",
                "B12",
                "A13",
                "B13",
                "A14",
                "B14",
                "A15",
                "B15",
                "B16",
                "B17",
                "A18",
                "B18",
                "B19",
                "A19",
                "A20",
                "B20",
                "B22",
                "A23",
                "B23",
                "A24",
                "B24",
                "A25",
                "B26",
                "A27",
                "B27",
                "A28",
                "B28",
                "B29",
                "B30",
                "B31",
                "A32",
                "B32",
                "A33",
                "B34",
                "B35",
                "A35",
                "A36",
                "B36",
                "A37",
                "B37",
                "A38",
                "B38",
                "B39",
                "A39",
                "A40",
                "B40",
                "B41",
                "A41",
                "B42",
                "A43",
                "B43",
                "B44",
                "A44",
                "B45",
                "A45",
                "A46",
                "B46",
                "B47",
                "A47",
            },
        )


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
        plt.close("all")

    def test_plotting(self):
        self.pf.plot_analyzed_image()
        self.pf_updown.plot_analyzed_image()

    def test_saving_image(self):
        save_file(self.pf.save_analyzed_image)
        save_file(self.pf_updown.save_analyzed_image)

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

    def test_plot_histogram(self):
        self.pf.plot_histogram()

    def test_save_histogram(self):
        # to disk
        save_file(self.pf.save_histogram)
        # to binary stream
        save_file(self.pf.save_histogram, as_file_object='b')

    def test_plot_leaf_profile(self):
        self.pf.plot_leaf_profile(20, 3)

    def test_save_leaf_profile(self):
        # to disk
        save_file(self.pf.save_leaf_profile, 20, 3)
        # to binary stream
        save_file(self.pf.save_leaf_profile, 20, 3, as_file_object='b')


class PFTestMixin(CloudFileMixin):
    """Base Mixin for testing a picketfence image."""
    dir_path = ['picket_fence']
    picket_orientation = Orientation.UP_DOWN
    mlc = "Millennium"
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
    separate_leaves = False
    nominal_gap_mm = 1
    mlc_skew = 0
    max_error_picket = None
    max_error_leaf = None

    @classmethod
    def get_logfile(cls):
        """Return the canonical path to the log file."""
        if cls.log is not None:
            return osp.join(*cls.dir_path, *cls.log)

    @classmethod
    def setUpClass(cls):
        cls.pf = PicketFence(cls.get_filename(), log=cls.get_logfile())
        if cls.pass_num_pickets:
            cls.pf.analyze(
                sag_adjustment=cls.sag_adjustment,
                num_pickets=cls.num_pickets,
                invert=cls.invert,
                separate_leaves=cls.separate_leaves,
                nominal_gap_mm=cls.nominal_gap_mm,
            )
        else:
            cls.pf.analyze(
                sag_adjustment=cls.sag_adjustment,
                invert=cls.invert,
                separate_leaves=cls.separate_leaves,
                nominal_gap_mm=cls.nominal_gap_mm,
            )

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
        self.assertAlmostEqual(
            self.pf.abs_median_error, self.abs_median_error, delta=0.05
        )

    def test_picket_spacing(self):
        self.assertAlmostEqual(
            self.pf.mean_picket_spacing, self.mean_picket_spacing, delta=0.5
        )

    def test_mlc_skew(self):
        self.assertAlmostEqual(self.pf.mlc_skew(), self.mlc_skew, delta=0.3)

    def test_max_picket(self):
        if self.max_error_picket:
            self.assertEqual(self.pf.max_error_picket, self.max_error_picket)

    def test_max_leaf(self):
        if self.max_error_leaf:
            self.assertIn(self.max_error_leaf, self.pf.max_error_leaf)


class PFDemo(PFTestMixin, TestCase):
    """Tests specifically for the EPID demo image."""

    picket_orientation = Orientation.UP_DOWN
    max_error = 0.08
    abs_median_error = 0.06
    max_error_picket = 0
    max_error_leaf = 31

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
        self.assertAlmostEqual(pf.percent_passing, 100, delta=1)


class WideGapSimulation(PFTestMixin, TestCase):
    file_name = 'noisy-wide-gap-pf.dcm'
    max_error = 0.11
    abs_median_error = 0.06
    num_pickets = 7
    mean_picket_spacing = 30


class WideGapSimulationSeparate(WideGapSimulation):
    separate_leaves = True
    nominal_gap_mm = 16
    max_error = 0.3
    abs_median_error = 0.07
    percent_passing = 100


class FFFWideGapSimulation(PFTestMixin, TestCase):
    file_name = 'noisy-FFF-wide-gap-pf.dcm'
    max_error = 0.17
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
    invert = True


class ElektaCloseEdges(PFTestMixin, TestCase):
    file_name = 'PF,-Elekta,-pickets-near-edges.dcm'
    max_error = 0.23
    abs_median_error = 0.07
    num_pickets = 9
    mean_picket_spacing = 30
    mlc_skew = -0.7


class ElektaCloseEdgesRot90(PFTestMixin, TestCase):
    file_name = 'PF,-Elekta,-pickets-near-edges.dcm'
    max_error = 0.23
    abs_median_error = 0.07
    num_pickets = 9
    mean_picket_spacing = 30
    picket_orientation = Orientation.LEFT_RIGHT
    mlc_skew = 0.7

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
    delete_file = False

    @classmethod
    def setUpClass(cls):
        path1 = get_file_from_cloud_test_repo([TEST_DIR, 'combo-jaw.dcm'])
        path2 = get_file_from_cloud_test_repo([TEST_DIR, 'combo-mlc.dcm'])
        cls.pf = PicketFence.from_multiple_images([path1, path2], stretch_each=True)
        cls.pf.analyze(
            sag_adjustment=cls.sag_adjustment,
            orientation=Orientation.LEFT_RIGHT
        )


class AS500(PFTestMixin, TestCase):
    """Tests for the AS500 image."""
    file_name = 'AS500_PF.dcm'
    max_error = 0.15
    abs_median_error = 0.04


class AS5002(PFTestMixin, TestCase):
    """Tests for the AS500#2 image."""
    file_name = 'AS500#2.dcm'
    max_error = 0.12
    abs_median_error = 0.03
    mlc_skew = -0.3


class AS5003(PFTestMixin, TestCase):
    """Tests for the AS500#3 image."""
    file_name = 'AS500#3.dcm'
    max_error = 0.16
    abs_median_error = 0.03


class AS5004(PFTestMixin, TestCase):
    """Tests for the AS500#4 image."""
    file_name = 'AS500#4.dcm'
    max_error = 0.28
    abs_median_error = 0.06
    mlc_skew = -0.3


class AS5005(PFTestMixin, TestCase):
    """Tests for the AS500#4 image."""
    file_name = 'AS500#5.dcm'
    max_error = 0.23
    abs_median_error = 0.04


class AS5006(PFTestMixin, TestCase):
    """Tests for the AS500#4 image."""
    file_name = 'AS500#6.dcm'
    picket_orientation = Orientation.LEFT_RIGHT
    max_error = 0.23
    abs_median_error = 0.06


class AS5007(PFTestMixin, TestCase):
    """Tests for the AS500#4 image."""
    file_name = 'AS500#7.dcm'
    max_error = 0.24
    abs_median_error = 0.05
    mlc_skew = -0.3


class AS5008(PFTestMixin, TestCase):
    """Tests for the AS500#4 image."""
    file_name = 'AS500#8.dcm'
    max_error = 0.2
    abs_median_error = 0.04
    mlc_skew = -0.3


class AS5009(PFTestMixin, TestCase):
    """Tests for the AS500#4 image."""
    file_name = 'AS500#9.dcm'
    max_error = 0.24
    abs_median_error = 0.04
    mlc_skew = -0.3


class AS50010(PFTestMixin, TestCase):
    """Tests for the AS500#4 image."""
    file_name = 'AS500#10.dcm'
    picket_orientation = Orientation.LEFT_RIGHT
    max_error = 0.24
    abs_median_error = 0.05


class AS500error(PFTestMixin, TestCase):
    """Tests for the AS500#2 image."""
    file_name = 'AS500-error.dcm'
    num_pickets = 6
    percent_passing = 99
    max_error = 0.55
    abs_median_error = 0.07
    passes = False
    mean_picket_spacing = 20
    mlc_skew = -0.3


class AS1000(PFTestMixin, TestCase):
    """Tests for the AS1000 image."""
    file_name = 'AS1000_PF.dcm'
    max_error = 0.29
    abs_median_error = 0.06


class AS1000_2(PFTestMixin, TestCase):
    """Tests for the AS1000 image."""
    file_name = 'AS1000#2.dcm'
    max_error = 0.24
    abs_median_error = 0.07


class AS1000_3(PFTestMixin, TestCase):
    """Tests for the AS1000 image."""
    file_name = 'AS1000#3.dcm'
    max_error = 0.13
    abs_median_error = 0.05


class AS1000_4(PFTestMixin, TestCase):
    """Tests for the AS1000 image."""
    file_name = 'AS1000#4.dcm'
    picket_orientation = Orientation.LEFT_RIGHT
    max_error = 0.18
    abs_median_error = 0.05


class AS1000_90(PFTestMixin, TestCase):
    """Tests for the AS1000 image."""
    file_name = 'AS1000-90.dcm'
    picket_orientation = Orientation.LEFT_RIGHT
    max_error = 0.23
    abs_median_error = 0.05


class AS1000HDSmall(PFTestMixin, TestCase):
    """Tests for the AS1000 image."""
    file_name = 'AS1000-HD-small.dcm'
    mlc = 'HD'
    max_error = 0.05
    abs_median_error = 0.05


class AS1000HDFull(PFTestMixin, TestCase):
    """Tests for the AS1000 image with a smaller pattern (only inner leaves)."""
    file_name = 'AS1000-HD-full.dcm'
    mlc = 'HD'
    max_error = 0.2
    abs_median_error = 0.06


class AS1000HDFullVMAT(PFTestMixin, TestCase):
    """Tests for the AS1000 image with a smaller pattern (only inner leaves)."""
    file_name = 'AS1000-HD-full-VMAT.dcm'
    mlc = 'HD'
    max_error = 0.2
    abs_median_error = 0.08


@skip  # says file isn't real DICOM TODO: Figure out why not real DICOM
class AS1000HDFullError(PFTestMixin, TestCase):
    """Tests for the AS1000 image with a few errors introduced."""
    file_name = 'AS1000-HD-full-error.dcm'
    mlc = 'HD'
    num_pickets = 6
    abs_median_error = 0.03
    max_error = 0.39

    def test_lower_tolerance_fails(self):
        """This image has an introduced error; this should catch with a reasonable tolerance."""
        pf = PicketFence(self.file_path)
        pf.analyze(tolerance=0.3, hdmlc=self.hdmlc)
        self.assertFalse(pf.passed)


class AS1200(PFTestMixin, TestCase):
    """Tests for the AS1200 image."""
    file_name = 'AS1200.dcm'
    max_error = 0.08
    abs_median_error = 0.02


class AS1200Error(PFTestMixin, TestCase):
    """Tests for the AS1200 image."""
    file_name = 'AS1200-error.dcm'
    num_pickets = 6
    max_error = 0.48
    abs_median_error = 0.05
    sag_adjustment = -1.2
    mean_picket_spacing = 20


class AS1200ExtendedSID(PFTestMixin, TestCase):
    """Tests for the AS1200 image."""
    file_name = 'AS1200-ExtendedSID.dcm'
    max_error = 0.12
    abs_median_error = 0.04


class AS1200ExtendedSIDVMAT(PFTestMixin, TestCase):
    """Tests for the AS1200 image."""
    file_name = 'AS1200-ExtendedSID-VMAT.dcm'
    max_error = 0.18
    abs_median_error = 0.06


# @expectedFailure  # too dirty
# class AS1200HD(PFTestMixin, TestCase):
#     """Tests for the AS1200 image."""
#     file_name = 'AS1200-HD.dcm'
#     mlc = 'HD'
#     max_error = 0.05
#     abs_median_error = 0.02
#     num_pickets = 10
#     pass_num_pickets = True


# @expectedFailure  # terribly dirty image with artifacts all over.
# class AS1200HDTranslated(PFTestMixin, TestCase):
#     """Tests for the AS1200 image."""
#     file_name = 'AS1200-HD-translated.dcm'
#     mlc = 'HD'
#     max_error = 0.15
#     abs_median_error = 0.02
#     num_pickets = 10
#     pass_num_pickets = True


class ChicagoNoError(PFTestMixin, TestCase):
    dir_path = [TEST_DIR, 'Chicago']
    file_name = 'PF no error.dcm'
    # log = ['Chicago', 'PF no error tlog.bin']
    mlc = 'HD'
    max_error = 0.24


class ChicagoError(PFTestMixin, TestCase):
    dir_path = [TEST_DIR, 'Chicago']
    file_name = 'PF point2mm error.dcm'
    # log = ['Chicago', 'PF point2mm tlog.bin']
    mlc = 'HD'
    max_error = 0.3


@skip
class CharlestonRA(PFTestMixin, TestCase):
    file_name = ['Charleston', 'TB1', 'July2016', 'RA.dcm']
    max_error = 0.17


@skip
class CharlestonG0(PFTestMixin, TestCase):
    file_name = ['Charleston', 'TB1', 'July2016', 'G0.dcm']
    max_error = 0.1
