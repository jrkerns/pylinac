from unittest import TestCase

from pylinac.picketfence import PicketFence, osp, np, UP_DOWN, LEFT_RIGHT

from tests.utils import save_file, LoadingTestBase

TEST_DIR = osp.join(osp.dirname(__file__), 'test_files', 'Picket Fence')


class TestLoading(LoadingTestBase, TestCase):
    klass = PicketFence
    constructor_input = osp.join(osp.dirname(osp.dirname(__file__)), 'pylinac', 'demo_files', 'picket_fence',
                                 'EPID-PF-LR.dcm')
    url = 'AS500-UD.dcm'

    def test_filter_on_load(self):
        PicketFence(self.constructor_input, filter=3)  # shouldn't raise


class GeneralTests(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.pf = PicketFence.from_demo_image()
        cls.pf.analyze()

    def test_bad_tolerance_values(self):
        self.assertRaises(ValueError, self.pf.analyze, 0.2, 0.3)


class TestPlottingSaving(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.pf = PicketFence.from_demo_image()
        cls.pf.analyze()

    def test_plotting(self):
        self.pf.plot_analyzed_image()

    def test_saving_image(self):
        save_file(self.pf.save_analyzed_image)
        save_file(self.pf.save_analyzed_image, interactive=True)


class PFTestMixin:
    """Base Mixin for testing a picketfence image."""
    im_path = ''
    dir_location = TEST_DIR
    picket_orientation = UP_DOWN
    hdmlc = False
    num_pickets = 10
    percent_passing = 100
    max_error = 0
    abs_median_error = 0
    sag_adjustment = 0
    passes = True

    @classmethod
    def setUpClass(cls):
        cls.pf = PicketFence(cls.get_filename())
        cls.pf.analyze(hdmlc=cls.hdmlc, sag_adjustment=cls.sag_adjustment)

    @classmethod
    def get_filename(cls):
        return osp.join(cls.dir_location, cls.im_path)

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

    def test_all_orientations(self):
        median_errors = []
        for rotation in range(4):
            self.pf.image.rot90()
            self.pf.analyze(hdmlc=self.hdmlc, sag_adjustment=self.sag_adjustment)
            median_errors.append(self.pf.abs_median_error)

        for error in median_errors:
            self.assertAlmostEqual(error, np.mean(median_errors), delta=0.1)


class PFDemo(PFTestMixin, TestCase):
    """Tests specifically for the EPID demo image."""
    picket_orientation = LEFT_RIGHT
    max_error = 0.217
    abs_median_error = 0.06

    @classmethod
    def get_filename(cls):
        return osp.join(osp.dirname(osp.dirname(__file__)), 'pylinac', 'demo_files', 'picket_fence', 'EPID-PF-LR.dcm')

    def test_demo(self):
        self.pf.run_demo()

    def test_demo_lower_tolerance(self):
        pf = PicketFence.from_demo_image()
        pf.analyze(0.15, action_tolerance=0.05)
        pf.plot_analyzed_image()
        self.assertAlmostEqual(pf.percent_passing, 94, delta=1)


class MultipleImagesPF(PFTestMixin, TestCase):
    """Test of a multiple image picket fence; e.g. EPID images."""
    max_error = 0.112
    abs_median_error = 0.019
    picket_orientation = LEFT_RIGHT
    num_pickets = 5

    @classmethod
    def setUpClass(cls):
        path1 = osp.join(TEST_DIR, 'combo-jaw.dcm')
        path2 = osp.join(TEST_DIR, 'combo-mlc.dcm')
        cls.pf = PicketFence.from_multiple_images([path1, path2])
        cls.pf.analyze(hdmlc=cls.hdmlc, sag_adjustment=cls.sag_adjustment)
