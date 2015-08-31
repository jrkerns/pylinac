import unittest
import time
import os

from pylinac.picketfence import PicketFence, osp, np

test_file_dir = osp.join(osp.dirname(__file__), 'test_files', 'Picket Fence')


class PFTestMixin:
    """Base Mixin for testing a picketfence image."""
    im_path = ''
    picket_orientation = 'Up-Down'
    hdmlc = False
    num_pickets = 10
    percent_passing = 100
    max_error = 0
    abs_median_error = 0
    sag_adjustment = 0
    passes = True

    @classmethod
    def setUpClass(cls):
        cls.pf = PicketFence(cls.im_path)
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

    def test_all_orientations(self):
        median_errors = []
        for rotation in range(4):
            self.pf.image.rot90()
            self.pf.analyze(hdmlc=self.hdmlc, sag_adjustment=self.sag_adjustment)
            median_errors.append(self.pf.abs_median_error)

        for error in median_errors:
            self.assertAlmostEqual(error, np.mean(median_errors), delta=0.1)

    def test_plotting(self):
        self.pf.plot_analyzed_image()

    def test_saving_image(self):
        filename = 'tester.png'
        self.pf.save_analyzed_image(filename)

        time.sleep(0.1)  # sleep just to let OS work
        self.assertTrue(osp.isfile(filename), "Save file did not successfully save the image")

        # cleanup
        os.remove(filename)
        self.assertFalse(osp.isfile(filename), "Save file test did not clean up saved image")


class PFDemo(PFTestMixin, unittest.TestCase):
    """Tests specifically for the EPID demo image."""
    im_path = osp.join(osp.dirname(osp.dirname(__file__)), 'pylinac', 'demo_files', 'picket_fence', 'EPID-PF-LR.dcm')
    picket_orientation = 'Left-Right'
    max_error = 0.217
    abs_median_error = 0.06

    def test_demo(self):
        self.pf.run_demo()

    def test_demo_lower_tolerance(self):
        pf = PicketFence.from_demo_image()
        pf.analyze(0.15, action_tolerance=0.05)
        pf.plot_analyzed_image()
        self.assertAlmostEqual(pf.percent_passing, 93, delta=1)


class GeneralTests(unittest.TestCase):

    def test_filter_on_load(self):
        PicketFence(osp.join(osp.dirname(osp.dirname(__file__)), 'pylinac', 'demo_files', 'picket_fence',
                    'EPID-PF-LR.dcm'), filter=3)

    def test_bad_tolerance_values(self):
        pf = PicketFence.from_demo_image()
        self.assertRaises(ValueError, pf.analyze, 0.2, 0.3)

    def test_from_url(self):
        """Test getting a PF image from a URL."""
        url = 'https://s3.amazonaws.com/assuranceqa-staging/uploads/imgs/AS500-UD.dcm'
        pf = PicketFence.from_url(url)
        pf.analyze()

        bad_url = 'https://s3.amazonaws.com/assuranceqa-staging/uploads/imgs/AS500-UD_not_real.dcm'
        with self.assertRaises(ConnectionError):
            pf = PicketFence.from_url(bad_url)
