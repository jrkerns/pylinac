import unittest

from pylinac.picketfence import *


class PF_EPID_demo(unittest.TestCase):
    """Tests specifically for the EPID demo image."""
    test_im_path = osp.abspath(osp.join(osp.dirname(__file__), 'test_files', 'Picket Fence'))

    def setUp(self):
        self.pf = PicketFence()
        self.pf.load_demo_image()
        self.pf.analyze()

    def test_demo(self):
        self.pf.run_demo()

    def test_rotated_demo(self):
        self.pf = PicketFence()
        self.pf.load_demo_image()
        self.pf.image.rot90()
        self.pf.analyze()
        self.pf.return_results()
        self.pf.plot_analyzed_image()

    def test_demo_lower_tolerance(self):
        self.pf = PicketFence()
        self.pf.load_demo_image()
        self.pf.analyze(0.15, action_tolerance=0.05)
        self.pf.plot_analyzed_image()
        self.assertAlmostEqual(self.pf.percent_passing, 95, delta=0.5)

    def test_image_orientation(self):
        """Test image orientation."""
        # check original orientation
        self.assertEqual(self.pf.orientation, orientations['LR'])
        # check that 90 degree orientation is other way
        self.pf.image.rot90()
        self.pf._threshold()
        self.pf._find_orientation()
        self.assertEqual(self.pf.orientation, orientations['UD'])

    def test_number_pickets_found(self):
        # check number of strips found
        self.assertEqual(self.pf.num_pickets, 10, msg="MLC strips found not expected value")

    def test_error_values(self):
        # check error values
        self.assertAlmostEqual(self.pf.abs_median_error, 0.067, delta=0.02)
        self.assertAlmostEqual(self.pf.max_error, 0.213, delta=0.08)

    def test_filter_on_load(self):
        pf = PicketFence()
        pf.load_image(osp.join(self.test_im_path, 'EPID-PF.dcm'), filter=3)

    def test_passed(self):
        self.assertTrue(self.pf.passed)

    def test_analyze_tol_values(self):
        self.assertRaises(ValueError, self.pf.analyze, 0.2, 0.3)

    def test_percent_passing(self):
        self.assertEqual(self.pf.percent_passing, 100)

    def test_all_orientations_give_same_error(self):
        """Test that all orientations (0, 90, 180, 270) result in close agreement of error."""
        medians = np.zeros(4)
        maxs = np.zeros(4)

        medians[0] = self.pf.abs_median_error
        maxs[0] = self.pf.max_error

        for idx in range(1,4):
            new_pf = PicketFence()
            new_pf.load_demo_image()
            new_pf.image.rot90(idx)
            new_pf.analyze()
            medians[idx] = new_pf.abs_median_error
            maxs[idx] = new_pf.max_error
            self.assertAlmostEqual(medians[idx], medians[idx - 1], delta=0.05)
            self.assertAlmostEqual(maxs[idx], maxs[idx - 1], delta=0.1)