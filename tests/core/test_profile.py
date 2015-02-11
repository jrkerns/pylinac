import unittest

import scipy.signal as sps

from pylinac.core.profile import *


class Test_Profile(unittest.TestCase):
    # create reference profiles
    # ...create periodic xdata to calculate over; periodic because of sawtooth func
    xdata = np.arange(8 * np.pi, step=0.05 * np.pi)
    # ...create periodic sawtooth profile
    values = sps.sawtooth(xdata, width=0.7)

    def test_inputs(self):
        p = Profile()
        self.assertIsNone(p.y_values)

        # test that xdata is generated if not passed in
        p = Profile(self.values)
        xvals = np.arange(len(p.y_values))
        self.assertCountEqual(p.x_values, xvals)

        # test that x and y values are numpy arrays
        self.assertIsInstance(p.y_values, np.ndarray)
        self.assertIsInstance(p.x_values, np.ndarray)

    def test_ground_profile(self):
        """Test that the profile is properly grounded to 0."""
        p = Profile(self.values)
        # the minimum shouldn't be zero to start with
        self.assertFalse(p.y_values.min() == 0)
        # but it should be after grounding
        p.ground()
        self.assertTrue(p.y_values.min() == 0)

    def test_find_peaks(self):
        p = Profile(self.values)
        p.find_peaks()
        known_peak_locs = (28, 68, 108, 148)
        for found_peak, known_peak in zip(p.peaks, known_peak_locs):
            self.assertEqual(found_peak.idx, known_peak)

    def test_find_valleys(self):
        p = Profile(self.values)
        p.find_valleys()
        known_valley_locs = (40, 80, 120)
        for found_valley, known_valley in zip(p.valleys, known_valley_locs):
            self.assertEqual(found_valley.idx, known_valley)

    def test_find_FWHM_peaks(self):
        p = Profile(self.values)
        p.find_FWXM_peaks()
        known_peak_locs = (24, 64, 103, 143)
        for found_peak, known_peak in zip(p.peaks, known_peak_locs):
            self.assertAlmostEqual(found_peak.idx, known_peak, delta=1)

    def test_subdivide_profiles(self):
        p = Profile(self.values)
        p.find_peaks()
        subprofiles = p._subdivide_profiles()

        # test subprofile is a singleprofile instance
        self.assertIsInstance(subprofiles[0], SingleProfile)

        # test that x-values got propagated
        # I.e. profile1.x_values = [0,1,2] profile2.x_values = [3,4,5], ...
        known_end_points = (0, 28, 68, 108)
        found_left_ends = [subprofile.x_values[0] for subprofile in subprofiles]
        for found_left_end, known_left_end in zip(found_left_ends, known_end_points):
            self.assertEqual(found_left_end, known_left_end)


class Test_CircleProfile(unittest.TestCase):

    def test_inputs(self):
        cp = CircleProfile()
        self.assertEqual(cp.center.y, 0)
        self.assertIsNone(cp.y_values)

        center = Point(50, 50)
        radius = 30
        cp = CircleProfile(center, radius)
        self.assertEqual(cp.center.x, center.x)

    def test_get_profile(self):
        img_array = np.ones((101, 101))
        center = Point(50, 50)
        radius = 30

        cp = CircleProfile(center, radius)
        cp.get_profile(img_array)
        # test random sample is 1 (the whole matrix is 1's)
        self.assertEqual(cp.y_values[500], 1)
        # test x_locs set properly
        self.assertAlmostEqual(cp.x_locs[0], center.x+radius, delta=0.01)

        small_img_array = np.ones((20,20))
        cp = CircleProfile(center, radius)
        self.assertRaises(ValueError, cp.get_profile, small_img_array)


class Test_SingleProfile(unittest.TestCase):
    # create reference profiles
    # ...create periodic xdata to calculate over; periodic because of sawtooth func
    xdata = np.arange(1.9 * np.pi, step=0.02 * np.pi)
    # ...create periodic sawtooth profile
    values = sps.sawtooth(xdata, width=0.7)
    initial_peak = np.pi

    penums = (20, 40, 50, 70, 80, 90)

    def test_inputs(self):
        p = SingleProfile(self.values)
        # test that xdata is generated if not passed in
        xvals = np.arange(len(p.y_values))
        self.assertCountEqual(p.x_values, xvals)
        self.assertEqual(p.x_values[-1], len(xvals)-1)

        # test that the initial peak found is reasonable
        p = SingleProfile(self.values, self.xdata)
        known_peak = 4
        self.assertAlmostEqual(p.initial_peak, known_peak, delta=2)

        # test that a bad peak passed raises an error
        bad_peak_idx = 20
        self.assertRaises(IndexError, SingleProfile, self.values, self.xdata, initial_peak=bad_peak_idx)

        # test that bad profile raises error
        xdata = np.arange(1.5 * np.pi, step=0.02 * np.pi)
        bad_profile = sps.sawtooth(xdata, width=0.7)
        self.assertRaises(ValueError, SingleProfile, bad_profile)

    def test_get_X_penum_idx(self):
        p = SingleProfile(self.values)
        # test various penumbra points

        left_known_idxs = (13, 27, 34, 48, 55, 62)
        for penum, known_idx in zip(self.penums, left_known_idxs):
            self.assertEqual(p.get_X_penum_idx('left', penum), known_idx)

        right_known_idxs = (89, 84, 82, 77, 75, 73)
        for penum, known_idx in zip(self.penums, right_known_idxs):
            self.assertEqual(p.get_X_penum_idx('right', penum), known_idx)

    def test_getFWXM(self):
        p = SingleProfile(self.values)
        known_widths = (76, 57, 48, 29, 20, 11)
        for penum, known_width in zip(self.penums, known_widths):
            self.assertEqual(p.get_FWXM(penum), known_width)

    def test_get_FWXM_center(self):
        p = SingleProfile(self.values)
        known_idxs = (51, 55.5, 58, 62.5, 65, 67.5)
        for penum, known_width in zip(self.penums, known_idxs):
            self.assertAlmostEqual(p.get_FWXM_center(penum), known_width, delta=0.1)

            # test that values are rounded when asked for
        known_idxs = (51, 56, 58, 62, 65, 68)
        for penum, known_width in zip(self.penums, known_idxs):
            self.assertEqual(p.get_FWXM_center(penum, True), known_width)

    def test_get_penum_width(self):
        p = SingleProfile(self.values)
        # test that error raised when lower > upper
        self.assertRaises(ValueError, p.get_penum_width, upper=20, lower=80)

        sides = ('left', 'right', 'left')
        lower_penums = (10, 20, 30)
        known_lower_penums = (49, 14, 35)
        for side, penum, known_penum in zip(sides, lower_penums, known_lower_penums):
            self.assertEqual(p.get_penum_width(side, lower=penum), known_penum)

        upper_penums = (70, 80, 90)
        known_upper_penums = (35, 14, 49)
        for side, penum, known_penum in zip(sides, upper_penums, known_upper_penums):
            self.assertEqual(p.get_penum_width(side, upper=penum), known_penum)

    def test_get_field_value(self):
        p = SingleProfile(self.values)

        field_widths = (0.3, 0.6, 0.8, 0.9)
        values = ('mean', 'median', 'max', 'min')
        known_values = (0.63, 0.63, 1.0, 0.03)
        for width, value, known_value in zip(field_widths, values, known_values):
            self.assertAlmostEqual(p.get_field_calculation(width, value), known_value, delta=0.03)
