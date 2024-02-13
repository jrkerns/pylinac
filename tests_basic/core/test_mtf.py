import unittest

import numpy as np

from pylinac.core.mtf import MTF


class TestMTF(unittest.TestCase):
    def test_normal_mtf(self):
        pair_units = (0.1, 0.2, 0.3)
        maxs = (500, 300, 100)
        mins = (25, 50, 75)

        m = MTF(pair_units, maxs, mins)
        rm = m.relative_resolution(x=50)
        self.assertAlmostEqual(rm, 0.24, delta=0.03)
        rm = m.relative_resolution(x=90)
        self.assertAlmostEqual(rm, 0.15, delta=0.03)

    def test_mtf_lower_than_values(self):
        # should generate a warning
        pair_units = (0.1, 0.2, 0.3)
        maxs = (500, 300, 100)
        mins = (25, 50, 75)

        m = MTF(pair_units, maxs, mins)
        rm = m.relative_resolution(x=10)
        self.assertAlmostEqual(rm, 0.3, delta=0.03)

    def test_non_decreasing_mtf(self):
        # this will return the first occurrence where the condition is met.
        # should generate a warning
        pair_units = (0.1, 0.2, 0.3, 0.4)
        maxs = (500, 300, 500, 100)
        mins = (25, 50, 25, 75)

        MTF(pair_units, maxs, mins)

    def test_no_zero_division_of_line_pair_distances(self):
        old_settings = np.geterr()
        # set overflow to cause errors
        # shouldn't raise
        np.seterr(all="raise")
        pair_units = (0.1, 0.2, 0.3)
        maxs = (500, 300, 100)
        mins = (25, 50, 75)

        m = MTF(pair_units, maxs, mins)
        # shouldn't raise
        m.plot()
        # reset to old settings
        np.seterr(**old_settings)
