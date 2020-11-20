import unittest
import os.path as osp

from pylinac.core.utilities import *


class TestUtilities(unittest.TestCase):

    def test_isnumeric(self):
        # test numerics
        numerics = (5, 3.2, np.array((5,6))[0])
        for number in numerics:
            self.assertTrue(isnumeric(number))

        # test non-numerics
        notnumerics = ('5', np.array((5,6)))
        for notnumeric in notnumerics:
            self.assertFalse(isnumeric(notnumeric))

    def test_is_iterable(self):
        # test iterables
        iters = ((1,2,'t'), [4, 8, 'r'], np.array((5,6,7)))
        for iter in iters:
            self.assertTrue(is_iterable(iter))
        # test non-iterables
        noniters = (5,)
        for iter in noniters:
            self.assertFalse(is_iterable(iter))
