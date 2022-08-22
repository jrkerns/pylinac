import unittest

from pylinac import Interpolation
from pylinac.core.utilities import *


class TestUtilities(unittest.TestCase):

    def test_is_iterable(self):
        # test iterables
        iters = ((1,2,'t'), [4, 8, 'r'], np.array((5,6,7)))
        for iter in iters:
            self.assertTrue(is_iterable(iter))
        # test non-iterables
        noniters = (5,)
        for iter in noniters:
            self.assertFalse(is_iterable(iter))

    def test_convert_to_enum(self):
        self.assertEqual(Interpolation.LINEAR, convert_to_enum('Linear', Interpolation))
        self.assertEqual(Interpolation.LINEAR, convert_to_enum(Interpolation.LINEAR, Interpolation))
        with self.assertRaises(ValueError):
            convert_to_enum('baffled', Interpolation)
