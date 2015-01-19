"""Test suite for the pylinac.io module."""

import unittest

from pylinac.core.io import *

class Test_IO(unittest.TestCase):

    def test_is_valid_file(self):
        not_a_file = "file"
        self.assertFalse(is_valid_file(not_a_file))

        self.assertRaises(FileExistsError, is_valid_file, not_a_file, True)

        real_file = __file__
        self.assertTrue(is_valid_file(real_file))
