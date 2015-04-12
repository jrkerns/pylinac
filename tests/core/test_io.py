"""Test suite for the pylinac.io module."""

import unittest
import os.path as osp

from pylinac.core.io import is_valid_dir, is_valid_file

class Test_IO(unittest.TestCase):

    def test_is_valid_file(self):
        not_a_file = "file"
        self.assertFalse(is_valid_file(not_a_file, raise_error=False))

        self.assertRaises(FileExistsError, is_valid_file, not_a_file)

        real_file = __file__
        self.assertTrue(is_valid_file(real_file))

    def test_is_valid_dir(self):
        not_a_dir = "dir"
        self.assertRaises(NotADirectoryError, is_valid_dir, not_a_dir)
        self.assertFalse(is_valid_dir(not_a_dir, False))

        is_a_dir = osp.dirname(__file__)
        self.assertTrue(is_valid_dir(is_a_dir))
