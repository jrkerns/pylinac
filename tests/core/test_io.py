"""Test suite for the pylinac.io module."""

import unittest
import os.path as osp

from pylinac.core.io import is_valid_dir, is_valid_file, open_file

class Test_IO(unittest.TestCase):

    def test_invalid_file_str(self):
        not_a_file = "file"
        self.assertFalse(is_valid_file(not_a_file, raise_error=False))

        self.assertRaises(FileExistsError, is_valid_file, not_a_file)

    def test_valid_file(self):
        real_file = __file__
        self.assertTrue(is_valid_file(real_file))

        real_file_obj = open(__file__)
        self.assertTrue(is_valid_file(real_file_obj))

    def test_is_valid_dir(self):
        not_a_dir = "dir"
        self.assertRaises(NotADirectoryError, is_valid_dir, not_a_dir)
        self.assertFalse(is_valid_dir(not_a_dir, False))

        is_a_dir = osp.dirname(__file__)
        self.assertTrue(is_valid_dir(is_a_dir))

    def test_open_file_str(self):
        not_a_file = "file"
        self.assertRaises(FileNotFoundError, open_file, not_a_file)

        real_file = __file__
        f = open_file(real_file)  # shouldn't raise

    def test_open_file_obj(self):
        real_file_obj = open(__file__)
        f = open_file(real_file_obj)  # shouldn't raise

        f2 = open_file(f)  # also shouldn't raise on previously-read file

    def test_cursor_reset(self):
        """Test that opening a file object if it's already open resets the cursor to 0."""
        real_file_obj = open(__file__)
        self.assertEqual(real_file_obj.tell(), 0)
        f = open_file(real_file_obj)  # shouldn't raise
        self.assertEqual(f.tell(), 0)

        f2 = open_file(f)  # also shouldn't raise on previously-read file
        self.assertEqual(f2.tell(), 0)
