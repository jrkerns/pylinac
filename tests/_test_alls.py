"""Test the '_test_all_* test files. This is a much more complete test suite,
but accesses the larger, private test data."""
import os.path as osp

from tests._test_all import run


test_dir = osp.dirname(__file__)
run(test_dir, pattern='_test_all_*.py')
