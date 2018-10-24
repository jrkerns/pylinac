"""Test the entire package; an underscore precedes this file name
so it does not include itself in the test discovery."""
import os.path as osp
import unittest

# from tests import run_tests


def run_tests(directory, pattern='test*.py'):
    # import a runner to run tests
    runner = unittest.TextTestRunner()
    # run test discovery
    test_suite = unittest.defaultTestLoader.discover(directory, pattern)
    # run test runner
    runner.run(test_suite)

test_dir = osp.dirname(__file__)
run_tests(test_dir)
