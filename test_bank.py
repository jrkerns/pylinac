"""Test the entire package; an underscore precedes this file name
so it does not include itself in the test discovery."""
import os.path as osp
import unittest

import matplotlib


matplotlib.use('Agg')


def run_tests(directory, pattern='test*.py'):
    # import a runner to run tests_basic
    runner = unittest.TextTestRunner()
    # run test discovery
    test_suite = unittest.defaultTestLoader.discover(directory, pattern)
    # run test runner
    runner.run(test_suite)


test_dir = osp.join(osp.dirname(__file__), 'tests_bank')
run_tests(test_dir)
