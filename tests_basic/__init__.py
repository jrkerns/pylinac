import os.path as osp
import unittest

import matplotlib

TEST_FILES_DIR = osp.abspath(osp.join(osp.dirname(__file__), 'test_files'))
TEST_BANK_DIR = osp.abspath(osp.join(osp.dirname(__file__), '..', '..', 'pylinac test files'))

matplotlib.use('Agg')


def run_tests(directory, pattern='test*.py'):
    # import a runner to run tests_basic
    runner = unittest.TextTestRunner()
    # run test discovery
    test_suite = unittest.defaultTestLoader.discover(directory, pattern)
    # run test runner
    runner.run(test_suite)
