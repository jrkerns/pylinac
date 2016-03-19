"""Test the entire package; an underscore precedes this file name
so it does not include itself in the test discovery."""
import unittest
import os.path as osp


def run(directory, pattern='test*.py'):
    # import a runner to run tests
    runner = unittest.TextTestRunner()
    # run test discovery
    test_suite = unittest.defaultTestLoader.discover(directory, pattern)
    # run test runner
    runner.run(test_suite)

if __name__ == "__main__":
    test_dir = osp.dirname(__file__)
    run(test_dir)
