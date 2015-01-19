# Test the entire package; an underscore precedes this file name
# so it does not include itself in the test discovery.

import unittest
import os.path as osp


def full_test_suite_discovery():
    """Return a test suite that discovers test throughout the entire package."""
    # discover tests recursively starting from the top level module
    package_dir = osp.dirname(__file__)
    # suite = unittest.defaultTestLoader.loadTestsFromModule(osp.join(test_dir, 'test_cbct'))
    # loader.loadTestsFromModule(osp.join(test_dir, 'test_cbct'))
    suite = unittest.TestLoader().discover(package_dir)
    return suite


# ---------------------------
# Actual testing
# ---------------------------
# import a runner to run tests
runner = unittest.TextTestRunner()
# run test discovery
test_suite = full_test_suite_discovery()
# run test runner
runner.run(test_suite)
