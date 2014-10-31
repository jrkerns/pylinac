# Test the entire package

import unittest


def full_test_suite_discovery():
    """Return a test suite that discovers test throughout the entire package."""
    # discover tests recursively starting from the top level module
    suite = unittest.TestLoader().discover('pylinac')
    return suite


# ---------------------------
# Actual testing
# ---------------------------
# import a runner to run tests
runner = unittest.TextTestRunner()
# run test discovery
test_suite = full_test_suite_discovery()
# run it
runner.run(test_suite)