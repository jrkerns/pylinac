from __future__ import division, absolute_import, print_function

import unittest

from pylinac.cbctqa.cbct import CBCT


# Current pass/fail rates depend on tolerances. E.g. within 5 HU/%/mm might be better
class CBCT_end2end_tests(unittest.TestCase):
    """This is an end-to-end test of the CBCT module.
        The demo files are loaded and analyzed. The results should match what is expected. This is the most general
        test to run; similar to an RPC phantom test. If the end results are what we expect we're good, so goes
        the philosophy.
        We will test the 3 demo image packs included with pylinac.
        """

    def setUp(self):
        self.cbct = CBCT()

    def test_hihead(self):
        """End-to-end test of High-quality head."""
        self.cbct.load_demo_images('hi_head')
        self.analyze_pass()

    def test_lothorax(self):
        """End-to-end test of low dose thorax."""
        self.cbct.load_demo_images('thorax')
        self.analyze_pass()

    def test_pelvis(self):
        """End-to-end test of High-quality head."""
        self.cbct.load_demo_images('pelvis')
        self.analyze_pass()

    def analyze_pass(self):
        """Analyze the CBCT and assert that all parameters pass."""
        self.cbct.analyze()
        self.assertTrue((self.cbct.HU_passed, self.cbct.GEO_passed, self.cbct.UNIF_passed))
