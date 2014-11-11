from __future__ import division, absolute_import, print_function

import unittest

from pylinac.cbctqa.cbct import CBCT


# Current pass/fail rates depend on tolerances. E.g. within 5 HU/%/mm might be better
class CBCT_demos(unittest.TestCase):
    """This is an end-to-end test of the CBCT module.
        The demo files are loaded and analyzed. The results should match what is expected. This is the most general
        test to run; similar to an RPC phantom test. If the end results are what we expect we're good, so goes
        the philosophy.
        We will test the 3 demo image packs included with pylinac.
        """

    def setUp(self):
        self.cbct = CBCT()

    def test_head_all_pass(self):
        """Test that all parameters (HU, Unif, scaling, SR) for the High-Quality Head demo passes."""
        self.cbct.load_demo_images('hi_head')
        self.analyze_pass()

    def test_thorax_all_pass(self):
        """Test that all parameters (HU, Unif, scaling, SR) for the High-Quality Head demo passes."""
        self.cbct.load_demo_images('thorax')
        self.analyze_pass()

    def test_pelvis_all_pass(self):
        """Test that all parameters (HU, Unif, scaling, SR) for the High-Quality Head demo passes."""
        self.cbct.load_demo_images('pelvis')
        self.analyze_pass()

    def analyze_pass(self):
        """Analyze the CBCT and assert that all parameters pass."""
        self.cbct.analyze()
        self.assertTrue((self.cbct.HU.passed_test, self.cbct.GEO.passed_test, self.cbct.UN.passed_test))
