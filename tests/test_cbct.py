import unittest

from pylinac.cbct import CBCT


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


