import unittest

from pylinac.core.decorators import *


class TestDecorators(unittest.TestCase):

    def test_timethis(self):
        @timethis
        def dumb_function():
            return [item for item in range(1000000)]

        dumb_function()
