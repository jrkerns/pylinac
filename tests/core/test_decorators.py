import unittest

from pylinac.core.decorators import *


class TestDecorators(unittest.TestCase):

    def test_timethis(self):
        @timethis
        def dumb_function():
            return [item for item in range(1000000)]

        dumb_function()

    def test_type_accept(self):
        @type_accept(param1=(int, float), param2=(str, None))
        def dumb_function(param1, param2):
            pass
        @type_accept(param1=int)
        def dumb_func2(param1):
            pass

        # valid parameter types
        df = dumb_function(1, 'r')
        df = dumb_function(2.1, None)

        # invalid parameter types
        self.assertRaises(TypeError, dumb_function, 'r', 't')
        self.assertRaises(TypeError, dumb_function, 1, 2)
        self.assertRaises(TypeError, dumb_func2, 2.2)

    def test_value_accept(self):

        @value_accept(param1=(1,5), param2=('left', 'right'))
        def dumb_function(param1, param2):
            pass

        # valid parameter values
        df = dumb_function(2, 'right')
        df = dumb_function(4.999, 'left')

        # invalid parameter values
        self.assertRaises(ValueError, dumb_function, 2, 'light')
        self.assertRaises(ValueError, dumb_function, 6, 'right')
        self.assertRaises(ValueError, dumb_function, 9, 'bright')

    def test_value_and_type_accept(self):

        @type_accept(param1=(int, float))
        @value_accept(param1=(1,5))
        def dumb_function(param1):
            pass

        # valid use
        df = dumb_function(3)
        df = dumb_function(1.01)

        # invalid use
        self.assertRaises(ValueError, dumb_function, 6)  # correct type, wrong value
        self.assertRaises(TypeError, dumb_function, '3')  # incorrect type
