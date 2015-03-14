import unittest

from pylinac.core.decorators import *


class Test_Decorators(unittest.TestCase):

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

    def test_lazyproperty(self):

        class ExpensiveClass:
            @lazyproperty
            def expensive_property(self):
                time.sleep(0.2)
                return
        ec = ExpensiveClass()

        # run the expensive property for the first time
        start = time.time()
        _ = ec.expensive_property
        end = time.time()
        first_access_time = end - start

        # run it for the second time; should access cached property
        start = time.time()
        _ = ec.expensive_property
        end = time.time()
        cached_access_time = end - start

        self.assertLess(cached_access_time, first_access_time)

    def test_unwrap_func(self):
        #TODO: can't figure this one out
        class DumbClass:
            # a wrapped function (due to decorator that doesn't use @wraps)
            @lazyproperty
            @type_accept()
            def dumb_property(self):
                return 'result'

        d = DumbClass()
        # do something like assert dumb_property name isn't same as unwrapped name,
        # but dumb_property doesn't have __name__ attr.

