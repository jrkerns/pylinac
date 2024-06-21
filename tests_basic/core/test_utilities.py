import json
import unittest
from builtins import AttributeError
from unittest import TestCase

import numpy as np

from pylinac import Interpolation
from pylinac.core.scale import abs360, wrap360
from pylinac.core.utilities import (
    OptionListMixin,
    convert_to_enum,
    is_close_degrees,
    is_iterable,
    simple_round,
)


class TestIsClose(unittest.TestCase):
    def test_is_close_degrees_simple(self):
        self.assertTrue(is_close_degrees(0, 0))
        self.assertTrue(is_close_degrees(0, 360))
        self.assertTrue(is_close_degrees(360, 0))
        self.assertTrue(is_close_degrees(0, 1))
        self.assertTrue(is_close_degrees(0, 359))
        self.assertTrue(is_close_degrees(73, 73.2))
        self.assertTrue(is_close_degrees(-90, 270))

    def test_is_close_degrees_not_close(self):
        self.assertFalse(is_close_degrees(0, 180))
        self.assertFalse(is_close_degrees(0, 90))
        self.assertFalse(is_close_degrees(0, 270))
        self.assertFalse(is_close_degrees(90, -90))

    def test_delta(self):
        self.assertTrue(is_close_degrees(0, 4, delta=5))
        self.assertTrue(is_close_degrees(0, 356, delta=5))
        self.assertTrue(is_close_degrees(95, 110, delta=15))
        with self.assertRaises(ValueError):
            is_close_degrees(0, 1, delta=-1)


class TestUtilities(unittest.TestCase):
    def test_is_iterable(self):
        # test iterables
        iters = ((1, 2, "t"), [4, 8, "r"], np.array((5, 6, 7)))
        for iter in iters:
            self.assertTrue(is_iterable(iter))
        # test non-iterables
        noniters = (5,)
        for iter in noniters:
            self.assertFalse(is_iterable(iter))

    def test_convert_to_enum(self):
        self.assertEqual(Interpolation.LINEAR, convert_to_enum("Linear", Interpolation))
        self.assertEqual(
            Interpolation.LINEAR, convert_to_enum(Interpolation.LINEAR, Interpolation)
        )
        with self.assertRaises(ValueError):
            convert_to_enum("baffled", Interpolation)

    def test_absolute_360(self):
        self.assertEqual(abs360(-90), 270)
        self.assertEqual(abs360(-5), 355)
        self.assertEqual(abs360(12), 12)
        self.assertEqual(abs360(359), 359)

    def test_wrap_360_over(self):
        self.assertEqual(wrap360(361), 1)
        self.assertEqual(wrap360(360), 0)
        self.assertEqual(wrap360(359.6), 359.6)
        self.assertEqual(wrap360(359), 359)
        self.assertEqual(wrap360(180), 180)


class TestSimpleRound(unittest.TestCase):
    def test_precision(self):
        self.assertEqual(simple_round(0.12345, decimals=0), 0.0)
        self.assertEqual(simple_round(0.12345, decimals=1), 0.1)
        self.assertEqual(simple_round(0.12345, decimals=2), 0.12)
        self.assertEqual(simple_round(0.12345, decimals=3), 0.123)
        self.assertEqual(simple_round(0.12345, decimals=4), 0.1234)
        self.assertEqual(simple_round(0.12345, decimals=None), 0.12345)

    def test_0decimals_returns_int(self):
        self.assertIsInstance(simple_round(0.12345, decimals=0), int)

    def test_1_or_more_decimals_returns_float(self):
        self.assertIsInstance(simple_round(0.12345, decimals=1), float)
        self.assertIsInstance(simple_round(12, decimals=1), float)
        self.assertIsInstance(simple_round(12, decimals=2), float)

    def test_int_in_is_int_out(self):
        self.assertIsInstance(simple_round(12, decimals=None), int)


class TestOptionMixin(TestCase):
    def test_option_list(self):
        class MyOptions(OptionListMixin):
            APPLES = "aPpLes"
            ORANGES = "Oranges"

        self.assertIsInstance(MyOptions.options(), list)
        self.assertEqual(len(MyOptions.options()), 2)
        self.assertEqual(MyOptions.APPLES, "aPpLes")
        self.assertListEqual(MyOptions.options(), ["aPpLes", "Oranges"])

    def test_option_with_method(self):
        class MyOptions(OptionListMixin):
            APPLES = "aPpLes"
            ORANGES = "Oranges"

            def kill_me_now(self):
                pass

        self.assertIsInstance(MyOptions.options(), list)
        self.assertEqual(len(MyOptions.options()), 2)
        self.assertEqual(MyOptions.APPLES, "aPpLes")
        self.assertListEqual(MyOptions.options(), ["aPpLes", "Oranges"])


class ResultsDataBase:
    model: callable

    def construct_analyzed_instance(self):
        try:
            instance = self.model.from_demo_image()
        except AttributeError:
            instance = self.model.from_demo_images()
        instance.analyze()
        return instance

    def test_results(self):
        instance = self.construct_analyzed_instance()
        instance.results_data()
        data_dict = instance.results_data(as_dict=True)
        # ensure it's json serializable; shouldn't raise
        json.dumps(data_dict)

        data_json = instance.results_data(as_json=True)
        self.assertIsInstance(data_json, str)
        json.loads(data_json)
