import json
import tempfile
import unittest
from builtins import AttributeError
from pathlib import Path
from unittest import TestCase

import numpy as np
import quaac
from parameterized import parameterized
from quaac import Attachment, Equipment, User

from pylinac import Interpolation
from pylinac.core.scale import wrap360
from pylinac.core.utilities import (
    OptionListMixin,
    convert_to_enum,
    is_close_degrees,
    is_iterable,
    simple_round,
    uniquify,
)

performer = User(name="James Kerns", email="j@j.com")
linac = Equipment(
    name="Clinac",
    model="Clinac 21EX",
    manufacturer="Varian",
    type="linac",
    serial_number="1234",
)


class TestQuaacMixin(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        from pylinac.starshot import Starshot  # local so no circular import

        cls.instance = Starshot.from_demo_image()
        cls.instance.analyze()

    def test_cant_overwrite_by_default(self):
        # write to a file first
        with tempfile.NamedTemporaryFile(delete=False) as t:
            with open(t.name, "w") as f:
                f.write("blah")
        # try to write to the same file
        with self.assertRaises(FileExistsError):
            self.instance.to_quaac(t.name, performer=performer, primary_equipment=linac)

    def test_overwrite(self):
        # write to a file first
        with tempfile.NamedTemporaryFile(delete=False) as t:
            with open(t.name, "w") as f:
                f.write("blah")
        # try to write to the same file
        self.instance.to_quaac(
            t.name, performer=performer, primary_equipment=linac, overwrite=True
        )

    def test_to_json(self):
        p = Path(tempfile.gettempdir()) / "test.json"
        if p.exists():
            p.unlink()
        self.instance.to_quaac(
            p, format="json", performer=performer, primary_equipment=linac
        )

    def test_add_attachment(self):
        with tempfile.NamedTemporaryFile(delete=False) as t:
            a = Attachment.from_file(
                path=t.name, name="test", comment="test attachment"
            )
        with tempfile.NamedTemporaryFile(delete=False) as f:
            self.instance.to_quaac(
                f.name,
                performer=performer,
                primary_equipment=linac,
                attachments=[a],
                overwrite=True,
            )
        quaac_file = quaac.Document.from_yaml_file(f.name)
        self.assertEqual(len(quaac_file.attachments), 1)


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

    def test_wrap_360_over(self):
        self.assertEqual(wrap360(-90), 270)
        self.assertEqual(wrap360(-5), 355)
        self.assertEqual(wrap360(12), 12)
        self.assertEqual(wrap360(359), 359)
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


class QuaacTestBase:
    """Base class for testing quaac generation for pylinac classes."""

    def quaac_instance(self):
        raise NotImplementedError

    def test_write_quaac(self):
        phantom = self.quaac_instance()
        with tempfile.NamedTemporaryFile(delete=False) as t:
            phantom.to_quaac(
                path=t.name,
                performer=performer,
                primary_equipment=linac,
                overwrite=True,
            )
        quaac.Document.from_yaml_file(str(t.name))


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


class TestUniquifyFunction(unittest.TestCase):
    @parameterized.expand(
        [
            (
                "test_unique_name_not_in_existing",
                ["apples", "bananas"],
                "cherries",
                "cherries",
            ),
            (
                "test_unique_name_in_existing_no_suffixes",
                ["apples", "bananas"],
                "bananas",
                "bananas-1",
            ),
            (
                "test_unique_name_in_existing_with_suffixes",
                ["apples", "bananas", "bananas-1", "bananas-2"],
                "bananas",
                "bananas-3",
            ),
            (
                "test_unique_name_with_gaps_in_suffixes",
                ["file", "file-1", "file-2", "file-4"],
                "file",
                "file-3",
            ),
            (
                "test_unique_name_large_number_of_suffixes",
                ["item"] + [f"item-{i}" for i in range(1, 1000)],
                "item",
                "item-1000",
            ),
            ("test_empty_existing_list", [], "unique", "unique"),
            ("test_empty_string_name", ["", "-1", "-2"], "", "-3"),
            (
                "test_name_already_with_suffix",
                ["report", "report-1", "report-2"],
                "report-1",
                "report-1-1",
            ),
            (
                "test_name_with_multiple_dashes",
                ["data-set", "data-set-1"],
                "data-set",
                "data-set-2",
            ),
            ("test_case_sensitivity", ["Bananas", "bananas-1"], "bananas", "bananas"),
            (
                "test_numeric_suffix_collision",
                ["item", "item-1", "item-01"],
                "item",
                "item-2",
            ),
        ]
    )
    def test_uniquify(self, name, existing, input_name, expected):
        result = uniquify(existing, input_name)
        self.assertEqual(result, expected)
