import math
from unittest import TestCase

import numpy as np

from pylinac.core.array_utils import (
    bit_invert,
    convert_to_dtype,
    filter,
    geometric_center_idx,
    geometric_center_value,
    ground,
    invert,
    normalize,
    stretch,
)


class TestGeometricCenter(TestCase):
    def test_geometric_idx_odd(self):
        arr3 = np.array((1, 3, 4))
        assert geometric_center_idx(arr3) == 1

    def test_geometric_idx_even(self):
        arr3 = np.array((1, 3, 4, 5))
        assert geometric_center_idx(arr3) == 1.5

    def test_geometric_idx_empty(self):
        arr = np.array([])
        with self.assertRaises(ValueError):
            geometric_center_idx(arr)

    def test_geometric_value_odd(self):
        arr1 = np.array((1, 2, 3))
        assert geometric_center_value(arr1) == 2

    def test_geometric_value_even(self):
        arr2 = np.array((1, 1, 2, 2))
        assert geometric_center_value(arr2) == 1.5

    def test_geometric_value_empty(self):
        arr = np.array([])
        with self.assertRaises(ValueError):
            geometric_center_value(arr)

    def test_geometric_idx_multidim_fails(self):
        arr1 = np.random.randn(2, 2)
        with self.assertRaises(ValueError):
            geometric_center_idx(arr1)

    def test_geometric_value_multidim_fails(self):
        arr1 = np.random.randn(2, 2)
        with self.assertRaises(ValueError):
            geometric_center_value(arr1)


class TestNormalize(TestCase):
    def test_normalize_to_max(self):
        arr1 = np.array((1, 2, 3, 4))
        norm_arr = normalize(arr1)
        assert norm_arr.max() == 1.0
        assert norm_arr[-1] == 1.0
        assert math.isclose(norm_arr[0], 0.25)

    def test_normalize_to_value(self):
        arr1 = np.array((1, 2, 3, 4), dtype=float)
        norm_arr = normalize(arr1, 2)
        assert math.isclose(norm_arr.max(), 2.0)
        assert math.isclose(norm_arr[1], 1)


class TestInvert(TestCase):
    def test_invert(self):
        arr1 = np.array([0, 10])
        inv_arr = invert(arr1)
        assert np.array_equal(inv_arr, [10, 0])

    def test_invert_negative(self):
        arr1 = np.array([-5, -1])
        inv_arry = invert(arr1)
        assert np.array_equal(inv_arry, [-1, -5])


class TestGround(TestCase):
    def test_ground(self):
        arr = np.array([3, 4, 5])
        grounded_arr = ground(arr)
        assert np.array_equal(grounded_arr, np.array([0, 1, 2]))

    def test_ground_negative(self):
        arr = np.array([-3, -4, -5])
        grounded_arr = ground(arr)  # moves min from -5->0
        assert np.array_equal(grounded_arr, np.array([2, 1, 0]))

    def test_ground_value(self):
        arr = np.array([3, 4, 5])
        grounded_arr = ground(arr, value=10)
        assert np.array_equal(grounded_arr, np.array([10, 11, 12]))


class TestFilter(TestCase):
    def test_median_filter_size1(self):
        arr = np.array([0, 0, 0, 3, 0, 0, 0])
        f_arr = filter(arr, size=1, kind="median")  # size=1 is a no-op
        assert np.array_equal(f_arr, np.array([0, 0, 0, 3, 0, 0, 0]))

    def test_median_filter_float(self):
        arr = np.array([0, 0, 0, 3, 0, 0, 0])
        f_arr = filter(
            arr, size=0.1, kind="median"
        )  # this float will resolve to size=1
        assert np.array_equal(f_arr, [0, 0, 0, 3, 0, 0, 0])

    def test_median_filter_size3(self):
        arr = np.array([0, 0, 0, 3, 0, 0, 0])
        f_arr = filter(arr, size=3, kind="median")
        # median with x surrounded by 0's will be 0
        assert np.array_equal(f_arr, [0, 0, 0, 0, 0, 0, 0])

        arr = np.array([0, 0, 3, 3, 0, 0, 0])
        f_arr = filter(arr, size=3, kind="median")
        # with two adjacent values, the median will be the value
        assert np.array_equal(f_arr, [0, 0, 3, 3, 0, 0, 0])

    def test_gaussian_filter_size1(self):
        arr = np.array([0, 0, 0, 3, 0, 0, 0])
        f_arr = filter(arr, size=1, kind="gaussian")
        # even a gaussian w/ size 1 will smooth out a peak like this
        assert np.array_equal(f_arr, [0, 0, 0, 1, 0, 0, 0])

    def test_invalid_float(self):
        arr = np.array([0, 0, 0, 3, 0, 0, 0])
        with self.assertRaises(ValueError) as exc:
            filter(arr, size=2.3, kind="gaussian")
            self.assertIn("was not between", str(exc.exception))

    def test_invalid_filter(self):
        arr = np.array([0, 0, 0, 3, 0, 0, 0])
        with self.assertRaises(ValueError) as exc:
            filter(arr, size=1, kind="filterthis")
            self.assertIn("unsupported", str(exc.exception))


class TestComplement(TestCase):
    def test_complement_uint8(self):
        arr1 = np.array([0, 10], dtype=np.uint8)
        inv_arr = bit_invert(arr1)
        assert np.array_equal(inv_arr, [255, 245])

    def test_complement_uint16(self):
        arr1 = np.array([0, 10], dtype=np.uint16)
        inv_arr = bit_invert(arr1)
        assert np.array_equal(inv_arr, [65535, 65525])

    def test_complement_int8(self):
        # int8 has a sign, so this will actually revolve about 0 minus 1
        arr1 = np.array([0, 10], dtype=np.int8)
        inv_arr = bit_invert(arr1)
        assert np.array_equal(inv_arr, [-1, -11])

    def test_complement_float_fails(self):
        arr1 = np.array([0, 10], dtype=float)
        with self.assertRaises(ValueError) as ve:
            bit_invert(arr1)
            self.assertIn("could not be safely", str(ve.exception))


class TestStretch(TestCase):
    def test_simple_stretch(self):
        """Default is equivalent to ground/normalize"""
        arr1 = np.array([5, 6, 7])
        s_arr = stretch(arr1)
        assert np.array_equal(s_arr, [0, 0.5, 1])

    def test_stretch_new_max(self):
        arr1 = np.array([5, 6, 7])
        s_arr = stretch(arr1, max=10)
        assert np.array_equal(s_arr, [0, 5, 10])

    def test_stretch_new_min(self):
        arr1 = np.array([5, 6, 7])
        s_arr = stretch(arr1, min=8, max=10)
        assert np.array_equal(s_arr, [8, 9, 10])

    def test_stretch_a_lot(self):
        arr1 = np.array([5, 20, 30])
        s_arr = stretch(arr1, min=8, max=50)
        assert np.array_equal(s_arr, [8, 33.2, 50])

    def test_bad_max(self):
        arr1 = np.array([5, 6, 7])
        with self.assertRaises(ValueError):
            stretch(arr1, min=8)  # default max is 1

    def test_max_outside_dtype(self):
        arr1 = np.array([5, 6, 7], dtype=np.uint8)
        with self.assertRaises(ValueError) as ve:
            stretch(arr1, max=30000)
            self.assertIn("was larger", str(ve.exception))

    def test_min_outside_dtype(self):
        arr1 = np.array([5, 6, 7], dtype=np.uint8)
        with self.assertRaises(ValueError) as ve:
            stretch(arr1, min=-1)
            self.assertIn("was smaller", str(ve.exception))


class TestConvertDataType(TestCase):
    def test_upward_conversion(self):
        arr1 = np.array([5, 6, 7], dtype=np.uint8)
        c_arr = convert_to_dtype(arr1, dtype=np.uint16)
        # going from uint8 to uint16 will increase by 65535/255 = 257x
        assert np.array_equal(c_arr, [1285, 1542, 1799])
        self.assertEqual(c_arr.dtype, np.uint16)

    def test_downward_conversion(self):
        arr1 = np.array([0, 100, 1000, 10000, 65535], dtype=np.uint16)
        c_arr = convert_to_dtype(arr1, dtype=np.uint8)
        # going from uint16 to uint8 will decrease by 255/65535 = 1/257
        assert np.array_equal(c_arr, [0, 1, 4, 39, 255])
        self.assertEqual(c_arr.dtype, np.uint8)

    def test_unsigned_to_signed_conversion(self):
        arr1 = np.array([0, 255], dtype=np.uint8)
        c_arr = convert_to_dtype(arr1, dtype=np.int8)
        # values will be negative because we are near the bottom of the range of values
        assert np.array_equal(c_arr, [-128, 127])
        self.assertEqual(c_arr.dtype, np.int8)

    def test_float_becomes_normalized(self):
        arr1 = np.array([0, 255.2], dtype=float)
        c_arr = convert_to_dtype(arr1, dtype=np.uint16)
        assert np.array_equal(c_arr, [0, 65535])
        self.assertEqual(c_arr.dtype, np.uint16)
