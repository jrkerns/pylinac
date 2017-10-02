import unittest
import os.path as osp

from pylinac.core.utilities import *


class TestUtilities(unittest.TestCase):

    def test_isnumeric(self):
        # test numerics
        numerics = (5, 3.2, np.array((5,6))[0])
        for number in numerics:
            self.assertTrue(isnumeric(number))

        # test non-numerics
        notnumerics = ('5', np.array((5,6)))
        for notnumeric in notnumerics:
            self.assertFalse(isnumeric(notnumeric))

    def test_is_iterable(self):
        # test iterables
        iters = ((1,2,'t'), [4, 8, 'r'], np.array((5,6,7)))
        for iter in iters:
            self.assertTrue(is_iterable(iter))
        # test non-iterables
        noniters = (5,)
        for iter in noniters:
            self.assertFalse(is_iterable(iter))

    def test_is_dicom(self):
        """Test the is_dicom function."""
    
        test_file = osp.join(osp.dirname(osp.dirname(__file__)), 'test_files', 'VMAT', 'DRGSdmlc-105-example.dcm')
        invalid_file = test_file.replace('DR', 'DR_')
        notdicom_file = osp.abspath(__file__)

        # valid file returns True
        self.assertTrue(is_dicom(test_file))

        # return false for real file but not dicom
        self.assertFalse(is_dicom(notdicom_file))

        # test invalid path
        self.assertRaises(IOError, is_dicom, invalid_file)

    def test_typed_property(self):

        class DumbClass:
            intprop = typed_property('dumbprop', int)
            floatstrprop = typed_property('floatstrprop', (float, str))

        dc = DumbClass()

        # test the intprop
        self.assertIsNone(dc.intprop)
        dc.intprop = 3
        self.assertEqual(dc.intprop, 3)
        self.assertRaises(TypeError, setattr, dc, 'intprop', 1.0)

        # test the intstrprop
        dc.floatstrprop = 3.3
        self.assertEqual(dc.floatstrprop, 3.3)
        dc.floatstrprop = 'mystring'
        self.assertEqual(dc.floatstrprop, dc._floatstrprop)
        self.assertRaises(TypeError, setattr, dc, 'floatstrprop', 3)

