import copy
from unittest import TestCase
import os.path as osp

import numpy as np

from pylinac.core.geometry import Point
from pylinac.core.image import Image, DicomImage, ArrayImage, FileImage


tif_path = osp.join(osp.dirname(osp.dirname(__file__)), 'test_files', 'Starshot', 'Starshot#1.tif')
png_path = osp.join(osp.dirname(osp.dirname(__file__)), 'test_files', 'Starshot', 'Starshot#1.png')
dcm_path = osp.join(osp.dirname(osp.dirname(__file__)), 'test_files', 'VMAT', 'DRGSdmlc-105-example.dcm')
dcm_url = 'https://github.com/jrkerns/pylinac/blob/master/pylinac/demo_files/picket_fence/EPID-PF-LR.dcm?raw=true'


class TestImage(TestCase):
    """Test the Image class."""

    def test_url(self):
        img = Image.load_url(dcm_url)
        self.assertIsInstance(img, DicomImage)

    def test_dicom(self):
        img = Image.load(dcm_path)
        self.assertIsInstance(img, DicomImage)

    def test_file(self):
        img = Image.load(tif_path)
        self.assertIsInstance(img, FileImage)

    def test_array(self):
        arr = np.arange(36).reshape(6, 6)
        img = Image.load(arr)
        self.assertIsInstance(img, ArrayImage)

    def test_multiples(self):
        paths = [dcm_path, dcm_path, dcm_path]
        img = Image.load_multiples(paths)
        self.assertIsInstance(img, DicomImage)

        # test non-superimposable images
        paths = [dcm_path, tif_path]
        with self.assertRaises(ValueError):
            Image.load_multiples(paths)

    def test_nonsense(self):
        with self.assertRaises(TypeError):
            Image.load('blahblah')


class TestBaseImage(TestCase):
    """Test the methods of BaseImage. Since it's a semi-abstract class, its subclasses (DicomImage,
    ArrayImage, and FileImage) are tested."""

    def setUp(self):
        self.img = Image.load(tif_path)
        self.dcm = Image.load(dcm_path)
        array = np.arange(42).reshape(6, 7)
        self.arr = Image.load(array)

    def test_remove_edges(self):
        """Remove the edges from a pixel array."""
        crop = 15
        orig_shape = self.img.shape
        orig_dpi = self.img.dpi
        self.img.remove_edges(crop)
        new_shape = self.img.shape
        new_dpi = self.img.dpi
        self.assertEqual(new_shape[0]+crop*2, orig_shape[0])
        # ensure original metadata is still the same
        self.assertEqual(new_dpi, orig_dpi)

    def test_filter(self):
        # test integer filter size
        filter_size = 3
        self.arr.filter(filter_size)
        self.assertEqual(self.arr.array[0, 0], 1)
        # test float filter size
        filter_size = 0.03
        self.arr.filter(filter_size)
        # test using invalid float value
        self.assertRaises(TypeError, self.img.filter, 1.1)
        # test using a gaussian filter
        self.arr.filter(kind='gaussian')

    def test_ground(self):
        old_min_val = copy.copy(self.dcm.array.min())
        ground_val = self.dcm.ground()
        self.assertEqual(old_min_val, ground_val)
        # test that array was also changed
        self.assertAlmostEqual(self.dcm.array.min(), 0)

    def test_resize(self):
        new_size = (200, 300)
        self.img.resize(new_size)
        self.assertEqual(self.img.shape, new_size)

    def test_invert(self):
        self.img.invert()
        self.dcm.invert()
        self.arr.invert()

    def test_dist2edge_min(self):
        dist = self.arr.dist2edge_min(Point(1,3))
        self.assertEqual(dist, 1)

        dist = self.arr.dist2edge_min((1,3))
        self.assertEqual(dist, 1)

    def test_center(self):
        self.assertIsInstance(self.img.center, Point)
        img_known_center = Point(512, 1702)
        dcm_known_center = Point(512, 384)
        self.assertEqual(self.img.center.x, img_known_center.x)
        self.assertEqual(self.dcm.center.y, dcm_known_center.y)

    def test_plot(self):
        self.img.plot()  # shouldn't raise

    def test_roll(self):
        orig_val = self.arr[0, 0]
        self.arr.roll()
        shifted_val = self.arr[0, 1]
        self.assertEqual(orig_val, shifted_val)

    def test_rot90(self):
        orig_val = self.dcm[0, 0]
        self.dcm.rot90()
        shifted_val = self.dcm[-1, 0]
        self.assertEqual(orig_val, shifted_val)

    def test_check_inversion(self):
        # original dicom needs inversion
        orig_val = self.dcm[0, 0]
        self.dcm.check_inversion()
        inverted_val = self.dcm[0, 0]
        self.assertNotEqual(orig_val, inverted_val)

    def test_threshold(self):
        # apply high-pass threshold
        orig_val = self.arr[0, 4]
        self.arr.threshold(threshold=10)
        zeroed_val = self.arr[0, 4]
        self.assertNotEqual(orig_val, zeroed_val)
        self.assertEqual(zeroed_val, 0)

        # apply low-pass threshold
        orig_val = self.arr[-1, -1]
        self.arr.threshold(threshold=20, kind='low')
        zeroed_val = self.arr[-1, -1]
        self.assertNotEqual(orig_val, zeroed_val)
        self.assertEqual(zeroed_val, 0)

    def test_gamma(self):
        array = np.arange(49).reshape((7,7))
        ref_img = Image.load(array, dpi=1)
        comp_img = Image.load(array, dpi=1)
        comp_img.roll(amount=1)
        g_map = ref_img.gamma(comp_img)

        num_passing_pixels = np.nansum(g_map < 1)
        num_calced_pixels = np.nansum(g_map >= 0)
        pass_pct = num_passing_pixels / num_calced_pixels * 100
        average_gamma = np.nanmean(g_map)
        expected_pass_pct = 86
        expected_avg_gamma = 0.78
        self.assertAlmostEqual(pass_pct, expected_pass_pct, delta=1)
        self.assertAlmostEqual(average_gamma, expected_avg_gamma, delta=0.02)


class TestDicomImage(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dcm = Image.load(dcm_path)

    def test_sid(self):
        self.assertEqual(self.dcm.sid, 1050)

    def test_dpi(self):
        self.assertAlmostEqual(self.dcm.dpi, 68, delta=0.1)

    def test_dpmm(self):
        self.assertAlmostEqual(self.dcm.dpmm, 2.68, delta=0.01)

    def test_cax(self):
        self.assertLessEqual(self.dcm.cax.distance_to(self.dcm.center), 1.5)


class TestFileImage(TestCase):

    def test_sid(self):
        # default sid is 1000
        fi = FileImage(tif_path)
        self.assertEqual(fi.sid, 1000)

        # SID can be set though
        fi2 = FileImage(tif_path, sid=1500)
        self.assertEqual(fi2.sid, 1500)

        # SID also affects the dpi
        orig_dpi = fi.dpi
        scaled_dpi = fi2.dpi
        self.assertEqual(orig_dpi, scaled_dpi*2/3)

    def test_dpi_dpmm(self):
        # DPI is usually in TIF files
        fi = FileImage(tif_path)
        # shouldn't raise
        fi.dpi
        fi.dpmm

        # not in certain other files
        fi_jpg = FileImage(png_path)
        with self.assertRaises(AttributeError):
            fi_jpg.dpi
        # but DPI can be set though
        fi_jpg2 = FileImage(png_path, dpi=100)
        # shouldn't raise
        fi_jpg2.dpi
        fi_jpg2.dpmm


class TestArrayImage(TestCase):

    def test_dpmm(self):
        arr = np.arange(42).reshape(6, 7)
        ai = ArrayImage(arr)
        # errors if not passed in
        with self.assertRaises(AttributeError):
            ai.dpi

        ai2 = ArrayImage(arr, dpi=20)
        self.assertEqual(ai2.dpi, 20)
        self.assertEqual(ai2.dpmm, 20/25.4)
