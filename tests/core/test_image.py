import unittest
import os.path as osp

import dicom
import numpy as np

from pylinac.core.geometry import Point
from pylinac.core.image import Image, IMAGE, DICOM, ARRAY


img_path = osp.join(osp.dirname(osp.dirname(__file__)), 'test_files', 'Starshot', '6XCollStar.tif')
dcm_path = osp.join(osp.dirname(osp.dirname(__file__)), 'test_files', 'VMAT', 'DRGSmlc-105-example.dcm')

class Test_Image_Load(unittest.TestCase):

    def test_open(self):
        """Test the open class method."""
        # load a tif file
        img = Image(img_path)
        self.assertEqual(img.im_type, IMAGE)

        # load a dicom file
        img2 = Image(dcm_path)
        self.assertEqual(img2.im_type, DICOM)

        # try loading a bad file
        bad_file = osp.abspath(__file__)
        self.assertRaises(IOError, Image, bad_file)

        # load an array
        dcm = dicom.read_file(dcm_path)
        img = Image(dcm.pixel_array)
        self.assertEqual(img.im_type, ARRAY)


class Test_Image_Methods(unittest.TestCase):

    def setUp(self):
        self.img = Image(img_path)
        self.dcm = Image(dcm_path)
        small_array = np.arange(42).reshape(6,7)
        self.sm_arr = Image(small_array)

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

    def test_median_filter(self):
        filter_size = 3
        self.sm_arr.median_filter(filter_size)
        self.assertEqual(self.sm_arr.array[0, 0], 1)
        filter_size = 0.03
        self.sm_arr.median_filter(filter_size)

    # def test_rotate(self):
    #     self.sm_arr.rotate(90)
    #     self.assertEqual(self.sm_arr.array[0,0], 5)

    def test_resize(self):
        new_size = (200, 300)
        self.img.resize(new_size)
        self.assertEqual(self.img.shape, new_size)

    def test_invert(self):
        self.img.invert()

    def test_dpi_change(self):
        """Change the DPI and see if DPMM also changes."""
        new_dpi = 100
        self.img.dpi = new_dpi
        self.assertEqual(self.img.dpmm, new_dpi/25.4)

        new_dpmm = 3
        self.img.dpmm = new_dpmm
        self.assertEqual(self.img.dpi, new_dpmm*25.4)

    def test_dist2edge_min(self):
        dist = self.sm_arr.dist2edge_min(Point(1,3))
        self.assertEqual(dist, 1)

        dist = self.sm_arr.dist2edge_min((1,3))
        self.assertEqual(dist, 1)