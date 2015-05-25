import copy
import unittest
import os.path as osp
from tkinter.filedialog import askopenfilename

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
        self.assertRaises(TypeError, Image, bad_file)

        # not a valid parameter
        bad_input = 3.5
        self.assertRaises(TypeError, Image, bad_input)

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

        self.assertRaises(ValueError, self.img.median_filter, 1.1)

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

    def test_center(self):
        self.assertIsInstance(self.img.center, Point)
        img_known_center = Point(1420, 1702)
        dcm_known_center = Point(512, 384)
        self.assertEqual(self.img.center.x, img_known_center.x)
        self.assertEqual(self.dcm.center.y, dcm_known_center.y)

    def test_SID(self):
        self.assertEqual(self.dcm.SID, 105)
        self.assertRaises(ValueError, setattr, self.dcm, 'SID', '105')
        # test that the SID and DPMM change when set
        old_dpmm = copy.copy(self.dcm.dpmm)
        self.dcm.SID = 150
        self.assertEqual(self.dcm.SID, 150)
        self.assertEqual(self.dcm.dpmm, 150/100 * old_dpmm)

    def test_combine_multiples(self):
        bad_img_path = [dcm_path, img_path]
        self.assertRaises(AttributeError, Image.combine_multiples, bad_img_path)

        good_img_path = [img_path, img_path]
        combined_img = Image.combine_multiples(good_img_path)
        self.assertIsInstance(combined_img, Image)
