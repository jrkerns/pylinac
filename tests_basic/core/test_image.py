import copy
import io
import json
import shutil
import tempfile
import unittest
from builtins import ValueError
from pathlib import Path
from unittest import TestCase

import numpy as np
import PIL.Image
import pydicom
from numpy.testing import assert_array_almost_equal
from pydicom.uid import UID

from pylinac.core import image
from pylinac.core.geometry import Point
from pylinac.core.image import (
    XIM,
    ArrayImage,
    DicomImage,
    DicomImageStack,
    FileImage,
    LazyDicomImageStack,
    LinacDicomImage,
    _rescale_dicom_values,
    _unscale_dicom_values,
    equate_images,
    gamma_2d,
    load,
    tiff_to_dicom,
)
from pylinac.core.io import TemporaryZipDirectory
from tests_basic.utils import (
    get_file_from_cloud_test_repo,
    get_folder_from_cloud_test_repo,
    save_file,
)

bad_tif_path = get_file_from_cloud_test_repo(["Winston-Lutz", "AQA_A_03082023.tif"])
tif_path = get_file_from_cloud_test_repo(["Starshot", "Starshot-1.tif"])
png_path = get_file_from_cloud_test_repo(["Starshot", "Starshot-1.png"])
dcm_path = get_file_from_cloud_test_repo(["VMAT", "DRGSdmlc-105-example.dcm"])
as500_path = get_file_from_cloud_test_repo(["picket_fence", "AS500#5.dcm"])
xim_path = get_file_from_cloud_test_repo(["ximdcmtest.xim"])
xim_dcm_path = get_file_from_cloud_test_repo(["ximdcmtest.dcm"])
dcm_url = "https://storage.googleapis.com/pylinac_demo_files/EPID-PF-LR.dcm"


class TestFromMultiples(TestCase):
    def test_from_multiples(self):
        """Test that loading multiple images works"""
        paths = [dcm_path, dcm_path, dcm_path]
        img = image.load_multiples(paths)
        self.assertIsInstance(img, DicomImage)

    def test_different_sizes_fails(self):
        """Different-sized images can't be superimposed"""
        paths = [dcm_path, tif_path]
        with self.assertRaises(ValueError):
            image.load_multiples(paths)

    def test_round_trip_via_temp_file(self):
        """Test that loading multiple images works"""
        paths = [dcm_path, dcm_path]
        img = image.load_multiples(paths)
        with tempfile.NamedTemporaryFile(delete=False) as tf:
            img.save(tf.name)
            # shouldn't raise
            img_cropped = image.load(tf.name)
        self.assertEqual(img.shape, img_cropped.shape)

    def test_round_trip_via_stream(self):
        """Test that loading multiple images works"""
        paths = [dcm_path, dcm_path]
        img = image.load_multiples(paths)
        with io.BytesIO() as stream:
            img.save(stream)
            stream.seek(0)
            # shouldn't raise
            img_cropped = image.load(stream)
        self.assertEqual(img.shape, img_cropped.shape)

    def test_max_stays_same(self):
        """Test the max value of the image stays the same when using max method"""
        paths = [dcm_path, dcm_path]
        img = image.load_multiples(paths, method="max", stretch_each=False)
        self.assertEqual(np.max(image.load(dcm_path).array), np.max(img.array))

    def test_max_goes_to_1_with_stretch(self):
        """Test the max value of the image goes to 1 when using max method and stretch"""
        paths = [dcm_path, dcm_path]
        img = image.load_multiples(paths, method="max", stretch_each=True)
        self.assertEqual(np.max(img.array), 1)

    def test_rescale_to_int16_max_stretch_true(self):
        """Test that loading multiple images works"""
        paths = [dcm_path, dcm_path]
        img = image.load_multiples(
            paths, stretch_each=True, method="max"
        )  # will cause a value >1
        with io.BytesIO() as stream:
            img.save(stream)  # converts back to uint16, unscaled
            stream.seek(0)
            img_cropped = image.load(stream, raw_pixels=True)  # load unscaled
        self.assertEqual(np.max(img_cropped.array), 65535)  # max of uint16

    def test_rescale_to_int16_max_stretch_false(self):
        """Test that loading multiple images works"""
        paths = [dcm_path, dcm_path]
        img = image.load_multiples(
            paths, stretch_each=False, method="max"
        )  # will still be a float array.
        with io.BytesIO() as stream:
            img.save(stream)  # converts back to uint16, unscaled
            stream.seek(0)
            img_cropped = image.load(stream, raw_pixels=True)  # load unscaled
        self.assertEqual(np.max(img_cropped.array), 65535)  # max of uint16


class TestDICOMScaling(TestCase):
    def test_raw_pixels_doesnt_change_array(self):
        """Test that loading a dicom with raw_pixels=True doesn't change the array"""
        ds = pydicom.dcmread(dcm_path)
        array = _rescale_dicom_values(ds.pixel_array, ds, raw_pixels=True)
        assert np.array_equal(ds.pixel_array, array)

    def test_mr_storage_doesnt_rescale(self):
        """Test that loading a dicom with MR storage doesn't rescale"""
        dcm_path = get_file_from_cloud_test_repo(["ACR", "MRI", "GE - 3T", "IM_0001"])
        ds = pydicom.dcmread(dcm_path)
        array = _rescale_dicom_values(ds.pixel_array, ds, raw_pixels=False)
        assert np.array_equal(ds.pixel_array, array)

    def test_ct_image(self):
        """Test an older CT Image where we are guaranteed to have rescale slope and intercept tags"""
        dcm_path = get_file_from_cloud_test_repo(
            ["CBCT", "CatPhan_504", "Case3_Philips_1mm", "1mm", "EE035381"]
        )
        ds = pydicom.dcmread(dcm_path)
        array = _rescale_dicom_values(ds.pixel_array, ds, raw_pixels=False)
        assert not np.array_equal(ds.pixel_array, array)
        assert np.isclose(array.max(), 194)

    def test_no_tags_inverts(self):
        ds = pydicom.Dataset()
        ds.SOPClassUID = UID("1.2.840.122332")  # junk UID; not a real image type
        array = np.zeros((3, 3))
        array[0, 0] = 100
        scaled_array = _rescale_dicom_values(array, ds, raw_pixels=False)
        assert not np.array_equal(array, scaled_array)
        assert np.isclose(
            scaled_array[0, 0], 0
        )  # test we inverted and the value is now opposite
        assert np.isclose(scaled_array.max(), 100)

    def test_all_tags_negative_sign(self):
        ds = pydicom.Dataset()
        ds.RescaleSlope = 1
        ds.RescaleIntercept = -1000
        ds.PixelIntensityRelationshipSign = -1
        ds.SOPClassUID = UID("1.2.840.122332")  # junk UID; not a real image type
        array = np.ones((3, 3))
        scaled_array = _rescale_dicom_values(array, ds, raw_pixels=False)
        # with all tags, this will be negative and offset
        assert not np.array_equal(array, scaled_array)
        assert np.isclose(scaled_array.max(), 999)

    def test_all_tags_positive_sign(self):
        ds = pydicom.Dataset()
        ds.RescaleSlope = 1
        ds.RescaleIntercept = -1000
        ds.PixelIntensityRelationshipSign = 1
        ds.SOPClassUID = UID("1.2.840.122332")  # junk UID; not a real image type
        array = np.ones((3, 3))
        scaled_array = _rescale_dicom_values(array, ds, raw_pixels=False)
        # with all tags, this will be negative and offset
        assert not np.array_equal(array, scaled_array)
        assert np.isclose(scaled_array.max(), -999)


class TestDICOMUnscaling(TestCase):
    def test_unscale_raw_pixels(self):
        """Test when we unscale the image that the values are the same"""
        ds = pydicom.dcmread(dcm_path)
        array = _unscale_dicom_values(ds.pixel_array, ds, raw_pixels=True)
        assert np.array_equal(ds.pixel_array, array)

    def test_unscale_mr_storage(self):
        dcm_path = get_file_from_cloud_test_repo(["ACR", "MRI", "GE - 3T", "IM_0001"])
        ds = pydicom.dcmread(dcm_path)
        array = _unscale_dicom_values(ds.pixel_array, ds, raw_pixels=False)
        assert np.array_equal(ds.pixel_array, array)

    def test_unscale_ct_image(self):
        """Test an older CT Image where we are guaranteed to have rescale slope and intercept tags"""
        dcm_path = get_file_from_cloud_test_repo(
            ["CBCT", "CatPhan_504", "Case3_Philips_1mm", "1mm", "EE035381"]
        )
        ds = pydicom.dcmread(dcm_path)
        original_array = ds.pixel_array
        scaled_array = _rescale_dicom_values(ds.pixel_array, ds, raw_pixels=False)
        unscaled_array = _unscale_dicom_values(scaled_array, ds, raw_pixels=False)
        assert np.array_equal(original_array, unscaled_array)

    def test_no_tags_inverts_back(self):
        ds = pydicom.Dataset()
        ds.SOPClassUID = UID("1.2.840.122332")  # junk UID; not a real image type
        original_array = np.zeros((3, 3), dtype=np.uint16)
        original_array[0, 0] = 100
        scaled_array = _rescale_dicom_values(original_array, ds, raw_pixels=False)
        assert not np.array_equal(original_array, scaled_array)
        unscaled_array = _unscale_dicom_values(scaled_array, ds, raw_pixels=False)
        assert np.array_equal(original_array, unscaled_array)

    def test_all_tags_negative_sign(self):
        ds = pydicom.Dataset()
        ds.RescaleSlope = 1
        ds.RescaleIntercept = -1000
        ds.PixelIntensityRelationshipSign = -1
        ds.SOPClassUID = UID("1.2.840.122332")  # junk UID; not a real image type
        original_array = np.ones((3, 3))
        scaled_array = _rescale_dicom_values(original_array, ds, raw_pixels=False)
        assert not np.array_equal(scaled_array, original_array)
        unscaled_array = _unscale_dicom_values(scaled_array, ds, raw_pixels=False)
        assert np.array_equal(unscaled_array, original_array)


class TestEquateImages(TestCase):
    def test_same_sized_images_work(self):
        """As found here: https://github.com/jrkerns/pylinac/issues/446"""

        gen = np.random.default_rng()
        image1 = load(gen.random((20, 20)), dpi=10)
        image2 = load(gen.random((10, 10)), dpi=5)

        img1, img2 = equate_images(image1, image2)
        self.assertEqual(img1.shape, img2.shape)


class TestLoaders(TestCase):
    """Test the image loading functions."""

    def test_load_url(self):
        img = image.load_url(dcm_url)
        self.assertIsInstance(img, DicomImage)

    def test_load_dicom(self):
        img = image.load(dcm_path)
        self.assertIsInstance(img, DicomImage)

    def test_load_dicom_from_stream(self):
        img_ref = image.load(dcm_path)
        with open(dcm_path, "rb") as f:
            p = io.BytesIO(f.read())
            img = image.load(p)
        self.assertIsInstance(img, DicomImage)
        self.assertEqual(img.dpi, img_ref.dpi)

    def test_load_dicom_from_file_object(self):
        img_ref = image.load(dcm_path)
        with open(dcm_path, "rb") as f:
            img = image.load(f)
        self.assertIsInstance(img, DicomImage)
        self.assertEqual(img.dpi, img_ref.dpi)

    def test_load_file(self):
        img = image.load(tif_path)
        self.assertIsInstance(img, FileImage)

    def test_load_file_from_stream(self):
        img_ref = image.load(tif_path)
        with open(tif_path, "rb") as f:
            p = io.BytesIO(f.read())
            img = image.load(p)
        self.assertIsInstance(img, FileImage)
        self.assertEqual(img.path, "")
        self.assertEqual(img.center, img_ref.center)

    def test_load_file_from_temp_file(self):
        img_ref = image.load(tif_path)
        tmp = tempfile.NamedTemporaryFile(delete=False)
        tmp.write(open(tif_path, "rb").read())
        img = image.load(tmp)
        self.assertIsInstance(img, FileImage)
        self.assertIsInstance(img.path, str)
        self.assertGreater(len(img.path), 10)  # has temp name
        self.assertEqual(img.center, img_ref.center)

    def test_load_file_from_file_object(self):
        img_ref = image.load(tif_path)
        with open(tif_path, "rb") as f:
            img = image.load(f)
        self.assertIsInstance(img, FileImage)
        self.assertEqual(img.center, img_ref.center)

    def test_load_array(self):
        arr = np.arange(36).reshape(6, 6)
        img = image.load(arr)
        self.assertIsInstance(img, ArrayImage)

    def test_load_multiples(self):
        paths = [dcm_path, dcm_path, dcm_path]
        img = image.load_multiples(paths)
        self.assertIsInstance(img, DicomImage)

        # test non-superimposable images
        paths = [dcm_path, tif_path]
        with self.assertRaises(ValueError):
            image.load_multiples(paths)

    def test_load_multiples_custom_loader(self):
        """Use a custom loader to load multiple images"""
        paths = [dcm_path, dcm_path, dcm_path]
        img = image.load_multiples(paths, loader=image.LinacDicomImage)
        self.assertIsInstance(img, image.LinacDicomImage)

    def test_nonsense(self):
        with self.assertRaises(FileNotFoundError):
            image.load("blahblah")

    def test_is_image(self):
        self.assertTrue(image.is_image(dcm_path))
        # not an image
        self.assertFalse(
            image.is_image(
                get_file_from_cloud_test_repo(["mlc_logs", "dlogs", "Adlog1.dlg"])
            )
        )


class TestBaseImage(TestCase):
    """Test the basic methods of BaseImage. Since it's a semi-abstract class, its subclasses (DicomImage,
    ArrayImage, and FileImage) are tested."""

    def setUp(self):
        self.img = image.load(tif_path)
        self.dcm = image.load(dcm_path)
        array = np.arange(42).reshape(6, 7)
        self.arr = image.load(array)

    def test_crop(self):
        """Remove the edges from a pixel array."""
        crop = 15
        orig_shape = self.img.shape
        orig_dpi = self.img.dpi
        self.img.crop(crop)
        new_shape = self.img.shape
        new_dpi = self.img.dpi
        self.assertEqual(new_shape[0] + crop * 2, orig_shape[0])
        # ensure original metadata is still the same
        self.assertEqual(new_dpi, orig_dpi)

    def test_crop_must_be_positive(self):
        """Crop must be manifestly positive"""
        crop = 0
        with self.assertRaises(ValueError):
            self.img.crop(crop)

    def test_filter(self):
        # test integer filter size
        filter_size = 3
        self.arr.filter(filter_size)
        self.assertEqual(self.arr.array[0, 0], 1)
        # test float filter size
        filter_size = 0.03
        self.arr.filter(filter_size)
        # test using invalid float value
        self.assertRaises(ValueError, self.img.filter, 1.1)
        # test using a gaussian filter
        self.arr.filter(kind="gaussian")

    def test_ground(self):
        old_min_val = copy.copy(self.dcm.array.min())
        ground_val = self.dcm.ground()
        self.assertEqual(old_min_val, ground_val)
        # test that array was also changed
        self.assertAlmostEqual(self.dcm.array.min(), 0)

    def test_invert(self):
        self.img.invert()
        self.dcm.invert()
        self.arr.invert()

    def test_dist2edge_min(self):
        dist = self.arr.dist2edge_min(Point(1, 3))
        self.assertEqual(dist, 1)

        dist = self.arr.dist2edge_min((1, 3))
        self.assertEqual(dist, 1)

    def test_center(self):
        self.assertIsInstance(self.img.center, Point)
        img_known_center = Point(511.5, 1702)
        dcm_known_center = Point(511.5, 383.5)
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

    def test_rotate(self):
        # rotate by 90 because it's easy to test
        # we test where the value was relatively high and
        # ensure it's now low
        orig_val = self.dcm[150, 500]
        self.dcm.rotate(angle=90)
        shifted_val = self.dcm[150, 500]
        self.assertNotAlmostEqual(orig_val, shifted_val, delta=0.2)

    def test_threshold(self):
        # apply high-pass threshold
        orig_val = self.arr[0, 4]
        self.arr.threshold(threshold=10)
        zeroed_val = self.arr[0, 4]
        self.assertNotEqual(orig_val, zeroed_val)
        self.assertEqual(zeroed_val, 0)

        # apply low-pass threshold
        orig_val = self.arr[-1, -1]
        self.arr.threshold(threshold=20, kind="low")
        zeroed_val = self.arr[-1, -1]
        self.assertNotEqual(orig_val, zeroed_val)
        self.assertEqual(zeroed_val, 0)

    @unittest.skip("Skip until overhaul of gamma method to true gamma")
    def test_gamma(self):
        array = np.arange(49).reshape((7, 7))
        ref_img = image.load(array, dpi=1)
        comp_img = image.load(array, dpi=1)
        comp_img.roll(amount=1)
        g_map = ref_img.gamma(comp_img)

        num_passing_pixels = np.nansum(g_map < 1)
        num_calced_pixels = np.nansum(g_map >= 0)
        pass_pct = num_passing_pixels / num_calced_pixels * 100
        average_gamma = np.nanmean(g_map)
        expected_pass_pct = 81
        expected_avg_gamma = 0.87
        self.assertAlmostEqual(pass_pct, expected_pass_pct, delta=1)
        self.assertAlmostEqual(average_gamma, expected_avg_gamma, delta=0.02)


class TestDicomImage(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dcm: DicomImage = image.load(dcm_path)

    def test_sid(self):
        self.assertEqual(self.dcm.sid, 1050)

    def test_dpi(self):
        self.assertAlmostEqual(self.dcm.dpi, 68, delta=0.1)

    def test_dpmm(self):
        self.assertAlmostEqual(self.dcm.dpmm, 2.68, delta=0.01)

    def test_cax(self):
        self.assertLessEqual(self.dcm.cax.distance_to(self.dcm.center), 1.5)

    def test_save(self):
        save_file(self.dcm.save)

    def test_save_round_trip_has_same_pixel_values(self):
        original_dcm = pydicom.dcmread(dcm_path)
        dcm_img = DicomImage(dcm_path)
        with io.BytesIO() as stream:
            dcm_img.save(stream)
            stream.seek(0)
            reloaded_dcm = pydicom.dcmread(stream)
        # full array comparison fails due to some rounding,
        # but testing min/max and argmax should be sufficient
        self.assertAlmostEqual(
            reloaded_dcm.pixel_array.max(), original_dcm.pixel_array.max(), places=-1
        )
        self.assertAlmostEqual(
            reloaded_dcm.pixel_array.min(), original_dcm.pixel_array.min(), places=-1
        )
        self.assertEqual(
            reloaded_dcm.pixel_array.argmax(), original_dcm.pixel_array.argmax()
        )

    def test_save_out_of_bounds_values_normalizes(self):
        # this occurs if we add multiple images together and the values go out of bounds
        # causing a bit overflow
        original_dcm = pydicom.dcmread(dcm_path)
        dcm_img = DicomImage(dcm_path)
        dcm_img.array *= 10e6  # values out of normal bounds of uint16
        with io.BytesIO() as stream, self.assertWarns(UserWarning):
            dcm_img.save(stream)
            stream.seek(0)
            reloaded_dcm = pydicom.dcmread(stream)
        # because values are out of bounds, the values will get stretched
        # and equality cannot be assumed
        # max will be uint16 (original datatype) max
        self.assertAlmostEqual(reloaded_dcm.pixel_array.max(), 65535, places=-1)
        # indices should be the same however.
        self.assertEqual(
            reloaded_dcm.pixel_array.argmax(), original_dcm.pixel_array.argmax()
        )

    def test_manipulation_still_saves_correctly(self):
        dcm = image.load(dcm_path)
        original_shape = dcm.shape
        dcm.crop(15)
        with tempfile.NamedTemporaryFile(delete=False) as tf:
            dcm.save(tf.name)
            # shouldn't raise
            dcm_cropped = image.load(tf.name)
        self.assertEqual(
            original_shape[0] - 30, dcm_cropped.shape[0]
        )  # 15 from each side = 2 * 15 = 30
        self.assertEqual(
            original_shape[1] - 30, dcm_cropped.shape[1]
        )  # 15 from each side = 2 * 15 = 30

    def test_raw_pixels_via_load(self):
        dcm_ds = pydicom.dcmread(dcm_path)
        dcm_raw = image.load(dcm_path, raw_pixels=True)
        assert np.array_equal(dcm_raw.array, dcm_ds.pixel_array)
        # test that the pixels are corrected otherwise
        dcm_corr = image.load(dcm_path)
        assert not np.array_equal(dcm_corr.array, dcm_ds.pixel_array)

    def test_z_location_mri(self):
        path = get_file_from_cloud_test_repo(["ACR", "MRI", "GE - 3T", "IM_0008"])
        img: DicomImage = image.load(path)
        self.assertAlmostEqual(img.z_position, -20.2, places=1)

    def test_z_location_ct(self):
        path = get_file_from_cloud_test_repo(
            ["CBCT", "CatPhan_504", "Case3_Philips_1mm", "1mm", "EE035381"]
        )
        img: DicomImage = image.load(path)
        self.assertAlmostEqual(img.z_position, -500.5, places=1)

    def test_z_location_none_for_epid(self):
        """EPID has no z location"""
        with self.assertRaises(AttributeError):
            self.dcm.z_position

    def test_slice_spacing_mri(self):
        path = get_file_from_cloud_test_repo(["ACR", "MRI", "GE - 3T", "IM_0008"])
        img: DicomImage = image.load(path)
        self.assertAlmostEqual(img.slice_spacing, 10, places=1)

    def test_slice_spacing_ct(self):
        path = get_file_from_cloud_test_repo(
            ["CBCT", "CatPhan_504", "Case3_Philips_1mm", "1mm", "EE035381"]
        )
        img: DicomImage = image.load(path)
        self.assertAlmostEqual(img.slice_spacing, 0.5, places=1)

    def test_negative_slice_spacing_ct(self):
        path = get_file_from_cloud_test_repo(
            ["CBCT", "CatPhan_604", "negative_spacing.zip"]
        )
        stack = DicomImageStack.from_zip(path)
        self.assertAlmostEqual(stack[1].slice_spacing, 2, places=1)


class TestXIMImage(TestCase):
    def test_normal_load(self):
        xim = XIM(xim_path)
        self.assertIsInstance(xim.array, np.ndarray)
        self.assertEqual(xim.array.shape, (1280, 1280))
        self.assertIsInstance(xim.properties, dict)

    def test_high_diff_files(self):
        """RAM-2414; 2-byte 'diff' sizes were causing XIM to choke. See ticket/PR for more"""
        # shouldn't choke; this file has 1-element switches betwen byte sizes
        XIM(get_file_from_cloud_test_repo(["xim_1_element_diff_switch.xim"]))

        XIM(get_file_from_cloud_test_repo(["xim_2byte_diff.xim"]))

    def test_dont_read_pixels(self):
        xim = XIM(xim_path, read_pixels=False)
        with self.assertRaises(AttributeError):
            xim.array
        self.assertIsInstance(xim.properties, dict)

    def test_equivalent_to_dcm(self):
        """The pixel info should be the same between dicom and xim"""
        dcm_img = DicomImage(xim_dcm_path)
        xim_img = XIM(xim_path)
        assert_array_almost_equal(dcm_img.array, xim_img.array)

    def test_save_png_stream(self):
        xim = XIM(xim_path)
        s = io.BytesIO()
        xim.save_as(s, format="png")
        s.seek(0)
        pimg = PIL.Image.open(s)
        png_array = np.asarray(pimg)
        assert_array_almost_equal(png_array, xim.array)

    def test_save_png(self):
        xim = XIM(xim_path)
        with tempfile.NamedTemporaryFile(delete=False) as tf:
            xim.save_as(tf, format="png")
            pimg = PIL.Image.open(tf.name)
        png_array = np.asarray(pimg)
        assert_array_almost_equal(png_array, xim.array)
        # make sure the properties were saved
        assert (
            xim.properties["AcquisitionSystemVersion"]
            == pimg.info["AcquisitionSystemVersion"]
        )
        mlc_a = json.loads(pimg.info["MLCLeafsA"])
        self.assertIsInstance(mlc_a, list)

    def test_save_tiff(self):
        xim = XIM(xim_path)
        with tempfile.NamedTemporaryFile(delete=False) as tf:
            xim.save_as(tf, format="tiff")
            pimg = PIL.Image.open(tf.name)
        png_array = np.asarray(pimg)
        assert_array_almost_equal(png_array, xim.array)

    def test_overflow(self):
        old_settings = np.geterr()
        # set overflow to cause errors
        # shouldn't raise
        np.seterr(all="raise")
        xim_path = get_file_from_cloud_test_repo(["IsoCal-kV-08.xim"])
        XIM(xim_path)
        # reset to old settings
        np.seterr(**old_settings)


class TestLinacDicomImage(TestCase):
    def test_normal_image(self):
        img = LinacDicomImage(as500_path, use_filenames=False)
        self.assertEqual(img.gantry_angle, 0)

    def test_passing_axis_info_through_filename(self):
        new_name = as500_path.replace("AS500#5", "AS500Gantry78Coll13Couch44")
        shutil.copy(as500_path, new_name)
        img = LinacDicomImage(new_name, use_filenames=True)
        self.assertEqual(img.gantry_angle, 78)
        self.assertEqual(img.collimator_angle, 13)
        self.assertEqual(img.couch_angle, 44)

    def test_passing_axis_info_directly(self):
        img = LinacDicomImage(
            as500_path, use_filenames=False, gantry=24, coll=60, couch=8
        )
        self.assertEqual(img.gantry_angle, 24)
        self.assertEqual(img.collimator_angle, 60)
        self.assertEqual(img.couch_angle, 8)

    def test_using_axes_no_precision_doesnt_round(self):
        img = LinacDicomImage(
            as500_path,
            use_filenames=False,
            axes_precision=None,
            gantry=0.12345,
            coll=89.483,
            couch=359.6145,  # None is default
        )
        self.assertEqual(img.gantry_angle, 0.12345)
        self.assertEqual(img.collimator_angle, 89.483)
        self.assertEqual(img.couch_angle, 359.6145)

    def test_using_axes_2precision(self):
        img = LinacDicomImage(
            as500_path,
            use_filenames=False,
            axes_precision=2,
            gantry=0.12345,
            coll=89.483,
            couch=359.6145,
        )
        self.assertEqual(img.gantry_angle, 0.12)
        self.assertEqual(img.collimator_angle, 89.48)
        self.assertEqual(img.couch_angle, 359.61)

    def test_using_axes_0precision(self):
        img = LinacDicomImage(
            as500_path,
            use_filenames=False,
            axes_precision=0,
            gantry=0.12345,
            coll=89.48,
            couch=359.6145,
        )
        self.assertEqual(img.gantry_angle, 0)
        self.assertEqual(img.collimator_angle, 89)
        self.assertEqual(img.couch_angle, 0)


class TestFileImage(TestCase):
    def test_sid(self):
        # default sid is None
        fi = FileImage(tif_path)
        self.assertIsNone(fi.sid)

        # SID can be set though
        fi2 = FileImage(tif_path, sid=1500)
        self.assertEqual(fi2.sid, 1500)

        # SID also affects the dpi
        orig_dpi = fi.dpi
        scaled_dpi = fi2.dpi
        self.assertEqual(orig_dpi, scaled_dpi * 2 / 3)

    def test_dpi_dpmm(self):
        # DPI is usually in TIF files
        fi = FileImage(tif_path)
        # shouldn't raise
        fi.dpi
        fi.dpmm

        # not in certain other files
        fi_jpg = FileImage(png_path)
        self.assertIsNone(fi_jpg.dpi)

        # but DPI can be set though
        fi_jpg2 = FileImage(png_path, dpi=100)
        # shouldn't raise
        fi_jpg2.dpi
        fi_jpg2.dpmm

    def test_dpi_abnormal(self):
        # has DPI of 1. Nonsensical
        fimg = FileImage(bad_tif_path)
        with self.assertRaises(ValueError):
            fimg.dpi


class TestTiff(TestCase):
    """A special case of the FileImage"""

    def test_all_tiffs_have_tags_and_are_2d(self):
        """Test all tiffs will load. Just raw ingestion"""
        all_starshot_files = get_folder_from_cloud_test_repo(["Starshot"])
        for img in Path(all_starshot_files).iterdir():
            if img.suffix in (".tif", ".tiff"):
                fimg = FileImage(img)
                self.assertTrue(fimg, "tags")
                self.assertEqual(len(fimg.array.shape), 2)


class TestArrayImage(TestCase):
    def test_dpmm(self):
        arr = np.arange(42).reshape(6, 7)
        ai = ArrayImage(arr)

        self.assertIsNone(ai.dpi)

        ai2 = ArrayImage(arr, dpi=20)
        self.assertEqual(ai2.dpi, 20)
        self.assertEqual(ai2.dpmm, 20 / 25.4)


class TestDicomStack(TestCase):
    stack_location = get_file_from_cloud_test_repo(["CBCT", "CBCT_4.zip"])

    def test_loading(self):
        # test normal construction
        with TemporaryZipDirectory(self.stack_location) as tmpzip:
            dstack = DicomImageStack(tmpzip)
            self.assertEqual(len(dstack), 64)
        # test zip
        dstack = DicomImageStack.from_zip(self.stack_location)

    def test_lazy_and_normal_have_same_length(self):
        dstack = DicomImageStack.from_zip(self.stack_location)
        dstack_lazy = LazyDicomImageStack.from_zip(self.stack_location)
        self.assertEqual(len(dstack), len(dstack_lazy))

    def test_images_are_the_same(self):
        dstack = DicomImageStack.from_zip(self.stack_location)
        dstack_lazy = LazyDicomImageStack.from_zip(self.stack_location)
        for img, img_lazy in zip(dstack, dstack_lazy):
            assert_array_almost_equal(img.array, img_lazy.array)

    @unittest.skip("Wait until better error checking is implemented")
    def test_mixed_studies(self):
        mixed_study_zip = get_file_from_cloud_test_repo(["CBCT", "mixed_studies.zip"])
        with self.assertRaises(ValueError):
            DicomImageStack.from_zip(mixed_study_zip)

    def test_writing_back_to_lazy_stack(self):
        dstack_lazy = LazyDicomImageStack.from_zip(self.stack_location)
        original_offset = np.copy(dstack_lazy[0].array) + 50
        for idx, img in enumerate(dstack_lazy):
            img.array += 50
            dstack_lazy[idx] = img
        manipulated_offset = np.copy(dstack_lazy[0].array)
        # assert that writing back to the stack works
        assert_array_almost_equal(manipulated_offset, original_offset)

    def test_slice_spacing_ct(self):
        dstack = DicomImageStack.from_zip(self.stack_location)
        self.assertAlmostEqual(dstack.slice_spacing, 2.5, delta=0.001)

    def test_slice_spacing_mri(self):
        path = get_file_from_cloud_test_repo(["ACR", "MRI", "T1-Single.zip"])
        dstack = DicomImageStack.from_zip(path, min_number=5)
        self.assertAlmostEqual(dstack.slice_spacing, 10, delta=0.001)


class TestGamma2D(TestCase):
    def test_perfect_match_is_0(self):
        ref = eval = np.ones((5, 5))
        gamma = gamma_2d(reference=ref, evaluation=eval)
        self.assertEqual(gamma.max(), 0)
        self.assertEqual(gamma.min(), 0)
        self.assertEqual(gamma.size, 25)

        # test a high measurement value
        ref = eval = np.ones((5, 5)) * 50
        gamma = gamma_2d(reference=ref, evaluation=eval)
        self.assertEqual(gamma.max(), 0)
        self.assertEqual(gamma.min(), 0)
        self.assertEqual(gamma.size, 25)

    def test_gamma_perfectly_at_1(self):
        # offset a profile exactly by the dose to agreement
        ref = np.ones((5, 5))
        eval = np.ones((5, 5)) * 1.01
        gamma = gamma_2d(reference=ref, evaluation=eval, dose_to_agreement=1)
        self.assertAlmostEqual(gamma.max(), 1, delta=0.001)
        self.assertAlmostEqual(gamma.min(), 1, delta=0.001)

        # test same but eval is LOWER than ref
        ref = np.ones((5, 5))
        eval = np.ones((5, 5)) * 0.99
        gamma = gamma_2d(reference=ref, evaluation=eval, dose_to_agreement=1)
        self.assertAlmostEqual(gamma.max(), 1, delta=0.001)
        self.assertAlmostEqual(gamma.min(), 1, delta=0.001)

    def test_gamma_some_on_some_off(self):
        ref = np.ones((5, 5))
        eval = np.ones((5, 5))
        eval[(0, 0, 1, 1), (0, 1, 1, 0)] = 1.03  # set top left corner to 3% off
        gamma = gamma_2d(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=1,
            distance_to_agreement=1,
            gamma_cap_value=5,
        )
        self.assertAlmostEqual(gamma[0, 0], 3, delta=0.01)  # fully off by 3
        self.assertAlmostEqual(
            gamma[0, 1], 1, delta=0.01
        )  # dose at next pixel matches (dose=0, dist=1)
        self.assertAlmostEqual(gamma[-1, -1], 0, delta=0.01)  # gamma at end is perfect

        # check inverted pattern is mirrored (checks off-by-one errors)
        ref = np.ones((5, 5))
        eval = np.ones((5, 5))
        eval[
            (-1, -1, -2, -2), (-1, -2, -2, -1)
        ] = 1.03  # set bottom right corner to 3% off
        gamma = gamma_2d(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=1,
            distance_to_agreement=1,
            gamma_cap_value=5,
        )
        self.assertAlmostEqual(gamma[0, 0], 0, delta=0.01)
        self.assertAlmostEqual(gamma[-1, -2], 1, delta=0.01)
        self.assertAlmostEqual(gamma[-1, -1], 3, delta=0.01)

    def test_localized_dose(self):
        ref = np.ones((5, 5))
        ref[0, 0] = 100
        eval = np.ones((5, 5))
        eval[0, 0] = 103
        eval[0, 1] = 1.03
        # with global, element 2 is easily under gamma 1 since DTA there is 3
        gamma = gamma_2d(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=3,
            distance_to_agreement=1,
            gamma_cap_value=5,
            global_dose=False,
            dose_threshold=0,
        )
        self.assertAlmostEqual(gamma[0, 0], 1, delta=0.01)  # fully off by 3
        self.assertAlmostEqual(
            gamma[0, 1], 1, delta=0.01
        )  # dose here is also off by 3% relative dose
        self.assertAlmostEqual(gamma[-1, -1], 0, delta=0.01)  # gamma at end is perfect

    def test_threshold(self):
        ref = np.zeros((5, 5))
        ref[0, 0] = 1
        eval = ref
        # only one point should be computed as rest are under default threshold
        gamma = gamma_2d(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=3,
            distance_to_agreement=1,
            gamma_cap_value=5,
            global_dose=False,
            dose_threshold=5,
        )
        self.assertAlmostEqual(gamma[0, 0], 0, delta=0.01)
        self.assertTrue(np.isnan(gamma[0, 1]))
        self.assertTrue(np.isnan(gamma[-1, -1]))

    def test_fill_value(self):
        ref = np.zeros((5, 5))
        ref[0, 0] = 1
        eval = ref
        # only one point should be computed as rest are under default threshold
        gamma = gamma_2d(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=3,
            distance_to_agreement=1,
            gamma_cap_value=5,
            global_dose=False,
            dose_threshold=5,
            fill_value=0.666,
        )
        self.assertAlmostEqual(gamma[0, 0], 0, delta=0.01)
        self.assertAlmostEqual(
            gamma[0, 1], 0.666, delta=0.01
        )  # dose here is also off by 3% relative dose
        self.assertAlmostEqual(gamma[-1, -1], 0.666, delta=0.01)

    def test_gamma_half(self):
        # offset a profile by half the dose to agreement to ensure it's 0.5
        ref = np.ones((5, 5))
        eval = np.ones((5, 5)) / 1.005
        gamma = gamma_2d(reference=ref, evaluation=eval, dose_to_agreement=1)
        self.assertAlmostEqual(gamma.max(), 0.5, delta=0.01)
        self.assertAlmostEqual(gamma.min(), 0.5, delta=0.01)

    def test_gamma_cap(self):
        # cap to the value
        ref = np.ones((5, 5))
        eval = np.ones((5, 5)) * 10
        gamma = gamma_2d(
            reference=ref, evaluation=eval, dose_to_agreement=1, gamma_cap_value=2
        )
        self.assertEqual(gamma.max(), 2)
        self.assertEqual(gamma.min(), 2)

    def test_non_2d_array(self):
        ref = np.ones(5)
        eval = np.ones((5, 5))
        with self.assertRaises(ValueError):
            gamma_2d(reference=ref, evaluation=eval)

        ref = np.ones((5, 5))
        eval = np.ones(5)
        with self.assertRaises(ValueError):
            gamma_2d(reference=ref, evaluation=eval)


class TestTiffToDicom(TestCase):
    def test_conversion_can_be_loaded_as_dicom(self):
        ds = tiff_to_dicom(
            tif_path,
            sid=1000,
            dpi=200,
            gantry=10,
            coll=22,
            couch=33,
        )
        # shouldn't raise
        LinacDicomImage.from_dataset(ds)

    def test_conversion_captures_axes(self):
        ds = tiff_to_dicom(
            tif_path,
            sid=1000,
            dpi=200,
            gantry=10,
            coll=22,
            couch=33,
        )
        dicom_img = LinacDicomImage.from_dataset(ds)
        self.assertEqual(dicom_img.gantry_angle, 10)
        self.assertEqual(dicom_img.collimator_angle, 22)
        self.assertEqual(dicom_img.couch_angle, 33)

    def test_conversion_of_dpmm(self):
        ds = tiff_to_dicom(tif_path, sid=1000, gantry=10, coll=22, couch=33)
        dicom_img = LinacDicomImage.from_dataset(ds)
        self.assertEqual(dicom_img.dpi, 150)
        self.assertEqual(dicom_img.dpmm, 150 / 25.4)

    def test_conversion_goes_to_uint16(self):
        tiff_img = FileImage(tif_path)
        ds = tiff_to_dicom(
            tif_path,
            sid=1000,
            dpi=200,
            gantry=10,
            coll=22,
            couch=33,
        )
        self.assertEqual(tiff_img.array.dtype, np.uint8)
        dicom_img = LinacDicomImage.from_dataset(ds)
        self.assertEqual(dicom_img.array.dtype, np.uint16)

    def test_mass_conversion(self):
        """Mass conversion; shouldn't fail. All images have dpi tag"""
        all_starshot_files = get_folder_from_cloud_test_repo(["Starshot"])
        for img in Path(all_starshot_files).iterdir():
            if img.suffix in (".tif", ".tiff"):
                tiff_to_dicom(img, sid=1000, gantry=10, coll=11, couch=12)
