import io
import os
from pathlib import Path
from unittest import TestCase

from matplotlib import pyplot as plt
from scipy import ndimage

from pylinac.acr import CT464, ACRCTResult
from pylinac.core.geometry import Point
from pylinac.core.io import TemporaryZipDirectory
from tests_basic.utils import InitTesterMixin, FromZipTesterMixin, get_file_from_cloud_test_repo, save_file, \
    CloudFileMixin

TEST_DIR_CT = ['ACR', 'CT']


class TestACRCT(TestCase, FromZipTesterMixin, InitTesterMixin):
    klass = CT464
    zip = [*TEST_DIR_CT, 'Philips.zip']
    init_file = [*TEST_DIR_CT, 'Philips_CT']
    is_folder = True

    def test_load_from_list_of_paths(self):
        # shouldn't raise
        with TemporaryZipDirectory(get_file_from_cloud_test_repo(self.zip)) as zfolder:
            paths = [Path(zfolder) / f for f in os.listdir(zfolder)]
            CT464(paths)

    def test_load_from_list_of_streams(self):
        # shouldn't raise
        with TemporaryZipDirectory(get_file_from_cloud_test_repo(self.zip)) as zfolder:
            paths = [Path(zfolder, f) for f in os.listdir(zfolder)]
            paths = [io.BytesIO(open(p, 'rb').read()) for p in paths]
            CT464(paths)


class TestGeneral(TestCase):

    def setUp(self):
        path = get_file_from_cloud_test_repo([*TEST_DIR_CT, 'Philips.zip'])
        self.ct = CT464.from_zip(path)

    def test_phan_center(self):
        """Test locations of the phantom center."""
        known_phan_center = Point(258, 253)
        self.ct.analyze()
        self.assertAlmostEqual(self.ct.ct_calibration_module.phan_center.x, known_phan_center.x, delta=0.7)
        self.assertAlmostEqual(self.ct.ct_calibration_module.phan_center.y, known_phan_center.y, delta=0.7)

    def test_results_data(self):
        self.ct.analyze()
        data = self.ct.results_data()
        self.assertIsInstance(data, ACRCTResult)
        self.assertEqual(data.num_images, self.ct.num_images)

        # check the additional modules got added
        self.assertIsInstance(data.ct_module.rois, dict)


class TestPlottingSaving(TestCase):

    @classmethod
    def setUpClass(cls):
        path = get_file_from_cloud_test_repo([*TEST_DIR_CT, 'Philips.zip'])
        cls.ct = CT464.from_zip(path)
        cls.ct.analyze()

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def test_plot_images(self):
        """Test that saving an image does something."""
        save_file(self.ct.plot_images)

    def test_save_images(self):
        """Test that saving an image does something."""
        save_file(self.ct.save_images, to_single_file=False)

    def test_subimages_errors(self):
        """We don't use subimages here. easier to pass as a list of figs"""
        with self.assertRaises(NotImplementedError):
            self.ct.plot_analyzed_subimage('sr')
        with self.assertRaises(NotImplementedError):
            self.ct.save_analyzed_subimage('sr')

    def test_set_figure_size(self):
        self.ct.plot_analyzed_image(figsize=(8, 13))
        fig = plt.gcf()
        self.assertEqual(fig.bbox_inches.height, 13)
        self.assertEqual(fig.bbox_inches.width, 8)


class ACRCTMixin(CloudFileMixin):
    dir_path = ['ACR', 'CT']
    origin_slice: int
    phantom_roll: float = 0
    mtf_50: float
    slice_thickness: float
    hu_values: dict
    unif_values: dict

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        cls.ct = CT464.from_zip(filename)
        cls.ct.analyze()

    def test_roll(self):
        self.assertAlmostEqual(self.ct.catphan_roll, self.phantom_roll, delta=0.3)

    def test_mtf(self):
        self.assertAlmostEqual(self.ct.spatial_resolution_module.mtf.relative_resolution(50), self.mtf_50, delta=0.1)

    def test_HU_values(self):
        """Test HU values."""
        for key, roi in self.ct.ct_calibration_module.rois.items():
            exp_val = self.hu_values[key]
            meas_val = roi.pixel_value
            self.assertAlmostEqual(exp_val, meas_val, delta=5)

    def test_uniformity_values(self):
        """Test Uniformity HU values."""
        for key, exp_val in self.unif_values.items():
            meas_val = self.ct.uniformity_module.rois[key].pixel_value
            self.assertAlmostEqual(exp_val, meas_val, delta=5)


class ACRPhilips(ACRCTMixin, TestCase):
    file_name = 'Philips.zip'
    mtf_50 = 0.54
    phantom_roll = -0.3
    hu_values = {'Poly': -87, 'Acrylic': 126, 'Bone': 904, 'Air': -987, 'Water': 4}
    unif_values = {'Center': 1, 'Left': 1, 'Right': 1, 'Top': 1, 'Bottom': 1}


class ACRPhilipsOffset(ACRCTMixin, TestCase):
    """Shift the phantom over by several pixels to ensure no row/col algorithm issues

    Unfortunately, I can't move them that far because the FOV is very tight
    """
    file_name = 'Philips.zip'
    mtf_50 = 0.54
    phantom_roll = -0.3
    hu_values = {'Poly': -87, 'Acrylic': 126, 'Bone': 904, 'Air': -987, 'Water': 4}
    unif_values = {'Center': 1, 'Left': 1, 'Right': 1, 'Top': 1, 'Bottom': 1}

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        cls.ct = CT464.from_zip(filename)
        for img in cls.ct.dicom_stack:
            img.roll(direction='x', amount=5)
        cls.ct.localize()
        cls.ct.analyze()


class ACRPhilipsRotated(ACRCTMixin, TestCase):
    """Rotate the phantom over by several pixels to ensure no row/col algorithm issues

    Unfortunately, I can't move them that far because the FOV is very tight
    """
    file_name = 'Philips.zip'
    mtf_50 = 0.54
    phantom_roll = -3.3
    hu_values = {'Poly': -87, 'Acrylic': 126, 'Bone': 904, 'Air': -987, 'Water': 4}
    unif_values = {'Center': 1, 'Left': 1, 'Right': 1, 'Top': 1, 'Bottom': 1}

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        cls.ct = CT464.from_zip(filename)
        for img in cls.ct.dicom_stack:
            img.array = ndimage.rotate(img.array, angle=3, mode='nearest')
        cls.ct.localize()
        cls.ct.analyze()
