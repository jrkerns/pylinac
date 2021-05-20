"""Test suite for the pylinac.io module."""
import unittest
import os
import os.path as osp

from pylinac.core.io import TemporaryZipDirectory, get_url, URLError, is_dicom, SNCProfiler
from pylinac.core.profile import SingleProfile
from tests_basic.utils import get_file_from_cloud_test_repo


class TestIO(unittest.TestCase):

    def test_temp_zip_dir(self):
        """Test the TemporaryZipDirectory."""
        zfile = get_file_from_cloud_test_repo(['VMAT', 'DRMLC.zip'])

        # test context manager use; shouldn't raise
        with TemporaryZipDirectory(zfile) as tmpzip:
            files = [osp.join(tmpzip, file) for file in os.listdir(tmpzip)]
            # test that they are real files
            self.assertTrue(osp.isfile(files[0]))
            # test that both images were unpacked
            self.assertEqual(len(files), 2, msg="There were not 2 files found")

    def test_get_url(self):
        """Test the URL retreiver."""
        # test webpage
        webpage_url = 'http://google.com'
        get_url(webpage_url)  # shouldn't raise
        # test file
        file_url = 'https://s3.amazonaws.com/pylinac/winston_lutz.zip'
        local_file = get_url(file_url)
        osp.isfile(local_file)
        # bad URL
        with self.assertRaises(URLError):
            get_url('http://asdfasdfasdfasdfasdfasdfasdfasdf.org')

    def test_is_dicom(self):
        """Test the is_dicom function."""

        test_file = get_file_from_cloud_test_repo(['VMAT', 'DRGSdmlc-105-example.dcm'])
        invalid_file = test_file.replace('DR', 'DR_')
        notdicom_file = osp.abspath(__file__)

        # valid file returns True
        self.assertTrue(is_dicom(test_file))

        # return false for real file but not dicom
        self.assertFalse(is_dicom(notdicom_file))

        # test invalid path
        self.assertRaises(IOError, is_dicom, invalid_file)


class TestSNCProfiler(unittest.TestCase):

    def test_loading(self):
        path = get_file_from_cloud_test_repo(['9E-GA0.prs'])
        prof = SNCProfiler(path)
        self.assertEqual(len(prof.detectors), 254)

    def test_to_profiles(self):
        path = get_file_from_cloud_test_repo(['9E-GA0.prs'])
        prof = SNCProfiler(path)
        profs = prof.to_profiles()
        self.assertEqual(len(profs), 4)
        self.assertIsInstance(profs[0], SingleProfile)
