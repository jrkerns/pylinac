"""Test suite for the pylinac.io module."""
import unittest
import os
import os.path as osp

from pylinac.core.io import TemporaryZipDirectory, get_url, URLError


class TestIO(unittest.TestCase):

    def test_temp_zip_dir(self):
        """Test the TemporaryZipDirectory."""
        zfile = osp.join(osp.dirname(__file__), '..', 'test_files', 'VMAT', 'DRMLC.zip')

        # test context manager
        with TemporaryZipDirectory(zfile) as tmpzip:
            files = os.listdir(tmpzip)
            # test that they are real files
            osp.isfile(files[0])
            # test that both images were unpacked
            self.assertEqual(len(files), 2, msg="There were not 2 files found")

    def test_get_url(self):
        """Test the URL retreiver."""
        # test webpage
        webpage_url = 'http://google.com'
        get_url(webpage_url)  # shouldn't raise
        # test file
        file_url = 'https://s3.amazonaws.com/assuranceqa-staging/uploads/imgs/winston_lutz.zip'
        local_file = get_url(file_url)
        osp.isfile(local_file)
        # bad URL
        with self.assertRaises(URLError):
            get_url('http://asdfasdfasdfasdfasdfasdfasdfasdf.org')
