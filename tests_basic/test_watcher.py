"""Test the watcher script by starting the script, moving files into the directory, and ensuring it processes the files correctly."""
import os
import os.path as osp
import shutil
import unittest

from pylinac.watcher import process
from pylinac.core.io import retrieve_filenames
from tests_basic import TEST_FILES_DIR
WATCHER_DIR = osp.join(TEST_FILES_DIR, 'watcher')
# WATCHER_SCRIPT = osp.abspath(osp.join(osp.dirname(__file__), '../', 'pylinac', 'scripts.py'))
# PYLINAC_DIR = osp.abspath(osp.join(osp.dirname(__file__), '../', 'pylinac'))


class TestProcess(unittest.TestCase):
    """Class for moving a set of files from a source dir to destination dir for the purposes
    of processing the destination files, but 1) retaining the original files and 2) cleaing the
    destination directory at the end of the tests_basic.
    """

    @classmethod
    def copy_files(cls):
        files = [
            osp.join(TEST_FILES_DIR, 'Starshot', 'Starshot#1.tif'),
            osp.join(TEST_FILES_DIR, 'Picket Fence', 'AS500_PF.dcm'),
            osp.join(TEST_FILES_DIR, 'VMAT', 'DRMLC.zip'),
            osp.join(TEST_FILES_DIR, 'CBCT', 'CBCT_4.zip'),
            osp.join(TEST_FILES_DIR, 'MLC logs', 'tlogs', 'Anonymous_4DC Treatment_JS0_TX_20140712095629.bin'),
            osp.join(TEST_FILES_DIR, 'Winston-Lutz', 'lutz.zip'),
            osp.join(TEST_FILES_DIR, 'Planar imaging', 'Leeds_ccw.dcm'),
            osp.join(TEST_FILES_DIR, 'Planar imaging', 'QC3 2.5MV.dcm'),
        ]
        files += retrieve_filenames(osp.join(TEST_FILES_DIR, 'CBCT', 'Pelvis'))
        for file in files:
            shutil.copy(file, WATCHER_DIR)
        print("Files copied over")

    @staticmethod
    def clean_dir():
        """Clean the watcher directory."""
        files = os.listdir(WATCHER_DIR)
        files.remove('dummy.txt')
        for fyle in files:
            fyle = osp.join(WATCHER_DIR, fyle)
            os.remove(fyle)
        print("Files cleaned up")

    def setUp(self):
        """Clear out ."""
        self.clean_dir()
        self.copy_files()

    def tearDown(self):
        """Remove all files from the destination directory."""
        self.clean_dir()

    def test_process_files(self):
        process(WATCHER_DIR, copy_new_files=True)
