"""Test the watcher script by starting the script, moving files into the directory, and ensuring it processes the files correctly."""
import os
import os.path as osp
import shutil
import unittest

from pylinac.watcher import process
from tests import TEST_FILES_DIR
WATCHER_DIR = osp.join(TEST_FILES_DIR, 'watcher')
WATCHER_SCRIPT = osp.abspath(osp.join(osp.dirname(__file__), '../', 'pylinac', 'scripts.py'))
PYLINAC_DIR = osp.abspath(osp.join(osp.dirname(__file__), '../', 'pylinac'))


class TestProcess(unittest.TestCase):
    """Class for moving a set of files from a source dir to destination dir for the purposes
    of processing the destination files, but 1) retaining the original files and 2) cleaing the
    destination directory at the end of the tests.
    """
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
    DEST_DIR = WATCHER_DIR

    @classmethod
    def copy_files(cls):
        for file in cls.files:
            shutil.copy(file, cls.DEST_DIR)

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
        process(self.DEST_DIR)

    def test_process_files_with_classifier(self):
        config = osp.join(osp.dirname(__file__), 'test_config.yaml')
        process(self.DEST_DIR, config_file=config)



# class WatcherTest(unittest.TestCase):
#     """Test the watcher script by adding a file from every type of module."""
#     files2analyze = [
#         osp.join(TEST_FILES_DIR, 'Starshot', 'Starshot#1.tif'),
#         osp.join(TEST_FILES_DIR, 'Picket Fence', 'AS500_PF.dcm'),
#         osp.join(TEST_FILES_DIR, 'VMAT', 'DRMLC.zip'),
#         osp.join(TEST_FILES_DIR, 'CBCT', 'Varian', 'CBCT_4.zip'),
#         osp.join(TEST_FILES_DIR, 'MLC logs', 'tlogs', 'Anonymous_4DC Treatment_JS0_TX_20140712095629.bin'),
#         osp.join(TEST_FILES_DIR, 'Winston-Lutz', 'lutz.zip'),
#         osp.join(TEST_FILES_DIR, 'Planar imaging', 'Leeds_ccw.dcm'),
#         osp.join(TEST_FILES_DIR, 'Planar imaging', 'PIPSpro 2.5MV.dcm'),
#     ]
#
#     @classmethod
#     def setUpClass(cls):
#         """Setup the watcher."""
#         cls.clean_dir()
#         cls.watching_process = subprocess.Popen([sys.executable, WATCHER_SCRIPT, WATCHER_DIR])
#         time.sleep(8)  # process needs some time before it starts reading file changes
#
#     @classmethod
#     def tearDownClass(cls):
#         """Kill watcher and remove all files from the watcher director."""
#         cls.watching_process.kill()
#         time.sleep(1)
#         cls.clean_dir()
#
#     @staticmethod
#     def clean_dir():
#         """Clean the watcher directory."""
#         files = os.listdir(WATCHER_DIR)
#         files.remove('dummy.txt')
#         for fyle in files:
#             fyle = osp.join(WATCHER_DIR, fyle)
#             os.remove(fyle)
#         print("Files cleaned up")
#
#     def test_all(self):
#         """Test all the files."""
#         for fyle in self.files2analyze:
#             self.process_file(fyle)
#
#     def process_file(self, filepath):
#         # copy the file over to the watcher directory
#         shutil.copy(filepath, WATCHER_DIR)
#         print("{} moved into {}".format(osp.basename(filepath), WATCHER_DIR))
#         # construct new expected file name
#         new_path = osp.join(WATCHER_DIR, osp.basename(filepath))
#         filename_stub = osp.splitext(osp.basename(new_path))[0]
#         new_path = new_path.replace(filename_stub, filename_stub + '_analysis')
#         _, file_ext = osp.splitext(new_path)
#         new_path = new_path.replace(file_ext, '.png')
#
#         # wait for processing to finish by checking when the .txt file is generated
#         safety_timer = 0
#         while not osp.isfile(new_path):
#             time.sleep(0.3)
#             safety_timer += 0.3
#             if safety_timer > 30:
#                 print("Process timed out after {} seconds and was not completed.".format(int(safety_timer)))
#                 break
#         time.sleep(1)
