"""Test the watcher script by starting the script, moving files into the directory, and ensuring it processes the files correctly."""
import os
import os.path as osp
import subprocess
import shutil
import sys
import time
import unittest

TEST_DIR = osp.join(osp.dirname(osp.abspath(__file__)), 'test_files')
WATCHER_DIR = osp.join(TEST_DIR, 'watcher')

PYLINAC_DIR = osp.join(osp.dirname(osp.dirname(osp.abspath(__file__))), 'pylinac')
WATCHER_SCRIPT = osp.join(PYLINAC_DIR, 'watcher.py')


class WatcherTest(unittest.TestCase):
    """Test the watcher script by adding a file from every type of module."""
    files2analyze = [
        osp.join(TEST_DIR, 'Starshot', 'Starshot#1.tif'),
        osp.join(TEST_DIR, 'Picket Fence', 'AS500_PF.dcm'),
        osp.join(TEST_DIR, 'VMAT', 'DRMLC.zip'),
        osp.join(TEST_DIR, 'CBCT', 'Varian', 'CBCT_4.zip'),
        osp.join(TEST_DIR, 'MLC logs', 'tlogs', 'Anonymous_4DC Treatment_JS0_TX_20140712095629.bin'),
        osp.join(PYLINAC_DIR, 'demo_files', 'winston_lutz', 'winston_lutz.zip')
    ]

    @classmethod
    def setUpClass(cls):
        """Setup the watcher."""
        cls.watching_process = subprocess.Popen([sys.executable, WATCHER_SCRIPT, WATCHER_DIR])
        time.sleep(3)  # process needs some time before it starts reading file changes

    @classmethod
    def tearDownClass(cls):
        """Kill watcher and remove all files from the watcher director."""
        cls.watching_process.kill()
        time.sleep(1)
        files = os.listdir(WATCHER_DIR)
        files.remove('dummy.txt')
        for fyle in files:
            fyle = osp.join(WATCHER_DIR, fyle)
            os.remove(fyle)

    def test_all(self):
        """Test all the files."""
        for fyle in self.files2analyze:
            self.process_file(fyle)

    def process_file(self, filepath):
        # copy the file over to the watcher directory
        shutil.copy(filepath, WATCHER_DIR)
        print("{} moved into watcher folder".format(filepath))
        new_path = osp.join(WATCHER_DIR, osp.basename(filepath))

        # wait for processing to finish by checking when the .txt file is generated
        _, file_ext = osp.splitext(new_path)
        safety_timer = 0
        while not osp.isfile(new_path.replace(file_ext, '.txt')):
            time.sleep(0.5)
            safety_timer += 0.5
            if safety_timer > 10:
                print("Process timed out after 30 seconds and was not completed.")
                break

