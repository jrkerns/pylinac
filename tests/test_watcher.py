import unittest
import os
import subprocess
import time
import os.path as osp
import shutil
import sys

test_files_dir = osp.join(osp.dirname(osp.abspath(__file__)), 'test_files')
watcher_dir = osp.join(test_files_dir, 'watcher')
watcher_script = osp.join(osp.dirname(osp.dirname(osp.abspath(__file__))), 'pylinac', 'watcher.py')


class WatcherTest(unittest.TestCase):
    """Test the watcher script by adding a file from every type of module."""
    files2analyze = [
        osp.join(test_files_dir, 'Starshot', 'Starshot#1.tif'),
        osp.join(test_files_dir, 'Picket Fence', 'AS500_PF.dcm'),
        osp.join(test_files_dir, 'VMAT', 'DRMLC.zip'),
        osp.join(test_files_dir, 'CBCT', 'GE_CT.zip'),
        osp.join(test_files_dir, 'MLC logs', 'tlogs', 'qqq2106_4DC Treatment_JS0_TX_20140712095629.bin')
    ]

    @classmethod
    def setUpClass(cls):
        """Setup the watcher."""
        cls.watching_process = subprocess.Popen([sys.executable, watcher_script, watcher_dir])
        time.sleep(1)

    @classmethod
    def tearDownClass(cls):
        """Kill watcher and remove all files from the watcher director."""
        cls.watching_process.kill()
        time.sleep(1)
        files = os.listdir(watcher_dir)
        files.remove('dummy.txt')
        for fyle in files:
            fyle = osp.join(watcher_dir, fyle)
            os.remove(fyle)

    def test_all(self):
        """Test all the files."""
        for fyle in self.files2analyze:
            self.process_file(fyle)

    def process_file(self, filepath):
        # copy the file over to the watcher directory
        shutil.copy(filepath, watcher_dir)
        new_path = osp.join(watcher_dir, osp.basename(filepath))

        # wait for processing to finish by checking when the .txt file is generated
        _, file_ext = osp.splitext(new_path)
        safety_timer = 0
        while not osp.isfile(new_path.replace(file_ext, '.txt')) and safety_timer < 30:
            time.sleep(0.5)
            safety_timer += 0.5
