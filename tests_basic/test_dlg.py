import unittest
from pathlib import Path

from pylinac.dlg import DLG
from pylinac.picketfence import MLC


class TestDLG(unittest.TestCase):
    file_path = Path(__file__).parent / Path('test_files/DLG 1.5+0.2.dcm')

    def test_measured_dlg(self):
        dlg = DLG(self.file_path)
        dlg.analyze(gaps=(-0.9, -1.1, -1.3, -1.5, -1.7, -1.9), mlc=MLC.MILLENNIUM)
        self.assertAlmostEqual(dlg.measured_dlg, 1.503, delta=0.001)
