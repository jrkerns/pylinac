import unittest

from pylinac.dlg import DLG
from pylinac.picketfence import MLC
from tests_basic.utils import requires_cloud_data


class TestDLG(unittest.TestCase):
    @requires_cloud_data(files={"file_path": ["DLG_1.5_0.2.dcm"]})
    def test_measured_dlg(self, file_path: str):
        dlg = DLG(file_path)
        dlg.analyze(gaps=(-0.9, -1.1, -1.3, -1.5, -1.7, -1.9), mlc=MLC.MILLENNIUM)
        self.assertAlmostEqual(dlg.measured_dlg, 1.503, delta=0.001)
