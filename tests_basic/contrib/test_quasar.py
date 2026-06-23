from pathlib import Path
from unittest import TestCase

from pylinac.contrib.quasar import QuasarLightRadScaling
from tests_basic.utils import requires_cloud_data


class TestQuasar(TestCase):
    @classmethod
    @requires_cloud_data(folders={"cloud_dir": ["planar_imaging", "Quasar"]})
    def setUpClass(cls, cloud_dir: str) -> None:
        cls.dir = Path(cloud_dir)

    def test_analyze(self):
        for file in self.dir.iterdir():
            if file.is_file():
                q = QuasarLightRadScaling(file)
                q.analyze()
                q.plot_analyzed_image()
