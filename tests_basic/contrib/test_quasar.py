from pathlib import Path
from unittest import TestCase

from pylinac.contrib.quasar import QuasarLightRadScaling
from tests_basic.utils import get_folder_from_cloud_test_repo


class TestQuasar(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.dir = Path(get_folder_from_cloud_test_repo(["planar_imaging", "Quasar"]))

    def test_analyze(self):
        for file in self.dir.iterdir():
            if file.is_file():
                q = QuasarLightRadScaling(file)
                q.analyze()
                q.plot_analyzed_image()
