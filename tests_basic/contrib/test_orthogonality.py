import pytest
from pathlib import Path
from unittest import TestCase

from pylinac.contrib.orthogonality import JawOrthogonality
from tests_basic.utils import get_folder_from_cloud_repo


@pytest.mark.proprietary
class TestOrthogonality(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.dir = Path(get_folder_from_cloud_repo(["planar_imaging", "Orthogonality"]))

    def test_analyze(self):
        for file in self.dir.iterdir():
            if file.is_file():
                q = JawOrthogonality(file)
                q.analyze()
                q.plot_analyzed_image()
