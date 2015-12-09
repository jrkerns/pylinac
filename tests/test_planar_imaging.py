import os.path as osp
from unittest import TestCase

from pylinac.planar_imaging import LeedsTOR

TEST_DIR = osp.join(osp.dirname(__file__), 'test_files', 'kV')


class LeedsTORTestMixin:
    filename = ''

    @classmethod
    def setUpClass(cls):
        cls.leeds = LeedsTOR(osp.join(TEST_DIR, cls.filename))

    def test_analyze(self):
        # test that the phantom runs
        self.leeds.analyze()
        self.leeds.plot_analyzed_image()


class LeedsCW(LeedsTORTestMixin, TestCase):
    filename = 'Leeds_cw.dcm'


class LeedsCCW(LeedsTORTestMixin, TestCase):
    filename = 'Leeds_ccw.dcm'


class Leeds1(LeedsTORTestMixin, TestCase):
    filename = 'Leeds_1.dcm'


class Leeds2(LeedsTORTestMixin, TestCase):
    filename = 'Leeds_2.dcm'
