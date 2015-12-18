import os.path as osp
from unittest import TestCase

import matplotlib.pyplot as plt

from pylinac import LeedsTOR, PipsPro

TEST_DIR = osp.join(osp.dirname(__file__), 'test_files', 'Planar imaging')


class LeedsTORTestMixin:
    filename = ''

    @classmethod
    def setUpClass(cls):
        cls.leeds = LeedsTOR(osp.join(TEST_DIR, cls.filename))

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def test_analyze(self):
        # test that the phantom runs
        self.leeds.analyze()
        self.leeds.plot_analyzed_image()


class LeedsCCW(LeedsTORTestMixin, TestCase):
    filename = 'Leeds_ccw.dcm'


class PipsProTestMixin:
    filename = ''

    @classmethod
    def setUpClass(cls):
        cls.pipspro = PipsPro(osp.join(TEST_DIR, cls.filename))

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def test_analyze(self):
        # test that the phantom runs
        self.pipspro.analyze()
        self.pipspro.plot_analyzed_image()


class PipsPro1(PipsProTestMixin, TestCase):
    filename = 'PIPSpro 2.5MV.dcm'