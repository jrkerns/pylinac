import os.path as osp
from unittest import TestCase

import matplotlib.pyplot as plt

from pylinac import LeedsTOR, PipsProQC3
from tests.utils import save_file

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
        self.leeds.analyze()
        self.leeds.analyze(invert=True)

    def test_plotting(self):
        self.leeds.plot_analyzed_image()
        self.leeds.plot_analyzed_image(low_contrast=False, high_contrast=False)
        self.leeds.plot_analyzed_image(image=False, low_contrast=False, high_contrast=False)


class LeedsDemo(LeedsTORTestMixin, TestCase):

    @classmethod
    def setUpClass(cls):
        cls.leeds = LeedsTOR.from_demo_image()

    def test_demo(self):
        LeedsTOR.run_demo()  # shouldn't raise

    def test_url(self):
        LeedsTOR.from_url('https://s3.amazonaws.com/assuranceqa-staging/uploads/imgs/leeds.dcm')


class LeedsCCW(LeedsTORTestMixin, TestCase):
    filename = 'Leeds_ccw.dcm'


class PipsProTestMixin:
    filename = ''

    @classmethod
    def setUpClass(cls):
        cls.pipspro = PipsProQC3(osp.join(TEST_DIR, cls.filename))

    def test_analyze(self):
        self.pipspro.analyze()
        self.pipspro.analyze(invert=True)

    def test_plotting(self):
        self.pipspro.plot_analyzed_image()
        self.pipspro.plot_analyzed_image(low_contrast=False, high_contrast=False)
        self.pipspro.plot_analyzed_image(image=False, low_contrast=False, high_contrast=False)

    def test_saving(self):
        self.pipspro.plot_analyzed_image()
        save_file(self.pipspro.save_analyzed_image)

class PipsProDemo(PipsProTestMixin, TestCase):

    @classmethod
    def setUpClass(cls):
        cls.pipspro = PipsProQC3.from_demo_image()

    def test_demo(self):
        PipsProQC3.run_demo()  # shouldn't raise

    def test_url(self):
        PipsProQC3.from_url('https://s3.amazonaws.com/assuranceqa-staging/uploads/imgs/pipspro.dcm')


class PipsPro1(PipsProTestMixin, TestCase):
    filename = 'PIPSpro 2.5MV.dcm'
