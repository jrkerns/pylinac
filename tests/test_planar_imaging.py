import os.path as osp
from unittest import TestCase

from pylinac import LeedsTOR, PipsProQC3
from tests.utils import save_file

TEST_DIR = osp.join(osp.dirname(__file__), 'test_files', 'Planar imaging')


class PlanarPhantomMixin:
    klass = object
    filename = None

    @classmethod
    def setUpClass(cls):
        if cls.filename is None:
            cls.instance = cls.klass.from_demo_image()
        else:
            cls.instance = cls.klass(osp.join(TEST_DIR, cls.filename))

    def test_analyze(self):
        self.instance.analyze()
        self.instance.analyze(invert=True)

    def test_plotting(self):
        self.instance.plot_analyzed_image()
        self.instance.plot_analyzed_image(low_contrast=False, high_contrast=False)
        self.instance.plot_analyzed_image(image=False, low_contrast=False, high_contrast=False)

    def test_saving(self):
        self.instance.plot_analyzed_image()
        save_file(self.instance.save_analyzed_image)


class LeedsTORTestMixin(PlanarPhantomMixin):
    klass = LeedsTOR


class LeedsDemo(LeedsTORTestMixin, TestCase):

    def test_demo(self):
        LeedsTOR.run_demo()  # shouldn't raise

    def test_url(self):
        LeedsTOR.from_url('https://s3.amazonaws.com/assuranceqa-staging/uploads/imgs/leeds.dcm')


class LeedsCCW(LeedsTORTestMixin, TestCase):
    filename = 'Leeds_ccw.dcm'


class PipsProTestMixin(PlanarPhantomMixin):
    klass = PipsProQC3


class PipsProDemo(PipsProTestMixin, TestCase):

    def test_demo(self):
        PipsProQC3.run_demo()  # shouldn't raise

    def test_url(self):
        PipsProQC3.from_url('https://s3.amazonaws.com/assuranceqa-staging/uploads/imgs/pipspro.dcm')


class PipsPro1(PipsProTestMixin, TestCase):
    filename = 'PIPSpro 2.5MV.dcm'
