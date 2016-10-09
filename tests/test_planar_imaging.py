import os.path as osp
from unittest import TestCase

import matplotlib.pyplot as plt

from pylinac import LeedsTOR, StandardImagingQC3
from tests.utils import save_file, LocationMixin

TEST_DIR = osp.join(osp.dirname(__file__), 'test_files', 'Planar imaging')


class PlanarPhantomMixin(LocationMixin):
    klass = object
    dir_location = TEST_DIR

    @classmethod
    def setUpClass(cls):
        if not cls.file_path:
            cls.instance = cls.klass.from_demo_image()
        else:
            cls.instance = cls.klass(cls.get_filename())

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

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


class LeedsCCW(LeedsTORTestMixin, TestCase):
    file_path = ['Leeds_ccw.dcm']


class SIQC3TestMixin(PlanarPhantomMixin):
    klass = StandardImagingQC3


class SIQC3Demo(SIQC3TestMixin, TestCase):

    def test_demo(self):
        StandardImagingQC3.run_demo()  # shouldn't raise


class SIQC3_1(SIQC3TestMixin, TestCase):
    file_path = ['PIPSpro 2.5MV.dcm']
