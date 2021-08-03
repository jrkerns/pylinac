import io
import unittest
from typing import Callable
from unittest import TestCase

import matplotlib.pyplot as plt
import numpy as np

from pylinac import LeedsTOR, StandardImagingQC3, LasVegas, DoselabMC2kV, DoselabMC2MV
from pylinac.planar_imaging import PlanarResult, SNCkV, SNCMV, StandardImagingQCkV, PTWEPIDQC
from tests_basic.utils import save_file, CloudFileMixin, get_file_from_cloud_test_repo

TEST_DIR = 'planar_imaging'


class GeneralTests(TestCase):

    def test_from_file_object(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, 'Leeds_ccw.dcm'])
        with open(path, 'rb') as f:
            phan = LeedsTOR(f)
            phan.analyze()
        self.assertIsInstance(phan, LeedsTOR)

    def test_from_stream(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, 'Leeds_ccw.dcm'])
        with open(path, 'rb') as f:
            s = io.BytesIO(f.read())
            phan = LeedsTOR(s)
            phan.analyze()
        self.assertIsInstance(phan, LeedsTOR)

    def test_overrides(self):
        phan = DoselabMC2kV.from_demo_image()
        phan.analyze()

    def test_results_data(self):
        phan = LeedsTOR.from_demo_image()
        phan.analyze()
        data = phan.results_data()
        self.assertIsInstance(data, PlanarResult)
        self.assertEqual(data.median_contrast, np.median([roi.contrast for roi in phan.low_contrast_rois]))

        data_dict = phan.results_data(as_dict=True)
        self.assertIsInstance(data_dict, dict)
        self.assertEqual(len(data_dict), 8)
        self.assertIn('pylinac_version', data_dict)

    def test_results_data_no_mtf(self):
        phan = LasVegas.from_demo_image()
        phan.analyze()

        data_dict = phan.results_data(as_dict=True)
        self.assertEqual(len(data_dict), 8)


class PlanarPhantomMixin(CloudFileMixin):
    klass: Callable
    dir_path = ['planar_imaging']
    mtf_50 = None
    invert = False
    ssd = 1000
    file_name = None

    @classmethod
    def setUpClass(cls):
        if not cls.file_name:
            cls.instance = cls.klass.from_demo_image()
        else:
            cls.instance = cls.klass(cls.get_filename())
        cls.instance.analyze(ssd=cls.ssd, invert=cls.invert)

    @classmethod
    def tearDownClass(cls):
        plt.close('all')
        del cls.instance

    def test_plotting(self):
        self.instance.plot_analyzed_image()
        self.instance.plot_analyzed_image(low_contrast=False, high_contrast=False)
        self.instance.plot_analyzed_image(image=False, low_contrast=False, high_contrast=False)

    def test_saving(self):
        self.instance.plot_analyzed_image()
        save_file(self.instance.save_analyzed_image)

    def test_pdf(self):
        save_file(self.instance.publish_pdf)

    def test_mtf(self):
        if self.mtf_50 is not None:
            self.assertAlmostEqual(self.mtf_50, self.instance.mtf.relative_resolution(50), delta=0.3)

    def test_results(self):
        self.assertIsInstance(self.instance.results(), str)


class LeedsDemo(PlanarPhantomMixin, TestCase):
    klass = LeedsTOR
    mtf_50 = 1.5

    def test_demo(self):
        LeedsTOR.run_demo()  # shouldn't raise


class LeedsCCW(PlanarPhantomMixin, TestCase):
    klass = LeedsTOR
    mtf_50 = 1.5
    file_name = 'Leeds_ccw.dcm'


class Leeds45Deg(PlanarPhantomMixin, TestCase):
    klass = LeedsTOR
    invert = True  # inverted in v3.0 due to changed default inversion behavior
    mtf_50 = 1.9
    ssd = 1500
    file_name = 'Leeds-45deg.dcm'


class LeedsDirtyEdges(PlanarPhantomMixin, TestCase):
    klass = LeedsTOR
    mtf_50 = 1.3
    ssd = 1000
    file_name = 'Leeds-dirty-edges.dcm'


@unittest.skip("Phantom appears distorted. MTF locations are different than other phantoms")
class LeedsClosedBlades(PlanarPhantomMixin, TestCase):
    klass = LeedsTOR
    mtf_50 = 1.3
    ssd = 1500
    file_name = 'Leeds-closed-blades.dcm'


class SIQC3Demo(PlanarPhantomMixin, TestCase):
    klass = StandardImagingQC3
    mtf_50 = 0.53

    def test_demo(self):
        StandardImagingQC3.run_demo()  # shouldn't raise


class SIQC3_1(PlanarPhantomMixin, TestCase):
    klass = StandardImagingQC3
    file_name = 'QC3-2.5MV.dcm'
    mtf_50 = 1.19


class SIQC3_2(PlanarPhantomMixin, TestCase):
    klass = StandardImagingQC3
    file_name = 'QC3-2.5MV-2.dcm'
    mtf_50 = 1.16

    def test_wrong_ssd_fails(self):
        self.instance = self.klass(self.get_filename())
        with self.assertRaises(ValueError):
            self.instance.analyze(ssd=1500)  # really at 1000


class LasVegasTestMixin(PlanarPhantomMixin):
    klass = LasVegas
    phantom_angle = 0

    def test_plotting(self):
        self.instance.plot_analyzed_image()
        self.instance.plot_analyzed_image(low_contrast=False)
        self.instance.plot_analyzed_image(image=False, low_contrast=False)

    def test_angle(self):
        self.assertAlmostEqual(self.instance.phantom_angle, self.phantom_angle, delta=1)


class LasVegasDemo(PlanarPhantomMixin, TestCase):
    klass = LasVegas
    phantom_angle = 0

    def test_demo(self):
        LasVegas.run_demo()  # shouldn't raise


class DoselabMVDemo(PlanarPhantomMixin, TestCase):
    klass = DoselabMC2MV
    mtf_50 = 0.54

    def test_demo(self):
        DoselabMC2MV.run_demo()


class DoselabkVDemo(PlanarPhantomMixin, TestCase):
    klass = DoselabMC2kV
    mtf_50 = 2.16

    def test_demo(self):
        DoselabMC2kV.run_demo()


class SNCkVDemo(PlanarPhantomMixin, TestCase):
    klass = SNCkV
    mtf_50 = 1.76

    def test_demo(self):
        SNCkV.run_demo()


class SNCMVDemo(PlanarPhantomMixin, TestCase):
    klass = SNCMV
    mtf_50 = 0.43

    def test_demo(self):
        SNCkV.run_demo()


class SIQCkVDemo(PlanarPhantomMixin, TestCase):
    klass = StandardImagingQCkV
    mtf_50 = 1.81

    def test_demo(self):
        StandardImagingQCkV.run_demo()


class PTWEPIDDemo(PlanarPhantomMixin, TestCase):
    klass = PTWEPIDQC
    mtf_50 = 0.79

    def test_demo(self):
        PTWEPIDQC.run_demo()
