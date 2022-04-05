import io
import os.path as osp
import unittest
from typing import Callable
from unittest import TestCase

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pylinac import LeedsTOR, StandardImagingQC3, LasVegas, DoselabMC2kV, DoselabMC2MV
from pylinac.core import image
from pylinac.planar_imaging import PlanarResult, SNCkV, SNCMV, StandardImagingQCkV, PTWEPIDQC, StandardImagingFC2
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

    def test_multiple_plots(self):
        phan = LeedsTOR.from_demo_image()
        phan.analyze()
        figs, names = phan.plot_analyzed_image(split_plots=True)
        self.assertEqual(len(figs), 3)
        files = phan.save_analyzed_image(filename='a.png', split_plots=True)
        names = ('a_image.png', 'a_low_contrast.png', 'a_high_contrast.png')
        for name in names:
            self.assertIn(name, files)

        # regular single plot produces one image/file
        figs, names = phan.plot_analyzed_image()
        self.assertEqual(len(figs), 0)
        name = 'b.png'
        phan.save_analyzed_image('b.png')
        self.assertTrue(osp.isfile(name))

        # stream buffer shouldn't fail
        with io.BytesIO() as tmp:
            phan.save_analyzed_image(tmp)

        # to streams should return streams
        streams = phan.save_analyzed_image(split_plots=True, to_streams=True)
        self.assertEqual(len(streams.keys()), 3)
        with self.assertRaises(ValueError):
            phan.save_analyzed_image()  # no filename and no streams is an error

    def test_passing_image_kwargs(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, 'Leeds_ccw.dcm'])

        # do normal analysis
        phan = LeedsTOR(path)
        phan.analyze()
        x = phan.results_data().phantom_center_x_y[0]

        # pass kwarg; use same dpi as image; results should be the same.
        img = image.load(path)
        phan = LeedsTOR(path, image_kwargs={'dpi': img.dpi})
        phan.analyze()
        x_manual_dpi = phan.results_data().phantom_center_x_y[0]

        self.assertEqual(x, x_manual_dpi)


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


class FC2Mixin(PlanarPhantomMixin):
    klass = StandardImagingFC2
    dir_path = ['planar_imaging', 'SI FC2']
    field_size_x_mm = 150
    field_size_y_mm = 150
    field_epid_offset_x_mm = 0
    field_epid_offset_y_mm = 0
    field_bb_offset_x_mm = 0
    field_bb_offset_y_mm = 0
    fwxm = 50

    @classmethod
    def setUpClass(cls):
        if not cls.file_name:
            cls.instance = cls.klass.from_demo_image()
        else:
            cls.instance = cls.klass(cls.get_filename())
        cls.instance.analyze(invert=cls.invert, fwxm=cls.fwxm)

    def test_plotting(self):
        self.instance.plot_analyzed_image()

    def test_field_size(self):
        results_data = self.instance.results_data()
        assert results_data.field_size_x_mm == pytest.approx(self.field_size_x_mm, abs=0.3)
        assert results_data.field_size_y_mm == pytest.approx(self.field_size_y_mm, abs=0.3)
        assert results_data.field_epid_offset_x_mm == pytest.approx(self.field_epid_offset_x_mm, abs=0.2)
        assert results_data.field_epid_offset_y_mm == pytest.approx(self.field_epid_offset_y_mm, abs=0.2)
        assert results_data.field_bb_offset_x_mm == pytest.approx(self.field_bb_offset_x_mm, abs=0.2)
        assert results_data.field_bb_offset_y_mm == pytest.approx(self.field_bb_offset_y_mm, abs=0.2)


class FC2Demo(FC2Mixin, TestCase):
    klass = StandardImagingFC2
    field_size_x_mm = 148.5
    field_size_y_mm = 149.1
    field_epid_offset_x_mm = -0.7
    field_epid_offset_y_mm = 0.3
    field_bb_offset_x_mm = -1.2
    field_bb_offset_y_mm = 1.2

    def test_demo(self):
        StandardImagingFC2.run_demo()


class FC210x10_10FFF(FC2Mixin, TestCase):
    file_name = 'FC-2-10x10-10fff.dcm'
    klass = StandardImagingFC2
    field_size_x_mm = 98.7
    field_size_y_mm = 99.3
    field_epid_offset_x_mm = 0.2
    field_bb_offset_y_mm = 0.8


class FC210x10_10X(FC2Mixin, TestCase):
    file_name = 'FC-2-10x10-10x.dcm'
    klass = StandardImagingFC2
    field_size_x_mm = 99.3
    field_size_y_mm = 99.6
    field_epid_offset_y_mm = 0.2
    field_epid_offset_x_mm = -0.5
    field_bb_offset_y_mm = 1.1
    field_bb_offset_x_mm = -0.8


class FC210x10_15X(FC2Mixin, TestCase):
    file_name = 'FC-2-10x10-15x.dcm'
    klass = StandardImagingFC2
    field_size_x_mm = 99.3
    field_size_y_mm = 99.6
    field_epid_offset_y_mm = 0.1
    field_epid_offset_x_mm = -0.5
    field_bb_offset_y_mm = 1.1
    field_bb_offset_x_mm = -0.8


class FC215x15_10X(FC2Mixin, TestCase):
    file_name = 'FC-2-15x15-10X.dcm'
    klass = StandardImagingFC2
    field_size_y_mm = 149.2
    field_size_x_mm = 149.2
    field_epid_offset_y_mm = 0.1
    field_epid_offset_x_mm = -0.5
    field_bb_offset_y_mm = 1.1
    field_bb_offset_x_mm = -0.8


class FC215x15_10FFF(FC2Mixin, TestCase):
    file_name = 'FC-2-15x15-10XFFF.dcm'
    fwxm = 30
    klass = StandardImagingFC2
    field_size_y_mm = 149.5
    field_size_x_mm = 149.6
    field_epid_offset_y_mm = -0.1
    field_epid_offset_x_mm = 0.2
    field_bb_offset_y_mm = 0.8
    field_bb_offset_x_mm = -0.1
