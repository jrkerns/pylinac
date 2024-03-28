import io
import json
import os.path as osp
from typing import Callable
from unittest import TestCase, skip

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import rotate

from pylinac import (
    DoselabMC2kV,
    DoselabMC2MV,
    DoselabRLf,
    ElektaLasVegas,
    IBAPrimusA,
    LasVegas,
    LeedsTOR,
    StandardImagingQC3,
)
from pylinac.core import image
from pylinac.planar_imaging import (
    PTWEPIDQC,
    SNCFSQA,
    SNCMV,
    SNCMV12510,
    IMTLRad,
    IsoAlign,
    LeedsTORBlue,
    LightRadResult,
    PlanarResult,
    SNCkV,
    StandardImagingFC2,
    StandardImagingQCkV,
    percent_integral_uniformity,
)
from tests_basic.utils import CloudFileMixin, get_file_from_cloud_test_repo, save_file

TEST_DIR = "planar_imaging"


class TestPercentIntegralUniformity(TestCase):
    def test_normal(self):
        self.assertAlmostEqual(
            percent_integral_uniformity(max=1000, min=900), 94.73, delta=0.1
        )

    def test_perfect(self):
        self.assertAlmostEqual(
            percent_integral_uniformity(max=1000, min=1000), 100, delta=0.1
        )

    def test_min_0(self):
        self.assertAlmostEqual(
            percent_integral_uniformity(max=1000, min=0), 0, delta=0.1
        )


class GeneralTests(TestCase):
    def test_from_file_object(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, "Leeds", "Leeds_ccw.dcm"])
        with open(path, "rb") as f:
            phan = LeedsTOR(f)
            phan.analyze()
        self.assertIsInstance(phan, LeedsTOR)

    def test_from_stream(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, "Leeds", "Leeds_ccw.dcm"])
        with open(path, "rb") as f:
            s = io.BytesIO(f.read())
            phan = LeedsTOR(s)
            phan.analyze()
        self.assertIsInstance(phan, LeedsTOR)

    def test_overrides(self):
        phan = DoselabMC2kV.from_demo_image()
        phan.analyze()

    def test_results(self):
        phan = LeedsTOR.from_demo_image()
        phan.analyze()
        data = phan.results()
        self.assertIsInstance(data, str)

        data_list = phan.results(as_list=True)
        self.assertIsInstance(data_list, list)
        self.assertEqual(len(data_list), 9)

    def test_results_data(self):
        phan = LeedsTOR.from_demo_image()
        phan.analyze()
        data = phan.results_data()
        self.assertIsInstance(data, PlanarResult)
        self.assertEqual(
            data.median_contrast,
            np.median([roi.contrast for roi in phan.low_contrast_rois]),
        )

        data_dict = phan.results_data(as_dict=True)
        self.assertIsInstance(data_dict, dict)
        self.assertEqual(len(data_dict), 11)
        self.assertIn("pylinac_version", data_dict)

        data_dict = phan.results_data(as_dict=True)
        self.assertIsInstance(data_dict, dict)

        data_json = phan.results_data(as_json=True)
        self.assertIsInstance(data_json, str)
        # shouldn't raise
        json.loads(data_json)

    def test_results_data_no_mtf(self):
        phan = LasVegas.from_demo_image()
        phan.analyze()

        data_dict = phan.results_data(as_dict=True)
        self.assertEqual(len(data_dict), 11)

    def test_set_figure_size(self):
        phan = LeedsTOR.from_demo_image()
        phan.analyze()
        phan.plot_analyzed_image(figsize=(7, 11))
        fig = plt.gcf()
        self.assertEqual(fig.bbox_inches.height, 11)
        self.assertEqual(fig.bbox_inches.width, 7)

    def test_set_figure_size_splot_plots(self):
        phan = LeedsTOR.from_demo_image()
        phan.analyze()
        figs, _ = phan.plot_analyzed_image(figsize=(7, 11), split_plots=True)
        self.assertEqual(figs[0].bbox_inches.height, 11)
        self.assertEqual(figs[0].bbox_inches.width, 7)

    def test_multiple_plots(self):
        phan = LeedsTOR.from_demo_image()
        phan.analyze()
        figs, names = phan.plot_analyzed_image(split_plots=True)
        self.assertEqual(len(figs), 3)
        files = phan.save_analyzed_image(filename="a.png", split_plots=True)
        names = ("a_image.png", "a_low_contrast.png", "a_high_contrast.png")
        for name in names:
            self.assertIn(name, files)

        # regular single plot produces one image/file
        figs, names = phan.plot_analyzed_image()
        self.assertEqual(len(figs), 0)
        name = "b.png"
        phan.save_analyzed_image("b.png")
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
        path = get_file_from_cloud_test_repo([TEST_DIR, "Leeds", "Leeds_ccw.dcm"])

        # do normal analysis
        phan = LeedsTOR(path)
        phan.analyze()
        x = phan.results_data().phantom_center_x_y[0]

        # pass kwarg; use same dpi as image; results should be the same.
        img = image.load(path)
        phan = LeedsTOR(path, image_kwargs={"dpi": img.dpi})
        phan.analyze()
        x_manual_dpi = phan.results_data().phantom_center_x_y[0]

        self.assertEqual(x, x_manual_dpi)

    def test_ssd_values(self):
        """Test various SSD values"""
        phan = LeedsTOR.from_demo_image()
        phan.analyze(ssd="auto")  # shouldn't raise
        phan = LeedsTOR.from_demo_image()
        phan.analyze(ssd=1000)  # shouldn't raise
        with self.assertRaises(ValueError):
            phan = LeedsTOR.from_demo_image()
            phan.analyze(ssd=1500)  # really at 1000

    def test_scaling_area(self):
        """Test various scaling area values"""
        phan = LeedsTOR.from_demo_image()
        phan.analyze()
        self.assertAlmostEqual(phan.results_data().phantom_area, 17760.9, delta=0.3)


class PlanarPhantomMixin(CloudFileMixin):
    klass: Callable
    dir_path = ["planar_imaging"]
    mtf_50 = None
    invert = False
    ssd = "auto"
    median_contrast = None
    median_cnr = None
    file_name = None
    rois_seen = None
    piu = None

    @classmethod
    def setUpClass(cls):
        cls.instance = cls.create_instance()
        cls.preprocess(cls.instance)
        cls.instance.analyze(ssd=cls.ssd, invert=cls.invert)

    @classmethod
    def create_instance(cls):
        if not cls.file_name:
            return cls.klass.from_demo_image()
        else:
            return cls.klass(cls.get_filename())

    @classmethod
    def preprocess(cls, instance):
        pass

    @classmethod
    def tearDownClass(cls):
        plt.close("all")
        del cls.instance

    def test_bad_inversion_recovers(self):
        instance = self.create_instance()
        instance.image.invert()
        instance.analyze(ssd=self.ssd, invert=self.invert)
        # check that the MTF is the expected value. This is a surrogate for the angle being wrong
        if self.mtf_50:
            self.assertAlmostEqual(
                self.mtf_50, instance.mtf.relative_resolution(50), delta=0.2
            )

    def test_percent_integral_uniformity(self):
        """Test the PIU if it exists in the results data"""
        if self.piu is not None:
            self.assertAlmostEqual(
                self.piu,
                self.instance.results_data().percent_integral_uniformity,
                delta=0.1,
            )
        elif not isinstance(self.instance.results_data(), LightRadResult):
            self.assertIsNone(self.instance.results_data().percent_integral_uniformity)

    def test_plotting(self):
        self.instance.plot_analyzed_image()
        self.instance.plot_analyzed_image(low_contrast=False, high_contrast=False)
        self.instance.plot_analyzed_image(
            image=False, low_contrast=False, high_contrast=False
        )

    def test_saving(self):
        self.instance.plot_analyzed_image()
        save_file(self.instance.save_analyzed_image)

    def test_pdf(self):
        save_file(self.instance.publish_pdf)

    def test_mtf(self):
        if self.mtf_50 is not None:
            self.assertAlmostEqual(
                self.mtf_50, self.instance.mtf.relative_resolution(50), delta=0.3
            )

    def test_rois_seen(self):
        if self.rois_seen is not None:
            self.assertEqual(
                self.rois_seen, self.instance.results_data().num_contrast_rois_seen
            )

    def test_median_contrast(self):
        if self.median_contrast is not None:
            self.assertAlmostEqual(
                self.median_contrast,
                self.instance.results_data().median_contrast,
                delta=0.03,
            )

    def test_median_cnr(self):
        if self.median_cnr is not None:
            self.assertAlmostEqual(
                self.median_cnr,
                self.instance.results_data().median_cnr,
                delta=0.01 * self.median_cnr,
            )

    def test_results(self):
        self.assertIsInstance(self.instance.results(), str)


class LeedsMixin(PlanarPhantomMixin):
    klass = LeedsTOR
    dir_path = ["planar_imaging", "Leeds"]


class LeedsDemo(LeedsMixin, TestCase):
    mtf_50 = 1.5
    piu = 91.8

    def test_demo(self):
        LeedsTOR.run_demo()  # shouldn't raise


class LeedsCCW(LeedsMixin, TestCase):
    mtf_50 = 1.5
    file_name = "Leeds_ccw.dcm"
    piu = 93.5


class Leeds45Deg(LeedsMixin, TestCase):
    mtf_50 = 1.9
    ssd = "auto"
    file_name = "Leeds-45deg.dcm"
    piu = 95.5


class LeedsDirtyEdges(LeedsMixin, TestCase):
    mtf_50 = 1.53
    ssd = "auto"
    file_name = "Leeds-dirty-edges.dcm"
    piu = 96.8


class LeedsOffsetHighRes(LeedsMixin, TestCase):
    mtf_50 = 1.85
    ssd = "auto"
    file_name = "Leeds_offset_high_res_rois.dcm"
    piu = 89.3


class LeedsBlue(LeedsMixin, TestCase):
    klass = LeedsTORBlue
    mtf_50 = 1.5
    ssd = "auto"
    file_name = "Leeds_Blue.dcm"
    piu = 97


class LeedsBlueRotated(LeedsMixin, TestCase):
    klass = LeedsTORBlue
    mtf_50 = 1.5
    ssd = "auto"
    file_name = "Leeds_Blue.dcm"
    piu = 97

    @classmethod
    def preprocess(cls, instance):
        instance.image.array = rotate(
            instance.image.array, angle=180, mode="mirror", order=0
        )


@skip("Phantom appears distorted. MTF locations are different than other phantoms")
class LeedsClosedBlades(LeedsMixin, TestCase):
    mtf_50 = 1.3
    ssd = "auto"
    file_name = "Leeds-closed-blades.dcm"


class LeedsACB1(LeedsMixin, TestCase):
    dir_path = ["planar_imaging", "Leeds", "ACB 1"]
    file_path = "1.dcm"
    mtf_50 = 1.69
    piu = 91.8


class LeedsBadInversion(LeedsMixin, TestCase):
    """Radmachine image where inversion was bad. pylinac should be able to correct"""

    file_path = "Leeds bad inversion.dcm"
    mtf_50 = 1.69
    piu = 91.8


class SIQC3Demo(PlanarPhantomMixin, TestCase):
    klass = StandardImagingQC3
    mtf_50 = 0.53
    rois_seen = 5
    piu = 98

    def test_demo(self):
        StandardImagingQC3.run_demo()  # shouldn't raise


class SIQC3_1(PlanarPhantomMixin, TestCase):
    klass = StandardImagingQC3
    file_name = "QC3-2.5MV.dcm"
    mtf_50 = 1.19
    rois_seen = 5
    piu = 91.8


class SIQC3_2(PlanarPhantomMixin, TestCase):
    klass = StandardImagingQC3
    file_name = "QC3-2.5MV-2.dcm"
    mtf_50 = 1.16
    ssd = 1000
    rois_seen = 5
    piu = 91.8

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


class LasVegasDemo(LasVegasTestMixin, TestCase):
    rois_seen = 12
    piu = 98.4

    def test_demo(self):
        LasVegas.run_demo()  # shouldn't raise


@skip(
    "Non-cardinal angles no longer supported. If support is re-added these can be reactivated"
)
class LasVegas10deg(LasVegasTestMixin, TestCase):
    file_path = [
        "TrueBeam 1 - 2364",
        "2.5MV LV HQ 10deg - ImageRT_2016-10-6 20-12-58.dcm",
    ]
    phantom_angle = 290


@skip(
    "Non-cardinal angles no longer supported. If support is re-added these can be reactivated"
)
class LasVegasrotated(LasVegasTestMixin, TestCase):
    file_path = [
        "TrueBeam 1 - 2364",
        "2.5MV LV HQ side2 - ImageRT_2016-10-6 20-43-3.dcm",
    ]
    phantom_angle = 284


@skip(
    "Non-cardinal angles no longer supported. If support is re-added these can be reactivated"
)
class LasVegasTB1(LasVegasTestMixin, TestCase):
    file_path = [
        "TrueBeam 1 - 2364",
        "6MV LasVegas HQ 0deg - ImageRT_2016-10-6 20-10-17.dcm",
    ]
    phantom_angle = 284.5


class ElektaLasVegasMixin(LasVegasTestMixin):
    dir_path = ["planar_imaging", "Elekta Las Vegas"]
    klass = ElektaLasVegas

    @classmethod
    def setUpClass(cls):
        cls.instance = cls.create_instance()
        cls.preprocess(cls.instance)
        cls.instance.image.rot90(n=3)
        cls.instance.analyze(ssd=cls.ssd, invert=cls.invert)


class ElektaDemo(ElektaLasVegasMixin, TestCase):
    rois_seen = 17
    piu = 99.5

    def test_demo(self):
        ElektaLasVegas.run_demo()  # shouldn't raise


class Elekta2MU(ElektaLasVegasMixin, TestCase):
    file_name = "LasVegas_2MU.dcm"
    rois_seen = 12
    piu = 99.11


class Elekta10MU(ElektaLasVegasMixin, TestCase):
    file_name = "LasVegas_10MU.dcm"
    rois_seen = 17
    piu = 99.4


class DoselabMVDemo(PlanarPhantomMixin, TestCase):
    klass = DoselabMC2MV
    mtf_50 = 0.54
    piu = 51.8

    def test_demo(self):
        DoselabMC2MV.run_demo()


class DoselabkVDemo(PlanarPhantomMixin, TestCase):
    klass = DoselabMC2kV
    mtf_50 = 2.0
    piu = 4.2

    def test_demo(self):
        DoselabMC2kV.run_demo()


class DoselabkV70kVp(PlanarPhantomMixin, TestCase):
    klass = DoselabMC2kV
    dir_path = ["planar_imaging", "Doselab MC2"]
    file_name = "DL kV 70kVp.dcm"
    mtf_50 = 1.14
    piu = 0


class SNCkVDemo(PlanarPhantomMixin, TestCase):
    klass = SNCkV
    mtf_50 = 1.76
    median_contrast = 0.17
    median_cnr = 69.4
    piu = 98.7

    def test_demo(self):
        SNCkV.run_demo()


class SNCMVDemo(PlanarPhantomMixin, TestCase):
    klass = SNCMV
    median_cnr = 81
    median_contrast = 0.21
    mtf_50 = 0.43
    piu = 98.4

    def test_demo(self):
        SNCMV.run_demo()


class SNCMV12510_6MV1(PlanarPhantomMixin, TestCase):
    klass = SNCMV12510
    mtf_50 = 0.91
    median_contrast = 0.254
    median_cnr = 65.34
    dir_path = ["planar_imaging", "SNC MV Old"]
    file_name = "SNC_MV_Old1.dcm"
    piu = 96.3

    def test_demo(self):
        SNCMV12510.run_demo()


class SNCMV12510_6MV2(PlanarPhantomMixin, TestCase):
    klass = SNCMV12510
    mtf_50 = 0.85
    median_contrast = 0.255
    median_cnr = 66.43
    dir_path = ["planar_imaging", "SNC MV Old"]
    file_name = "SNC_MV_Old2.dcm"
    piu = 96.4


class SNCMV12510_Jig(PlanarPhantomMixin, TestCase):
    """Phantom where the jig is touching and gets in the way of analysis"""

    klass = SNCMV12510
    mtf_50 = 0.92
    median_contrast = 0.23
    median_cnr = 58.6
    dir_path = ["planar_imaging", "SNC MV Old"]
    file_name = "SNC_MV_jig.dcm"
    piu = 96.7


class IBAPrimusDemo(PlanarPhantomMixin, TestCase):
    klass = IBAPrimusA
    dir_path = ["planar_imaging", "PrimusL"]
    file_name = "Demo.dcm"
    mtf_50 = 1.66
    ssd = 1395
    median_cnr = 1084.4
    median_contrast = 0.62
    piu = 90.1

    def test_demo(self):
        IBAPrimusA.run_demo()


class IBAPrimusBasic(IBAPrimusDemo):
    # same as demo but no test_demo method; this is inherited so no need to call test_demo a lot

    def test_demo(self):
        pass


class IBAPrimusDemo0(IBAPrimusBasic):
    """Rotate image to 0 (pointing towards gun) to ensure it still analyzes and results are similar"""

    @classmethod
    def preprocess(cls, instance):
        instance.image.rot90()


class IBAPrimusShifted(IBAPrimusBasic):
    """Shift the image slightly to ensure we can handle slightly offset phantom placements"""

    @classmethod
    def preprocess(cls, instance):
        instance.image.array = np.roll(instance.image.array, shift=50)


class IBAPrimusDemoMinus90(IBAPrimusBasic):
    """Rotate image to -90 (pointing left in BEV) to ensure it still analyzes and results are similar"""

    @classmethod
    def preprocess(cls, instance):
        instance.image.rot90(2)


class IBAPrimusDemoBadInversion(IBAPrimusBasic):
    """Force a bad inversion and ensure recovery"""

    @classmethod
    def preprocess(cls, instance):
        instance.image.invert()


class IBAPrimusFarSSD(PlanarPhantomMixin, TestCase):
    klass = IBAPrimusA
    dir_path = ["planar_imaging", "PrimusL"]
    file_name = "Primus_farSSD.dcm"
    mtf_50 = 2.33
    ssd = 2790
    median_cnr = 3990
    median_contrast = 0.6
    piu = 89.8


class SIQCkVDemo(PlanarPhantomMixin, TestCase):
    klass = StandardImagingQCkV
    mtf_50 = 1.81
    rois_seen = 5
    piu = 96.7

    def test_demo(self):
        StandardImagingQCkV.run_demo()


class PTWEPIDDemo(PlanarPhantomMixin, TestCase):
    klass = PTWEPIDQC
    mtf_50 = 0.79
    piu = 95.4

    def test_demo(self):
        PTWEPIDQC.run_demo()


class PTWEPIDQC1(PlanarPhantomMixin, TestCase):
    klass = PTWEPIDQC
    mtf_50 = 0.79
    rois_seen = 9
    median_contrast = 0.26
    median_cnr = 40.9
    dir_path = ["planar_imaging", "PTW-EPID"]
    file_name = "PTW EPID QC Phantom.dcm"
    piu = 94.7


class PTWEPID15MV(PlanarPhantomMixin, TestCase):
    klass = PTWEPIDQC
    mtf_50 = 0.5
    rois_seen = 9
    median_contrast = 0.17
    median_cnr = 26.7
    dir_path = ["planar_imaging", "PTW-EPID"]
    file_name = "PTW-EPID-15MV.dcm"
    piu = 96.3


class PTWEPID6xHigh(PlanarPhantomMixin, TestCase):
    klass = PTWEPIDQC
    mtf_50 = 0.79
    rois_seen = 9
    median_contrast = 0.28
    median_cnr = 72.1
    dir_path = ["planar_imaging", "PTW-EPID"]
    file_name = "QA EPI 6x High.dcm"
    piu = 94.3


class PTWEPID6xHighQuality(PlanarPhantomMixin, TestCase):
    klass = PTWEPIDQC
    mtf_50 = 0.79
    rois_seen = 9
    median_contrast = 0.254
    median_cnr = 37.9
    dir_path = ["planar_imaging", "PTW-EPID"]
    file_name = "TB1_PTW_EPID_Phan 6x_HighQuality.dcm"
    piu = 94.6


class PTWEPIDTB3(PlanarPhantomMixin, TestCase):
    klass = PTWEPIDQC
    mtf_50 = 0.79
    rois_seen = 9
    median_contrast = 0.31
    median_cnr = 43.2
    dir_path = ["planar_imaging", "PTW-EPID"]
    file_name = "TB3_PTW_EPID_Phan.dcm"
    piu = 92.3


class PTWEPIDTB4(PlanarPhantomMixin, TestCase):
    klass = PTWEPIDQC
    mtf_50 = 0.79
    rois_seen = 9
    median_contrast = 0.30
    median_cnr = 39.1
    dir_path = ["planar_imaging", "PTW-EPID"]
    file_name = "TB4 PTW_EPID_phan.dcm"
    piu = 92.2


class FC2Mixin(PlanarPhantomMixin):
    klass = StandardImagingFC2
    dir_path = ["planar_imaging", "SI FC2"]
    field_size_x_mm = 150
    field_size_y_mm = 150
    field_epid_offset_x_mm = 0
    field_epid_offset_y_mm = 0
    field_bb_offset_x_mm = 0
    field_bb_offset_y_mm = 0
    fwxm = 50

    @classmethod
    def setUpClass(cls):
        cls.instance = cls.create_instance()
        cls.instance.analyze(invert=cls.invert, fwxm=cls.fwxm)

    def test_bad_inversion_recovers(self):
        # no inversion issues w/ this phantom
        pass

    def test_plotting(self):
        self.instance.plot_analyzed_image()

    def test_field_size(self):
        results_data = self.instance.results_data()
        self.assertAlmostEqual(
            results_data.field_size_x_mm, self.field_size_x_mm, delta=0.3
        )
        self.assertAlmostEqual(
            results_data.field_size_y_mm, self.field_size_y_mm, delta=0.3
        )

    def test_epid_offset(self):
        results_data = self.instance.results_data()
        self.assertAlmostEqual(
            results_data.field_epid_offset_x_mm, self.field_epid_offset_x_mm, delta=0.2
        )
        self.assertAlmostEqual(
            results_data.field_epid_offset_y_mm, self.field_epid_offset_y_mm, delta=0.2
        )

    def test_bb_offset(self):
        results_data = self.instance.results_data()
        self.assertAlmostEqual(
            results_data.field_bb_offset_x_mm, self.field_bb_offset_x_mm, delta=0.2
        )
        self.assertAlmostEqual(
            results_data.field_bb_offset_y_mm, self.field_bb_offset_y_mm, delta=0.2
        )


class FC2Demo(FC2Mixin, TestCase):
    field_size_y_mm = 148.5
    field_size_x_mm = 149.1
    field_epid_offset_y_mm = -0.7
    field_epid_offset_x_mm = 0.3
    field_bb_offset_x_mm = 0.1
    field_bb_offset_y_mm = -0.2

    def test_demo(self):
        StandardImagingFC2.run_demo()


class FC210x10_10FFF(FC2Mixin, TestCase):
    file_name = "FC-2-10x10-10fff.dcm"
    field_size_y_mm = 98.7
    field_size_x_mm = 99.3
    field_epid_offset_x_mm = 0
    field_epid_offset_y_mm = 0.3
    field_bb_offset_y_mm = 1
    field_bb_offset_x_mm = -0.3


class FC210x10_10X(FC2Mixin, TestCase):
    file_name = "FC-2-10x10-10x.dcm"
    field_size_y_mm = 99.3
    field_size_x_mm = 99.6
    field_epid_offset_x_mm = 0.2
    field_epid_offset_y_mm = -0.5
    field_bb_offset_y_mm = 0.4
    field_bb_offset_x_mm = -0.1


class FC210x10_15X(FC2Mixin, TestCase):
    file_name = "FC-2-10x10-15x.dcm"
    field_size_y_mm = 99.3
    field_size_x_mm = 99.6
    field_epid_offset_x_mm = 0.1
    field_epid_offset_y_mm = -0.5
    field_bb_offset_y_mm = 0.5
    field_bb_offset_x_mm = -0.2


class FC215x15_10X(FC2Mixin, TestCase):
    file_name = "FC-2-15x15-10X.dcm"
    field_size_x_mm = 149.2
    field_size_y_mm = 149.2
    field_epid_offset_x_mm = 0.1
    field_epid_offset_y_mm = -0.5
    field_bb_offset_y_mm = 0.25
    field_bb_offset_x_mm = 0


class FC215x15_10FFF(FC2Mixin, TestCase):
    file_name = "FC-2-15x15-10XFFF.dcm"
    fwxm = 30
    field_size_x_mm = 149.5
    field_size_y_mm = 149.6
    field_epid_offset_x_mm = -0.1
    field_epid_offset_y_mm = 0.2
    field_bb_offset_y_mm = 1.0
    field_bb_offset_x_mm = -0.3


class FC2Yoda(FC2Mixin, TestCase):
    file_name = "FC-2-Yoda.dcm"
    field_size_y_mm = 148.2
    field_size_x_mm = 149.2
    field_epid_offset_y_mm = -1
    field_epid_offset_x_mm = 0.2
    field_bb_offset_y_mm = -1
    field_bb_offset_x_mm = 0.5


class FC2Perfect(FC2Mixin, TestCase):
    file_name = "fc2-perfect.dcm"
    field_size_y_mm = 120.3
    field_size_x_mm = 120.3
    field_epid_offset_y_mm = 0
    field_epid_offset_x_mm = 0
    field_bb_offset_y_mm = 0
    field_bb_offset_x_mm = 0


class FC2FieldDown1mm(FC2Mixin, TestCase):
    file_name = "fc2-down1mm.dcm"
    field_size_y_mm = 120.3
    field_size_x_mm = 120.3
    field_epid_offset_y_mm = -1.0
    field_epid_offset_x_mm = 0
    field_bb_offset_y_mm = -1.0
    field_bb_offset_x_mm = 0


class FC2BBDownRight1mm(FC2Mixin, TestCase):
    file_name = "fc2-bbdownright1mm.dcm"
    field_size_y_mm = 120.3
    field_size_x_mm = 120.3
    field_epid_offset_y_mm = 0
    field_epid_offset_x_mm = 0
    field_bb_offset_y_mm = 1
    field_bb_offset_x_mm = 1


class DoselabRLfMixin(FC2Mixin):
    klass = DoselabRLf
    dir_path = ["planar_imaging", "Doselab RLf"]


class DoselabRLfDemo(DoselabRLfMixin, TestCase):
    field_size_y_mm = 148.2
    field_size_x_mm = 149.3
    field_epid_offset_x_mm = 0.1
    field_epid_offset_y_mm = 0.6
    field_bb_offset_y_mm = 0.8
    field_bb_offset_x_mm = 0.2

    def test_demo(self):
        DoselabRLf.run_demo()


class DoselabRLf10x10(DoselabRLfMixin, TestCase):
    file_name = "FS 10x10.dcm"
    field_size_y_mm = 98.2
    field_size_x_mm = 99.2
    field_epid_offset_x_mm = 0.2
    field_epid_offset_y_mm = 0.8
    field_bb_offset_y_mm = 1.1
    field_bb_offset_x_mm = 0.2


class DoselabRLfKB(DoselabRLfMixin, TestCase):
    """A failed dataset from KB."""

    file_name = "KB_RLf.dcm"
    field_size_y_mm = 98.9
    field_size_x_mm = 100.8
    field_epid_offset_x_mm = -0.2
    field_epid_offset_y_mm = 0.8
    field_bb_offset_y_mm = -0.1
    field_bb_offset_x_mm = -0.3


class DoselabRLfKB2(DoselabRLfMixin, TestCase):
    """A failed dataset from KB."""

    file_name = "KB2_RLf.dcm"
    field_size_y_mm = 98.2
    field_size_x_mm = 99.2
    field_epid_offset_x_mm = 0.2
    field_epid_offset_y_mm = 0.8
    field_bb_offset_y_mm = 1
    field_bb_offset_x_mm = 0.2


class IsoAlignMixin(FC2Mixin):
    klass = IsoAlign
    dir_path = ["planar_imaging", "Doselab RLf"]


class IsoAlignDemo(IsoAlignMixin, TestCase):
    field_size_y_mm = 149.6
    field_size_x_mm = 150.2
    field_epid_offset_y_mm = 0.4
    field_epid_offset_x_mm = -0.1
    field_bb_offset_y_mm = 0.0
    field_bb_offset_x_mm = 0.2

    def test_demo(self):
        IsoAlign.run_demo()


class IMTLRadMixin(FC2Mixin):
    klass = IMTLRad
    dir_path = ["planar_imaging", "IMT_L-Rad"]


class IMTLRadDemo(IMTLRadMixin, TestCase):
    field_size_y_mm = 210.8
    field_size_x_mm = 210.5
    field_epid_offset_x_mm = 1.3
    field_epid_offset_y_mm = 0
    field_bb_offset_y_mm = 0.8
    field_bb_offset_x_mm = -0.6

    def test_demo(self):
        IMTLRad.run_demo()


class IMTLRad21x21(IMTLRadMixin, TestCase):
    file_name = "RTIMAGE_11_1.dcm"
    field_size_y_mm = 210.8
    field_size_x_mm = 210.5
    field_epid_offset_x_mm = 1.3
    field_epid_offset_y_mm = 0
    field_bb_offset_y_mm = 0.8
    field_bb_offset_x_mm = -0.6


class IMTLRad21x21_2(IMTLRadMixin, TestCase):
    file_name = "RTIMAGE_11_2.dcm"
    field_size_y_mm = 210.9
    field_size_x_mm = 210.7
    field_epid_offset_x_mm = 1.4
    field_epid_offset_y_mm = -0.3
    field_bb_offset_y_mm = 0.4
    field_bb_offset_x_mm = -0.7


class IMTLRadPerfect(IMTLRadMixin, TestCase):
    file_name = "perfect_imt.dcm"
    field_size_y_mm = 150
    field_size_x_mm = 150
    field_epid_offset_x_mm = 0
    field_epid_offset_y_mm = 0
    field_bb_offset_y_mm = 0
    field_bb_offset_x_mm = 0


class IMTLRadOffset(IMTLRadMixin, TestCase):
    file_name = "offset_imt.dcm"
    field_size_y_mm = 150
    field_size_x_mm = 150
    field_epid_offset_x_mm = 0
    field_epid_offset_y_mm = 0
    field_bb_offset_y_mm = 2
    field_bb_offset_x_mm = 0


class SNCFSQAMixin(FC2Mixin):
    klass = SNCFSQA
    dir_path = ["planar_imaging", "SNC FSQA"]


class SNCFSQADemo(SNCFSQAMixin, TestCase):
    # The 6x 15x15 file in the test repo is the same as the demo
    field_size_y_mm = 149
    field_size_x_mm = 150.4
    field_epid_offset_x_mm = 0.7
    field_epid_offset_y_mm = 0.3
    field_bb_offset_y_mm = -0.1
    field_bb_offset_x_mm = -0.4

    def test_demo(self):
        SNCFSQA.run_demo()


class SNCFSQA10x10(SNCFSQAMixin, TestCase):
    file_name = "6x_FSQA_10x10.dcm"
    field_size_y_mm = 99.6
    field_size_x_mm = 100.2
    field_epid_offset_x_mm = 0.7
    field_epid_offset_y_mm = -0.1
    field_bb_offset_y_mm = -0.5
    field_bb_offset_x_mm = -0.5
