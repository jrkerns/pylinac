import io
import json
import tempfile
from functools import partial
from pathlib import Path
from typing import Iterable, Type, Union
from unittest import TestCase

import numpy as np
from matplotlib import pyplot as plt

from pylinac import DRGS, DRMLC
from pylinac.core.geometry import Point
from pylinac.core.image_generator import (
    AS1200Image,
    FilterFreeFieldLayer,
    GaussianFilterLayer,
    RandomNoiseLayer,
)
from pylinac.vmat import VMATResult
from tests_basic.utils import (
    FromDemoImageTesterMixin,
    FromURLTesterMixin,
    get_file_from_cloud_test_repo,
    save_file,
)

TEST_DIR = "VMAT"

within_5 = partial(TestCase().assertAlmostEqual, delta=5)
within_1 = partial(TestCase().assertAlmostEqual, delta=1)


class LoadingBase(FromURLTesterMixin, FromDemoImageTesterMixin):
    demo_load_method = "from_demo_images"
    klass: Union[Type[DRGS], Type[DRMLC]]

    def test_normal_instantiation(self):
        one = get_file_from_cloud_test_repo([TEST_DIR, "no_test_or_image_type_1.dcm"])
        two = get_file_from_cloud_test_repo([TEST_DIR, "no_test_or_image_type_2.dcm"])
        instance = self.klass(image_paths=(one, two))
        self.assertIsInstance(instance, self.klass)

    def test_from_stream(self):
        one = get_file_from_cloud_test_repo([TEST_DIR, "no_test_or_image_type_1.dcm"])
        two = get_file_from_cloud_test_repo([TEST_DIR, "no_test_or_image_type_2.dcm"])
        with open(one, "rb") as s1, open(two, "rb") as s2:
            s11 = io.BytesIO(s1.read())
            s22 = io.BytesIO(s2.read())
            instance = self.klass(image_paths=(s11, s22))
            instance.analyze()
        self.assertIsInstance(instance, self.klass)

    def test_from_file_object(self):
        one = get_file_from_cloud_test_repo([TEST_DIR, "no_test_or_image_type_1.dcm"])
        two = get_file_from_cloud_test_repo([TEST_DIR, "no_test_or_image_type_2.dcm"])
        with open(one, "rb") as s1, open(two, "rb") as s2:
            instance = self.klass(image_paths=(s1, s2))
            instance.analyze()
        self.assertIsInstance(instance, self.klass)

    def test_passing_3_images_fails(self):
        """Test passing the wrong number of images."""
        with self.assertRaises(ValueError):
            self.klass(image_paths=("", "", ""))

    def test_print_results(self):
        instance = self.klass.from_demo_images()
        instance.analyze()
        self.assertIsInstance(instance.results(), str)

    def test_plot_analyzed_image(self):
        instance = self.klass.from_demo_images()
        instance.analyze()
        instance.plot_analyzed_image()  # shouldn't raise

    def test_publish_pdf(self):
        instance = self.klass.from_demo_images()
        instance.analyze()
        save_file(instance.publish_pdf)

    def test_results_data(self):
        instance = self.klass.from_demo_images()
        instance.analyze()
        data = instance.results_data()
        self.assertIsInstance(data, VMATResult)
        self.assertEqual(data.test_type, instance._result_header)
        data_dict = instance.results_data(as_dict=True)
        self.assertEqual(len(data_dict), 9)

        data_dict = instance.results_data(as_dict=True)
        self.assertIsInstance(data_dict, dict)
        self.assertIn("pylinac_version", data_dict)
        self.assertEqual(data_dict["max_deviation_percent"], instance.max_r_deviation)

        data_str = instance.results_data(as_json=True)
        self.assertIsInstance(data_str, str)
        # shouldn't raise
        json.loads(data_str)

    def test_custom_roi_config(self):
        my_drgs = DRGS.from_demo_images()
        my_drgs.analyze(roi_config={"DR: 150 MU/min": {"offset_mm": 0}})
        self.assertEqual(len(my_drgs.segments), 1)
        results_data = my_drgs.results_data()
        self.assertIn("DR: 150 MU/min", results_data.named_segment_data.keys())

    def test_custom_num_rois_and_spacing(self):
        """Kraken minimizes the number of inputs; accepts the # of ROIs and spacing. Essentially the same implementation as Kraken"""
        num_roi = 5
        spacing_mm = 30
        offsets = np.arange(0, num_roi * spacing_mm, 30)
        centered_offsets = offsets - np.mean(offsets)
        roi_config = {
            f"ROI {idx + 1}": {"offset_mm": offset}
            for idx, offset in enumerate(centered_offsets)
        }
        my_drgs = DRGS.from_demo_images()
        my_drgs.analyze(roi_config=roi_config)
        self.assertEqual(len(my_drgs.segments), 5)


class TestDRGSLoading(LoadingBase, TestCase):
    url = "drgs.zip"
    klass = DRGS


class TestDRMLCLoading(LoadingBase, TestCase):
    url = "drmlc.zip"
    klass = DRMLC


class VMATMixin:
    klass = object
    filepaths = Iterable[str]
    is_zip = False
    segment_positions = {1: Point(100, 200)}
    segment_values = {
        0: {"r_dev": 0, "r_corr": 100},
        4: {"r_dev": 0, "r_corr": 100},
    }
    init_kwargs = {}
    analyze_kwargs = {}
    avg_abs_r_deviation = 0
    avg_r_deviation = 0
    max_r_deviation = 0
    passes = True
    print_debug = False

    @classmethod
    def absolute_path(cls):
        if cls.is_zip:
            path = get_file_from_cloud_test_repo([TEST_DIR, *cls.filepaths])
        else:
            path = [
                get_file_from_cloud_test_repo([TEST_DIR, path])
                for path in cls.filepaths
            ]
        return path

    def setUp(self):
        if self.is_zip:
            self.vmat = self.klass.from_zip(self.absolute_path(), **self.init_kwargs)
        else:
            self.vmat = self.klass(self.absolute_path(), **self.init_kwargs)
        self.vmat.analyze(**self.analyze_kwargs)
        if self.print_debug:
            print(self.vmat.results())
            print(
                f"Segment 0: rdev {self.vmat.segments[0].r_dev:2.3f}, rcorr {self.vmat.segments[0].r_corr:2.3f}"
            )
            if self.klass == DRGS:
                print(
                    f"Segment 4: rdev {self.vmat.segments[4].r_dev:2.3f}, rcorr {self.vmat.segments[4].r_corr:2.3f}"
                )
            else:
                print(
                    f"Segment 2: rdev {self.vmat.segments[2].r_dev:2.3f}, rcorr {self.vmat.segments[2].r_corr:2.3f}"
                )
            print("Max dev", self.vmat.max_r_deviation)

    def test_overall_passed(self):
        self.assertEqual(self.vmat.passed, self.passes)

    def test_fail_with_tight_tolerance(self):
        self.vmat.analyze(tolerance=0.001)
        self.assertFalse(self.vmat.passed)

    def test_segment_positions(self):
        for key, value in self.segment_positions.items():
            within_5(self.vmat.segments[key].center.x, value.x)
            within_5(self.vmat.segments[key].center.y, value.y)

    def test_segment_values(self):
        for key, value in self.segment_values.items():
            within_1(self.vmat.segments[key].r_dev, value["r_dev"])
            within_1(self.vmat.segments[key].r_corr, value["r_corr"])
            if "stdev" in value:
                self.assertAlmostEqual(
                    self.vmat.segments[key].stdev, value["stdev"], places=2
                )

    def test_deviations(self):
        self.assertAlmostEqual(
            self.vmat.avg_abs_r_deviation, self.avg_abs_r_deviation, delta=0.05
        )
        self.assertAlmostEqual(
            self.vmat.avg_r_deviation, self.avg_r_deviation, delta=0.02
        )
        self.assertAlmostEqual(
            self.vmat.max_r_deviation, self.max_r_deviation, delta=0.1
        )

    def test_different_segment_size_is_nearly_the_same(self):
        segment_width_mm = 10
        segment_height_mm = 50
        self.vmat.analyze(segment_size_mm=(segment_width_mm, segment_height_mm))
        self.assertTrue(self.vmat.segments[0]._nominal_width_mm, segment_width_mm)
        self.assertTrue(self.vmat.segments[0]._nominal_height_mm, segment_height_mm)


class TestDRGSDemo(VMATMixin, TestCase):
    """Tests of the result values of the DRGS demo images."""

    segment_positions = {0: Point(161, 192), 4: Point(314, 192)}
    segment_values = {
        0: {"r_dev": 0.965, "r_corr": 6.2, "stdev": 0.0008},
        4: {"r_dev": -0.459, "r_corr": 6, "stdev": 0.0007},
    }
    avg_abs_r_deviation = 0.74
    max_r_deviation = 1.8
    passes = False

    def setUp(self):
        self.vmat = DRGS.from_demo_images()
        self.vmat.analyze()

    def test_demo(self):
        """Run the demo; no errors should arise."""
        self.vmat.run_demo()

    def test_set_figure_size(self):
        self.vmat.plot_analyzed_image(figsize=(7, 11))
        fig = plt.gcf()
        self.assertEqual(fig.bbox_inches.height, 11)
        self.assertEqual(fig.bbox_inches.width, 7)


class TestDRMLCDemo(VMATMixin, TestCase):
    """Tests of the result values of the DRMLC demo images."""

    segment_positions = {0: Point(170, 192), 2: Point(285, 192)}
    segment_values = {
        0: {"r_dev": -0.7, "r_corr": 5.7, "stdev": 0.00086},
        2: {"r_dev": -0.405, "r_corr": 5.8, "stdev": 0.00085},
    }
    avg_abs_r_deviation = 0.44
    max_r_deviation = 0.89

    def setUp(self):
        self.vmat = DRMLC.from_demo_images(**self.init_kwargs)
        self.vmat.analyze()

    def test_demo(self):
        self.vmat.run_demo()


class TestDRMLCDemoRawPixels(TestDRMLCDemo):
    """Use raw DICOM pixel values, like doselab does."""

    init_kwargs = {"raw_pixels": True, "ground": False, "check_inversion": False}
    segment_values = {
        0: {"r_dev": -0.55, "r_corr": 138.55},
        2: {"r_dev": 0.56, "r_corr": 140},
    }
    avg_abs_r_deviation = 0.54
    max_r_deviation = 0.56


class TestDRMLC105(VMATMixin, TestCase):
    """Tests of the result values of MLCS images at 105cm SID."""

    klass = DRMLC
    filepaths = ("DRMLCopen-105-example.dcm", "DRMLCdmlc-105-example.dcm")
    segment_positions = {0: Point(391, 384), 2: Point(552, 384)}
    segment_values = {
        0: {
            "r_dev": -2.1,
            "r_corr": 13.66,
        },  # r_corr changed in v3.0 due to difference in default inversion/scaling
        2: {"r_dev": 0.22, "r_corr": 14},
    }
    avg_abs_r_deviation = 1.06
    max_r_deviation = 2.11
    passes = False


class TestDRGS105(VMATMixin, TestCase):
    """Tests of the result values of DRMLC images at 105cm SID."""

    filepaths = ("DRGSopen-105-example.dcm", "DRGSdmlc-105-example.dcm")
    klass = DRGS
    segment_positions = {0: Point(378, 384), 2: Point(485, 384)}
    segment_values = {
        0: {
            "r_dev": 1.385,
            "r_corr": 15.11,
        },  # r_corr changed in v3.0 due to difference in default inversion/scaling
        4: {"r_dev": -0.8, "r_corr": 14.8},
    }
    avg_abs_r_deviation = 0.68
    max_r_deviation = 1.38


class TestDRMLC2(VMATMixin, TestCase):
    """Tests of the result values of MLCS images at 105cm SID."""

    filepaths = ("DRMLC-2_open.dcm", "DRMLC-2_dmlc.dcm")
    klass = DRMLC
    segment_positions = {0: Point(199, 192), 2: Point(275, 192)}
    segment_values = {
        0: {"r_dev": 0.77, "r_corr": 6.1},
        2: {"r_dev": -1.1, "r_corr": 6},
    }
    avg_abs_r_deviation = 1.4
    max_r_deviation = 2.11
    passes = False


class TestDRGS2(VMATMixin, TestCase):
    """Tests of the result values of DRMLC images at 105cm SID."""

    filepaths = ("DRGS-2_open.dcm", "DRGS-2_dmlc.dcm")
    klass = DRGS
    segment_positions = {0: Point(191, 192), 2: Point(242, 192)}
    segment_values = {
        0: {"r_dev": 1.5, "r_corr": 6.4},
        4: {"r_dev": -0.7, "r_corr": 6.3},
    }
    avg_abs_r_deviation = 0.7
    max_r_deviation = 1.5


class TestDRMLCWideGaps(VMATMixin, TestCase):
    """Tests of the result values of a perfect DRMLC but with very wide gaps."""

    filepaths = ("vmat-drgs-open-wide-gaps.dcm", "vmat-drgs-open-wide-gaps.dcm")
    klass = DRMLC
    segment_positions = {0: Point(439, 640), 2: Point(707, 640)}
    segment_values = {
        0: {"r_dev": 0, "r_corr": 100},
        2: {"r_dev": 0, "r_corr": 100},
    }
    avg_abs_r_deviation = 0
    max_r_deviation = 0.0
    passes = True

    def test_fail_with_tight_tolerance(self):
        pass


class TestDRMLCOverlapGaps(VMATMixin, TestCase):
    """Tests of the result values of a perfect DRMLC but with gaps that are overlapping (e.g. from a poor DLG)."""

    filepaths = ("vmat-drgs-open-overlap.dcm", "vmat-drgs-open-overlap.dcm")
    klass = DRMLC
    segment_positions = {0: Point(439, 640), 2: Point(707, 640)}
    segment_values = {
        0: {"r_dev": 0, "r_corr": 100},
        2: {"r_dev": 0, "r_corr": 100},
    }
    avg_abs_r_deviation = 0
    max_r_deviation = 0.0
    passes = True

    def test_fail_with_tight_tolerance(self):
        pass


class TestHalcyonDRGS(VMATMixin, TestCase):
    """A Halcyon image is FFF and goes to the edge of the EPID. Causes bad inversion w/o FWXM profile type."""

    klass = DRGS
    filepaths = ("HalcyonDRGS.zip",)
    analyze_kwargs = {
        "roi_config": {
            "ROI 1": {"offset_mm": -120},
            "ROI 2": {"offset_mm": -80},
            "ROI 3": {"offset_mm": -40},
            "ROI 4": {"offset_mm": 0},
            "ROI 5": {"offset_mm": 40},
            "ROI 6": {"offset_mm": 80},
            "ROI 7": {"offset_mm": 120},
        }
    }
    is_zip = True
    segment_positions = {0: Point(89, 640), 2: Point(456, 640)}
    segment_values = {
        0: {"r_dev": 0.583, "r_corr": 13.62},
        2: {"r_dev": -0.30, "r_corr": 13.5},
    }
    avg_abs_r_deviation = 0.266
    max_r_deviation = 0.58


class TestHalcyonDRMLC(VMATMixin, TestCase):
    """A Halcyon image is FFF and goes to the edge of the EPID. Causes bad inversion w/o FWXM profile type."""

    klass = DRMLC
    filepaths = ("HalcyonDRMLC.zip",)
    is_zip = True
    analyze_kwargs = {
        "roi_config": {
            "ROI 1": {"offset_mm": -115},
            "ROI 2": {"offset_mm": -57.5},
            "ROI 3": {"offset_mm": 0},
            "ROI 4": {"offset_mm": 57.5},
            "ROI 5": {"offset_mm": 115},
        }
    }
    segment_positions = {0: Point(112, 640), 2: Point(639, 640)}
    segment_values = {
        0: {"r_dev": -0.34, "r_corr": 2.8},
        2: {"r_dev": 1.15, "r_corr": 2.85},
    }
    avg_abs_r_deviation = 0.66
    max_r_deviation = 1.15
    passes = True


class TestHalcyonDRGS2(VMATMixin, TestCase):
    """A Hal image w/ deep gaps between the ROIs. Causes a shift in the ROIs from RAM-3483"""

    klass = DRGS
    filepaths = ("DRGS_Halcyon2.zip",)
    is_zip = True
    segment_positions = {0: Point(364, 640), 2: Point(543, 640)}
    segment_values = {
        0: {"r_dev": 1.17, "r_corr": 13.73},
        2: {"r_dev": -0.206, "r_corr": 13.56},
    }
    avg_abs_r_deviation = 0.44
    max_r_deviation = 0.803
    passes = True


class TestHalcyonDRGS3(VMATMixin, TestCase):
    """A TB image w/ deep gaps between the ROIs. Causes a shift in the ROIs from RAM-3483"""

    klass = DRGS
    filepaths = ("DRGS_example_PM.zip",)
    is_zip = True
    segment_positions = {0: Point(280, 384), 2: Point(433, 384)}
    segment_values = {
        0: {"r_dev": -0.37, "r_corr": 13.73},
        2: {"r_dev": -0.206, "r_corr": 13.56},
    }
    avg_abs_r_deviation = 0.89
    max_r_deviation = 1.87
    passes = False


class TestContrivedWideGapTest(VMATMixin, TestCase):
    """A contrived test with a wide gap between the segments."""

    klass = DRMLC
    is_zip = False
    segment_positions = {0: Point(506, 640), 2: Point(685, 640)}
    segment_values = {
        0: {"r_dev": 0, "r_corr": 100},
        2: {"r_dev": 0, "r_corr": 100},
    }
    avg_abs_r_deviation = 0
    max_r_deviation = 0.0
    passes = True

    def create_synthetic_images(self):
        tmp_dir = Path(tempfile.gettempdir())
        as1200_open = AS1200Image(1000)
        as1200_open.add_layer(FilterFreeFieldLayer(field_size_mm=(110, 110)))
        as1200_open.add_layer(GaussianFilterLayer())
        open_path = tmp_dir / "contrived_wide_gap_open.dcm"
        as1200_open.generate_dicom(open_path)
        # generate the DMLC image
        as1200_dmlc = AS1200Image(1000)
        as1200_dmlc.add_layer(
            FilterFreeFieldLayer(field_size_mm=(150, 20), cax_offset_mm=(0, 45))
        )
        as1200_dmlc.add_layer(
            FilterFreeFieldLayer(field_size_mm=(150, 20), cax_offset_mm=(0, 15))
        )
        as1200_dmlc.add_layer(
            FilterFreeFieldLayer(field_size_mm=(150, 20), cax_offset_mm=(0, -15))
        )
        as1200_dmlc.add_layer(
            FilterFreeFieldLayer(field_size_mm=(150, 20), cax_offset_mm=(0, -45))
        )
        as1200_dmlc.add_layer(GaussianFilterLayer())
        as1200_dmlc.add_layer(RandomNoiseLayer(sigma=0.005))
        dmlc_path = tmp_dir / "contrived_wide_gap_dmlc.dcm"
        as1200_dmlc.generate_dicom(dmlc_path)
        return open_path, dmlc_path

    def setUp(self):
        open_path, dmlc_path = self.create_synthetic_images()
        self.vmat = self.klass(image_paths=(open_path, dmlc_path), **self.init_kwargs)
        self.vmat.analyze(**self.analyze_kwargs)
