from __future__ import annotations

import io
import json
import os.path as osp
from dataclasses import dataclass
from unittest import TestCase

from matplotlib import pyplot as plt

from pylinac import Centering, Edge, Interpolation, Normalization, Protocol
from pylinac.core import image
from pylinac.core.exceptions import NotAnalyzed
from pylinac.core.io import retrieve_demo_file
from pylinac.field_analysis import FieldAnalysis as FieldAnalysisV1
from pylinac.field_profile_analysis import FieldProfileAnalysis
from pylinac.metrics.profile import (
    CAXToLeftEdgeMetric,
    CAXToRightEdgeMetric,
    FlatnessDifferenceMetric,
    PenumbraLeftMetric,
    PenumbraRightMetric,
    ProfileMetric,
    SymmetryPointDifferenceMetric,
)
from tests_basic.utils import (
    CloudFileMixin,
    get_file_from_cloud_test_repo,
    has_www_connection,
    save_file,
)

TEST_DIR = "flatness_symmetry"  # trigger


def create_instance(model=FieldProfileAnalysis):
    fs = model.from_demo_image()
    fs.analyze()
    return fs


@dataclass
class MetricTest:
    metric: ProfileMetric
    value: float | dict[str, float]
    delta: float = 0.3


class FieldProfileAnalysisBase(CloudFileMixin):
    dir_path = ["flatness_symmetry"]
    invert = False
    edge_detection_method = Edge.INFLECTION_DERIVATIVE
    centering = Centering.BEAM_CENTER
    normalization_method = Normalization.GEOMETRIC_CENTER
    print_results = False
    position: tuple[float, float] = (0.5, 0.5)
    x_width: float = 0.03
    y_width: float = 0.03
    x_field_size: float
    y_field_size: float
    metric_tests: list[MetricTest]

    @classmethod
    def metrics(cls) -> list[ProfileMetric]:
        return [test.metric for test in cls.metric_tests]

    @classmethod
    def setUpClass(cls):
        cls.fa = FieldProfileAnalysis(cls.get_filename())
        cls.fa.analyze(
            centering=cls.centering,
            invert=cls.invert,
            edge_type=cls.edge_detection_method,
            metrics=cls.metrics(),
            position=cls.position,
            x_width=cls.x_width,
            y_width=cls.y_width,
        )
        if cls.print_results:
            print(cls.fa.results())
        super().setUpClass()

    @classmethod
    def tearDownClass(cls) -> None:
        plt.close("all")
        del cls.fa
        super().tearDownClass()

    def test_x_field_size(self):
        self.assertAlmostEqual(
            self.fa.results_data().x_metrics["Field Width (mm)"],
            self.x_field_size,
            delta=1,
        )

    def test_y_field_size(self):
        self.assertAlmostEqual(
            self.fa.results_data().y_metrics["Field Width (mm)"],
            self.y_field_size,
            delta=1,
        )

    def test_metrics(self):
        data = self.fa.results_data()
        for metric in self.metric_tests:
            if isinstance(metric.value, (float, int)):
                self.assertAlmostEqual(
                    data.y_metrics[metric.metric.full_name],
                    metric.value,
                    delta=metric.delta,
                    msg=metric.metric.full_name,
                )
                self.assertAlmostEqual(
                    data.x_metrics[metric.metric.full_name],
                    metric.value,
                    delta=metric.delta,
                    msg=metric.metric.full_name,
                )
            else:
                for p, v in metric.value.items():
                    if "x" == p:
                        self.assertAlmostEqual(
                            data.x_metrics[metric.metric.full_name],
                            v,
                            delta=metric.delta,
                            msg=f"{metric.metric.full_name} [x]",
                        )
                    if "y" == p:
                        self.assertAlmostEqual(
                            data.y_metrics[metric.metric.full_name],
                            v,
                            delta=metric.delta,
                            msg=f"{metric.metric.full_name} [y]",
                        )


class FieldProfileAnalysisPlusV1Comparison(FieldProfileAnalysisBase):
    """This is a V1 analysis test base used for comparing v1 to v2"""

    fs: FieldAnalysisV1
    is_FFF = False
    penumbra = (20, 80)
    protocol = Protocol.VARIAN
    interpolation_method = Interpolation.LINEAR
    in_field_ratio = 0.8
    slope_exclusion_ratio = 0.2

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.fs = FieldAnalysisV1(cls.get_filename())
        cls.fs.analyze(
            protocol=cls.protocol,
            centering=cls.centering,
            vert_position=cls.position[1],
            horiz_position=cls.position[0],
            vert_width=cls.y_width,
            horiz_width=cls.x_width,
            in_field_ratio=cls.in_field_ratio,
            slope_exclusion_ratio=cls.slope_exclusion_ratio,
            invert=cls.invert,
            is_FFF=cls.is_FFF,
            penumbra=cls.penumbra,
            interpolation=cls.interpolation_method,
            edge_detection_method=cls.edge_detection_method,
        )
        if cls.print_results:
            print(cls.fs.results())

    def test_symmetry(self):
        r1 = self.fs.results_data()
        r2 = self.fa.results_data()
        # x
        dv1 = r1.protocol_results["symmetry_horizontal"]
        dv2 = r2.x_metrics["Point Difference Symmetry (%)"]
        self.assertAlmostEqual(dv2, dv1, delta=0.15)
        # y
        dv1 = r1.protocol_results["symmetry_vertical"]
        dv2 = r2.y_metrics["Point Difference Symmetry (%)"]
        self.assertAlmostEqual(dv2, dv1, delta=0.15)

    def test_flatness(self):
        r1 = self.fs.results_data()
        r2 = self.fa.results_data()
        # x
        dv1 = r1.protocol_results["flatness_horizontal"]
        dv2 = r2.x_metrics["Flatness (Difference) (%)"]
        self.assertAlmostEqual(dv2, dv1, delta=0.1)
        # y
        dv1 = r1.protocol_results["flatness_vertical"]
        dv2 = r2.y_metrics["Flatness (Difference) (%)"]
        self.assertAlmostEqual(dv2, dv1, delta=0.1)

    def test_field_sizes(self):
        r1 = self.fs.results_data()
        r2 = self.fa.results_data()
        # x
        dv1 = r1.field_size_horizontal_mm
        dv2 = r2.x_metrics["Field Width (mm)"]
        self.assertAlmostEqual(dv2, dv1, delta=1)
        # y
        dv1 = r1.field_size_vertical_mm
        dv2 = r2.y_metrics["Field Width (mm)"]
        self.assertAlmostEqual(dv2, dv1, delta=1)

    def test_cax_to_edge(self):
        r1 = self.fs.results_data()
        r2 = self.fa.results_data()
        # left
        dv1 = r1.cax_to_left_mm
        dv2 = r2.x_metrics["CAX to Left Beam Edge (mm)"]
        self.assertAlmostEqual(dv2, dv1, delta=0.1)
        # right
        dv1 = r1.cax_to_right_mm
        dv2 = r2.x_metrics["CAX to Right Beam Edge (mm)"]
        self.assertAlmostEqual(dv2, dv1, delta=0.1)
        # top
        dv1 = r1.cax_to_top_mm
        dv2 = r2.y_metrics["CAX to Left Beam Edge (mm)"]
        self.assertAlmostEqual(dv2, dv1, delta=0.1)
        # bottom
        dv1 = r1.cax_to_bottom_mm
        dv2 = r2.y_metrics["CAX to Right Beam Edge (mm)"]
        self.assertAlmostEqual(dv2, dv1, delta=0.1)

    def test_penumbra(self):
        r1 = self.fs.results_data()
        r2 = self.fa.results_data()
        # left
        dv1 = r1.left_penumbra_mm
        dv2 = r2.x_metrics["Left Penumbra (mm)"]
        self.assertAlmostEqual(dv2, dv1, delta=0.5)
        # right
        dv1 = r1.right_penumbra_mm
        dv2 = r2.x_metrics["Right Penumbra (mm)"]
        self.assertAlmostEqual(dv2, dv1, delta=0.5)
        # top
        dv1 = r1.top_penumbra_mm
        dv2 = r2.y_metrics["Left Penumbra (mm)"]
        self.assertAlmostEqual(dv2, dv1, delta=0.5)
        # bottom
        dv1 = r1.bottom_penumbra_mm
        dv2 = r2.y_metrics["Right Penumbra (mm)"]
        self.assertAlmostEqual(dv2, dv1, delta=0.5)


class FieldAnalysisTests(TestCase):
    def tearDown(self) -> None:
        plt.close("all")

    def test_load_from_file_object(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, "6x-auto-bulb-2.dcm"])
        ref_fa = FieldProfileAnalysis(path)
        ref_fa.analyze()
        with open(path, "rb") as f:
            fa = FieldProfileAnalysis(f)
            fa.analyze()
        self.assertIsInstance(fa, FieldProfileAnalysis)
        self.assertEqual(fa.image.shape, ref_fa.image.shape)

    def test_load_from_stream(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, "6x-auto-bulb-2.dcm"])
        ref_fa = FieldProfileAnalysis(path)
        ref_fa.analyze()
        with open(path, "rb") as f:
            s = io.BytesIO(f.read())
            fa = FieldProfileAnalysis(s)
            fa.analyze()
        self.assertIsInstance(fa, FieldProfileAnalysis)
        self.assertEqual(fa.image.shape, ref_fa.image.shape)

    def test_demo_is_reachable(self):
        if has_www_connection():
            file = retrieve_demo_file(name="flatsym_demo.dcm", force=True)
            self.assertTrue(osp.isfile(file))

    def test_demo_loads_properly(self):
        """Loading the demo shouldn't raise an error"""
        FieldProfileAnalysis.from_demo_image()  # shouldn't raise

    def test_profile_limits(self):
        """Extreme profile limits should not raise an error"""
        fs = FieldProfileAnalysis.from_demo_image()
        fs.analyze()
        fs.analyze()
        fs.analyze()

    def test_analyze_sets_analyzed_flag(self):
        fs = create_instance()
        self.assertTrue(fs._is_analyzed)

    def test_results(self):
        fs = create_instance()
        self.assertIsInstance(fs.results(), str)

    def test_results_data(self):
        fs = create_instance()
        fs.analyze()
        # shouldn't raise
        fs.results_data()

        data_dict = fs.results_data(as_dict=True)
        self.assertIsInstance(data_dict, dict)

        data_json = fs.results_data(as_json=True)
        self.assertIsInstance(data_json, str)
        # shouldn't raise
        json.loads(data_json)

    def test_results_warnings(self):
        fs = create_instance()
        fs.analyze()
        # shouldn't raise
        data = fs.results_data()
        self.assertEqual(len(data.warnings), 0)

    def test_results_fails_if_not_analyzed(self):
        fs = FieldProfileAnalysis.from_demo_image()
        with self.assertRaises(NotAnalyzed):
            fs.results()

    def test_plot_works(self):
        fs = create_instance()
        fs.plot_analyzed_images()
        fs.plot_analyzed_images(show=True)

    def test_plot_fails_if_not_analysed(self):
        fs = FieldProfileAnalysis.from_demo_image()
        with self.assertRaises(NotAnalyzed):
            fs.plot_analyzed_images()

    def test_set_figure_size(self):
        fs = create_instance()
        figs = fs.plot_analyzed_images(figsize=(7, 11))
        self.assertEqual(figs[0].bbox_inches.height, 11)
        self.assertEqual(figs[0].bbox_inches.width, 7)

    def test_pdf_gets_generated(self):
        fs = create_instance()
        save_file(fs.publish_pdf)

    def test_pdf_fails_if_not_analyzed(self):
        fs = FieldProfileAnalysis.from_demo_image()
        with self.assertRaises(NotAnalyzed):
            fs.publish_pdf("dummy.pdf")

    def test_string_type_works_for_centering_interpolation_normalization_edge(self):
        fa = FieldProfileAnalysis.from_demo_image()
        fa.analyze(
            centering="Beam center",
            normalization="Beam center",
            edge_type="FWHM",
        )
        fa2 = FieldProfileAnalysis.from_demo_image()
        fa2.analyze(
            centering=Centering.BEAM_CENTER,
            normalization=Normalization.BEAM_CENTER,
            edge_type=Edge.FWHM,
        )
        self.assertEqual(
            fa.results_data().centering,
            fa2.results_data().centering,
        )
        self.assertEqual(
            fa.results_data().y_metrics["Field Width (mm)"],
            fa2.results_data().y_metrics["Field Width (mm)"],
        )

    def test_invalid_string_of_enum_fails(self):
        fa = FieldProfileAnalysis.from_demo_image()
        with self.assertRaises(ValueError):
            fa.analyze(normalization="limmerick")

    def test_image_kwargs(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, "6x-auto-bulb-2.dcm"])

        ref_fa = FieldProfileAnalysis(path)
        ref_fa.analyze()
        vert_ref = ref_fa.results_data().y_metrics["Field Width (mm)"]

        # pass kwarg; use same dpi as image; CAX offset should be the same, would be different with different DPI
        img = image.load(path)
        fa = FieldProfileAnalysis(path, dpi=img.dpi)
        fa.analyze()
        vert_manual = fa.results_data().y_metrics["Field Width (mm)"]

        self.assertEqual(vert_ref, vert_manual)


class FlatSymDemo(FieldProfileAnalysisPlusV1Comparison, TestCase):
    x_field_size = 141
    y_field_size = 200.2
    centering = Centering.BEAM_CENTER
    width = 0.03
    metric_tests = [
        MetricTest(
            metric=PenumbraRightMetric(),
            value={"x": 3.0, "y": 2.6},
        ),
        MetricTest(
            metric=PenumbraLeftMetric(),
            value={"x": 2.8, "y": 3.2},
        ),
        MetricTest(
            metric=FlatnessDifferenceMetric(),
            value=1.7,
        ),
        MetricTest(
            metric=SymmetryPointDifferenceMetric(),
            value={"y": -2.65, "x": -2.99},
        ),
        MetricTest(
            metric=CAXToLeftEdgeMetric(),
            value={"x": 60.4, "y": 99.8},
        ),
        MetricTest(
            metric=CAXToRightEdgeMetric(),
            value={"x": 80.5, "y": 100.5},
        ),
    ]

    @classmethod
    def get_filename(cls):
        return retrieve_demo_file(name="flatsym_demo.dcm")


class FlatSymWideDemo(FlatSymDemo, TestCase):
    """Use a wider slit to generate the profile. Results are the same"""

    width = 0.1


class NormalOpenField(FieldProfileAnalysisPlusV1Comparison, TestCase):
    """Typical field w/ horns"""

    file_name = "flat_open_15x15.dcm"
    edge_detection_method = Edge.FWHM
    x_field_size = 150
    y_field_size = 150
    metric_tests = [
        MetricTest(
            metric=PenumbraRightMetric(),
            value=3.3,
        ),
        MetricTest(
            metric=PenumbraLeftMetric(),
            value=3.3,
        ),
        MetricTest(
            metric=CAXToLeftEdgeMetric(),
            value=75,
        ),
        MetricTest(
            metric=CAXToRightEdgeMetric(),
            value=75,
        ),
        MetricTest(metric=FlatnessDifferenceMetric(), value=1.25, delta=0.1),
        MetricTest(metric=SymmetryPointDifferenceMetric(), value=0, delta=0.1),
    ]


class PerfectOpenField(FieldProfileAnalysisPlusV1Comparison, TestCase):
    """Completely flat field"""

    file_name = "perfect_open_15x15.dcm"
    edge_detection_method = Edge.FWHM
    x_field_size = 150
    y_field_size = 150
    metric_tests = [
        MetricTest(
            metric=PenumbraRightMetric(),
            value=3.3,
        ),
        MetricTest(
            metric=PenumbraLeftMetric(),
            value=3.3,
        ),
        MetricTest(
            metric=CAXToLeftEdgeMetric(),
            value=75,
        ),
        MetricTest(
            metric=CAXToRightEdgeMetric(),
            value=75,
        ),
        MetricTest(
            metric=FlatnessDifferenceMetric(),
            value=0,
        ),
        MetricTest(
            metric=SymmetryPointDifferenceMetric(),
            value=0,
        ),
    ]


class FFFOpenField(FieldProfileAnalysisPlusV1Comparison, TestCase):
    """FFF field. Note the same field size and penumbra as a flat beam"""

    file_name = "fff_open_15x15.dcm"
    edge_detection_method = Edge.INFLECTION_DERIVATIVE
    x_field_size = 150
    y_field_size = 150
    metric_tests = [
        MetricTest(
            metric=PenumbraRightMetric(),
            value=3.4,
        ),
        MetricTest(
            metric=PenumbraLeftMetric(),
            value=3.4,
        ),
        MetricTest(
            metric=CAXToLeftEdgeMetric(),
            value=75,
        ),
        MetricTest(
            metric=CAXToRightEdgeMetric(),
            value=75,
        ),
        MetricTest(metric=FlatnessDifferenceMetric(), value=5.2, delta=0.05),
        MetricTest(metric=SymmetryPointDifferenceMetric(), value=-0.11, delta=0.05),
    ]


class FFFOpenFieldHill(FFFOpenField, TestCase):
    """FFF field using Hill inflection. All metrics are the same except the y symmetry.
    This is because of an extra pixel (due to rounding) being included. The worst symmetry
    is at the extreme edge of the field."""

    edge_detection_method = Edge.INFLECTION_HILL
    metric_tests = [
        MetricTest(
            metric=PenumbraRightMetric(),
            value=3.4,
        ),
        MetricTest(
            metric=PenumbraLeftMetric(),
            value=3.4,
        ),
        MetricTest(
            metric=CAXToLeftEdgeMetric(),
            value=75,
        ),
        MetricTest(
            metric=CAXToRightEdgeMetric(),
            value=75,
        ),
        MetricTest(metric=FlatnessDifferenceMetric(), value=5.2, delta=0.05),
        MetricTest(
            metric=SymmetryPointDifferenceMetric(),
            value={"x": -0.11, "y": -0.22},  # only difference
            delta=0.05,
        ),
    ]


class FlatSym6X(FieldProfileAnalysisPlusV1Comparison, TestCase):
    file_name = "6x-auto-bulb-2.dcm"
    x_field_size = 99.4
    y_field_size = 99.4
    metric_tests = [
        MetricTest(
            metric=PenumbraRightMetric(),
            value={"x": 2.5, "y": 3.0},
        ),
        MetricTest(
            metric=PenumbraLeftMetric(),
            value={"x": 2.8, "y": 3.0},
        ),
        MetricTest(
            metric=FlatnessDifferenceMetric(),
            value=1.5,
        ),
        MetricTest(
            metric=SymmetryPointDifferenceMetric(),
            value={"x": -0.5, "y": -0.22},
        ),
        MetricTest(
            metric=CAXToLeftEdgeMetric(),
            value={"x": 49.2, "y": 49.9},
        ),
        MetricTest(
            metric=CAXToRightEdgeMetric(),
            value={"x": 50.3, "y": 49.5},
        ),
    ]


class FlatSym18X(FieldProfileAnalysisPlusV1Comparison, TestCase):
    file_name = "18x-auto-bulb2.dcm"
    x_field_size = 99.4
    y_field_size = 99.4
    metric_tests = [
        MetricTest(
            metric=PenumbraRightMetric(),
            value={"x": 3.0, "y": 3.4},
        ),
        MetricTest(
            metric=PenumbraLeftMetric(),
            value={"x": 3.4, "y": 3.7},
        ),
        MetricTest(
            metric=FlatnessDifferenceMetric(),
            value=1.5,
        ),
        MetricTest(
            metric=SymmetryPointDifferenceMetric(),
            value={"x": -0.5, "y": 0.4},
        ),
        MetricTest(
            metric=CAXToLeftEdgeMetric(),
            value={"x": 49.2, "y": 49.9},
        ),
        MetricTest(
            metric=CAXToRightEdgeMetric(),
            value={"x": 50.3, "y": 49.5},
        ),
    ]
