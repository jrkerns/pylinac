import json
from pathlib import Path
from unittest import TestCase

import numpy as np
from matplotlib import pyplot as plt

from pylinac.nuclear import (
    CenterOfRotation,
    CenterOfRotationResults,
    FourBarResolution,
    MaxCountRate,
    Nuclide,
    PlanarUniformity,
    PlanarUniformityResults,
    QuadrantResolution,
    QuadrantResolutionResults,
    SimpleSensitivity,
    SimpleSensitivityResults,
    TomographicContrast,
    TomographicContrastResults,
    TomographicResolution,
    TomographicUniformity,
    TomographicUniformityResults,
    determine_binning,
    integral_uniformity,
)
from tests_basic.utils import get_file_from_cloud_test_repo

TEST_DIR = Path("nuclear")


class TestMaxCountRate(TestCase):
    def setUp(self) -> None:
        plt.close("all")
        p = get_file_from_cloud_test_repo([TEST_DIR, "MaxCountRate.dcm"])
        self.m = MaxCountRate(p)

    def test_max_count_rate(self):
        self.m.analyze()
        # known value from NMQC toolkit
        self.assertAlmostEqual(self.m.max_countrate, 358437, delta=1)

    def test_max_frame(self):
        self.m.analyze()
        self.assertEqual(self.m.max_frame, 15)

    def test_frame_duration_is_accounted_for(self):
        # the same frame in half the time should have twice the countrate
        self.m.analyze(frame_duration=0.5)
        self.assertAlmostEqual(self.m.max_countrate, 358437 * 2, delta=1)

    def test_plot(self):
        # should not raise
        self.m.analyze()
        self.m.plot()
        self.m.plot(show=False)


class TestPlanarUniformityUtils(TestCase):
    def test_binning(self):
        binning = determine_binning(2)
        self.assertEqual(binning, 4)

        binning = determine_binning(4.1)
        self.assertEqual(binning, 2)

        binning = determine_binning(4.9)
        self.assertEqual(binning, 1)

    def test_integral_uniformity(self):
        # double check of michelson
        u = integral_uniformity(np.asarray([1, 1]))
        self.assertEqual(u, 0)

        u = integral_uniformity(np.asarray([10, 20]))
        self.assertAlmostEqual(u, 10 / 30 * 100)


class TestPlanarUniformity(TestCase):
    def setUp(self) -> None:
        p = get_file_from_cloud_test_repo([TEST_DIR, "PlanarUniformity.dcm"])
        self.m = PlanarUniformity(p)

    def test_results(self):
        self.m.analyze()
        self.assertEqual(len(self.m.frame_results), 1)
        self.assertAlmostEqual(
            self.m.frame_results["1"]["ufov"].differential_uniformity, 1.728, delta=0.15
        )  # value from NMQC
        self.assertAlmostEqual(
            self.m.frame_results["1"]["ufov"].integral_uniformity, 2.38, delta=0.01
        )
        self.assertAlmostEqual(
            self.m.frame_results["1"]["cfov"].differential_uniformity, 1.728, delta=0.01
        )
        self.assertAlmostEqual(
            self.m.frame_results["1"]["cfov"].integral_uniformity, 2.38, delta=0.01
        )
        self.assertEqual(self.m.frame_results["1"]["binned_frame"].shape, (64, 64))

    def test_results_data(self):
        self.m.analyze()
        data = self.m.results_data()
        self.assertIsInstance(data, dict)
        self.assertIsInstance(data["Frame 1"], PlanarUniformityResults)
        d = self.m.results_data(as_dict=True)
        self.assertIsInstance(d, dict)
        self.assertIsInstance(d["Frame 1"], dict)

    def test_changing_window_size_affects_results(self):
        self.m.analyze(window_size=5)
        self.assertAlmostEqual(
            self.m.frame_results["1"]["ufov"].differential_uniformity, 1.728, delta=0.15
        )  # default
        self.m.analyze(window_size=2)
        self.assertAlmostEqual(
            self.m.frame_results["1"]["ufov"].differential_uniformity, 1.01, delta=0.01
        )


class TestTwoFramePlanarUniformity(TestCase):
    def setUp(self) -> None:
        plt.close("all")
        p = get_file_from_cloud_test_repo([TEST_DIR, "TwoFrameUniformity.dcm"])
        self.m = PlanarUniformity(p)

    def test_results(self):
        self.m.analyze()
        self.assertEqual(len(self.m.frame_results), 2)
        # the NMQC toolkit gives a strange circular ROI for this, so values are just from pylinac
        # and quite different from NMQC
        self.assertAlmostEqual(
            self.m.frame_results["1"]["ufov"].differential_uniformity, 4.23, delta=0.01
        )
        self.assertAlmostEqual(
            self.m.frame_results["1"]["ufov"].integral_uniformity, 27.6, delta=0.01
        )
        self.assertAlmostEqual(
            self.m.frame_results["1"]["cfov"].differential_uniformity, 3.11, delta=0.01
        )
        self.assertAlmostEqual(
            self.m.frame_results["1"]["cfov"].integral_uniformity, 17.15, delta=0.01
        )
        self.assertEqual(self.m.frame_results["1"]["binned_frame"].shape, (128, 128))

        self.assertAlmostEqual(
            self.m.frame_results["2"]["ufov"].differential_uniformity, 4.05, delta=0.01
        )
        self.assertAlmostEqual(
            self.m.frame_results["2"]["ufov"].integral_uniformity, 28.16, delta=0.01
        )
        self.assertAlmostEqual(
            self.m.frame_results["2"]["cfov"].differential_uniformity, 3.33, delta=0.01
        )
        self.assertAlmostEqual(
            self.m.frame_results["2"]["cfov"].integral_uniformity, 17.67, delta=0.01
        )
        self.assertEqual(self.m.frame_results["2"]["binned_frame"].shape, (128, 128))


class TestCenterOfRotation102(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        plt.close("all")
        p = get_file_from_cloud_test_repo([TEST_DIR, "COR_102.dcm"])
        cls.m = CenterOfRotation(p)
        cls.m.analyze()

    def test_x_cor_deviation(self):
        self.assertAlmostEqual(self.m.x_cor_deviation_mm, 0.138, delta=0.01)

    def test_y_cor_deviation(self):
        self.assertAlmostEqual(self.m.y_cor_deviation_mm, 0.129, delta=0.01)

    def test_results_data(self):
        data = self.m.results_data()
        self.assertIsInstance(data, CenterOfRotationResults)
        self.assertAlmostEqual(data.x_deviation_mm, 0.138, delta=0.01)
        self.assertAlmostEqual(data.y_deviation_mm, 0.129, delta=0.01)

    def test_plot(self):
        # shouldn't raise
        figs, axes = self.m.plot()
        self.assertEqual(len(figs), 3)
        self.assertEqual(len(axes), 3)


class TestCenterOfRotation180(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        plt.close("all")
        p = get_file_from_cloud_test_repo([TEST_DIR, "COR_180.dcm"])
        cls.m = CenterOfRotation(p)
        cls.m.analyze()

    def test_x_cor_deviation(self):
        self.assertAlmostEqual(self.m.x_cor_deviation_mm, 0.170, delta=0.01)

    def test_y_cor_deviation(self):
        self.assertAlmostEqual(self.m.y_cor_deviation_mm, 0.207, delta=0.01)


class TestTomographicResolution(TestCase):
    def setUp(self) -> None:
        plt.close("all")
        p = get_file_from_cloud_test_repo([TEST_DIR, "TomoResolution.dcm"])
        self.m = TomographicResolution(p)
        self.m.analyze()

    def test_x_fwhm(self):
        self.assertAlmostEqual(self.m.x_axis.fwhm, 22.31, delta=0.01)

    def test_x_fwtm(self):
        self.assertAlmostEqual(self.m.x_axis.fwtm, 40.66, delta=0.01)

    def test_y_fwhm(self):
        self.assertAlmostEqual(self.m.y_axis.fwhm, 22.33, delta=0.01)

    def test_y_fwtm(self):
        self.assertAlmostEqual(self.m.y_axis.fwtm, 40.69, delta=0.01)

    def test_z_fwhm(self):
        self.assertAlmostEqual(self.m.z_axis.fwhm, 22.47, delta=0.01)

    def test_z_fwtm(self):
        self.assertAlmostEqual(self.m.z_axis.fwtm, 40.96, delta=0.01)


class TestSimpleSensitivityNoBackground(TestCase):
    def setUp(self) -> None:
        plt.close("all")
        p = get_file_from_cloud_test_repo([TEST_DIR, "Sensitivity_Phantom.dcm"])
        self.m = SimpleSensitivity(p)
        self.m.analyze(activity_mbq=10, nuclide=Nuclide.Tc99m)

    def test_duration(self):
        self.assertAlmostEqual(self.m.duration_s, 300, delta=1)

    def test_sensitivity_mbq(self):
        self.assertAlmostEqual(
            self.m.sensitivity_mbq, 245.83, delta=0.05
        )  # value from NMQC

    def test_sensitivity_uci(self):
        self.assertAlmostEqual(
            self.m.sensitivity_uci, 545.817, delta=0.05
        )  # value from NMQC

    def test_results_data(self):
        data = self.m.results_data()
        self.assertIsInstance(data, SimpleSensitivityResults)
        self.assertAlmostEqual(data.duration_s, 300, delta=0.01)
        self.assertAlmostEqual(data.sensitivity_mbq, 245.83, delta=0.05)

    def test_results(self):
        r = self.m.results()
        self.assertIsInstance(r, str)
        self.assertIn("Duration", r)


class TestSimpleSensitivityWithBackground(TestCase):
    @classmethod
    def setUp(cls) -> None:
        plt.close("all")
        p = get_file_from_cloud_test_repo([TEST_DIR, "Sensitivity_Phantom.dcm"])
        b = get_file_from_cloud_test_repo([TEST_DIR, "Sensitivity_Background.dcm"])
        cls.m = SimpleSensitivity(p, background_path=b)
        cls.m.analyze(activity_mbq=10, nuclide=Nuclide.Tc99m)

    def test_duration(self):
        self.assertAlmostEqual(self.m.duration_s, 300, delta=1)

    def test_sensitivity_mbq(self):
        self.assertAlmostEqual(self.m.sensitivity_mbq, 243.912, delta=0.2)

    def test_sensitivity_uci(self):
        self.assertAlmostEqual(self.m.sensitivity_uci, 541.484, delta=0.5)


class TestFourBar(TestCase):
    def setUp(self) -> None:
        p = get_file_from_cloud_test_repo([TEST_DIR, "FourBar.dcm"])
        self.m = FourBarResolution(p)
        self.m.analyze(separation_mm=100, roi_width_mm=10)

    # all reference values from NMQC

    def test_x_fwhm(self):
        self.assertAlmostEqual(self.m.x_axis.fwhm, 7.567, delta=0.1)

    def test_y_fwhm(self):
        self.assertAlmostEqual(self.m.y_axis.fwhm, 7.53, delta=0.1)

    def test_x_measured_pixel_size(self):
        self.assertAlmostEqual(self.m.x_axis.measured_pixel_size, 1.040, delta=0.01)

    def test_y_measured_pixel_size(self):
        self.assertAlmostEqual(self.m.y_axis.measured_pixel_size, 1.043, delta=0.01)

    def test_x_measured_difference(self):
        self.assertAlmostEqual(self.m.x_axis.pixel_size_difference, 1.065, delta=0.02)

    def test_y_measured_difference(self):
        self.assertAlmostEqual(self.m.y_axis.pixel_size_difference, 0.6899, delta=0.03)


class TestQuadrantResolution(TestCase):
    def setUp(self) -> None:
        p = get_file_from_cloud_test_repo([TEST_DIR, "QuadrantBar.dcm"])
        self.m = QuadrantResolution(p)
        self.m.analyze(bar_widths=(4.23, 3.18, 2.54, 2.12))

    def test_wrong_bar_lengths(self):
        with self.assertRaises(ValueError):
            self.m.analyze(bar_widths=(4.23, 3.18, 2.54))

    def test_mtfs(self):
        known_mtfs = (0.64, 0.462, 0.30, 0.202)  # from NMQC
        for known_mtf, (lpmm, mtf) in zip(known_mtfs, self.m.mtf.mtfs.items()):
            self.assertAlmostEqual(mtf, known_mtf, delta=0.01)

    def test_fwhms(self):
        known_fwhms = (2.995, 2.962, 2.954, 2.843)
        for known_fwhm, (lpmm, fwhm) in zip(known_fwhms, self.m.mtf.fwhms.items()):
            self.assertAlmostEqual(fwhm, known_fwhm, delta=0.05)

    def test_results_data(self):
        data = self.m.results_data(as_dict=True)
        self.assertIsInstance(data, dict)
        data = self.m.results_data()
        self.assertIsInstance(data, QuadrantResolutionResults)
        self.assertEqual(len(data.quadrants), 4)
        self.assertEqual(data.quadrants["1"]["spacing"], 4.23)
        self.assertAlmostEqual(data.quadrants["1"]["lpmm"], 0.118, delta=0.001)

    def test_plot(self):
        # shouldn't raise
        self.m.plot()
        self.m.plot(show=False)


class TestTomographicUniformity(TestCase):
    @classmethod
    def setUp(cls) -> None:
        p = get_file_from_cloud_test_repo([TEST_DIR, "Jaszczak.dcm"])
        cls.m = TomographicUniformity(p)

    def test_results(self):
        self.m.analyze()
        results = self.m.results()
        self.assertIsInstance(results, str)

    def test_plot(self):
        self.m.analyze()
        # shouldn't raise
        self.m.plot()
        self.m.plot(show=False)

    def test_results_data(self):
        self.m.analyze()
        data = self.m.results_data(as_dict=True)
        self.assertIsInstance(data, dict)
        data = self.m.results_data()
        self.assertIsInstance(data, TomographicUniformityResults)

    def test_first_frame_lower_than_0_raises(self):
        with self.assertRaises(ValueError):
            self.m.analyze(first_frame=-3)

    def test_last_frame_rolls_over_to_positive(self):
        # 45 frames + 1-index based - 3 = 43
        self.m.analyze(last_frame=-3)
        self.assertEqual(self.m.last_frame, 43)

    def test_last_frame_outside_bounds(self):
        with self.assertRaises(ValueError):
            self.m.analyze(last_frame=100)

    def test_last_frame_less_than_first_frame(self):
        with self.assertRaises(ValueError):
            self.m.analyze(first_frame=10, last_frame=5)

    def test_center_ratio(self):
        self.m.analyze()
        self.assertAlmostEqual(
            self.m.center_ratio, 1.004, delta=0.005
        )  # value from NMQC

    def test_cfov_integral_uniformity(self):
        self.m.analyze()
        r = self.m.results_data()
        self.assertAlmostEqual(r.cfov_integral_uniformity, 13.031, delta=0.3)

    def test_cfov_differential_uniformity(self):
        self.m.analyze()
        r = self.m.results_data()
        self.assertAlmostEqual(r.cfov_differential_uniformity, 9.605, delta=0.2)


class TestTomographicContrast(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        p = get_file_from_cloud_test_repo([TEST_DIR, "Jaszczak.dcm"])
        cls.m = TomographicContrast(p)
        cls.m.analyze()

    def test_results(self):
        results = self.m.results()
        self.assertIsInstance(results, str)
        self.assertIn("Contrast", results)

    def test_results_data(self):
        data = self.m.results_data()
        self.assertIsInstance(data, TomographicContrastResults)
        self.assertAlmostEqual(
            data.uniformity_baseline, 1331, delta=12
        )  # value from NMQC
        self.assertAlmostEqual(data.spheres["1"]["max_contrast"], 72.1, delta=0.2)
        self.assertAlmostEqual(data.spheres["1"]["mean_contrast"], 34.2, delta=1.7)
        self.assertAlmostEqual(data.spheres["6"]["max_contrast"], 18.3, delta=0.4)
        self.assertAlmostEqual(
            data.spheres["6"]["mean_contrast"], 13.65, delta=4
        )  # the notable exception in difference

        data_dict = self.m.results_data(as_dict=True)
        self.assertIsInstance(data_dict, dict)

        data_str = self.m.results_data(as_json=True)
        self.assertIsInstance(data_str, str)
        # shouldn't raise
        json.loads(data_str)

    def plot(self):
        # shouldn't raise
        self.m.plot()
        self.m.plot(show=False)

    def test_different_diameter_vs_angle_lengths_error(self):
        with self.assertRaises(ValueError):
            self.m.analyze(sphere_diameters_mm=(10, 13), sphere_angles=(0, 90, 180))
