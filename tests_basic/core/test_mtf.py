import unittest

import numpy as np
from matplotlib import pyplot as plt
from parameterized import parameterized

from pylinac.core.mtf import MTF, EdgeSpreadFunctionMTF


class TestMTF(unittest.TestCase):
    def test_normal_mtf(self):
        pair_units = (0.1, 0.2, 0.3)
        maxs = (500, 300, 100)
        mins = (25, 50, 75)

        m = MTF(pair_units, maxs, mins)
        rm = m.relative_resolution(x=50)
        self.assertAlmostEqual(rm, 0.24, delta=0.03)
        rm = m.relative_resolution(x=90)
        self.assertAlmostEqual(rm, 0.15, delta=0.03)

    def test_mtf_lower_than_values(self):
        # should generate a warning
        pair_units = (0.1, 0.2, 0.3)
        maxs = (500, 300, 100)
        mins = (25, 50, 75)

        m = MTF(pair_units, maxs, mins)
        rm = m.relative_resolution(x=10)
        self.assertAlmostEqual(rm, 0.3, delta=0.03)

    def test_non_decreasing_mtf(self):
        # this will return the first occurrence where the condition is met.
        # should generate a warning
        pair_units = (0.1, 0.2, 0.3, 0.4)
        maxs = (500, 300, 500, 100)
        mins = (25, 50, 25, 75)

        MTF(pair_units, maxs, mins)

    def test_no_zero_division_of_line_pair_distances(self):
        old_settings = np.geterr()
        # set overflow to cause errors
        # shouldn't raise
        np.seterr(all="raise")
        pair_units = (0.1, 0.2, 0.3)
        maxs = (500, 300, 100)
        mins = (25, 50, 75)

        m = MTF(pair_units, maxs, mins)
        # shouldn't raise
        m.plot()
        # reset to old settings
        np.seterr(**old_settings)


class TestEdgeSpreadFunctionMTF(unittest.TestCase):
    def _test_mtf_of_step_function(self, mtf: EdgeSpreadFunctionMTF, sample_spacing=1):
        # In the ideal case: ESF ~ ..0001111.., LSF ~ ..0011000.., MTF = cos(pi * f)

        # Frequency
        n_profile_points = 2 * len(mtf.mtf)
        unscaled_sampling_frequency = 1 / n_profile_points
        unscaled_nyquist_frequency = 0.5
        unscaled_frequency_array = np.arange(
            0, unscaled_nyquist_frequency, unscaled_sampling_frequency
        )
        freq_nominal = unscaled_frequency_array / sample_spacing
        freq_actual = mtf.freq
        self.assertTrue(np.allclose(freq_nominal, freq_actual))

        # MTF (note: MTF = cos(pi*unscaled_freq), ie it does not depend on the sampling_freq)
        mtf_nominal = np.cos(np.pi * freq_nominal * sample_spacing)
        mtf_actual = mtf.mtf
        self.assertTrue(np.allclose(mtf_nominal, mtf_actual))

        # Relative resolution
        targets = np.array([30, 50, 80])
        resolution_nominal = np.arccos(targets / 100) / np.pi / sample_spacing
        resolution_actual = [mtf.relative_resolution(target) for target in targets]
        self.assertTrue(np.allclose(resolution_nominal, resolution_actual))

    def test_single_esf(self):
        samples_size = [8]
        esf = [np.append(np.zeros(x // 2), np.ones(x // 2)) for x in samples_size]
        mtf = EdgeSpreadFunctionMTF(esf)
        self._test_mtf_of_step_function(mtf)

    def test_multiple_esf(self):
        samples_size = [8, 6]
        esf = [np.append(np.zeros(x // 2), np.ones(x // 2)) for x in samples_size]
        mtf = EdgeSpreadFunctionMTF(esf)
        self._test_mtf_of_step_function(mtf)

    def test_sample_spacing(self):
        samples_size = [8, 6]
        esf = [np.append(np.zeros(x // 2), np.ones(x // 2)) for x in samples_size]
        sample_spacing = 10
        mtf = EdgeSpreadFunctionMTF(esf, sample_spacing)
        self._test_mtf_of_step_function(mtf, sample_spacing)

    def test_windowing(self):
        from scipy.signal import windows

        samples_size = [8, 6]
        esf = [np.append(np.zeros(x // 2), np.ones(x // 2)) for x in samples_size]
        mtf = EdgeSpreadFunctionMTF(esf, windowing=windows.kaiser, beta=0.5)
        self._test_mtf_of_step_function(mtf)

    def test_uncentered_esf(self):
        # The windowing function is symmetric so if the edge is not centered, windowing will skew the LSF.
        # This bias will mostly affect the higher frequencies (e.g. when testing agreement for 20% MTF)

        from scipy.signal import windows

        sample_size = 256
        sample_shift_from_center = 100
        esf = np.zeros(sample_size)
        esf[sample_size // 2 + sample_shift_from_center :] = 1

        # Without windowing there is no skew introduced in the step function so we get the ideal MTF
        mtf = EdgeSpreadFunctionMTF([esf], windowing=None)
        self._test_mtf_of_step_function(mtf)

        # With default windowing there is skew introduced in LSF
        mtf = EdgeSpreadFunctionMTF([esf])
        with self.assertRaises(AssertionError):
            self._test_mtf_of_step_function(mtf)

        # Adjusting the windowing can help correct the situation
        mtf = EdgeSpreadFunctionMTF([esf], windowing=windows.tukey, alpha=0.2)
        self._test_mtf_of_step_function(mtf)

    def test_padding_mode_none(self):
        samples_size = [256, 256]
        esf = [np.append(np.zeros(x // 2), np.ones(x // 2)) for x in samples_size]
        mtf = EdgeSpreadFunctionMTF(esf, padding_mode="none")
        self.assertEqual(len(mtf.mtf), samples_size[0] // 2)

    def test_padding_mode_fixed(self):
        samples_size = [256, 257]
        esf = [np.append(np.zeros(x // 2), np.ones(x // 2)) for x in samples_size]
        padding_size = 1000
        mtf = EdgeSpreadFunctionMTF(esf, padding_mode="fixed", num_samples=padding_size)
        self.assertEqual(len(mtf.mtf), padding_size // 2)

    def test_padding_mode_auto_default(self):
        samples_size = [256, 257]
        esf = [np.append(np.zeros(x // 2), np.ones(x // 2)) for x in samples_size]
        num_samples_default = 1024
        mtf = EdgeSpreadFunctionMTF(esf, padding_mode="auto")
        self.assertEqual(len(mtf.mtf), num_samples_default // 2)

    def test_padding_mode_auto_next_power_of_two(self):
        samples_size = [256, 1026]
        esf = [np.append(np.zeros(x // 2), np.ones(x // 2)) for x in samples_size]
        next_power_of_two = 2048
        mtf = EdgeSpreadFunctionMTF(esf, padding_mode="auto")
        self.assertEqual(len(mtf.mtf), next_power_of_two // 2)

    def test_error_if_none_padding_and_esf_have_different_size(self):
        samples_size = [256, 258]
        esf = [np.append(np.zeros(x // 2), np.ones(x // 2)) for x in samples_size]
        with self.assertRaises(ValueError):
            EdgeSpreadFunctionMTF(esf, padding_mode="none")

    def test_error_if_fixed_padding_size_smaller_than_samples_size(self):
        samples_size = [256, 258]
        esf = [np.append(np.zeros(x // 2), np.ones(x // 2)) for x in samples_size]
        padding_size = 256
        with self.assertRaises(ValueError):
            EdgeSpreadFunctionMTF(esf, padding_mode="fixed", num_samples=padding_size)

    @parameterized.expand([(None, "Cycles / sample"), (10, "Line pairs / mm")])
    def test_plot(self, sample_spacing, x_label):
        samples_size = 255
        esf = np.zeros(samples_size)
        esf[samples_size // 2 :] = 1
        mtf = EdgeSpreadFunctionMTF([esf], sample_spacing=sample_spacing)

        fig, axis = plt.subplots()
        mtf.plot(axis)
        plt.show()

        self.assertEqual(axis.get_xlabel(), x_label)

    @unittest.skip("For debugging only")
    def test_plot_debug(self):
        samples_size = [100, 256]
        esf = [np.append(np.zeros(x // 2), np.ones(x // 2)) for x in samples_size]
        sample_spacing = 10

        mtf = EdgeSpreadFunctionMTF([esf[0]], sample_spacing=sample_spacing)
        mtf._plot_debug()

        mtf = EdgeSpreadFunctionMTF(esf, sample_spacing=sample_spacing)
        mtf._plot_debug()

        mtf = EdgeSpreadFunctionMTF(esf, sample_spacing=sample_spacing)
        mtf._plot_debug(plot_together=True)
