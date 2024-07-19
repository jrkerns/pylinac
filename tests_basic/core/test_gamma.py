from unittest import TestCase, skip

import numpy as np

from pylinac.core.gamma import gamma_1d, gamma_2d
from pylinac.core.image import DicomImage
from tests_basic.utils import get_file_from_cloud_test_repo


class TestAgnewMcGarry(TestCase):
    """Tests from the Agnew & McGarry paper. https://www.sciencedirect.com/science/article/abs/pii/S0167814015006660"""

    GEOM_1mm = ("Gamma", "Agnew", "Geometric", "1mm")
    GEOM_025mm = ("Gamma", "Agnew", "Geometric", "0.25mm")

    def test_dta_11_1mm(self):
        PASS_RATE = 92.3  # from Agnew's paper
        ref_file = get_file_from_cloud_test_repo(
            (*self.GEOM_1mm, "GeometricSquare_Reference_1mmPx.dcm")
        )
        eval_file = get_file_from_cloud_test_repo(
            (*self.GEOM_1mm, "GeometricSquare_Evaluated_1mmPx_DTA_Test.dcm")
        )
        ref_array = DicomImage(ref_file, raw_pixels=True).array
        eval_array = DicomImage(eval_file, raw_pixels=True).array
        g = gamma_2d(
            reference=ref_array,
            evaluation=eval_array,
            dose_to_agreement=1,
            distance_to_agreement=1,
            global_dose=True,
        )
        meas_pass_rate = np.sum(g < 1) / g.size * 100  # 92%
        self.assertAlmostEqual(meas_pass_rate, PASS_RATE, delta=0.31)

    def test_dd_11_1mm(self):
        PASS_RATE = 96  # from Agnew's paper; note that in the supplementary material the math is wrong 60000/62500 = 96%, not 92%.
        ref_file = get_file_from_cloud_test_repo(
            (*self.GEOM_1mm, "GeometricSquare_Reference_1mmPx.dcm")
        )
        eval_file = get_file_from_cloud_test_repo(
            (*self.GEOM_1mm, "GeometricSquare_Evaluated_1mmPx_DD_Test.dcm")
        )
        ref_array = DicomImage(ref_file, raw_pixels=True).array
        eval_array = DicomImage(eval_file, raw_pixels=True).array
        g = gamma_2d(
            reference=ref_array,
            evaluation=eval_array,
            dose_to_agreement=1,
            distance_to_agreement=1,
            global_dose=True,
        )
        meas_pass_rate = np.sum(g < 1) / g.size * 100
        self.assertAlmostEqual(meas_pass_rate, PASS_RATE, delta=0.05)

    @skip("Very slow; performed locally")
    def test_dta_11_025mm(self):
        PASS_RATE = 92.3  # from Agnew's paper
        ref_file = get_file_from_cloud_test_repo(
            (*self.GEOM_025mm, "GeometricSquare_Reference_0_25mmPx.dcm")
        )
        eval_file = get_file_from_cloud_test_repo(
            (*self.GEOM_025mm, "GeometricSquare_Evaluated_0_25mmPx_DTA_Test.dcm")
        )
        ref_array = DicomImage(ref_file, raw_pixels=True).array
        eval_array = DicomImage(eval_file, raw_pixels=True).array
        g = gamma_2d(
            reference=ref_array,
            evaluation=eval_array,
            dose_to_agreement=1,
            distance_to_agreement=4,  # 0.25mm pixel pitch; 4*0.25 = 1mm
            global_dose=True,
        )
        meas_pass_rate = np.sum(g < 1) / g.size * 100  # 92.19%
        self.assertAlmostEqual(meas_pass_rate, PASS_RATE, delta=0.15)

    @skip("Very slow; performed locally")
    def test_dd_11_025mm(self):
        PASS_RATE = 96  # from Agnew's paper
        ref_file = get_file_from_cloud_test_repo(
            (*self.GEOM_025mm, "GeometricSquare_Reference_0_25mmPx.dcm")
        )
        eval_file = get_file_from_cloud_test_repo(
            (*self.GEOM_025mm, "GeometricSquare_Evaluated_0_25mmPx_DD_Test.dcm")
        )
        ref_array = DicomImage(ref_file, raw_pixels=True).array
        eval_array = DicomImage(eval_file, raw_pixels=True).array
        g = gamma_2d(
            reference=ref_array,
            evaluation=eval_array,
            dose_to_agreement=1,
            distance_to_agreement=4,  # 0.25mm pixel pitch; 4*0.25 = 1mm
            global_dose=True,
        )
        meas_pass_rate = np.sum(g < 1) / g.size * 100  # 95.98%
        self.assertAlmostEqual(meas_pass_rate, PASS_RATE, delta=0.05)


class TestGamma2D(TestCase):
    def test_perfect_match_is_0(self):
        ref = eval = np.ones((5, 5))
        gamma = gamma_2d(reference=ref, evaluation=eval)
        self.assertEqual(gamma.max(), 0)
        self.assertEqual(gamma.min(), 0)
        self.assertEqual(gamma.size, 25)

        # test a high measurement value
        ref = eval = np.ones((5, 5)) * 50
        gamma = gamma_2d(reference=ref, evaluation=eval)
        self.assertEqual(gamma.max(), 0)
        self.assertEqual(gamma.min(), 0)
        self.assertEqual(gamma.size, 25)

    def test_gamma_perfectly_at_1(self):
        # offset a profile exactly by the dose to agreement
        ref = np.ones((5, 5))
        eval = np.ones((5, 5)) * 1.01
        gamma = gamma_2d(reference=ref, evaluation=eval, dose_to_agreement=1)
        self.assertAlmostEqual(gamma.max(), 1, delta=0.001)
        self.assertAlmostEqual(gamma.min(), 1, delta=0.001)

        # test same but eval is LOWER than ref
        ref = np.ones((5, 5))
        eval = np.ones((5, 5)) * 0.99
        gamma = gamma_2d(reference=ref, evaluation=eval, dose_to_agreement=1)
        self.assertAlmostEqual(gamma.max(), 1, delta=0.001)
        self.assertAlmostEqual(gamma.min(), 1, delta=0.001)

    def test_gamma_some_on_some_off(self):
        ref = np.ones((5, 5))
        eval = np.ones((5, 5))
        eval[(0, 0, 1, 1), (0, 1, 1, 0)] = 1.03  # set top left corner to 3% off
        gamma = gamma_2d(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=1,
            distance_to_agreement=1,
            gamma_cap_value=5,
        )
        self.assertAlmostEqual(gamma[0, 0], 3, delta=0.01)  # fully off by 3
        self.assertAlmostEqual(
            gamma[0, 1], 1, delta=0.01
        )  # dose at next pixel matches (dose=0, dist=1)
        self.assertAlmostEqual(gamma[-1, -1], 0, delta=0.01)  # gamma at end is perfect

        # check inverted pattern is mirrored (checks off-by-one errors)
        ref = np.ones((5, 5))
        eval = np.ones((5, 5))
        eval[
            (-1, -1, -2, -2), (-1, -2, -2, -1)
        ] = 1.03  # set bottom right corner to 3% off
        gamma = gamma_2d(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=1,
            distance_to_agreement=1,
            gamma_cap_value=5,
        )
        self.assertAlmostEqual(gamma[0, 0], 0, delta=0.01)
        self.assertAlmostEqual(gamma[-1, -2], 1, delta=0.01)
        self.assertAlmostEqual(gamma[-1, -1], 3, delta=0.01)

    def test_localized_dose(self):
        ref = np.ones((5, 5))
        ref[0, 0] = 100
        eval = np.ones((5, 5))
        eval[0, 0] = 103
        eval[0, 1] = 1.03
        # with global, element 2 is easily under gamma 1 since DTA there is 3
        gamma = gamma_2d(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=3,
            distance_to_agreement=1,
            gamma_cap_value=5,
            global_dose=False,
            dose_threshold=0,
        )
        self.assertAlmostEqual(gamma[0, 0], 1, delta=0.01)  # fully off by 3
        self.assertAlmostEqual(
            gamma[0, 1], 1, delta=0.01
        )  # dose here is also off by 3% relative dose
        self.assertAlmostEqual(gamma[-1, -1], 0, delta=0.01)  # gamma at end is perfect

    def test_threshold(self):
        ref = np.zeros((5, 5))
        ref[0, 0] = 1
        eval = ref
        # only one point should be computed as rest are under default threshold
        gamma = gamma_2d(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=3,
            distance_to_agreement=1,
            gamma_cap_value=5,
            global_dose=False,
            dose_threshold=5,
        )
        self.assertAlmostEqual(gamma[0, 0], 0, delta=0.01)
        self.assertTrue(np.isnan(gamma[0, 1]))
        self.assertTrue(np.isnan(gamma[-1, -1]))

    def test_fill_value(self):
        ref = np.zeros((5, 5))
        ref[0, 0] = 1
        eval = ref
        # only one point should be computed as rest are under default threshold
        gamma = gamma_2d(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=3,
            distance_to_agreement=1,
            gamma_cap_value=5,
            global_dose=False,
            dose_threshold=5,
            fill_value=0.666,
        )
        self.assertAlmostEqual(gamma[0, 0], 0, delta=0.01)
        self.assertAlmostEqual(
            gamma[0, 1], 0.666, delta=0.01
        )  # dose here is also off by 3% relative dose
        self.assertAlmostEqual(gamma[-1, -1], 0.666, delta=0.01)

    def test_gamma_half(self):
        # offset a profile by half the dose to agreement to ensure it's 0.5
        ref = np.ones((5, 5))
        eval = np.ones((5, 5)) / 1.005
        gamma = gamma_2d(reference=ref, evaluation=eval, dose_to_agreement=1)
        self.assertAlmostEqual(gamma.max(), 0.5, delta=0.01)
        self.assertAlmostEqual(gamma.min(), 0.5, delta=0.01)

    def test_gamma_cap(self):
        # cap to the value
        ref = np.ones((5, 5))
        eval = np.ones((5, 5)) * 10
        gamma = gamma_2d(
            reference=ref, evaluation=eval, dose_to_agreement=1, gamma_cap_value=2
        )
        self.assertEqual(gamma.max(), 2)
        self.assertEqual(gamma.min(), 2)

    def test_non_2d_array(self):
        ref = np.ones(5)
        eval = np.ones((5, 5))
        with self.assertRaises(ValueError):
            gamma_2d(reference=ref, evaluation=eval)

        ref = np.ones((5, 5))
        eval = np.ones(5)
        with self.assertRaises(ValueError):
            gamma_2d(reference=ref, evaluation=eval)


class TestGamma1D(TestCase):
    def test_same_profile_is_0_gamma(self):
        ref = eval = np.ones(5)
        gamma, _, _ = gamma_1d(reference=ref, evaluation=eval)
        self.assertEqual(max(gamma), 0)
        self.assertEqual(min(gamma), 0)
        self.assertEqual(len(gamma), 5)

        # test a high measurement value
        ref = eval = np.ones(5) * 50
        gamma, _, _ = gamma_1d(reference=ref, evaluation=eval)
        self.assertEqual(max(gamma), 0)
        self.assertEqual(min(gamma), 0)
        self.assertEqual(len(gamma), 5)

    def test_gamma_perfectly_at_1(self):
        # offset a profile exactly by the dose to agreement
        ref = np.ones(5)
        eval = np.ones(5) * 1.01
        gamma, _, _ = gamma_1d(reference=ref, evaluation=eval, dose_to_agreement=1)
        self.assertAlmostEqual(max(gamma), 1, delta=0.001)
        self.assertAlmostEqual(min(gamma), 1, delta=0.001)

        # test same but eval is LOWER than ref
        ref = np.ones(5)
        eval = np.ones(5) * 0.99
        gamma, _, _ = gamma_1d(reference=ref, evaluation=eval, dose_to_agreement=1)
        self.assertAlmostEqual(max(gamma), 1, delta=0.001)
        self.assertAlmostEqual(min(gamma), 1, delta=0.001)

    def test_gamma_half(self):
        # offset a profile by half the dose to agreement to ensure it's 0.5
        ref = np.ones(5)
        eval = np.ones(5) / 1.005
        gamma, _, _ = gamma_1d(reference=ref, evaluation=eval, dose_to_agreement=1)
        self.assertAlmostEqual(max(gamma), 0.5, delta=0.01)
        self.assertAlmostEqual(min(gamma), 0.5, delta=0.01)

    def test_gamma_some_on_some_off(self):
        ref = np.ones(5)
        eval = np.asarray((1.03, 1.03, 1, 1, 1))
        gamma, _, _ = gamma_1d(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=1,
            distance_to_agreement=1,
            gamma_cap_value=5,
        )
        self.assertAlmostEqual(gamma[0], 3, delta=0.01)  # fully off by 3
        self.assertAlmostEqual(gamma[1], 3, delta=0.01)
        self.assertAlmostEqual(gamma[-1], 0, delta=0.01)  # gamma at end is perfect

        # check inverted pattern is mirrored (checks off-by-one errors)
        ref = np.ones(5)
        eval = np.asarray((1, 1, 1, 1.03, 1.03))
        gamma, _, _ = gamma_1d(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=1,
            distance_to_agreement=1,
            gamma_cap_value=5,
        )
        self.assertAlmostEqual(gamma[0], 0, delta=0.01)
        self.assertAlmostEqual(gamma[-2], 3, delta=0.01)
        self.assertAlmostEqual(gamma[-1], 3, delta=0.01)

    def test_localized_dose(self):
        ref = eval = np.array((100, 1, 1, 1, 1))
        gamma, _, _ = gamma_1d(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=3,
            distance_to_agreement=1,
            gamma_cap_value=5,
            global_dose=False,
        )
        self.assertAlmostEqual(gamma[0], 0, delta=0.01)
        self.assertTrue(np.isnan(gamma[-2]))
        self.assertTrue(np.isnan(gamma[-1]))

    def test_threshold(self):
        ref = np.zeros(5)
        ref[0] = 1
        eval = ref
        # only one point should be computed as rest are under default threshold
        gamma, _, _ = gamma_1d(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=3,
            distance_to_agreement=1,
            gamma_cap_value=5,
            global_dose=False,
            dose_threshold_percent=5,
        )
        self.assertAlmostEqual(gamma[0], 0, delta=0.01)
        self.assertTrue(np.isnan(gamma[-2]))
        self.assertTrue(np.isnan(gamma[-1]))

    def test_fill_value(self):
        ref = np.zeros(5)
        ref[0] = 1
        eval = ref
        # only one point should be computed as rest are under default threshold
        gamma, _, _ = gamma_1d(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=3,
            distance_to_agreement=1,
            gamma_cap_value=5,
            global_dose=False,
            dose_threshold_percent=5,
            fill_value=0.666,
        )
        self.assertAlmostEqual(gamma[0], 0, delta=0.01)
        self.assertAlmostEqual(gamma[-2], 0.666, delta=0.01)
        self.assertAlmostEqual(gamma[-1], 0.666, delta=0.01)

    def test_gamma_cap(self):
        # cap to the value
        ref = np.ones(5)
        eval = np.ones(5) * 10
        gamma, _, _ = gamma_1d(
            reference=ref, evaluation=eval, dose_to_agreement=1, gamma_cap_value=2
        )
        self.assertEqual(max(gamma), 2)
        self.assertEqual(min(gamma), 2)

    def test_non_1d_array(self):
        ref = np.ones(5)
        eval = np.ones((5, 5))
        with self.assertRaises(ValueError):
            gamma_1d(reference=ref, evaluation=eval)
