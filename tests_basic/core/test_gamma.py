import math
from unittest import TestCase, skip

import numpy as np

from pylinac import Centering
from pylinac.core.gamma import _compute_distance, gamma_1d, gamma_2d, gamma_geometric
from pylinac.core.image import DicomImage
from pylinac.core.image_generator import AS1000Image, AS1200Image, PerfectFieldLayer
from pylinac.core.profile import FWXMProfile, FWXMProfilePhysical
from tests_basic.core.test_profile import generate_open_field
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


class TestGammaGeometric(TestCase):
    def test_point_projection_far_from_simplex(self):
        # test a point that is far from the simplex
        # goes to the closest point on the simplex (for 1D this is okay; not for 2D)
        point = np.array([2, 2])
        vertices = [np.array([0, 0]), np.array([1, 1])]
        distance = _compute_distance(point, vertices)
        self.assertAlmostEqual(distance, math.sqrt(2), delta=0.01)

    def test_point_along_simplex(self):
        # test a point that is along the simplex
        point = np.array([0.5, 0.5])
        vertices = [np.array([0, 0]), np.array([1, 1])]
        distance = _compute_distance(point, vertices)
        self.assertAlmostEqual(distance, 0, delta=0.01)

    def test_simple_distance(self):
        point = np.array([0, 1])
        vertices = [np.array([0, 0]), np.array([1, 1])]
        distance = _compute_distance(point, vertices)
        self.assertAlmostEqual(distance, math.sqrt(2) / 2, delta=0.01)

    def test_ref_x_values_not_same_length_as_ref(self):
        ref = eval = np.ones(5)
        with self.assertRaises(ValueError):
            gamma_geometric(
                reference=ref, evaluation=eval, reference_coordinates=np.arange(6)
            )

    def test_eval_x_values_not_same_length_as_eval(self):
        ref = eval = np.ones(5)
        with self.assertRaises(ValueError):
            gamma_geometric(
                reference=ref, evaluation=eval, evaluation_coordinates=np.arange(6)
            )

    def test_min_eval_x_lower_than_min_ref_x(self):
        ref = eval = np.ones(5)
        with self.assertRaises(ValueError):
            gamma_geometric(
                reference=ref,
                evaluation=eval,
                evaluation_coordinates=(np.arange(5) - 1),
                reference_coordinates=np.arange(5),
            )

    def test_max_eval_x_higher_than_max_ref_x(self):
        ref = eval = np.ones(5)
        with self.assertRaises(ValueError):
            gamma_geometric(
                reference=ref,
                evaluation=eval,
                evaluation_coordinates=(np.arange(5) + 1),
                reference_coordinates=np.arange(5),
            )

    def test_same_profile_is_0_gamma(self):
        ref = eval = np.ones(5)
        gamma = gamma_geometric(reference=ref, evaluation=eval)
        self.assertEqual(max(gamma), 0)
        self.assertEqual(min(gamma), 0)
        self.assertEqual(len(gamma), 5)

        # test a high measurement value
        ref = eval = np.ones(5) * 50
        gamma = gamma_geometric(reference=ref, evaluation=eval)
        self.assertEqual(max(gamma), 0)
        self.assertEqual(min(gamma), 0)
        self.assertEqual(len(gamma), 5)

    def test_gamma_perfectly_at_1(self):
        # offset a profile exactly by the dose to agreement
        ref = np.ones(5)
        eval = np.ones(5) * 1.01
        gamma = gamma_geometric(reference=ref, evaluation=eval, dose_to_agreement=1)
        self.assertAlmostEqual(max(gamma), 1, delta=0.001)
        self.assertAlmostEqual(min(gamma), 1, delta=0.001)

        # test same but eval is LOWER than ref
        ref = np.ones(5)
        eval = np.ones(5) * 0.99
        gamma = gamma_geometric(reference=ref, evaluation=eval, dose_to_agreement=1)
        self.assertAlmostEqual(max(gamma), 1, delta=0.001)
        self.assertAlmostEqual(min(gamma), 1, delta=0.001)

    def test_gamma_half(self):
        # offset a profile by half the dose to agreement to ensure it's 0.5
        ref = np.ones(5)
        eval = np.ones(5) / 1.005
        gamma = gamma_geometric(reference=ref, evaluation=eval, dose_to_agreement=1)
        self.assertAlmostEqual(max(gamma), 0.5, delta=0.01)
        self.assertAlmostEqual(min(gamma), 0.5, delta=0.01)

    def test_gamma_some_on_some_off(self):
        ref = np.ones(5)
        eval = np.asarray((1.03, 1.03, 1, 1, 1))
        gamma = gamma_geometric(
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
        gamma = gamma_geometric(
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
        gamma = gamma_geometric(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=3,
            distance_to_agreement=1,
            gamma_cap_value=5,
        )
        self.assertAlmostEqual(gamma[0], 0, delta=0.01)
        self.assertTrue(np.isnan(gamma[-2]))
        self.assertTrue(np.isnan(gamma[-1]))

    def test_threshold(self):
        ref = np.zeros(5)
        ref[0] = 1
        eval = ref
        # only one point should be computed as rest are under default threshold
        gamma = gamma_geometric(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=3,
            distance_to_agreement=1,
            gamma_cap_value=5,
            dose_threshold=5,
        )
        self.assertAlmostEqual(gamma[0], 0, delta=0.01)
        self.assertTrue(np.isnan(gamma[-2]))
        self.assertTrue(np.isnan(gamma[-1]))

    def test_fill_value(self):
        ref = np.zeros(5)
        ref[0] = 1
        eval = ref
        # only one point should be computed as rest are under default threshold
        gamma = gamma_geometric(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=3,
            distance_to_agreement=1,
            gamma_cap_value=5,
            dose_threshold=5,
            fill_value=0.666,
        )
        self.assertAlmostEqual(gamma[0], 0, delta=0.01)
        self.assertAlmostEqual(gamma[-2], 0.666, delta=0.01)
        self.assertAlmostEqual(gamma[-1], 0.666, delta=0.01)

    def test_gamma_cap(self):
        # cap to the value
        ref = np.ones(5)
        eval = np.ones(5) * 10
        gamma = gamma_geometric(
            reference=ref, evaluation=eval, dose_to_agreement=1, gamma_cap_value=2
        )
        self.assertEqual(max(gamma), 2)
        self.assertEqual(min(gamma), 2)

    def test_non_1d_array(self):
        ref = np.ones(5)
        eval = np.ones((5, 5))
        with self.assertRaises(ValueError):
            gamma_geometric(reference=ref, evaluation=eval)


class TestGamma1D(TestCase):
    def test_resolution_below_1mm(self):
        ref = eval = np.ones(5)
        with self.assertRaises(ValueError):
            gamma_1d(reference=ref, evaluation=eval, resolution_factor=0)

    def test_resolution_not_int(self):
        ref = eval = np.ones(5)
        with self.assertRaises(ValueError):
            gamma_1d(reference=ref, evaluation=eval, resolution_factor=1.5)

    def test_ref_x_values_not_same_length_as_ref(self):
        ref = eval = np.ones(5)
        with self.assertRaises(ValueError):
            gamma_1d(reference=ref, evaluation=eval, reference_coordinates=np.arange(6))

    def test_eval_x_values_not_same_length_as_eval(self):
        ref = eval = np.ones(5)
        with self.assertRaises(ValueError):
            gamma_1d(
                reference=ref, evaluation=eval, evaluation_coordinates=np.arange(6)
            )

    def test_min_eval_x_lower_than_min_ref_x(self):
        ref = eval = np.ones(5)
        with self.assertRaises(ValueError):
            gamma_1d(
                reference=ref,
                evaluation=eval,
                evaluation_coordinates=(np.arange(5) - 3),
                reference_coordinates=np.arange(5),
            )

    def test_max_eval_x_higher_than_max_ref_x(self):
        ref = eval = np.ones(5)
        with self.assertRaises(ValueError):
            gamma_1d(
                reference=ref,
                evaluation=eval,
                evaluation_coordinates=(np.arange(5) + 3),
                reference_coordinates=np.arange(5),
            )

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
            dose_threshold=5,
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
            dose_threshold=5,
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


class TestGammaFromProfile(TestCase):
    """Test the gamma method from the physical profile class"""

    def test_other_profile_not_physical(self):
        p = FWXMProfilePhysical(values=np.ones(5), dpmm=1)
        p2 = FWXMProfile(values=np.ones(5))
        with self.assertRaises(ValueError):
            p.gamma(reference_profile=p2)

    def test_gamma_on_self_is_zero(self):
        p = FWXMProfilePhysical(values=np.ones(5), dpmm=1)
        gamma = p.gamma(reference_profile=p)
        self.assertTrue(np.allclose(gamma, 0))

    def test_dose_gamma(self):
        vals = np.ones(5) * 1.03
        ref_vals = np.ones(5)
        p = FWXMProfilePhysical(values=vals, dpmm=1)
        ref_p = FWXMProfilePhysical(values=ref_vals, dpmm=1)
        gamma = p.gamma(reference_profile=ref_p, dose_to_agreement=3, gamma_cap_value=2)
        # all right at 1 because dose is 3% high
        self.assertTrue(np.allclose(gamma, 1))

    def test_low_density_eval(self):
        # the eval is 1/2 the spacing, but same values, as the reference
        vals = np.arange(1, 11)[::2]
        ref_vals = np.arange(11)
        p = FWXMProfilePhysical(values=vals, x_values=np.arange(1, 11, 2))
        ref_p = FWXMProfilePhysical(
            values=ref_vals,
        )
        gamma = p.gamma(reference_profile=ref_p, dose_to_agreement=1, gamma_cap_value=2)
        # gamma is perfect
        self.assertTrue(np.allclose(gamma, 0))

    def test_different_epids(self):
        """This test the same profile but with different EPIDs (i.e. pixel size)"""
        # we offset the reference by 1% to ensure we have a realistic gamma value
        img1200 = generate_open_field(
            field_size=(100, 100),
            imager=AS1200Image,
            alpha=0.99,
            field=PerfectFieldLayer,
        )
        img1000 = generate_open_field(
            field_size=(100, 100), imager=AS1000Image, field=PerfectFieldLayer
        )
        p1200 = img1200.image[640, :]
        p1000 = img1000.image[384, :]
        p1200_prof = FWXMProfilePhysical(values=p1200, dpmm=1 / img1200.pixel_size)
        p1000_prof = FWXMProfilePhysical(values=p1000, dpmm=1 / img1000.pixel_size)
        gamma = p1000_prof.gamma(
            reference_profile=p1200_prof, dose_to_agreement=1, gamma_cap_value=2
        )
        # gamma is very low; just pixel noise from the image generator
        self.assertAlmostEqual(np.nanmean(gamma), 0.948, delta=0.01)

    def test_beam_centering_offset_profile(self):
        # similar test as above, but profiles are offset and no dose difference
        # we expect gamma to be very low
        # note the 5mm offset below which would otherwise cause a large gamma range
        img1200 = generate_open_field(
            field_size=(100, 100),
            imager=AS1200Image,
            center=(0, 5),
            field=PerfectFieldLayer,
        )
        img1000 = generate_open_field(
            field_size=(100, 100), imager=AS1000Image, field=PerfectFieldLayer
        )
        p1200 = img1200.image[640, :]
        p1000 = img1000.image[384, :]
        p1200_prof = FWXMProfilePhysical(values=p1200, dpmm=1 / img1200.pixel_size)
        p1000_prof = FWXMProfilePhysical(values=p1000, dpmm=1 / img1000.pixel_size)
        gamma = p1000_prof.gamma(
            reference_profile=p1200_prof,
            dose_to_agreement=1,
            gamma_cap_value=2,
            centering=Centering.BEAM_CENTER,
        )
        # gamma is very low; just pixel noise from the image generator
        self.assertAlmostEqual(np.nanmean(gamma), 0.003, delta=0.01)

    def test_beam_centering_manual(self):
        # 5mm offset but we don't correct for it so expect a high gamma
        img1200 = generate_open_field(
            field_size=(100, 100),
            imager=AS1200Image,
            center=(0, 5),
            field=PerfectFieldLayer,
        )
        img1000 = generate_open_field(
            field_size=(100, 100), imager=AS1000Image, field=PerfectFieldLayer
        )
        p1200 = img1200.image[640, :]
        p1000 = img1000.image[384, :]
        p1200_prof = FWXMProfilePhysical(values=p1200, dpmm=1 / img1200.pixel_size)
        p1000_prof = FWXMProfilePhysical(values=p1000, dpmm=1 / img1000.pixel_size)
        gamma = p1000_prof.gamma(
            reference_profile=p1200_prof,
            dose_to_agreement=1,
            gamma_cap_value=2,
            centering=Centering.MANUAL,
        )
        # gamma is relatively high because of the 5mm offset
        self.assertAlmostEqual(np.nanmean(gamma), 0.327, delta=0.01)
        self.assertEqual(np.nanmax(gamma), 2.0)
