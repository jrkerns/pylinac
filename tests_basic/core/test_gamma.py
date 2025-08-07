import math
from unittest import TestCase

import numpy as np
from parameterized import parameterized

from pylinac.core.gamma import _compute_distance, gamma_1d, gamma_2d, gamma_geometric
from pylinac.core.image import DicomImage
from pylinac.core.image_generator import AS1000Image, AS1200Image
from pylinac.core.profile import FWXMProfile, FWXMProfilePhysical
from tests_basic.core.test_profile import generate_open_field
from tests_basic.utils import get_file_from_cloud_test_repo


class TestAgnewMcGarry(TestCase):
    """Tests from the Agnew & McGarry paper. https://www.sciencedirect.com/science/article/abs/pii/S0167814015006660

    .. warning::
        Agnew/McGarry used the opposite notation as Low's, hence this test class requires swapping the data (ref -> eval, eval -> ref).
        Since publication, Agnew/McGarry's code has been deprecated and replaced by pymedphys.gamma which uses the conventional notation.
    """

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
            reference=eval_array,
            evaluation=ref_array,
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
            reference=eval_array,
            evaluation=ref_array,
            dose_to_agreement=1,
            distance_to_agreement=1,
            global_dose=True,
        )
        meas_pass_rate = np.sum(g < 1) / g.size * 100
        self.assertAlmostEqual(meas_pass_rate, PASS_RATE, delta=0.05)

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
            reference=eval_array,
            evaluation=ref_array,
            dose_to_agreement=1,
            distance_to_agreement=4,  # 0.25mm pixel pitch; 4*0.25 = 1mm
            global_dose=True,
        )
        meas_pass_rate = np.sum(g < 1) / g.size * 100  # 92.19%
        self.assertAlmostEqual(meas_pass_rate, PASS_RATE, delta=0.15)

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
            reference=eval_array,
            evaluation=ref_array,
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
        eval[(-1, -1, -2, -2), (-1, -2, -2, -1)] = (
            1.03  # set bottom right corner to 3% off
        )
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

    def test_different_sizes(self):
        ref = np.ones((5, 5))
        eval = np.ones((6, 6))
        gamma = gamma_2d(reference=ref, evaluation=eval)
        self.assertEqual(gamma.shape, ref.shape)


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
        eval = np.asarray((1.03, 1.03, 1.03, 1.03, 1))
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
        eval = np.asarray((1, 1.03, 1.03, 1.03, 1.03))
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

    @parameterized.expand(
        [
            (1, 50),
            (2, 50),
            (3, 100),
            (4, 100),
            (5, 100),
        ]
    )
    def test_distance_scaling(self, dist: float, pass_rate: float):
        # see RAM-3754
        eval = np.asarray([0, 0.5, 1, 1, 1, 1])
        # ref is 3mm off. We want to test that gamma is <100 <3mm, then goes to 100 when we reach >=3mm
        ref = np.asarray([0, 0, 0, 0, 0.5, 1])
        # for dist in range(1, 5):
        gamma = gamma_geometric(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=0.1,  # large dose to ensure it's not the limiting factor
            distance_to_agreement=dist,
        )
        valid_elements = gamma[~np.isnan(gamma)]
        measured_pass_rate = 100 * (np.sum(valid_elements < 1) / len(valid_elements))
        self.assertAlmostEqual(measured_pass_rate, pass_rate, delta=1)

    def test_positive_distance_to_agreement(self):
        ref = np.ones(5)
        eval = np.ones(5)
        with self.assertRaises(ValueError):
            gamma_geometric(
                reference=ref,
                evaluation=eval,
                dose_to_agreement=1,
                distance_to_agreement=-1,
            )

    def test_postive_dose_to_agreement(self):
        ref = np.ones(5)
        eval = np.ones(5)
        with self.assertRaises(ValueError):
            gamma_geometric(
                reference=ref,
                evaluation=eval,
                dose_to_agreement=-1,
                distance_to_agreement=1,
            )

    def test_very_far_spacings(self):
        # if the measurement spacings are very far apart compared to the distance to
        # agreement, the indices for computing the simplex can end up being a single point which will error oout
        ref = np.asarray([0, 1, 3, 1, 0])
        eval = np.asarray([0, 1, 3, 1, 0])
        ref_x = np.arange(5)
        eval_x = np.arange(5)
        g = gamma_geometric(
            reference=ref,
            reference_coordinates=ref_x,
            evaluation=eval,
            evaluation_coordinates=eval_x,
            distance_to_agreement=0.1,
        )
        self.assertEqual(np.nanmax(g), 0)

    def test_reference_x_domain_smaller_than_eval(self):
        """Even if the reference x-domain is too small we can still
        evaluate the gamma."""
        ref_vals = [0, 0, 0, 0, 1, 2, 5, 8, 10, 10, 10, 10, 10, 8, 5, 2, 1, 0, 0, 0, 0]
        x_ref_vals = np.arange(len(ref_vals))
        vals = ref_vals[3:-3]  # we short-change the reference in low-dose areas
        x_vals = x_ref_vals[3:-3]
        gamma = gamma_geometric(
            reference=np.array(ref_vals),
            reference_coordinates=np.array(x_ref_vals),
            evaluation=np.array(vals),
            evaluation_coordinates=np.array(x_vals),
            distance_to_agreement=1,
            gamma_cap_value=2,
            dose_threshold=0,  # 0 threshold is important to ensure we calculate the gamma at the edges
        )
        self.assertEqual(np.nanmax(gamma), 2)
        self.assertAlmostEqual(np.nanmean(gamma), 0.476, places=2)

    @parameterized.expand(
        [
            (np.arange(5), np.arange(5), [np.nan, 0, 0, 0, 0]),
            (np.arange(4, -1, -1), np.arange(4, -1, -1), [0, 0, 0, 0, np.nan]),
        ],
    )
    def test_at_left_edge(self, ref, eval, expected):
        # test when the evaluation is at the edge of the reference
        # ensures no bounds error
        g = gamma_geometric(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=1,
            distance_to_agreement=1,
        )
        self.assertTrue(np.array_equal(g, expected, equal_nan=True))

    @parameterized.expand([(np.arange(4, -5, -1),), (np.arange(-4, 5, 1),)])
    def test_reversed_x_values_have_same_result(self, x):
        # test when the x values are reversed
        eval = np.asarray([0, 1, 3, 4, 5, 4, 4, 1, 0])
        ref = np.asarray([0, 1, 3, 4, 5, 4, 3, 1, 0])
        g = gamma_geometric(
            reference=ref,
            reference_coordinates=x,
            evaluation=eval,
            evaluation_coordinates=x,
            dose_to_agreement=1,
            distance_to_agreement=1,
        )
        self.assertTrue(
            np.allclose(
                g,
                [np.nan, 0, 0, 0, 0, 0, 0.3332, 0, np.nan],
                equal_nan=True,
                atol=0.001,
            )
        )

    @parameterized.expand([(np.arange(4, -5, -1),), (np.arange(-4, 5, 1),)])
    def test_reversed_x_values_for_one_profile(self, x):
        # test that when one profile is reversed, the gamma is still the same as if they were both in the same order
        eval = np.asarray([0, 1, 3, 4, 5, 4, 4, 1, 0])
        ref = np.asarray([0, 1, 3, 4, 5, 4, 3, 1, 0])
        g = gamma_geometric(
            reference=ref,
            reference_coordinates=x,
            evaluation=eval,
            evaluation_coordinates=np.flip(x),
            dose_to_agreement=1,
            distance_to_agreement=1,
        )
        self.assertTrue(
            np.allclose(
                g,
                [np.nan, 0, 0.3332, 0, 0, 0, 0, 0, np.nan],
                equal_nan=True,
                atol=0.001,
            )
        )

    @parameterized.expand(
        [
            ("bad eval", [0, 1, 2, 3, 4], [1, 3, 4, 3, 6]),
            ("bad ref", [1, 3, 4, 3, 6], [0, 1, 2, 3, 4]),
            ("bad both", [1, 3, 4, 3, 6], [1, 3, 4, 3, 6]),
        ]
    )
    def test_non_monotonic_fails(self, _, ref_x, eval_x):
        ref = eval = np.ones(5)
        with self.assertRaises(ValueError) as e:
            gamma_geometric(
                reference=ref,
                reference_coordinates=np.array(ref_x),
                evaluation=eval,
                evaluation_coordinates=np.array(eval_x),
            )
        self.assertIn("monotonically", str(e.exception))


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
        eval = np.asarray((1.03, 1.03, 1.03, 1, 1))
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
        eval = np.asarray((1, 1, 1.03, 1.03, 1.03))
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

    def test_different_sizes(self):
        ref = np.ones(5)
        eval = np.ones(6)
        gamma, _, _ = gamma_1d(
            reference=ref, evaluation=eval, dose_to_agreement=1, gamma_cap_value=2
        )
        self.assertEqual(len(gamma), len(ref))


class TestGammaFromProfile(TestCase):
    """Test the gamma method from the physical profile class"""

    def test_other_profile_not_physical(self):
        p = FWXMProfilePhysical(values=np.ones(5), dpmm=1)
        p2 = FWXMProfile(values=np.ones(5))
        with self.assertRaises(ValueError):
            p.gamma(evaluation_profile=p2)

    def test_gamma_on_self_is_zero(self):
        p = FWXMProfilePhysical(values=np.ones(5), dpmm=1)
        gamma = p.gamma(evaluation_profile=p)
        self.assertTrue(np.allclose(gamma, 0))

    def test_dose_gamma(self):
        ref_vals = np.ones(5) * 1.03
        vals = np.ones(5)
        p = FWXMProfilePhysical(values=vals, dpmm=1)
        ref_p = FWXMProfilePhysical(values=ref_vals, dpmm=1)
        gamma = p.gamma(
            evaluation_profile=ref_p, dose_to_agreement=3, gamma_cap_value=2
        )
        # all right at 1 because dose is 3% high
        self.assertTrue(np.allclose(gamma, 1))

    def test_low_density_eval(self):
        # the eval is 1/2 the spacing, but same values, as the reference
        ref_vals = np.arange(1, 11)[::2]
        vals = np.arange(11)
        p = FWXMProfilePhysical(values=ref_vals, x_values=np.arange(1, 11, 2))
        eval_p = FWXMProfilePhysical(values=vals)
        gamma = p.gamma(
            evaluation_profile=eval_p, dose_to_agreement=1, gamma_cap_value=2
        )
        # gamma is perfect
        self.assertTrue(np.allclose(gamma, 0))

    def test_different_epids(self):
        """This test the same profile but with different EPIDs (i.e. pixel size)"""
        img1200 = generate_open_field(field_size=(100, 100), imager=AS1200Image)
        img1000 = generate_open_field(field_size=(100, 100), imager=AS1000Image)
        p1200 = img1200.image[640, :]
        p1000 = img1000.image[384, :]
        p1200_prof = FWXMProfilePhysical(values=p1200, dpmm=1 / img1200.pixel_size)
        p1000_prof = FWXMProfilePhysical(values=p1000, dpmm=1 / img1000.pixel_size)
        gamma = p1000_prof.gamma(
            evaluation_profile=p1200_prof, dose_to_agreement=1, gamma_cap_value=2
        )
        # gamma is very low; just pixel noise from the image generator
        self.assertLessEqual(np.nanmean(gamma), 0.005)
