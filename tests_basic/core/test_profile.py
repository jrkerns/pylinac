import warnings
from typing import Type
from unittest import TestCase

import numpy as np
import scipy.signal as sps

from pylinac.core import image
from pylinac.core.array_utils import normalize
from pylinac.core.image_generator import (
    FilteredFieldLayer,
    FilterFreeFieldLayer,
    GaussianFilterLayer,
    PerfectFieldLayer,
)
from pylinac.core.image_generator.layers import Layer
from pylinac.core.image_generator.simulators import Simulator
from pylinac.core.profile import (
    CircleProfile,
    CollapsedCircleProfile,
    FWXMProfile,
    FWXMProfilePhysical,
    HillProfile,
    HillProfilePhysical,
    InflectionDerivativeProfile,
    InflectionDerivativeProfilePhysical,
    Interpolation,
    MultiProfile,
    Normalization,
    SingleProfile,
    gamma_1d,
    stretch,
)
from pylinac.metrics.profile import (
    PDD,
    Dmax,
    FlatnessDifferenceMetric,
    FlatnessRatioMetric,
    PenumbraLeftMetric,
    PenumbraRightMetric,
    SlopeMetric,
    SymmetryAreaMetric,
    SymmetryPointDifferenceMetric,
    SymmetryPointDifferenceQuotientMetric,
    TopDistanceMetric,
)
from tests_basic.utils import get_file_from_cloud_test_repo


def generate_open_field(
    field_size=(100, 100),
    sigma=2,
    center=(0, 0),
    field: Type[Layer] = FilteredFieldLayer,
) -> Simulator:
    from pylinac.core.image_generator import AS1000Image

    as1000 = AS1000Image()  # this will set the pixel size and shape automatically
    as1000.add_layer(
        field(field_size_mm=field_size, cax_offset_mm=center)
    )  # create a 50x50mm square field
    as1000.add_layer(
        GaussianFilterLayer(sigma_mm=sigma)
    )  # add an image-wide gaussian to simulate penumbra/scatter
    return as1000


def generate_profile(
    field_size=100, sigma=2, center=0, field: Type[Layer] = FilteredFieldLayer
) -> np.ndarray:
    img = generate_open_field(
        field_size=(field_size, field_size),
        sigma=sigma,
        center=(center, center),
        field=field,
    ).image
    arr = normalize(img[:, img.shape[1] // 2])
    return arr


class TestGamma1D(TestCase):
    def test_same_profile_is_0_gamma(self):
        ref = eval = np.ones(5)
        gamma = gamma_1d(reference=ref, evaluation=eval)
        self.assertEqual(max(gamma), 0)
        self.assertEqual(min(gamma), 0)
        self.assertEqual(len(gamma), 5)

        # test a high measurement value
        ref = eval = np.ones(5) * 50
        gamma = gamma_1d(reference=ref, evaluation=eval)
        self.assertEqual(max(gamma), 0)
        self.assertEqual(min(gamma), 0)
        self.assertEqual(len(gamma), 5)

    def test_gamma_perfectly_at_1(self):
        # offset a profile exactly by the dose to agreement
        ref = np.ones(5)
        eval = np.ones(5) * 1.01
        gamma = gamma_1d(reference=ref, evaluation=eval, dose_to_agreement=1)
        self.assertAlmostEqual(max(gamma), 1, delta=0.001)
        self.assertAlmostEqual(min(gamma), 1, delta=0.001)

        # test same but eval is LOWER than ref
        ref = np.ones(5)
        eval = np.ones(5) * 0.99
        gamma = gamma_1d(reference=ref, evaluation=eval, dose_to_agreement=1)
        self.assertAlmostEqual(max(gamma), 1, delta=0.001)
        self.assertAlmostEqual(min(gamma), 1, delta=0.001)

    def test_gamma_half(self):
        # offset a profile by half the dose to agreement to ensure it's 0.5
        ref = np.ones(5)
        eval = np.ones(5) / 1.005
        gamma = gamma_1d(reference=ref, evaluation=eval, dose_to_agreement=1)
        self.assertAlmostEqual(max(gamma), 0.5, delta=0.01)
        self.assertAlmostEqual(min(gamma), 0.5, delta=0.01)

    def test_gamma_some_on_some_off(self):
        ref = np.ones(5)
        eval = np.asarray((1.03, 1.03, 1, 1, 1))
        gamma = gamma_1d(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=1,
            distance_to_agreement=1,
            gamma_cap_value=5,
        )
        self.assertAlmostEqual(gamma[0], 3, delta=0.01)  # fully off by 3
        self.assertAlmostEqual(
            gamma[1], 1, delta=0.01
        )  # dose at next pixel matches (dose=0, dist=1)
        self.assertAlmostEqual(gamma[-1], 0, delta=0.01)  # gamma at end is perfect

        # check inverted pattern is mirrored (checks off-by-one errors)
        ref = np.ones(5)
        eval = np.asarray((1, 1, 1, 1.03, 1.03))
        gamma = gamma_1d(
            reference=ref,
            evaluation=eval,
            dose_to_agreement=1,
            distance_to_agreement=1,
            gamma_cap_value=5,
        )
        self.assertAlmostEqual(gamma[0], 0, delta=0.01)
        self.assertAlmostEqual(gamma[-2], 1, delta=0.01)
        self.assertAlmostEqual(gamma[-1], 3, delta=0.01)

    def test_localized_dose(self):
        ref = eval = np.array((100, 1, 1, 1, 1))
        gamma = gamma_1d(
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
        gamma = gamma_1d(
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
        gamma = gamma_1d(
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
        gamma = gamma_1d(
            reference=ref, evaluation=eval, dose_to_agreement=1, gamma_cap_value=2
        )
        self.assertEqual(max(gamma), 2)
        self.assertEqual(min(gamma), 2)

    def test_non_1d_array(self):
        ref = np.ones(5)
        eval = np.ones((5, 5))
        with self.assertRaises(ValueError):
            gamma_1d(reference=ref, evaluation=eval)


def create_simple_9_profile() -> np.array:
    # length of 9
    return np.array([0, 1, 2, 3, 4, 3, 2, 1, 0], dtype=float)


def create_simple_8_profile() -> np.array:
    # length of 8
    return np.array([0, 1, 2, 3, 3, 2, 1, 0], dtype=float)


def create_long_23_profile() -> np.array:
    # length of 23
    return np.array(
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        dtype=float,
    )


def create_long_22_profile() -> np.array:
    # length of 22
    return np.array(
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        dtype=float,
    )


def skewed_19_profile() -> np.array:
    """A profile where the peak is skewed to the right."""
    # length of 19
    return np.array(
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 8, 6, 4, 2, 0], dtype=float
    )


def symmetrical_sigmoidal_21_profile() -> np.array:
    """A curve with sigmoid shape on either side of the center"""
    # length of 21
    return np.array(
        [0, 1, 2, 4, 6, 8, 9, 10, 10, 10, 10, 10, 10, 10, 9, 8, 6, 4, 2, 1, 0],
        dtype=float,
    )


def symmetrical_sigmoidal_20_profile() -> np.array:
    """A curve with sigmoid shape on either side of the center"""
    # length of 20
    return np.array(
        [0, 1, 2, 4, 6, 8, 9, 10, 10, 10, 10, 10, 10, 9, 8, 6, 4, 2, 1, 0], dtype=float
    )


def symmetrical_sharp_sigmoidal_21_profile() -> np.array:
    """A curve with sharper sigmoid shape on either side of the center"""
    # length of 21
    return np.array(
        [0, 1, 1, 2, 5, 8, 9, 10, 10, 10, 10, 10, 10, 10, 9, 8, 5, 2, 1, 1, 0],
        dtype=float,
    )


class TestProfileGeneric(TestCase):
    def test_ground(self) -> None:
        offset_array = create_simple_9_profile() + 1
        self.assertEqual(offset_array.min(), 1)
        profile = FWXMProfile(offset_array, ground=True)
        # grounding makes min value 0
        self.assertEqual(profile.values.min(), 0)

    def test_max_normalization(self):
        array = create_simple_9_profile()
        self.assertNotEqual(array.max(), 1)
        max_prof = FWXMProfile(array, normalization=Normalization.MAX)
        self.assertEqual(max_prof.values.max(), 1)

    def test_no_normalization(self):
        array = create_simple_9_profile()
        no_prof = FWXMProfile(array, normalization=Normalization.NONE)
        self.assertEqual(no_prof.values.max(), array.max())

    def test_beam_normalization(self):
        array = create_long_23_profile()
        # dip the center so we don't conflate w/ max normalization
        array[11] = 8
        prof = FWXMProfile(array, normalization=Normalization.BEAM_CENTER)
        self.assertEqual(prof.values.max(), array.max() / 8)

    def test_geometric_normalization(self):
        array = create_long_23_profile()
        array = np.roll(array, 3)  # shift array so beam center is not at geom center
        prof = FWXMProfile(array, normalization=Normalization.GEOMETRIC_CENTER)
        self.assertNotEqual(prof.values.max(), array.max())

    def test_resample_with_ints_and_small_range_raises_warning(self):
        array = create_long_23_profile().astype(np.uint8)
        prof = FWXMProfile(array)
        with warnings.catch_warnings(record=True) as w:
            prof.as_resampled(interpolation_factor=2)
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[-1].category, UserWarning))
            self.assertIn("small", str(w[-1].message))

    def test_physical_resample_with_ints_and_small_range_raises_warning(self):
        array = create_long_23_profile().astype(np.uint8)
        prof = FWXMProfilePhysical(array, 1)
        with warnings.catch_warnings(record=True) as w:
            prof.as_resampled(interpolation_resolution_mm=0.1)
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[-1].category, UserWarning))
            self.assertIn("small", str(w[-1].message))


class TestFWXMProfile(TestCase):
    def test_x_values_decrease_is_okay(self):
        array = create_simple_9_profile()
        x_values = np.arange(len(array))[::-1]
        f = FWXMProfile(array, x_values=x_values)
        # both width and indices should respect the x-values
        self.assertEqual(f.field_width_px, 4)
        (left, right, width) = f.field_indices(in_field_ratio=0.8)
        self.assertEqual((left, right, width), (2, 6, 4))
        assert np.allclose(
            f.field_x_values(in_field_ratio=0.8), np.array([2, 3, 4, 5, 6])
        )

    def test_x_values_about_zero(self):
        """A profile that's centered on x should be correct"""
        array = create_long_23_profile()
        x_values = np.arange(len(array)) - 22 / 2
        f = FWXMProfile(array, x_values=x_values)
        (left, right, width) = f.field_indices(in_field_ratio=0.8)
        # the edges are symmetric about 0 because we offset the x-values above.
        self.assertEqual((left, right, width), (-5, 5, 10))
        # the x-values are symmetric about 0 because we offset the x-values above.
        assert np.allclose(
            f.field_x_values(in_field_ratio=0.8), np.arange(-5, 5.1), atol=0.01
        )

    def test_non_integer_increment_x_values_has_right_num_indices(self):
        """Arrays with x-values that are not integers used to cause problems when calculating the field array"""
        array = create_simple_9_profile()
        # x-values aren't simple integers
        x_values = (np.arange(len(array)) - 4.5) / 0.75
        f = FWXMProfile(array, x_values=x_values)
        self.assertAlmostEqual(f.field_width_px, 5.33, delta=0.01)
        (left, right, center) = f.field_indices(in_field_ratio=0.8)
        self.assertEqual((left, right, center), (-2, 2, 4))
        # the x-values are symmetric about 0 because we offset the x-values above.
        assert np.allclose(
            f.field_x_values(in_field_ratio=0.8),
            np.array([-2, -0.667, 0.6667, 2]),
            atol=0.01,
        )
        self.assertEqual(len(f.field_values(in_field_ratio=0.8)), 4)

    def test_not_monotonically_increasing_raises_error(self):
        array = create_simple_9_profile()
        x_values = np.arange(len(array))
        x_values[2] = 5
        with self.assertRaises(ValueError):
            FWXMProfile(array, x_values=x_values)

    def test_not_monotonically_decreasing_raises_error(self):
        array = create_simple_9_profile()
        x_values = np.arange(len(array))[::-1]
        x_values[5] = 5
        with self.assertRaises(ValueError):
            FWXMProfile(array, x_values=x_values)

    def test_center_idx(self):
        array = create_simple_9_profile()
        profile = FWXMProfile(array)
        self.assertEqual(profile.center_idx, 4)

    def test_center_even_idx(self):
        array = create_simple_8_profile()
        profile = FWXMProfile(array)
        self.assertEqual(profile.center_idx, 3.5)

    def test_field_edge_idx_50(self):
        array = create_simple_9_profile()
        profile = FWXMProfile(array, fwxm_height=50)
        self.assertEqual(profile.field_edge_idx("left"), 2)
        self.assertEqual(profile.field_edge_idx("right"), 6)

    def test_field_edge_idx_50_even(self):
        array = create_simple_8_profile()
        profile = FWXMProfile(array, fwxm_height=50)
        self.assertEqual(profile.field_edge_idx("left"), 1.5)
        self.assertEqual(profile.field_edge_idx("right"), 5.5)

    def test_field_edge_idx_25(self):
        array = create_simple_9_profile()
        profile = FWXMProfile(array, fwxm_height=25)
        self.assertEqual(profile.field_edge_idx("left"), 1)
        self.assertEqual(profile.field_edge_idx("right"), 7)

    def test_field_edge_idx_75(self):
        array = create_simple_9_profile()
        profile = FWXMProfile(array, fwxm_height=75)
        self.assertEqual(profile.field_edge_idx("left"), 3)
        self.assertEqual(profile.field_edge_idx("right"), 5)

    def test_field_edge_skewed(self):
        array = skewed_19_profile()
        profile = FWXMProfile(array, fwxm_height=50)
        self.assertEqual(profile.field_edge_idx("left"), 5)
        self.assertEqual(profile.field_edge_idx("right"), 14.5)

    def test_field_width_50(self):
        array = create_simple_9_profile()
        profile = FWXMProfile(array, fwxm_height=50)
        self.assertEqual(profile.field_width_px, 4)

    def test_field_width_25(self):
        array = create_simple_9_profile()
        profile = FWXMProfile(array, fwxm_height=25)
        self.assertEqual(profile.field_width_px, 6)

    def test_field_width_skewed(self):
        array = skewed_19_profile()
        profile = FWXMProfile(array, fwxm_height=50)
        self.assertEqual(profile.field_width_px, 9.5)

    def test_field_values_full_width(self):
        # ratio of 1 is FWHM
        array = create_simple_9_profile()
        profile = FWXMProfile(array, fwxm_height=50)
        field_values = profile.field_values(in_field_ratio=1)
        self.assertIsInstance(field_values, np.ndarray)
        self.assertEqual(len(field_values), 5)

    def test_field_values_half_width(self):
        # ratio of 0.5 is FHWM
        array = create_simple_9_profile()
        profile = FWXMProfile(array, fwxm_height=50)
        field_values = profile.field_values(in_field_ratio=0.5)
        self.assertIsInstance(field_values, np.ndarray)
        self.assertEqual(len(field_values), 3)

    def test_resample(self):
        array = create_long_23_profile()
        profile = FWXMProfile(array, fwxm_height=50)
        resampled_profile = profile.as_resampled(interpolation_factor=2)
        self.assertEqual(len(resampled_profile), len(profile) * 2)
        self.assertIsInstance(resampled_profile, FWXMProfile)
        # ensure x-values are the same; i.e. that we didn't just multiple x-values
        self.assertEqual(resampled_profile.x_values.max(), profile.x_values.max())
        # y values should be similar.
        self.assertAlmostEqual(
            resampled_profile.values.max(), profile.values.max(), delta=0.1
        )

    def test_resample_10(self):
        array = create_long_23_profile()
        profile = FWXMProfile(array, fwxm_height=50)
        resampled_profile = profile.as_resampled(interpolation_factor=10)
        self.assertEqual(len(resampled_profile), len(profile) * 10)
        self.assertIsInstance(resampled_profile, FWXMProfile)
        # ensure x-values are the same; i.e. that we didn't just multiple x-values
        self.assertEqual(resampled_profile.x_values.max(), profile.x_values.max())
        # y values should be similar.
        self.assertAlmostEqual(
            resampled_profile.values.max(), profile.values.max(), delta=0.1
        )

    def test_resample_in_half(self):
        array = create_long_23_profile()
        profile = FWXMProfile(array, fwxm_height=50)
        resampled_profile = profile.as_resampled(interpolation_factor=0.5)
        self.assertEqual(len(resampled_profile), 12)
        self.assertIsInstance(resampled_profile, FWXMProfile)
        # ensure x-values are the same; i.e. that we didn't just multiple x-values
        self.assertEqual(resampled_profile.x_values.max(), profile.x_values.max())
        # y values should be similar.
        self.assertAlmostEqual(
            resampled_profile.values.max(), profile.values.max(), delta=0.1
        )


class TestInflectionDerivativeProfile(TestCase):
    def test_center_idx(self):
        array = symmetrical_sigmoidal_21_profile()
        profile = InflectionDerivativeProfile(array)
        self.assertAlmostEqual(profile.center_idx, 10, delta=0.01)

    def test_center_idx_even(self):
        array = symmetrical_sigmoidal_20_profile()
        profile = InflectionDerivativeProfile(array)
        self.assertAlmostEqual(profile.center_idx, 9.5, delta=0.01)

    def test_field_edge_idx(self):
        array = symmetrical_sigmoidal_21_profile()
        profile = InflectionDerivativeProfile(array)
        self.assertAlmostEqual(profile.field_edge_idx("left"), 3.5, delta=0.01)
        self.assertAlmostEqual(profile.field_edge_idx("right"), 16.5, delta=0.01)

    def test_field_width(self):
        array = symmetrical_sigmoidal_21_profile()
        profile = InflectionDerivativeProfile(array)
        self.assertAlmostEqual(profile.field_width_px, 13, delta=0.01)

    def test_resample_10(self):
        array = create_long_23_profile()
        profile = InflectionDerivativeProfile(array)
        resampled_profile = profile.as_resampled(interpolation_factor=10)
        self.assertEqual(len(resampled_profile), len(profile) * 10)
        self.assertIsInstance(resampled_profile, InflectionDerivativeProfile)
        # ensure x-values are the same; i.e. that we didn't just multiple x-values
        self.assertEqual(resampled_profile.x_values.max(), profile.x_values.max())
        # y values should be similar.
        self.assertAlmostEqual(
            resampled_profile.values.max(), profile.values.max(), delta=0.1
        )

    def test_resampled_half(self):
        array = create_long_23_profile()
        profile = InflectionDerivativeProfile(array)
        resampled_profile = profile.as_resampled(interpolation_factor=0.5)
        self.assertEqual(len(resampled_profile), 12)
        self.assertIsInstance(resampled_profile, InflectionDerivativeProfile)
        # ensure x-values are the same; i.e. that we didn't just multiple x-values
        self.assertEqual(resampled_profile.x_values.max(), profile.x_values.max())
        # y values should be similar.
        self.assertAlmostEqual(
            resampled_profile.values.max(), profile.values.max(), delta=0.1
        )


class TestHillProfile(TestCase):
    def test_center_idx(self):
        array = symmetrical_sharp_sigmoidal_21_profile()
        profile = HillProfile(array, hill_window_ratio=0.2)
        self.assertAlmostEqual(profile.center_idx, 10, delta=0.1)

    def test_field_edge_idx(self):
        array = symmetrical_sharp_sigmoidal_21_profile()
        profile = HillProfile(array, hill_window_ratio=0.2)
        self.assertAlmostEqual(profile.field_edge_idx("left"), 3.8, delta=0.1)
        self.assertAlmostEqual(profile.field_edge_idx("right"), 15.9, delta=0.1)

    def test_field_width(self):
        array = symmetrical_sharp_sigmoidal_21_profile()
        profile = HillProfile(array, hill_window_ratio=0.2)
        self.assertAlmostEqual(profile.field_width_px, 12, delta=0.1)

    def test_resample_10(self):
        array = create_long_23_profile()
        profile = HillProfile(array, hill_window_ratio=0.2)
        resampled_profile = profile.as_resampled(interpolation_factor=10)
        self.assertEqual(len(resampled_profile), len(profile) * 10)
        self.assertIsInstance(resampled_profile, InflectionDerivativeProfile)
        # ensure x-values are the same; i.e. that we didn't just multiple x-values
        self.assertEqual(resampled_profile.x_values.max(), profile.x_values.max())
        # y values should be similar.
        self.assertAlmostEqual(
            resampled_profile.values.max(), profile.values.max(), delta=0.1
        )

    def test_resampled_half(self):
        array = create_long_23_profile()
        profile = HillProfile(array)
        resampled_profile = profile.as_resampled(interpolation_factor=0.5)
        self.assertEqual(len(resampled_profile), 12)
        self.assertIsInstance(resampled_profile, InflectionDerivativeProfile)
        # ensure x-values are the same; i.e. that we didn't just multiple x-values
        self.assertEqual(resampled_profile.x_values.max(), profile.x_values.max())
        # y values should be similar.
        self.assertAlmostEqual(
            resampled_profile.values.max(), profile.values.max(), delta=0.1
        )


class TestFWXMProfilePhysical(TestCase):
    def test_field_width_mm(self):
        array = create_long_23_profile()
        profile = FWXMProfilePhysical(array, fwxm_height=50, dpmm=2)
        width_px = profile.field_width_px
        self.assertEqual(profile.field_width_mm, width_px / 2)

    def test_resample_same_size(self):
        array = create_long_23_profile()
        profile = FWXMProfilePhysical(array, fwxm_height=50, dpmm=2)
        resampled_profile = profile.as_resampled(interpolation_resolution_mm=0.5)
        # 2 dpmm == 0.5mm resolution, so no change
        self.assertEqual(len(resampled_profile), len(profile))
        self.assertIsInstance(resampled_profile, FWXMProfilePhysical)
        # ensure x-values are the same; i.e. that we didn't just multiple x-values
        self.assertEqual(resampled_profile.x_values.max(), profile.x_values.max())
        # y values should be similar.
        self.assertAlmostEqual(
            resampled_profile.values.max(), profile.values.max(), delta=0.1
        )

    def test_resample_doubled(self):
        array = create_long_23_profile()
        profile = FWXMProfilePhysical(array, fwxm_height=50, dpmm=1)
        resampled_profile = profile.as_resampled(interpolation_resolution_mm=0.1)
        # will be 10x as large
        self.assertEqual(len(resampled_profile), len(profile) * 10)
        self.assertIsInstance(resampled_profile, FWXMProfilePhysical)
        # ensure x-values are properly offset
        self.assertEqual(
            resampled_profile.x_values.max(), profile.x_values.max() + 0.45
        )
        # y values should be similar.
        self.assertAlmostEqual(
            resampled_profile.values.max(), profile.values.max(), delta=0.1
        )
        self.assertEqual(resampled_profile.dpmm, 10)

    def test_resample_adds_half_pixel(self):
        """With phyiscal profiles, when interpolating,
        we have to account for the 'half pixel' offset that must be accounted for.
        See the grid_mode parameter of scikit-image zoom"""
        array = create_long_23_profile()
        profile = FWXMProfilePhysical(array, fwxm_height=50, dpmm=1)
        resampled_profile = profile.as_resampled(interpolation_resolution_mm=0.1)
        # Add half pixel to either side to account for physical size
        self.assertAlmostEqual(resampled_profile.x_values[0], -0.45, delta=0.01)
        self.assertAlmostEqual(resampled_profile.x_values[-1], 22.45, delta=0.01)

    def test_resample_of_resample(self):
        """Test that resampling twice **correctly** is the same as a 1-step resample"""
        array = create_long_23_profile()
        profile = FWXMProfilePhysical(array, fwxm_height=50, dpmm=1)
        prof10 = profile.as_resampled(interpolation_resolution_mm=0.1)
        prof100 = profile.as_resampled(interpolation_resolution_mm=0.01)
        # we need grid = false because we have already resampled the profile once
        # doing it again is not appropriate.
        prof100_2 = prof10.as_resampled(interpolation_resolution_mm=0.01, grid=False)
        self.assertEqual(len(prof100.values), len(prof100_2.values))
        self.assertEqual(prof100.center_idx, prof100_2.center_idx)
        self.assertAlmostEqual(
            prof100.field_width_px, prof100_2.field_width_px, delta=0.01
        )

    def test_physical_x_values(self):
        """Test converting dpmm to absolute x-values works and also includes the half-pixel on either end."""
        array = create_long_23_profile()
        dpmm = 3
        half_pixel = 1 / (dpmm * 2)
        profile = FWXMProfilePhysical(array, fwxm_height=50, dpmm=dpmm)
        phys_x_values = profile.physical_x_values
        self.assertAlmostEqual(min(phys_x_values), half_pixel)
        self.assertAlmostEqual(
            max(phys_x_values), 23 / dpmm - half_pixel
        )  # 23 is len of x-values

    def test_converting_to_simple_profile(self):
        pixel_size = 0.390625  # as1000 pixel size
        profile = generate_profile()
        fwxm = FWXMProfilePhysical(values=profile, dpmm=1 / pixel_size)
        abs_fwxm = fwxm.as_simple_profile()
        # ensure the field widths are the same
        self.assertAlmostEqual(fwxm.field_width_mm, abs_fwxm.field_width_px, delta=0.01)
        # ensure the CAXs are the same to ensure proper pixel offset
        self.assertAlmostEqual(
            fwxm.center_idx * pixel_size + pixel_size / 2,
            abs_fwxm.center_idx,
            delta=0.01,
        )
        self.assertIsInstance(abs_fwxm, FWXMProfile)


class TestInflectionProfilePhysical(TestCase):
    def test_field_width_mm(self):
        array = create_long_23_profile()
        profile = InflectionDerivativeProfilePhysical(array, dpmm=2)
        width_px = profile.field_width_px
        self.assertEqual(profile.field_width_mm, width_px / 2)

    def test_resample_same_size(self):
        array = create_long_23_profile()
        profile = InflectionDerivativeProfilePhysical(array, dpmm=2)
        resampled_profile = profile.as_resampled(interpolation_resolution_mm=0.5)
        # 2 dpmm == 0.5mm resolution, so no change
        self.assertEqual(len(resampled_profile), len(profile))
        self.assertIsInstance(resampled_profile, InflectionDerivativeProfilePhysical)
        # ensure x-values are the same; i.e. that we didn't just multiple x-values
        self.assertEqual(resampled_profile.x_values.max(), profile.x_values.max())
        # y values should be similar.
        self.assertAlmostEqual(
            resampled_profile.values.max(), profile.values.max(), delta=0.1
        )

    def test_resample_doubled(self):
        array = create_long_23_profile()
        profile = InflectionDerivativeProfilePhysical(array, dpmm=1)
        resampled_profile = profile.as_resampled(interpolation_resolution_mm=0.1)
        # will be 10x as large
        self.assertEqual(len(resampled_profile), len(profile) * 10)
        self.assertIsInstance(resampled_profile, InflectionDerivativeProfilePhysical)
        # ensure x-values are the nearly the same; i.e. that we didn't just multiple x-values
        # they aren't exactly the same due to physical pixel size
        self.assertEqual(
            resampled_profile.x_values.max(), profile.x_values.max() + 0.45
        )
        # y values should be similar.
        self.assertAlmostEqual(
            resampled_profile.values.max(), profile.values.max(), delta=0.1
        )

    def test_resample_adds_half_pixel(self):
        """With phyiscal profiles, when interpolating,
        we have to account for the 'half pixel' offset that must be accounted for.
        See the grid_mode parameter of scipy's zoom function"""
        array = create_long_23_profile()
        profile = InflectionDerivativeProfilePhysical(array, dpmm=1)
        resampled_profile = profile.as_resampled(interpolation_resolution_mm=0.1)
        # Add half pixel to either side to account for physical size
        self.assertAlmostEqual(resampled_profile.x_values[0], -0.45, delta=0.01)
        self.assertAlmostEqual(resampled_profile.x_values[-1], 22.45, delta=0.01)

    def test_resample_of_resample(self):
        """Test that resampling twice **correctly** is the same as a 1-step resample"""
        array = symmetrical_sharp_sigmoidal_21_profile()
        profile = InflectionDerivativeProfilePhysical(array, dpmm=1)
        prof10 = profile.as_resampled(interpolation_resolution_mm=0.1)
        prof100 = profile.as_resampled(interpolation_resolution_mm=0.01)
        # we need grid = false because we have already resampled the profile once
        # doing it again is not appropriate.
        prof100_2 = prof10.as_resampled(interpolation_resolution_mm=0.01, grid=False)
        self.assertEqual(len(prof100.values), len(prof100_2.values))
        self.assertAlmostEqual(prof100.center_idx, prof100_2.center_idx, delta=0.001)
        self.assertAlmostEqual(
            prof100.field_width_px, prof100_2.field_width_px, delta=0.001
        )


class TestHillProfilePhysical(TestCase):
    def test_field_width_mm(self):
        array = symmetrical_sharp_sigmoidal_21_profile()
        profile = HillProfilePhysical(array, dpmm=2, hill_window_ratio=0.2)
        width_px = profile.field_width_px
        self.assertEqual(profile.field_width_mm, width_px / 2)

    def test_resample_same_size(self):
        array = symmetrical_sharp_sigmoidal_21_profile()
        profile = HillProfilePhysical(array, dpmm=2, hill_window_ratio=0.2)
        resampled_profile = profile.as_resampled(interpolation_resolution_mm=0.5)
        # 2 dpmm == 0.5mm resolution, so no change
        self.assertEqual(len(resampled_profile), len(profile))
        self.assertIsInstance(resampled_profile, HillProfilePhysical)
        # ensure x-values are the same; i.e. that we didn't just multiple x-values
        self.assertEqual(resampled_profile.x_values.max(), profile.x_values.max())
        # y values should be similar.
        self.assertAlmostEqual(
            resampled_profile.values.max(), profile.values.max(), delta=0.1
        )

    def test_resample_doubled(self):
        array = symmetrical_sharp_sigmoidal_21_profile()
        profile = HillProfilePhysical(array, dpmm=1, hill_window_ratio=0.2)
        resampled_profile = profile.as_resampled(interpolation_resolution_mm=0.1)
        # will be 10x as large
        self.assertEqual(len(resampled_profile), len(profile) * 10)
        self.assertIsInstance(resampled_profile, HillProfilePhysical)
        # ensure x-values are the same plus offset; i.e. that we didn't just multiple x-values
        self.assertEqual(
            resampled_profile.x_values.max(), profile.x_values.max() + 0.45
        )
        # y values should be similar.
        self.assertAlmostEqual(
            resampled_profile.values.max(), profile.values.max(), delta=0.1
        )

    def test_resample_adds_half_pixel(self):
        """With phyiscal profiles, when interpolating,
        we have to account for the 'half pixel' offset that must be accounted for.
        See the grid_mode parameter of scikit-image zoom"""
        array = symmetrical_sharp_sigmoidal_21_profile()
        profile = HillProfilePhysical(array, dpmm=1)
        resampled_profile = profile.as_resampled(interpolation_resolution_mm=0.1)
        # Add half pixel to either side to account for physical size
        self.assertAlmostEqual(resampled_profile.x_values[0], -0.45, delta=0.01)
        self.assertAlmostEqual(resampled_profile.x_values[-1], 20.45, delta=0.01)

    def test_resample_of_resample(self):
        """Test that resampling twice **correctly** is the same as a 1-step resample"""
        array = symmetrical_sharp_sigmoidal_21_profile()
        profile = InflectionDerivativeProfilePhysical(array, dpmm=1)
        prof10 = profile.as_resampled(interpolation_resolution_mm=0.1)
        prof100 = profile.as_resampled(interpolation_resolution_mm=0.01)
        # we need grid = false because we have already resampled the profile once
        # doing it again is not appropriate.
        prof100_2 = prof10.as_resampled(interpolation_resolution_mm=0.01, grid=False)
        self.assertEqual(len(prof100.values), len(prof100_2.values))
        self.assertAlmostEqual(prof100.center_idx, prof100_2.center_idx, delta=0.001)
        self.assertAlmostEqual(
            prof100.field_width_px, prof100_2.field_width_px, delta=0.001
        )


class TestResampleTo(TestCase):
    def test_resample_same(self):
        array = symmetrical_sharp_sigmoidal_21_profile()
        profile = InflectionDerivativeProfilePhysical(array, dpmm=1)
        a2 = symmetrical_sigmoidal_21_profile()
        p2 = InflectionDerivativeProfilePhysical(a2, dpmm=1)
        new_p = profile.resample_to(p2)
        self.assertIsInstance(new_p, InflectionDerivativeProfile)
        self.assertEqual(len(new_p.x_values), len(profile.x_values))

    def test_extrapolate_not_allowed(self):
        array = create_simple_9_profile()
        source = InflectionDerivativeProfilePhysical(array, dpmm=1)
        a2 = symmetrical_sigmoidal_21_profile()
        target = InflectionDerivativeProfilePhysical(a2, dpmm=1)
        with self.assertRaises(ValueError):
            source.resample_to(target)

    def test_normal_subrange_physical(self):
        array = create_simple_9_profile()
        target = InflectionDerivativeProfilePhysical(array, dpmm=1)
        a2 = symmetrical_sigmoidal_21_profile()
        source = InflectionDerivativeProfilePhysical(a2, dpmm=1)
        new_profile = source.resample_to(target)
        self.assertIsInstance(new_profile, InflectionDerivativeProfile)
        self.assertEqual(len(new_profile.x_values), len(target.x_values))
        self.assertEqual(new_profile.x_values[3], target.physical_x_values[3])

    def test_normal_subrange_simple(self):
        array = create_simple_9_profile()
        target = InflectionDerivativeProfile(array)
        a2 = symmetrical_sigmoidal_21_profile()
        source = InflectionDerivativeProfile(a2)
        new_profile = source.resample_to(target)
        self.assertIsInstance(new_profile, InflectionDerivativeProfile)
        self.assertEqual(len(new_profile.x_values), len(target.x_values))
        self.assertEqual(new_profile.x_values[3], target.x_values[3])


def create_pdd_x_y() -> (np.array, np.array):
    x = [
        -1.65,
        -0.65,
        0.35,
        1.35,
        2.35,
        3.35,
        4.35,
        5.35,
        6.35,
        7.35,
        8.35,
        9.35,
        10.35,
        11.35,
        12.35,
        13.35,
        14.35,
        15.35,
        16.35,
        17.35,
        18.35,
        19.35,
        20.35,
        21.35,
        22.35,
        23.35,
        24.35,
        25.35,
        26.35,
        27.35,
        28.35,
        29.35,
        30.35,
        31.35,
        32.35,
        33.35,
        34.35,
        35.35,
        36.35,
        37.35,
        38.35,
        39.35,
        40.35,
        41.35,
        42.35,
        43.35,
        44.35,
        45.35,
        46.35,
        47.35,
        48.35,
        49.35,
        50.35,
        51.35,
        52.35,
        53.35,
        54.35,
        55.35,
        56.35,
        57.35,
        58.35,
        59.35,
        60.35,
        61.35,
        62.35,
        63.35,
        64.35,
        65.35,
        66.35,
        67.35,
        68.35,
        69.35,
        70.35,
        71.35,
        72.35,
        73.35,
        74.35,
        75.35,
        76.35,
        77.35,
        78.35,
        79.35,
        80.35,
        81.35,
        82.35,
        83.35,
        84.35,
        85.35,
        86.35,
        87.35,
        88.35,
        89.35,
        90.35,
        91.35,
        92.35,
        93.35,
        94.35,
        95.35,
        96.35,
        97.35,
        98.35,
        99.35,
        100.35,
        101.35,
        102.35,
        103.35,
        104.35,
        105.35,
        106.35,
        107.35,
        108.35,
        109.35,
        110.35,
        111.35,
        112.35,
        113.35,
        114.35,
        115.35,
        116.35,
        117.35,
        118.35,
        119.35,
        120.35,
        121.35,
        122.35,
        123.35,
        124.35,
        125.35,
        126.35,
        127.35,
        128.35,
        129.35,
        130.35,
        131.35,
        132.35,
        133.35,
        134.35,
        135.35,
        136.35,
        137.35,
        138.35,
        139.35,
        140.35,
        141.35,
        142.35,
        143.35,
        144.35,
        145.35,
        146.35,
        147.35,
        148.35,
        149.35,
        150.35,
        151.35,
        152.35,
        153.35,
        154.35,
        155.35,
        156.35,
        157.35,
        158.35,
        159.35,
        160.35,
        161.35,
        162.35,
        163.35,
        164.35,
        165.35,
        166.35,
        167.35,
        168.35,
        169.35,
        170.35,
        171.35,
        172.35,
        173.35,
        174.35,
        175.35,
        176.35,
        177.35,
        178.35,
        179.35,
        180.35,
        181.35,
        182.35,
        183.35,
        184.35,
        185.35,
        186.35,
        187.35,
        188.35,
        189.35,
        190.35,
        191.35,
        192.35,
        193.35,
        194.35,
        195.35,
        196.35,
        197.35,
        198.35,
        199.35,
        200.35,
        201.35,
        202.35,
        203.35,
        204.35,
        205.35,
        206.35,
        207.35,
        208.35,
        209.35,
        210.35,
        211.35,
        212.35,
        213.35,
        214.35,
        215.35,
        216.35,
        217.35,
        218.35,
        219.35,
        220.35,
        221.35,
        222.35,
        223.35,
        224.35,
        225.35,
        226.35,
        227.35,
        228.35,
        229.35,
        230.35,
        231.35,
        232.35,
        233.35,
        234.35,
        235.35,
        236.35,
        237.35,
        238.35,
        239.35,
        240.35,
        241.35,
        242.35,
        243.35,
        244.35,
        245.35,
        246.35,
        247.35,
        248.35,
        249.35,
        250.35,
        251.35,
        252.35,
        253.35,
        254.35,
        255.35,
        256.35,
        257.35,
        258.35,
        259.35,
        260.35,
        261.35,
        262.35,
        263.35,
        264.35,
        265.35,
        266.35,
        267.35,
        268.35,
        269.35,
        270.35,
        271.35,
        272.35,
        273.35,
        274.35,
        275.35,
        276.35,
        277.35,
        278.35,
        279.35,
        280.35,
        281.35,
        282.35,
        283.35,
        284.35,
        285.35,
        286.35,
        287.35,
        288.35,
        289.35,
        290.35,
        291.35,
        292.35,
        293.35,
        294.35,
        295.35,
        296.35,
        297.35,
        298.35,
        299.35,
        300.35,
        301.35,
        302.35,
        303.35,
        304.35,
        305.35,
        306.35,
        307.35,
        308.35,
    ]

    y = [
        1.30e00,
        1.35e00,
        1.44e00,
        1.57e00,
        1.83e00,
        2.10e00,
        2.29e00,
        2.44e00,
        2.56e00,
        2.68e00,
        2.75e00,
        2.83e00,
        2.89e00,
        2.93e00,
        2.97e00,
        3.01e00,
        3.03e00,
        3.06e00,
        3.08e00,
        3.10e00,
        3.10e00,
        3.11e00,
        3.13e00,
        3.11e00,
        3.12e00,
        3.12e00,
        3.12e00,
        3.12e00,
        3.10e00,
        3.10e00,
        3.07e00,
        3.08e00,
        3.06e00,
        3.06e00,
        3.04e00,
        3.02e00,
        3.01e00,
        3.00e00,
        2.99e00,
        2.98e00,
        2.96e00,
        2.95e00,
        2.94e00,
        2.92e00,
        2.92e00,
        2.90e00,
        2.89e00,
        2.87e00,
        2.87e00,
        2.85e00,
        2.84e00,
        2.83e00,
        2.81e00,
        2.80e00,
        2.78e00,
        2.76e00,
        2.77e00,
        2.75e00,
        2.74e00,
        2.71e00,
        2.69e00,
        2.69e00,
        2.68e00,
        2.67e00,
        2.65e00,
        2.65e00,
        2.63e00,
        2.61e00,
        2.60e00,
        2.59e00,
        2.59e00,
        2.56e00,
        2.55e00,
        2.54e00,
        2.54e00,
        2.52e00,
        2.51e00,
        2.49e00,
        2.48e00,
        2.47e00,
        2.45e00,
        2.45e00,
        2.43e00,
        2.42e00,
        2.40e00,
        2.40e00,
        2.39e00,
        2.37e00,
        2.36e00,
        2.34e00,
        2.34e00,
        2.32e00,
        2.32e00,
        2.30e00,
        2.28e00,
        2.28e00,
        2.26e00,
        2.26e00,
        2.24e00,
        2.24e00,
        2.23e00,
        2.21e00,
        2.21e00,
        2.18e00,
        2.18e00,
        2.16e00,
        2.15e00,
        2.14e00,
        2.13e00,
        2.13e00,
        2.10e00,
        2.10e00,
        2.09e00,
        2.07e00,
        2.07e00,
        2.06e00,
        2.05e00,
        2.04e00,
        2.02e00,
        2.02e00,
        2.01e00,
        2.00e00,
        1.99e00,
        1.98e00,
        1.97e00,
        1.96e00,
        1.95e00,
        1.94e00,
        1.93e00,
        1.92e00,
        1.91e00,
        1.89e00,
        1.89e00,
        1.88e00,
        1.87e00,
        1.87e00,
        1.85e00,
        1.84e00,
        1.83e00,
        1.83e00,
        1.82e00,
        1.81e00,
        1.80e00,
        1.78e00,
        1.78e00,
        1.77e00,
        1.75e00,
        1.75e00,
        1.75e00,
        1.73e00,
        1.72e00,
        1.71e00,
        1.71e00,
        1.70e00,
        1.69e00,
        1.69e00,
        1.67e00,
        1.66e00,
        1.65e00,
        1.65e00,
        1.64e00,
        1.62e00,
        1.63e00,
        1.61e00,
        1.60e00,
        1.60e00,
        1.59e00,
        1.58e00,
        1.58e00,
        1.56e00,
        1.56e00,
        1.56e00,
        1.54e00,
        1.53e00,
        1.53e00,
        1.52e00,
        1.51e00,
        1.51e00,
        1.50e00,
        1.49e00,
        1.48e00,
        1.48e00,
        1.47e00,
        1.46e00,
        1.45e00,
        1.45e00,
        1.44e00,
        1.43e00,
        1.43e00,
        1.42e00,
        1.41e00,
        1.40e00,
        1.39e00,
        1.39e00,
        1.38e00,
        1.38e00,
        1.36e00,
        1.36e00,
        1.35e00,
        1.35e00,
        1.34e00,
        1.33e00,
        1.33e00,
        1.31e00,
        1.31e00,
        1.31e00,
        1.30e00,
        1.29e00,
        1.29e00,
        1.28e00,
        1.27e00,
        1.27e00,
        1.26e00,
        1.26e00,
        1.25e00,
        1.24e00,
        1.24e00,
        1.23e00,
        1.22e00,
        1.22e00,
        1.22e00,
        1.21e00,
        1.20e00,
        1.19e00,
        1.19e00,
        1.18e00,
        1.18e00,
        1.17e00,
        1.17e00,
        1.16e00,
        1.16e00,
        1.15e00,
        1.15e00,
        1.14e00,
        1.13e00,
        1.12e00,
        1.12e00,
        1.11e00,
        1.11e00,
        1.10e00,
        1.10e00,
        1.09e00,
        1.09e00,
        1.08e00,
        1.08e00,
        1.07e00,
        1.07e00,
        1.06e00,
        1.05e00,
        1.05e00,
        1.05e00,
        1.04e00,
        1.03e00,
        1.03e00,
        1.02e00,
        1.02e00,
        1.02e00,
        1.01e00,
        1.01e00,
        1.00e00,
        9.95e-01,
        9.92e-01,
        9.82e-01,
        9.83e-01,
        9.74e-01,
        9.71e-01,
        9.66e-01,
        9.66e-01,
        9.60e-01,
        9.52e-01,
        9.49e-01,
        9.41e-01,
        9.38e-01,
        9.35e-01,
        9.33e-01,
        9.25e-01,
        9.19e-01,
        9.17e-01,
        9.09e-01,
        9.07e-01,
        9.05e-01,
        8.95e-01,
        8.93e-01,
        8.91e-01,
        8.85e-01,
        8.83e-01,
        8.79e-01,
        8.74e-01,
        8.69e-01,
        8.67e-01,
        8.64e-01,
        8.55e-01,
        8.52e-01,
        8.50e-01,
        8.50e-01,
        8.44e-01,
        8.35e-01,
        8.32e-01,
        8.30e-01,
        8.24e-01,
        8.21e-01,
        8.17e-01,
        8.11e-01,
        8.05e-01,
        8.07e-01,
        7.98e-01,
        8.00e-01,
        7.95e-01,
        7.90e-01,
        7.88e-01,
        7.81e-01,
    ]
    return np.asarray(x), np.asarray(y)


def create_electron_pdd_x_y() -> (np.ndarray, np.ndarray):
    ex = [
        -1,
        -0.5,
        0,
        0.5,
        1,
        1.5,
        2,
        2.5,
        3,
        3.5,
        4,
        4.5,
        5,
        5.5,
        6,
        6.5,
        7,
        7.5,
        8,
        8.5,
        9,
        9.5,
        10,
        10.5,
        11,
        11.5,
        12,
        12.5,
        13,
        13.5,
        14,
        14.5,
        15,
        15.5,
        16,
        16.5,
        17,
        17.5,
        18,
        18.5,
        19,
        19.5,
        20,
        20.5,
        21,
        21.5,
        22,
        22.5,
        23,
        23.5,
        24,
        24.5,
        25,
        25.5,
        26,
        26.5,
        27,
        27.5,
        28,
        28.5,
        29,
        29.5,
        30,
        30.5,
        31,
        31.5,
        32,
        32.5,
        33,
        33.5,
        34,
        34.5,
        35,
        35.5,
        36,
        36.5,
        37,
        37.5,
        38,
        38.5,
        39,
        39.5,
        40,
        40.5,
        41,
        41.5,
        42,
        42.5,
        43,
        43.5,
        44,
        44.5,
        45,
        45.5,
        46,
        46.5,
        47,
        47.5,
        48,
        48.5,
        49,
        49.5,
        50,
        50.5,
        51,
        51.5,
        52,
        52.5,
        53,
        53.1,
        54.1,
        55.1,
        56.1,
        57.1,
        58.1,
        59.1,
        60.1,
        61.1,
        62.1,
        63.1,
        64.1,
        65.1,
        66.1,
        67.1,
        68.1,
        69.1,
        70.1,
        71.1,
        72.1,
        73.1,
        73.2,
        75.2,
        77.2,
        79.2,
        81.2,
        83.2,
        85.2,
        85.6,
    ]

    ey = [
        0.669,
        0.67,
        0.67117,
        0.67228,
        0.67319,
        0.67554,
        0.67926,
        0.68323,
        0.68858,
        0.69343,
        0.69892,
        0.70448,
        0.70894,
        0.71315,
        0.71755,
        0.72171,
        0.72581,
        0.72951,
        0.73387,
        0.7379,
        0.74101,
        0.74494,
        0.74863,
        0.75244,
        0.75567,
        0.76046,
        0.76399,
        0.7675,
        0.77162,
        0.77551,
        0.77973,
        0.78298,
        0.78742,
        0.79076,
        0.79392,
        0.79776,
        0.80171,
        0.8038,
        0.80718,
        0.80833,
        0.81018,
        0.81278,
        0.81371,
        0.81423,
        0.81479,
        0.81546,
        0.81349,
        0.81206,
        0.8105,
        0.80642,
        0.80239,
        0.79823,
        0.79271,
        0.78586,
        0.7787,
        0.77018,
        0.76085,
        0.75013,
        0.73949,
        0.72737,
        0.71492,
        0.70111,
        0.68574,
        0.66938,
        0.65188,
        0.63321,
        0.615,
        0.59442,
        0.574,
        0.55177,
        0.52885,
        0.50554,
        0.48185,
        0.45682,
        0.43291,
        0.40886,
        0.3839,
        0.35889,
        0.33471,
        0.30897,
        0.28472,
        0.26287,
        0.24043,
        0.21955,
        0.19884,
        0.17902,
        0.15896,
        0.14222,
        0.12569,
        0.11078,
        0.096744,
        0.083955,
        0.071901,
        0.061399,
        0.051575,
        0.043986,
        0.036873,
        0.031042,
        0.026497,
        0.022253,
        0.019259,
        0.016458,
        0.014621,
        0.01305,
        0.011997,
        0.011171,
        0.010505,
        0.010195,
        0.0098217,
        0.00981,
        0.0094065,
        0.009293,
        0.0092355,
        0.009114,
        0.0090332,
        0.0089828,
        0.008935,
        0.0088745,
        0.008789,
        0.0088073,
        0.0087009,
        0.0086413,
        0.0085368,
        0.0085075,
        0.0084537,
        0.0084732,
        0.0083848,
        0.0082986,
        0.008268,
        0.0081769,
        0.0082267,
        0.008086,
        0.0080325,
        0.0078916,
        0.0077547,
        0.0077245,
        0.007652,
        0.0075915,
    ]
    return np.asarray(ex), np.asarray(ey)


class TestPDDMetric(TestCase):
    def test_normal(self):
        x, y = create_pdd_x_y()
        profile = FWXMProfile(values=y, x_values=x)
        pdd50 = profile.compute(metrics=[PDD(depth_mm=50)])
        self.assertAlmostEqual(pdd50, 90.16, delta=0.03)

    def test_passing_in_reverse_order_wont_change(self):
        # if the x and y data is passed in reverse order, the result should be the same.
        # i.e. x=[300, 299, 298, ...] should be the same as x=[0, 1, 2, ...]
        x, y = create_pdd_x_y()
        profile = FWXMProfile(values=y[::-1], x_values=x[::-1])
        pdd50 = profile.compute(metrics=[PDD(depth_mm=50)])
        self.assertAlmostEqual(pdd50, 90.16, delta=0.03)

    def test_pdd_normalized_to_max(self):
        x, y = create_pdd_x_y()
        profile = FWXMProfile(values=y, x_values=x)
        dmax = profile.compute(metrics=[PDD(depth_mm=50, normalize_to="max")])
        self.assertAlmostEqual(dmax, 89.96, delta=0.01)

    def test_cant_extrapolate(self):
        # extrapolation is not allowed
        x, y = create_pdd_x_y()
        profile = FWXMProfile(values=y, x_values=x)
        with self.assertRaises(ValueError):
            profile.compute(metrics=[PDD(depth_mm=500)])
        with self.assertRaises(ValueError):
            profile.compute(metrics=[PDD(depth_mm=-500)])

    def test_edge_of_profile_isnt_allowed(self):
        x, y = create_pdd_x_y()
        profile = FWXMProfile(values=y, x_values=x)
        with self.assertRaises(ValueError):
            profile.compute(metrics=[PDD(depth_mm=0)])

    def test_dmax_parameters_available(self):
        x, y = create_pdd_x_y()
        profile = FWXMProfile(values=y, x_values=x)
        pdd = profile.compute(
            metrics=[
                PDD(
                    depth_mm=50,
                    normalize_to="fit",
                    dmax_window_mm=10,
                    dmax_poly_order=3,
                )
            ]
        )
        self.assertAlmostEqual(pdd, 90.22, delta=0.01)

    def test_dmax_max(self):
        x, y = create_pdd_x_y()
        profile = FWXMProfile(values=y, x_values=x)
        pdd = profile.compute(metrics=[PDD(depth_mm=50, normalize_to="max")])
        self.assertAlmostEqual(pdd, 89.96, delta=0.01)


class TestSlopeMetric(TestCase):
    def test_normal(self):
        x, y = create_pdd_x_y()
        profile = FWXMProfilePhysical(values=y, x_values=x, dpmm=1)
        slope = profile.compute(metrics=[SlopeMetric()])
        self.assertAlmostEqual(slope, 0.00136, delta=0.001)

    def test_passing_inverted_ratios_fails(self):
        x, y = create_pdd_x_y()
        profile = FWXMProfilePhysical(values=y, x_values=x, dpmm=-1)
        with self.assertRaises(ValueError):
            profile.compute(metrics=[SlopeMetric(ratio_edges=(0.8, 0.2))])

    def test_passing_three_ratios_fails(self):
        x, y = create_pdd_x_y()
        profile = FWXMProfilePhysical(values=y, x_values=x, dpmm=1)
        with self.assertRaises(ValueError):
            profile.compute(metrics=[SlopeMetric(ratio_edges=(0.8, 0.2, 0.1))])


class TestDmaxMetric(TestCase):
    def test_normal(self):
        x, y = create_pdd_x_y()
        profile = FWXMProfile(values=y, x_values=x)
        dmax = profile.compute(metrics=[Dmax()])
        self.assertAlmostEqual(dmax, 22.04, delta=0.01)

    def test_passing_in_reverse_order_wont_change(self):
        # if the x and y data is passed in reverse order, the result should be the same.
        # i.e. x=[300, 299, 298, ...] should be the same as x=[0, 1, 2, ...]
        x, y = create_pdd_x_y()
        profile = FWXMProfile(values=y[::-1], x_values=x[::-1])
        dmax = profile.compute(metrics=[Dmax()])
        self.assertAlmostEqual(dmax, 22.04, delta=0.01)

    def test_large_window(self):
        x, y = create_pdd_x_y()
        profile = FWXMProfile(values=y, x_values=x)
        dmax = profile.compute(metrics=[Dmax(window_mm=50)])
        self.assertAlmostEqual(dmax, 21.3, delta=0.01)

    def test_poly_order(self):
        x, y = create_pdd_x_y()
        profile = FWXMProfile(values=y, x_values=x)
        dmax = profile.compute(metrics=[Dmax(poly_order=3)])
        self.assertAlmostEqual(dmax, 21.85, delta=0.01)

    def test_electron_pdd(self):
        x, y = create_electron_pdd_x_y()
        profile = FWXMProfile(values=y, x_values=x)
        dmax = profile.compute(metrics=[Dmax()])
        self.assertAlmostEqual(dmax, 20.92, delta=0.01)


class TestProfilePlugins(TestCase):
    def test_plot_without_metric_is_fine(self):
        array = generate_profile()
        profile = FWXMProfile(array, fwxm_height=50)
        profile.plot()

    def test_analyze_method(self):
        # tests the .analyze method, not the plugin itself
        array = generate_profile()
        profile = FWXMProfile(array, fwxm_height=50)
        profile.compute(metrics=[SymmetryPointDifferenceMetric()])
        self.assertIsInstance(profile.metric_values, dict)
        self.assertEqual(profile.metric_values["Point Difference Symmetry"], 0)

    def test_symmetry_point_difference_perfect(self):
        array = generate_profile()
        profile = FWXMProfile(array)
        profile.compute(metrics=[SymmetryPointDifferenceMetric()])
        self.assertEqual(profile.metric_values["Point Difference Symmetry"], 0)

    def test_symmetry_point_difference_right_negative(self):
        """When the profile skews higher on the right, the symmetry should be negative"""
        array = generate_profile(center=5)
        profile = FWXMProfile(array)
        profile.compute(metrics=[SymmetryPointDifferenceMetric()])
        self.assertAlmostEqual(
            profile.metric_values["Point Difference Symmetry"], -0.85, delta=0.01
        )

    def test_symmetry_point_difference_left_positive(self):
        """When the profile skews higher on the left, the symmetry should be positive"""
        array = generate_profile(center=-5)
        profile = FWXMProfile(array)
        profile.compute(metrics=[SymmetryPointDifferenceMetric()])
        self.assertAlmostEqual(
            profile.metric_values["Point Difference Symmetry"], 0.85, delta=0.01
        )

    def test_top_distance_perfect(self):
        """A perfect profile should have the top position at 0 for FFF"""
        array = generate_profile(field=FilterFreeFieldLayer)
        profile = FWXMProfilePhysical(array, dpmm=1)
        profile.compute(metrics=[TopDistanceMetric()])
        self.assertEqual(profile.metric_values["Top Distance"], 0)

    def test_top_distance_left(self):
        array = generate_profile(field=FilterFreeFieldLayer, center=5)
        profile = FWXMProfilePhysical(array, dpmm=1)
        profile.compute(metrics=[TopDistanceMetric()])
        self.assertAlmostEqual(profile.metric_values["Top Distance"], -18.8, delta=0.1)

    def test_top_distance_right(self):
        """A perfect profile should have the top position at 0 for FFF"""
        array = generate_profile(field=FilterFreeFieldLayer, center=-5)
        profile = FWXMProfilePhysical(array, dpmm=1)
        profile.compute(metrics=[TopDistanceMetric()])
        self.assertAlmostEqual(profile.metric_values["Top Distance"], 18.8, delta=0.1)

    def test_symmetry_quotient_perfect(self):
        """A perfectly symmetric profile should have a symmetry quotient of 100"""
        array = generate_profile()
        profile = FWXMProfilePhysical(array, dpmm=1)
        profile.compute(metrics=[SymmetryPointDifferenceQuotientMetric()])
        self.assertEqual(
            profile.metric_values["Point Difference Quotient Symmetry"], 100
        )

    def test_symmetry_quotient_offset(self):
        """The quotient will always be 100 or above"""
        array = generate_profile(center=5)
        profile = FWXMProfilePhysical(array, dpmm=1)
        profile.compute(metrics=[SymmetryPointDifferenceQuotientMetric()])
        self.assertAlmostEqual(
            profile.metric_values["Point Difference Quotient Symmetry"],
            100.84,
            delta=0.01,
        )

    def test_flatness_ratio_perfect(self):
        """A perfectly flat profile should have a flatness ratio of 1"""
        array = generate_profile(field=PerfectFieldLayer)
        profile = FWXMProfile(array)
        profile.compute(metrics=[FlatnessRatioMetric()])
        self.assertEqual(profile.metric_values["Flatness (Ratio)"], 100)

    def test_flatness_ratio_normal(self):
        array = generate_profile()
        profile = FWXMProfile(array)
        profile.compute(metrics=[FlatnessRatioMetric()])
        self.assertAlmostEqual(
            profile.metric_values["Flatness (Ratio)"], 103.02, delta=0.01
        )

    def test_flatness_difference_perfect(self):
        """A perfectly flat profile should have a flatness ratio of 1"""
        array = generate_profile(field=PerfectFieldLayer)
        profile = FWXMProfile(array)
        profile.compute(metrics=[FlatnessDifferenceMetric()])
        self.assertEqual(profile.metric_values["Flatness (Difference)"], 0)

    def test_flatness_difference_normal(self):
        """A perfectly flat profile should have a flatness ratio of 1"""
        array = generate_profile()
        profile = FWXMProfile(array)
        profile.compute(metrics=[FlatnessDifferenceMetric()])
        self.assertAlmostEqual(
            profile.metric_values["Flatness (Difference)"], 1.49, delta=0.01
        )

    def test_symmetry_area_perfect(self):
        """A perfectly symmetric profile should have a symmetry area of 0"""
        array = generate_profile(field=PerfectFieldLayer)
        profile = FWXMProfile(array)
        profile.compute(metrics=[SymmetryAreaMetric()])
        self.assertEqual(profile.metric_values["Symmetry (Area)"], 0)

    def test_symmetry_area_right_higher(self):
        array = generate_profile(center=5)
        profile = FWXMProfile(array)
        profile.compute(metrics=[SymmetryAreaMetric()])
        self.assertAlmostEqual(
            profile.metric_values["Symmetry (Area)"], -0.24, delta=0.01
        )

    def test_symmetry_area_left_higher(self):
        array = generate_profile(center=-5)
        profile = FWXMProfile(array)
        profile.compute(metrics=[SymmetryAreaMetric()])
        self.assertAlmostEqual(
            profile.metric_values["Symmetry (Area)"], 0.24, delta=0.01
        )

    def test_penumbra_left(self):
        array = generate_profile(field=PerfectFieldLayer)
        profile = FWXMProfilePhysical(array, dpmm=1)
        profile.compute(metrics=[PenumbraLeftMetric()])
        self.assertAlmostEqual(profile.metric_values["Left Penumbra"], 8.63, delta=0.01)

    def test_penumbra_right(self):
        array = generate_profile(field=PerfectFieldLayer)
        profile = FWXMProfilePhysical(array, dpmm=1)
        profile.compute(metrics=[PenumbraRightMetric()])
        self.assertAlmostEqual(
            profile.metric_values["Right Penumbra"], 8.63, delta=0.01
        )

    def test_penumbra_is_based_on_field_height(self):
        """The penumbra should be based on the field height, not the profile height"""
        array = generate_profile(field=PerfectFieldLayer)
        profile = FWXMProfilePhysical(array, dpmm=1, fwxm_height=30)
        profile.compute(metrics=[PenumbraLeftMetric()])
        self.assertAlmostEqual(profile.metric_values["Left Penumbra"], 5.78, delta=0.01)


class SingleProfileTests(TestCase):
    def test_normalization_max(self):
        """changed default parameter value to None in 3.10. 'max' should still work"""
        rng = np.random.default_rng()
        array = rng.random(100)

        # don't apply normalization initially, do it later
        p = SingleProfile(
            array,
            normalization_method=Normalization.NONE,
            interpolation=Interpolation.NONE,
            ground=False,
        )
        p.normalize("max")
        self.assertEqual(p.values.max(), 1)

    def test_normalization(self):
        rng = np.random.default_rng()
        array = rng.random(100)

        # don't apply normalization
        max_v = array.max()
        p = SingleProfile(
            array,
            normalization_method=Normalization.NONE,
            interpolation=Interpolation.NONE,
            ground=False,
        )
        self.assertEqual(max_v, p.values.max())

        # apply max norm
        p = SingleProfile(
            array,
            normalization_method=Normalization.MAX,
            interpolation=Interpolation.NONE,
        )
        self.assertEqual(1.0, p.values.max())

        # passing parameter as str
        p = SingleProfile(
            array, normalization_method="Max", interpolation=Interpolation.NONE
        )
        self.assertEqual(1.0, p.values.max())

        # make sure interpolation doesn't affect the norm
        p = SingleProfile(
            array,
            normalization_method=Normalization.MAX,
            interpolation=Interpolation.LINEAR,
        )
        self.assertEqual(1.0, p.values.max())

        # test out a real field
        field = generate_open_field()
        p = SingleProfile(field.image[:, 500], normalization_method=Normalization.MAX)
        self.assertEqual(1.0, p.values.max())

        # filtered beam center is less than max value
        p = SingleProfile(
            field.image[:, 500], normalization_method=Normalization.BEAM_CENTER
        )
        self.assertGreaterEqual(p.values.max(), 1.0)

    def test_x_values_are_monotonically_increasing(self):
        rng = np.random.default_rng()
        array = rng.random(100)
        with self.assertRaises(ValueError):
            SingleProfile(array, x_values=array, interpolation=Interpolation.NONE)

    def test_beam_center(self):
        # centered field
        field = generate_open_field()
        p = SingleProfile(
            field.image[:, int(field.shape[1] / 2)], interpolation=Interpolation.NONE
        )
        self.assertAlmostEqual(
            p.beam_center()["index (exact)"], field.shape[0] / 2, delta=1
        )

        # offset field
        field = generate_open_field(center=(10, 10))
        p = SingleProfile(
            field.image[:, int(field.shape[1] / 2)], interpolation=Interpolation.NONE
        )
        self.assertAlmostEqual(p.beam_center()["index (exact)"], 422, delta=1)

    def test_field_values_length(self):
        field = generate_open_field()
        p = SingleProfile(
            field.image[:, int(field.shape[1] / 2)], interpolation=Interpolation.NONE
        )
        field_data = p.field_data()
        width = field_data["right index (rounded)"] - (
            field_data["left index (rounded)"] - 1
        )  # we subtract one because the values include the left index element
        self.assertEqual(len(field_data["field values"]), width)

    def test_cax_odd_sized_array(self):
        arr = np.array(
            [0, 1, 1, 1, 2, 3, 5, 7, 9, 10, 10, 10, 9, 7, 5, 3, 2, 1, 1, 1, 0]
        )
        p = SingleProfile(
            arr,
            interpolation=Interpolation.NONE,
            normalization_method=Normalization.NONE,
        )
        field_data = p.field_data()
        self.assertEqual(field_data["cax index (rounded)"], 10)
        self.assertEqual(field_data["cax index (exact)"], 10)

    def test_cax_even_sized_array(self):
        arr = np.array([0, 1, 1, 1, 2, 3, 5, 7, 9, 10, 10, 9, 7, 5, 3, 2, 1, 1, 1, 0])
        p = SingleProfile(
            arr,
            interpolation=Interpolation.NONE,
            normalization_method=Normalization.NONE,
        )
        field_data = p.field_data()
        self.assertEqual(field_data["cax index (rounded)"], 10)
        self.assertEqual(field_data["cax index (exact)"], 9.5)

    def test_field_data_is_symmetric_odd(self):
        """For a symmetric array, the field values should be centered and symmetric."""
        arr = np.array(
            [0, 1, 1, 1, 2, 3, 5, 7, 9, 10, 10, 10, 9, 7, 5, 3, 2, 1, 1, 1, 0]
        )
        p = SingleProfile(
            arr,
            interpolation=Interpolation.NONE,
            normalization_method=Normalization.NONE,
        )
        field_data = p.field_data()
        self.assertEqual(field_data["field values"][0], field_data["field values"][-1])

    def test_field_data_is_symmetric_even(self):
        """For a symmetric array, the field values should be centered and symmetric."""
        arr = np.array([0, 1, 1, 1, 2, 3, 5, 7, 9, 10, 10, 9, 7, 5, 3, 2, 1, 1, 1, 0])
        p = SingleProfile(
            arr,
            interpolation=Interpolation.NONE,
            normalization_method=Normalization.NONE,
        )
        field_data = p.field_data()
        self.assertEqual(field_data["field values"][0], field_data["field values"][-1])

    def test_geometric_center(self):
        # centered field
        field = generate_open_field()
        p = SingleProfile(
            field.image[:, int(field.shape[1] / 2)], interpolation=Interpolation.NONE
        )
        self.assertAlmostEqual(
            p.geometric_center()["index (exact)"], field.shape[0] / 2, delta=1
        )

        # offset field should still be centered
        field = generate_open_field(center=(20, 20))
        p = SingleProfile(
            field.image[:, int(field.shape[1] / 2)], interpolation=Interpolation.NONE
        )
        self.assertAlmostEqual(
            p.geometric_center()["index (exact)"], field.shape[0] / 2, delta=1
        )

    def test_interpolation(self):
        # centered field
        field = generate_open_field()
        # no interp
        p = SingleProfile(
            field.image[:, int(field.shape[1] / 2)], interpolation=Interpolation.NONE
        )
        self.assertEqual(len(p.values), len(field.image[:, int(field.shape[1] / 2)]))

        # interp as str
        p = SingleProfile(field.image[:, int(field.shape[1] / 2)], interpolation=None)
        self.assertEqual(len(p.values), len(field.image[:, int(field.shape[1] / 2)]))

        # linear interp
        p = SingleProfile(
            field.image[:, int(field.shape[1] / 2)],
            interpolation=Interpolation.LINEAR,
            interpolation_factor=10,
        )
        self.assertEqual(
            len(p.values), len(field.image[:, int(field.shape[1] / 2)]) * 10
        )

        p = SingleProfile(
            field.image[:, int(field.shape[1] / 2)],
            interpolation=Interpolation.LINEAR,
            dpmm=1 / field.pixel_size,
            interpolation_resolution_mm=0.1,
        )
        # right length
        self.assertEqual(len(p.values), 3000)
        # dpmm is still the same
        self.assertEqual(p.dpmm, 1 / field.pixel_size)

        # spline interp
        p = SingleProfile(
            field.image[:, int(field.shape[1] / 2)],
            interpolation=Interpolation.SPLINE,
            interpolation_factor=10,
        )
        self.assertEqual(
            len(p.values), len(field.image[:, int(field.shape[1] / 2)]) * 10
        )

        p = SingleProfile(
            field.image[:, int(field.shape[1] / 2)],
            interpolation=Interpolation.SPLINE,
            dpmm=1 / field.pixel_size,
            interpolation_resolution_mm=0.1,
        )
        # right length
        self.assertEqual(len(p.values), 3000)
        # dpmm is still the same
        self.assertEqual(p.dpmm, 1 / field.pixel_size)

    def test_resample(self):
        # test resampling doubles the length of the values
        field = generate_open_field()
        profile = SingleProfile(
            field.image[:, int(field.shape[1] / 2)],
            interpolation=Interpolation.LINEAR,
            dpmm=1 / field.pixel_size,
            interpolation_resolution_mm=0.1,
        )
        resampled_profile = profile.resample(interpolation_resolution_mm=0.05)
        self.assertEqual(len(profile.values), len(resampled_profile.values) / 2)
        self.assertIsInstance(resampled_profile, SingleProfile)
        original_data = profile.fwxm_data()
        resampled_data = resampled_profile.fwxm_data()
        self.assertAlmostEqual(
            original_data["width (exact)"], resampled_data["width (exact)"], delta=0.001
        )

        # test resampling via factor
        field = generate_open_field()
        profile = SingleProfile(
            field.image[:, int(field.shape[1] / 2)],
            interpolation=Interpolation.LINEAR,
            dpmm=None,
            interpolation_factor=10,
        )
        resampled_profile = profile.resample(
            interpolation_factor=5
        )  # 5x resolution as first profile
        self.assertEqual(len(profile), len(resampled_profile) / 5)
        self.assertIsInstance(resampled_profile, SingleProfile)
        original_data = profile.fwxm_data()
        resampled_data = resampled_profile.fwxm_data()
        self.assertAlmostEqual(
            original_data["width (exact)"], resampled_data["width (exact)"], delta=0.001
        )

    def test_gamma_perfect(self):
        # centered field
        field = generate_open_field()
        reference = evaluation = SingleProfile(
            field.image[:, int(field.shape[1] / 2)],
            interpolation=Interpolation.LINEAR,
            dpmm=1 / field.pixel_size,
            interpolation_resolution_mm=0.1,
        )
        gamma = reference.gamma(evaluation)
        self.assertIsInstance(gamma, np.ndarray)
        self.assertTrue(np.isnan(gamma[0]))  # first element is under threshold
        self.assertAlmostEqual(
            np.nanmean(gamma), 0, delta=0.0003
        )  # it's 0 because it's the same profile +/- presicion error

    def test_gamma_perfect_downsampled(self):
        # centered field
        field = generate_open_field()
        reference = SingleProfile(
            field.image[:, int(field.shape[1] / 2)],
            interpolation=Interpolation.LINEAR,
            dpmm=1 / field.pixel_size,
            interpolation_resolution_mm=0.1,
        )
        # downsample eval profile. We will upsample inside gamma and it should be pretty close to the original
        evaluation = reference.resample(interpolation_resolution_mm=0.2)
        gamma = reference.gamma(evaluation)
        self.assertIsInstance(gamma, np.ndarray)
        self.assertTrue(np.isnan(gamma[0]))  # first element is under threshold
        self.assertLess(
            np.nanmean(gamma), 0.002
        )  # there are small errors introduced by down/up sampled but they're very small

    def test_gamma_perfect_upsampled(self):
        # centered field
        field = generate_open_field()
        reference = SingleProfile(
            field.image[:, int(field.shape[1] / 2)],
            interpolation=Interpolation.LINEAR,
            dpmm=1 / field.pixel_size,
            interpolation_resolution_mm=0.1,
        )
        # upsample eval profile. We will downsample inside gamma and it should be pretty close to the original
        evaluation = reference.resample(interpolation_resolution_mm=0.05)
        gamma = reference.gamma(evaluation)
        self.assertIsInstance(gamma, np.ndarray)
        self.assertTrue(np.isnan(gamma[0]))  # first element is under threshold
        self.assertLess(
            np.nanmean(gamma), 0.002
        )  # there are small errors introduced by down/up sampled but they're very small

    def test_gamma_unequal_samples(self):
        # centered field
        field = generate_open_field()
        reference = SingleProfile(
            field.image[:, int(field.shape[1] / 2)],
            interpolation=Interpolation.LINEAR,
            dpmm=1 / field.pixel_size,
            interpolation_resolution_mm=0.1,
        )
        evaluation = SingleProfile(
            field.image[:-50, int(field.shape[1] / 2)],
            interpolation=Interpolation.LINEAR,
            dpmm=1 / field.pixel_size,
            interpolation_resolution_mm=0.1,
        )
        gamma = reference.gamma(evaluation)
        self.assertIsInstance(gamma, np.ndarray)
        self.assertEqual(len(gamma), len(reference.values))
        self.assertTrue(np.isnan(gamma[0]))  # first element is under threshold
        self.assertAlmostEqual(
            np.nanmean(gamma), 0, delta=0.001
        )  # it's 0 because it's the same profile +/- presicion error

    def test_gamma_no_dpmm(self):
        field = generate_open_field()
        reference = evaluation = SingleProfile(
            field.image[:, int(field.shape[1] / 2)], dpmm=None
        )
        with self.assertRaises(ValueError):
            reference.gamma(evaluation)


class MultiProfileTestMixin:
    values = np.array
    peak_max_idxs = (0,)
    valley_max_idxs = (0,)
    peak_fwxm_idxs = (0,)
    subdivide_fwxm_centers = (0,)

    @classmethod
    def setUpClass(cls):
        cls.profile = MultiProfile(cls.values)

    def test_find_peaks(self):
        peaks, _ = self.profile.find_peaks()
        for peak, known_peak in zip(peaks, self.peak_max_idxs):
            self.assertAlmostEqual(peak, known_peak, delta=1)

    def test_find_fwxm_peaks(self):
        peak_idxs, _ = self.profile.find_fwxm_peaks()
        for peak, known_peak in zip(peak_idxs, self.peak_fwxm_idxs):
            self.assertAlmostEqual(peak, known_peak, delta=1)

    def test_find_valleys(self):
        valleys, _ = self.profile.find_valleys()
        for valley, known_valley in zip(valleys, self.valley_max_idxs):
            self.assertAlmostEqual(valley, known_valley, delta=1)


class MultiProfileTriangle(MultiProfileTestMixin, TestCase):
    x_values = np.linspace(0, 8 * np.pi, num=200)
    values = sps.sawtooth(x_values, width=0.5)
    valley_max_idxs = (50, 100, 150)
    peak_max_idxs = (25, 75, 125, 175)
    peak_fwxm_idxs = (25, 75, 125, 175)

    def test_ground_profile(self):
        """Test that the profile is properly grounded to 0."""
        p = MultiProfile(self.values)
        # the minimum shouldn't be zero to start with
        self.assertFalse(p.values.min() == 0)
        # but it should be after grounding
        p.ground()
        self.assertTrue(p.values.min() == 0)


class CircleProfileTestMixin:
    klass = CircleProfile
    image_file_location = get_file_from_cloud_test_repo(["Starshot", "Starshot-1.tif"])
    radius = 300
    peak_idxs = (0,)
    valley_idxs = (0,)
    fwxm_peak_idxs = (0,)
    center_point = (507, 650)

    @classmethod
    def setUpClass(cls):
        img = image.load(cls.image_file_location)
        cls.profile = cls.klass(cls.center_point, cls.radius, img.array)
        cls.profile.filter(size=0.01, kind="gaussian")

    def test_locations(self):
        first_x_location = self.profile.radius + self.profile.center.x
        self.assertAlmostEqual(first_x_location, self.profile.x_locations[0], delta=1)

    def test_peak_idxs(self):
        for known, meas in zip(self.peak_idxs, self.profile.find_peaks()[0]):
            self.assertAlmostEqual(known, meas, delta=1)

    def test_valley_idxs(self):
        for known, meas in zip(
            self.valley_idxs, self.profile.find_valleys(min_distance=0.08)[0]
        ):
            self.assertAlmostEqual(known, meas, delta=1)

    def test_fwxm_peak_idxs(self):
        for known, meas in zip(self.fwxm_peak_idxs, self.profile.find_fwxm_peaks()[0]):
            self.assertAlmostEqual(known, meas, delta=1)

    def test_add_to_axes(self):
        # shouldn't raise
        self.profile.plot2axes()


class CircleProfileStarshot(CircleProfileTestMixin, TestCase):
    peak_idxs = [219, 480, 738, 984, 1209, 1421, 1633, 1864]
    valley_idxs = [95, 348, 607, 860, 1098, 1316, 1527, 1743]
    fwxm_peak_idxs = [218, 480, 738, 984, 1209, 1421, 1633, 1864]


class CollapsedCircleProfileStarshot(CircleProfileTestMixin, TestCase):
    klass = CollapsedCircleProfile
    peak_idxs = [241, 529, 812, 1084, 1331, 1563, 1797, 2051]
    valley_idxs = [104, 397, 667, 946, 1210, 1451, 1680, 1916]
    fwxm_peak_idxs = [241, 529, 812, 1084, 1331, 1563, 1797, 2052]


class TestStretchDeprecation(TestCase):
    def test_stretch_deprecation(self):
        arr = np.array([0, 2, 5, 10], dtype=np.uint8)
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger a warning.
            stretch(arr)
            # Verify some things
            assert len(w) >= 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert "Using stretch from the profile module is deprecated" in str(
                w[0].message
            )
