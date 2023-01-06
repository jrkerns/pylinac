from unittest import TestCase

import numpy as np
import scipy.signal as sps

from pylinac.core import image
from pylinac.core.image_generator.simulators import Simulator
from pylinac.core.profile import (
    SingleProfile,
    MultiProfile,
    CircleProfile,
    CollapsedCircleProfile,
    Normalization,
    Interpolation,
    gamma_1d,
)
from tests_basic.utils import get_file_from_cloud_test_repo


def generate_open_field(field_size=(100, 100), sigma=2, center=(0, 0)) -> Simulator:
    from pylinac.core.image_generator import AS1000Image
    from pylinac.core.image_generator.layers import (
        FilteredFieldLayer,
        GaussianFilterLayer,
    )

    as1000 = AS1000Image()  # this will set the pixel size and shape automatically
    as1000.add_layer(
        FilteredFieldLayer(field_size_mm=field_size, cax_offset_mm=center)
    )  # create a 50x50mm square field
    as1000.add_layer(
        GaussianFilterLayer(sigma_mm=sigma)
    )  # add an image-wide gaussian to simulate penumbra/scatter
    return as1000


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


class SingleProfileTests(TestCase):
    def test_normalization(self):
        array = np.random.rand(1, 100).squeeze()

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

    values = np.ndarray
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
