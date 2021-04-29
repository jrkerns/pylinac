from unittest import TestCase
import os.path as osp

import numpy as np
import scipy.signal as sps

from pylinac.core import image
from pylinac.core.image_generator.simulators import Simulator
from pylinac.core.profile import SingleProfile, MultiProfile, CircleProfile, CollapsedCircleProfile, Normalization, \
    Interpolation


def generate_open_field(field_size=(100, 100), sigma=2, center=(0, 0)) -> Simulator:
    from pylinac.core.image_generator import AS1000Image
    from pylinac.core.image_generator.layers import FilteredFieldLayer, GaussianFilterLayer

    as1000 = AS1000Image()  # this will set the pixel size and shape automatically
    as1000.add_layer(FilteredFieldLayer(field_size_mm=field_size, cax_offset_mm=center))  # create a 50x50mm square field
    as1000.add_layer(GaussianFilterLayer(sigma_mm=sigma))  # add an image-wide gaussian to simulate penumbra/scatter
    return as1000


class SingleProfileTests(TestCase):

    def test_normalization(self):
        array = np.random.rand(1, 100).squeeze()

        # don't apply normalization
        max_v = array.max()
        p = SingleProfile(array, normalization_method=Normalization.NONE, interpolation=Interpolation.NONE, ground=False)
        self.assertEqual(max_v, p.values.max())

        # apply max norm
        p = SingleProfile(array, normalization_method=Normalization.MAX, interpolation=Interpolation.NONE)
        self.assertEqual(1.0, p.values.max())

        # make sure interpolation doesn't affect the norm
        p = SingleProfile(array, normalization_method=Normalization.MAX, interpolation=Interpolation.LINEAR)
        self.assertEqual(1.0, p.values.max())

        # test out a real field
        field = generate_open_field()
        p = SingleProfile(field.image[:, 500], normalization_method=Normalization.MAX)
        self.assertEqual(1.0, p.values.max())

        # filtered beam center is less than max value
        p = SingleProfile(field.image[:, 500], normalization_method=Normalization.BEAM_CENTER)
        self.assertGreaterEqual(p.values.max(), 1.0)

    def test_beam_center(self):
        # centered field
        field = generate_open_field()
        p = SingleProfile(field.image[:, int(field.shape[1]/2)], interpolation=Interpolation.NONE)
        self.assertAlmostEqual(p.beam_center()['index (exact)'], field.shape[0]/2, delta=1)

        # offset field
        field = generate_open_field(center=(10, 10))
        p = SingleProfile(field.image[:, int(field.shape[1]/2)], interpolation=Interpolation.NONE)
        self.assertAlmostEqual(p.beam_center()['index (exact)'], 422, delta=1)

    def test_geometric_center(self):
        # centered field
        field = generate_open_field()
        p = SingleProfile(field.image[:, int(field.shape[1]/2)], interpolation=Interpolation.NONE)
        self.assertAlmostEqual(p.geometric_center()['index (exact)'], field.shape[0]/2, delta=1)

        # offset field should still be centered
        field = generate_open_field(center=(20, 20))
        p = SingleProfile(field.image[:, int(field.shape[1]/2)], interpolation=Interpolation.NONE)
        self.assertAlmostEqual(p.geometric_center()['index (exact)'], field.shape[0]/2, delta=1)

    def test_interpolation(self):
        # centered field
        field = generate_open_field()
        # no interp
        p = SingleProfile(field.image[:, int(field.shape[1] / 2)], interpolation=Interpolation.NONE)
        self.assertEqual(len(p.values), len(field.image[:, int(field.shape[1] / 2)]))

        # linear interp
        p = SingleProfile(field.image[:, int(field.shape[1] / 2)], interpolation=Interpolation.LINEAR, interpolation_factor=10)
        self.assertEqual(len(p.values), len(field.image[:, int(field.shape[1] / 2)])*10)

        p = SingleProfile(field.image[:, int(field.shape[1] / 2)], interpolation=Interpolation.LINEAR,
                          dpmm=1/field.pixel_size, interpolation_resolution_mm=0.1)
        # right length
        self.assertEqual(len(p.values), len(field.image[:, int(field.shape[1] / 2)])*field.pixel_size/0.1)
        # right dpmm
        self.assertEqual(p.dpmm, 10)

        # spline interp
        p = SingleProfile(field.image[:, int(field.shape[1] / 2)], interpolation=Interpolation.SPLINE, interpolation_factor=10)
        self.assertEqual(len(p.values), len(field.image[:, int(field.shape[1] / 2)])*10)

        p = SingleProfile(field.image[:, int(field.shape[1] / 2)], interpolation=Interpolation.SPLINE,
                          dpmm=1/field.pixel_size, interpolation_resolution_mm=0.1)
        # right length
        self.assertEqual(len(p.values), len(field.image[:, int(field.shape[1] / 2)])*field.pixel_size/0.1)
        # right dpmm
        self.assertEqual(p.dpmm, 10)




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

    x_values = np.linspace(0, 8*np.pi, num=200)
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
    image_file_location = osp.join(osp.dirname(osp.dirname(osp.abspath(__file__))), 'test_files', 'Starshot',
                                   'Starshot#1.tif')
    radius = 300
    peak_idxs = (0,)
    valley_idxs = (0,)
    fwxm_peak_idxs = (0,)
    center_point = (507, 650)

    @classmethod
    def setUpClass(cls):
        img = image.load(cls.image_file_location)
        cls.profile = cls.klass(cls.center_point, cls.radius, img.array)
        cls.profile.filter(size=0.01, kind='gaussian')

    def test_locations(self):
        first_x_location = self.profile.radius + self.profile.center.x
        self.assertAlmostEqual(first_x_location, self.profile.x_locations[0], delta=1)

    def test_peak_idxs(self):
        for known, meas in zip(self.peak_idxs, self.profile.find_peaks()[0]):
            self.assertAlmostEqual(known, meas, delta=1)

    def test_valley_idxs(self):
        for known, meas in zip(self.valley_idxs, self.profile.find_valleys(min_distance=0.08)[0]):
            self.assertAlmostEqual(known, meas, delta=1)

    def test_fwxm_peak_idxs(self):
        for known, meas in zip(self.fwxm_peak_idxs, self.profile.find_fwxm_peaks()[0]):
            self.assertAlmostEqual(known, meas, delta=1)

    def test_add_to_axes(self):
        # shouldn't raise
        self.profile.plot2axes()


class CircleProfileStarshot(CircleProfileTestMixin, TestCase):

    peak_idxs = [219,  480,  738,  984, 1209, 1421, 1633, 1864]
    valley_idxs = [95,  348,  607,  860, 1098, 1316, 1527, 1743]
    fwxm_peak_idxs = [218,  480,  738,  984, 1209, 1421, 1633, 1864]


class CollapsedCircleProfileStarshot(CircleProfileTestMixin, TestCase):

    klass = CollapsedCircleProfile
    peak_idxs = [241,  529,  812, 1084, 1331, 1563, 1797, 2051]
    valley_idxs = [104,  397,  667,  946, 1210, 1451, 1680, 1916]
    fwxm_peak_idxs = [241,  529,  812, 1084, 1331, 1563, 1797, 2052]
