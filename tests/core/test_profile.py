from unittest import TestCase
import os.path as osp

import numpy as np
import scipy.signal as sps

from pylinac.core import image
from pylinac.core.profile import SingleProfile, MultiProfile, CircleProfile, CollapsedCircleProfile


class SingleProfileMixin:

    ydata = np.ndarray
    normalize_sides = True
    fwxm_indices = {30: 0, 50: 0, 80: 0}
    fwxm_center_values = {40: 0, 60: 0, 80: 0}
    fwxm_center_indices = {40: 0, 60: 0, 80: 0}
    field_value_length = 0
    penumbra_widths_8020 = {'left': 0, 'right': 0, 'both': 0}
    penumbra_widths_9010 = {'left': 0, 'right': 0, 'both': 0}
    field_edge_indices = (0, 0)
    field_calculations = {'max': 0, 'mean': 0, 'min': 0}
    peak_idx = 0

    @classmethod
    def setUpClass(cls):
        cls.profile = SingleProfile(cls.ydata, normalize_sides=cls.normalize_sides)

    def test_fwxms(self):
        for fwxm, fwhm_idx in self.fwxm_indices.items():
            self.assertAlmostEqual(self.profile.fwxm(fwxm), fwhm_idx, delta=1)
        for fwxm, fwhm_idx in self.fwxm_indices.items():
            self.assertAlmostEqual(self.profile.fwxm(fwxm, interpolate=True), fwhm_idx, delta=1)

    def test_fwxm_centers(self):
        # test indices, interpolated and not interpolated
        for fwxm, fwhm_val in self.fwxm_center_values.items():
            self.assertAlmostEqual(self.profile.fwxm_center(fwxm, kind='value'), fwhm_val, delta=0.1)
        for fwxm, fwhm_val in self.fwxm_center_values.items():
            self.assertAlmostEqual(self.profile.fwxm_center(fwxm, kind='value', interpolate=True), fwhm_val, delta=0.1)

        # test indices, interpolated and not interpolated
        for fwxm, fwhm_idx in self.fwxm_center_indices.items():
            self.assertAlmostEqual(self.profile.fwxm_center(fwxm), fwhm_idx, delta=1)
        for fwxm, fwhm_idx in self.fwxm_center_indices.items():
            self.assertAlmostEqual(self.profile.fwxm_center(fwxm, interpolate=True), fwhm_idx, delta=1)

    def test_penum_widths(self):
        # test 80/20, interp and non-interp
        for side, val in self.penumbra_widths_8020.items():
            self.assertAlmostEqual(self.profile.penumbra_width(side, lower=20, upper=80), val, delta=0.1)
        for side, val in self.penumbra_widths_8020.items():
            self.assertAlmostEqual(self.profile.penumbra_width(side, lower=20, upper=80, interpolate=True), val, delta=1)

        # test 90/10
        for side, val in self.penumbra_widths_9010.items():
            self.assertAlmostEqual(self.profile.penumbra_width(side, lower=10, upper=90), val, delta=0.1)
        for side, val in self.penumbra_widths_9010.items():
            self.assertAlmostEqual(self.profile.penumbra_width(side, lower=10, upper=90, interpolate=True), val, delta=1)

    def test_field_value_length(self):
        field_values = self.profile.field_values()
        self.assertAlmostEqual(len(field_values), self.field_value_length, delta=2)

    def test_field_edges(self):
        for meas, known in zip(self.field_edge_indices, self.profile.field_edges()):
            self.assertAlmostEqual(meas, known, delta=0.1)

    def test_field_calculations(self):
        for calc, val in self.field_calculations.items():
            self.assertAlmostEqual(self.profile.field_calculation(calculation=calc), val, delta=0.1)

    def test_initial_peak(self):
        detected_initial_peak_idx = self.profile._initial_peak_idx
        self.assertAlmostEqual(detected_initial_peak_idx, self.peak_idx, delta=1)

    def test_unnormalized_peaks(self):
        pass


class SingleProfileTriangle(SingleProfileMixin, TestCase):

    xdata = np.linspace(0, 2*np.pi, num=200)
    ydata = sps.sawtooth(xdata, width=0.5)
    fwxm_indices = {30: 140, 50: 101, 80: 41}
    fwxm_center_values = {40: 1, 60: 1, 80: 1}
    fwxm_center_indices = {40: 100, 60: 100, 80: 100}
    penumbra_widths_8020 = {'left': 60, 'right': 60, 'both': 60}
    penumbra_widths_9010 = {'left': 80, 'right': 80, 'both': 80}
    field_edge_indices = (60, 140)
    field_calculations = {'max': 0.99, 'mean': 0.60, 'min': 0.21}
    field_value_length = 80
    peak_idx = 100


class SingleProfileCutoffTriangle(SingleProfileMixin, TestCase):
    """A triangle cut short on the right side. Can effectively test the normalization of each side."""
    xdata = np.linspace(0, 1.7 * np.pi, num=200)
    ydata = sps.sawtooth(xdata, width=0.5)
    fwxm_indices = {30: 139, 50: 100, 80: 40}
    fwxm_center_values = {40: 0.83, 60: 0.88, 80: 0.95}
    fwxm_center_indices = {40: 107, 60: 110.5, 80: 114}
    penumbra_widths_8020 = {'left': 70, 'right': 49, 'both': 59.5}
    penumbra_widths_9010 = {'left': 94, 'right': 65, 'both': 79.5}
    field_edge_indices = (68, 148)
    field_calculations = {'max': 0.99, 'mean': 0.64, 'min': 0.18}
    field_value_length = 80
    peak_idx = 117


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
        peaks = self.profile.find_peaks()
        for peak, known_peak in zip(peaks, self.peak_max_idxs):
            self.assertAlmostEqual(peak, known_peak, delta=1)

    def test_find_fwxm_peaks(self):
        peakidxs = self.profile.find_fwxm_peaks()
        for peak, known_peak in zip(peakidxs, self.peak_fwxm_idxs):
            self.assertAlmostEqual(peak, known_peak, delta=1)

    def test_find_valleys(self):
        valleys = self.profile.find_valleys()
        for valley, known_valley in zip(valleys, self.valley_max_idxs):
            self.assertAlmostEqual(valley, known_valley, delta=1)

    def test_subdivide(self):
        self.profile.find_peaks()
        profiles = self.profile.subdivide()
        for profile, known_fwxm_center in zip(profiles, self.subdivide_fwxm_centers):
            fwxm_center = profile.fwxm_center()
            self.assertAlmostEqual(fwxm_center, known_fwxm_center, delta=1)


class MultiProfileTriangle(MultiProfileTestMixin, TestCase):

    x_values = np.linspace(0, 8*np.pi, num=200)
    values = sps.sawtooth(x_values, width=0.5)
    valley_max_idxs = (50, 100, 150)
    peak_max_idxs = (25, 75, 125, 175)
    peak_fwxm_idxs = (25, 75, 125, 175)
    subdivide_fwxm_centers = (25, 50, 50, 50)

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

    def test_locations(self):
        first_x_location = self.profile.radius + self.profile.center.x
        self.assertAlmostEqual(first_x_location, self.profile.x_locations[0], delta=1)

    def test_peak_idxs(self):
        for known, meas in zip(self.peak_idxs, self.profile.find_peaks()):
            self.assertAlmostEqual(known, meas, delta=1)

    def test_valley_idxs(self):
        for known, meas in zip(self.valley_idxs, self.profile.find_valleys()):
            self.assertAlmostEqual(known, meas, delta=1)

    def test_fwxm_peak_idxs(self):
        for known, meas in zip(self.fwxm_peak_idxs, self.profile.find_fwxm_peaks()):
            self.assertAlmostEqual(known, meas, delta=1)

    def test_add_to_axes(self):
        # shouldn't raise
        self.profile.plot2axes()


class CircleProfileStarshot(CircleProfileTestMixin, TestCase):

    peak_idxs = [218., 480., 738., 985., 1209., 1420., 1633., 1857.]
    valley_idxs = [118., 338., 606., 911., 1138., 1364., 1529., 1799.]
    fwxm_peak_idxs = [219.5, 479.5, 738.0, 984.5, 1209.0, 1421.0, 1633.5, 1857.5]


class CollapsedCircleProfileStarshot(CircleProfileTestMixin, TestCase):

    klass = CollapsedCircleProfile
    peak_idxs = [241., 529., 812., 1083., 1330., 1563., 1796., 2044.]
    valley_idxs = [100., 405., 673., 960., 1241., 1481., 1714., 1916.]
    fwxm_peak_idxs = [241.0, 529.5, 812.5, 1084.0, 1330.5, 1563.0, 1797.0, 2043.5]
