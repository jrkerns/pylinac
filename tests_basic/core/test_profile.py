from unittest import TestCase
import os.path as osp

import numpy as np
import scipy.signal as sps

from pylinac.core import image
from pylinac.core.profile import SingleProfile, MultiProfile, CircleProfile, CollapsedCircleProfile


class SingleProfileMixin:

    ydata = np.ndarray
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
        cls.profile = SingleProfile(cls.ydata)

    def test_fwxms(self):
        for fwxm, fwhm_idx in self.fwxm_indices.items():
            self.assertAlmostEqual(self.profile.fwxm(fwxm), fwhm_idx, delta=1)

    def test_fwxm_centers(self):
        # test values, interpolated and not interpolated
        for fwxm, fwhm_val in self.fwxm_center_values.items():
            self.assertAlmostEqual(self.profile.fwxm_center(fwxm)[1], fwhm_val, delta=0.1)
        for fwxm, fwhm_val in self.fwxm_center_values.items():
            self.assertAlmostEqual(self.profile.fwxm_center(fwxm, interpolate=True)[1], fwhm_val, delta=0.1)

        # test indices, interpolated and not interpolated
        for fwxm, fwhm_idx in self.fwxm_center_indices.items():
            self.assertAlmostEqual(self.profile.fwxm_center(fwxm)[0], fwhm_idx, delta=1)
        for fwxm, fwhm_idx in self.fwxm_center_indices.items():
            self.assertAlmostEqual(self.profile.fwxm_center(fwxm, interpolate=True)[0], fwhm_idx, delta=1)

    def test_penum_widths(self):
        # test 80/20, interp and non-interp
        lt_penum, rt_penum = self.profile.penumbra_width(lower=20, upper=80)
        self.assertAlmostEqual(lt_penum, self.penumbra_widths_8020['left'], delta=1)
        self.assertAlmostEqual(rt_penum, self.penumbra_widths_8020['right'], delta=1)
        self.assertAlmostEqual(np.mean([lt_penum, rt_penum]), self.penumbra_widths_8020['both'], delta=1)
        # test 90/10
        lt_penum, rt_penum = self.profile.penumbra_width(lower=10, upper=90)
        self.assertAlmostEqual(lt_penum, self.penumbra_widths_9010['left'], delta=1)
        self.assertAlmostEqual(rt_penum, self.penumbra_widths_9010['right'], delta=1)
        self.assertAlmostEqual(np.mean([lt_penum, rt_penum]), self.penumbra_widths_9010['both'], delta=1)

    def test_field_value_length(self):
        field_values = self.profile.field_values()
        self.assertAlmostEqual(len(field_values), self.field_value_length, delta=2)

    def test_field_edges(self):
        for meas, known in zip(self.field_edge_indices, self.profile.field_edges()):
            self.assertAlmostEqual(meas, known, delta=0.1)

    def test_field_calculations(self):
        for calc, val in self.field_calculations.items():
            self.assertAlmostEqual(self.profile.field_calculation(calculation=calc), val, delta=0.1)

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
    fwxm_indices = {30: 115, 50: 82, 80: 33}
    fwxm_center_values = {40: 1, 60: 1, 80: 1}
    fwxm_center_indices = {40: 117, 60: 117, 80: 117}
    penumbra_widths_8020 = {'left': 49, 'right': 49, 'both': 49}
    penumbra_widths_9010 = {'left': 65, 'right': 65, 'both': 65}
    field_edge_indices = (84, 150)
    field_calculations = {'max': 0.99, 'mean': 0.64, 'min': 0.43}
    field_value_length = 66
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
