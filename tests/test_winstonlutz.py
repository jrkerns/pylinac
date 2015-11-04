"""Tests for the Winston-Lutz module."""
import os.path as osp
from unittest import TestCase

from pylinac import WinstonLutz
from pylinac.core.geometry import Vector, vector_is_close


image_bank_dir = osp.abspath(osp.join('..', '..', 'unorganized linac data', 'Winston-Lutz', 'Chicago', 'WL-Final_C&G&C_Final'))


class GeneralTests(TestCase):
    image_dir = osp.abspath(osp.join('..', '..', 'unorganized linac data', 'Winston-Lutz', 'Chicago', 'WL-Final_C&G&C_Final'))

    @classmethod
    def setUpClass(cls):
        cls.wl = WinstonLutz(cls.image_dir)

    def test_plot(self):
        self.wl.plot_images()  # shouldn't raise
        self.wl.plot_gantry_sag()

    def test_results(self):
        print(self.wl.results())  # shouldn't raise


class WinstonLutzMixin:
    image_dir = ''
    num_images = 0
    gantry_iso_size = 0
    gantry_iso2bb_vector = Vector
    gantry_sag = 0
    collimator_iso_size = 0
    collimator_iso2bb_vector = Vector
    couch_iso_size = 0
    couch_iso2bb_vector = Vector
    cax2bb_max_distance = 0
    cax2bb_median_distance = 0
    variable_axes = {0: 'Gantry'}  # fill with as many {image: axis} pairs as desired

    @classmethod
    def setUpClass(cls):
        cls.wl = WinstonLutz(cls.image_dir)

    def test_number_of_images(self):
        self.assertEqual(len(self.wl.images), self.num_images)

    def test_variable_axes(self):
        for idx, axis in self.variable_axes.items():
            self.assertEqual(self.wl.images[idx].variable_axis, axis)

    def test_gantry_iso(self):
        # test iso size
        self.assertAlmostEqual(self.wl.gantry_iso_size, self.gantry_iso_size, delta=0.2)
        # test iso vector
        self.assertTrue(vector_is_close(self.wl.gantry_iso2bb_vector, self.gantry_iso2bb_vector))

    def test_gantry_sage(self):
        self.assertAlmostEqual(self.wl.gantry_sag(), self.gantry_sag, delta=0.15)

    def test_collimator_iso(self):
        # test iso size
        self.assertAlmostEqual(self.wl.collimator_iso_size, self.collimator_iso_size, delta=0.2)
        # test iso vector
        self.assertTrue(vector_is_close(self.wl.collimator_iso2bb_vector, self.collimator_iso2bb_vector))

    def test_couch_iso(self):
        # test iso size
        self.assertAlmostEqual(self.wl.couch_iso_size, self.couch_iso_size, delta=0.2)
        # test iso vector
        self.assertTrue(vector_is_close(self.wl.couch_iso2bb_vector, self.couch_iso2bb_vector))


class WLKiX(WinstonLutzMixin, TestCase):
    image_dir = image_bank_dir = osp.abspath(osp.join('..', '..', 'unorganized linac data', 'Winston-Lutz', 'Katy iX'))
    num_images = 17
    gantry_iso_size = 0.72
    gantry_iso2bb_vector = Vector(-0.3, 0.1, -0.2)
    gantry_sag = 1
    collimator_iso_size = 0.5
    collimator_iso2bb_vector = Vector(-0.6, -0.3)
    couch_iso_size = 0.7
    couch_iso2bb_vector = Vector(-0.2)
    cax2bb_max_distance = 1.07
    cax2bb_median_distance = 0.8
    variable_axes = {0: 'Reference'}


class WLChTB(WinstonLutzMixin, TestCase):
    image_dir = osp.abspath(osp.join('..', '..', 'unorganized linac data', 'Winston-Lutz', 'Chicago', 'WL-Final_C&G&C_Final'))
    num_images = 17
    gantry_iso_size = 0.37
    gantry_iso2bb_vector = Vector(-0.2)
    gantry_sag = 0.54
    collimator_iso_size = 0.1
    collimator_iso2bb_vector = Vector(-0.2, -0.3)
    couch_iso_size = 0.1
    couch_iso2bb_vector = Vector(-0.25, -0.2)
    cax2bb_max_distance = 0.5
    cax2bb_median_distance = 0.3
    variable_axes = {0: 'Collimator'}


