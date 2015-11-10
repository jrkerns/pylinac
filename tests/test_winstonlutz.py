"""Tests for the Winston-Lutz module."""
import os.path as osp
from unittest import TestCase

from pylinac import WinstonLutz
from pylinac.core.geometry import Vector, vector_is_close


class GeneralTests(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.wl = WinstonLutz.from_demo_images()

    def test_loading_demo_files(self):
        wl = WinstonLutz.from_demo_images()  # shouldn't raise

    def test_run_demo(self):
        WinstonLutz().run_demo()  # shouldn't raise

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

    def test_gantry_sag(self):
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


class WLDemo(WinstonLutzMixin, TestCase):
    num_images = 17
    gantry_sag = 1
    gantry_iso_size = 0.7
    gantry_iso2bb_vector = Vector(0.2, 0.1, 0.4)
    collimator_iso_size = 0.5
    collimator_iso2bb_vector = Vector(0.2, 0.4, 0)
    couch_iso_size = 1.2
    couch_iso2bb_vector = Vector(-1, -0.4, 0)
    variable_axes = {0: 'Reference'}

    @classmethod
    def setUpClass(cls):
        cls.wl = WinstonLutz.from_demo_images()
