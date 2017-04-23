import os.path as osp
from unittest import TestCase

import matplotlib.pyplot as plt

from pylinac import WinstonLutz
from pylinac.winston_lutz import GANTRY, COLLIMATOR, COUCH
from pylinac.core.geometry import Vector, vector_is_close
from tests import TEST_BANK_DIR
from tests.utils import save_file, LoadingTestBase, LocationMixin


class TestWLLoading(LoadingTestBase, TestCase):
    klass = WinstonLutz
    demo_load_method = 'from_demo_images'
    url = 'winston_lutz.zip'


class GeneralTests(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.wl = WinstonLutz.from_demo_images()

    def test_run_demo(self):
        WinstonLutz.run_demo()  # shouldn't raise

    def test_results(self):
        print(self.wl.results())  # shouldn't raise


class TestPlottingSaving(TestCase):

    def setUp(self):
        self.wl = WinstonLutz.from_demo_images()

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def test_plot(self):
        self.wl.plot_images()  # shouldn't raise

    def test_save(self):
        save_file(self.wl.save_summary)
        save_file(self.wl.save_images)


class WinstonLutzMixin(LocationMixin):
    dir_location = osp.join(TEST_BANK_DIR, 'Winston-Lutz')
    num_images = 0
    zip = True
    gantry_deviation = None
    gantry_iso_size = 0
    gantry_iso2bb_vector = Vector
    collimator_iso_size = 0
    collimator_iso2bb_vector = Vector
    couch_iso_size = 0
    couch_iso2bb_vector = Vector
    cax2bb_max_distance = 0
    cax2bb_median_distance = 0
    variable_axes = {0: 'Reference'}  # fill with as many {image: axis} pairs as desired

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        if cls.zip:
            cls.wl = WinstonLutz.from_zip(filename)
        else:
            cls.wl = WinstonLutz(filename)

    def test_number_of_images(self):
        self.assertEqual(len(self.wl.images), self.num_images)

    def test_variable_axes(self):
        for idx, axis in self.variable_axes.items():
            self.assertEqual(axis, self.wl.images[idx].variable_axis)

    def test_gantry_iso(self):
        # test iso size
        self.assertAlmostEqual(self.wl.gantry_iso_size, self.gantry_iso_size, delta=0.2)
        # test iso vector
        self.assertTrue(vector_is_close(self.wl.gantry_iso2bb_vector, self.gantry_iso2bb_vector),
                        msg="{} was != {} within {}".format(self.wl.gantry_iso2bb_vector, self.gantry_iso2bb_vector, 0.2))

    def test_collimator_iso(self):
        # test iso size
        if self.collimator_iso_size is not None:
            self.assertAlmostEqual(self.wl.collimator_iso_size, self.collimator_iso_size, delta=0.2)
        # test iso vector
        if self.collimator_iso2bb_vector is not None:
            self.assertTrue(vector_is_close(self.wl.collimator_iso2bb_vector, self.collimator_iso2bb_vector),
                            msg="{} was != {}".format(self.wl.collimator_iso2bb_vector, self.collimator_iso2bb_vector))

    def test_couch_iso(self):
        # test iso size
        if self.couch_iso_size is not None:
            self.assertAlmostEqual(self.wl.couch_iso_size, self.couch_iso_size, delta=0.2)
        # test iso vector
        if self.couch_iso2bb_vector is not None:
            self.assertTrue(vector_is_close(self.wl.couch_iso2bb_vector, self.couch_iso2bb_vector),
                            msg="{} was != {}".format(self.wl.couch_iso2bb_vector, self.couch_iso2bb_vector))

    def test_axis_deviation(self):
        if self.gantry_deviation is not None:
            self.assertAlmostEqual(self.wl.axis_rms_deviation(GANTRY, value='range'), self.gantry_deviation, delta=0.2)


class WLDemo(WinstonLutzMixin, TestCase):
    num_images = 17
    gantry_iso_size = 1
    gantry_iso2bb_vector = Vector(-0.4, 0.2, -0.5)
    collimator_iso_size = 1.2
    collimator_iso2bb_vector = Vector(-0.15, 0.4, 0)
    couch_iso_size = 3.5
    couch_iso2bb_vector = Vector(1.3, -0.7, 0)
    variable_axes = {0: 'Reference'}
    gantry_deviation = 0.9

    @classmethod
    def setUpClass(cls):
        cls.wl = WinstonLutz.from_demo_images()
