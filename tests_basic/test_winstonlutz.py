import os.path as osp
import tempfile
from unittest import TestCase

import matplotlib.pyplot as plt
import numpy as np

from pylinac import WinstonLutz
from pylinac.winston_lutz import GANTRY, COLLIMATOR, COUCH, REFERENCE, GB_COMBO, GBP_COMBO, EPID
from pylinac.core.geometry import Vector, vector_is_close
from tests_basic import TEST_BANK_DIR, TEST_FILES_DIR
from tests_basic.utils import save_file, LoadingTestBase, LocationMixin


TEST_DIR = osp.join(TEST_FILES_DIR, 'Winston-Lutz')


class TestWLLoading(LoadingTestBase, TestCase):
    klass = WinstonLutz
    demo_load_method = 'from_demo_images'
    url = 'winston_lutz.zip'

    def test_loading_1_image_fails(self):
        with self.assertRaises(ValueError):
            WinstonLutz(osp.join(TEST_DIR, 'lutz', '1_image'))

    def test_invalid_dir(self):
        with self.assertRaises(ValueError):
            WinstonLutz(r'nonexistant/dir')


class GeneralTests(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.wl = WinstonLutz.from_demo_images()

    def test_run_demo(self):
        WinstonLutz.run_demo()  # shouldn't raise

    def test_results(self):
        print(self.wl.results())  # shouldn't raise

    def test_bb_shift_instructions(self):
        move = self.wl.bb_shift_instructions()
        self.assertTrue("RIGHT" in move)

        move = self.wl.bb_shift_instructions(couch_vrt=-2, couch_lat=1, couch_lng=100)
        self.assertTrue("RIGHT" in move)
        self.assertTrue("VRT" in move)


class TestPublishPDF(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.wl = WinstonLutz.from_demo_images()

    def test_publish_pdf(self):
        # normal publish; shouldn't raise
        with tempfile.TemporaryFile() as t:
            self.wl.publish_pdf(t)

    def test_publish_w_metadat_and_notes(self):
        with tempfile.TemporaryFile() as t:
            self.wl.publish_pdf(t, notes='stuff', metadata={'Unit': 'TB1'})


class TestPlottingSaving(TestCase):

    def setUp(self):
        self.wl = WinstonLutz.from_demo_images()

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def test_plot(self):
        self.wl.plot_images()  # shouldn't raise
        self.wl.plot_images(axis=GANTRY)
        self.wl.plot_images(axis=COLLIMATOR)
        self.wl.plot_images(axis=COUCH)
        self.wl.plot_images(axis=GB_COMBO)
        self.wl.plot_images(axis=GBP_COMBO)

    def test_save(self):
        save_file(self.wl.save_summary)
        save_file(self.wl.save_images)

    def test_plot_wo_all_axes(self):
        # test that analyzing images w/o gantry images doesn't fail
        wl_zip = osp.join(TEST_FILES_DIR, 'Winston-Lutz', 'Naming.zip')
        wl = WinstonLutz.from_zip(wl_zip, use_filenames=True)
        wl.plot_summary()  # shouldn't raise


class WinstonLutzMixin(LocationMixin):
    dir_location = osp.join(TEST_BANK_DIR, 'Winston-Lutz')
    num_images = 0
    zip = True
    gantry_iso_size = 0
    collimator_iso_size = 0
    couch_iso_size = 0
    cax2bb_max_distance = 0
    cax2bb_median_distance = 0
    epid_deviation = None
    bb_shift_vector = Vector()  # vector to place BB at iso
    axis_of_rotation = {0: 'Reference'}  # fill with as many {image#: known_axis_of_rotation} pairs as desired
    print_results = False
    use_filenames = False

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        if cls.zip:
            cls.wl = WinstonLutz.from_zip(filename, use_filenames=cls.use_filenames)
        else:
            cls.wl = WinstonLutz(filename, use_filenames=cls.use_filenames)
        if cls.print_results:
            print(cls.wl.results())
            print(cls.wl.bb_shift_vector)

    def test_number_of_images(self):
        self.assertEqual(len(self.wl.images), self.num_images)

    def test_gantry_iso(self):
        # test iso size
        self.assertAlmostEqual(self.wl.gantry_iso_size, self.gantry_iso_size, delta=0.15)

    def test_collimator_iso(self):
        # test iso size
        if self.collimator_iso_size is not None:
            self.assertAlmostEqual(self.wl.collimator_iso_size, self.collimator_iso_size, delta=0.15)

    def test_couch_iso(self):
        # test iso size
        if self.couch_iso_size is not None:
            self.assertAlmostEqual(self.wl.couch_iso_size, self.couch_iso_size, delta=0.15)

    def test_epid_deviation(self):
        if self.epid_deviation is not None:
            self.assertAlmostEqual(max(self.wl.axis_rms_deviation(EPID)), self.epid_deviation, delta=0.15)

    def test_bb_max_distance(self):
        self.assertAlmostEqual(self.wl.cax2bb_distance(metric='max'), self.cax2bb_max_distance, delta=0.15)

    def test_bb_median_distance(self):
        self.assertAlmostEqual(self.wl.cax2bb_distance(metric='median'), self.cax2bb_median_distance, delta=0.1)

    def test_bb_shift_vector(self):
        self.assertTrue(vector_is_close(self.wl.bb_shift_vector, self.bb_shift_vector, delta=0.15), msg="The vector {} is not sufficiently close to vector {}".format(self.wl.bb_shift_vector, self.bb_shift_vector))

    def test_known_axis_of_rotation(self):
        for idx, axis in self.axis_of_rotation.items():
            self.assertEqual(axis, self.wl.images[idx].variable_axis)


class WLDemo(WinstonLutzMixin, TestCase):
    num_images = 17
    gantry_iso_size = 1
    collimator_iso_size = 1.2
    couch_iso_size = 2.3
    cax2bb_max_distance = 1.2
    cax2bb_median_distance = 0.7
    epid_deviation = 1.3
    axis_of_rotation = {0: 'Reference'}
    bb_shift_vector = Vector(x=0.4, y=-0.4, z=-0.2)

    @classmethod
    def setUpClass(cls):
        cls.wl = WinstonLutz.from_demo_images()


class WLPerfect30x8(WinstonLutzMixin, TestCase):
    """30x30mm field, 8mm BB"""
    dir_location = TEST_FILES_DIR
    file_path = ['Winston-Lutz', 'perfect_WL_30x8.zip']
    num_images = 10
    gantry_iso_size = 0
    collimator_iso_size = 0
    couch_iso_size = 0
    cax2bb_max_distance = 0
    cax2bb_median_distance = 0
    epid_deviation = 0
    bb_shift_vector = Vector()


class WLPerfect10x4(WinstonLutzMixin, TestCase):
    """10x10mm field, 4mm BB"""
    dir_location = TEST_FILES_DIR
    file_path = ['Winston-Lutz', 'perfect_WL_10x4.zip']
    num_images = 10
    gantry_iso_size = 0
    collimator_iso_size = 0
    couch_iso_size = 0
    cax2bb_max_distance = 0
    cax2bb_median_distance = 0
    epid_deviation = 0
    bb_shift_vector = Vector()


class WLNoisy30x5(WinstonLutzMixin, TestCase):
    """30x30mm field, 5mm BB. S&P noise added"""
    dir_location = TEST_FILES_DIR
    file_path = ['Winston-Lutz', 'noisy_WL_30x5.zip']
    num_images = 10
    gantry_iso_size = 0.08
    collimator_iso_size = 0
    couch_iso_size = 0
    cax2bb_max_distance = 0
    cax2bb_median_distance = 0
    epid_deviation = 0
    bb_shift_vector = Vector()


class WLLateral3mm(WinstonLutzMixin, TestCase):
    # verified independently
    dir_location = TEST_FILES_DIR
    file_path = ['Winston-Lutz', 'lat3mm.zip']
    num_images = 4
    gantry_iso_size = 0.5
    cax2bb_max_distance = 3.8
    cax2bb_median_distance = 2.3
    bb_shift_vector = Vector(x=-3.6, y=0.5, z=0.6)


class WLLongitudinal3mm(WinstonLutzMixin, TestCase):
    # verified independently
    dir_location = TEST_FILES_DIR
    file_path = ['Winston-Lutz', 'lng3mm.zip']
    num_images = 4
    gantry_iso_size = 0.5
    cax2bb_max_distance = 3.9
    cax2bb_median_distance = 3.7
    bb_shift_vector = Vector(x=-0.63, y=3.6, z=0.6)


class WLVertical3mm(WinstonLutzMixin, TestCase):
    dir_location = TEST_FILES_DIR
    file_path = ['Winston-Lutz', 'vrt3mm.zip']
    num_images = 4
    gantry_iso_size = 0.5
    cax2bb_max_distance = 3.8
    cax2bb_median_distance = 2.3
    bb_shift_vector = Vector(x=-0.5, y=0.5, z=3.6)
    print_results = True


class WLDontUseFileNames(WinstonLutzMixin, TestCase):
    dir_location = TEST_FILES_DIR
    file_path = ['Winston-Lutz', 'Naming.zip']
    num_images = 4
    gantry_iso_size = 0.3
    cax2bb_max_distance = 0.9
    cax2bb_median_distance = 0.8
    bb_shift_vector = Vector(x=-0.4, y=0.6, z=0.6)
    axis_of_rotation = {0: REFERENCE, 1: GANTRY, 2: GANTRY, 3: GANTRY}


class WLUseFileNames(WinstonLutzMixin, TestCase):
    dir_location = TEST_FILES_DIR
    file_path = ['Winston-Lutz', 'Naming.zip']
    use_filenames = True
    num_images = 4
    collimator_iso_size = 1.2
    cax2bb_max_distance = 0.9
    cax2bb_median_distance = 0.8
    bb_shift_vector = Vector(y=0.6)
    axis_of_rotation = {0: COLLIMATOR, 1: COLLIMATOR, 2: COLLIMATOR, 3: COLLIMATOR}


class WLBadFilenames(TestCase):

    def test_bad_filenames(self):
        # tests_basic that using filenames with incorrect syntax will fail
        wl_dir = osp.join(TEST_FILES_DIR, 'Winston-Lutz', 'Bad Names.zip')
        with self.assertRaises(ValueError):
            wl = WinstonLutz.from_zip(wl_dir, use_filenames=True)
