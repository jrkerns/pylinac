import io
import tempfile
from unittest import TestCase

import matplotlib.pyplot as plt

import pylinac
from pylinac import WinstonLutz
from pylinac.core.geometry import Vector, vector_is_close
from pylinac.winston_lutz import Axis, WinstonLutzResult
from tests_basic.utils import save_file, CloudFileMixin, get_folder_from_cloud_test_repo, \
    get_file_from_cloud_test_repo, FromDemoImageTesterMixin, FromURLTesterMixin

TEST_DIR = 'Winston-Lutz'


class TestWLLoading(TestCase, FromDemoImageTesterMixin, FromURLTesterMixin):
    klass = WinstonLutz
    demo_load_method = 'from_demo_images'
    url = 'winston_lutz.zip'

    def test_loading_1_image_fails(self):
        with self.assertRaises(ValueError):
            folder = get_folder_from_cloud_test_repo(['Winston-Lutz', 'lutz', '1_image'])
            WinstonLutz(folder)

    def test_invalid_dir(self):
        with self.assertRaises(ValueError):
            WinstonLutz(r'nonexistant/dir')

    def test_load_from_file_object(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, 'noisy_WL_30x5.zip'])
        ref_w = WinstonLutz.from_zip(path)
        ref_w.analyze()
        with open(path, 'rb') as f:
            w = WinstonLutz.from_zip(f)
            w.analyze()
        self.assertIsInstance(w, WinstonLutz)
        self.assertEqual(w.gantry_iso_size, ref_w.gantry_iso_size)

    def test_load_from_stream(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, 'noisy_WL_30x5.zip'])
        ref_w = WinstonLutz.from_zip(path)
        ref_w.analyze()
        with open(path, 'rb') as f:
            s = io.BytesIO(f.read())
            w = WinstonLutz.from_zip(s)
            w.analyze()
        self.assertIsInstance(w, WinstonLutz)
        self.assertEqual(w.gantry_iso_size, ref_w.gantry_iso_size)


class GeneralTests(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.wl = WinstonLutz.from_demo_images()
        cls.wl.analyze()

    def test_run_demo(self):
        WinstonLutz.run_demo()  # shouldn't raise

    def test_results(self):
        print(self.wl.results())  # shouldn't raise

    def test_not_yet_analyzed(self):
        wl = WinstonLutz.from_demo_images()
        with self.assertRaises(ValueError):
            wl.results()  # not yet analyzed

        with self.assertRaises(ValueError):
            wl.plot_images()

        with self.assertRaises(ValueError):
            wl.plot_summary()

    def test_str_or_enum(self):
        # shouldn't raise
        self.wl.plot_images('Gantry')
        self.wl.plot_images(Axis.GANTRY)

        self.wl.plot_axis_images('Gantry')
        self.wl.plot_axis_images(Axis.GANTRY)

    def test_bb_override(self):
        with self.assertRaises(ValueError):
            wl = pylinac.WinstonLutz.from_demo_images()
            wl.analyze(bb_size_mm=8)

    def test_bb_shift_instructions(self):
        move = self.wl.bb_shift_instructions()
        self.assertTrue("RIGHT" in move)

        move = self.wl.bb_shift_instructions(couch_vrt=-2, couch_lat=1, couch_lng=100)
        self.assertTrue("RIGHT" in move)
        self.assertTrue("VRT" in move)

    def test_results_data(self):
        data = self.wl.results_data()
        self.assertIsInstance(data, WinstonLutzResult)
        self.assertEqual(data.num_couch_images, self.wl._get_images(axis=(Axis.COUCH, Axis.REFERENCE))[0])
        self.assertEqual(data.max_2d_cax_to_epid_mm, self.wl.cax2epid_distance('max'))
        self.assertEqual(data.median_2d_cax_to_epid_mm, self.wl.cax2epid_distance('median'))

        data_dict = self.wl.results_data(as_dict=True)
        self.assertIn('pylinac_version', data_dict)
        self.assertEqual(data_dict['gantry_3d_iso_diameter_mm'], self.wl.gantry_iso_size)

    def test_bb_too_far_away_fails(self):
        """BB is >20mm from CAX"""
        file = get_file_from_cloud_test_repo([TEST_DIR, 'bb_too_far_away.zip'])
        wl = WinstonLutz.from_zip(file)
        with self.assertRaises(ValueError):
            wl.analyze()


class TestPublishPDF(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.wl = WinstonLutz.from_demo_images()
        cls.wl.analyze()

    def test_publish_pdf(self):
        # normal publish; shouldn't raise
        with tempfile.TemporaryFile() as t:
            self.wl.publish_pdf(t)

    def test_publish_w_metadata_and_notes(self):
        with tempfile.TemporaryFile() as t:
            self.wl.publish_pdf(t, notes='stuff', metadata={'Unit': 'TB1'})


class TestPlottingSaving(TestCase):

    def setUp(self):
        self.wl = WinstonLutz.from_demo_images()
        self.wl.analyze()

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def test_plot(self):
        self.wl.plot_images()  # shouldn't raise
        self.wl.plot_images(axis=Axis.GANTRY)
        self.wl.plot_images(axis=Axis.COLLIMATOR)
        self.wl.plot_images(axis=Axis.COUCH)
        self.wl.plot_images(axis=Axis.GB_COMBO)
        self.wl.plot_images(axis=Axis.GBP_COMBO)

    def test_save(self):
        save_file(self.wl.save_summary)
        save_file(self.wl.save_images)

    def test_plot_wo_all_axes(self):
        # test that analyzing images w/o gantry images doesn't fail
        wl_zip = get_file_from_cloud_test_repo([TEST_DIR, 'Naming.zip'])
        wl = WinstonLutz.from_zip(wl_zip, use_filenames=True)
        wl.analyze()
        wl.plot_summary()  # shouldn't raise


class WinstonLutzMixin(CloudFileMixin):
    dir_path = ['Winston-Lutz']
    num_images = 0
    zip = True
    bb_size = 5
    gantry_iso_size = 0
    collimator_iso_size = 0
    couch_iso_size = 0
    cax2bb_max_distance = 0
    cax2bb_median_distance = 0
    epid_deviation = None
    bb_shift_vector = Vector()  # vector to place BB at iso
    axis_of_rotation = {0: Axis.REFERENCE}  # fill with as many {image#: known_axis_of_rotation} pairs as desired
    print_results = False
    use_filenames = False

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        if cls.zip:
            cls.wl = WinstonLutz.from_zip(filename, use_filenames=cls.use_filenames)
        else:
            cls.wl = WinstonLutz(filename, use_filenames=cls.use_filenames)
        cls.wl.analyze(bb_size_mm=cls.bb_size)
        if cls.print_results:
            print(cls.wl.results())
            print(cls.wl.bb_shift_vector)

    def test_number_of_images(self):
        self.assertEqual(self.num_images, len(self.wl.images))

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
            self.assertAlmostEqual(max(self.wl.axis_rms_deviation(Axis.EPID)), self.epid_deviation, delta=0.15)

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
    axis_of_rotation = {0: Axis.REFERENCE}
    bb_shift_vector = Vector(x=0.4, y=-0.4, z=-0.2)
    delete_file = False

    @classmethod
    def setUpClass(cls):
        cls.wl = WinstonLutz.from_demo_images()
        cls.wl.analyze()


class WLPerfect30x8(WinstonLutzMixin, TestCase):
    """30x30mm field, 8mm BB"""
    file_name = 'perfect_WL_30x8.zip'
    num_images = 4
    gantry_iso_size = 0
    collimator_iso_size = 0
    couch_iso_size = 0
    cax2bb_max_distance = 0
    cax2bb_median_distance = 0
    epid_deviation = 0
    bb_shift_vector = Vector()


class WLPerfect30x2(WinstonLutzMixin, TestCase):
    """30x30mm field, 2mm BB"""
    file_name = 'perfect_WL_30x2mm.zip'
    num_images = 4
    gantry_iso_size = 0
    collimator_iso_size = 0
    couch_iso_size = 0
    cax2bb_max_distance = 0
    cax2bb_median_distance = 0
    epid_deviation = 0
    bb_shift_vector = Vector()
    bb_size = 2


class WLPerfect10x4(WinstonLutzMixin, TestCase):
    """10x10mm field, 4mm BB"""
    file_name = 'perfect_WL_10x4.zip'
    num_images = 4
    gantry_iso_size = 0
    collimator_iso_size = 0
    couch_iso_size = 0
    cax2bb_max_distance = 0
    cax2bb_median_distance = 0
    epid_deviation = 0
    bb_shift_vector = Vector()


class WLNoisy30x5(WinstonLutzMixin, TestCase):
    """30x30mm field, 5mm BB. S&P noise added"""
    file_name = 'noisy_WL_30x5.zip'
    num_images = 4
    gantry_iso_size = 0.08
    collimator_iso_size = 0
    couch_iso_size = 0
    cax2bb_max_distance = 0
    cax2bb_median_distance = 0
    epid_deviation = 0
    bb_shift_vector = Vector()


class WLLateral3mm(WinstonLutzMixin, TestCase):
    # verified independently
    file_name = 'lat3mm.zip'
    num_images = 4
    gantry_iso_size = 0.5
    cax2bb_max_distance = 3.8
    cax2bb_median_distance = 2.3
    bb_shift_vector = Vector(x=-3.6, y=0.5, z=0.6)


class WLLongitudinal3mm(WinstonLutzMixin, TestCase):
    # verified independently
    file_name = 'lng3mm.zip'
    num_images = 4
    gantry_iso_size = 0.5
    cax2bb_max_distance = 3.9
    cax2bb_median_distance = 3.7
    bb_shift_vector = Vector(x=-0.63, y=3.6, z=0.6)


class WLVertical3mm(WinstonLutzMixin, TestCase):
    file_name = 'vrt3mm.zip'
    num_images = 4
    gantry_iso_size = 0.5
    cax2bb_max_distance = 3.8
    cax2bb_median_distance = 2.3
    bb_shift_vector = Vector(x=-0.5, y=0.5, z=3.6)
    print_results = True


class WLDontUseFileNames(WinstonLutzMixin, TestCase):
    file_name = 'Naming.zip'
    num_images = 4
    gantry_iso_size = 0.3
    cax2bb_max_distance = 0.9
    cax2bb_median_distance = 0.8
    bb_shift_vector = Vector(x=-0.4, y=0.6, z=0.6)
    axis_of_rotation = {0: Axis.REFERENCE, 1: Axis.GANTRY, 2: Axis.GANTRY, 3: Axis.GANTRY}


class WLUseFileNames(WinstonLutzMixin, TestCase):
    file_name = 'Naming.zip'
    use_filenames = True
    num_images = 4
    collimator_iso_size = 1.2
    cax2bb_max_distance = 0.9
    cax2bb_median_distance = 0.8
    bb_shift_vector = Vector(y=0.6)
    axis_of_rotation = {0: Axis.COLLIMATOR, 1: Axis.COLLIMATOR, 2: Axis.COLLIMATOR, 3: Axis.COLLIMATOR}


class WLBadFilenames(TestCase):

    def test_bad_filenames(self):
        # tests_basic that using filenames with incorrect syntax will fail
        wl_dir = get_file_from_cloud_test_repo([TEST_DIR, 'Bad-Names.zip'])
        with self.assertRaises(ValueError):
            wl = WinstonLutz.from_zip(wl_dir, use_filenames=True)
            wl.analyze()
