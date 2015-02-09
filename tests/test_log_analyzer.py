from unittest import TestCase
import time

from pylinac.log_analyzer import *


class Test_Log_Loading(TestCase):
    """Tests of dynalog files, mostly using the demo file."""
    test_dir = osp.join(osp.dirname(__file__), 'test_files', 'MLC logs')

    def test_loading(self):
        """Test that loading the badly-named dynalog still loads and identifies properly."""
        test_tlog = osp.join(self.test_dir, 'tlogs', "qqq2106_4DC Treatment_JS0_TX_20140712095629.bin")
        # should produce no errors
        # load method 1
        MachineLog(test_tlog)
        # load method 2
        log = MachineLog()
        self.assertFalse(log.log_is_loaded)
        log.load(test_tlog)
        self.assertTrue(log.log_is_loaded)

        # throw an error for files that aren't logs
        not_a_file = test_tlog.replace(".bin", 'blahblah.bin')
        self.assertRaises(FileExistsError, MachineLog, not_a_file)
        not_a_log = osp.join(osp.dirname(__file__), 'test_files', 'VMAT', 'DRGSmlc-105-example.dcm')
        self.assertRaises(IOError, MachineLog, not_a_log)

    def test_dynalog_loading(self):
        a_file = osp.join(self.test_dir, 'dlogs', 'Adlog1.dlg')
        MachineLog(a_file)

        b_file = osp.join(self.test_dir, 'dlogs', 'Bdlog1.dlg')
        MachineLog(b_file)

        a_but_not_b_dir = osp.join(self.test_dir, 'a_no_b_dir', 'Adlog1.dlg')
        self.assertRaises(FileNotFoundError, MachineLog, a_but_not_b_dir)
        b_but_not_a_dir = osp.join(self.test_dir, 'b_no_a_dir', 'Bdlog1.dlg')
        self.assertRaises(FileNotFoundError, MachineLog, b_but_not_a_dir)

        bad_name_dlg = osp.join(self.test_dir, 'bad_names', 'bad_name_dlg.dlg')
        self.assertRaises(ValueError, MachineLog, bad_name_dlg)


class Test_Dynalog_Demo(TestCase):
    """Tests having to do with trajectory log files."""
    @classmethod
    def setUpClass(cls):
        cls.log = MachineLog()
        cls.log.load_demo_dynalog()

    def test_type(self):
        """Test all kinds of things about the dynalog demo."""
        self.assertTrue(self.log.log_type == log_types['dlog'])

    def test_header(self):
        """Test header info of the dynalog; ensures data integrity."""
        header = self.log.header
        self.assertEqual(header.version, ['B'])
        self.assertEqual(header.patient_name, ['Clinac4 QA', '', 'Clinac4 QA'])
        self.assertEqual(header.plan_filename, ['1.2.246.352.71.5.1399119341.107477.20110923193623', '21'])
        self.assertEqual(header.tolerance, 102)
        self.assertEqual(header.num_mlc_leaves, 120)
        self.assertEqual(header.clinac_scale, [' 1'])

    def test_axis_data(self):
        """Sample a few points from the axis data to ensure integrity."""
        axis_data = self.log.axis_data
        # properties
        self.assertEqual(axis_data.num_beamholds, 20)
        self.assertEqual(axis_data.num_snapshots, 99)
        self.assertEqual(len(axis_data.beam_hold.actual), 99)
        # MU data
        self.assertEqual(axis_data.mu.actual[0], 0)
        self.assertEqual(axis_data.mu.actual[-1], 25000)
        self.assertEqual(axis_data.mu.expected[0], 0)
        self.assertEqual(axis_data.mu.expected[-1], 25000)
        # jaws
        self.assertEqual(axis_data.jaws.x1.actual[0], 8)
        self.assertEqual(axis_data.jaws.y1.actual[-1], 20)
        self.assertRaises(AttributeError, axis_data.jaws.x2.plot_expected)

    def test_mlc(self):
        """Test integrity of MLC data & methods."""
        mlc = self.log.axis_data.mlc
        self.assertEqual(mlc.num_leaves, 120)
        self.assertEqual(mlc.num_pairs, 60)
        self.assertEqual(mlc.num_snapshots, 21)  # snapshots where beam was on
        self.assertEqual(len(mlc.moving_leaves), 60)
        self.assertFalse(mlc.hdmlc)

    def test_mlc_positions(self):
        """Test some MLC positions."""
        mlc = self.log.axis_data.mlc
        self.assertAlmostEqual(mlc.leaf_axes[1].actual[0], 7.564, delta=0.001)
        self.assertAlmostEqual(mlc.leaf_axes[120].expected[-1], -4.994, delta=0.001)
        self.assertAlmostEqual(mlc.leaf_axes[1].difference[0], 0, delta=0.1)

    def test_mlc_leafpair_moved(self):
        mlc = self.log.axis_data.mlc
        self.assertTrue(mlc.leaf_moved(9))
        self.assertFalse(mlc.leaf_moved(8))
        self.assertTrue(mlc.pair_moved(3))

    def test_RMS_error(self):
        mlc = self.log.axis_data.mlc
        self.assertAlmostEqual(mlc.get_RMS_avg(), 0.0373, delta=0.001)
        self.assertAlmostEqual(mlc.get_RMS_avg(bank='a'), 0.0375, delta=0.001)
        self.assertAlmostEqual(mlc.get_RMS_avg(only_moving_leaves=True), 0.074, delta=0.01)
        self.assertAlmostEqual(mlc.get_RMS_max(), 0.0756, delta=0.001)
        self.assertAlmostEqual(mlc.get_RMS_percentile(), 0.0754, delta=0.001)
        self.assertAlmostEqual(len(mlc.get_RMS('b')), 60)
        self.assertAlmostEqual(mlc.get_RMS((1,3)).mean(), 0.0717, delta=0.001)
        self.assertAlmostEqual(mlc.create_error_array((2,3), False).mean(), 0.034, delta=0.001)

    def test_under_jaws(self):
        mlc = self.log.axis_data.mlc
        self.assertFalse(mlc.leaf_under_y_jaw(4))


class Test_Dlog_Fluence(TestCase):
    def setUp(self):
        self.log = MachineLog()
        self.log.load_demo_dynalog()

    def test_demo(self):
        self.log.run_dlog_demo()

    def test_fluence(self):
        fluence = self.log.fluence

        self.assertFalse(fluence.actual.map_calced)
        self.assertFalse(fluence.expected.map_calced)
        self.assertFalse(fluence.gamma.map_calced)
        self.assertRaises(AttributeError, fluence.actual.plot_map)

        # do repeating fluence calcs; ensure semi-lazy property
        start = time.time()
        fluence.actual.calc_map()
        end = time.time()
        first_calc_time = end - start
        start = time.time()
        fluence.actual.calc_map()
        end = time.time()
        second_calc_time = end - start
        self.assertLess(second_calc_time, first_calc_time)

        # same for gamma
        start = time.time()
        fluence.gamma.calc_map(resolution=0.15)
        end = time.time()
        first_calc_time = end - start
        start = time.time()
        fluence.gamma.calc_map(resolution=0.15)
        end = time.time()
        second_calc_time = end - start
        self.assertLess(second_calc_time, first_calc_time)

        self.assertAlmostEqual(fluence.gamma.pass_prcnt, 99.85, delta=0.1)
        self.assertAlmostEqual(fluence.gamma.avg_gamma, 0.019, delta=0.005)
        self.assertAlmostEqual(fluence.gamma.histogram[0], 155804, delta=100)


class Test_Tlog_Demo(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.log = MachineLog()
        cls.log.load_demo_trajectorylog()

    def test_demo(self):
        MachineLog().run_tlog_demo()

    def test_type(self):
        self.assertTrue(self.log.log_type, log_types['tlog'])

    def test_header(self):
        header = self.log.header
        self.assertEqual(header.header, 'VOSTL')
        self.assertEqual(header.version, 2.1)
        self.assertEqual(header.header_size, 1024)
        self.assertEqual(header.sampling_interval, 20)
        self.assertEqual(header.num_axes, 14)
        self.assertEqual(header.axis_enum.size, 14)
        self.assertEqual(header.samples_per_axis.size, 14)
        self.assertEqual(header.samples_per_axis[-1], 122)
        self.assertEqual(header.num_mlc_leaves, 120)
        self.assertEqual(header.clinac_scale, 1)
        self.assertEqual(header.num_subbeams, 2)
        self.assertEqual(header.is_truncated, 0)
        self.assertEqual(header.num_snapshots, 5200)
        self.assertEqual(header.mlc_model, 3)

    def test_axis_data(self):
        axis_data = self.log.axis_data
        self.assertAlmostEqual(axis_data.collimator.actual[0], 180, delta=0.1)
        self.assertAlmostEqual(axis_data.mu.difference[1], 0.0000337, delta=0.01)
        self.assertAlmostEqual(axis_data.gantry.expected[2], 310, delta=0.1)

    def test_mlc(self):
        mlc = self.log.axis_data.mlc
        self.assertTrue(mlc.hdmlc)
        self.assertEqual(mlc.num_leaves, 120)
        self.assertEqual(mlc.num_pairs, 60)
        self.assertEqual(mlc.num_snapshots, 1021)

        log_no_exclusion = MachineLog()
        log_no_exclusion.load_demo_trajectorylog(exclude_beam_off=False)
        self.assertEqual(log_no_exclusion.axis_data.mlc.num_snapshots, 5200)

    def test_under_jaws(self):
        mlc = self.log.axis_data.mlc
        self.assertTrue(mlc.leaf_under_y_jaw(4))

    def test_RMS_error(self):
        mlc = self.log.axis_data.mlc
        self.assertAlmostEqual(mlc.get_RMS_avg(), 0.001, delta=0.001)
        self.assertAlmostEqual(mlc.get_RMS_max(), 0.00216, delta=0.001)
        self.assertAlmostEqual(mlc.get_RMS_percentile(), 0.00196, delta=0.001)


class Test_Tlog_Fluence(TestCase):
    def setUp(self):
        self.log = MachineLog()
        self.log.load_demo_trajectorylog()

    def test_fluence(self):
        fluence = self.log.fluence
        fluence.gamma.calc_map()
        self.assertAlmostEqual(fluence.gamma.pass_prcnt, 100, delta=0.1)
        self.assertAlmostEqual(fluence.gamma.avg_gamma, 0.001, delta=0.005)
        self.assertAlmostEqual(fluence.gamma.histogram[0], 240000, delta=100)