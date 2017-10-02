import os.path as osp
import os
from unittest import TestCase
import shutil


from pylinac.log_analyzer import MachineLogs, STATIC_IMRT, DYNAMIC_IMRT, \
    VMAT, anonymize, TrajectoryLog, Dynalog, load_log, DynalogMatchError, NotADynalogError, IMAGING
from tests.utils import save_file, LoadingTestBase, LocationMixin

TEST_DIR = osp.join(osp.dirname(__file__), 'test_files', 'MLC logs')
ANONYMOUS_SOURCE_FOLDER = osp.join(TEST_DIR, '_anonbase')
ANONYMOUS_DEST_FOLDER = osp.join(TEST_DIR, 'anonymous')


class TestAnonymizeFunction(TestCase):
    """Test the anonymization method."""

    def setUp(self):
        # move over files from other directory, since the filenames get overridden
        for file in os.listdir(ANONYMOUS_SOURCE_FOLDER):
            basefile = osp.join(ANONYMOUS_SOURCE_FOLDER, file)
            destfile = osp.join(ANONYMOUS_DEST_FOLDER, file)
            if not osp.isfile(destfile):
                shutil.copy(basefile, ANONYMOUS_DEST_FOLDER)

    @classmethod
    def tearDownClass(cls):
        # remove files from anonymous folder
        files = os.listdir(ANONYMOUS_DEST_FOLDER)
        files.remove('dummy.txt')
        for file in files:
            file = osp.join(ANONYMOUS_DEST_FOLDER, file)
            os.remove(file)

    def test_anonymize_function(self):
        # shouldn't raise
        anonymize(osp.join(ANONYMOUS_DEST_FOLDER, 'A1234_patientid.dlg'))
        anonymize(ANONYMOUS_DEST_FOLDER, inplace=False)
        anonymize(ANONYMOUS_DEST_FOLDER, recursive=False)

    def test_dynalog(self):
        # test making an anonymized copy
        dlog_file = osp.join(ANONYMOUS_DEST_FOLDER, 'A1234_patientid.dlg')
        dlog = Dynalog(dlog_file)
        dlog.anonymize()

        # test doing inplace anonymization
        files = dlog.anonymize(inplace=True, suffix='inplace')
        for file in files:
            self.assertTrue('inplace' in file)

    def test_destination(self):
        tlog_file = osp.join(ANONYMOUS_DEST_FOLDER, 'PatientID_4DC Treatment_JST90_TX_20140712094246.bin')
        tlog = TrajectoryLog(tlog_file)
        tlog.anonymize(destination=ANONYMOUS_DEST_FOLDER)  # shouldn't raise

    def test_bad_name(self):
        """Test that a log with a bad name (no underscore) fails gracefully."""
        dlog_file = osp.join(ANONYMOUS_DEST_FOLDER, 'A1234patientid.dlg')
        dlog = Dynalog(dlog_file)
        with self.assertRaises(NameError):
            dlog.anonymize()


class TestLogPlottingSavingMixin:
    """Test the plotting methods and plot saving methods."""

    def test_plot_axes(self):
        for methodname in ('plot_actual', 'plot_expected', 'plot_difference'):
            method = getattr(self.log.axis_data.mlc.leaf_axes[10], methodname)
            method()  # shouldn't raise

    def test_save_axes(self):
        for methodname in ('save_plot_actual', 'save_plot_expected', 'save_plot_difference'):
            # save matplotlib figures
            method = getattr(self.log.axis_data.mlc.leaf_axes[10], methodname)
            save_file(method)

            # save MPLD3 HTML
            save_file(method, interactive=True, as_file_object='str')

    def test_fluence_plotting(self):
        # raise error if map hasn't yet been calc'ed.
        with self.assertRaises(AttributeError):
            self.log.fluence.actual.plot_map()

        self.log.fluence.actual.calc_map()
        self.log.fluence.actual.plot_map()
        self.log.fluence.gamma.calc_map()
        self.log.fluence.gamma.histogram()
        self.log.fluence.gamma.plot_histogram()
        self.log.fluence.gamma.plot_passfail_map()

    def test_saving_fluence_plots(self):
        self.log.fluence.gamma.calc_map()
        save_file(self.log.fluence.gamma.save_map)
        save_file(self.log.fluence.gamma.save_histogram)

    def test_save_summary(self):
        self.log.fluence.gamma.calc_map()
        save_file(self.log.save_summary)


class TestLogBase:
    klass = object

    def setUp(self):
        self.log = self.klass.from_demo()
        # move over files from other directory, since the filenames get overridden
        for file in os.listdir(ANONYMOUS_SOURCE_FOLDER):
            basefile = osp.join(ANONYMOUS_SOURCE_FOLDER, file)
            destfile = osp.join(ANONYMOUS_DEST_FOLDER, file)
            if not osp.isfile(destfile):
                shutil.copy(basefile, ANONYMOUS_DEST_FOLDER)

    @classmethod
    def tearDownClass(cls):
        # remove files from anonymous folder
        files = os.listdir(ANONYMOUS_DEST_FOLDER)
        files.remove('dummy.txt')
        for file in files:
            file = osp.join(ANONYMOUS_DEST_FOLDER, file)
            os.remove(file)

    def test_run_demo(self):
        self.log.run_demo()

    def test_anonymize(self):
        log_file = osp.join(ANONYMOUS_DEST_FOLDER, self.anon_file)
        log = self.klass(log_file)

        files = log.anonymize(inplace=True, suffix='inplace')
        # self.assertIsInstance(files, list)
        for file in files:
            self.assertTrue('inplace' in file)


class TestTrajectoryLog(TestLogPlottingSavingMixin, LoadingTestBase, TestLogBase, TestCase):
    klass = TrajectoryLog
    demo_load_method = 'from_demo'
    url = 'Tlog.bin'
    anon_file = 'PatientID_4DC Treatment_JST90_TX_20140712094246.bin'

    def test_not_logs(self):
        # throw an error for files that aren't logs
        test_tlog = osp.join(TEST_DIR, 'tlogs', "Anonymous_4DC Treatment_JS0_TX_20140712095629.bin")
        not_a_file = test_tlog.replace(".bin", 'blahblah.bin')
        self.assertRaises(IOError, TrajectoryLog, not_a_file)
        not_a_log = osp.join(osp.dirname(__file__), 'test_files', 'VMAT', 'DRGSmlc-105-example.dcm')
        self.assertRaises(IOError, TrajectoryLog, not_a_log)

    def test_save_to_csv(self):
        save_file(self.log.to_csv)

    def test_txt_file_also_loads_if_around(self):
        # has a .txt file
        log_with_txt = osp.join(TEST_DIR, 'mixed_types', "Anonymous_4DC Treatment_JST90_TX_20140712094246.bin")

        log = TrajectoryLog(log_with_txt)
        self.assertIsNotNone(log.txt)
        self.assertIsInstance(log.txt, dict)
        self.assertEqual(log.txt['Patient ID'], 'Anonymous')

        # DOESN'T have a txt file
        log_no_txt = osp.join(TEST_DIR, 'tlogs', "Anonymous_4DC Treatment_JS0_TX_20140712095629.bin")

        log = TrajectoryLog(log_no_txt)
        self.assertIsNone(log.txt)


class TestDynalog(TestLogPlottingSavingMixin, LoadingTestBase, TestLogBase, TestCase):
    klass = Dynalog
    demo_load_method = 'from_demo'
    anon_file = 'A1234_patientid.dlg'

    def test_loading_can_find_paired_file(self):
        # shouldn't raise since it can find B-file
        a_file = osp.join(TEST_DIR, 'dlogs', 'Adlog1.dlg')
        Dynalog(a_file)

        # ditto for A-file
        b_file = osp.join(TEST_DIR, 'dlogs', 'Bdlog1.dlg')
        Dynalog(b_file)

    def test_loading_bad_names(self):
        a_but_not_b_dir = osp.join(TEST_DIR, 'a_no_b_dir', 'Adlog1.dlg')
        self.assertRaises(DynalogMatchError, Dynalog, a_but_not_b_dir)

        b_but_not_a_dir = osp.join(TEST_DIR, 'b_no_a_dir', 'Bdlog1.dlg')
        self.assertRaises(DynalogMatchError, Dynalog, b_but_not_a_dir)

        bad_name_dlg = osp.join(TEST_DIR, 'bad_names', 'bad_name_dlg.dlg')
        self.assertRaises(ValueError, Dynalog, bad_name_dlg)


class TestIndividualLogBase(LocationMixin):
    """Mixin to use when testing a single machine log; must be mixed with unittest.TestCase."""
    num_mlc_leaves = 120
    num_snapshots = 0
    num_beamholds = 0
    num_moving_leaves = 0
    treatment_type = ''
    static_axes = []
    moving_axes = []
    leaf_move_status = {'moving': tuple(), 'static': tuple()}
    average_rms = 0
    maximum_rms = 0
    average_gamma = 0
    percent_pass_gamma = 100
    mu_delivered = 0

    @classmethod
    def setUpClass(cls):
        cls.log = load_log(cls.get_filename())
        if cls.log.treatment_type != IMAGING:
            cls.log.fluence.gamma.calc_map()

    def test_num_leaves(self):
        """Test the number of MLC leaves and pairs."""
        self.assertEqual(self.log.header.num_mlc_leaves, self.num_mlc_leaves)

    def test_treatment_type(self):
        """Test the treatment type."""
        self.assertEqual(self.treatment_type, self.log.treatment_type)

    def test_rms_error(self):
        """Test the average and maximum RMS errors."""
        self.assertAlmostEqual(self.log.axis_data.mlc.get_RMS_avg(), self.average_rms, delta=0.01)
        self.assertAlmostEqual(self.log.axis_data.mlc.get_RMS_max(), self.maximum_rms, delta=0.01)

    def test_fluence_gamma(self):
        """Test gamma results for fluences."""
        if self.log.treatment_type != IMAGING:
            self.assertAlmostEqual(self.log.fluence.gamma.avg_gamma, self.average_gamma, delta=0.02)
            self.assertAlmostEqual(self.log.fluence.gamma.pass_prcnt, self.percent_pass_gamma, delta=0.1)

    def test_mu_delivered(self):
        """Test the number of MU delivered during the log."""
        self.assertAlmostEqual(self.log.axis_data.mu.actual[-1], self.mu_delivered, delta=1)

    def test_static_axes(self):
        """Test that certain axes did not move during treatment."""
        for axis_name in self.static_axes:
            axis = getattr(self.log.axis_data, axis_name)
            self.assertFalse(axis.moved)

    def test_leaf_moved_status(self):
        """Test that the given leaves either moved or did not move."""
        moving_leaves = self.leaf_move_status['moving']
        for leaf in moving_leaves:
            self.assertTrue(self.log.axis_data.mlc.leaf_moved(leaf))

        static_leaves = self.leaf_move_status['static']
        for leaf in static_leaves:
            self.assertFalse(self.log.axis_data.mlc.leaf_moved(leaf))


class TestIndividualTrajectoryLog(TestIndividualLogBase):
    version = 2.1  # or 3.0
    header = 'VOSTL'
    header_size = 1024
    sampling_interval = 20
    num_axes = 14
    axis_scale = 1
    num_subbeams = 0
    is_truncated = 0
    mlc_model = 2
    first_subbeam_data = {'gantry_angle': 0, 'collimator_angle': 0, 'jaw_x1': 0, 'jaw_x2': 0, 'jaw_y1': 0, 'jaw_y2': 0}

    def test_first_subbeam_data(self):
        """Test the first subbeam data."""
        first_subbeam = self.log.subbeams[0]
        for key, known_value in self.first_subbeam_data.items():
            axis = getattr(first_subbeam, key)
            self.assertAlmostEqual(known_value, axis.actual, delta=0.1)

    def test_header(self):
        """Test a few header values; depends on log type."""
        header = self.log.header
        self.assertEqual(header.version, self.version)
        self.assertEqual(header.header, self.header)
        self.assertEqual(header.header_size, self.header_size)
        self.assertEqual(header.sampling_interval, self.sampling_interval)
        self.assertEqual(header.num_axes, self.num_axes)
        self.assertEqual(header.axis_scale, self.axis_scale)
        self.assertEqual(header.num_subbeams, self.num_subbeams)
        self.assertEqual(header.is_truncated, self.is_truncated)
        self.assertEqual(header.mlc_model, self.mlc_model)

    def test_num_snapshots(self):
        """Test the number of snapshots in the log."""
        self.assertEqual(self.log.header.num_snapshots, self.num_snapshots)

    def test_num_beamholds(self):
        """Test the number of times the beam was held in the log."""
        self.assertEqual(self.log.num_beamholds, self.num_beamholds)


class TestIndividualDynalog(TestIndividualLogBase):
    tolerance = 102
    clinac_scale = 1
    mu_delivered = 25000
    version = 'B'

    def test_num_snapshots(self):
        """Test the number of snapshots in the log."""
        self.assertEqual(self.log.axis_data.num_snapshots, self.num_snapshots)

    def test_num_beamholds(self):
        """Test the number of times the beam was held in the log."""
        self.assertEqual(self.log.num_beamholds, self.num_beamholds)


class TestDynalogDemo(TestIndividualDynalog, TestCase):
    """Tests of the dynalog demo."""
    treatment_type = DYNAMIC_IMRT
    num_beamholds = 20
    num_snapshots = 99
    average_rms = 0.037
    maximum_rms = 0.076
    average_gamma = 0.47
    percent_pass_gamma = 91
    leaf_move_status = {'moving': (9, 3), 'static': (8, )}

    @classmethod
    def setUpClass(cls):
        cls.log = Dynalog.from_demo()
        cls.log.fluence.gamma.calc_map()


class TestTrajectoryLogDemo(TestIndividualTrajectoryLog, TestCase):
    """Tests for the demo trajectory log."""
    num_snapshots = 5200  # excluded: 1021
    num_subbeams = 2
    num_beamholds = 19
    mlc_model = 3
    treatment_type = DYNAMIC_IMRT
    static_axes = ['collimator']
    moving_axes = ['gantry']
    average_rms = 0.001
    maximum_rms = 0.002
    percent_pass_gamma = 100
    mu_delivered = 183
    first_subbeam_data = {'gantry_angle': 310, 'collimator_angle': 180, 'jaw_x1': 3.7, 'jaw_x2': 3.4, 'jaw_y1': 3.8,
                          'jaw_y2': 3.9}

    @classmethod
    def setUpClass(cls):
        cls.log = TrajectoryLog.from_demo()
        cls.log.fluence.gamma.calc_map()


class TestMachineLogs(TestCase):
    _logs_dir = osp.abspath(osp.join(osp.dirname(__file__), '.', 'test_files', 'MLC logs'))
    logs_dir = osp.join(_logs_dir, 'mixed_types')
    logs_altdir = osp.join(_logs_dir, 'altdir')
    mix_type_dir = osp.join(_logs_dir, 'mixed_types')

    def test_loading(self):
        # test root level directory
        logs = MachineLogs(self.logs_dir, recursive=False)
        self.assertEqual(logs.num_logs, 3)
        # test recursive
        logs = MachineLogs(self.logs_dir)
        self.assertEqual(logs.num_logs, 3)
        # test using zip file
        zfile = osp.join(self._logs_dir, 'mixed_types.zip')
        logs = MachineLogs.from_zip(zfile)
        self.assertEqual(logs.num_logs, 3)

    def test_basic_parameters(self):
        # no real test other than to make sure it works
        logs = MachineLogs(self.logs_dir)
        logs.report_basic_parameters()

    def test_num_logs(self):
        logs = MachineLogs(self.logs_dir, recursive=False)
        self.assertEqual(logs.num_logs, 3)
        self.assertEqual(logs.num_tlogs, 2)
        self.assertEqual(logs.num_dlogs, 1)

        logs = MachineLogs(self.mix_type_dir)
        self.assertEqual(logs.num_dlogs, 1)
        self.assertEqual(logs.num_tlogs, 2)

    def test_empty_dir(self):
        empty_dir = osp.join(self._logs_dir, 'empty_dir')
        logs = MachineLogs(empty_dir)
        self.assertEqual(logs.num_logs, 0)

    def test_mixed_types(self):
        """test mixed directory (tlogs & dlogs)"""
        log_dir = osp.join(self._logs_dir, 'mixed_types')
        logs = MachineLogs(log_dir)
        self.assertEqual(logs.num_logs, 3)

    def test_dlog_matches_missing(self):
        """Test that Dlogs without a match are skipped."""
        log_dir = osp.join(self._logs_dir, 'some_matches_missing')
        logs = MachineLogs(log_dir)
        self.assertEqual(logs.num_logs, 1)

    def test_append(self):
        # append a directory
        logs = MachineLogs(self.logs_altdir)
        logs.append(self.logs_altdir)
        self.assertEqual(logs.num_logs, 8)
        # append a file string
        single_file = osp.join(self.logs_altdir, 'Anonymous_4DC Treatment_JST90_TX_20140712094246.bin')
        logs.append(single_file)
        # append a MachineLog
        single_log = load_log(single_file)
        logs.append(single_log)

        # try to append something that's not a Log
        log = None
        with self.assertRaises(TypeError):
            logs.append(log)

    def test_avg_gamma(self):
        logs = MachineLogs(self.logs_dir, recursive=False)
        gamma = logs.avg_gamma()
        self.assertAlmostEqual(gamma, 0, delta=0.002)

    def test_avg_gamma_pct(self):
        logs = MachineLogs(self.logs_dir, recursive=False)
        gamma = logs.avg_gamma_pct()
        self.assertAlmostEqual(gamma, 100, delta=0.01)

    def test_writing_to_csv(self):
        logs = MachineLogs(self.logs_dir, recursive=False)
        files = logs.to_csv()
        self.assertIsInstance(files, list)
        # clean up by deleting files
        for file in files:
            os.remove(file)
