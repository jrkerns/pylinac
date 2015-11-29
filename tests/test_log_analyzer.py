import os.path as osp
import os
from unittest import TestCase
import shutil

import matplotlib.pyplot as plt

from pylinac.log_analyzer import MachineLog, MachineLogs, DYNALOG, TRAJECTORY_LOG, STATIC_IMRT, DYNAMIC_IMRT, VMAT
from tests.utils import save_file

test_dir = osp.join(osp.dirname(__file__), 'test_files', 'MLC logs')


class TestAnonymize(TestCase):
    """Test the anonymization method."""
    anon_originals_folder = osp.join(test_dir, '_anonbase')
    anon_folder = osp.join(test_dir, 'anonymous')
    files = []

    def setUp(self):
        # move over files from other directory, since the filenames get overridden
        for file in os.listdir(self.anon_originals_folder):
            basefile = osp.join(self.anon_originals_folder, file)
            destfile = osp.join(self.anon_folder, file)
            if not osp.isfile(destfile):
                shutil.copy(basefile, self.anon_folder)

    @classmethod
    def tearDownClass(cls):
        # remove files from anonymous folder
        files = os.listdir(cls.anon_folder)
        files.remove('dummy.txt')
        for file in files:
            file = osp.join(cls.anon_folder, file)
            os.remove(file)
        plt.close('all')

    def test_dynalog(self):
        # test making an anonymized copy
        dlog_file = osp.join(self.anon_folder, 'A1234_patientid.dlg')
        dlog = MachineLog(dlog_file)
        dlog.anonymize()

        # test doing inplace anonymization
        files = dlog.anonymize(inplace=True, suffix='inplace')
        self.assertIsInstance(files, list)
        for file in files:
            self.assertTrue('inplace' in file)

    def test_trajectorylog(self):
        tlog_file = osp.join(self.anon_folder, 'PatientID_4DC Treatment_JST90_TX_20140712094246.bin')
        tlog = MachineLog(tlog_file)
        tlog.anonymize()

        files = tlog.anonymize(inplace=True, suffix='inplace')
        self.assertIsInstance(files, list)
        for file in files:
            self.assertTrue('inplace' in file)

    def test_destination(self):
        tlog_file = osp.join(self.anon_folder, 'PatientID_4DC Treatment_JST90_TX_20140712094246.bin')
        tlog = MachineLog(tlog_file)
        tlog.anonymize(destination=self.anon_folder)  # shouldn't raise

    def test_from_stream(self):
        """Anonymizing a log that was loaded from a stream should fail (no filename to replace)."""
        url = 'https://s3.amazonaws.com/assuranceqa-staging/uploads/imgs/Tlog2.bin'
        tlog = MachineLog.from_url(url)
        with self.assertRaises(IOError):
            tlog.anonymize()

    def test_bad_name(self):
        """Test that a log with a bad name (no underscore) fails gracefully."""
        dlog_file = osp.join(self.anon_folder, 'Adlog1.dlg')
        dlog = MachineLog(dlog_file)
        with self.assertRaises(NameError):
            dlog.anonymize()

    def test_machinelogs(self):
        """Test anonymization of multiple machine logs."""
        logs = MachineLogs()
        logs.append(osp.join(self.anon_folder, 'A1234_patientid.dlg'))
        logs.append(osp.join(self.anon_folder, 'PatientID_4DC Treatment_JST90_TX_20140712094246.bin'))

        logs.anonymize()  # shouldn't raise


class TestLogLoading(TestCase):
    """Tests of dynalog files, mostly using the demo file."""
    test_dir = osp.join(osp.dirname(__file__), 'test_files', 'MLC logs')

    def test_loading(self):
        """Test that loading the badly-named dynalog still loads and identifies properly."""
        test_tlog = osp.join(self.test_dir, 'tlogs', "Anonymous_4DC Treatment_JS0_TX_20140712095629.bin")
        # should produce no errors
        # load method 1
        MachineLog(test_tlog)
        # load method 2
        log = MachineLog()
        self.assertFalse(log.is_loaded)
        log.load(test_tlog)
        self.assertTrue(log.is_loaded)

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

    def test_txt_file_also_loads_if_around(self):
        # has a .txt file
        log_with_txt = osp.join(self.test_dir, 'mixed_types', "Anonymous_4DC Treatment_JST90_TX_20140712094246.bin")

        log = MachineLog(log_with_txt)
        self.assertTrue(hasattr(log, 'txt'))
        self.assertIsInstance(log.txt, dict)
        self.assertEqual(log.txt['Patient ID'], 'Anonymous')

        # DOESN'T have a txt file
        log_no_txt = osp.join(self.test_dir, 'tlogs', "Anonymous_4DC Treatment_JS0_TX_20140712095629.bin")

        log = MachineLog(log_no_txt)
        self.assertFalse(hasattr(log, 'txt'))

    def test_from_url(self):
        url = 'https://s3.amazonaws.com/assuranceqa-staging/uploads/imgs/Tlog2.bin'
        MachineLog.from_url(url)  # shouldn't raise


class TestLogPlottingSaving(TestCase):
    """Test the plotting methods and plot saving methods."""
    tlog = MachineLog.from_demo_trajectorylog()
    dlog = MachineLog.from_demo_dynalog()

    def test_save_tlog_to_csv(self):
        save_file(self.tlog.to_csv)

    def test_plot_axes(self):
        for methodname in ('plot_actual', 'plot_expected', 'plot_difference'):
            method = getattr(self.tlog.axis_data.gantry, methodname)
            method()  # shouldn't raise

    def test_save_axes(self):
        for methodname in ('save_plot_actual', 'save_plot_expected', 'save_plot_difference'):
            # save matplotlib figures
            method = getattr(self.tlog.axis_data.gantry, methodname)
            save_file(method)

            # save MPLD3 HTML
            save_file(method, interactive=True, as_file_object='str')

    def test_fluence_plotting(self):
        # raise error if map hasn't yet been calc'ed.
        with self.assertRaises(AttributeError):
            self.dlog.fluence.actual.plot_map()

        self.dlog.fluence.actual.calc_map()
        self.dlog.fluence.actual.plot_map()

    def test_saving_fluence_plots(self):
        self.dlog.fluence.gamma.calc_map()
        save_file(self.dlog.fluence.gamma.save_map)

    def test_save_summary(self):
        self.tlog.fluence.gamma.calc_map()
        save_file(self.tlog.save_summary)


class TestLogMixin:
    """Mixin to use when testing a single machine log; must be mixed with unittest.TestCase."""
    log_path = ''
    log_type = TRAJECTORY_LOG
    num_mlc_leaves = 120
    num_snapshots = 0
    num_beamholds = 0
    num_moving_leaves = 0
    treatment_type = STATIC_IMRT
    static_axes = []
    moving_axes = []
    leaf_move_status = {'moving': tuple(), 'static': tuple()}
    average_rms = 0
    maximum_rms = 0
    average_gamma = 0
    percent_pass_gamma = 100
    mu_delivered = 25000  # always 25,000 for dynalogs, variable for trajectory logs
    version = 2.1  # or 3.0, or 'B' if dynalog

    # Trajectory log-specific data
    header = 'VOSTL'
    header_size = 1024
    sampling_interval = 20
    num_axes = 14
    axis_scale = 1
    num_subbeams = 0
    is_truncated = 0
    mlc_model = 3
    first_subbeam_data = {'gantry_angle': 0, 'collimator_angle': 0, 'jaw_x1': 0, 'jaw_x2': 0, 'jaw_y1': 0, 'jaw_y2': 0}

    # Dynalog-specific data
    tolerance = 102
    clinac_scale = 1

    @classmethod
    def setUpClass(cls):
        cls.log = MachineLog(cls.log_path)
        cls.log.fluence.gamma.calc_map()

    def test_log_type(self):
        """Test all kinds of things about the dynalog demo."""
        self.assertEqual(self.log.log_type, self.log_type)

    def test_num_leaves(self):
        """Test the number of MLC leaves and pairs."""
        self.assertEqual(self.log.header.num_mlc_leaves, self.num_mlc_leaves)

    def test_treatment_type(self):
        """Test the treatment type."""
        self.assertEqual(self.log.treatment_type, self.treatment_type)

    def test_num_snapshots(self):
        """Test the number of snapshots in the log."""
        if self.log.log_type == DYNALOG:
            self.assertEqual(self.log.axis_data.num_snapshots, self.num_snapshots)
        else:
            self.assertEqual(self.log.header.num_snapshots, self.num_snapshots)

    def test_num_beamholds(self):
        """Test the number of times the beam was held in the log."""
        self.assertEqual(self.log.axis_data.num_beamholds, self.num_beamholds)

    def test_rms_error(self):
        """Test the average and maximum RMS errors."""
        self.assertAlmostEqual(self.log.axis_data.mlc.get_RMS_avg(), self.average_rms, delta=0.01)
        self.assertAlmostEqual(self.log.axis_data.mlc.get_RMS_max(), self.maximum_rms, delta=0.01)

    def test_fluence_gamma(self):
        """Test gamma results for fluences."""
        self.assertAlmostEqual(self.log.fluence.gamma.avg_gamma, self.average_gamma, delta=0.01)
        self.assertAlmostEqual(self.log.fluence.gamma.pass_prcnt, self.percent_pass_gamma, delta=0.1)

    def test_header(self):
        """Test a few header values; depends on log type."""
        header = self.log.header
        self.assertEqual(header.version, self.version)
        if self.log_type == DYNALOG:
            self.assertEqual(header.tolerance, self.tolerance)
            self.assertEqual(header.clinac_scale, self.clinac_scale)
        else:
            self.assertEqual(header.header, self.header)
            self.assertEqual(header.header_size, self.header_size)
            self.assertEqual(header.sampling_interval, self.sampling_interval)
            self.assertEqual(header.num_axes, self.num_axes)
            self.assertEqual(header.axis_scale, self.axis_scale)
            self.assertEqual(header.num_subbeams, self.num_subbeams)
            self.assertEqual(header.is_truncated, self.is_truncated)
            self.assertEqual(header.mlc_model, self.mlc_model)

    def test_mu_delivered(self):
        """Test the number of MU delivered during the log."""
        self.assertAlmostEqual(self.log.axis_data.mu.actual[-1], self.mu_delivered, delta=0.01)

    def test_static_axes(self):
        """Test that certain axes did not move during treatment."""
        for axis_name in self.static_axes:
            axis = getattr(self.log.axis_data, axis_name)
            self.assertFalse(axis.moved)

    def test_subbeam_data(self):
        """Test the first subbeam data."""
        if self.log_type == TRAJECTORY_LOG:
            first_subbeam = self.log.subbeams[0]
            for key, known_value in self.first_subbeam_data.items():
                axis = getattr(first_subbeam, key)
                self.assertAlmostEqual(known_value, axis.actual, delta=0.1)

    def test_leaf_moved_status(self):
        """Test that the given leaves either moved or did not move."""
        moving_leaves = self.leaf_move_status['moving']
        for leaf in moving_leaves:
            self.assertTrue(self.log.axis_data.mlc.leaf_moved(leaf))

        static_leaves = self.leaf_move_status['static']
        for leaf in static_leaves:
            self.assertFalse(self.log.axis_data.mlc.leaf_moved(leaf))


class DynalogDemo(TestLogMixin, TestCase):
    """Tests of the dynalog demo."""
    log_type = DYNALOG
    treatment_type = DYNAMIC_IMRT
    num_beamholds = 20
    num_snapshots = 99
    average_rms = 0.037
    maximum_rms = 0.076
    average_gamma = 0.019
    percent_pass_gamma = 99.85
    version = 'B'
    tolerance = 102
    clinac_scale = 1
    leaf_move_status = {'moving': (9, 3), 'static': (8, )}

    @classmethod
    def setUpClass(cls):
        cls.log = MachineLog.from_demo_dynalog()
        cls.log.fluence.gamma.calc_map()

    def test_demo(self):
        """Test the run demo method."""
        # shouldn't raise
        MachineLog().run_dlog_demo()


class TrajectoryLogDemo(TestLogMixin, TestCase):
    """Tests for the demo trajectory log."""
    log_type = TRAJECTORY_LOG
    num_snapshots = 5200  # excluded: 1021
    num_subbeams = 2
    version = 2.1
    static_axes = ['collimator']
    moving_axes = ['gantry']
    average_rms = 0.001
    maximum_rms = 0.002
    average_gamma = 0.001
    percent_pass_gamma = 100
    mu_delivered = 183
    first_subbeam_data = {'gantry_angle': 310, 'collimator_angle': 180, 'jaw_x1': 3.7, 'jaw_x2': 3.4, 'jaw_y1': 3.8,
                          'jaw_y2': 3.9}

    @classmethod
    def setUpClass(cls):
        cls.log = MachineLog.from_demo_trajectorylog()
        cls.log.fluence.gamma.calc_map()

    def test_demo(self):
        """Test the run demo method."""
        # shouldn't raise
        MachineLog().run_tlog_demo()


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
        # test using method
        logs = MachineLogs()
        logs.load_folder(self.logs_dir)
        self.assertEqual(logs.num_logs, 3)
        # test using zip file
        zfile = osp.join(self._logs_dir, 'mixed_types.zip')
        logs = MachineLogs.from_zip(zfile)
        self.assertEqual(logs.num_logs, 4)

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
        logs = MachineLogs()
        logs.append(self.logs_altdir)
        self.assertEqual(logs.num_logs, 4)
        # append a file string
        logs = MachineLogs()
        single_file = osp.join(self.logs_altdir, 'Anonymous_4DC Treatment_JST90_TX_20140712094246.bin')
        logs.append(single_file)
        # append a MachineLog
        logs = MachineLogs()
        single_log = MachineLog(single_file)
        logs.append(single_log)

        # try to append something that's not a Log
        log = None
        with self.assertRaises(TypeError):
            logs.append(log)

    def test_empty_op(self):
        """Test that error is raised if trying to do op with no logs."""
        logs = MachineLogs()
        self.assertRaises(ValueError, logs.avg_gamma)

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
