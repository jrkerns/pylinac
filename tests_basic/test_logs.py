import io
import os.path as osp
import os
from unittest import TestCase
import shutil
import tempfile

import numpy as np

from pylinac.log_analyzer import MachineLogs, TreatmentType, \
    anonymize, TrajectoryLog, Dynalog, load_log, DynalogMatchError, NotADynalogError, NotALogError
from tests_basic.utils import save_file, CloudFileMixin, get_file_from_cloud_test_repo, \
    get_folder_from_cloud_test_repo, FromDemoImageTesterMixin, FromURLTesterMixin

TEST_DIR = 'mlc_logs'
ANONYMOUS_SOURCE_FOLDER = get_folder_from_cloud_test_repo(['mlc_logs', '_anonbase'])
ANONYMOUS_DEST_FOLDER = get_folder_from_cloud_test_repo(['mlc_logs', 'anonymous'])


class TestAnonymizeFunction(TestCase):
    """Test the anonymization method."""

    def setUp(self):
        anon_source = get_folder_from_cloud_test_repo(['mlc_logs', '_anonbase'])
        anon_dest = get_folder_from_cloud_test_repo(['mlc_logs', 'anonymous'])
        # move over files from other directory, since the filenames get overridden
        for file in os.listdir(anon_source):
            basefile = osp.join(anon_source, file)
            destfile = osp.join(anon_dest, file)
            if not osp.isfile(destfile):
                shutil.copy(basefile, anon_dest)

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

    def test_invalid(self):
        invalid_path = r'nonexistant/path'
        with self.assertRaises(NotALogError):
            anonymize(invalid_path)


class TestPublishPDF(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tlog = TrajectoryLog.from_demo()
        cls.dlog = Dynalog.from_demo()

    def test_publish_pdf(self):
        # normal publish; shouldn't raise
        with tempfile.TemporaryFile() as t:
            self.dlog.publish_pdf(t)
        with tempfile.TemporaryFile() as t:
            self.tlog.publish_pdf(t)

    def test_publish_pdf_w_imaging_log(self):
        imaging_tlog = TrajectoryLog(get_file_from_cloud_test_repo(['mlc_logs', 'tlogs', 'imaging.bin']))
        with self.assertRaises(ValueError), tempfile.NamedTemporaryFile() as t:
            imaging_tlog.publish_pdf(t.name)

    def test_publish_pdf_w_metadata_and_notes(self):
        with tempfile.TemporaryFile() as t:
            self.dlog.publish_pdf(t, metadata={'unit': 'TB1'}, notes='extra string')

        with tempfile.TemporaryFile() as t:
            self.tlog.publish_pdf(t, notes=['stuff', 'to', 'list'])


class LogPlottingSavingMixin:
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

    def test_fluence_plotting(self):
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


class TestTrajectoryTreatmentTypes(TestCase):

    def test_imaging_log(self):
        tlog = TrajectoryLog(get_file_from_cloud_test_repo(['mlc_logs', 'tlogs', 'imaging.bin']))
        self.assertTrue(tlog.treatment_type, TreatmentType.IMAGING.value)

    def test_vmat_log(self):
        tlog = TrajectoryLog(get_file_from_cloud_test_repo(['mlc_logs', 'tlogs', 'vmat.bin']))
        self.assertTrue(tlog.treatment_type, TreatmentType.VMAT.value)

    def test_static_imrt_log(self):
        tlog = TrajectoryLog(get_file_from_cloud_test_repo(['mlc_logs', 'tlogs', 'static_imrt.bin']))
        self.assertTrue(tlog.treatment_type, TreatmentType.STATIC_IMRT.value)

    def test_dynamic_imrt_log(self):
        tlog = TrajectoryLog(get_file_from_cloud_test_repo(['mlc_logs', 'tlogs', 'dynamic_imrt.bin']))
        self.assertTrue(tlog.treatment_type, TreatmentType.DYNAMIC_IMRT.value)


class TestDynalogTreatmentTypes(TestCase):

    def test_vmat_log(self):
        get_folder_from_cloud_test_repo(['mlc_logs', 'dlogs'])
        dlog = Dynalog(get_file_from_cloud_test_repo(['mlc_logs', 'dlogs', 'A_vmat.dlg']))
        self.assertTrue(dlog.treatment_type, TreatmentType.VMAT)

    def test_static_imrt_log(self):
        get_folder_from_cloud_test_repo(['mlc_logs', 'dlogs'])
        dlog = Dynalog(get_file_from_cloud_test_repo(['mlc_logs', 'dlogs', 'A_static_imrt.dlg']))
        self.assertTrue(dlog.treatment_type, TreatmentType.STATIC_IMRT)

    def test_dynamic_imrt_log(self):
        pass  # need to find one


class TestLoadLog(TestCase):

    def test_load_trajectory_log_from_file_object(self):
        path = get_file_from_cloud_test_repo(['mlc_logs', 'tlogs', 'dynamic_imrt.bin'])
        ref_log = TrajectoryLog(path)
        with open(path, 'rb') as f:
            t = TrajectoryLog(f)
        self.assertIsInstance(t, TrajectoryLog)
        self.assertEqual(t.num_beamholds, ref_log.num_beamholds)

    def test_dynalog_file(self):
        dynalog = get_file_from_cloud_test_repo(['mlc_logs', 'dlogs', 'A_static_imrt.dlg'])
        self.assertIsInstance(load_log(dynalog), Dynalog)

    def test_tlog_file(self):
        tlog = get_file_from_cloud_test_repo(['mlc_logs', 'tlogs', 'dynamic_imrt.bin'])
        self.assertIsInstance(load_log(tlog), TrajectoryLog)

    def test_url(self):
        url = r'https://s3.amazonaws.com/pylinac/Tlog.bin'
        self.assertIsInstance(load_log(url), TrajectoryLog)

    def test_dir(self):
        dlog_dir = get_folder_from_cloud_test_repo(['mlc_logs', 'dlogs'])
        self.assertIsInstance(load_log(dlog_dir), MachineLogs)

    def test_zip(self):
        zip_file = get_file_from_cloud_test_repo(['mlc_logs', 'mixed_types.zip'])
        self.assertIsInstance(load_log(zip_file), MachineLogs)

    def test_invalid_file(self):
        invalid_file = get_file_from_cloud_test_repo(['mlc_logs', 'Demo-subbeam-0-actual-fluence.npy'])
        with self.assertRaises(NotALogError):
            load_log(invalid_file)

    def test_invalid_path(self):
        invalid_path = r'nonexistant/path'
        with self.assertRaises(NotALogError):
            load_log(invalid_path)


class LogBase:
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


class TestTrajectoryLog(LogPlottingSavingMixin, LogBase, TestCase, FromDemoImageTesterMixin, FromURLTesterMixin):
    klass = TrajectoryLog
    demo_load_method = 'from_demo'
    url = 'Tlog.bin'
    anon_file = 'PatientID_4DC Treatment_JST90_TX_20140712094246.bin'

    def test_not_logs(self):
        # throw an error for files that aren't logs
        test_tlog = get_file_from_cloud_test_repo(['mlc_logs', 'tlogs', "Anonymous_4DC_Treatment_JS0_TX_20140712095629.bin"])
        not_a_file = test_tlog.replace(".bin", 'blahblah.bin')
        self.assertRaises(IOError, TrajectoryLog, not_a_file)
        not_a_log = get_file_from_cloud_test_repo(['VMAT', 'DRGSdmlc-105-example.dcm'])
        self.assertRaises(IOError, TrajectoryLog, not_a_log)

    def test_save_to_csv(self):
        save_file(self.log.to_csv)

    def test_txt_file_also_loads_if_around(self):
        # has a .txt file
        _ = get_folder_from_cloud_test_repo(['mlc_logs', 'mixed_types'])
        log_with_txt = get_file_from_cloud_test_repo(['mlc_logs', 'mixed_types', "Anonymous_4DC Treatment_JST90_TX_20140712094246.bin"])

        log = TrajectoryLog(log_with_txt)
        self.assertIsNotNone(log.txt)
        self.assertIsInstance(log.txt, dict)
        self.assertEqual(log.txt['Patient ID'], 'Anonymous')

        # DOESN'T have a txt file
        _ = get_folder_from_cloud_test_repo(['mlc_logs', 'tlogs'])
        log_no_txt = get_file_from_cloud_test_repo(['mlc_logs', 'tlogs', "Anonymous_4DC_Treatment_JS0_TX_20140712095629.bin"])

        log = TrajectoryLog(log_no_txt)
        self.assertIsNone(log.txt)


class TestDynalog(LogPlottingSavingMixin, LogBase, TestCase, FromDemoImageTesterMixin):
    klass = Dynalog
    demo_load_method = 'from_demo'
    anon_file = 'A1234_patientid.dlg'

    def test_loading_can_find_paired_file(self):
        # get all the test files
        get_folder_from_cloud_test_repo(['mlc_logs', 'dlogs'])

        # shouldn't raise since it can find B-file
        a_file = get_file_from_cloud_test_repo(['mlc_logs', 'dlogs', 'Adlog1.dlg'])
        Dynalog(a_file)

        # ditto for A-file
        b_file = get_file_from_cloud_test_repo(['mlc_logs', 'dlogs', 'Bdlog1.dlg'])
        Dynalog(b_file)

    def test_loading_bad_names(self):
        a_but_not_b_dir = get_file_from_cloud_test_repo(['mlc_logs', 'a_no_b_dir', 'Adlog1.dlg'])
        self.assertRaises(DynalogMatchError, Dynalog, a_but_not_b_dir)

        b_but_not_a_dir = get_file_from_cloud_test_repo(['mlc_logs', 'b_no_a_dir', 'Bdlog1.dlg'])
        self.assertRaises(DynalogMatchError, Dynalog, b_but_not_a_dir)

        bad_name_dlg = get_file_from_cloud_test_repo(['mlc_logs', 'bad_names', 'bad_name_dlg.dlg'])
        self.assertRaises(ValueError, Dynalog, bad_name_dlg)


class IndividualLogBase(CloudFileMixin):
    """Mixin to use when testing a single machine log; must be mixed with unittest.TestCase."""
    num_mlc_leaves = 120
    num_snapshots = 0
    num_beamholds = 0
    num_moving_leaves = 0
    treatment_type = ''
    dir_path = ['mlc_logs']
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
        if cls.log.treatment_type != TreatmentType.IMAGING.value:
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
        if self.log.treatment_type != TreatmentType.IMAGING.value:
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

    def test_publish_pdf(self):
        with io.BytesIO() as temp:
            self.log.publish_pdf(temp)


class IndividualTrajectoryLog(IndividualLogBase):
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

    def test_subbeam_fluences_unequal_to_cumulative(self):
        # as raised in #154
        if self.num_subbeams > 1:
            cumulative_fluence = self.log.fluence.actual.calc_map()
            subbeam_fluences = [subbeam.fluence.actual.calc_map() for subbeam in self.log.subbeams]
            if len(self.log.subbeams) > 0:
                for subbeam_fluence in subbeam_fluences:
                    self.assertFalse(np.array_equal(subbeam_fluence, cumulative_fluence))

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


class TestTrajectoryLogV4(IndividualTrajectoryLog, TestCase):
    version = 4.0
    dir_path = ['mlc_logs', 'tlogs']
    file_name = 'v4_log.bin'
    header = 'VOSTL'
    header_size = 1024
    sampling_interval = 20
    num_axes = 16
    mu_delivered = 100
    num_snapshots = 506
    axis_scale = 1
    num_subbeams = 1
    treatment_type = TreatmentType.STATIC_IMRT.value
    is_truncated = 0
    mlc_model = 2
    first_subbeam_data = {'gantry_angle': 180, 'collimator_angle': 270, 'jaw_x1': 10, 'jaw_x2': 10, 'jaw_y1': 10, 'jaw_y2': 10}
    plan_name = '4DC Treatment'

    def test_metadata(self):
        self.assertEqual(self.log.header.metadata.plan_name, self.plan_name)


class IndividualDynalog(IndividualLogBase):
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


class TestDynalogDemo(IndividualDynalog, TestCase):
    """Tests of the dynalog demo."""
    treatment_type = TreatmentType.DYNAMIC_IMRT.value
    num_beamholds = 20
    num_snapshots = 99
    average_rms = 0.04
    maximum_rms = 0.07
    average_gamma = 0.47
    percent_pass_gamma = 91
    leaf_move_status = {'moving': (9, 3), 'static': (8, )}
    delete_file = False

    @classmethod
    def setUpClass(cls):
        cls.log = Dynalog.from_demo()
        cls.log.fluence.gamma.calc_map()

    def test_fluences(self):
        reference_fluence = np.load(get_file_from_cloud_test_repo(['mlc_logs', 'Dynalog-demo-actual-fluence.npy']))
        self.log.fluence.actual.calc_map()
        demo_fluence = self.log.fluence.actual.array
        self.assertTrue(np.array_equal(demo_fluence, reference_fluence))


class TestTrajectoryLogDemo(IndividualTrajectoryLog, TestCase):
    """Tests for the demo trajectory log."""
    num_snapshots = 5200  # excluded: 1021
    num_subbeams = 2
    num_beamholds = 19
    mlc_model = 3
    treatment_type = TreatmentType.DYNAMIC_IMRT.value
    static_axes = ['collimator']
    moving_axes = ['gantry']
    average_rms = 0.001
    maximum_rms = 0.002
    percent_pass_gamma = 100
    mu_delivered = 183
    first_subbeam_data = {'gantry_angle': 310, 'collimator_angle': 180, 'jaw_x1': 3.7, 'jaw_x2': 3.4, 'jaw_y1': 3.8,
                          'jaw_y2': 3.9}
    delete_file = False

    @classmethod
    def setUpClass(cls):
        cls.log = TrajectoryLog.from_demo()
        cls.log.fluence.gamma.calc_map()

    def test_subbeam_fluences(self):
        # subbeam 0
        reference_fluence_0 = np.load(get_file_from_cloud_test_repo(['mlc_logs', 'Demo-subbeam-0-actual-fluence.npy']))
        self.log.subbeams[0].fluence.actual.calc_map()
        demo_fluence_0 = self.log.subbeams[0].fluence.actual.array
        self.assertTrue(np.array_equal(demo_fluence_0, reference_fluence_0))

        # subbeam 1
        reference_fluence_1 = np.load(get_file_from_cloud_test_repo(['mlc_logs', 'Demo-subbeam-1-actual-fluence.npy']))
        self.log.subbeams[1].fluence.actual.calc_map()
        demo_fluence_1 = self.log.subbeams[1].fluence.actual.array
        self.assertTrue(np.array_equal(demo_fluence_1, reference_fluence_1))

    def test_calc_gamma_early_fails(self):
        log = TrajectoryLog.from_demo()
        with self.assertRaises(ValueError):
            log.fluence.gamma.plot_map()


class TestMachineLogs(TestCase):

    @property
    def logs_dir(self):
        return get_folder_from_cloud_test_repo(['mlc_logs', 'mixed_types'])

    def test_loading(self):
        # test root level directory
        logs = MachineLogs(self.logs_dir, recursive=False)
        self.assertEqual(logs.num_logs, 3)
        # test recursive
        logs = MachineLogs(self.logs_dir)
        self.assertEqual(logs.num_logs, 3)
        # test using zip file
        zfile = get_file_from_cloud_test_repo(['mlc_logs', 'mixed_types.zip'])
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

    def test_empty_dir(self):
        empty_dir = get_folder_from_cloud_test_repo(['mlc_logs', 'empty_dir'])
        logs = MachineLogs(empty_dir)
        self.assertEqual(logs.num_logs, 0)
        with self.assertRaises(ValueError):
            logs.avg_gamma()

    def test_mixed_types(self):
        """test mixed directory (tlogs & dlogs)"""
        log_dir = get_folder_from_cloud_test_repo(['mlc_logs', 'mixed_types'])
        logs = MachineLogs(log_dir)
        self.assertEqual(logs.num_logs, 3)

    def test_dlog_matches_missing(self):
        """Test that Dlogs without a match are skipped."""
        log_dir = get_folder_from_cloud_test_repo(['mlc_logs', 'some_matches_missing'])
        logs = MachineLogs(log_dir)
        self.assertEqual(logs.num_logs, 1)

    def test_append(self):
        # append a directory
        logs = MachineLogs(get_folder_from_cloud_test_repo(['mlc_logs', 'altdir']))
        logs.append(get_folder_from_cloud_test_repo(['mlc_logs', 'altdir']))
        self.assertEqual(logs.num_logs, 8)
        # append a file string
        single_file = get_file_from_cloud_test_repo(['mlc_logs', 'altdir', 'Anonymous_4DC Treatment_JST90_TX_20140712094246.bin'])
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

    def test_writing_csv_with_no_logs(self):
        empty_dir = get_folder_from_cloud_test_repo(['mlc_logs', 'empty_dir'])
        logs = MachineLogs(empty_dir)
        logs.to_csv()  # shouldn't raise but will print a statement

    def test_anonymize(self):
        logs = MachineLogs(self.logs_dir, recursive=False)
        files = logs.anonymize(inplace=False, suffix='_suffixed')
        self.assertIsInstance(files, list)
        # cleanup
        for pdir, sdir, files in os.walk(self.logs_dir):
            to_remove = [file for file in files if 'suffixed' in file]
            for file in to_remove:
                os.remove(osp.join(pdir, file))
