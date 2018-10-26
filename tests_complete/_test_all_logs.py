import os.path as osp
from unittest import TestCase

from pylinac.log_analyzer import VMAT, IMAGING, STATIC_IMRT, DYNAMIC_IMRT
from tests_basic.test_logs import TestIndividualLogBase
from tests_basic import TEST_BANK_DIR


class LogBankBase(TestIndividualLogBase):
    dir_location = osp.join(TEST_BANK_DIR, 'Machine logs')


class One(LogBankBase, TestCase):
    file_path = ['Anonymous_4DC Treatment_A_TX_20120928131920.bin']
    treatment_type = STATIC_IMRT
    num_subbeams = 1
    mu_delivered = 209
    num_snapshots = 1098
    first_subbeam_data = {'gantry_angle': 185, 'collimator_angle': 180, 'jaw_y1': 10.5}


class Two(LogBankBase, TestCase):
    file_path = ['Anonymous_4DC Treatment_A1_TX_20120928132942.bin']
    treatment_type = DYNAMIC_IMRT
    num_subbeams = 10
    mu_delivered = 681
    num_beamholds = 142
    num_snapshots = 28268
    first_subbeam_data = {'gantry_angle': 340, 'collimator_angle': 180, 'jaw_y1': 10.8}


class DynalogArc(LogBankBase, TestCase):
    file_path = ['Katy iX', 'A20120712122417_Anonymous.dlg']
    treatment_type = VMAT
    version = 'B'
    tolerance = 510
    average_gamma = 0.06
    mu_delivered = 2696
    num_snapshots = 1151
    average_rms = 0.16
    maximum_rms = 0.205


class Four(LogBankBase, TestCase):
    file_path = ['Chicago', 'T-Log HDMLC', 'anonymized_4DC Treatment_1.1_TX_20151015093202.bin']
    num_snapshots = 6356
    version = 3
    treatment_type = DYNAMIC_IMRT
    num_subbeams = 2
    num_axes = 16
    mu_delivered = 535
    num_beamholds = 2
    mlc_model = 3
    first_subbeam_data = {'gantry_angle': 178.9, 'jaw_x2': 5.2}


class CBCTSetup(LogBankBase, TestCase):
    file_path = ['Chicago', 'T-Log HDMLC', 'anonymized_4DC Treatment_KVCBCT_Setup_20151015090308.bin']
    num_snapshots = 1238
    version = 3
    treatment_type = IMAGING
    num_subbeams = 1
    num_axes = 16
    mu_delivered = 0
    num_beamholds = 0
    mlc_model = 3
    first_subbeam_data = {}


class KVSetup(LogBankBase, TestCase):
    file_path = ['Chicago', 'T-Log HDMLC', 'anonymized_4DC Treatment_KVKV_SetupPair_20151015130741.bin']
    num_snapshots = 185
    version = 3
    treatment_type = IMAGING
    mlc_model = 3
    num_subbeams = 1
    num_axes = 16
    mu_delivered = 0
    num_beamholds = 0
    first_subbeam_data = {}


class DoubleExposure(LogBankBase, TestCase):
    file_path = ['Chicago', 'T-Log HDMLC', 'anonymized_4DC Treatment_Planned_Double_Exposure_ADHOC_20151015140943.bin']
    num_snapshots = 750
    version = 3
    treatment_type = IMAGING
    mlc_model = 3
    average_gamma = 0
    num_subbeams = 2
    num_axes = 16
    mu_delivered = 2
    num_beamholds = 4
    first_subbeam_data = {'gantry_angle': 180, 'jaw_x2': 4}


class Five(LogBankBase, TestCase):
    file_path = ['Chicago', 'T-Log Mil120', 'anonymized_4DC Treatment_1.1_TX_20151015102651.bin']
    num_snapshots = 10728
    version = 3
    num_subbeams = 3
    treatment_type = DYNAMIC_IMRT
    num_axes = 16
    mu_delivered = 428
    num_beamholds = 3
    first_subbeam_data = {'gantry_angle': 176.7, 'jaw_x2': 8.2}


class OpenPort(LogBankBase, TestCase):
    file_path = ['Chicago', 'T-Log Mil120', 'anonymized_4DC Treatment_Planned_Open_Port_Image_ADHOC_20151015131101.bin']
    num_snapshots = 72
    version = 3
    treatment_type = IMAGING
    num_subbeams = 1
    num_axes = 16
    mu_delivered = 1
    num_beamholds = 3
    first_subbeam_data = {'gantry_angle': 180, 'jaw_x2': 6}


class Six(LogBankBase, TestCase):
    file_path = ['Bay Area iX', 'A20121212123129_Anonymous.dlg']
    treatment_type = VMAT
    version = 'B'
    tolerance = 510
    num_snapshots = 1150
    average_rms = 0.11
    maximum_rms = 0.14
    num_subbeams = 1
    num_axes = 16
    mu_delivered = 2696
    average_gamma = 0.03
