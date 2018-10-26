import os.path as osp
from unittest import TestCase

from tests_basic import TEST_BANK_DIR
from tests_basic.test_planar_imaging import LeedsTORTestMixin, SIQC3TestMixin, LasVegasTestMixin


# mixins
class LeedsBankBase(LeedsTORTestMixin):
    dir_location = osp.join(TEST_BANK_DIR, '2D Image quality phantoms', 'Leeds')


class QC3BankBase(SIQC3TestMixin):
    dir_location = osp.join(TEST_BANK_DIR, '2D Image quality phantoms', 'QC-3')


class LasVegasBankBase(LasVegasTestMixin):
    dir_location = osp.join(TEST_BANK_DIR, '2D Image quality phantoms', 'Las Vegas')


# actual test cases
class LeedsACB1(LeedsBankBase, TestCase):
    file_path = ['ACB 1', '1.dcm']


class LasVegas10deg(LasVegasBankBase, TestCase):
    file_path = ['TrueBeam 1 - 2364', '2.5MV LV HQ 10deg - ImageRT_2016-10-6 20-12-58.dcm']
    phantom_angle = 290


class LasVegasrotated(LasVegasBankBase, TestCase):
    file_path = ['TrueBeam 1 - 2364', '2.5MV LV HQ side2 - ImageRT_2016-10-6 20-43-3.dcm']
    phantom_angle = 284


class LasVegasTB1(LasVegasBankBase, TestCase):
    file_path = ['TrueBeam 1 - 2364', '6MV LasVegas HQ 0deg - ImageRT_2016-10-6 20-10-17.dcm']
    phantom_angle = 284.5
