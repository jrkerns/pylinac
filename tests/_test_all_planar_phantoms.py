import os.path as osp
from unittest import TestCase

from tests import TEST_BANK_DIR
from tests.test_planar_imaging import LeedsTORTestMixin, PipsProTestMixin


class LeedsBankBase(LeedsTORTestMixin):
    dir_location = osp.join(TEST_BANK_DIR, '2D Image quality phantoms', 'Leeds')


class PipsProBankBase(PipsProTestMixin):
    dir_location = osp.join(TEST_BANK_DIR, '2D Image quality phantoms', 'PipsPro')


class LeedsACB1(LeedsBankBase, TestCase):
    file_path = ['ACB 1', '1.dcm']
