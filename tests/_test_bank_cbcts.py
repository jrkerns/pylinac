"""Run through the CBCT image bank. Since the criteria and analysis are different than the default DataBankMixin, we must create our own."""
import os.path as osp
from unittest import TestCase
from zipfile import BadZipfile

from pylinac import CatPhan600, CatPhan503, CatPhan504
from tests.utils import DataBankMixin


def run_catphan504(path):
    """Function to pass to the process pool executor to process cbct images."""
    try:
        mypf = CatPhan504.from_zip(path)
        mypf.analyze()
        return 'Success'
    except (ValueError, FileNotFoundError, BadZipfile) as e:
        return 'Failure: {} @ {}'.format(e, path)


def run_catphan503(path):
    try:
        mypf = CatPhan503.from_zip(path)
        mypf.analyze()
        return 'Success'
    except (ValueError, FileNotFoundError, BadZipfile) as e:
        return 'Failure: {} @ {}'.format(e, path)


def run_catphan600(path):
    try:
        mypf = CatPhan600.from_zip(path)
        mypf.analyze()
        return 'Success'
    except (ValueError, FileNotFoundError, BadZipfile) as e:
        return 'Failure: {} @ {}'.format(e, path)


class CBCTTestBank(DataBankMixin, TestCase):
    DATA_DIR = ['CBCTs']
    print_success_path = False

    def file_should_be_processed(self, filepath):
        return filepath.endswith('.zip')

    def test_all(self):
        orig_dir = self.DATA_DIR
        self.DATA_DIR = osp.join(orig_dir, 'CatPhan 503')
        super().test_all(run_catphan503)
        self.DATA_DIR = osp.join(orig_dir, 'CatPhan 600')
        super().test_all(run_catphan600)
        self.DATA_DIR = osp.join(orig_dir, 'CatPhan 504')
        super().test_all(run_catphan504)
