"""Run through the CBCT image bank. Since the criteria and analysis are different than the default DataBankMixin, we must create our own."""
from unittest import TestCase
from zipfile import BadZipfile

from pylinac import CBCT
from tests.utils import DataBankMixin


def run_cbct(path):
    """Function to pass to the process pool executor to process cbct images."""
    try:
        mypf = CBCT.from_zip(path)
        mypf.analyze()
        return 'Success'
    except (ValueError, FileNotFoundError, BadZipfile) as e:
        return 'Failure: {} @ {}'.format(e, path)


class CBCTTestBank(DataBankMixin, TestCase):
    DATA_DIR = ['CBCTs']
    print_success_path = True

    def file_should_be_processed(self, filepath):
        return filepath.endswith('.zip')

    def test_all(self):
        super().test_all(run_cbct)
