"""Run through the Machine Log bank."""
from unittest import TestCase

from pylinac import MachineLog
from tests.utils import DataBankMixin


def run_log(path):
    """Function to pass to the process pool executor to process machine logs."""
    try:
        log = MachineLog(path)
        log.fluence.gamma.calc_map()
        if log.fluence.gamma.pass_prcnt < 90:
            raise Exception("Gamma pass % < 90")
        return 'Success'
    except Exception as e:
        return 'Failure: {} @ {}'.format(e, path)


class TestLogBank(DataBankMixin, TestCase):
    DATA_DIR = ['Machine logs']

    def file_should_be_processed(self, filepath):
        return filepath.endswith('.dlg') or filepath.endswith('.bin')

    def test_all(self):
        super().test_all(run_log)
