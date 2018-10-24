"""Run through the Machine Log bank."""
from unittest import TestCase

from pylinac import load_log
from pylinac.log_analyzer import IMAGING
from tests_basic.utils import DataBankMixin


def run_log(path):
    """Function to pass to the process pool executor to process machine logs."""
    try:
        log = load_log(path)
        if not log.treatment_type == IMAGING:
            log.fluence.gamma.calc_map()
            if log.fluence.gamma.pass_prcnt < 90:
                raise Exception("Gamma pass % < 90")
        ret = "Success"
    except Exception as e:
        ret = 'Failure: {} @ {}'.format(e, path)
    return ret


class TestLogBank(DataBankMixin, TestCase):
    DATA_DIR = ['Machine logs']
    print_success_path = True

    def file_should_be_processed(self, filepath):
        return filepath.endswith('.dlg') or filepath.endswith('.bin')

    def test_all(self):
        super().test_all(run_log)
