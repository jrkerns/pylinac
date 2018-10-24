"""Run through the Winston-Lutz image bank."""
from unittest import TestCase

from pylinac import WinstonLutz
from tests_basic.utils import DataBankMixin


def run_wl(path):
    """Function to pass to the process pool executor to process picket fence images."""
    try:
        my_wl = WinstonLutz.from_zip(path)
        if my_wl.gantry_iso_size > 3:
            raise ValueError
        return 'Success'
    except Exception as e:
        return 'Failure: {} @ {}'.format(e, path)


class WinstonLutzTestBank(DataBankMixin, TestCase):
    DATA_DIR = ['Winston-Lutz']

    def file_should_be_processed(self, filepath):
        return filepath.endswith('.zip')

    def test_all(self):
        super().test_all(run_wl)
