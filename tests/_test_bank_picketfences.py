"""Run through the Picket Fence image bank."""
from unittest import TestCase

from pylinac import PicketFence
from tests.utils import DataBankMixin


def run_pf(path):
    """Function to pass to the process pool executor to process picket fence images."""
    try:
        mypf = PicketFence(path)
        mypf.analyze()
        return 'Success'
    except ValueError:
        try:
            mypf = PicketFence(path, filter=3)
            mypf.analyze()
            return 'Success'
        except:
            return 'Failure at {}'.format(path)
    except:
        return 'Failure at {}'.format(path)


class PicketFenceTestBank(DataBankMixin, TestCase):
    DATA_DIR = ['Picket Fences']

    def file_should_be_processed(self, filepath):
        return filepath.endswith('.dcm')

    def test_all(self):
        super().test_all(run_pf)
