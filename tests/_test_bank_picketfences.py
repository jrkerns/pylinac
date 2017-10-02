"""Run through the Picket Fence image bank."""
from unittest import TestCase

from pylinac import PicketFence
from tests.utils import DataBankMixin


def run_pf(path):
    """Function to pass to the process pool executor to process picket fence images."""
    try:
        mypf = PicketFence(path)
        mypf.analyze()
        if mypf.max_error > 1.2:
            raise Exception("Max MLC peak error > 1.2mm")
        return 'Success'
    except ValueError:
        try:
            mypf = PicketFence(path, filter=3)
            mypf.analyze()
            if mypf.max_error > 1.2:
                raise Exception("Max MLC peak error > 1.2mm")
            return 'Success'
        except (ValueError,) as e:
            return 'Failure: {} @ {}'.format(e, path)
    except Exception as e:
        return 'Failure: {} @ {}'.format(e, path)


class PicketFenceTestBank(DataBankMixin, TestCase):
    DATA_DIR = ['Picket Fences']
    write_failures_to_file = True

    def file_should_be_processed(self, filepath):
        return filepath.endswith('.dcm')

    def test_all(self):
        super().test_all(run_pf)
