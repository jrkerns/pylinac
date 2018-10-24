"""Run through the Starshot image bank."""
from unittest import TestCase

from pylinac import Starshot
from tests_basic.utils import DataBankMixin


def run_star(path):
    """Function to pass to the process pool executor to process starshot images."""
    try:
        mystar = Starshot(path, sid=1000)
        mystar.analyze()
        if mystar.wobble.diameter_mm > 3:
            raise Exception("Diamater was > 3mm.")
        return 'Success'
    except Exception as e:
        return 'Failure: {} @ {}'.format(e, path)


class StarshotBank(DataBankMixin, TestCase):
    DATA_DIR = ['Starshots']

    def test_all(self):
        super().test_all(run_star)
