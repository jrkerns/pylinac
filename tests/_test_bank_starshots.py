"""Run through the Starshot image bank."""
from unittest import TestCase

from pylinac import Starshot
from tests.utils import DataBankMixin


def run_star(path):
    """Function to pass to the process pool executor to process starshot images."""
    try:
        mystar = Starshot(path)
        mystar.analyze()
        return 'Success'
    except:
        return 'Failure at {}'.format(path)


class StarshotBank(DataBankMixin, TestCase):
    DATA_DIR = 'Starshots'

    def test_all(self):
        super().test_all(run_star)
