"""Run through the 2D planar phantom image bank."""
from unittest import TestCase

from pylinac import LeedsTOR, PipsProQC3
from tests.utils import DataBankMixin


def process_phantom(phantom, path):
    try:
        phantom.analyze()
        phantom.plot_analyzed_image()
    except:
        return 'Failure at {}'.format(path)
    else:
        return 'Success'


def run_leeds(path):
    """Function to pass to the process pool executor to process Leeds phantom images."""
    leeds = LeedsTOR(path)
    return process_phantom(leeds, path)


def run_pipspro(path):
    """Function to pass to the process pool executor to process PipsPro-QC3 images."""
    pp = PipsProQC3(path)
    return process_phantom(pp, path)


class TestLeedsImageBank(DataBankMixin, TestCase):
    DATA_DIR = '2D Image quality phantoms/Leeds'

    def test_all(self):
        super().test_all(run_leeds)


class TestPipsProImageBank(DataBankMixin, TestCase):
    DATA_DIR = '2D Image quality phantoms/PipsPro'

    def test_all(self):
        super().test_all(run_pipspro)
