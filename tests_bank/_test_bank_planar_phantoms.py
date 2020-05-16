"""Run through the 2D planar phantom image bank."""
from unittest import TestCase

# from tests_basic import prep_mpl_testing
import matplotlib.pyplot as plt

from pylinac import LeedsTOR, StandardImagingQC3, LasVegas
from tests_basic.utils import DataBankMixin


def process_phantom(phantom, path):
    try:
        phantom.analyze()
        phantom.plot_analyzed_image(show=False)
        plt.close('all')
    except Exception as e:
        return 'Failure: {} @ {}'.format(e, path)
    else:
        return 'Success'


def run_leeds(path):
    """Function to pass to the process pool executor to process Leeds phantom images."""
    leeds = LeedsTOR(path)
    return process_phantom(leeds, path)


def run_qc3(path):
    """Function to pass to the process pool executor to process PipsPro-QC3 images."""
    pp = StandardImagingQC3(path)
    return process_phantom(pp, path)


def run_lasvegas(path):
    lv = LasVegas(path)
    return process_phantom(lv, path)


class TestLeedsImageBank(DataBankMixin, TestCase):
    DATA_DIR = ['2D Image quality phantoms', 'Leeds']

    def test_all(self):
        super().test_all(run_leeds)


class TestQC3ImageBank(DataBankMixin, TestCase):
    DATA_DIR = ['2D Image quality phantoms', 'QC-3']
    print_success_path = True
    write_failures_to_file = True

    def test_all(self):
        super().test_all(run_qc3)


class TestLasVegasImageBank(DataBankMixin, TestCase):
    DATA_DIR = ['2D Image quality phantoms', 'Las Vegas']
    write_failures_to_file = True

    def test_all(self):
        super().test_all(run_lasvegas)
