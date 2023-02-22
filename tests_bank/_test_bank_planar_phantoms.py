"""Run through the 2D planar phantom image bank."""
from unittest import TestCase

# from tests_basic import prep_mpl_testing
import matplotlib.pyplot as plt

from pylinac import LasVegas, LeedsTOR, StandardImagingQC3
from pylinac.planar_imaging import PTWEPIDQC, SNCMV, SNCkV, StandardImagingQCkV
from tests_basic.utils import DataBankMixin


def process_phantom(phantom, path):
    try:
        phantom.analyze()
        phantom.plot_analyzed_image(show=False)
        plt.close("all")
    except Exception as e:
        return f"Failure: {e} @ {path}"
    else:
        return "Success"


def run_leeds(path):
    """Function to pass to the process pool executor to process Leeds phantom images."""
    leeds = LeedsTOR(path)
    return process_phantom(leeds, path)


class TestLeedsImageBank(DataBankMixin, TestCase):
    DATA_DIR = ["2D Image quality phantoms", "Leeds"]

    def test_all(self):
        super().test_all(run_leeds)


def run_qc3(path):
    """Function to pass to the process pool executor to process PipsPro-QC3 images."""
    pp = StandardImagingQC3(path)
    return process_phantom(pp, path)


class TestQC3ImageBank(DataBankMixin, TestCase):
    DATA_DIR = ["2D Image quality phantoms", "SI QC-3"]
    print_success_path = True
    write_failures_to_file = True

    def test_all(self):
        super().test_all(run_qc3)


def run_si_qckv(path):
    """Function to pass to the process pool executor to process PipsPro-QC3 images."""
    pp = StandardImagingQCkV(path)
    return process_phantom(pp, path)


class TestQCkVImageBank(DataBankMixin, TestCase):
    DATA_DIR = ["2D Image quality phantoms", "SI QC-kV"]
    print_success_path = True
    write_failures_to_file = True

    def test_all(self):
        super().test_all(run_si_qckv)


def run_lasvegas(path):
    lv = LasVegas(path)
    return process_phantom(lv, path)


class TestLasVegasImageBank(DataBankMixin, TestCase):
    DATA_DIR = ["2D Image quality phantoms", "Las Vegas"]
    write_failures_to_file = True

    def test_all(self):
        super().test_all(run_lasvegas)


def run_ptwepid(path):
    lv = PTWEPIDQC(path)
    return process_phantom(lv, path)


class TestPTWEPIDBank(DataBankMixin, TestCase):
    DATA_DIR = ["2D Image quality phantoms", "PTW-EPID"]
    write_failures_to_file = True

    def test_all(self):
        super().test_all(run_ptwepid)


def run_snckv(path):
    lv = SNCkV(path)
    return process_phantom(lv, path)


class TestSNCkVBank(DataBankMixin, TestCase):
    DATA_DIR = ["2D Image quality phantoms", "Sun Nuclear", "kV"]
    write_failures_to_file = True

    def test_all(self):
        super().test_all(run_snckv)


def run_sncmv(path):
    phan = SNCMV(path)
    return process_phantom(phan, path)


class TestSNCMVBank(DataBankMixin, TestCase):
    DATA_DIR = ["2D Image quality phantoms", "Sun Nuclear", "MV"]
    write_failures_to_file = True

    def test_all(self):
        super().test_all(run_sncmv)
