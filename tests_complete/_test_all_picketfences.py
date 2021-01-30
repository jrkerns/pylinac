"""Travis CI memory can't handle all the picketfences; thus only test them when explicitly asked to."""
import os.path as osp
from unittest import TestCase, skip

from tests_basic import TEST_BANK_DIR
from tests_basic.test_picketfence import PFTestMixin, PicketFence, LEFT_RIGHT


class PFBankMixin(PFTestMixin):
    """Base Picket Fence class for analyzing images in the bank directory."""
    dir_location = osp.join(TEST_BANK_DIR, 'Picket Fences')


class AS500(PFBankMixin, TestCase):
    """Tests for the AS500 image."""
    file_path = ['AS500_PF.dcm']
    max_error = 0.15
    abs_median_error = 0.04


class AS5002(PFBankMixin, TestCase):
    """Tests for the AS500#2 image."""
    file_path = ['AS500#2.dcm']
    max_error = 0.12
    abs_median_error = 0.03


class AS5003(PFBankMixin, TestCase):
    """Tests for the AS500#3 image."""
    file_path = ['AS500#3.dcm']
    max_error = 0.16
    abs_median_error = 0.03


class AS5004(PFBankMixin, TestCase):
    """Tests for the AS500#4 image."""
    file_path = ['AS500#4.dcm']
    max_error = 0.28
    abs_median_error = 0.06


class AS5005(PFBankMixin, TestCase):
    """Tests for the AS500#4 image."""
    file_path = ['AS500#5.dcm']
    max_error = 0.23
    abs_median_error = 0.04


class AS5006(PFBankMixin, TestCase):
    """Tests for the AS500#4 image."""
    file_path = ['AS500#6.dcm']
    picket_orientation = LEFT_RIGHT
    max_error = 0.23
    abs_median_error = 0.06


class AS5007(PFBankMixin, TestCase):
    """Tests for the AS500#4 image."""
    file_path = ['AS500#7.dcm']
    max_error = 0.24
    abs_median_error = 0.05


class AS5008(PFBankMixin, TestCase):
    """Tests for the AS500#4 image."""
    file_path = ['AS500#8.dcm']
    max_error = 0.2
    abs_median_error = 0.04


class AS5009(PFBankMixin, TestCase):
    """Tests for the AS500#4 image."""
    file_path = ['AS500#9.dcm']
    max_error = 0.24
    abs_median_error = 0.04


class AS50010(PFBankMixin, TestCase):
    """Tests for the AS500#4 image."""
    file_path = ['AS500#10.dcm']
    picket_orientation = LEFT_RIGHT
    max_error = 0.24
    abs_median_error = 0.05


class AS500error(PFBankMixin, TestCase):
    """Tests for the AS500#2 image."""
    file_path = ['AS500-error.dcm']
    num_pickets = 6
    percent_passing = 99
    max_error = 0.55
    abs_median_error = 0.07
    passes = False
    mean_picket_spacing = 20


class AS1000(PFBankMixin, TestCase):
    """Tests for the AS1000 image."""
    file_path = ['AS1000_PF.dcm']
    max_error = 0.29
    abs_median_error = 0.06


class AS1000_2(PFBankMixin, TestCase):
    """Tests for the AS1000 image."""
    file_path = ['AS1000#2.dcm']
    max_error = 0.24
    abs_median_error = 0.07


class AS1000_3(PFBankMixin, TestCase):
    """Tests for the AS1000 image."""
    file_path = ['AS1000#3.dcm']
    max_error = 0.13
    abs_median_error = 0.05


class AS1000_4(PFBankMixin, TestCase):
    """Tests for the AS1000 image."""
    file_path = ['AS1000#4.dcm']
    picket_orientation = LEFT_RIGHT
    max_error = 0.18
    abs_median_error = 0.05


class AS1000_90(PFBankMixin, TestCase):
    """Tests for the AS1000 image."""
    file_path = ['AS1000-90.dcm']
    picket_orientation = LEFT_RIGHT
    max_error = 0.23
    abs_median_error = 0.05


class AS1000HDSmall(PFBankMixin, TestCase):
    """Tests for the AS1000 image."""
    file_path = ['AS1000-HD-small.dcm']
    mlc = 'HD'
    max_error = 0.05
    abs_median_error = 0.05


class AS1000HDFull(PFBankMixin, TestCase):
    """Tests for the AS1000 image with a smaller pattern (only inner leaves)."""
    file_path = ['AS1000-HD-full.dcm']
    mlc = 'HD'
    max_error = 0.2
    abs_median_error = 0.06


class AS1000HDFullVMAT(PFBankMixin, TestCase):
    """Tests for the AS1000 image with a smaller pattern (only inner leaves)."""
    file_path = ['AS1000-HD-full-VMAT.dcm']
    mlc = 'HD'
    max_error = 0.2
    abs_median_error = 0.08


@skip  # says file isn't real DICOM TODO: Figure out why not real DICOM
class AS1000HDFullError(PFBankMixin, TestCase):
    """Tests for the AS1000 image with a few errors introduced."""
    file_path = ['AS1000-HD-full-error.dcm']
    mlc = 'HD'
    num_pickets = 6
    abs_median_error = 0.03
    max_error = 0.39

    def test_lower_tolerance_fails(self):
        """This image has an introduced error; this should catch with a reasonable tolerance."""
        pf = PicketFence(self.file_path)
        pf.analyze(tolerance=0.3, hdmlc=self.hdmlc)
        self.assertFalse(pf.passed)


class AS1200(PFBankMixin, TestCase):
    """Tests for the AS1200 image."""
    file_path = ['AS1200.dcm']
    max_error = 0.08
    abs_median_error = 0.02


class AS1200Error(PFBankMixin, TestCase):
    """Tests for the AS1200 image."""
    file_path = ['AS1200-error.dcm']
    num_pickets = 6
    max_error = 0.48
    abs_median_error = 0.05
    sag_adjustment = -1.2
    mean_picket_spacing = 20


class AS1200ExtendedSID(PFBankMixin, TestCase):
    """Tests for the AS1200 image."""
    file_path = ['AS1200-ExtendedSID.dcm']
    max_error = 0.12
    abs_median_error = 0.04


class AS1200ExtendedSIDVMAT(PFBankMixin, TestCase):
    """Tests for the AS1200 image."""
    file_path = ['AS1200-ExtendedSID-VMAT.dcm']
    max_error = 0.18
    abs_median_error = 0.06


class AS1200HD(PFBankMixin, TestCase):
    """Tests for the AS1200 image."""
    file_path = ['AS1200-HD.dcm']
    mlc = 'HD'
    max_error = 0.05
    abs_median_error = 0.02
    num_pickets = 10
    pass_num_pickets = True


class AS1200HDTranslated(PFBankMixin, TestCase):
    """Tests for the AS1200 image."""
    file_path = ['AS1200-HD-translated.dcm']
    mlc = 'HD'
    max_error = 0.15
    abs_median_error = 0.02
    num_pickets = 10
    pass_num_pickets = True


class ChicagoNoError(PFBankMixin, TestCase):
    file_path = ['Chicago', 'PF no error.dcm']
    # log = ['Chicago', 'PF no error tlog.bin']
    mlc = 'HD'
    max_error = 0.24


class ChicagoError(PFBankMixin, TestCase):
    file_path = ['Chicago', 'PF point2mm error.dcm']
    # log = ['Chicago', 'PF point2mm tlog.bin']
    mlc = 'HD'
    max_error = 0.2


@skip
class CharlestonRA(PFBankMixin, TestCase):
    file_path = ['Charleston', 'TB1', 'July2016', 'RA.dcm']
    max_error = 0.17


@skip
class CharlestonG0(PFBankMixin, TestCase):
    file_path = ['Charleston', 'TB1', 'July2016', 'G0.dcm']
    max_error = 0.1
