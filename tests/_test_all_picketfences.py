"""Travis CI memory can't handle all the picketfences; thus only test them when explicitly asked to."""
import concurrent.futures
import os
import os.path as osp
import time
from unittest import TestCase, skip

from tests.test_picketfence import PFTestMixin, TEST_DIR, PicketFence, LEFT_RIGHT


class AS500(PFTestMixin, TestCase):
    """Tests for the AS500 image."""
    im_path = osp.join(TEST_DIR, 'AS500_PF.dcm')
    max_error = 0.15
    abs_median_error = 0.04


class AS5002(PFTestMixin, TestCase):
    """Tests for the AS500#2 image."""
    im_path = osp.join(TEST_DIR, 'AS500#2.dcm')
    max_error = 0.12
    abs_median_error = 0.03


class AS5003(PFTestMixin, TestCase):
    """Tests for the AS500#3 image."""
    im_path = osp.join(TEST_DIR, 'AS500#3.dcm')
    max_error = 0.16
    abs_median_error = 0.03


class AS5004(PFTestMixin, TestCase):
    """Tests for the AS500#4 image."""
    im_path = osp.join(TEST_DIR, 'AS500#4.dcm')
    max_error = 0.28
    abs_median_error = 0.06


class AS5005(PFTestMixin, TestCase):
    """Tests for the AS500#4 image."""
    im_path = osp.join(TEST_DIR, 'AS500#5.dcm')
    max_error = 0.23
    abs_median_error = 0.04


class AS5006(PFTestMixin, TestCase):
    """Tests for the AS500#4 image."""
    im_path = osp.join(TEST_DIR, 'AS500#6.dcm')
    picket_orientation = LEFT_RIGHT
    max_error = 0.23
    abs_median_error = 0.06


class AS5007(PFTestMixin, TestCase):
    """Tests for the AS500#4 image."""
    im_path = osp.join(TEST_DIR, 'AS500#7.dcm')
    max_error = 0.24
    abs_median_error = 0.05


class AS5008(PFTestMixin, TestCase):
    """Tests for the AS500#4 image."""
    im_path = osp.join(TEST_DIR, 'AS500#8.dcm')
    max_error = 0.2
    abs_median_error = 0.04


class AS5009(PFTestMixin, TestCase):
    """Tests for the AS500#4 image."""
    im_path = osp.join(TEST_DIR, 'AS500#9.dcm')
    max_error = 0.24
    abs_median_error = 0.04


class AS50010(PFTestMixin, TestCase):
    """Tests for the AS500#4 image."""
    im_path = osp.join(TEST_DIR, 'AS500#10.dcm')
    picket_orientation = LEFT_RIGHT
    max_error = 0.24
    abs_median_error = 0.05


class AS500error(PFTestMixin, TestCase):
    """Tests for the AS500#2 image."""
    im_path = osp.join(TEST_DIR, 'AS500-error.dcm')
    num_pickets = 6
    percent_passing = 97.5
    max_error = 0.55
    abs_median_error = 0.07
    passes = False


class AS1000(PFTestMixin, TestCase):
    """Tests for the AS1000 image."""
    im_path = osp.join(TEST_DIR, 'AS1000_PF.dcm')
    max_error = 0.29
    abs_median_error = 0.06


class AS1000_2(PFTestMixin, TestCase):
    """Tests for the AS1000 image."""
    im_path = osp.join(TEST_DIR, 'AS1000#2.dcm')
    max_error = 0.24
    abs_median_error = 0.07


class AS1000_3(PFTestMixin, TestCase):
    """Tests for the AS1000 image."""
    im_path = osp.join(TEST_DIR, 'AS1000#3.dcm')
    max_error = 0.13
    abs_median_error = 0.05


class AS1000_4(PFTestMixin, TestCase):
    """Tests for the AS1000 image."""
    im_path = osp.join(TEST_DIR, 'AS1000#4.dcm')
    picket_orientation = LEFT_RIGHT
    max_error = 0.18
    abs_median_error = 0.05


class AS1000_90(PFTestMixin, TestCase):
    """Tests for the AS1000 image."""
    im_path = osp.join(TEST_DIR, 'AS1000-90.dcm')
    picket_orientation = LEFT_RIGHT
    max_error = 0.23
    abs_median_error = 0.05


class AS1000HDSmall(PFTestMixin, TestCase):
    """Tests for the AS1000 image."""
    im_path = osp.join(TEST_DIR, 'AS1000-HD-small.dcm')
    hdmlc = True
    max_error = 0.18
    abs_median_error = 0.05


class AS1000HDFull(PFTestMixin, TestCase):
    """Tests for the AS1000 image with a smaller pattern (only inner leaves)."""
    im_path = osp.join(TEST_DIR, 'AS1000-HD-full.dcm')
    hdmlc = True
    max_error = 0.2
    abs_median_error = 0.06


class AS1000HDFullVMAT(PFTestMixin, TestCase):
    """Tests for the AS1000 image with a smaller pattern (only inner leaves)."""
    im_path = osp.join(TEST_DIR, 'AS1000-HD-full-VMAT.dcm')
    hdmlc = True
    max_error = 0.2
    abs_median_error = 0.08


class AS1000HDFullError(PFTestMixin, TestCase):
    """Tests for the AS1000 image with a few errors introduced."""
    im_path = osp.join(TEST_DIR, 'AS1000-HD-full-error.dcm')
    hdmlc = True
    num_pickets = 6
    abs_median_error = 0.03
    max_error = 0.39

    def test_lower_tolerance_fails(self):
        """This image has an introduced error; this should catch with a reasonable tolerance."""
        pf = PicketFence(self.im_path)
        pf.analyze(tolerance=0.3, hdmlc=self.hdmlc)
        self.assertFalse(pf.passed)


class AS1200(PFTestMixin, TestCase):
    """Tests for the AS1200 image."""
    im_path = osp.join(TEST_DIR, 'AS1200.dcm')
    max_error = 0.08
    abs_median_error = 0.02


class AS1200Error(PFTestMixin, TestCase):
    """Tests for the AS1200 image."""
    im_path = osp.join(TEST_DIR, 'AS1200-error.dcm')
    num_pickets = 6
    max_error = 0.48
    abs_median_error = 0.05
    sag_adjustment = -1.2


class AS1200ExtendedSID(PFTestMixin, TestCase):
    """Tests for the AS1200 image."""
    im_path = osp.join(TEST_DIR, 'AS1200-ExtendedSID.dcm')
    max_error = 0.12
    abs_median_error = 0.04


class AS1200ExtendedSIDVMAT(PFTestMixin, TestCase):
    """Tests for the AS1200 image."""
    im_path = osp.join(TEST_DIR, 'AS1200-ExtendedSID-VMAT.dcm')
    max_error = 0.18
    abs_median_error = 0.06


@skip
class AS1200HD(PFTestMixin, TestCase):
    """Tests for the AS1200 image."""
    im_path = osp.join(TEST_DIR, 'AS1200-HD.dcm')
    hdmlc = True
    max_error = 0.05
    abs_median_error = 0.02

    @classmethod
    def setUpClass(cls):
        cls.pf = PicketFence(cls.im_path)
        cls.pf.analyze(hdmlc=cls.hdmlc, num_pickets=cls.num_pickets)


@skip
class AS1200HDTranslated(PFTestMixin, TestCase):
    """Tests for the AS1200 image."""
    im_path = osp.join(TEST_DIR, 'AS1200-HD-translated.dcm')
    hdmlc = True
    max_error = 0.05
    abs_median_error = 0.02
