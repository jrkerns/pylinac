"""Travis CI memory can't handle all the CBCTs; thus only test them when explicitly asked to."""
import os.path as osp
from unittest import TestCase, skip

from tests import TEST_BANK_DIR
from tests.test_cbct import CBCTMixin


class CBCTBankMixin(CBCTMixin):
    dir_location = osp.join(TEST_BANK_DIR, 'CBCTs')
    

class VarianPelvis(CBCTBankMixin, TestCase):
    """Test the Varian Pelvis protocol CBCT."""
    file_path = ['Pelvis.zip']
    expected_roll = -0.24
    slice_locations = {'HU': 32, 'UN': 7, 'SR': 44, 'LC': 20}
    hu_values = {'Poly': -36, 'Acrylic': 114, 'Delrin': 342, 'Air': -993, 'Teflon': 992, 'PMP': -188, 'LDPE': -95}
    unif_values = {'Center': 17, 'Left': 5, 'Right': 4, 'Top': 4, 'Bottom': 4}
    mtf_values = {60: 0.65, 70: 0.560, 80: 0.48, 90: 0.40, 95: 0.30}
    avg_line_length = 49.8
    lowcon_visible = 3


class VarianPelvisSpotlight(CBCTBankMixin, TestCase):
    """Test the Varian Pelvis Spotlight protocol CBCT."""
    file_path = ['Pelvis spotlight.zip']
    expected_roll = -0.26
    slice_locations = {'HU': 32, 'UN': 7, 'SR': 44, 'LC': 20}
    hu_values = {'Poly': -43, 'Acrylic': 118, 'Delrin': 341, 'Air': -998, 'Teflon': 967, 'PMP': -198, 'LDPE': -100}
    unif_values = {'Center': 15, 'Left': 3, 'Right': -1, 'Top': -1, 'Bottom': 0}
    mtf_values = {60: 0.97, 70: 0.87, 80: 0.75, 90: 0.59, 95: 0.45}
    avg_line_length = 49.94
    lowcon_visible = 5


class VarianLowDoseThorax(CBCTBankMixin, TestCase):
    """Test the Varian Low-Dose Thorax protocol CBCT."""
    file_path = ['Low dose thorax.zip']
    expected_roll = -0.29
    slice_locations = {'HU': 32, 'UN': 7, 'SR': 44, 'LC': 20}
    hu_values = {'Poly': -46, 'Acrylic': 119, 'Delrin': 341, 'Air': -998, 'Teflon': 992, 'PMP': -193, 'LDPE': -97}
    unif_values = {'Center': 19, 'Left': 7, 'Right': -1, 'Top': 3, 'Bottom': 2}
    mtf_values = {60: 0.56, 70: 0.50, 80: 0.43, 90: 0.34, 95: 0.27}
    avg_line_length = 49.76
    lowcon_visible = 2


class VarianStandardHead(CBCTBankMixin, TestCase):
    """Test the Varian Standard Head protocol CBCT."""
    file_path = ['Standard head.zip']
    expected_roll = -0.19
    slice_locations = {'HU': 31, 'UN': 6, 'SR': 43, 'LC': 19}
    hu_values = {'Poly': -43, 'Acrylic': 124, 'Delrin': 345, 'Air': -991, 'Teflon': 997, 'PMP': -199, 'LDPE': -101}
    unif_values = {'Center': 17, 'Left': 15, 'Right': 4, 'Top': 9, 'Bottom': 9}
    mtf_values = {60: 0.95, 70: 0.85, 80: 0.72, 90: 0.36, 95: 0.32}
    avg_line_length = 49.94
    lowcon_visible = 1
    thickness_passed = False


class VarianLowDoseHead(CBCTBankMixin, TestCase):
    """Test the Varian Low-Dose Head protocol CBCT."""
    file_path = ['Low dose head.zip']
    expected_roll = -0.4
    slice_locations = {'HU': 32, 'UN': 6, 'SR': 44, 'LC': 20}
    hu_values = {'Poly': -41, 'Acrylic': 123, 'Delrin': 350, 'Air': -990, 'Teflon': 998, 'PMP': -200, 'LDPE': -103}
    unif_values = {'Center': 16, 'Left': 11, 'Right': 3, 'Top': 7, 'Bottom': 6}
    mtf_values = {60: 0.87, 70: 0.84, 80: 0.66, 90: 0.49, 95: 0.41}
    lowcon_visible = 1
    avg_line_length = 49.93


@skip  # SR slice fails; TODO: make slices modular, so if 1 fails, the others still run
class GEMonthlyCT(CBCTBankMixin, TestCase):
    """Test a monthly CT scan from GE."""
    file_path = ['GE_CT.zip']
    hu_tolerance = 90
    slice_locations = {'HU': 143, 'UN': 91, 'SR': 167, 'LC': 119}
    hu_values = {'Poly': -32, 'Acrylic': 119, 'Delrin': 333, 'Air': -944, 'Teflon': 909, 'PMP': -173, 'LDPE': -87}
    unif_values = {'Center': 11, 'Left': 11, 'Right': 11, 'Top': 11, 'Bottom': 11}
    mtf_values = {60: 0.51, 70: 0.45, 80: 0.39, 90: 0.30, 95: 0.25}
    lowcon_visible = 4


class ToshibaMonthlyCT(CBCTBankMixin, TestCase):
    """Test a monthly CT scan from Toshiba."""
    file_path = ['Toshiba.zip']
    hu_tolerance = 240
    slice_locations = {'HU': 36, 'UN': 14, 'SR': 46, 'LC': 26}
    hu_values = {'Poly': -32, 'Acrylic': 106, 'Delrin': 467, 'Air': -994, 'Teflon': 1214, 'PMP': -165, 'LDPE': -85}
    unif_values = {'Center': 8, 'Left': 7, 'Right': 7, 'Top': 7, 'Bottom': 6}
    mtf_values = {60: 0.48, 70: 0.40, 80: 0.33, 90: 0.28, 95: 0.24}
    lowcon_visible = 5


class CBCT1(CBCTBankMixin, TestCase):
    """A Varian CBCT dataset"""
    file_path = ['CBCT_1.zip']
    expected_roll = -0.53
    slice_locations = {'HU': 31, 'UN': 6, 'SR': 43, 'LC': 19}
    hu_values = {'Poly': -39, 'Acrylic': 130, 'Delrin': 347, 'Air': -986, 'Teflon': 1002, 'PMP': -189, 'LDPE': -90}
    unif_values = {'Center': 13, 'Left': 17, 'Right': 5, 'Top': 13, 'Bottom': 11}
    mtf_values = {60: 0.97, 70: 0.88, 80: 0.75, 90: 0.54, 95: 0.45}
    avg_line_length = 49.9
    lowcon_visible = 1
    thickness_passed = False


class CBCT2(CBCTBankMixin, TestCase):
    """A Varian CBCT dataset"""
    file_path = ['CBCT_2.zip']
    expected_roll = -0.2
    hu_tolerance = 50
    slice_locations = {'HU': 34, 'UN': 8, 'SR': 46, 'LC': 22}
    hu_values = {'Poly': -16, 'Acrylic': 135, 'Delrin': 367, 'Air': -967, 'Teflon': 1017, 'PMP': -163, 'LDPE': -71}
    unif_values = {'Center': 47, 'Left': 35, 'Right': 41, 'Top': 39, 'Bottom': 37}
    mtf_values = {60: 0.68, 70: 0.57, 80: 0.49, 90: 0.40, 95: 0.31}
    lowcon_visible = 3
    avg_line_length = 49.82


class CBCT3(CBCTBankMixin, TestCase):
    """A Varian CBCT dataset"""
    file_path = ['CBCT_3.zip']
    expected_roll = -2.66
    hu_tolerance = 50
    slice_locations = {'HU': 35, 'UN': 10, 'SR': 47, 'LC': 23}
    hu_values = {'Poly': -44, 'Acrylic': 110, 'Delrin': 325, 'Air': -979, 'Teflon': 949, 'PMP': -194, 'LDPE': -107}
    unif_values = {'Center': 2, 'Left': -1, 'Right': 11, 'Top': 9, 'Bottom': 2}
    mtf_values = {60: 0.64, 70: 0.56, 80: 0.48, 90: 0.40, 95: 0.31}
    avg_line_length = 49.87
    lowcon_visible = 4
    thickness_passed = False

# CBCT4 is in the regular test_cbct.py file


class CBCT5(CBCTBankMixin, TestCase):
    """A Varian CBCT dataset"""
    file_path = ['CBCT_5.zip']
    slice_locations = {'HU': 35, 'UN': 8, 'SR': 47, 'LC': 23}
    hu_values = {'Poly': -54, 'Acrylic': 101, 'Delrin': 328, 'Air': -999, 'Teflon': 975, 'PMP': -203, 'LDPE': -110}
    unif_values = {'Center': 19, 'Left': -8, 'Right': -5, 'Top': -7, 'Bottom': -6}
    mtf_values = {80: 0.46, 90: 0.38, 60: 0.59, 70: 0.52, 95: 0.29}
    avg_line_length = 49.55
    lowcon_visible = 5


class CBCT6(CBCTBankMixin, TestCase):
    """A Varian CBCT dataset"""
    file_path = ['CBCT_6.zip']
    expected_roll = 0.2
    slice_locations = {'HU': 38, 'UN': 12, 'SR': 50, 'LC': 26}
    hu_values = {'Poly': -42, 'Acrylic': 107, 'Delrin': 327, 'Air': -994, 'Teflon': 972, 'PMP': -192, 'LDPE': -100}
    unif_values = {'Center': -5, 'Left': 0, 'Right': -13, 'Top': -7, 'Bottom': -6}
    mtf_values = {80: 0.68, 90: 0.46, 60: 0.89, 70: 0.79, 95: 0.31}
    avg_line_length = 49.94
    lowcon_visible = 5


class CBCT7(CBCTBankMixin, TestCase):
    """A Varian CBCT dataset"""
    file_path = ['CBCT_7.zip']
    expected_roll = 0.5
    slice_locations = {'HU': 36, 'UN': 10, 'SR': 48, 'LC': 24}
    hu_values = {'Poly': -50, 'Acrylic': 108, 'Delrin': 334, 'Air': -999, 'Teflon': 982, 'PMP': -200, 'LDPE': -107}
    unif_values = {'Center': 14, 'Left': -5, 'Right': -5, 'Top': -5, 'Bottom': -5}
    mtf_values = {80: 0.46, 90: 0.37, 60: 0.6, 70: 0.53, 95: 0.28}
    avg_line_length = 49.65
    lowcon_visible = 3


class CBCT8(CBCTBankMixin, TestCase):
    """A Varian CBCT dataset"""
    file_path = ['CBCT_8.zip']
    expected_roll = 0.55
    slice_locations = {'HU': 39, 'UN': 14, 'SR': 51, 'LC': 27}
    hu_values = {'Poly': -37, 'Acrylic': 114, 'Delrin': 334, 'Air': -994, 'Teflon': 982, 'PMP': -186, 'LDPE': -97}
    unif_values = {'Center': -4, 'Left': 2, 'Right': -5, 'Top': 0, 'Bottom': -1}
    mtf_values = {80: 0.65, 90: 0.36, 60: 0.89, 70: 0.77, 95: 0.28}
    avg_line_length = 49.95
    lowcon_visible = 5


class CBCT9(CBCTBankMixin, TestCase):
    """A Varian CBCT dataset"""
    file_path = ['CBCT_9.zip']
    expected_roll = 0.4
    slice_locations = {'HU': 35, 'UN': 9, 'SR': 47, 'LC': 23}
    hu_values = {'Poly': -53, 'Acrylic': 105, 'Delrin': 330, 'Air': -999, 'Teflon': 978, 'PMP': -199, 'LDPE': -107}
    unif_values = {'Center': 10, 'Left': -5, 'Right': -4, 'Top': -6, 'Bottom': -5}
    mtf_values = {80: 0.49, 90: 0.41, 60: 0.65, 70: 0.56, 95: 0.32}
    avg_line_length = 49.64
    lowcon_visible = 3


class CBCT10(CBCTBankMixin, TestCase):
    """A Varian CBCT dataset"""
    file_path = ['CBCT_10.zip']
    expected_roll = 0.4
    slice_locations = {'HU': 38, 'UN': 12, 'SR': 50, 'LC': 26}
    hu_values = {'Poly': -37, 'Acrylic': 109, 'Delrin': 334, 'Air': -992, 'Teflon': 985, 'PMP': -186, 'LDPE': -93}
    unif_values = {'Center': -4, 'Left': 4, 'Right': -5, 'Top': 0, 'Bottom': -1}
    mtf_values = {80: 0.73, 90: 0.58, 60: 0.94, 70: 0.85, 95: 0.34}
    lowcon_visible = 4


class CBCT11(CBCTBankMixin, TestCase):
    """A Varian CBCT dataset"""
    file_path = ['CBCT_11.zip']
    expected_roll = 0.5
    slice_locations = {'HU': 38, 'UN': 12, 'SR': 50, 'LC': 26}
    hu_values = {'Poly': -37, 'Acrylic': 108, 'Delrin': 332, 'Air': -992, 'Teflon': 982, 'PMP': -189, 'LDPE': -95}
    unif_values = {'Center': -6, 'Left': 2, 'Right': -7, 'Top': -2, 'Bottom': -1}
    mtf_values = {80: 0.69, 90: 0.38, 60: 0.9, 70: 0.79, 95: 0.29}
    avg_line_length = 49.94
    lowcon_visible = 4


class CBCT12(CBCTBankMixin, TestCase):
    """A Varian CBCT dataset"""
    file_path = ['CBCT_12.zip']
    slice_locations = {'HU': 35, 'UN': 10, 'SR': 47, 'LC': 23}
    hu_values = {'Poly': -55, 'Acrylic': 112, 'Delrin': 335, 'Air': -999, 'Teflon': 982, 'PMP': -201, 'LDPE': -107}
    unif_values = {'Center': 5, 'Left': -5, 'Right': -9, 'Top': -7, 'Bottom': -6}
    mtf_values = {80: 0.42, 90: 0.32, 60: 0.54, 70: 0.48, 95: 0.26}
    avg_line_length = 49.59
    lowcon_visible = 2


class CBCT13(CBCTBankMixin, TestCase):
    """A Varian CBCT dataset"""
    file_path = ['CBCT_13.zip']
    expected_roll = 0.2
    slice_locations = {'HU': 36, 'UN': 10, 'SR': 48, 'LC': 24}
    hu_values = {'Poly': -53, 'Acrylic': 106, 'Delrin': 329, 'Air': -999, 'Teflon': 976, 'PMP': -200, 'LDPE': -107}
    unif_values = {'Center': 3, 'Left': -7, 'Right': -6, 'Top': -8, 'Bottom': -6}
    mtf_values = {80: 0.46, 90: 0.38, 60: 0.61, 70: 0.53, 95: 0.29}
    avg_line_length = 49.66
    lowcon_visible = 3


class CBCT14(CBCTBankMixin, TestCase):
    """A Varian CBCT dataset"""
    file_path = ['CBCT_14.zip']
    expected_roll = 0.84
    slice_locations = {'HU': 32, 'UN': 6, 'SR': 44, 'LC': 20}
    hu_values = {'Poly': -41, 'Acrylic': 125, 'Delrin': 334, 'Air': -995, 'Teflon': 986, 'PMP': -184, 'LDPE': -89}
    unif_values = {'Center': 18, 'Left': 13, 'Right': 15, 'Top': 14, 'Bottom': 14}
    mtf_values = {80: 0.42, 90: 0.33, 60: 0.54, 70: 0.48, 95: 0.26}
    avg_line_length = 49.4
    lowcon_visible = 3


class CBCT15(CBCTBankMixin, TestCase):
    """A Varian CBCT dataset."""
    file_path = ['CBCT_15.zip']
    hu_tolerance = 50
    slice_locations = {'HU': 60, 'UN': 28, 'SR': 75, 'LC': 45}
    hu_values = {'Poly': -32, 'Acrylic': 121, 'Delrin': 353, 'Air': -998, 'Teflon': 945, 'PMP': -186, 'LDPE': -93}
    unif_values = {'Center': -2, 'Left': 6, 'Right': 5, 'Top': 11, 'Bottom': 3}
    mtf_values = {80: 0.53, 90: 0.46, 60: 0.7, 70: 0.61, 95: 0.43}
    lowcon_visible = 6


class CBCT16(CBCTBankMixin, TestCase):
    """A Varian CBCT dataset"""
    file_path = ['CBCT_16.zip']
    expected_roll = 0.2
    slice_locations = {'HU': 32, 'UN': 6, 'SR': 44, 'LC': 20}
    hu_values = {'Poly': -37, 'Acrylic': 128, 'Delrin': 342, 'Air': -995, 'Teflon': 1000, 'PMP': -181, 'LDPE': -87}
    unif_values = {'Center': 17, 'Left': 20, 'Right': 18, 'Top': 19, 'Bottom': 19}
    mtf_values = {80: 0.42, 90: 0.33, 60: 0.53, 70: 0.48, 95: 0.26}
    avg_line_length = 49.6
    lowcon_visible = 3


class CBCT17(CBCTBankMixin, TestCase):
    """A Varian CBCT dataset"""
    file_path = ['CBCT_17.zip']
    expected_roll = 0.45
    slice_locations = {'HU': 34, 'UN': 8, 'SR': 46, 'LC': 22}
    hu_values = {'Poly': -46, 'Acrylic': 117, 'Delrin': 344, 'Air': -989, 'Teflon': 989, 'PMP': -197, 'LDPE': -101}
    unif_values = {'Center': 5, 'Left': 0, 'Right': -7, 'Top': -6, 'Bottom': -2}
    mtf_values = {80: 0.68, 90: 0.34, 60: 0.93, 70: 0.8, 95: 0.27}
    avg_line_length = 49.94
    lowcon_visible = 1


class Katy1(CBCTBankMixin, TestCase):
    """CBCT with very high HU values."""
    file_path = ['Katy iX', 'Monday, March 10, 2014 1-05-47 PM.zip']
    expected_roll = -0.3
    slice_locations = {'HU': 44, 'UN': 18, 'SR': 56, 'LC': 32}
    hu_values = {'Poly': 584, 'Acrylic': 797, 'Delrin': 1046, 'Air': -404, 'Teflon': 1720, 'PMP': 424, 'LDPE': 522}
    hu_passed = False
    unif_values = {'Center': 612, 'Left': 628, 'Right': 645, 'Top': 642, 'Bottom': 631}
    unif_passed = False
    mtf_values = {80: 0.66, 90: 0.41, 60: 0.9, 70: 0.8, 95: 0.3}
    lowcon_visible = 3


class Elekta4(CBCTBankMixin, TestCase):
    """An Elekta CBCT dataset"""
    file_path = ['Elekta_4.zip']
    slice_locations = {'HU': 130, 'UN': 20, 'SR': 100, 'LC': 100}
    hu_values = {'Poly': -319, 'Acrylic': -224, 'Delrin': -91, 'Air': -863, 'Teflon': 249, 'PMP': -401, 'LDPE': -352}
    hu_passed = False
    unif_values = {'Center': -273, 'Left': -267, 'Right': -266, 'Top': -267, 'Bottom': -267}
    unif_passed = False
    mtf_values = {80: 0.52, 90: 0.45, 60: 0.69, 70: 0.59, 95: 0.41}
    lowcon_visible = 1


class Elekta7(CBCTBankMixin, TestCase):
    """An Elekta CBCT dataset"""
    file_path = ['Elekta_7.zip']
    slice_locations = {'HU': 159, 'UN': 49, 'SR': 129, 'LC': 129}
    hu_values = {'Poly': -128, 'Acrylic': -16, 'Delrin': 141, 'Air': -782, 'Teflon': 541, 'PMP': -226, 'LDPE': -164}
    hu_passed = False
    unif_values = {'Center': -81, 'Left': -73, 'Right': -72, 'Top': -73, 'Bottom': -73}
    unif_passed = False
    mtf_values = {80: 0.5, 90: 0.43, 60: 0.68, 70: 0.58, 95: 0.36}


class Elekta8(CBCTBankMixin, TestCase):
    """An Elekta CBCT dataset"""
    file_path = ['Elekta_8.zip']
    slice_locations = {'HU': 161, 'UN': 51, 'SR': 131, 'LC': 131}
    hu_values = {'Poly': -324, 'Acrylic': -229, 'Delrin': -97, 'Air': -868, 'Teflon': 245, 'PMP': -406, 'LDPE': -356}
    hu_passed = False
    unif_values = {'Center': -293, 'Left': -286, 'Right': -285, 'Top': -286, 'Bottom': -286}
    unif_passed = False
    mtf_values = {80: 0.42, 90: 0.42, 60: 0.7, 70: 0.51, 95: 0.25}
