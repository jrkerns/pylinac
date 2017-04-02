import unittest
from unittest import TestCase, skip

from tests.test_winstonlutz import WinstonLutzMixin, Vector


class KatyiX0(WinstonLutzMixin, TestCase):
    file_path = ['Katy iX', '0.zip']
    num_images = 17
    gantry_iso_size = 1
    gantry_iso2bb_vector = Vector(0.45, 0.7, 0.2)
    gantry_sag = 0.8
    epid_sag = 1.1
    collimator_iso_size = 1.1
    collimator_iso2bb_vector = Vector(0.6, -0.3)
    couch_iso_size = 1.3
    couch_iso2bb_vector = Vector(0.2, -0.1)
    cax2bb_max_distance = 1.07
    cax2bb_median_distance = 0.8


class KatyiX1(WinstonLutzMixin, TestCase):
    file_path = ['Katy iX', '1.zip']
    num_images = 17
    gantry_iso_size = 1.1
    gantry_iso2bb_vector = Vector(-0.5, -0.5, -0.4)
    gantry_sag = 1
    epid_sag = 1.2
    collimator_iso_size = 0.7
    collimator_iso2bb_vector = Vector(-0.3, 0.2, 0)
    couch_iso_size = 0.6
    couch_iso2bb_vector = Vector(0.05, -0.1, 0)
    cax2bb_max_distance = 1.2
    cax2bb_median_distance = 0.3


class KatyiX2(WinstonLutzMixin, TestCase):
    file_path = ['Katy iX', '2.zip']
    num_images = 17
    gantry_iso_size = 0.9
    gantry_iso2bb_vector = Vector(-0.2, -0.2, -0.3)
    gantry_sag = 0.8
    epid_sag = 1.2
    collimator_iso_size = 0.8
    collimator_iso2bb_vector = Vector(-0.5, 0, 0)
    couch_iso_size = 1
    couch_iso2bb_vector = Vector(0, 0.3, 0)
    cax2bb_max_distance = 1.1
    cax2bb_median_distance = 0.5


class KatyiX3(WinstonLutzMixin, TestCase):
    file_path = ['Katy iX', '3 (with crosshair).zip']
    num_images = 17
    gantry_iso_size = 1.3
    gantry_iso2bb_vector = Vector(0.4, 0.3, 0.3)
    gantry_sag = 1
    epid_sag = 1.5
    collimator_iso_size = 1.5
    collimator_iso2bb_vector = Vector(0.45, -0.65)
    couch_iso_size = 1.9
    couch_iso2bb_vector = Vector(0.2, 0.2)
    cax2bb_max_distance = 1.25
    cax2bb_median_distance = 0.8


@unittest.skip  # TODO: Investigate why set is giving DICOM error
class KatyTB0(WinstonLutzMixin, TestCase):
    file_path = ['Katy TB', '0.zip']
    num_images = 17
    gantry_iso_size = 0.72
    gantry_iso2bb_vector = Vector(-0.3, 0.1, -0.2)
    gantry_sag = 1.25
    collimator_iso_size = 0.5
    collimator_iso2bb_vector = Vector(-0.6, -0.3)
    couch_iso_size = 0.9
    couch_iso2bb_vector = Vector(-0.2)
    cax2bb_max_distance = 1.07
    cax2bb_median_distance = 0.8


class KatyTB1(WinstonLutzMixin, TestCase):
    file_path = ['Katy TB', '1.zip']
    num_images = 16
    gantry_iso_size = 0.9
    gantry_iso2bb_vector = Vector(0.5, 0.1, -0.4)
    gantry_sag = 0.6
    epid_sag = 0.6
    collimator_iso_size = 0.8
    collimator_iso2bb_vector = Vector(0.6, -0.1)
    couch_iso_size = 0.9
    couch_iso2bb_vector = Vector(0.6, 0.2)
    cax2bb_max_distance = 1
    cax2bb_median_distance = 0.7


class KatyTB2(WinstonLutzMixin, TestCase):
    file_path = ['Katy TB', '2.zip']
    num_images = 17
    gantry_iso_size = 0.9
    gantry_iso2bb_vector = Vector(0, 0.7, -0.5)
    gantry_sag = 0.6
    epid_sag = 0.64
    collimator_iso_size = 0.8
    collimator_iso2bb_vector = Vector(-0.1, 0)
    couch_iso_size = 0.6
    couch_iso2bb_vector = Vector(0.4, 0.1)
    cax2bb_max_distance = 1.1
    cax2bb_median_distance = 0.4


class ChicagoTBFinal(WinstonLutzMixin, TestCase):
    file_path = ['Chicago', 'WL-Final_C&G&C_Final.zip']
    num_images = 17
    gantry_iso_size = 0.91
    gantry_iso2bb_vector = Vector(-0.1, 0, 0)
    gantry_sag = 0.5
    epid_sag = 0.43
    collimator_iso_size = 0.1
    collimator_iso2bb_vector = Vector(0.2, -0.3)
    couch_iso_size = 0.2
    couch_iso2bb_vector = Vector(0.25, -0.2)
    cax2bb_max_distance = 0.5
    cax2bb_median_distance = 0.3


class ChicagoTB52915(WinstonLutzMixin, TestCase):
    file_path = ['Chicago', 'WL_05-29-15_Final.zip']
    num_images = 16
    gantry_iso_size = 0.37
    gantry_iso2bb_vector = Vector(0, -0.2, 0.1)
    gantry_sag = 0.5
    epid_sag = 0.25
    collimator_iso_size = 0.15
    collimator_iso2bb_vector = Vector(0.1, -0.1)
    couch_iso_size = 0.2
    couch_iso2bb_vector = Vector(0.3, -0.2)
    cax2bb_max_distance = 0.5
    cax2bb_median_distance = 0.3


class TrueBeam3120213(WinstonLutzMixin, TestCase):
    file_path = ['TrueBeam 3', '120213.zip']
    num_images = 26
    cax2bb_max_distance = 0.9
    cax2bb_median_distance = 0.35
    gantry_sag = 0.95
    epid_sag = 0.6
    gantry_iso_size = 1.1
    gantry_iso2bb_vector = Vector(0.3, -0.1, -0.2)
    collimator_iso_size = 0.7
    collimator_iso2bb_vector = Vector(0.1, 0.1)
    couch_iso_size = 0.7
    couch_iso2bb_vector = Vector(0.2, 0.2)


class SugarLandiX2015(WinstonLutzMixin, TestCase):
    file_path = ['Sugarland iX', '2015', 'Lutz2.zip']
    num_images = 17
    gantry_iso_size = 1.25
    gantry_iso2bb_vector = Vector(-0.4, -0.1, -1.0)
    gantry_sag = 0.8
    epid_sag = 1.15
    collimator_iso_size = 0.5
    collimator_iso2bb_vector = Vector(-0.6, 0.6)
    couch_iso_size = 1.5
    couch_iso2bb_vector = Vector(-1.1)
    cax2bb_max_distance = 1.67
    cax2bb_median_distance = 1.05


class BayAreaiX0(WinstonLutzMixin, TestCase):
    file_path = ['Bay Area iX', '0.zip']
    num_images = 17
    gantry_iso_size = 1
    gantry_iso2bb_vector = Vector(-0.4, 0.2, -0.5)
    gantry_sag = 1.0
    epid_sag = 1.4
    collimator_iso_size = 1.1
    collimator_iso2bb_vector = Vector(-0.2, 0.4)
    couch_iso_size = 3.5
    couch_iso2bb_vector = Vector(1.3, -0.7)
    cax2bb_max_distance = 1.25
    cax2bb_median_distance = 0.6


@skip  # not a full dataset, but important to test that dead pixels are accounted for.
class DeadPixel(WinstonLutzMixin, TestCase):
    file_path = ['DeadPixelWL.zip']
    num_images = 8
    gantry_iso_size = 0.2
    gantry_iso2bb_vector = Vector(0.3, 0.1, 0.4)
    gantry_sag = 0.4
    collimator_iso_size = 0.5
    collimator_iso2bb_vector = Vector(0.2, 0.4)
    couch_iso_size = 1.2
    couch_iso2bb_vector = Vector(-1.0, -0.4)
    cax2bb_max_distance = 1.25
    cax2bb_median_distance = 0.6


class DAmoursElektaOffset(WinstonLutzMixin, TestCase):
    """An Elekta dataset, with the BB centered."""
    file_path = ['Michel DAmours - WLGantry_Offset_x=-1cm,y=+1cm,z=-1cm.zip']
    num_images = 8
    gantry_iso_size = 1.1
    gantry_iso2bb_vector = Vector(-10.3, 11, -9.3)
    gantry_sag = 1
    epid_sag = 2.6
    collimator_iso_size = None
    collimator_iso2bb_vector = None
    couch_iso_size = None
    couch_iso2bb_vector = None
    cax2bb_max_distance = 1.25
    cax2bb_median_distance = 0.6


class DAmoursElektaXOffset(WinstonLutzMixin, TestCase):
    """An Elekta dataset, with the BB centered."""
    file_path = ['Michel D\'Amours - WL_Shift_x=+1cm.zip']
    num_images = 8
    gantry_iso_size = 1.1
    gantry_iso2bb_vector = Vector(9.5, 0, 0.2)
    gantry_sag = 1
    epid_sag = 2.4
    collimator_iso_size = None
    collimator_iso2bb_vector = None
    couch_iso_size = None
    couch_iso2bb_vector = None
    cax2bb_max_distance = 1.25
    cax2bb_median_distance = 0.6


class DAmoursElektaCentered(WinstonLutzMixin, TestCase):
    """An Elekta dataset, with the BB centered."""
    file_path = ['Michel D\'Amours - GantryWL_BBCentered.zip']
    num_images = 8
    gantry_iso_size = 1.1
    gantry_iso2bb_vector = Vector(0, 0, 0.3)
    gantry_sag = 1
    epid_sag = 2.25
    collimator_iso_size = None
    collimator_iso2bb_vector = None
    couch_iso_size = None
    couch_iso2bb_vector = None
    cax2bb_max_distance = 1.25
    cax2bb_median_distance = 0.6


class DeBr6XElekta(WinstonLutzMixin, TestCase):
    """An Elekta dataset, with the BB centered."""
    file_path = ['DeBr', '6X_Elekta_Ball_Bearing.zip']
    num_images = 8
    gantry_iso_size = 1.6
    gantry_iso2bb_vector = Vector(-0.3, 0, 0)
    gantry_sag = 1.2
    epid_sag = 2.6
    collimator_iso_size = 0.5
    collimator_iso2bb_vector = Vector(-0.7, -0.4, 0)
    couch_iso_size = None
    couch_iso2bb_vector = None
    cax2bb_max_distance = 1.0
    cax2bb_median_distance = 0.7