import unittest
from unittest import TestCase

from tests.test_winstonlutz import WinstonLutzMixin, Vector


class KatyiX0(WinstonLutzMixin, TestCase):
    file_path = ['Katy iX', '0.zip']
    num_images = 17
    gantry_iso_size = 0.72
    gantry_iso2bb_vector = Vector(-0.3, 0, -0.2)
    gantry_sag = 1
    collimator_iso_size = 0.5
    collimator_iso2bb_vector = Vector(-0.6, -0.3)
    couch_iso_size = 0.7
    couch_iso2bb_vector = Vector(-0.2)
    cax2bb_max_distance = 1.07
    cax2bb_median_distance = 0.8
    variable_axes = {0: 'Reference'}


class KatyiX1(WinstonLutzMixin, TestCase):
    file_path = ['Katy iX', '1.zip']
    num_images = 17
    gantry_iso_size = 0.5
    gantry_iso2bb_vector = Vector(0.1, -0.2, 0.4)
    gantry_sag = 1
    collimator_iso_size = 0.35
    collimator_iso2bb_vector = Vector(0.3, 0.2, 0)
    couch_iso_size = 0.25
    couch_iso2bb_vector = Vector(0.1, 0, 0)
    cax2bb_max_distance = 1.2
    cax2bb_median_distance = 0.3
    variable_axes = {0: 'Collimator'}


class KatyiX2(WinstonLutzMixin, TestCase):
    file_path = ['Katy iX', '2.zip']
    num_images = 17
    gantry_iso_size = 0.5
    gantry_iso2bb_vector = Vector(0.25, -0.1, 0.3)
    gantry_sag = 1.1
    collimator_iso_size = 0.4
    collimator_iso2bb_vector = Vector(0.5, 0, 0)
    couch_iso_size = 0.9
    couch_iso2bb_vector = Vector(-0.4, 0.3, 0)
    cax2bb_max_distance = 1.1
    cax2bb_median_distance = 0.5
    variable_axes = {0: 'Gantry'}


class KatyiX3(WinstonLutzMixin, TestCase):
    file_path = ['Katy iX', '3 (with crosshair).zip']
    num_images = 17
    gantry_iso_size = 0.65
    gantry_iso2bb_vector = Vector(-0.4, 0, -0.2)
    gantry_sag = 1.25
    collimator_iso_size = 0.6
    collimator_iso2bb_vector = Vector(-0.4, -0.5)
    couch_iso_size = 0.9
    couch_iso2bb_vector = Vector(-0.2)
    cax2bb_max_distance = 1.25
    cax2bb_median_distance = 0.8
    variable_axes = {0: 'Collimator'}


@unittest.skip
class KatyTB0(WinstonLutzMixin, TestCase):
    # TODO: Investigate why set is giving DICOM error
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
    variable_axes = {0: 'Collimator'}


class KatyTB1(WinstonLutzMixin, TestCase):
    file_path = ['Katy TB', '1.zip']
    num_images = 16
    gantry_iso_size = 0.4
    gantry_iso2bb_vector = Vector(-0.5, -0.1, 0.4)
    gantry_sag = 1.2
    collimator_iso_size = 0.4
    collimator_iso2bb_vector = Vector(-0.5, -0.1)
    couch_iso_size = 0.5
    couch_iso2bb_vector = Vector(-0.3)
    cax2bb_max_distance = 1
    cax2bb_median_distance = 0.7
    variable_axes = {0: 'Reference'}


class KatyTB2(WinstonLutzMixin, TestCase):
    file_path = ['Katy TB', '2.zip']
    num_images = 17
    gantry_iso_size = 0.5
    gantry_iso2bb_vector = Vector(-0.2, -0.3, 0.5)
    gantry_sag = 1.1
    collimator_iso_size = 0.4
    collimator_iso2bb_vector = Vector(0.1, 0)
    couch_iso_size = 0.3
    couch_iso2bb_vector = Vector(-0.2, 0.3)
    cax2bb_max_distance = 1.1
    cax2bb_median_distance = 0.4
    variable_axes = {0: 'Couch'}


class ChicagoTBFinal(WinstonLutzMixin, TestCase):
    file_path = ['Chicago', 'WL-Final_C&G&C_Final.zip']
    num_images = 25
    gantry_iso_size = 0.37
    gantry_iso2bb_vector = Vector(-0.2)
    gantry_sag = 0.54
    collimator_iso_size = 0.1
    collimator_iso2bb_vector = Vector(-0.2, -0.3)
    couch_iso_size = 0.1
    couch_iso2bb_vector = Vector(-0.25, -0.2)
    cax2bb_max_distance = 0.5
    cax2bb_median_distance = 0.3
    variable_axes = {0: 'Gantry'}


class ChicagoTB52915(WinstonLutzMixin, TestCase):
    file_path = ['Chicago', 'WL_05-29-15_Final.zip']
    num_images = 16
    gantry_iso_size = 0.37
    gantry_iso2bb_vector = Vector(-0.1, -0.2, -0.1)
    gantry_sag = 0.5
    collimator_iso_size = 0.15
    collimator_iso2bb_vector = Vector(-0.05, -0.1)
    couch_iso_size = 0.2
    couch_iso2bb_vector = Vector(-0.1, -0.4)
    cax2bb_max_distance = 0.5
    cax2bb_median_distance = 0.3
    variable_axes = {0: 'Gantry'}


class TrueBeam3120213(WinstonLutzMixin, TestCase):
    file_path = ['TrueBeam 3', '120213.zip']
    num_images = 26
    cax2bb_max_distance = 0.9
    cax2bb_median_distance = 0.35
    gantry_sag = 0.95
    gantry_iso_size = 0.65
    gantry_iso2bb_vector = Vector(0.1, 0.1, 0.25)
    collimator_iso_size = 0.4
    collimator_iso2bb_vector = Vector(-0.1, 0.1)
    couch_iso_size = 0.3
    couch_iso2bb_vector = Vector(-0.1, 0.1)


class SugarLandiX1(WinstonLutzMixin, TestCase):
    file_path = ['Sugarland iX', '1.zip']
    num_images = 17
    gantry_iso_size = 0.55
    gantry_iso2bb_vector = Vector(0.8, 0.1, 1.0)
    gantry_sag = 0.85
    collimator_iso_size = 0.25
    collimator_iso2bb_vector = Vector(0.7, 0.6)
    couch_iso_size = 0.7
    couch_iso2bb_vector = Vector(1.15)
    cax2bb_max_distance = 1.67
    cax2bb_median_distance = 1.05
    variable_axes = {0: 'Gantry'}


class BayAreaiX0(WinstonLutzMixin, TestCase):
    file_path = ['Bay Area iX', '0.zip']
    num_images = 17
    gantry_iso_size = 0.7
    gantry_iso2bb_vector = Vector(0.3, 0.1, 0.4)
    gantry_sag = 1.0
    collimator_iso_size = 0.5
    collimator_iso2bb_vector = Vector(0.2, 0.4)
    couch_iso_size = 1.2
    couch_iso2bb_vector = Vector(-1.0, -0.4)
    cax2bb_max_distance = 1.25
    cax2bb_median_distance = 0.6
    variable_axes = {0: 'Reference'}