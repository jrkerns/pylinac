from unittest import TestCase

import numpy as np

from tests_basic.test_winstonlutz import WinstonLutzMixin, Vector


class KatyiX0(WinstonLutzMixin, TestCase):
    # independently verified
    file_path = ['Katy iX', '0.zip']
    num_images = 17
    gantry_iso_size = 1
    collimator_iso_size = 1
    couch_iso_size = 1.3
    cax2bb_max_distance = 1.2
    cax2bb_median_distance = 0.8
    bb_shift_vector = Vector(x=-0.5, y=0.4, z=-0.5)
    print_results = True


class KatyiX1(WinstonLutzMixin, TestCase):
    file_path = ['Katy iX', '1.zip']
    num_images = 17
    gantry_iso_size = 1.1
    collimator_iso_size = 0.7
    couch_iso_size = 0.6
    cax2bb_max_distance = 1.2
    cax2bb_median_distance = 0.3
    bb_shift_vector = Vector(x=0.3, y=-0.2, z=0.3)


class KatyiX2(WinstonLutzMixin, TestCase):
    file_path = ['Katy iX', '2.zip']
    num_images = 17
    gantry_iso_size = 0.9
    collimator_iso_size = 0.8
    couch_iso_size = 1.5
    cax2bb_max_distance = 1.1
    cax2bb_median_distance = 0.5
    bb_shift_vector = Vector(x=0.4, y=-0.1, z=0.1)


class KatyiX3(WinstonLutzMixin, TestCase):
    file_path = ['Katy iX', '3 (with crosshair).zip']
    num_images = 17
    gantry_iso_size = 1.1
    collimator_iso_size = 1.3
    couch_iso_size = 1.8
    cax2bb_max_distance = 1.25
    cax2bb_median_distance = 0.8
    bb_shift_vector = Vector(x=-0.3, y=0.4, z=-0.5)


class KatyTB0(WinstonLutzMixin, TestCase):
    file_path = ['Katy TB', '0.zip']
    num_images = 17
    gantry_iso_size = 0.9
    collimator_iso_size = 0.8
    couch_iso_size = 1.3
    cax2bb_max_distance = 1.07
    cax2bb_median_distance = 0.8
    bb_shift_vector = Vector(x=-0.7, y=-0.1, z=-0.2)


class KatyTB1(WinstonLutzMixin, TestCase):
    file_path = ['Katy TB', '1.zip']
    num_images = 16
    gantry_iso_size = 0.9
    collimator_iso_size = 0.8
    couch_iso_size = 1.1
    cax2bb_max_distance = 1
    cax2bb_median_distance = 0.7
    bb_shift_vector = Vector(x=-0.6, y=-0.2)


class KatyTB2(WinstonLutzMixin, TestCase):
    file_path = ['Katy TB', '2.zip']
    num_images = 17
    gantry_iso_size = 1
    collimator_iso_size = 0.7
    couch_iso_size = 0.7
    cax2bb_max_distance = 1.1
    cax2bb_median_distance = 0.4
    bb_shift_vector = Vector(x=0.0, y=-0.2, z=-0.6)


class ChicagoTBFinal(WinstonLutzMixin, TestCase):
    # verified independently
    file_path = ['Chicago', 'WL-Final_C&G&C_Final.zip']
    num_images = 17
    gantry_iso_size = 0.91
    collimator_iso_size = 0.1
    couch_iso_size = 0.3
    cax2bb_max_distance = 0.5
    cax2bb_median_distance = 0.3
    bb_shift_vector = Vector(y=0.1)


class ChicagoTB52915(WinstonLutzMixin, TestCase):
    file_path = ['Chicago', 'WL_05-29-15_Final.zip']
    num_images = 16
    gantry_iso_size = 0.6
    collimator_iso_size = 0.3
    couch_iso_size = 0.3
    cax2bb_max_distance = 0.5
    cax2bb_median_distance = 0.3
    bb_shift_vector = Vector(z=0.2)


class TrueBeam3120213(WinstonLutzMixin, TestCase):
    file_path = ['TrueBeam 3', '120213.zip']
    num_images = 26
    cax2bb_max_distance = 0.9
    cax2bb_median_distance = 0.35
    gantry_iso_size = 1.1
    collimator_iso_size = 0.7
    couch_iso_size = 0.7
    bb_shift_vector = Vector(x=-0.1, y=-0.2, z=0.2)


class SugarLandiX2015(WinstonLutzMixin, TestCase):
    file_path = ['Sugarland iX', '2015', 'Lutz2.zip']
    num_images = 17
    gantry_iso_size = 1.3
    collimator_iso_size = 0.5
    couch_iso_size = 1.3
    cax2bb_max_distance = 1.67
    cax2bb_median_distance = 1.05
    bb_shift_vector = Vector(x=0.4, y=-0.7, z=0.1)


class BayAreaiX0(WinstonLutzMixin, TestCase):
    file_path = ['Bay Area iX', '0.zip']
    num_images = 17
    gantry_iso_size = 1
    collimator_iso_size = 1.1
    couch_iso_size = 2.3
    cax2bb_max_distance = 1.25
    cax2bb_median_distance = 0.6
    bb_shift_vector = Vector(x=0.3, y=-0.4, z=-0.2)


class DAmoursElektaOffset(WinstonLutzMixin, TestCase):
    """An Elekta dataset, with the BB centered."""
    file_path = ['Michel DAmours - WLGantry_Offset_x=-1cm,y=+1cm,z=-1cm.zip']
    num_images = 8
    gantry_iso_size = 1.1
    cax2bb_max_distance = 17.5
    cax2bb_median_distance = 14.3
    bb_shift_vector = Vector(x=10.2, y=-9.2, z=-11.1)  # independently verified


class DAmoursElektaXOffset(WinstonLutzMixin, TestCase):
    """An Elekta dataset, with the BB centered."""
    file_path = ['Michel D\'Amours - WL_Shift_x=+1cm.zip']
    num_images = 8
    gantry_iso_size = 1.1
    cax2bb_max_distance = 9.5
    cax2bb_median_distance = 6.9
    bb_shift_vector = Vector(x=-9.5, y=0.3, z=0.1)  # independently verified


class DAmoursElektaCentered(WinstonLutzMixin, TestCase):
    """An Elekta dataset, with the BB centered."""
    file_path = ['Michel D\'Amours - GantryWL_BBCentered.zip']
    num_images = 8
    gantry_iso_size = 1.1
    collimator_iso_size = None
    couch_iso_size = None
    cax2bb_max_distance = 0.8
    cax2bb_median_distance = 0.5
    bb_shift_vector = Vector(y=0.4)


class DeBr6XElekta(WinstonLutzMixin, TestCase):
    """An Elekta dataset, with the BB centered."""
    file_path = ['DeBr', '6X_Elekta_Ball_Bearing.zip']
    num_images = 8
    gantry_iso_size = 1.9
    collimator_iso_size = 1.9
    couch_iso_size = None
    cax2bb_max_distance = 1.0
    cax2bb_median_distance = 0.7
    bb_shift_vector = Vector(x=0.4, y=-0.2)
