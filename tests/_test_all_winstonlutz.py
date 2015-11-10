from unittest import TestCase

from tests.test_winstonlutz import WinstonLutzMixin, Vector, osp


class WLKiX(WinstonLutzMixin, TestCase):
    image_dir = osp.abspath(osp.join('..', '..', 'unorganized linac data', 'Winston-Lutz', 'Katy iX'))
    num_images = 17
    gantry_iso_size = 0.72
    gantry_iso2bb_vector = Vector(-0.3, 0.1, -0.2)
    gantry_sag = 1
    collimator_iso_size = 0.5
    collimator_iso2bb_vector = Vector(-0.6, -0.3)
    couch_iso_size = 0.7
    couch_iso2bb_vector = Vector(-0.2)
    cax2bb_max_distance = 1.07
    cax2bb_median_distance = 0.8
    variable_axes = {0: 'Reference'}


class WLChTB(WinstonLutzMixin, TestCase):
    image_dir = osp.abspath(osp.join('..', '..', 'unorganized linac data', 'Winston-Lutz', 'Chicago', 'WL-Final_C&G&C_Final'))
    num_images = 17
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