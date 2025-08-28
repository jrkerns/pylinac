from unittest import TestCase

import numpy as np

from pylinac.core.scale import (
    MachineScale,
    _mirror_360,
    _noop,
    _shift_and_mirror_360,
    convert,
    wrap180,
    wrap360,
)


class TestHelperMethods(TestCase):
    def test_wrap360(self):
        arr = np.array([-365, -270, -180, -5, 0, 5, 180, 270, 365])
        nominal = np.array([355, 90, 180, 355, 0, 5, 180, 270, 5])
        actual = wrap360(arr)
        assert np.all(actual == nominal)

    def test_wrap180(self):
        arr = np.array([-365, -270, -180, -5, 0, 5, 180, 270, 365])
        nominal = np.array([-5, 90, -180, -5, 0, 5, -180, -90, 5])
        actual = wrap180(arr)
        assert np.all(actual == nominal)


class TestPrivateMethods(TestCase):
    def test_noop(self):
        assert 5 == _noop(5)
        assert -5.3 == _noop(-5.3)

    def test_mirror_360(self):
        assert _mirror_360(5) == 355
        assert _mirror_360(355) == 5

    def test_shift_and_mirror_360(self):
        assert 355 == _shift_and_mirror_360(185)
        assert 185 == _shift_and_mirror_360(355)
        assert 185 == _shift_and_mirror_360(-5)

    def test_noop_multiple(self):
        arr = np.array([-270, -180, -5, 0, 5, 180, 270, 365])
        nominal = np.array([-270, -180, -5, 0, 5, 180, 270, 365])
        actual = _noop(arr)
        assert np.all(actual == nominal)

    def test_mirror_360_multiple(self):
        arr = np.array([-270, -180, -5, 0, 5, 180, 270, 365])
        nominal = np.array([270, 180, 5, 0, 355, 180, 90, 355])
        actual = _mirror_360(arr)
        assert np.all(actual == nominal)

    def test_shift_and_mirror_360_multiple(self):
        arr = np.array([-270, -180, -5, 0, 5, 180, 270, 365])
        nominal = np.array([90, 0, 185, 180, 175, 0, 270, 175])
        actual = _shift_and_mirror_360(arr)
        assert np.all(actual == nominal)


class TestScaleConversion(TestCase):
    def test_iec_to_iec(self):
        # should return the same values
        g = 5
        c = 5
        r = 5
        go, co, ro = convert(
            input_scale=MachineScale.IEC61217,
            output_scale=MachineScale.IEC61217,
            gantry=g,
            collimator=c,
            rotation=r,
        )
        assert g == go
        assert c == co
        assert r == ro

        g = 355
        c = 355
        r = 355
        go, co, ro = convert(
            input_scale=MachineScale.IEC61217,
            output_scale=MachineScale.IEC61217,
            gantry=g,
            collimator=c,
            rotation=r,
        )
        assert g == go
        assert c == co
        assert r == ro

    def test_iec_to_varian_iec(self):
        g = 5
        c = 5
        r = 5
        go, co, ro = convert(
            input_scale=MachineScale.IEC61217,
            output_scale=MachineScale.VARIAN_IEC,
            gantry=g,
            collimator=c,
            rotation=r,
        )
        assert g == go
        assert c == co
        assert ro == 355

        g = 355
        c = 355
        r = 355
        go, co, ro = convert(
            input_scale=MachineScale.IEC61217,
            output_scale=MachineScale.VARIAN_IEC,
            gantry=g,
            collimator=c,
            rotation=r,
        )
        assert g == go
        assert c == co
        assert ro == 5

    def test_varian_iec_to_iec(self):
        g = 5
        c = 5
        r = 5
        go, co, ro = convert(
            input_scale=MachineScale.VARIAN_IEC,
            output_scale=MachineScale.IEC61217,
            gantry=g,
            collimator=c,
            rotation=r,
        )
        assert g == go
        assert c == co
        assert ro == 355

        g = 355
        c = 355
        r = 355
        go, co, ro = convert(
            input_scale=MachineScale.VARIAN_IEC,
            output_scale=MachineScale.IEC61217,
            gantry=g,
            collimator=c,
            rotation=r,
        )
        assert g == go
        assert c == co
        assert ro == 5

    def test_iec_to_varian_standard(self):
        g = 5
        c = 5
        r = 5
        go, co, ro = convert(
            input_scale=MachineScale.IEC61217,
            output_scale=MachineScale.VARIAN_STANDARD,
            gantry=g,
            collimator=c,
            rotation=r,
        )
        assert go == 175
        assert co == 175
        assert ro == 175

        g = 355
        c = 355
        r = 355
        go, co, ro = convert(
            input_scale=MachineScale.IEC61217,
            output_scale=MachineScale.VARIAN_STANDARD,
            gantry=g,
            collimator=c,
            rotation=r,
        )
        assert go == 185
        assert co == 185
        assert ro == 185

    def test_varian_standard_to_iec(self):
        g = 175
        c = 175
        r = 175
        go, co, ro = convert(
            input_scale=MachineScale.VARIAN_STANDARD,
            output_scale=MachineScale.IEC61217,
            gantry=g,
            collimator=c,
            rotation=r,
        )
        assert go == 5
        assert co == 5
        assert ro == 5

        g = 185
        c = 185
        r = 185
        go, co, ro = convert(
            input_scale=MachineScale.VARIAN_STANDARD,
            output_scale=MachineScale.IEC61217,
            gantry=g,
            collimator=c,
            rotation=r,
        )
        assert go == 355
        assert co == 355
        assert ro == 355

    def test_iec_to_elekta_iec(self):
        g = 5
        c = 5
        r = 5
        go, co, ro = convert(
            input_scale=MachineScale.IEC61217,
            output_scale=MachineScale.ELEKTA_IEC,
            gantry=g,
            collimator=c,
            rotation=r,
        )
        assert g == go
        assert c == co
        assert ro == 355

        g = 355
        c = 355
        r = 355
        go, co, ro = convert(
            input_scale=MachineScale.IEC61217,
            output_scale=MachineScale.ELEKTA_IEC,
            gantry=g,
            collimator=c,
            rotation=r,
        )
        assert g == go
        assert c == co
        assert ro == 5

    def test_elekta_iec_to_iec(self):
        g = 5
        c = 5
        r = 5
        go, co, ro = convert(
            input_scale=MachineScale.ELEKTA_IEC,
            output_scale=MachineScale.IEC61217,
            gantry=g,
            collimator=c,
            rotation=r,
        )
        assert g == go
        assert c == co
        assert ro == 355

        g = 355
        c = 355
        r = 355
        go, co, ro = convert(
            input_scale=MachineScale.ELEKTA_IEC,
            output_scale=MachineScale.IEC61217,
            gantry=g,
            collimator=c,
            rotation=r,
        )
        assert g == go
        assert c == co
        assert ro == 5
