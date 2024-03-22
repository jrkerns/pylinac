from pylinac.core.scale import (
    MachineScale,
    convert,
    inv_shift_and_mirror_360,
    mirror_360,
    noop,
    shift_and_mirror_360,
)


def test_noop():
    assert 5 == noop(5)
    assert -5.3 == noop(-5.3)


def test_mirror_360():
    assert mirror_360(5) == 355
    assert mirror_360(355) == 5


def test_shift_and_mirror_360():
    assert 355 == shift_and_mirror_360(185)
    assert 185 == shift_and_mirror_360(355)
    assert 185 == shift_and_mirror_360(-5)


def test_inv_shift_and_mirror_360():
    assert 180 == inv_shift_and_mirror_360(0)
    assert 175 == inv_shift_and_mirror_360(5)


def test_iec_to_iec():
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


def test_iec_to_varian_iec():
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


def test_varian_iec_to_iec():
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


def test_iec_to_varian_standard():
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


def test_varian_standard_to_iec():
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


def test_iec_to_elekta_iec():
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


def test_elekta_iec_to_iec():
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
