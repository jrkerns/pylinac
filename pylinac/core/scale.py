from __future__ import annotations

from enum import Enum

import argue


def noop(value: float) -> float:
    """Don't do anything."""
    return value


def mirror_360(value: float) -> float:
    """Mirror about 0"""
    argue.verify_bounds(value, argue.POSITIVE)
    return 360 - value


def shift_and_mirror_360(value: float) -> float:
    """Shift by 180 degrees and then mirror about 0"""
    argue.verify_bounds(value, argue.POSITIVE)
    v = value - 180
    if v > 0:
        return mirror_360(v)
    else:
        return abs(v)


def inv_shift_and_mirror_360(value: float) -> float:
    """Inverse shift and mirror"""
    v = 180 - value
    if v < 0:
        return mirror_360(abs(v))
    else:
        return v


class MachineScale(Enum):
    """Possible machine scales. Used for specifying input and output scales for conversion.
    The enum keys are conversion functions for each axis relative to IEC 61217"""

    IEC61217 = {
        "gantry_to_iec": noop,
        "collimator_to_iec": noop,
        "rotation_to_iec": noop,
        "gantry_from_iec": noop,
        "collimator_from_iec": noop,
        "rotation_from_iec": noop,
    }
    ELEKTA_IEC = {
        "gantry_to_iec": noop,
        "collimator_to_iec": noop,
        "rotation_to_iec": mirror_360,
        "gantry_from_iec": noop,
        "collimator_from_iec": noop,
        "rotation_from_iec": mirror_360,
    }
    VARIAN_IEC = {
        "gantry_to_iec": noop,
        "collimator_to_iec": noop,
        "rotation_to_iec": mirror_360,
        "gantry_from_iec": noop,
        "collimator_from_iec": noop,
        "rotation_from_iec": mirror_360,
    }
    VARIAN_STANDARD = {
        "gantry_to_iec": shift_and_mirror_360,
        "collimator_to_iec": shift_and_mirror_360,
        "rotation_to_iec": shift_and_mirror_360,
        "gantry_from_iec": inv_shift_and_mirror_360,
        "collimator_from_iec": inv_shift_and_mirror_360,
        "rotation_from_iec": inv_shift_and_mirror_360,
    }


def convert(
    input_scale: MachineScale,
    output_scale: MachineScale,
    gantry: float,
    collimator: float,
    rotation: float,
) -> (float, float, float):
    """Convert from one coordinate scale to another. Returns gantry, collimator, rotation."""
    # convert to IEC61217 since everything is defined relative to it
    g = input_scale.value["gantry_to_iec"](gantry)
    c = input_scale.value["collimator_to_iec"](collimator)
    r = input_scale.value["rotation_to_iec"](rotation)
    # now apply the inverse to go from 61217 to output scale
    g_out = output_scale.value["gantry_from_iec"](g)
    c_out = output_scale.value["collimator_from_iec"](c)
    r_out = output_scale.value["rotation_from_iec"](r)
    return g_out, c_out, r_out


def abs360(value: float) -> float:
    """Convert angles to always be positive. E.g. -90 -> 270"""
    return (360 + value) % 360


def wrap360(value: float) -> float:
    """Wrap the input value around 360. E.g. 361 -> 1"""
    return value % 360
