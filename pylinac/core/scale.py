from __future__ import annotations

from enum import Enum

from numpy import ndarray


def _noop(value: float | ndarray) -> float | ndarray:
    """Don't do anything."""
    return value


def _mirror_360(value: float | ndarray) -> float | ndarray:
    """Mirror about 0"""
    return wrap360(-value)


def _shift_and_mirror_360(value: float | ndarray) -> float | ndarray:
    """Shift by 180 degrees and then mirror about 0"""
    return wrap360(180 - value)


def wrap360(value: float | ndarray) -> float | ndarray:
    """Wrap the input values to the interval [0, 360)"""
    return value % 360


def wrap180(value: float | ndarray) -> float | ndarray:
    """Wrap the input values to the interval [-180, 180)"""
    return wrap360(value + 180) - 180


class MachineScale(Enum):
    """Possible machine scales. Used for specifying input and output scales for conversion.
    The enum keys are conversion functions for each axis relative to IEC 61217"""

    IEC61217 = {
        "gantry_to_iec": _noop,
        "collimator_to_iec": _noop,
        "rotation_to_iec": _noop,
        "gantry_from_iec": _noop,
        "collimator_from_iec": _noop,
        "rotation_from_iec": _noop,
    }
    ELEKTA_IEC = {
        "gantry_to_iec": _noop,
        "collimator_to_iec": _noop,
        "rotation_to_iec": _mirror_360,
        "gantry_from_iec": _noop,
        "collimator_from_iec": _noop,
        "rotation_from_iec": _mirror_360,
    }
    VARIAN_IEC = {
        "gantry_to_iec": _noop,
        "collimator_to_iec": _noop,
        "rotation_to_iec": _mirror_360,
        "gantry_from_iec": _noop,
        "collimator_from_iec": _noop,
        "rotation_from_iec": _mirror_360,
    }
    VARIAN_STANDARD = {
        "gantry_to_iec": _shift_and_mirror_360,
        "collimator_to_iec": _shift_and_mirror_360,
        "rotation_to_iec": _shift_and_mirror_360,
        "gantry_from_iec": _shift_and_mirror_360,
        "collimator_from_iec": _shift_and_mirror_360,
        "rotation_from_iec": _shift_and_mirror_360,
    }


def convert(
    input_scale: MachineScale,
    output_scale: MachineScale,
    gantry: float | ndarray,
    collimator: float | ndarray,
    rotation: float | ndarray,
) -> tuple[float | ndarray, float | ndarray, float | ndarray]:
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
