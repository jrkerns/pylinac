from enum import Enum
from typing import Tuple


def noop(value: float) -> float:
    return value

# def wrap_360(value: float) -> float:
#     return value % 360

def mirror_360(value: float) -> float:
    return 360 - value

def shift_and_mirror_360(value: float) -> float:
    return mirror_360(180 - value)

def wrap_1000(value: float) -> float:
    return value % 1000


class MachineScale(Enum):
    """Possible machine scales. Used for specifying input and and output scales for conversion.
    The enum keys are conversion functions for each axis relative to IEC 61217"""
    IEC61217 = {'gantry': noop, 'collimator': noop, 'rotation': noop}
    ELEKTA = "Elekta"
    ELEKTA_IEC = "Elekta IEC"
    VARIAN_IEC = {'gantry': noop, 'collimator': noop, 'rotation': mirror_360}
    VARIAN_STANDARD = {'gantry': shift_and_mirror_360, 'collimator': shift_and_mirror_360, 'rotation': shift_and_mirror_360}


def convert(input_scale: MachineScale, output_scale: MachineScale, gantry: float, collimator: float, couch: float) -> Tuple[float, float, float]:
    return
