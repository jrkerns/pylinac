"""Utility functions for pylinac."""
from __future__ import annotations

import os
import os.path as osp
import struct
from abc import abstractmethod
from collections.abc import Iterable
from datetime import datetime
from enum import Enum
from typing import BinaryIO, Generic, Sequence, TypeVar

import numpy as np
import pydicom
from pydantic import BaseModel, ConfigDict, Field

from .. import __version__


def convert_to_enum(value: str | Enum | None, enum: type[Enum]) -> Enum:
    """Convert a value to an enum representation from an enum value if needed"""
    if isinstance(value, enum):
        return value
    else:
        return enum(value)


class OptionListMixin:
    """A mixin class that will create a list of the class attributes.
    Used for enum-like classes"""

    @classmethod
    def options(cls) -> list[str]:
        return [
            option
            for attr, option in cls.__dict__.items()
            if not callable(option) and not attr.startswith("__")
        ]


class ResultBase(BaseModel):
    model_config = ConfigDict(
        arbitrary_types_allowed=True
    )  # https://docs.pydantic.dev/latest/api/config/#pydantic.config.ConfigDict.arbitrary_types_allowed
    pylinac_version: str = __version__  #:
    date_of_analysis: datetime = Field(default_factory=datetime.today)  #:


T = TypeVar("T")


class ResultsDataMixin(Generic[T]):
    """A mixin for classes that generate results data. This mixin is used to generate the results data and present it in different formats.
    The generic types allow correct type hinting of the results data."""

    @abstractmethod
    def _generate_results_data(self) -> T:
        pass

    def results_data(
        self, as_dict: bool = False, as_json: bool = False
    ) -> T | dict | str:
        """Present the results data and metadata as a dataclass, dict, or tuple.
        The default return type is a dataclass.

        Parameters
        ----------
        as_dict : bool
            If True, return the results as a dictionary.
        as_json : bool
            If True, return the results as a JSON string. Cannot be True if as_dict is True.
        """
        if as_dict and as_json:
            raise ValueError("Cannot return as both dict and JSON. Pick one.")
        data = self._generate_results_data()
        if as_dict:
            return data.model_dump()
        if as_json:
            return data.model_dump_json()
        return data


def clear_data_files():
    """Delete all demo files, image classifiers, etc from the demo folder"""
    demo_folder = osp.join(osp.dirname(osp.dirname(__file__)), "demo_files")
    if osp.isdir(demo_folder):
        for file in os.listdir(demo_folder):
            full_file = osp.join(demo_folder, file)
            if osp.isfile(full_file):
                os.remove(full_file)
    print("Pylinac data files cleared.")


def assign2machine(source_file: str, machine_file: str):
    """Assign a DICOM RT Plan file to a specific machine. The source file is overwritten to contain
    the machine of the machine file.

    Parameters
    ----------
    source_file : str
        Path to the DICOM RTPlan file that contains the fields/plan desired
        (e.g. a Winston Lutz set of fields or Varian's default PF files).
    machine_file : str
        Path to a DICOM RTPlan file that has the desired machine. This is easily obtained from pushing a plan from the TPS
        for that specific machine. The file must contain at least one valid field.
    """
    dcm_source = pydicom.dcmread(source_file)
    dcm_machine = pydicom.dcmread(machine_file)
    for beam in dcm_source.BeamSequence:
        beam.TreatmentMachineName = dcm_machine.BeamSequence[0].TreatmentMachineName
    dcm_source.save_as(source_file)


def is_close(val: float, target: float | Sequence, delta: float = 1):
    """Return whether the value is near the target value(s).

    Parameters
    ----------
    val : number
        The value being compared against.
    target : number, iterable
        If a number, the values are simply evaluated.
        If a sequence, each target is compared to ``val``.
        If any values of ``target`` are close, the comparison is considered True.

    Returns
    -------
    bool
    """
    try:
        targets = (value for value in target)
    except (AttributeError, TypeError):
        targets = [target]
    for target in targets:
        if target - delta < val < target + delta:
            return True
    return False


def simple_round(number: float | int, decimals: int | None = 0) -> float | int:
    """Round a number to the given number of decimals. Fixes small floating number errors. If decimals is None, no rounding is performed"""
    if decimals is None:
        return number
    num = int(round(number * 10**decimals))
    if decimals >= 1:
        num /= 10**decimals
    return num


def is_iterable(object) -> bool:
    """Determine if an object is iterable."""
    return isinstance(object, Iterable)


class Structure:
    """A simple structure that assigns the arguments to the object."""

    def __init__(self, **kwargs):
        self.__dict__.update(**kwargs)

    def update(self, **kwargs):
        self.__dict__.update(**kwargs)


def decode_binary(
    file: BinaryIO,
    dtype: type[int] | type[float] | type[str] | str | np.dtype,
    num_values: int = 1,
    cursor_shift: int = 0,
    strip_empty: bool = True,
) -> int | float | str | np.ndarray | list:
    """Read in a raw binary file and convert it to given data types.

    Parameters
    ----------
    file
        The open file object.
    dtype
        The expected data type to return. If int or float and num_values > 1, will return numpy array.
    num_values
        The expected number of dtype to return

        .. note:: This is not the same as the number of bytes.

    cursor_shift : int
        The number of bytes to move the cursor forward after decoding. This is used if there is a
        reserved section after the read-in segment.
    strip_empty : bool
        Whether to strip trailing empty/null values for strings.
    """
    f = file

    if isinstance(dtype, str):
        s = struct.calcsize(dtype) * num_values
        output = struct.unpack(dtype * num_values, f.read(s))
        if len(output) == 1:
            output = output[0]
    elif dtype == str:  # if string
        ssize = struct.calcsize("c") * num_values
        output = struct.unpack("c" * num_values, f.read(ssize))
        if strip_empty:
            output = "".join(o.decode() for o in output if o != b"\x00")
        else:
            output = "".join(o.decode() for o in output)
    elif dtype == int:
        ssize = struct.calcsize("i") * num_values
        output = np.asarray(struct.unpack("i" * num_values, f.read(ssize)))
        if len(output) == 1:
            output = int(np.squeeze(output))
    elif dtype == float:
        ssize = struct.calcsize("f") * num_values
        output = np.asarray(struct.unpack("f" * num_values, f.read(ssize)))
        if len(output) == 1:
            output = float(np.squeeze(output))
    else:
        raise TypeError(f"datatype '{dtype}' was not valid")

    # shift cursor if need be (e.g. if a reserved section follows)
    if cursor_shift:
        f.seek(cursor_shift, 1)
    return output
