"""Utility functions for pylinac."""

from __future__ import annotations

import json
import os
import os.path as osp
import struct
from abc import abstractmethod
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import BinaryIO, Generic, Literal, TypeVar

import numpy as np
import pydicom
from pydantic import BaseModel, ConfigDict, Field
from quaac import Attachment, DataPoint, Document, Equipment, User

from .. import __version__, version
from .scale import wrap360
from .warnings import WarningCollectorMixin


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
    pylinac_version: str = Field(
        default=__version__,
        title="Pylinac version",
        description="The version of Pylinac used for the analysis.",
    )
    date_of_analysis: datetime = Field(
        default_factory=datetime.today,
        title="Date of Analysis",
        description="The date the analysis was performed.",
    )
    warnings: list[dict] = Field(
        title="Warnings",
        description="Code warnings that occurred during the analysis.",
        default_factory=list,
    )


T = TypeVar("T")


class ResultsDataMixin(Generic[T], WarningCollectorMixin):
    """A mixin for classes that generate results data. This mixin is used to generate the results data and present it in different formats.
    The generic types allow correct type hinting of the results data."""

    @abstractmethod
    def _generate_results_data(self) -> T:
        pass

    def results_data(
        self,
        as_dict: bool = False,
        as_json: bool = False,
        by_alias: bool = False,
        exclude: set[str] | None = None,
    ) -> T | dict | str:
        """Present the results data and metadata as a dataclass, dict, or tuple.
        The default return type is a dataclass.

        Parameters
        ----------
        as_dict : bool
            If True, return the results as a dictionary.
        as_json : bool
            If True, return the results as a JSON string. Cannot be True if as_dict is True.
        by_alias : bool
            If True, use the alias names of the dataclass fields. These are generally the more human-readable names.
        exclude : set
            A set of fields to exclude from the results data.
        """
        if as_dict and as_json:
            raise ValueError("Cannot return as both dict and JSON. Pick one.")
        data = self._generate_results_data()
        if as_dict:
            return json.loads(data.model_dump_json(by_alias=by_alias, exclude=exclude))
        if as_json:
            return data.model_dump_json(by_alias=by_alias, exclude=exclude)
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


def is_close_degrees(angle1: float, angle2: float, delta: float = 1) -> bool:
    """A sister function to is_close that takes into account the circular nature of degrees.

    Parameters
    ----------
    angle1 : float
        The first angle in degrees.
    angle2 : float
        The second angle in degrees.
    delta : float
        The maximum difference allowed between the angles in degrees
    """
    if delta < 0:
        raise ValueError("Delta must be positive")
    angle1 = wrap360(angle1)
    angle2 = wrap360(angle2)
    simple_diff = abs(angle1 - angle2)
    other_side_of_circle = 360 - simple_diff
    return min(simple_diff, other_side_of_circle) <= delta


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


class TemporaryAttribute:
    """Context manager to temporarily set a class attribute."""

    def __init__(self, cls, attribute_name, temporary_value):
        self.cls = cls
        self.attribute_name = attribute_name
        self.temporary_value = temporary_value
        self.original_value = getattr(cls, attribute_name)

    def __enter__(self):
        setattr(self.cls, self.attribute_name, self.temporary_value)

    def __exit__(self, exc_type, exc_value, traceback):
        setattr(self.cls, self.attribute_name, self.original_value)


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
    elif dtype is str:
        ssize = struct.calcsize("c") * num_values
        output = struct.unpack("c" * num_values, f.read(ssize))
        if strip_empty:
            output = "".join(o.decode() for o in output if o != b"\x00")
        else:
            output = "".join(o.decode() for o in output)
    elif dtype is int:
        ssize = struct.calcsize("i") * num_values
        output = np.asarray(struct.unpack("i" * num_values, f.read(ssize)))
        if len(output) == 1:
            output = int(np.squeeze(output))
    elif dtype is float:
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


@dataclass  # dataclasses can have default values; typed dicts cannot
class QuaacDatum:
    """Individual data point to be saved to a QuAAC file."""

    value: str | float | int
    unit: str = ""
    description: str = ""
    reference_value: str | float | int | None = None


class QuaacMixin:
    """A mixin for pylinac analysis classes to save results to a QuAAC file."""

    @abstractmethod
    def _quaac_datapoints(self) -> dict[str, QuaacDatum]:
        """Return the data points to be saved to the QuAAC file. The tuple is in the format of
        (name, value, unit, description)."""
        raise NotImplementedError

    def to_quaac(
        self,
        path: str | Path,
        performer: User,
        primary_equipment: Equipment,
        format: Literal["json", "yaml"] = "yaml",
        attachments: list[Attachment] | None = None,
        overwrite: bool = False,
        **kwargs,
    ) -> None:
        """Write an analysis to a QuAAC file. This will include the items
        from results_data() and the PDF report.

        Parameters
        ----------
        path : str, Path
            The file to write the results to.
        performer : User
            The user who performed the analysis.
        primary_equipment : Equipment
            The equipment used in the analysis.
        format : {'json', 'yaml'}
            The format to write the file in.
        attachments : list of Attachment
            Additional attachments to include in the QuAAC file.
        overwrite : bool
            Whether to overwrite the file if it already exists.
        kwargs
            Additional keyword arguments to pass to the Document instantiation.
        """
        attachments = attachments or []
        if Path(path).exists() and not overwrite:
            raise FileExistsError(
                f"{path} already exists. Pass 'overwrite=True' to overwrite."
            )
        datapoints = []
        data_values = self._quaac_datapoints()
        for name, datum in data_values.items():
            dp = DataPoint(
                performer=performer,
                perform_datetime=datetime.now(),
                primary_equipment=primary_equipment,
                name=name,
                measurement_value=datum.value,
                measurement_unit=datum.unit,
                description=datum.description,
                reference_value=datum.reference_value,
                attachments=attachments,
                parameters={"pylinac version": version.__version__},
            )
            datapoints.append(dp)
        d = Document(datapoints=datapoints, **kwargs)
        if format == "json":
            d.to_json_file(path)
        elif format == "yaml":
            d.to_yaml_file(path)


def uniquify(seq: list[str] | tuple[str], value: str) -> str:
    """Create a unique string if it already exists in the sequence"""
    if value not in seq:
        return value
    i = 1
    while True:
        new_value = f"{value}-{i}"
        if new_value not in seq:
            return new_value
        i += 1
