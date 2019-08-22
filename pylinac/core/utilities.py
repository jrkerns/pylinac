"""Utility functions for pylinac."""
from collections import Iterable
import decimal
import os
import os.path as osp
import subprocess
import struct
from typing import Union, Sequence
from datetime import datetime

import pydicom
import numpy as np
from sklearn.preprocessing import MinMaxScaler

from .typing import NumberLike


def clear_data_files():
    """Delete all demo files, image classifiers, etc from the demo folder"""
    demo_folder = osp.join(osp.dirname(osp.dirname(__file__)), 'demo_files')
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


def is_close(val: NumberLike, target: Union[NumberLike, Sequence], delta: NumberLike=1):
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


def typed_property(name, expected_type_or_tuple_of_types):
    """Type-enforced property. Python Cookbook 9.21 (3rd ed)."""
    storage_name = '_' + name

    @property
    def prop(self):
        return getattr(self, storage_name, None)

    @prop.setter
    def prop(self, value):
        if not isinstance(value, expected_type_or_tuple_of_types):
            raise TypeError(f"{name} must be a {expected_type_or_tuple_of_types}. Got: {type(value)}")
        setattr(self, storage_name, value)

    return prop


def simple_round(number: NumberLike, decimals: int=0):
    """Round a number to the given number of decimals. Fixes small floating number errors."""
    num = int(round(number * 10 ** decimals))
    num /= 10 ** decimals
    return num


def isnumeric(object):
    """Check whether the passed object is numeric in any sense."""
    return isinstance(object, (int, float, decimal.Decimal, np.number))


def is_float_like(number):
    return isinstance(number, (float, np.float, np.float16, np.float32, np.float64))


def is_int_like(number):
    return isinstance(number, (int, np.int, np.int16, np.int32, np.int64, np.int8))


def is_iterable(object):
    """Determine if an object is iterable."""
    return isinstance(object, Iterable)


def minmax_scale(array, feature_range=(0, 1), axis=0, copy=True):
    """Copy of scikit-learn's minmax_scale function. Reproduced here for backwards compatibility."""
    original_ndim = array.ndim

    if original_ndim == 1:
        array = array.reshape(array.shape[0], 1)

    s = MinMaxScaler(feature_range=feature_range, copy=copy)
    if axis == 0:
        array = s.fit_transform(array)
    else:
        array = s.fit_transform(array.T).T

    if original_ndim == 1:
        array = array.ravel()

    return array


class Structure:
    """A simple structure that assigns the arguments to the object."""
    def __init__(self, **kwargs):
        self.__dict__.update(**kwargs)

    def update(self, **kwargs):
        self.__dict__.update(**kwargs)


def decode_binary(file, dtype, num_values=1, cursor_shift=0):
    """Read in a raw binary file and convert it to given data types.

    Parameters
    ----------
    file : file object
        The open file object.
    dtype : int, float, str
        The expected data type to return. If int or float, will return numpy array.
    num_values : int
        The expected number of dtype to return

        .. note:: This is not the same as the number of bytes.

    cursor_shift : int
        The number of bytes to move the cursor forward after decoding. This is used if there is a
        reserved section after the read-in segment.
    """
    f = file

    if dtype == str:  # if string
        output = f.read(num_values)
        if type(f) is not str:  # in py3 fc will be bytes
            output = output.decode()
        # strip the padding ("\x00")
        output = output.strip('\x00')
    elif dtype == int:
        ssize = struct.calcsize('i') * num_values
        output = np.asarray(struct.unpack('i' * num_values, f.read(ssize)))
        if len(output) == 1:
            output = int(output)
    elif dtype == float:
        ssize = struct.calcsize('f') * num_values
        output = np.asarray(struct.unpack('f' * num_values, f.read(ssize)))
        if len(output) == 1:
            output = float(output)
    else:
        raise TypeError(f"datatype '{dtype}' was not valid")

    # shift cursor if need be (e.g. if a reserved section follows)
    if cursor_shift:
        f.seek(cursor_shift, 1)
    return output


def open_path(path: str):
    """Open the specified path in the system default viewer."""

    if os.name == 'darwin':
        launcher = "open"
    elif os.name == 'posix':
        launcher = "xdg-open"
    elif os.name == 'nt':
        launcher = "explorer"
    subprocess.call([launcher, path])


def file_exists(filename: str):
    """Check if the file exists and if it does add a timestamp"""
    if osp.exists(filename):
        filename, ext = osp.splitext(filename)
        mytime = datetime.now().strftime("%Y%m%d%H%M%S")
        filename = filename + mytime + ext
    return filename
