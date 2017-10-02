"""Utility functions for pylinac."""
from collections import Iterable
import decimal
import os
import os.path as osp
import struct

import dicom
import numpy as np
from sklearn.preprocessing import MinMaxScaler


def clear_data_files():
    """Delete all demo files, image classifiers, etc from the demo folder"""
    demo_folder = osp.join(osp.dirname(osp.dirname(__file__)), 'demo_files')
    if osp.isdir(demo_folder):
        for file in os.listdir(demo_folder):
            full_file = osp.join(demo_folder, file)
            if osp.isfile(full_file):
                os.remove(full_file)
    print("Pylinac data files cleared.")


def assign2machine(source_file, machine_file):
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
    dcm_source = dicom.read_file(source_file)
    dcm_machine = dicom.read_file(machine_file)
    for beam in dcm_source.BeamSequence:
        beam.TreatmentMachineName = dcm_machine.BeamSequence[0].TreatmentMachineName
    dcm_source.save_as(source_file)


def is_close(val, target, delta=1):
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


def import_mpld3():
    """Try importing MPLD3. Raises error if not installed. Returns the MPLD3 library."""
    try:
        import mpld3
    except ImportError:
        raise ImportError("The MPLD3 library must be installed to make interactive plots. See http://mpld3.github.io/index.html for info.")
    return mpld3


def typed_property(name, expected_type_or_tuple_of_types):
    """Type-enforced property. Python Cookbook 9.21 (3rd ed)."""
    storage_name = '_' + name

    @property
    def prop(self):
        return getattr(self, storage_name, None)

    @prop.setter
    def prop(self, value):
        if not isinstance(value, expected_type_or_tuple_of_types):
            raise TypeError("{0} must be a {1}. Got: {2}".format(name, expected_type_or_tuple_of_types, type(value)))
        setattr(self, storage_name, value)

    return prop


def simple_round(number, decimals=0):
    """Round a number to the given number of decimals. Fixes small floating number errors."""
    num = int(round(number * 10 ** decimals))
    num /= 10 ** decimals
    return num


def is_dicom(file):
    """Boolean specifying if file is a proper DICOM file.

    This function is a pared down version of read_preamble meant for a fast return.
    The file is read for a proper preamble ('DICM'), returning True if so,
    and False otherwise. This is a conservative approach.

    Parameters
    ----------
    file : str
        The path to the file.

    See Also
    --------
    pydicom.filereader.read_preamble
    pydicom.filereader.read_partial
    """
    fp = open(file, 'rb')
    preamble = fp.read(0x80)
    prefix = fp.read(4)
    return prefix == b"DICM"


def is_dicom_image(file):
    """Boolean specifying if file is a proper DICOM file with a image

    Parameters
    ----------
    file : str
        The path to the file.

    See Also
    --------
    pydicom.filereader.read_preamble
    pydicom.filereader.read_partial
    """
    result = False
    try:
        img = dicom.read_file(file, force=True)
        img.pixel_array
        result = True
    except:
        pass
    return result


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
        raise TypeError("datatype '{}' was not valid".format(dtype))

    # shift cursor if need be (e.g. if a reserved section follows)
    if cursor_shift:
        f.seek(cursor_shift, 1)
    return output
