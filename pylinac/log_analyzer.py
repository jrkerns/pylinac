from abc import ABCMeta, abstractproperty
import struct
import os.path as osp
import csv
import copy

import numpy as np
import scipy.ndimage.filters as spf
import matplotlib.pyplot as plt

from pylinac.core import io
from pylinac.core.decorators import type_accept, lazyproperty
from pylinac.core.io import is_valid_file


np.seterr(invalid='ignore')  # ignore warnings for invalid numpy operations. Used for np.where() operations on partially-NaN arrays.

"""Named Constants"""
log_types = {'dlog': 'Dynalog', 'tlog': 'Trajectory log'}  # the two current log types
dlg_file_exts = ('.dlg',)  # add more if need be
tlg_file_exts = ('.bin',)  # ditto
dynalog_leaf_conversion = 1.96614  # MLC physical plane scaling factor to iso (100cm SAD) plane


class Axis:
    """Represents an 'Axis' of a Trajectory log or dynalog file, holding actual and possibly expected values."""
    def __init__(self, actual, expected=None):
        """
        Parameters
        ----------
        actual : numpy.ndarray
            The array of actual position values.
        expected : numpy.ndarray, optional
            The array of expected position values. Not applicable for dynalog axes other than MLCs.
        """
        self.actual = actual
        if expected is not None:
            if len(actual) != len(expected):
                raise ValueError("Actual and expected Axis parameters are not equal length")
            self.expected = expected

    @lazyproperty
    def difference(self):
        """Return an array of the difference between actual and expected positions.

        Returns
        -------
        numpy.ndarray
            Array the same length as actual/expected.
        """
        if self.expected is not None:
            return self.actual - self.expected
        else:
            raise AttributeError("Expected positions not passed to Axis")

    def plot_actual(self):
        """Plot the actual positions as a matplotlib figure."""
        plt.plot(self.actual)
        plt.show()

    def plot_expected(self):
        """Plot the expected positions as a matplotlib figure."""
        plt.plot(self.expected)
        plt.show()

    def plot_difference(self):
        """Plot the difference of positions as a matplotlib figure."""
        plt.plot(self.difference)
        plt.show()


class LeafAxis(Axis):
    """Axis holding leaf information."""
    def __init__(self, actual, expected):
        # force expected argument to be supplied
        super().__init__(actual, expected)


class GantryAxis(Axis):
    """Axis holding gantry information."""
    pass


class HeadAxis(Axis):
    """Axis holding head information (e.g. jaw positions, collimator)."""
    pass


class CouchAxis(Axis):
    """Axis holding couch information."""
    pass


class BeamAxis(Axis):
    """Axis holding beam information (e.g. MU, beam hold status)."""
    pass


class Fluence(metaclass=ABCMeta):
    """An abstract class to be used for the actual and expected fluences."""
    pixel_map = np.ndarray
    resolution = -1
    fluence_type = ''  # must be specified by subclass

    def __init__(self, mlc_struct=None, mu_axis=None):
        """
        Parameters
        ----------
        mlc_struct : MLC_Struct
        mu_axis : BeamAxis
        """
        self.mlc = mlc_struct
        self.mu = mu_axis

    @property
    def map_calced(self):
        """Return a boolean specifying whether the fluence has been calculated."""
        if self.pixel_map.size != 0:
            return True
        else:
            return False

    def _same_conditions(self, resolution):
        """Return whether the conditions passed are the same as prior conditions (for semi-lazy operations)."""
        if self.resolution != resolution:
            return False
        else:
            return True

    def calc_map(self, resolution=1):
        """Calculate a fluence pixel map.

        Parameters
        ----------
        resolution : float
            The resolution in mm of the fluence calculation in the leaf-moving direction.

         Returns
         -------
         numpy.ndarray
             A numpy array reconstructing the actual fluence of the log. The size will
             be the number of MLC pairs by 400 / resolution since the MLCs can move anywhere within the
             40cm-wide linac head opening.
         """
        # return if map has already been calculated under the same conditions
        if self.map_calced and self._same_conditions(resolution):
            return self.pixel_map
        # preallocate arrays for expected and actual fluence of number of leaf pairs-x-4000 (40cm = 4000um, etc)
        fluence = np.zeros((self.mlc.num_pairs, 400 / resolution), dtype=float)

        # calculate the MU delivered in each snapshot. For Tlogs this is absolute; for dynalogs it's normalized.
        mu_matrix = getattr(self.mu, self.fluence_type)
        MU_differential = np.zeros(len(mu_matrix))
        MU_differential[0] = mu_matrix[0]
        MU_differential[1:] = np.diff(mu_matrix)
        MU_differential = MU_differential / mu_matrix[-1]
        MU_cumulative = mu_matrix[-1]

        # calculate each "line" of fluence (the fluence of an MLC leaf pair, e.g. 1 & 61, 2 & 62, etc),
        # and add each "line" to the total fluence matrix
        fluence_line = np.zeros((400 / resolution))
        leaf_offset = self.mlc.num_pairs
        pos_offset = int(np.round(200 / resolution))
        for leaf in range(1, self.mlc.num_pairs + 1):
            if self.mlc.leaf_is_under_y_jaw(leaf):
                continue
            else:
                fluence_line[:] = 0  # emtpy the line values on each new leaf pair
                left_leaf_data = getattr(self.mlc.leaf_axes[leaf], self.fluence_type)
                left_leaf_data = -np.round(left_leaf_data * 10 / resolution) + pos_offset
                right_leaf_data = getattr(self.mlc.leaf_axes[leaf + leaf_offset], self.fluence_type)
                right_leaf_data = np.round(right_leaf_data * 10 / resolution) + pos_offset
                if self.mlc.leaf_did_move(leaf):
                    for snapshot in self.mlc.snapshot_idx:
                        # TODO: incorporate numexpr and multithreading
                        left_pos = left_leaf_data[snapshot]
                        right_pos = right_leaf_data[snapshot]
                        fluence_line[left_pos:right_pos] += MU_differential[snapshot]
                else:  # leaf didn't move; no need to calc over every snapshot
                    first_snapshot = self.mlc.snapshot_idx[0]
                    left_pos = left_leaf_data[first_snapshot]
                    right_pos = right_leaf_data[first_snapshot]
                    fluence_line[left_pos:right_pos] = MU_cumulative
                fluence[leaf - 1, :] = fluence_line

        self.pixel_map = fluence
        self.resolution = resolution
        return fluence

    def plot_map(self):
        """Plot the fluence; the pixel map must have been calculated first."""
        if not self.map_calced:
            raise AttributeError("Map not yet calculated; use calc_map()")
        plt.imshow(self.pixel_map, aspect='auto')
        plt.show()


class ActualFluence(Fluence):
    """The actual fluence object."""
    fluence_type = 'actual'


class ExpectedFluence(Fluence):
    """The expected fluence object."""
    fluence_type = 'expected'


class GammaFluence(Fluence):
    """Gamma object, including pixel maps of gamma, binary pass/fail pixel map, and others."""
    distTA = -1
    doseTA = -1
    threshold = -1
    pass_prcnt = -1
    avg_gamma = -1
    doseTA_map = np.ndarray
    distTA_map = np.ndarray
    passfail_map = np.ndarray

    def __init__(self, actual_fluence, expected_fluence, mlc_struct):
        """
        Parameters
        ----------
        actual_fluence : ActualFluence
            The actual fluence object.
        expected_fluence : ExpectedFluence
            The expected fluence object.
        mlc_struct : MLC_Struct
            The MLC structure, so fluence can be calculated from leaf positions.
        """
        self.actual_fluence = actual_fluence
        self.expected_fluence = expected_fluence
        self.mlc = mlc_struct

    def _same_conditions(self, doseTA, distTA, threshold, resolution):
        """Determine if the passed conditions are the same as the existing ones.

        See calc_map for parameter info.
        """
        if (self.distTA == distTA & self.doseTA == doseTA
            & self.threshold == threshold & self.resolution == resolution):
            return True
        else:
            return False

    def calc_map(self, DoseTA=2, DistTA=1, threshold=10, resolution=0.1, calc_individual_maps=False):
        """Calculate the gamma from the actual and expected fluences.

        The gamma calculation is based on `Bakai et al
        <http://iopscience.iop.org/0031-9155/48/21/006/>`_ eq.6,
        which is a quicker alternative to the standard Low gamma equation.

        Parameters
        ----------
        DoseTA : int, float
            Dose-to-agreement in percent; e.g. 2 is 2%.
        DistTA : int, float
            Distance-to-agreement in mm
        threshold : int, float
            The dose threshold percentage of the maximum dose, below which is not analyzed.
        resolution : int, float
            The resolution in mm of the resulting gamma map in the leaf-movement direction

        Returns
        -------
        numpy.ndarray
            A num_mlc_leaves-x-400/resolution numpy array.
        """
        # if gamma has been calculated under same conditions, return it
        if self.map_calced and self._same_conditions(DistTA, DoseTA, threshold, resolution):
            return self.pixel_map

        # calc fluences if need be
        if not self.actual_fluence.map_calced or resolution != self.actual_fluence.resolution:
            self.actual_fluence.calc_map(resolution)
        if not self.expected_fluence.map_calced or resolution != self.expected_fluence.resolution:
            self.expected_fluence.calc_map(resolution)

        actual = copy.copy(self.actual_fluence.pixel_map)
        expected = copy.copy(self.expected_fluence.pixel_map)

        # set dose values below threshold to 0 so gamma doesn't calculate over it
        actual[actual < (threshold / 100) * np.max(actual)] = np.nan
        expected[expected < (threshold / 100) * np.max(expected)] = np.nan

        # preallocation
        gamma_map = np.zeros(expected.shape)

        # image gradient in x-direction (leaf movement direction) using sobel filter
        img_x = spf.sobel(actual, 1)

        # equation: (measurement - reference) / sqrt ( doseTA^2 + distTA^2 * image_gradient^2 )
        for leaf in range(self.mlc.num_pairs):
            gamma_map[leaf] = np.abs(actual[leaf, :] - expected[leaf, :]) / np.sqrt(
                (DoseTA / 100.0 ** 2) + ((DistTA / resolution ** 2) * (img_x[leaf, :] ** 2)))
        # construct binary pass/fail map
        self.passfail_map = np.array(np.where(gamma_map >= 1, 1, 0)[0], dtype=bool)

        if calc_individual_maps:
            # calculate DoseTA map (drops distTA calc from gamma eq)
            self.doseTA_map = np.zeros(gamma_map.shape)
            for leaf in range(self.mlc.num_pairs):
                self.doseTA_map[leaf] = (actual[leaf, :] - expected[leaf, :]) / np.sqrt(
                    (DoseTA / 100.0 ** 2))
            # calculate DistTA map (drops DoseTA calc from gamma eq)
            self.distTA_map = np.zeros(gamma_map.shape)
            for leaf in range(self.mlc.num_pairs):
                self.distTA_map[leaf] = (actual[leaf, :] - expected[leaf, :]) / np.sqrt(
                    ((DistTA / resolution ** 2) * (img_x[leaf, :] ** 2)))

        # calculate standard metrics
        self.pass_prcnt = (np.sum(gamma_map < 1) / np.sum(gamma_map >= 0)) * 100
        self.avg_gamma = np.nanmean(gamma_map)

        self.pixel_map = gamma_map
        return gamma_map

    def plot_passfail_map(self):
        """Plot the binary gamma map, only showing whether pixels passed or failed."""
        plt.imshow(self.passfail_map)
        plt.show()


class Fluence_Struct:
    """Structure for data and methods having to do with fluences."""
    def __init__(self, mlc_struct=None, mu_axis=None):
        self.actual = ActualFluence(mlc_struct, mu_axis)
        self.expected = ExpectedFluence(mlc_struct, mu_axis)
        self.gamma = GammaFluence(self.actual, self.expected, mlc_struct)


class MLC:
    """The MLC class holds MLC information and retreives relevant data about the MLCs and positions."""
    def __init__(self, snapshot_idx=None, jaw_struct=None, HDMLC=False):
        """
        Parameters
        ----------
        snapshot_idx : array, list
            The snapshots to be considered for RMS and error calculations (can be all snapshots or just when beam was on).
        jaw_struct : Jaw_Struct
        HDMLC : boolean
            If False (default), indicates a regular MLC model (e.g. Millenium 120).
            If True, indicates an HD MLC model (e.g. Millenium 120 HD).
        """
        super().__init__()
        self.leaf_axes = {}
        self.snapshot_idx = snapshot_idx
        self.jaws = jaw_struct
        self.hdmlc = HDMLC

    @type_accept(leaf_axis=LeafAxis, leaf_num=int)
    def add_leaf_axis(self, leaf_axis, leaf_num):
        """Add a leaf axis to the MLC data structure.

        Parameters
        ----------
        leaf_axis : LeafAxis
            The leaf axis to be added.
        leaf_num : int
            The leaf number; must be between 1 and 120. 1-60 corresponds to the A-bank; 61-120 to the B-bank.
        """
        self.leaf_axes[leaf_num] = leaf_axis

    @lazyproperty
    def num_pairs(self):
        """Return the number of MLC pairs."""
        return int(self.num_leaves/2)

    @lazyproperty
    def num_leaves(self):
        """Return the number of MLC leaves."""
        return len(self.leaf_axes)

    @lazyproperty
    def num_snapshots(self):
        return len(self.snapshot_idx)

    @lazyproperty
    def moving_leaves(self):
        """Return an tuple holding the leaf numbers that moved during treatment."""
        threshold = 0.003
        indices = ()
        for leaf_num, leafdata in self.leaf_axes.items():
            leaf_stdev = np.std(leafdata.actual[self.snapshot_idx])
            if leaf_stdev > threshold:
                indices += (leaf_num,)
        return np.array(indices)

    def leaf_did_move(self, leaf_num):
        """Return whether the leaf moved during treatment."""
        if leaf_num in self.moving_leaves:
            return True
        else:
            return False

    @lazyproperty
    def all_leaf_indices(self):
        """Return an array enumerated over the number of leaves."""
        return np.array(range(1, len(self.leaf_axes) + 1))

    def get_RMS_avg(self, bank=None, only_moving_leaves=False):
        """Return the overall average RMS of given leaves."""
        leaves = self.get_leaves(bank, only_moving_leaves)
        rms_array = self.create_RMS_array(leaves)
        return np.mean(rms_array)

    def get_RMS_max(self, bank=None, only_moving_leaves=False):
        """Return the overall maximum RMS of given leaves."""
        leaves = self.get_leaves(bank, only_moving_leaves)
        rms_array = self.create_RMS_array(leaves)
        return np.max(rms_array)

    def get_leaves(self, bank=None, only_moving_leaves=False):
        """Return a list of leaves that match the given conditions.

        Parameters
        ----------
        bank : {'A', 'B'}, None
            Specifies which bank is desired. If None, will return both banks.
        only_moving_leaves : boolean
            If False (default), include all the leaves.
            If True, will remove the leaves that were static during treatment.

            .. warning::
                The RMS and error will nearly always be lower if all leaves are included since non-moving leaves
                have an error of 0 and will drive down the average values. Convention would include all leaves,
                but prudence would use only the moving leaves to get a more accurate assessment of error/RMS.
        """
        # get all leaves or only the moving leaves
        if only_moving_leaves:
            leaves = self.moving_leaves
        else:
            leaves = self.all_leaf_indices

        # select leaves by bank if desired
        if bank is not None:
            if 'a' in bank.lower():
                leaves = leaves[leaves <= self.num_pairs]
            elif 'b' in bank.lower():
                leaves = leaves[leaves > self.num_pairs]

        return leaves

    def get_nth_perc_error(self, percentile=95, bank=None, only_moving_leaves=False):
        """Calculate the n-th percentile error of the leaf error."""
        leaves = self.get_leaves(bank, only_moving_leaves)
        error_array = self.create_error_array(leaves)

        abs_error = np.abs(error_array)
        return np.percentile(abs_error, percentile)

    def create_error_array(self, leaves):
        """Create an error array of only the leaves specified."""
        error_array = self._error_array_all_leaves
        return error_array[leaves, :]

    def create_RMS_array(self, leaves):
        """Create an RMS array of only the leaves specified."""
        rms_array = self._RMS_array_all_leaves
        leaves -= 1
        return rms_array[leaves]

    @lazyproperty
    def _error_array_all_leaves(self):
        # preallocate
        mlc_error = np.zeros((self.num_leaves, self.num_snapshots))
        # construct numpy array for easy array calculation
        for leaf in range(self.num_leaves):
            mlc_error[leaf, :] = self.leaf_axes[leaf+1].difference[self.snapshot_idx]
        return mlc_error

    @lazyproperty
    def _error_array_moving_leaves(self):
        """Return an error array only for leaves that moved."""
        moving_leaves = self.moving_leaves
        error_array = self._error_array_all_leaves
        return error_array[moving_leaves, :]

    @lazyproperty
    def _RMS_array_all_leaves(self):
        """Return the RMS of all leaves."""
        rms_array = np.array([np.sqrt(np.sum(leafdata.difference[self.snapshot_idx] ** 2) / self.num_snapshots) for leafdata in self.leaf_axes.values()])
        return rms_array

    @lazyproperty
    def _RMS_array_moving_leaves(self):
        """Return an array of RMS of only the moving leaves."""
        moving_leaves = self.moving_leaves
        rms_array = self._RMS_array_all_leaves
        return rms_array[moving_leaves, :]

    def leaf_is_under_y_jaw(self, leaf_num):
        """Return a boolean specifying if the given leaf is under one of the y jaws."""
        outer_leaf_thickness = 10  # mm
        inner_leaf_thickness = 5
        mlc_position = 0
        if self.hdmlc:
            outer_leaf_thickness /= 2
            inner_leaf_thickness /= 2
            mlc_position = 100
        for leaf in range(1, leaf_num+1):
            if 10 >= leaf or leaf >= 110:
                mlc_position += outer_leaf_thickness
            elif 50 >= leaf or leaf >= 70:
                mlc_position += inner_leaf_thickness
            else:
                mlc_position += outer_leaf_thickness

        y2_position = self.jaws.y2.actual[0]*10 + 200
        y1_position = 200 - self.jaws.y1.actual[0]*10
        if mlc_position < y1_position or mlc_position - outer_leaf_thickness > y2_position:
            return True
        else:
            return False

    def leaf_is_under_x_jaw(self, leaf_num):
        """Return a boolean specifying if the given leaf is under one of the x jaws."""
        raise NotImplementedError
        # most_retracted_leaf_position = np.max(self.leaf_axes[leaf_num].actual)
        # most_extended_leaf_position = np.min(self.leaf_axes[leaf_num].actual)
        #
        # x2_position = self.jaws.x2.actual[0] * 10 + 200
        # x1_position = 200 - self.jaws.x1.actual[0] * 10
        #
        #
        # if most_retracted_leaf_position < x1_position:
        #     return False
        # else:
        #     return True


class Jaw_Struct:
    """Jaw Axes data structure, holding jaw HeadAxes instances."""
    def __init__(self, x1, y1, x2, y2):
        if not all((
                isinstance(x1, HeadAxis),
                isinstance(y1, HeadAxis),
                isinstance(x2, HeadAxis),
                isinstance(y2, HeadAxis))):
            raise TypeError("HeadAxis not passed into Jaw structure")
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2

class Couch_Struct:
    """Couch Axes data structure, holding CouchAxis instances."""
    def __init__(self, vertical, longitudinal, lateral, rotational):
        if not all((isinstance(vertical, CouchAxis),
                    isinstance(longitudinal, CouchAxis),
                    isinstance(lateral, CouchAxis),
                    isinstance(rotational, CouchAxis))):
            raise TypeError("Couch Axes not passed into Couch structure")
        self.vert = vertical
        self.long = longitudinal
        self.latl = lateral
        self.rotn = rotational


class Log_Section(metaclass=ABCMeta):
    @abstractproperty
    def read(self):
        pass


class DLog_Section(Log_Section, metaclass=ABCMeta):
    def __init__(self, log_content):
        self.log_content = log_content


class TLog_Section(Log_Section, metaclass=ABCMeta):
    def __init__(self, log_content, cursor):
        self.log_content = log_content
        self._cursor = cursor

    def _decode_binary(self, filecontents, dtype, num_values=1, cursor_shift=0):
        """This method is the main "decoder" for reading in trajectory log binary data into human data types.

        Parameters
        ----------
        filecontents : file object
            The complete file having been read with .read().
        dtype : int, float, str
            The expected data type to return. If int or float, will return numpy array.
        num_values : int
            The expected number of dtype to return

            .. note:: This is not the same as the number of bytes.

        cursor_shift : int
            The number of bytes to move the cursor forward after decoding. This is used if there is a
            reserved section after the read-in segment.
        """
        fc = filecontents

        if dtype == str:  # if string
            output = fc[self._cursor:self._cursor + num_values]
            if type(fc) is not str:  # in py3 fc will be bytes
                output = output.decode()
            # Now, strip the padding ("\x00")
            output = output.strip('\x00')
            self._cursor += num_values
        elif dtype == int:
            ssize = struct.calcsize('i') * num_values
            output = np.asarray(struct.unpack('i' * num_values, fc[self._cursor:self._cursor + ssize]))
            self._cursor += ssize
        elif dtype == float:
            ssize = struct.calcsize('f') * num_values
            output = np.asarray(struct.unpack('f' * num_values, fc[self._cursor:self._cursor + ssize]))
            self._cursor += ssize
        else:
            raise TypeError("decode_binary datatype was not valid")

        self._cursor += cursor_shift  # shift cursor if need be (e.g. if a reserved section follows)
        return output

class Subbeam(TLog_Section):
    def __init__(self, log_content, cursor, log_version):
        super().__init__(log_content, cursor)
        self.log_version = log_version

    def read(self):
        self.control_point = self._decode_binary(self.log_content, int)
        self.mu_delivered = self._decode_binary(self.log_content, float)
        self.rad_time = self._decode_binary(self.log_content, float)
        self.sequence_num = self._decode_binary(self.log_content, int)
        # In Tlogs version 3.0 and up, beam names are 512 byte unicode strings, but in <3.0 they are 32 byte unicode strings
        if self.log_version >= 3:
            chars = 512
        else:
            chars = 32
        self.beam_name = self._decode_binary(self.log_content, str, chars, 32)

        return self, self._cursor


class Subbeamer:
    """One of 4 subsections of a trajectory log. Holds a list of subbeams; only applicable for auto-sequenced beams."""
    def __init__(self, log_content, cursor, header):
        self._log_content = log_content
        self._header = header
        self._cursor = cursor

    def read(self):
        subbeams = []
        if self._header.num_subbeams > 0:
            for subbeam_num in range(self._header.num_subbeams):
                subbeam, self._cursor = Subbeam(self._log_content, self._cursor, self._header.version).read()
                subbeams.append(subbeam)
        return subbeams, self._cursor


class Tlog_Header(TLog_Section):
    """A header object, one of 4 "sections" of a trajectory log. Holds sampling interval, version, etc."""

    def read(self):
        self.header = self._decode_binary(self.log_content, str, 16)  # for version 1.5 will be "VOSTL"
        self.version = float(self._decode_binary(self.log_content, str, 16))  # in the format of 2.x or 3.x
        self.header_size = self._decode_binary(self.log_content, int)  # fixed at 1024 in 1.5 specs
        self.sampling_interval = self._decode_binary(self.log_content, int)
        self.num_axes = self._decode_binary(self.log_content, int)
        self.axis_enum = self._decode_binary(self.log_content, int, self.num_axes)
        self.samples_per_axis = self._decode_binary(self.log_content, int, self.num_axes)
        self.num_mlc_leaves = self.samples_per_axis[-1] - 2  # subtract 2 (each carriage counts as an "axis" and must be removed)
        # self._cursor == self.num_axes * 4 # there is a reserved section after samples per axis. this moves it past it.
        self.clinac_scale = self._decode_binary(self.log_content, int)
        self.num_subbeams = self._decode_binary(self.log_content, int)
        self.is_truncated = self._decode_binary(self.log_content, int)
        self.num_snapshots = self._decode_binary(self.log_content, int)
        # the section after MLC model is reserved. Cursor is moved to the end of this reserved section.
        self.mlc_model = self._decode_binary(self.log_content, int, cursor_shift=1024 - (64 + self.num_axes * 8))

        return self, self._cursor

class Dlog_Header(DLog_Section):
    pass

class AxisData(TLog_Section):
    """One of four data structures outline in the Trajectory log file specification.
        Holds information on all Axes measured, like MU, gantry position, and MLC leaf positions."""
    def __init__(self, log_content, cursor, header):
        super().__init__(log_content, cursor)
        self._header = header

    def read(self, exclude_beam_off=True):
        # step size in bytes
        step_size = sum(self._header.samples_per_axis) * 2

        # read in all snapshot data at once, then assign
        snapshot_data = self._decode_binary(self.log_content, float, step_size * self._header.num_snapshots)

        # reshape snapshot data to be a x-by-num_snapshots matrix
        snapshot_data = snapshot_data.reshape(self._header.num_snapshots, -1)

        # collimator
        self.collimator = self._get_axis(snapshot_data, 0, HeadAxis)

        # gantry
        self.gantry = self._get_axis(snapshot_data, 2, GantryAxis)

        # jaws
        jaw_y1 = self._get_axis(snapshot_data, 4, HeadAxis)
        jaw_y2 = self._get_axis(snapshot_data, 6, HeadAxis)
        jaw_x1 = self._get_axis(snapshot_data, 8, HeadAxis)
        jaw_x2 = self._get_axis(snapshot_data, 10, HeadAxis)
        self.jaws = Jaw_Struct(jaw_x1, jaw_y1, jaw_x2, jaw_y2)

        # couch
        vrt = self._get_axis(snapshot_data, 12, CouchAxis)
        lng = self._get_axis(snapshot_data, 14, CouchAxis)
        lat = self._get_axis(snapshot_data, 16, CouchAxis)
        rtn = self._get_axis(snapshot_data, 18, CouchAxis)
        self.couch = Couch_Struct(vrt, lng, lat, rtn)

        # MU
        self.mu = self._get_axis(snapshot_data, 20, BeamAxis)

        # beam hold state
        self.beam_hold = self._get_axis(snapshot_data, 22, BeamAxis)

        # control point
        self.control_point = self._get_axis(snapshot_data, 24, BeamAxis)

        # carriages
        self.carraige_A = self._get_axis(snapshot_data, 26, HeadAxis)
        self.carriage_B = self._get_axis(snapshot_data, 28, HeadAxis)

        # remove snapshots where the beam wasn't on if flag passed
        if exclude_beam_off:
            snapshots = np.where(self.beam_hold.actual == 0)[0]
        else:
            snapshots = list(range(self._header.num_snapshots))

        if self._header.mlc_model == 2:
            hdmlc = False
        else:
            hdmlc = True

        self.mlc = MLC(snapshots, self.jaws, hdmlc)
        for leaf_num in range(1, self._header.num_mlc_leaves + 1):
            leaf_axis = self._get_axis(snapshot_data, 30 + 2 * (leaf_num - 1), LeafAxis)
            self.mlc.add_leaf_axis(leaf_axis, leaf_num)

        return self, self._cursor

    @lazyproperty
    def num_beamholds(self):
        """Return the number of times the beam was held."""
        diffmatrix = np.diff(self.beam_hold.actual)
        num_holds = int(np.sum(diffmatrix == 1))
        return num_holds


    def _get_axis(self, snapshot_data, column, axis_type):
        """Return column of data from snapshot data of the axis type passed.

        Parameters
        ----------
        snapshot_data : numpy.ndarray
            The data read in holding the axis data of the log.
        column : int
            The column of the desired data in snapshot_data
        axis_type : subclass of Axis
            The type of axis the data is.

        Returns
        -------
        axis_type
        """
        return axis_type(expected=snapshot_data[:, column],
                         actual=snapshot_data[:, column + 1])

class CRC(TLog_Section):
    """The last data section of a Trajectory log. Is a 2 byte cyclic redundancy check (CRC), specifically
        a CRC-16-CCITT. The seed is OxFFFF."""
    def __init__(self, log_content, cursor):
        super().__init__(log_content, cursor)

    def read(self):
        # crc = self._decode_binary(self.log_content, str, 2)
        # TODO: figure this out
        pass


class MachineLog:
    """Reads in and analyzes MLC log files from Varian linear accelerators.

    Attributes
    ----------
    log_type : str
        The log type loaded; either 'Dynalog' or 'Trajectory log'
    """
    @type_accept(filename=str)
    def __init__(self, filename=''):
        """
        Parameters
        ----------
        filename : str
            Path to the log file. For trajectory logs this is a single .bin file.
            For dynalog files this should be the A-file. The B-file will be automatically pulled when A is read in.
            The B-file must be in the same directory as the A-file or an error is thrown.
            If filename is not passed in on init, it will need to be loaded later.
        """
        self._filename = filename
        self._cursor = 0

        #------------------------
        # Generic Log attributes
        #------------------------
        # For dynalogs:   . For Tlogs: usually in 2.x TODO: find dynalog version info
        # self.version = ''
        # # the number of leaves in a log file. For dynalogs it is actually the # of PAIRS.
        # self.num_mlc_leaves = None
        # # For dynalogs: 0->Varian scale, 1->IEC scale. For Tlogs: 1-> machine scale, 2-> modified IEC
        # self.clinac_scale = None
        # # The number of "snapshots" taken during the treatment. Each snapshot captures numerous mechanical parameters
        # self.num_snapshots = None
        # # The sampling interval (i.e. time between snapshots) in ms. For dynalogs, this is 50ms. Tlogs are usu. ~20ms.
        # self.sampling_interval = None
        #
        # #----------------------------
        # # Dynalog-specific attributes
        # #----------------------------
        # # patient name, up to 25 characters
        # self.patient_name = ''
        # self.plan_filename = ''
        # # The mechanical tolerances. Machine will beam-hold if error is greater than tolerance. (mm)
        # self.tolerance = None
        # #----------------------------
        # # Trajectory log-specific attributes
        # #----------------------------
        # # 2->NDS 120 (Millennium), 3->NDS 120 HD
        # self.mlc_model = None
        # self.num_axes = None  # # of axes sampled
        # self.axis_enum = None
        # self.samples_per_axis = None
        # self.is_truncated = None  # 1-> truncated, 0 -> not truncated

        self.fluence = Fluence_Struct()

        # Read file if passed in
        if filename is not '':
            self.read_log()

    def run_tlog_demo(self):
        """Run the Trajectory log demo."""
        self.load_demo_trajectorylog()
        self.report_basic_parameters()
        self.fluence.gamma.plot_map()

    def run_dlog_demo(self):
        """Run the dynalog demo."""
        self.load_demo_dynalog()
        self.report_basic_parameters()
        self.fluence.gamma.plot_map()

    def load_demo_dynalog(self):
        """Load the demo dynalog file included with the package."""
        self._filename = osp.join(osp.dirname(__file__), 'demo_files', 'log_reader', 'AQA.dlg')
        self.read_log()

    def load_demo_trajectorylog(self):
        """Load the demo trajectory log included with the package."""
        self._filename = osp.join(osp.dirname(__file__), 'demo_files', 'log_reader', 'Tlog.bin')
        self.read_log()

    def load_logfile_UI(self):
        """Let user load a log file with a UI dialog box. """
        filename = io.get_filepath_UI()
        if filename: # if user didn't hit cancel...
            self._filename = filename
            self.read_log()

    @type_accept(filename=str)
    def load_logfile(self, filename):
        """Load the log file directly by passing the path to the file.

        Parameters
        ----------
        filename : str
            The path to the log file.
        """
        is_valid_file(filename, raise_error=True)
        self._filename = filename
        self.read_log()

    def report_basic_parameters(self):
        """Print the common parameters analyzed when investigating machine logs:

        -log type
        -average RMS
        -maximum RMS
        -95th percentile error
        -number of beam holdoffs
        -gamma pass percentage
        """
        print("MLC log type: " + self.log_type)
        print("Average RMS of all leaves: {:3.3f} mm".format(self.axis_data.mlc.get_RMS_avg(only_moving_leaves=True)))
        print("Max RMS error of all leaves: {:3.3f} mm".format(self.axis_data.mlc.get_RMS_max(only_moving_leaves=True)))
        print("95th percentile error: {:3.3f} mm".format(self.axis_data.mlc.get_nth_perc_error(95, only_moving_leaves=True)))
        print("Number of beam holdoffs: {:1.0f}".format(self.axis_data.num_beamholds))
        self.fluence.gamma.calc_map()
        print("Gamma pass %: {:2.2f}".format(self.fluence.gamma.pass_prcnt))
        print("Gamma average: {:2.3f}".format(self.fluence.gamma.avg_gamma))


    # @lazyproperty
    # def beam_on_snapshots(self):
    #     """Return the indices of the snapshots only when the beam was actually on.
    #
    #     For dynalogs this removes snapshots where the Beam On flag was 0 and where the Beam Hold was 0.
    #     For trajectory logs this removes the snapshots were the Beam Hold was 0 (there is no Beam On flag).
    #     """
    #     if self.log_type == log_types['tlog']:
    #         return np.where(self.beam_hold.actual == 0)[0]
    #     elif self.log_type == log_types['dlog']:
    #         holdidx = self.beam_hold.actual == 0
    #         beamonidx = self.beam_on.actual == 1
    #         return holdidx & beamonidx

    @lazyproperty
    def log_type(self):
        """Determine the MLC log type: Trajectory or Dynalog.

        First the file extension is examined and assumed. Failing that, it's opened and sample
        the first few bytes. If the sample matches what is found in standard dynalog or tlog files it will be
        assigned that log type.
        """
        __, ext = osp.splitext(self._filename)
        if ext == tlg_file_exts:
            log_type = log_types['tlog'] # trajectory log
        elif ext == dlg_file_exts:
            log_type = log_types['dlog'] # dynalog
        else:
            with open(self._filename, 'rb') as unknown_file:
                header_sample = unknown_file.read(5).decode()
                if 'B' in header_sample or 'A' in header_sample:
                    log_type = log_types['dlog']  # dynalog
                elif 'V' in header_sample:
                    log_type = log_types['tlog']  # trajectory log
                else:
                    raise ValueError("Log type unknown")

        return log_type

    @property
    def log_is_loaded(self):
        """Boolean specifying if a log has been loaded in yet."""
        if self._filename == '':
            return False
        else:
            return True

    def read_log(self):
        """Read in log based on what type of log it is: Trajectory or Dynalog."""
        if not self.log_is_loaded:
            raise AttributeError('Log file has not been specified. Use load_logfile_UI or load_logfile')

        # read log as appropriate to type
        if self.log_type == log_types['tlog']:
            self.read_Tlog()
        elif self.log_type == log_types['dlog']:
            self.read_Dlog()

    def _scale_dlog_mlc_pos(self):
        """Convert MLC leaf plane positions to isoplane positions.

        The function applies a scaling factor
        """
        # Units are initially in 100ths of mm, but are converted to mm.
        for leaf in range(1, self.mlc.num_leaves+1):
            self.mlc.leaf_axes[leaf].actual *= dynalog_leaf_conversion / 100
            self.mlc.leaf_axes[leaf].expected *= dynalog_leaf_conversion / 100

    def read_Dlog(self):
        """Read in Dynalog files from .dlg files (which are renamed CSV files).
        Formatting follows from the Dynalog File Viewer Reference Guide.
        """
        # if file is B-file, skip over file
        if self._filename[0] == 'B':
            return

        # Check that the B-file is in the same folder before getting too far
        bfile = return_B_file(self._filename)

        # create iterator object to read in lines
        with open(self._filename) as csvf:
            dlgdata = csv.reader(csvf, delimiter=',')
            self.version = next(dlgdata)
            self.patient_name = next(dlgdata)
            self.plan_filename = next(dlgdata)
            self.tolerance = int(next(dlgdata)[0])
            self.num_mlc_leaves = int(next(dlgdata)[0]) * 2  # the # of leaves in a dynalog is actually the # of *PAIRS*, hence the *2.
            self.clinac_scale = next(dlgdata)  # 0->Varian scale, 1->IEC scale

            # From here on out, each line is a "snapshot".
            # This reads in all snapshots, then assigns them. The reason this cannot be iterated over is because the
            # snapshots are column-order, not row order, so all rows must be read before any data can be assigned.
            matrix = np.array([line for line in dlgdata], dtype=float)

        self.num_snapshots = np.size(matrix, 0)

        # assignment of snapshot values
        # There is no "expected" MU in dynalogs, but for fluence calc purposes, it is set to that of the actual
        self.mu = Axis(matrix[:,0], matrix[:,0])
        self.DVA_segment = Axis(matrix[:,1])
        self.beam_hold = Axis(matrix[:,2])
        self.beam_on = Axis(matrix[:,3])
        self.prior_dose_idx = matrix[:, 4]  # currently not used for anything
        self.next_dose_idx = matrix[:, 5]  # ditto
        self.gantry = Axis(matrix[:, 6])
        self.collimator = Axis(matrix[:, 7])
        self.jaw_y1 = Axis(matrix[:, 8])
        self.jaw_y2 = Axis(matrix[:, 9])
        self.jaw_x1 = Axis(matrix[:, 10])
        self.jaw_x2 = Axis(matrix[:, 11])
        self.carriageA = Axis(matrix[:, 12])
        self.carriage_B = Axis(matrix[:, 13])
        # MLC positions are in hundredths of mm in the physical leaf plane. Positive is retracted, negative is extented.
        # Positions will be scaled to isocenter plane after bank B positions are added to matrix.

        self.mlc = MLC(self.beam_on.actual)
        for leaf in range(1, (self.num_mlc_leaves//2)+1):
            self.mlc[leaf] = Axis(expected=matrix[leaf-1, 14::4],
                                  actual=matrix[leaf-1, 15::4])
        # self._mlc_expected[:60,:] = matrix[:, 14::4].transpose()
        # self._mlc_actual[:60,:] = matrix[:, 15::4].transpose()

        # read in "B"-file to get bank B MLC positions. The file must be in the same folder as the "A"-file.
        # The header info is repeated but we already have that.
        with open(bfile) as csvf:
            dlgdata = csv.reader(csvf, delimiter=',')
            matrix = np.array([line for line in dlgdata if int(dlgdata.line_num) >= 7], dtype=float)

        # Add bank B MLC positions to mlc snapshot arrays
        for leaf in range(1, (self.num_mlc_leaves//2)+1):
            self.mlc[leaf+self.num_mlc_leaves//2] = Axis(expected=matrix[leaf-1, 14::4],
                                                         actual=matrix[leaf-1, 15::4])
            # self._mlc_expected[60:, :] = matrix[:, 14::4].transpose()
            # self._mlc_actual[60:, :] = matrix[:, 15::4].transpose()

        self._scale_dlog_mlc_pos()

    def read_Tlog(self, exclude_beam_off=True):
        """Read in Trajectory log from binary file according to TB 1.5/2.0 (i.e. Tlog v2.0/3.0) log file specifications.

        Parameters
        ----------
        exclude_beam_off : boolean
            If True (default), snapshots where the beam was not on will be removed.
            If False, all data will be included.

        .. warning::
            Including beam off data may affect fluence and gamma results. E.g. in a step-&-shoot IMRT
            delivery, leaves move between segments while the beam is held. If beam-off data is included,
            this may cause the RMS and fluence errors to rise unnecessarily. Be careful about changing
            this flag.
        """

        # read in trajectory log binary data
        fcontent = open(self._filename, 'rb').read()

        # Unpack the content according to respective section and data type (see log specification file).
        self.header, self._cursor = Tlog_Header(fcontent, self._cursor).read()

        # read in subbeam data. These are for auto-sequenced beams. If not autosequenced, separate logs are created.
        # Currently there is no good way of dealing with this data, but fortunately autosequencing is rare at this time.
        self.subbeams, self._cursor = Subbeamer(fcontent, self._cursor, self.header).read()

        self.axis_data, self._cursor = AxisData(fcontent, self._cursor, self.header).read(exclude_beam_off)

        self.crc = CRC(fcontent, self._cursor).read()

        self.fluence = Fluence_Struct(self.axis_data.mlc, self.axis_data.mu)


def return_B_file(a_filename):
    """Checks that the dynalog "B-file" for a given "A-file" exists within the same directory.

    Returns the absolute locations of the B-file.
    """
    bfile = "B" + osp.split(a_filename)[-1][1:]
    fullbfile = osp.abspath(osp.join(osp.split(a_filename)[0], bfile))
    bfileexists = osp.isfile(fullbfile)
    if not bfileexists:
        raise FileNotFoundError("B-file dynalog file not found; ensure B-file is in same directory as A-file")
    else:
        return fullbfile

# ------------------------
# MLC Log Viewer Example
# ------------------------
if __name__ == '__main__':
    # import cProfile
    # cProfile.run('MachineLog().run_tlog_demo()', sort=1)
    log = MachineLog()
    # log.load_demo_trajectorylog()
    log.run_tlog_demo()
