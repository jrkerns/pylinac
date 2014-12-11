from __future__ import division, absolute_import, print_function

import struct
import os.path as osp
import csv

from future.builtins import next
from future.builtins import range
from future.builtins import object
import numpy as np
import scipy.ndimage.filters as spf
import matplotlib.pyplot as plt

from pylinac.core.common_functions import get_filename
from decorators import type_accept, lazyproperty


"""Named Constants"""
log_types = ('dynalog', 'trajectory log')  # the two current log types
dlg_file_exts = ('.dlg')  # add more if need be
tlg_file_exts = ('.bin')  # ditto
dynalog_leaf_conversion = 1.96614  # MLC physical plane scaling factor to iso (100cm SAD) plane

class Axis(object):

    def __init__(self, actual, expected=None):
        self.actual = actual
        if expected is not None:
            assert len(actual) == len(expected)
            self.expected = expected

    @property
    def difference(self):
        if self.expected is not None:
            return self.actual - self.expected
        else:
            raise AttributeError("Expected positions not passed to Axis")

    def plot_actual(self):
        plt.plot(self.actual)
        plt.show()

    def plot_expected(self):
        plt.plot(self.expected)
        plt.show()

    def plot_difference(self):
        plt.plot(self.difference)
        plt.show()

class GantryAxis(Axis):
    pass

class LeafAxis(Axis):
    def __init__(self, actual, expected):
        super().__init__(actual, expected)

class HeadAxis(Axis):
    pass

class CouchAxis(Axis):
    pass

class BeamAxis(Axis):
    pass

class MLCBank(Axis):
    def __init__(self):
        super().__init__()


class MLC(dict):

    def __init__(self, beam_on_idx=None):
        super().__init__()
        self._cache = {}  # cache for holding calculated parameters
        self.beam_on_idx = beam_on_idx

    @lazyproperty
    def snapshots(self):
        if self.beam_on_idx is not None:
            return self.beam_on_idx
        else:
            return list(range(len(self[1].actual)))

    @lazyproperty
    def num_leaves(self):
        return len(self)

    @lazyproperty
    def num_snapshots(self):
        return len(self.snapshots)

    @lazyproperty
    def moving_leaf_indices(self):
        threshold = 0.003
        indices = ()
        for leaf_num, leafdata in self.items():
            leaf_stdev = np.std(leafdata.actual[self.snapshots])
            if leaf_stdev > threshold:
                indices += (leaf_num,)
        return indices

    @lazyproperty
    def all_leaf_indices(self):
        """Return an array enumerated over the number of leaves."""
        return np.array(list(range(1, len(self) + 1)))

    def get_RMS_per_leaf(self, bank=None, only_moving_leaves=False):
        """Calculate the root-mean-square (RMS) leaf error while the beam was on.

        return: average RMS, matrix of RMS per leaf (usually 1-x-120 numpy array)

        bank: Specify "A" or "B" if the RMS of the bank is desired. None will give both banks.
        only_moving_leaves: boolean specifying if the RMS should be calculated only on the leaves that moved.
            If False, the RMS will usually be lower since non-moving leaves have an RMS of 0 and will drive down
            the average value.
        """

        # get indices of leaves that moved if asked for
        if only_moving_leaves:
            leaves = self.moving_leaf_indices
        else:
            leaves = self.all_leaf_indices

        # get leaves of specific bank if asked for
        if bank is not None:
            leaves = get_bank_index(bank, leaves)

        # remove snapshots were beam was on if asked for
        if self.beam_on_idx is not None and len(self[1].actual) != self.num_snapshots:
            for leaf in self.keys():
                self[leaf].actual = self[leaf].actual[self.snapshots]
                self[leaf].expected = self[leaf].expected[self.snapshots]

        error_squared = self.mlc_error(leaves) ** 2
        RMS = np.sqrt(np.sum(error_squared, 1) / self.num_snapshots)
        return RMS
        # else:  # calculate both banks and attach them to self
        #     # bank A
        #     aidx = get_bank_index('A', idx)
        #     error_squared = self._get_mlc_error(aidx) ** 2
        #     RMS = np.sqrt(np.sum(error_squared, 1) / np.sum(self._beamon_idx))
        #     self.bankA_RMS = RMS
        #     # bank B
        #     bidx = get_bank_index('B', idx)
        #     error_squared = self._get_mlc_error(bidx) ** 2
        #     RMS = np.sqrt(np.sum(error_squared, 1) / np.sum(self._beamon_idx))
        #     self.bankB_RMS = RMS

    def get_overall_RMS_avg(self, bank=None, only_moving_leaves=False):
        """Return the overall average RMS of given leaves.
        """
        rms = self.get_RMS_per_leaf(bank=bank, only_moving_leaves=only_moving_leaves)
        return np.mean(rms)

    def get_overall_RMS_max(self, bank=None, only_moving_leaves=False):
        """Return the overall maximum RMS of given leaves.
        """
        rms = self.get_RMS_per_leaf(bank=bank, only_moving_leaves=only_moving_leaves)
        return np.max(rms)

    def mlc_error(self, leaves):
        """
        Calculate the difference between the planned and expected MLC positions for every leaf and every snapshot
        while the beam was on.
        """
        # preallocate
        mlc_error = np.zeros((len(leaves), self.num_snapshots))
        # construct numpy array for easy array calculation
        for idx, leaf in enumerate(leaves):
            mlc_error[idx, :] = self[leaf].difference
        return mlc_error

    def get_95th_perc_error(self, bank=None, only_moving_leaves=False):
        """Calculate the 95th percentile error of the leaves as asked for by TG-142.

        bank: Specify "A" or "B" if the 95th percentile error of the bank is desired. None will give both banks.
        only_moving_leaves: boolean specifying if the error should be calculated only on the leaves that moved.
            If False, the error will usually be lower since non-moving leaves have an error of 0 and will drive down
            the average value.
        return: 95th percentile error of all leaves, or of the given bank
        """
        # get indices of leaves that moved if asked for
        if only_moving_leaves:
            leaves = self.moving_leaf_indices
        else:
            leaves = self.all_leaf_indices

        # get bank indices if a specific bank is asked for.
        leaves = get_bank_index(bank, leaves)

        # calculate percentile
        error = np.abs(self.mlc_error(leaves))
        return np.percentile(error, 95)

    # def calc_RMS(self, only_moving_leaves=False):
    #     """Calculate the root mean square of the MLC leaves.
    #
    #     :param only_moving_leaves: boolean specifying whether to use only the leaves that moved during treatment
    #     """
    #     # get beam on index
    #     if not hasattr(self, '_beamon_idx'):
    #         self._get_beamon_idx()
    #
    #     # get indices of leaves that moved if asked for
    #     idx = self._get_leaf_idx(only_moving_leaves)
    #
    #     idx = get_bank_index('A', idx)
    #     error_squared = self._get_mlc_error(idx) ** 2
    #     RMS = np.sqrt(np.sum(error_squared, 1) / np.sum(self._beamon_idx))
    #     self.bankA_RMS = RMS
    #
    # def _get_leaf_idx(self, only_moving_leaves):
    #     """
    #     If only_moving_leaves is True, function gets the indices of the leaves that actually moved during the treatment,
    #     otherwise all the leaves are returned.
    #     """
    #     if only_moving_leaves:
    #         if not hasattr(self, '_leaves_that_moved'):
    #             idx = self.get_moving_leaves()
    #         else:
    #             idx = self._leaves_that_moved
    #     else:
    #         idx = np.array(list(range(self.num_mlc_leaves)))
    #     return idx


class Subbeam(object):
    """Holds subbeam information of Tlogs. Only applicable for auto-sequenced beams.
    """
    def __init__(self, log_content, cursor, log_version):
        self._cursor = cursor
        self.control_point = self._decode_binary(log_content, int)
        self.mu_delivered = self._decode_binary(log_content, float)
        self.rad_time = self._decode_binary(log_content, float)
        self.sequence_num = self._decode_binary(log_content, int)
        # In Tlogs version 3.0 and up, beam names are 512 byte unicode strings, but in <3.0 they are 32 byte unicode strings
        if log_version >= 3:
            chars = 512
        else:
            chars = 32
        self.beam_name = self._decode_binary(log_content, str, chars, 32)

    def _decode_binary(self, filecontents, dtype, num_values=1, cursor_shift=0):
        """This method is the main "decoder" for reading in trajectory log binary data into human data types.

        :param filecontents: the complete file having been read with .read().
        :param dtype: the expected data type to return. If int or float, will return numpy array.
        :type dtype: str, int, float
        :param num_values: the expected number of dtype to return; note that this is not the same as the # of bytes.
        :type num_values: int
        :param cursor_shift: the number of bytes to move the cursor forward after decoding. This is used if there is a
            reserved section after the read-in segment.
        :type cursor_shift: int
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

# def yield_axis(snapshot_data, axis_type):
#     column = 0
#     last_column = snapshot_data.shape[0]
#     while column < last_column:
#         yield axis_type(expected=snapshot_data[:, column],
#                    actual=snapshot_data[:, column+1])
#         column += 2


class MachineLog(object):
    """
    A class for reading in machine log files (both dynalogs and trajectory logs) from Varian machines
    and calculating various relevant parameters about them (RMS, 95th percentile error, etc).
    """

    @type_accept(filename=str, lightweight_mode=bool)
    def __init__(self, filename='', lightweight_mode=False):
        """
        :param filename: Path to the log file. For trajectory logs this is a single .bin file.
            For dynalog files this should be the A-file. The B-file will be automatically pulled when A is read in. The B-file must be in
            the same directory as the A-file or an error is thrown. If filename is not passed in on init, it will need to be loaded later.
        :type filename: string
        :param lightweight_mode: Specifies if fluences and gamma maps should be saved as attrs or if results should
            be returned on given calculation. lightweight mode is good for reading in lots of files so as to conserve
            memory, but requires on-the-fly calculations, which is slower then pre-calculating.
            Thus, when reading a handful of files, lightweight mode is really not necessary.
        :type lightweight_mode: bool
        """
        self._filename = filename

        self._cursor = 0

        self._in_lightweight_mode = lightweight_mode

        self._log_type = ''

        #------------------------
        # Generic Log attributes
        #------------------------
        # For dynalogs:   . For Tlogs: usually in 2.x TODO: find dynalog version info
        self.version = ''
        # the number of leaves in a log file. For dynalogs it is actually the # of PAIRS.
        self.num_mlc_leaves = None
        # For dynalogs: 0->Varian scale, 1->IEC scale. For Tlogs: 1-> machine scale, 2-> modified IEC
        self.clinac_scale = None
        # The number of "snapshots" taken during the treatment. Each snapshot captures numerous mechanical parameters
        self.num_snapshots = None
        # The sampling interval (i.e. time between snapshots) in ms. For dynalogs, this is 50ms. Tlogs are usu. ~20ms.
        self.sampling_interval = None
        # The actual locations of the MLCs in mm. Will be a num_leaves-x-num_snapshots numpy matrix.
        # self._mlc_actual = None
        # The expected locations of the MLCs in mm. Will be a num_leaves-x-num_snapshots numpy matrix.
        # self._mlc_expected = None

        #----------------------------
        # Dynalog-specific attributes
        #----------------------------
        # patient name, up to 25 characters
        self.patient_name = ''

        self.plan_filename = ''

        # The mechanical tolerances. Machine will beam-hold if error is greater than tolerance. (mm)
        self.tolerance = None

        #----------------------------
        # Trajectory log-specific attributes
        #----------------------------
        # 2->NDS 120 (Millennium), 3->NDS 120 HD
        self.mlc_model = None

        # A data list for auto-sequenced beams (currently not used)
        # self._subbeams = []

        self.num_axes = None  # # of axes sampled
        self.axis_enum = None
        self.samples_per_axis = None

        self.is_truncated = None  # 1-> truncated, 0 -> not truncated

        # Read file if passed in
        if filename is not '':
            self.read_log()

    def load_demo_dynalog(self):
        """Load the demo dynalog file included with the package."""
        self._filename = osp.join(osp.split(osp.abspath(__file__))[0], 'demo_files', 'AQA.dlg')
        self.read_log()

    def load_demo_trajectorylog(self):
        """Load the log file to the demo trajectory log included with the package."""
        self._filename = osp.join(osp.split(osp.abspath(__file__))[0], 'demo_files', 'Tlog.bin')
        self.read_log()

    def load_logfile_UI(self):
        """Let user load a log file with a UI dialog box. """
        filename = get_filename()
        if filename: # if user didn't hit cancel...
            self._filename = filename
            self.read_log()

    @type_accept(filename=str)
    def load_logfile(self, filename):
        """Load the log file directly by passing the path to the file.

        :param filename: The path to the log file.
        :type filename: str
        """
        if not osp.isfile(filename):
            try:
                raise FileExistsError("File does not exist")
            except:
                raise IOError("File does not exist")
        self._filename = filename
        self.read_log()

    def report_basic_parameters(self):
        """Print the common parameters analyzed when investigating machine logs:

        -log type
        -average RMS
        -95th percentile error
        -number of beam holdoffs
        -gamma pass percentage
        """
        print("MLC log type: " + self.log_type)
        print("Average RMS of all leaves: {:3.3f} mm".format(self.mlc.get_overall_RMS_avg(only_moving_leaves=True)))
        print("Max RMS error of all leaves: {:3.3f} mm".format(self.mlc.get_overall_RMS_max(only_moving_leaves=True)))
        print("95th percentile error: {:3.3f} mm".format(self.mlc.get_95th_perc_error(only_moving_leaves=True)))
        print("Number of beam holdoffs: {:1.0f}".format(self.num_beamholds))
        self.calc_gamma_stats()
        print("Gamma pass %: {:2.2f}".format(self._gamma_pass_pct))
        print("Gamma average: {:2.3f}".format(self._gamma_avg))

    @lazyproperty
    def num_beamholds(self):
        """The number of times the beam was held.
        """
        diffmatrix = np.diff(self.beam_hold.actual)
        num_holds = np.sum(diffmatrix == 1)
        return num_holds

    def calc_fluence(self, resolution=0.1):
        """Calculate the expected and actual fluences.

        :param resolution: The resolution in mm of the fluence calculation in the leaf-moving direction.
        :type resolution: float
        """
        # check for beam-on index
        if not hasattr(self,'_beamon_idx'):
            self._get_beamon_idx()

        # preallocate arrays for expected and actual fluence of number of leaf pairs-x-4000 (40cm = 4000um, etc)
        actual_fluence = np.zeros((self.num_mlc_leaves, 400/resolution), dtype=float)
        expected_fluence = np.zeros((self.num_mlc_leaves, 400/resolution), dtype=float)

        # calculate physical MLC positions based on the log readings. Must also be converted to resolution asked for.
        act_lft_mlc_pos = 200 / resolution - np.round(self._mlc_actual[:self.num_mlc_leaves/2, :] / 10)
        act_rt_mlc_pos = 200 / resolution + np.round(self._mlc_actual[self.num_mlc_leaves/2:, :] / 10)
        exp_lft_mlc_pos = 200 / resolution - np.round(self._mlc_expected[:self.num_mlc_leaves/2, :] / 10)
        exp_rt_mlc_pos = 200 / resolution + np.round(self._mlc_expected[self.num_mlc_leaves/2:, :] / 10)

        # calculate the MU delivered in each snapshot. For Tlogs this is absolute; for dynalogs it's normalized.
        act_MU_fx = np.zeros(len(self._MU_actual))
        act_MU_fx[0] = self._MU_actual[0]
        act_MU_fx[1:] = np.diff(self._MU_actual)
        act_MU_fx = act_MU_fx/self._MU_actual[-1]

        exp_MU_fx = np.zeros(len(self._MU_expected))
        exp_MU_fx[0] = self._MU_expected[0]
        exp_MU_fx[1:] = np.diff(self._MU_expected)
        exp_MU_fx = exp_MU_fx/self._MU_expected[-1]

        # calculate each "line" of fluence (the fluence of an MLC leaf pair, e.g. 1 & 61, 2 & 62, etc),
        # and add each "line" to the total fluence matrix
        plan_line = np.zeros((400/resolution), dtype=float)
        actual_line = np.zeros((400/resolution), dtype=float)
        for leaf in range(self.num_mlc_leaves//2):
            plan_line[:] = 0 # emtpy the line values on each new leaf pair
            actual_line[:] = 0
            for snapshot, __ in enumerate(self._beamon_idx):
                #TODO: incorporate numexpr and multithreading
                # neeval = ne.evaluate('plan_line[exp_lft_mlc_pos[leaf, snapshot]:exp_rt_mlc_pos[leaf,snapshot]] += exp_MU_fx[snapshot]')
                plan_line[exp_lft_mlc_pos[leaf, snapshot]:exp_rt_mlc_pos[leaf,snapshot]] += exp_MU_fx[snapshot]
                actual_line[act_lft_mlc_pos[leaf, snapshot]:act_rt_mlc_pos[leaf, snapshot]] += act_MU_fx[snapshot]
            expected_fluence[leaf, :] = plan_line
            actual_fluence[leaf, :] = actual_line

        if not self._in_lightweight_mode:
            self.fluence_actual = actual_fluence
            self.fluence_expected = expected_fluence
        else:
            return actual_fluence, expected_fluence

    def calc_gamma_map(self, DoseTA=2, DistTA=1, threshold=10, resolution=0.1, actual_fluence=None, expected_fluence=None):
        """
        Calculate the gamma from the actual and expected fluences. The calculation is based on `Bakai et al
        <http://iopscience.iop.org/0031-9155/48/21/006/>`_ eq.6,
        which is a quicker alternative to the standard gamma equation.


        :param DoseTA: dose-to-agreement in %
        :type DoseTA: float, int
        :param DistTA: distance-to-agreement in mm
        :type DistTA: float, int
        :param threshold: the dose threshold percentage of the maximum dose below which is not analyzed for gamma analysis
        :type threshold: float, int
        :param resolution: the resolution in mm of the resulting gamma map
        :type resolution: float
        :returns: a num_mlc_leaves-x-400/resolution numpy matrix
        """

        # if fluences are not passed in, assume they are attrs of self
        if actual_fluence is None or expected_fluence is None:
            try:
                actual_fluence = self.fluence_actual
                expected_fluence = self.fluence_expected
            except:
                # print "Fluence has not been calculated yet; calculating..."
                if self._in_lightweight_mode:
                    actual_fluence, expected_fluence = self.calc_fluence(resolution)
                else:
                    self.calc_fluence(resolution)
                    actual_fluence = self.fluence_actual
                    expected_fluence = self.fluence_expected

        # set dose values below threshold to 0 so gamma doesn't calculate over it
        actual_fluence[actual_fluence < (threshold/100)*np.max(actual_fluence)] = 0
        expected_fluence[expected_fluence < (threshold/100)*np.max(expected_fluence)] = 0

        # preallocation
        gamma_map = np.zeros((self.num_mlc_leaves,400/resolution))

        # image gradient in x-direction using sobel filter
        img_x = spf.sobel(actual_fluence,1)
        # flu = plt.imshow(img_x, aspect='auto')
        # plt.show()
        # img_y, img_xx = np.gradient(actual_fluence)
        # flu = plt.imshow(img_xx, aspect='auto')
        # plt.show()
        # flu = plt.imshow(img_y, aspect='auto')
        # plt.show()

        # equation: measurement - reference / sqrt ( doseTA^2 + distTA^2 * image_gradient^2 )
        for leaf in range(self.num_mlc_leaves):
            gamma_map[leaf] = (actual_fluence[leaf,:] - expected_fluence[leaf,:]) / np.sqrt(
                (DoseTA/100.0 ** 2) + ((DistTA/resolution ** 2) * (img_x[leaf,:] ** 2)))
        # flu = plt.imshow(gamma_map, aspect='auto')
        # plt.show()
        if not self._in_lightweight_mode:
            self._gamma_map = gamma_map
        else:
            return gamma_map

    def calc_gamma_stats(self, DoseTA=2, DistTA=1, threshold=10, resolution=0.1, only_nonzero_gamma=False, gamma_map=None):
        """
        Calculate common gamma parameters: pass percentage, average gamma

        gamma_map: gamma map matrix calculated with _calc_gamma() in lightweight mode
        only_nonzero_gamma: boolean of whether to calculate gamma on all pixels or only on the non-zero ones
        """

        if gamma_map is None:
            try:
                gamma_map = self._gamma_map
            except:
                if self._lightweight_mode:
                    gamma_map = self.calc_gamma_map(DoseTA=DoseTA, DistTA=DistTA, threshold=threshold,
                                                 resolution=resolution)
                else:
                    self.calc_gamma_map(DoseTA=DoseTA, DistTA=DistTA, threshold=threshold,
                                         resolution=resolution)
                    gamma_map = self._gamma_map

        if only_nonzero_gamma:
            pass_pct = (np.sum(gamma_map > 0) - np.sum(gamma_map > 1))/float(np.sum(gamma_map > 0))
            avg_gamma = np.mean(gamma_map[gamma_map > 0])
        else:
            pass_pct = np.sum(gamma_map < 1)/float(gamma_map.size)
            avg_gamma = np.mean(gamma_map)
        # turn pass_pct into percentage
        pass_pct *= 100

        if self._lightweight_mode:
            return pass_pct, avg_gamma
        else:
            self._gamma_pass_pct = pass_pct
            self._gamma_avg = avg_gamma

    # def get_moving_leaves(self, threshold=0.003):
    #     """
    #     Determine the leaves that actually moved during treatment. Useful for correcting stats when only a few leaves
    #     moved during the delivery.
    #
    #     :param threshold: the value that the standard deviation of the leaf position must be greater than in order to
    #         be considered a "moving" leaf. Rounding errors can result in nonzero standard deviation for standing leaves,
    #         thus, a small value is necessary.
    #     """
    #     stdevs = np.std(self._mlc_actual, 1)
    #     self._leaves_that_moved = np.squeeze(np.asarray(np.nonzero(stdevs > threshold), dtype=int))
    #     return self._leaves_that_moved

    def _get_beamon_idx(self, snapshot_data, beam_hold_column):
        """Return the indices of the snapshots that the beam was actually on.

        For dynalogs this removes snapshots where the Beam On flag was 0 and where the Beam Hold was 0.
        For trajectory logs this removes the snapshots were the Beam Hold was 0 (there is no Beam On flag).
        """
        if self.log_type == log_types[1]:  # trajectory log
            beamon_idx = np.squeeze(np.asarray(np.nonzero(snapshot_data[:,beam_hold_column] == 0), dtype=bool))
        elif self.log_type == log_types[0]:
            holdidx = self.beam_hold.actual == 0
            beamonidx = self.beam_hold.actual == 1
            beamon_idx = holdidx & beamonidx
        return beamon_idx

    def _get_logtype(self):
        """
        Determine the type of MLC log by first looking at its file extension. Failing that, open it and sample
        the first few bytes. If the sample matches what is found in standard dynalog or tlog files it will be
        assigned that log type.
        """
        __, ext = osp.splitext(self._filename)
        if ext == tlg_file_exts:
            self.log_type = log_types[1] # trajectory log
        elif ext == dlg_file_exts:
            self.log_type = log_types[0] # dynalog
        else:
            with open(self._filename, 'rb') as unknown_file:
                header_sample = unknown_file.read(5).decode()
                if 'B' in header_sample or 'A' in header_sample:
                    self.log_type = log_types[0]  # dynalog
                elif 'V' in header_sample:
                    self.log_type = log_types[1]  # trajectory log
                else:
                    raise ValueError("Log type unknown")

        return self.log_type

    def read_log(self):
        """Read in log based on what type of log it is: Trajectory or Dynalog."""
        if self._filename == '':
            raise AttributeError('Log file has not been specified. Use load_logfile_UI or load_logfile')

        # determine log type
        self._get_logtype()

        # read log as appropriate to type
        if self.log_type == log_types[1]: # if traj log
            self.read_Tlog()
        elif self.log_type == log_types[0]: # if dynalog
            self.read_Dlog()

    def _scale_dlog_mlc_pos(self):
        # convert MLC leaf plane positions to isoplane positions by applying scaling factor
        # Units are initially in 100ths of mm, but are converted to mm.
        for leaf in range(1, self.num_mlc_leaves+1):
            self.mlc[leaf].actual *= dynalog_leaf_conversion / 100
            self.mlc[leaf].expected *= dynalog_leaf_conversion / 100
        # self._mlc_expected = np.asarray(self._mlc_expected * self._dlg_leaf_conversion, dtype=float) / 100
        # self._mlc_actual = np.asarray(self._mlc_actual * self._dlg_leaf_conversion, dtype=float) / 100

    def read_Dlog(self):
        """
        Read in Dynalog files from .dlg files (which are renamed CSV files)
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
            self.num_mlc_leaves = int(next(dlgdata)[0]) * 2 # the # of leaves in a dynalog is actually the # of *PAIRS*, hence the *2.
            self.clinac_scale = next(dlgdata) # 0->Varian scale, 1->IEC scale

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
        self.carriageB = Axis(matrix[:, 13])
        # MLC positions are in hundredths of mm in the physical leaf plane. Positive is retracted, negative is extented.
        # Positions will be scaled to isocenter plane after bank B positions are added to matrix.
        self.mlc = {}
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

    def read_Tlog(self, exclude_beam_off=False):
        """Read in Trajectory log from binary file according to TB 1.5/2.0 (i.e. Tlog 2.0/3.0) log file specifications.

        :param exclude_beam_off: Flag specifiying whether to remove the snapshot data where the beam was off.
            If True (default), beam-off data will be removed.
            If False, all data will be included. Note that this will affect MLC and fluence calculations as error, fluence, etc
            are calculated from the data extracted here.
        :type exclude_beam_off: boolean
        """

        # read in trajectory log binary data
        fcontent = open(self._filename, 'rb').read()

        # Unpack the content according to respective section and data type (see log specification file).
        self.header = self._decode_binary(fcontent,str,16)  # for version 1.5 will be "VOSTL"
        self.version = float(self._decode_binary(fcontent, str, 16))  # in the format of 2.x or 3.x
        self.header_size = self._decode_binary(fcontent, int)  # fixed at 1024 in 1.5 specs
        self.sampling_interval = self._decode_binary(fcontent, int)
        self.num_axes = self._decode_binary(fcontent, int)
        self.axis_enum = self._decode_binary(fcontent, int, self.num_axes)
        self.samples_per_axis = self._decode_binary(fcontent, int, self.num_axes)
        self.num_mlc_leaves = self.samples_per_axis[-1]-2  # subtract 2 (each carriage counts as an "axis" and must be removed)
        # self._cursor == self.num_axes * 4 # there is a reserved section after samples per axis. this moves it past it.
        self.clinac_scale = self._decode_binary(fcontent, int)
        self.num_subbeams = self._decode_binary(fcontent, int)
        self.is_truncated = self._decode_binary(fcontent, int)
        self.num_snapshots = self._decode_binary(fcontent, int)
        # the section after MLC model is reserved. Cursor is moved to the end of this reserved section.
        self.mlc_model = self._decode_binary(fcontent, int, cursor_shift=1024 - (64 + self.num_axes * 8))

        # read in subbeam data. These are for auto-sequenced beams. If not autosequenced, separate logs are created.
        # Currently there is no good way of dealing with this data, but fortunately autosequencing is rare at this time.
        if self.num_subbeams:
            self.subbeams = []
            for beam in range(self.num_subbeams):
                self.subbeams.append(Subbeam(fcontent, self._cursor, self.version))
                # update cursor position to end of subbeam just analyzed.
                self._cursor = self.subbeams[beam]._cursor

        #----------------------------------------------------------------------
        # assignment of snapshot data (actual & expected of MLC, Jaw, Coll, etc)
        #----------------------------------------------------------------------



        # step size in bytes
        step_size = sum(self.samples_per_axis) * 2

        # read in all snapshot data at once, then assign
        snapshot_data = self._decode_binary(fcontent, float, step_size*self.num_snapshots)

        # reshape snapshot data to be a x-by-num_snapshots matrix
        snapshot_data = snapshot_data.reshape(self.num_snapshots,-1)



        # collimator
        self.collimator = self._get_axis(snapshot_data, 0, HeadAxis)

        # gantry
        self.gantry = self._get_axis(snapshot_data, 2, GantryAxis)

        # jaws
        self.jaw_y1 = self._get_axis(snapshot_data, 4, HeadAxis)
        self.jaw_y2 = self._get_axis(snapshot_data, 6, HeadAxis)
        self.jaw_x1 = self._get_axis(snapshot_data, 8, HeadAxis)
        self.jaw_x2 = self._get_axis(snapshot_data, 10, HeadAxis)

        # couch
        self.couch_vrt = self._get_axis(snapshot_data, 12, CouchAxis)
        self.couch_lng = self._get_axis(snapshot_data, 14, CouchAxis)
        self.couch_lat = self._get_axis(snapshot_data, 16, CouchAxis)
        self.couch_rtn = self._get_axis(snapshot_data, 18, CouchAxis)

        # MU
        self.mu = self._get_axis(snapshot_data, 20, BeamAxis)

        # beam hold state
        self.beam_hold = self._get_axis(snapshot_data, 22, BeamAxis)

        # control point
        self.control_point = self._get_axis(snapshot_data, 24, BeamAxis)

        # carriages
        self.carraigeA = self._get_axis(snapshot_data, 26, HeadAxis)
        self.carriageB = self._get_axis(snapshot_data, 28, HeadAxis)

        # preallocation
        # mlc_expected = np.zeros((self.num_mlc_leaves, self.num_snapshots)) # usually a 120-x-num_snapshots matrix
        # mlc_actual = np.zeros((self.num_mlc_leaves, self.num_snapshots))
        # assign MLC data for all leaves. Usually 1-60 is bank A, 61-120 is bank B
        # Units are initially in cm and are converted to mm by the *10.

        # for snapshot in self.num_snapshots:
            # mlc_expected[:, snapshot] = snapshot_data[30 + 2*snapshot::step_size].transpose() * 10
            # mlc_actual[:, snapshot] = snapshot_data[31 + 2*snapshot::step_size].transpose() * 10
            # mlc_expected[:, snapshot] = snapshot_data[30 + (idx * step_size):269 + (idx * step_size):2].transpose() * 10
            # mlc_actual[:, snapshot] = snapshot_data[31 + (idx * step_size):270 + (idx * step_size):2].transpose() * 10
        # remove snapshots where the beam wasn't on if flag passed
        if exclude_beam_off:
            beam_on_idx = self._get_beamon_idx(snapshot_data, 22)
        else:
            beam_on_idx = list(range(self.num_snapshots))
            # snapshot_data = snapshot_data[beam_on_idx, :]
        self.mlc = MLC(beam_on_idx)
        for leaf in range(1, self.num_mlc_leaves+1):
            self.mlc[leaf] = self._get_axis(snapshot_data, 30+2*(leaf-1), LeafAxis)
            # self.mlc[leaf] = (Axis(expected=snapshot_data[30+2*leaf::step_size] * 10,
            #                            actual=snapshot_data[31+2*leaf::step_size] * 10))


    def _get_axis(self, snapshot_data, column, axis_type):
        """Return column of data from snapshot data of the axis type passed."""
        return axis_type(expected=snapshot_data[:,column],
                         actual=snapshot_data[:,column+1])

    def _decode_binary(self, filecontents, dtype, num_values=1, cursor_shift=0):
        """This method is the main "decoder" for reading in trajectory log binary data into human data types.

        :param filecontents: the complete file having been read with .read().
        :param dtype: the expected data type to return. If int or float, will return numpy array.
        :type dtype: str, int, float
        :param num_values: the expected number of dtype to return; note that this is not the same as the # of bytes.
        :type num_values: int
        :param cursor_shift: the number of bytes to move the cursor forward after decoding. This is used if there is a
            reserved section after the read-in segment.
        :type cursor_shift: int
        """
        fc = filecontents

        if dtype == str: # if string
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

        self._cursor += cursor_shift # shift cursor if need be (e.g. if a reserved section follows)
        return output


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

def get_bank_index(bank, idx):
    """
    If a bank is specified, the function returns the leaf indices only from that bank, otherwise it is untouched
    """
    if bank is not None:
        if 'a' in bank.lower():
            idx = idx[idx < 60]
        elif 'b' in bank.lower():
            idx = idx[idx >= 60]
    return idx


def decode_binary(filecontents, dtype, cursor, num_values=1, cursor_shift=0):
    """This method is the main "decoder" for reading in trajectory log binary data into human data types.

    :param filecontents: the complete file having been read with .read().
    :param dtype: the expected data type to return. If int or float, will return numpy array.
    :type dtype: str, int, float
    :param num_values: the expected number of dtype to return; note that this is not the same as the # of bytes.
    :type num_values: int
    :param cursor_shift: the number of bytes to move the cursor forward after decoding. This is used if there is a
        reserved section after the read-in segment.
    :type cursor_shift: int
    """
    fc = filecontents

    if dtype == str:  # if string
        output = fc[cursor:cursor + num_values]
        if type(fc) is not str:  # in py3 fc will be bytes
            output = output.decode()
        # Strip the padding ("\x00")
        output = output.strip('\x00')
        cursor += num_values
    elif dtype == int:
        ssize = struct.calcsize('i') * num_values
        output = np.asarray(struct.unpack('i' * num_values, fc[cursor:cursor + ssize]))
        cursor += ssize
    elif dtype == float:
        ssize = struct.calcsize('f') * num_values
        output = np.asarray(struct.unpack('f' * num_values, fc[cursor:cursor + ssize]))
        cursor += ssize
    else:
        raise TypeError("decode_binary datatype was not valid")

    cursor += cursor_shift  # shift cursor if need be (e.g. if a reserved section follows)
    return output, cursor




# ------------------------
# MLC Log Viewer Example
# ------------------------
if __name__ == '__main__':
    mlc = MachineLog()
    # mlc.load_demo_trajectorylog()
    mlc.load_logfile_UI()
    # mlc.load_demo_dynalog()
    # pass
    # mlc.load_demo_file_trajectorylog()  # uncomment to use trajectory log file
    mlc.report_basic_parameters()

