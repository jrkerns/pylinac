from __future__ import division, absolute_import, print_function

import struct
import os.path as osp
import csv

from future.builtins import next
from future.builtins import range
from future.builtins import object
import numpy as np
import scipy.ndimage.filters as spf

from pylinac.common.common_functions import get_filename, open_PDF_file, go_up_dirlevel



# Set up some constants and put them at the top to make later adjustment easy
log_types = ('dynalog', 'trajectory log')  # the two current log types
dlg_file_exts = ('.dlg')  # add more if need be
tlg_file_exts = ('.bin')  # ditto
dynalog_leaf_conversion = (1.96614)  # MLC physical plane scaling factor to iso (100cm SAD) plane

class MachineLog(object):
    """
    A class for reading in machine log files (both dynalogs and trajectory logs) from Varian machines
    and calculating various relevant parameters about them (RMS, 95th percentile error, etc).
    """

    def __init__(self, filename='', lightweight_mode=False):
        """
        A class for reading in machine log files (both dynalogs and trajectory logs) from Varian machines
        and calculating various relevant parameters about them (RMS, 95th percentile error, etc)

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
        # File name/path to MLC log
        self._filename = filename

        # lightweight mode can be used when reading in lots of logs. Rather than setting large matrixes as class attrs,
        # it returns the values so that they can be later garbage-collected. Use is not recommended except when reading
        # many (~25+) files, although if you have the memory you can use lightweight mode anyway.
        self._lightweight_mode = lightweight_mode

        # The log type will be 'dynalog' or 'trajectory log', or as specified in the log_types set above.
        self.log_type = ''

        # An int to keep track of cursor place in file while reading Tlogs
        self._cursor = 0

        # A value to convert dynalog MLC snapshot data from the physical plane (at ~50cm) to the isoplace (100cm).
        self._dlg_leaf_conversion = dynalog_leaf_conversion

        #------------------------
        # Generic Log attributes
        #------------------------
        # For dynalogs:   . For Tlogs: usually in 2.x TODO: find dynalog version info
        self.version = ''
        # the number of leaves in a log file. For dynalogs it is actually the # of PAIRS.
        self.num_mlc_leaves = None
        # For dynalogs: 0->Varian scale, 1->IEC scale. For Tlogs: 1-> machine scale, 2-> modified IEC
        self._clinac_scale = None
        # The number of "snapshots" taken during the treatment. Each snapshot captures numerous mechanical parameters
        self.num_snapshots = None
        # The sampling interval (i.e. time between snapshots) in ms. For dynalogs, this is 50ms. Tlogs are usu. ~20ms.
        self.sampling_interval = None
        # The actual locations of the MLCs in mm. Will be a num_leaves-x-num_snapshots numpy matrix.
        self._mlc_actual = None
        # The expected locations of the MLCs in mm. Will be a num_leaves-x-num_snapshots numpy matrix.
        self._mlc_expected = None

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
        self._subbeams = []

        self._num_axes = None  # # of axes sampled
        self._axis_enum = None
        self._samples_per_axis = None

        self._is_truncated = None  # 1-> truncated, 0 -> not truncated


        # Read file if passed in
        if filename is not '':
            self.read_log()

    def load_demo_file_dynalog(self):
        """
        Set the log file to the demo dynalog file included with the package.
        """
        go_up_dirlevel()
        self._filename = osp.join(osp.split(osp.abspath(__file__))[0], 'demo files', 'AQA.dlg')
        self.read_log()

    def load_demo_file_trajectorylog(self):
        """
        Set the log file to the demo trajectory log included with the package.
        """
        self._filename = osp.join(osp.split(osp.abspath(__file__))[0], 'demo files', 'Tlog.bin')
        self.read_log()

    def load_logfile_UI(self):
        """
        Let user load a log file with a UI dialog box.
        """
        filename = get_filename()
        if filename is not None: # if user didn't hit cancel...
            self._filename = filename
            self.read_log()

    def load_logfile(self, filename):
        """
        Load the log non-interactively by passing path.

        :param filename: The path to the log file.
        :type filename: str
        """
        assert type(filename) is str, "Filename must be a valid string"
        #TODO: assert that file that string points to exists
        self._filename = filename
        self.read_log()

    def report_basic_parameters(self):
        """
        Print the common parameters analyzed when investigating machine logs:

        -log type
        -average RMS
        -95th percentile error
        -number of beam holdoffs
        -gamma pass percentage
        """
        print("MLC log type: " + self.log_type)
        print("Average RMS of all leaves: {:3.3f} mm".format(self.get_RMSavg(only_moving_leaves=False)))
        print("95th percentile error: {:3.3f} mm".format(self.get_95th_perc_error(only_moving_leaves=False)))
        print("Number of beam holdoffs: {:1.0f}".format(self.get_num_beamholds()))
        self.calc_gamma_stats()
        print("Gamma pass %: {:2.2f}".format(self._gamma_pass_pct))
        print("Gamma average: {:2.3f}".format(self._gamma_avg))

    def _get_leaf_idx(self, only_moving_leaves):
        """
        If only_moving_leaves is True, function gets the indices of the leaves that actually moved during the treatment,
        otherwise all the leaves are returned.
        """
        if only_moving_leaves:
            if not hasattr(self, '_leaves_that_moved'):
                idx = self.get_moving_leaves()
            else:
                idx = self._leaves_that_moved
        else:
            idx = np.array(list(range(self.num_mlc_leaves)))
        return idx

    def calc_RMS(self, only_moving_leaves=False):
        """Calculate the root mean square of the MLC leaves.

        :param only_moving_leaves: boolean specifying whether to use only the leaves that moved during treatment
        """
        # get beam on index
        if not hasattr(self, '_beamon_idx'):
            self._get_beamon_idx()

        # get indices of leaves that moved if asked for
        idx = self._get_leaf_idx(only_moving_leaves)

        idx = get_bank_index('A', idx)
        error_squared = self._get_mlc_error(idx) ** 2
        RMS = np.sqrt(np.sum(error_squared, 1) / np.sum(self._beamon_idx))
        self.bankA_RMS = RMS

    def get_RMS(self, bank=None, only_moving_leaves=True, return_it=True):
        """
        Calculate the root-mean-square (RMS) leaf error while the beam was on.

        return: average RMS, matrix of RMS per leaf (usually 1-x-120 numpy array)

        bank: Specify "A" or "B" if the RMS of the bank is desired. None will give both banks.
        only_moving_leaves: boolean specifying if the RMS should be calculated only on the leaves that moved.
            If False, the RMS will usually be lower since non-moving leaves have an RMS of 0 and will drive down
            the average value.
        """
        # get beam on index
        if not hasattr(self, '_beamon_idx'):
            self._get_beamon_idx()

        # get indices of leaves that moved if asked for
        idx = self._get_leaf_idx(only_moving_leaves)

        # get bank indices if a specific bank is asked for.
        if return_it:
            idx = get_bank_index(bank, idx)
            error_squared = self._get_mlc_error(idx) ** 2
            RMS = np.sqrt(np.sum(error_squared, 1) / np.sum(self._beamon_idx))
            return RMS
        else:  # calculate both banks and attach them to self
            # bank A
            aidx = get_bank_index('A', idx)
            error_squared = self._get_mlc_error(aidx) ** 2
            RMS = np.sqrt(np.sum(error_squared, 1) / np.sum(self._beamon_idx))
            self.bankA_RMS = RMS
            # bank B
            bidx = get_bank_index('B', idx)
            error_squared = self._get_mlc_error(bidx) ** 2
            RMS = np.sqrt(np.sum(error_squared, 1) / np.sum(self._beamon_idx))
            self.bankB_RMS = RMS

    def get_RMSavg(self, bank=None, only_moving_leaves=True):
        """
        Return the average RMS of given leaves.
        """
        rms = self.get_RMS(bank=bank, only_moving_leaves=only_moving_leaves)
        return np.mean(rms)

    def get_RMSmax(self, bank=None, only_moving_leaves=True):
        """
        Return the maximum RMS of given leaves.
        """
        rms = self.get_RMS(bank=bank, only_moving_leaves=only_moving_leaves)
        return np.max(rms)

    def get_num_beamholds(self):
        """
        Calculate the number of times the beam was held.
        """
        diffmatrix = np.diff(self._beamhold_actual)
        num_holds = np.sum(diffmatrix == 1)
        return num_holds

    def _get_mlc_error(self, idx):
        """
        Calculate the difference between the planned and expected MLC positions for every leaf and every snapshot
        while the beam was on.
        """
        # broadcasting limitations make indexing more syntactically convoluted than it should be
        return (self._mlc_actual[idx, :][:, self._beamon_idx] - self._mlc_expected[idx, :][:, self._beamon_idx])

    def get_95th_perc_error(self, bank=None, only_moving_leaves=True):
        """
        Calculate the 95th percentile error of the leaves as asked for by TG-142.

        return: 95th percentile error of all leaves, or of the given bank
        bank: Specify "A" or "B" if the 95th percentile error of the bank is desired. None will give both banks.
        only_moving_leaves: boolean specifying if the error should be calculated only on the leaves that moved.
            If False, the error will usually be lower since non-moving leaves have an error of 0 and will drive down
            the average value.
        """
        # get indices of leaves that moved if asked for
        idx = self._get_leaf_idx(only_moving_leaves)

        # get bank indices if a specific bank is asked for.
        idx = get_bank_index(bank, idx)

        # calculate percentile
        error = self._get_mlc_error(idx)
        return np.percentile(error,95)

    def calc_fluence(self, resolution=0.1):
        """
        Calculate the expected and actual fluence.

        :type resolution: object
        resolution: the resolution in mm of the fluence calculation in the leaf-moving direction.
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



        if not self._lightweight_mode:
            self.fluence_actual = actual_fluence
            self.fluence_expected = expected_fluence
        else:
            return actual_fluence, expected_fluence

    def calc_gamma_map(self, DoseTA=2, DistTA=1, threshold=10, resolution=0.1, actual_fluence=None, expected_fluence=None):
        """
        Calculate the gamma of the actual and expected fluences. The calculation is based on Bakai et al eq.6,
        which is a quicker alternative to the standard gamma equation.

        returns a num_mlc_leaves-x-400/resolution numpy matrix

        DoseTA: dose-to-agreement in %
        DistTA: distance-to-agreement in mm
        threshold: the dose threshold percentage of the maximum dose below which is not analyzed for gamma analysis
        resolution: the resolution in mm of the resulting gamma map
        """

        # if fluences are not passed in, assume they are attrs of self
        if actual_fluence is None or expected_fluence is None:
            try:
                actual_fluence = self.fluence_actual
                expected_fluence = self.fluence_expected
            except:
                # print "Fluence has not been calculated yet; calculating..."
                if self._lightweight_mode:
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
        if not self._lightweight_mode:
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

    def get_moving_leaves(self, threshold=0.003):
        """
        Determine the leaves that actually moved during treatment. Useful for correcting stats when only a few leaves
        moved during the delivery.

        :param threshold: the value that the standard deviation of the leaf position must be greater than in order to
            be considered a "moving" leaf. Rounding errors can result in nonzero standard deviation for standing leaves,
            thus, a small value is necessary.
        """
        stdevs = np.std(self._mlc_actual, 1)
        self._leaves_that_moved = np.squeeze(np.asarray(np.nonzero(stdevs > threshold), dtype=int))
        return self._leaves_that_moved

    def _get_beamon_idx(self):
        """
        Determine the snapshots that the beam was actually on.

        For dynalogs this removes snapshots where the Beam On flag was 0 and where the Beam Hold was 0.
        For trajectory logs this removes the snapshots were the Beam Hold was 0 (there is no Beam On flag).
        """
        if self.log_type == log_types[1]:  # trajectory log
            self._beamon_idx = np.squeeze(np.asarray(np.nonzero(self._beamhold_actual == 0), dtype=bool))
        elif self.log_type == log_types[0]:
            holdidx = self._beamhold_actual == 0
            beamonidx = self._beamon_actual == 1
            self._beamon_idx = holdidx & beamonidx
        return self._beamon_idx

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
                    raise TypeError("Log type unknown")

        return self.log_type

    def read_log(self):
        """
        Read in log based on what type of log it is: Trajectory or Dynalog.
        """
        assert hasattr(self,'_filename'), 'Log file has not been specified. Use load_logfile_UI or load_logfile'

        # determine log type
        self._get_logtype()

        # read log as appropriate to type
        if self.log_type == log_types[1]: # if traj log
            self.read_Tlog()
        elif self.log_type == log_types[0]: # if dynalog
            self.read_Dlog()

    def read_Dlog(self):
        """
        Read in Dynalog files from .dlg files (which are renamed CSV files)
        Formatting follows from the Dynalog File Viewer Reference Guide.
        """
        # if file is B-file, skip over file
        if self._filename[0] == 'B':
            return

        # Check that the B-file is in the same folder before getting too far
        bfile = check_B_file_exists(self._filename)

        # create iterator object to read in lines
        with open(self._filename) as csvf:
            dlgdata = csv.reader(csvf, delimiter=',')
            self.version = next(dlgdata)
            self.patient_name = next(dlgdata) # patient name, up to 25 characters
            self.plan_filename = next(dlgdata)
            self.tolerance = int(next(dlgdata)[0])
            self.num_mlc_leaves = int(next(dlgdata)[0]) * 2 # the # of leaves in a dynalog is actually the # of *PAIRS*
            self._clinac_scale = next(dlgdata) # 0->Varian scale, 1->IEC scale

            # From here on out, each line is a "snapshot".
            # This reads in all snapshots, then assigns them. The reason this cannot be iterated over is because the
            # snapshots are column-order, not row order, so all rows must be read before any data can be assigned.
            matrix = [line for line in dlgdata]

        # convert to numpy array
        matrix = np.array(matrix, dtype=int)

        self.num_snapshots = np.size(matrix, 0)

        # preallocation
        self._mlc_actual = np.zeros((self.num_mlc_leaves, self.num_snapshots), dtype=float)
        self._mlc_expected = np.zeros((self.num_mlc_leaves, self.num_snapshots), dtype=float)

        # assignment of values
        self._MU_actual = matrix[:,0]
        # There is no "expected" MU in dynalogs, but for fluence calc purposes, it is set to that of the actual
        self._MU_expected = self._MU_actual
        self._DVA_segment = matrix[:,1]
        self._beamhold_actual = matrix[:,2]
        self._beamon_actual = matrix[:,3]
        self._prior_dose_idx = matrix[:, 4]
        self._next_dose_idx = matrix[:, 5]
        self._gantry_actual = matrix[:, 6]
        self._coll_actual = matrix[:, 7]
        self._y1_actual = matrix[:, 8]
        self._y2_actual = matrix[:, 9]
        self._x1_actual = matrix[:, 10]
        self._x2_actual = matrix[:, 11]
        self._carriageA_actual = matrix[:, 12]
        self._carriageB_actual = matrix[:, 13]
        # MLC positions are in hundredths of mm in the physical leaf plane. Positive is retracted, negative is extented.
        # Positions will be scaled to isocenter plane after bank B positions are added to matrix.
        self._mlc_expected[:60,:] = matrix[:, 14::4].transpose()
        self._mlc_actual[:60,:] = matrix[:, 15::4].transpose()

        # read in "B"-file to get bank B MLC positions. The file must be in the same folder as the "A"-file.
        # The header info is repeated but we already have that.
        with open(bfile) as csvf:
            dlgdata = csv.reader(csvf, delimiter=',')
            matrix = [line for line in dlgdata if int(dlgdata.line_num) >= 7]
        matrix = np.array(matrix, dtype=int)

        # Add bank B MLC positions to mlc snapshot arrays
        self._mlc_expected[60:, :] = matrix[:, 14::4].transpose()
        self._mlc_actual[60:, :] = matrix[:, 15::4].transpose()

        # convert MLC leaf plane positions to isoplane positions by applying scaling factor
        # Units are initially in 100ths of mm, but are converted to mm.
        self._mlc_expected = np.asarray(self._mlc_expected * self._dlg_leaf_conversion, dtype=float) / 100
        self._mlc_actual = np.asarray(self._mlc_actual * self._dlg_leaf_conversion, dtype=float) / 100

    def read_Tlog(self):
        """
        Read in Trajectory log from binary file according to TB 1.5/2.0 log file specifications.
        See folder: log_viewer/Log file specifications for more info.
        """

        # read in trajectory log binary data
        fcontent = open(self._filename, 'rb').read()

        # Unpack the content according to respective section and data type (see log specification file).
        header = self._decode_binary(fcontent,str,16)  # for version 1.5 will be "VOSTL"
        self.version = float(self._decode_binary(fcontent, str, 16))  # in the format of 2.x or 3.x
        header_size = self._decode_binary(fcontent, int)  # fixed at 1024 in 1.5 specs
        self.sampling_interval = self._decode_binary(fcontent, int)
        self._num_axes = self._decode_binary(fcontent, int)
        self._axis_enum = self._decode_binary(fcontent, int, self._num_axes)
        self._samples_per_axis = self._decode_binary(fcontent,int,self._num_axes)
        self.num_mlc_leaves = self._samples_per_axis[-1]-2
        self._cursor == self._num_axes * 4 # there is a reserved section after samples per axis. this moves it past it.
        self._clinac_scale = self._decode_binary(fcontent, int)
        self.num_subbeams = self._decode_binary(fcontent, int)
        self._is_truncated = self._decode_binary(fcontent, int)
        self.num_snapshots = self._decode_binary(fcontent, int)
        self.mlc_model = self._decode_binary(fcontent, int, cursor_shift=1024 - (64 + self._num_axes * 8))
        # the next section is reserved. cursor is moved to the end of this reserved section

        # read in subbeam data. These are for auto-sequenced beams. If not autosequenced, separate logs are created.
        # Currently there is no good way of dealing with this data, but fortunately autosequencing is rare at this time.
        for idx in range(self.num_subbeams):
            cont_point = self._decode_binary(fcontent, int)
            mu = self._decode_binary(fcontent, float)
            exp_time = self._decode_binary(fcontent, float)
            seq_num = self._decode_binary(fcontent, float)
            # In Tlogs version 3.0 and up, beam names are 512 byte unicode strings, but in <3.0 they are 32 byte strings
            if self.version >= 3:
                beam_name = self._decode_binary(fcontent, str, 512, 32)
            else:
                beam_name = self._decode_binary(fcontent, str, 32, 32)
            self._subbeams.append({'control_point': cont_point, 'MU':mu, 'exp_time': exp_time,
                                   'sequence_number': seq_num, 'beam_name': beam_name})

        #----------------------------------------------------------------------
        # assignment of snapshot data (actual & expected of MLC, Jaw, Coll, etc)
        #----------------------------------------------------------------------

        # preallocation
        self._mlc_expected = np.zeros((self._samples_per_axis[-1] - 2, self.num_snapshots)) # usually a 120-x-num_snapshots matrix
        self._mlc_actual = np.zeros((self._samples_per_axis[-1] - 2, self.num_snapshots))

        # step size in bytes
        step_size = np.sum(self._samples_per_axis) * 2

        # read in all snapshot data at once, then assign
        snapshot_data = self._decode_binary(fcontent, float, step_size*self.num_snapshots)

        # collimator
        self._coll_expected = snapshot_data[0::step_size]
        self._coll_actual = snapshot_data[1::step_size]

        # gantry
        self._gantry_expected = snapshot_data[2::step_size]
        self._gantry_actual = snapshot_data[3::step_size]

        # jaws
        self._y1_expected = snapshot_data[4::step_size]
        self._y1_actual = snapshot_data[5::step_size]
        self._y2_expected = snapshot_data[6::step_size]
        self._y2_actual = snapshot_data[7::step_size]
        self._x1_expected = snapshot_data[8::step_size]
        self._x1_actual = snapshot_data[9::step_size]
        self._x2_expected = snapshot_data[10::step_size]
        self._x2_actual = snapshot_data[11::step_size]

        # couch
        self._couch_vrt_expected = snapshot_data[12::step_size]
        self._couch_vrt_actual = snapshot_data[13::step_size]
        self._couch_lng_expected = snapshot_data[14::step_size]
        self._couch_lng_actual = snapshot_data[15::step_size]
        self._couch_lat_expected = snapshot_data[16::step_size]
        self._couch_lat_actual = snapshot_data[17::step_size]
        self._couch_rtn_expected = snapshot_data[18::step_size]
        self._couch_rtn_actual = snapshot_data[19::step_size]

        # MU
        self._MU_expected = snapshot_data[20::step_size]
        self._MU_actual = snapshot_data[21::step_size]

        # beam hold state
        self._beamhold_expected = snapshot_data[22::step_size]
        self._beamhold_actual = snapshot_data[23::step_size]

        # control point
        self._controlpoint_expected = snapshot_data[24::step_size]
        self._controlpoint_actual = snapshot_data[25::step_size]

        # carriages
        self._carraigeA_expected = snapshot_data[26::step_size]
        self._carriageA_actual = snapshot_data[27::step_size]
        self._carriageB_expected = snapshot_data[28::step_size]
        self._carriageB_actual = snapshot_data[29::step_size]

        # assign MLC data for all leaves. Usually 1-60 is bank A, 61-120 is bank B
        # Units are initially in cm and are converted to mm.
        for idx in range(self.num_snapshots):
            self._mlc_expected[:, idx] = snapshot_data[30 + (idx * step_size):269 + (idx * step_size):2].transpose() * 10
            self._mlc_actual[:, idx] = snapshot_data[31 + (idx * step_size):270 + (idx * step_size):2].transpose() * 10


    def _decode_binary(self, filecontents, dtype, num_values=1, cursor_shift=0):
        """
        This method is the main "decoder" for reading in trajectory log binary data into another data type.

        filecontents: the complete file having been read with .read()
        dtype: the expected data type to return (usually str, int, or float). If int or float, will return numpy array
        num_values: the expected number of dtype to return; note that this is not the same as the # of bytes
        cursor_shift: the number of bytes to move the cursor forward after decoding. This is used if there is a
            reserved section after the read-in segment
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

    def open_dynalog_reference_PDF(self):
        """
         Open the Dynalog file reference PDF
        """
        open_PDF_file('Dynalog File Viewer & Reference Guide (2011).pdf')

    def open_trajectory_reference_PDF(self):
        """
         Open the Tlog (1.5) file reference PDF
        """
        open_PDF_file('TrueBeam Trajectory Log File Specification For TrueBeam 1.5 And Higher.pdf')


def check_B_file_exists(a_filename):
    """
    Checks that the dynalog "B-file" for a given "A-file" exists within the same directory.

    Returns the absolute locations of the B-file.
    """
    bfile = "B" + osp.split(a_filename)[-1][1:]
    fullbfile = osp.abspath(osp.join(osp.split(a_filename)[0], bfile))
    bfileexists = osp.isfile(fullbfile)
    if not bfileexists:
        raise IOError("B-file dynalog not found; ensure B-file is in same directory as A-file")
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

# ------------------------
# MLC Log Viewer Example
# ------------------------
if __name__ == '__main__':
    mlc = MachineLog()
    mlc.load_demo_file_dynalog()
    # mlc.load_logfile_UI()
    # pass
    # mlc.load_demo_file_trajectorylog()  # uncomment to use trajectory log file
    mlc.report_basic_parameters()

