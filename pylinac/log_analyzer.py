"""
The log analyzer module reads and parses Varian linear accelerator machine logs, both Dynalogs and Trajectory logs. The module also
calculates actual and expected fluences as well as performing gamma evaluations. Data is structured to be easily accessible and
easily plottable.

Unlike most other modules of pylinac, the log analyzer module has no end goal. Data is parsed from the logs, but what is done with that
info, and which info is analyzed is up to the user.

Features:

* **Analyze Dynalogs or Trajectory logs** - Either platform is supported. Tlog versions 2.1 and 3.0 supported.
* **Read in both .bin and .txt Trajectory log files** - Read in the machine data from both .bin and .txt files to get all the information recorded.
  See :attr:`~pylinac.log_analyzer.MachineLog.txt`.
* **Save Trajectory log data to CSV** - The Trajectory log binary data format does not allow for easy export of data. Pylinac lets you do
  that so you can use Excel or other software that you use with Dynalogs.
* **Plot or analyze any axis** - Every data axis can be accessed and plotted: the actual, expected, and even the difference.
* **Calculate fluences and gamma** - Besides reading in the MLC positions, pylinac calculates the actual and expected fluence
  as well as the gamma map; DTA and threshold values are adjustable.
"""
from abc import ABCMeta, abstractproperty
import copy
import csv
from functools import lru_cache
from io import BytesIO
import os
import os.path as osp
import struct
import warnings

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage.filters as spf

from pylinac.core.decorators import type_accept, value_accept
from pylinac.core.io import is_valid_file, is_valid_dir, get_folder_UI, get_filepath_UI, open_file
from pylinac.core.utilities import is_iterable, import_mpld3, get_url

np.seterr(invalid='ignore')  # ignore warnings for invalid numpy operations. Used for np.where() operations on partially-NaN arrays.

DYNALOG = 'Dynalog'
TRAJECTORY_LOG = 'Trajectory log'
STATIC_IMRT = 'Static IMRT'
DYNAMIC_IMRT = 'Dynamic IMRT'
VMAT = 'VMAT'


class MachineLogs(list):
    """
    .. versionadded:: 0.4.1

    Read in machine logs from a directory. Inherits from list. Batch methods are also provided."""
    @type_accept(folder=str)
    def __init__(self, folder=None, recursive=True, verbose=True):
        """
        Parameters
        ----------
        folder : str
            The directory of interest. Will walk through and process any logs, Trajectory or dynalog, it finds.
            Non-log files will be skipped.
        recursive : bool
            Whether to walk through subfolders of passed directory. Only used if ``folder`` is a valid log directory.
        verbose : bool
            If True (default), prints load status at each log. Only used if ``folder`` is a valid log directory.

        Examples
        --------
        Load a directory upon initialization::

            >>> log_folder = r'C:\path\log\directory'
            >>> logs = MachineLogs(log_folder)

        Or load them later directly::

            >>> logs = MachineLogs()
            >>> logs.load_folder(log_folder)

        Or even using a UI dialog box::

            >>> logs = MachineLogs.from_folder_UI()

        Batch methods include determining the average gamma and average gamma pass value::

            >>> logs.avg_gamma()
            >>> 0.05 # or whatever it is
            >>> logs.avg_gamma_pct()
            >>> 97.2
        """
        super().__init__()
        if folder is not None and is_valid_dir(folder):
            self.load_folder(folder, recursive, verbose)

    @property
    def num_logs(self):
        """Return the number of logs currently loaded."""
        return len(self)

    def _num_log_type(self, log_type):
        num = 0
        for log in self:
            if log.log_type == log_type:
                num += 1
        return num

    @property
    def num_tlogs(self):
        """Return the number of Trajectory logs currently loaded."""
        return self._num_log_type(TRAJECTORY_LOG)

    @property
    def num_dlogs(self):
        """Return the number of Trajectory logs currently loaded."""
        return self._num_log_type(DYNALOG)

    def load_folder(self, dir, recursive=True, verbose=True):
        """Load log files from a directory.

        Parameters
        ----------
        dir : str, None
            The directory of interest.
            If a string, will walk through and process any logs, Trajectory or dynalog, it finds.
            Non-log files will be skipped.
            If None, files must be loaded later using .load_dir() or .append().
        recursive : bool
            If True (default), will walk through subfolders of passed directory.
            If False, will only search root directory.
        verbose : bool
            If True (default), prints load status at each log.
        """
        # do initial walk to get file count
        num_logs = 0
        num_skipped = 0
        for root, dirs, files in os.walk(dir):
            _, num_logs, num_skipped = self._clean_log_filenames(files, root, num_logs, num_skipped)
            if not recursive:
                break
        if num_logs == 0:
            warnings.warn("No logs found.")
            return

        # actual log loading
        load_num = 1
        if verbose:
            print("{} logs found. \n{} logs skipped. \nLog loaded:".format(num_logs, num_skipped))
        for root, dirs, files in os.walk(dir):
            cleaned_files, _, _ = self._clean_log_filenames(files, root, num_logs, num_skipped)
            for name in cleaned_files:
                pth = osp.join(root, name)
                self.append(pth)
                if verbose:
                    print("{} of {}".format(load_num, num_logs))
                    load_num += 1
            if not recursive:
                break  # break out of for loop after top level search

    @classmethod
    def from_folder_UI(cls, recursive=True, verbose=True):
        """Construct a MachineLogs instance and load a folder from a UI dialog box.

        .. versionadded:: 0.6
        """
        obj = cls()
        obj.load_folder_UI(recursive, verbose)
        return obj

    def load_folder_UI(self, recursive=True, verbose=True):
        """Load a directory using a UI dialog box. See load_folder() for parameter info."""
        folder = get_folder_UI()
        if folder:
            self.load_folder(folder, recursive, verbose)

    def _clean_log_filenames(self, filenames, root, num_logs, num_skipped):
        """Extract the names of real log files from a list of files."""
        cleaned_filenames = []
        for idx, filename in enumerate(filenames):
            filename = osp.join(root, filename)
            if is_tlog(filename):
                num_logs += 1
                cleaned_filenames.append(filename)
            elif is_dlog(filename):
                opp_file = _return_other_dlg(filename, raise_find_error=False)
                if opp_file is not None:
                    num_logs += 1
                    opp_file = osp.basename(opp_file)
                    del filenames[filenames.index(opp_file)]
                    cleaned_filenames.append(filename)
                else:
                    num_skipped += 1
                    del filenames[idx]
        return cleaned_filenames, num_logs, num_skipped

    def _check_empty(self):
        """Check if any logs have been loaded."""
        if len(self) == 0:
            raise ValueError("No logs have been loaded yet.")

    def report_basic_parameters(self):
        """Report basic parameters of the logs.

        - Number of logs
        - Average gamma value of all logs
        - Average gamma pass percent of all logs
        """
        print("Number of logs: {}".format(self.num_logs))
        print("Average gamma: {:3.2f}".format(self.avg_gamma(verbose=False)))
        print("Average gamma pass percent: {:3.1f}".format(self.avg_gamma_pct(verbose=False)))

    def append(self, obj, recursive=True):
        """Append a log. Overloads list method.

        Parameters
        ----------
        obj : str, MachineLog
            If a string, must point to a log file.
            If a directory, must contain log files.
            If a MachineLog, then simply appends.
        recursive : bool
            Whether to walk through subfolders of passed directory. Only applicable if obj was a directory.
        """
        if isinstance(obj, str):
            if is_log(obj):
                log = MachineLog(obj)
                super().append(log)
            elif is_valid_dir(obj, raise_error=False):
                for root, dirs, files in os.walk(obj):
                    for name in files:
                        pth = osp.join(root, name)
                        self.append(pth)
        elif isinstance(obj, MachineLog):
            super().append(obj)
        else:
            raise TypeError("Can only append MachineLog or string pointing to a log or log directory.")

    def avg_gamma(self, doseTA=1, distTA=1, threshold=10, resolution=0.1, verbose=True):
        """Calculate and return the average gamma of all logs. See :meth:`~pylinac.log_analyzer.GammaFluence.calc_map()`
        for further parameter info."""
        self._check_empty()
        gamma_list = np.zeros(self.num_logs)
        if verbose:
            print("Calculating gammas:")
        for num, log in enumerate(self):
            log.fluence.gamma.calc_map(doseTA, distTA, threshold, resolution)
            gamma_list[num] = log.fluence.gamma.avg_gamma
            if verbose:
                print('{} of {}'.format(num+1, self.num_logs))
        return gamma_list.mean()

    def avg_gamma_pct(self, doseTA=1, distTA=1, threshold=10, resolution=0.1, verbose=True):
        """Calculate and return the average gamma pass percent of all logs. See :meth:`~pylinac.log_analyzer.GammaFluence.calc_map()`
        for further parameter info."""
        self._check_empty()
        gamma_list = np.zeros(self.num_logs)
        if verbose:
            print("Calculating gamma pass percent:")
        for num, log in enumerate(self):
            log.fluence.gamma.calc_map(doseTA, distTA, threshold, resolution)
            gamma_list[num] = log.fluence.gamma.pass_prcnt
            if verbose:
                print("{} of {}".format(num+1, self.num_logs))
        return gamma_list.mean()

    def to_csv(self):
        """Write trajectory logs to CSV. If there are both dynalogs and trajectory logs,
        only the trajectory logs will be written. File names will be the same as the original log file names."""
        num_written = 0
        for log in self:
            if is_tlog(log.filename):
                log.to_csv()
                num_written += 1
        if num_written:
            print('\n\nAll CSV files written!')
        else:
            print('\n\nNo files written')


class MachineLog:
    """Reads in and analyzes MLC log files, both dynalog and trajectory logs, from Varian linear accelerators.

    If reading Trajectory logs, the .txt file is also loaded if it's around.
    """
    def __init__(self, filename='', exclude_beam_off=True):
        """
        Parameters
        ----------
        filename : str
            Path to the log file. For trajectory logs this is a single .bin file.
            For dynalog files, load either the A*.dlg or B*.dlg file.

            .. warning:: Dynalogs must have names like "A*.dlg" and "B*.dlg" and be in the same folder, otherwise errors will ensue.

        exclude_beam_off : boolean
            If True (default), snapshots where the beam was not on will be removed.
            If False, all data will be included.

        Examples
        --------
        Load a file upon initialization::

            >>> mylogfile = "C:/path/to/log.bin"
            >>> log = MachineLog(mylogfile)

        Or::

            >>> mylogfile = "C:/path/to/A1.dlg"  # B must also be here
            >>> log = MachineLog()
            >>> log.load(mylogfile)

        Or load from a UI dialog::

            >>> log = MachineLog.from_UI()

        Run the demo::

            >>> MachineLog().run_dlog_demo()
            >>> MachineLog().run_tlog_demo()

        Attributes
        ----------
        filename : str
            The name of the file loaded in.
        header : A :class:`~pylinac.log_analyzer.Tlog_Header` or :class:`~pylinac.log_analyzer.Dlog_Header`, depending on the log type.
        axis_data : :class:`~pylinac.log_analyzer.Tlog_Axis_Data` or :class:`~pylinac.log_analyzer.Dlog_Axis_Data`, depending on the log type.
        subbeams : :class:`~pylinac.log_analyzer.SubbeamHandler`
            Will contain instances of :class:`~pylinac.log_analyzer.Subbeam`; will be empty if
            autosequencing was not done on v2.1.
        txt : dict
            If working with trajectory logs and the associated .txt file is in the same directory, the txt file data will be loaded as a dictionary.
        fluence : :class:`~pylinac.log_analyzer.Fluence_Struct`
            Contains actual and expected fluence data, including gamma.
        log_type : str
            The log type loaded; either 'Dynalog' or 'Trajectory log'
        log_is_loaded : bool
            Whether a log has yet been loaded.
        """
        self.filename = ''
        self.url = None
        self._cursor = 0
        self.fluence = FluenceStruct()

        # Read file if passed in
        if filename is not '':
            self.load(filename, exclude_beam_off)

    def run_tlog_demo(self):
        """Run the Trajectory log demo."""
        self.load_demo_trajectorylog()
        self.report_basic_parameters()
        self.plot_summary()

    def run_dlog_demo(self):
        """Run the dynalog demo."""
        self.load_demo_dynalog()
        self.report_basic_parameters()
        self.plot_summary()

    @classmethod
    def from_demo_dynalog(cls, exclude_beam_off=True):
        obj = cls()
        obj.load_demo_dynalog(exclude_beam_off)
        return obj

    def load_demo_dynalog(self, exclude_beam_off=True):
        """Load the demo dynalog file included with the package."""
        self.filename = osp.join(osp.dirname(__file__), 'demo_files', 'log_reader', 'AQA.dlg')
        self._read_log(exclude_beam_off)

    @classmethod
    def from_demo_trajectorylog(cls, exclude_beam_off=True):
        obj = cls()
        obj.load_demo_trajectorylog(exclude_beam_off)
        return obj

    def load_demo_trajectorylog(self, exclude_beam_off=True):
        """Load the demo trajectory log included with the package."""
        filename = osp.join(osp.dirname(__file__), 'demo_files', 'log_reader', 'Tlog.bin')
        self.load(filename, exclude_beam_off)

    @classmethod
    def from_url(cls, url, exclude_beam_off=True):
        """Instantiate a log from a URL.

        .. versionadded:: 0.7.1
        """
        obj = cls()
        obj.load_url(url, exclude_beam_off)
        return obj

    def load_url(self, url, exclude_beam_off=True):
        """Load a log from a URL.

        .. versionadded:: 0.7.1
        """
        response = get_url(url)
        stream = BytesIO(response.content)
        self.load(stream, exclude_beam_off)
        self.url = url

    @classmethod
    def from_UI(cls, exclude_beam_off=True):
        """Construct an instance and load a log from a UI dialog box.

        .. versionadded:: 0.6
        """
        obj = cls()
        obj.load_UI(exclude_beam_off)
        return obj

    def load_UI(self, exclude_beam_off=True):
        """Let user load a log file with a UI dialog box. """
        filename = get_filepath_UI()
        if filename: # if user didn't hit cancel...
            self.load(filename, exclude_beam_off)

    def load(self, filename, exclude_beam_off=True):
        """Load the log file directly by passing the path to the file.

        Parameters
        ----------
        filename : str
            The path to the log file.

        exclude_beam_off : boolean
            If True (default), snapshots where the beam was not on will be removed.
            If False, all data will be included.

            .. warning::
                Including beam off data may affect fluence and gamma results. E.g. in a step-&-shoot IMRT
                delivery, leaves move between segments while the beam is held. If beam-off data is included,
                the RMS and fluence errors may not correspond to what was delivered.
        """
        if is_valid_file(filename):
            if is_log(filename):
                self.filename = filename
                self._read_log(exclude_beam_off)
            else:
                raise IOError("File passed is not a valid log file")

    def plot_summary(self, show=True):
        """Plot actual & expected fluence, gamma map, gamma histogram,
            MLC error histogram, and MLC RMS histogram.
        """
        if self.fluence.gamma.map_calced:
            # plot the actual fluence
            ax = plt.subplot(2, 3, 1)
            ax.set_title('Actual Fluence', fontsize=10)
            ax.tick_params(axis='both', labelsize=8)
            ax.autoscale(tight=True)
            plt.imshow(self.fluence.actual.pixel_map, aspect='auto', interpolation='none')

            # plot the expected fluence
            ax = plt.subplot(2, 3, 2)
            ax.set_title("Expected Fluence", fontsize=10)
            ax.tick_params(axis='both', labelsize=8)
            plt.imshow(self.fluence.expected.pixel_map, aspect='auto', interpolation='none')

            # plot the gamma map
            gmma = self.fluence.gamma
            ax = plt.subplot(2, 3, 3)
            ax.set_title("Gamma Map ({:2.2f}% passing @ {}%/{}mm)".format(gmma.pass_prcnt, gmma.doseTA, gmma.distTA), fontsize=10)
            ax.tick_params(axis='both', labelsize=8)
            plt.imshow(self.fluence.gamma.pixel_map, aspect='auto', interpolation='none', vmax=1)
            plt.colorbar(ax=ax)

            # plot the gamma histogram
            ax = plt.subplot(2, 3, 4)
            ax.set_yscale('log')
            ax.set_title("Gamma Histogram (Avg: {:2.3f})".format(self.fluence.gamma.avg_gamma), fontsize=10)
            ax.tick_params(axis='both', labelsize=8)
            plt.hist(self.fluence.gamma.pixel_map.flatten(), self.fluence.gamma.bins)

            # plot the MLC error histogram
            ax = plt.subplot(2, 3, 5)
            p95error = self.axis_data.mlc.get_error_percentile()
            ax.set_title("Leaf Error Histogram (95th Perc: {:2.3f}cm)".format(p95error), fontsize=10)
            ax.tick_params(axis='both', labelsize=8)
            plt.hist(self.axis_data.mlc._abs_error_all_leaves.flatten())

            # plot the leaf RMSs
            ax = plt.subplot(2,3,6)
            ax.set_title("Leaf RMS Error (Max: {:2.3f}cm)".format(self.axis_data.mlc.get_RMS_max()), fontsize=10)
            ax.tick_params(axis='both', labelsize=8)
            ax.set_xlim([-0.5, self.axis_data.mlc.num_leaves+0.5])  # bit of padding since bar chart alignment is center
            plt.bar(np.arange(len(self.axis_data.mlc.get_RMS('both')))[::-1], self.axis_data.mlc.get_RMS('both'), align='center')

            if show:
                plt.show()
        else:
            raise AttributeError("Gamma map has not yet been calculated.")

    def save_summary(self, filename, **kwargs):
        """Save the summary image to file."""
        self.plot_summary(show=False)
        plt.savefig(filename, **kwargs)

    def report_basic_parameters(self, printout=True):
        """Print the common parameters analyzed when investigating machine logs:

        - Log type
        - Average MLC RMS
        - Maximum MLC RMS
        - 95th percentile MLC error
        - Number of beam holdoffs
        - Gamma pass percentage
        - Average gamma value
        """
        log_type = "MLC log type: {}\n".format(self.log_type)
        avg_rms = "Average RMS of all leaves: {:3.3f} cm\n".format(self.axis_data.mlc.get_RMS_avg(only_moving_leaves=False))
        max_rms = "Max RMS error of all leaves: {:3.3f} cm\n".format(self.axis_data.mlc.get_RMS_max())
        p95 = "95th percentile error: {:3.3f} cm\n".format(self.axis_data.mlc.get_error_percentile(95, only_moving_leaves=False))
        num_holdoffs = "Number of beam holdoffs: {:1.0f}\n".format(self.axis_data.num_beamholds)
        self.fluence.gamma.calc_map()
        gamma_pass = "Gamma pass %: {:2.2f}\n".format(self.fluence.gamma.pass_prcnt)
        gamma_avg = "Gamma average: {:2.3f}\n".format(self.fluence.gamma.avg_gamma)

        string = log_type + avg_rms + max_rms + p95 + num_holdoffs + gamma_pass + gamma_avg
        if printout:
            print(string)
        else:
            return string

    @property
    @lru_cache()
    def log_type(self):
        """Determine the MLC log type: Trajectory or Dynalog.

        The file is opened and sampled of
        the first few bytes. If the sample matches what is found in standard dynalog or tlog files it will be
        assigned that log type.

        Returns
        -------
        str : {'Dynalog', 'Trajectory Log'}

        Raises
        ------
        ValueError : If log type cannot be determined.
        """
        if is_dlog(self.filename):
            log_type = DYNALOG
        elif is_tlog(self.filename):
            log_type = TRAJECTORY_LOG
        return log_type

    @property
    def is_loaded(self):
        """Boolean specifying if a log has been loaded in yet."""
        if self.filename == '':
            return False
        else:
            return True

    @property
    def treatment_type(self):
        try:  # tlog
            gantry_std = self.subbeams[0].gantry_angle.actual.std()
        except AttributeError:  # dlog
            gantry_std = self.axis_data.gantry.actual.std()
        if gantry_std > 0.5:
            return 'VMAT'
        elif self.axis_data.mlc.get_RMS_avg(only_moving_leaves=True) > 0.05:
            return 'Dynamic IMRT'
        else:
            return 'Static IMRT'

    @property
    def _filename_str(self):
        if isinstance(self.filename, str):
            return self.filename
        else:
            return osp.basename(self.url)

    def to_csv(self, filename=None):
        """Write the log to a CSV file. Only applicable for Trajectory logs (Dynalogs are already similar to CSV).

        Parameters
        ----------
        filename : None, str
            If None (default), the CSV filename will be the same as the filename of the log.
            If a string, the filename will be named so.
        """
        is_file_object = False
        if not is_tlog(self.filename):
            raise TypeError("Writing to CSV is only applicable to trajectory logs.")
        if filename is None:
            filename = self.filename.replace('bin', 'csv')
        else:
            try:
                if not filename.endswith('.csv'):
                    filename = filename + '.csv'
            except AttributeError:
                is_file_object = True

        csv_file = open_file(filename, 'w')
        writer = csv.writer(csv_file, lineterminator='\n')
        # write header info
        header_titles = ('Tlog File:', 'Signature:', 'Version:', 'Header Size:', 'Sampling Inteval:',
                         'Number of Axes:', 'Axis Enumeration:', 'Samples per Axis:', 'Axis Scale:',
                         'Number of Subbeams:', 'Is Truncated?', 'Number of Snapshots:', 'MLC Model:')
        h = self.header
        header_values = (self._filename_str, h.header, h.version, h.header_size, h.sampling_interval,
                         h.num_axes, h.axis_enum, h.samples_per_axis, h.axis_scale, h.num_subbeams, h.is_truncated,
                         h.num_snapshots, h.mlc_model)
        for title, value in zip(header_titles, header_values):
            write_single_value(writer, title, value)

        # write axis data
        data_titles = ('Gantry', 'Collimator', 'Couch Lat', 'Couch Lng', 'Couch Rtn', 'MU',
                       'Beam Hold', 'Control Point', 'Carriage A', 'Carriage B')
        ad = self.axis_data
        data_values = (ad.gantry, ad.collimator, ad.couch.latl, ad.couch.long, ad.couch.rotn,
                       ad.mu, ad.beam_hold, ad.control_point, ad.carriage_A, ad.carriage_B)
        data_units = ('degrees', 'degrees', 'cm', 'cm', 'degrees', 'MU', None, None, 'cm',
                      'cm')
        for title, value, unit in zip(data_titles, data_values, data_units):
            write_array(writer, title, value, unit)

        # write leaf data
        for leaf_num, leaf in self.axis_data.mlc.leaf_axes.items():
            write_array(writer, 'Leaf ' + str(leaf_num), leaf, 'cm')

        if not is_file_object:
            print("CSV file written to: " + filename)

    def _read_log(self, exclude_beam_off):
        """Read in log based on what type of log it is: Trajectory or Dynalog."""
        if is_tlog(self.filename):
            self._read_tlog(exclude_beam_off)
        elif is_dlog(self.filename):
            self._read_dlog(exclude_beam_off)

    def _read_dlog(self, exclude_beam_off):
        """Read in Dynalog files from .dlg files (which are renamed CSV files).
        Formatting follows from the Dynalog File Viewer Reference Guide.
        """
        # if file is B*.dlg, replace with A*.dlg
        other_dlg_file = _return_other_dlg(self.filename)

        # create iterator object to read in lines
        with open(self.filename) as csvf:
            dlgdata = csv.reader(csvf, delimiter=',')
            self.header, dlgdata = DlogHeader(dlgdata)._read()
            self.axis_data = DlogAxisData(dlgdata, self.header, other_dlg_file)._read(exclude_beam_off)

        self.fluence = FluenceStruct(self.axis_data.mlc, self.axis_data.mu, self.axis_data.jaws)

    def _read_tlog(self, exclude_beam_off):
        """Read in Trajectory log from binary file according to TB 1.5/2.0 (i.e. Tlog v2.1/3.0) log file specifications."""
        # read in associated *.txt file if in the same directory.
        if is_tlog_txt_file_around(self.filename):
            self._read_txt_file()

        # read in trajectory log binary data
        fcontent = open_file(self.filename).read()

        # Unpack the content according to respective section and data type (see log specification file).
        self.header, self._cursor = TlogHeader(fcontent, self._cursor)._read()

        self.subbeams, self._cursor = SubbeamHandler(fcontent, self._cursor, self.header)._read()

        self.axis_data, self._cursor = TlogAxisData(fcontent, self._cursor, self.header)._read(exclude_beam_off)

        # self.crc = CRC(fcontent, self._cursor).read()

        self.fluence = FluenceStruct(self.axis_data.mlc, self.axis_data.mu, self.axis_data.jaws)

        self.subbeams.post_hoc_metadata(self.axis_data)

    def _read_txt_file(self):
        """Read a Tlog's associated .txt file and put in under the 'txt' attribute."""
        self.txt = {}
        txt_filename = self.filename.replace('.bin', '.txt')
        with open(txt_filename) as csvfile:
            txt_reader = csv.reader(csvfile, delimiter='\n')
            for row in txt_reader:
                if row and isinstance(row, list):
                    items = row[0].split(':', 1)
                    self.txt[items[0].strip()] = items[1].strip()


class Axis:
    """Represents an 'Axis' of a Trajectory log or dynalog file, holding actual and potentially expected and difference values.

    Attributes
    ----------
    Parameters are Attributes
    """
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
        self.expected = expected
        if expected is not None:
            try:
                if len(actual) != len(expected):
                    raise ValueError("Actual and expected Axis parameters are not equal length")
            except TypeError:
                pass
            self.expected = expected

    @property
    @lru_cache()
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

    def plot_actual(self, interactive=False):
        """Plot the actual positions as a matplotlib figure."""
        self._plot('actual', interactive=interactive)

    def save_plot_actual(self, filename, interactive=False, **kwargs):
        self._plot('actual', show=False)
        self._save(filename, interactive, **kwargs)

    def plot_expected(self, interactive=False):
        """Plot the expected positions as a matplotlib figure."""
        self._plot('expected', interactive=interactive)

    def save_plot_expected(self, filename, interactive=False, **kwargs):
        self._plot('expected', show=False)
        self._save(filename, interactive, **kwargs)

    def plot_difference(self, interactive=False):
        """Plot the difference of positions as a matplotlib figure."""
        self._plot('difference', interactive=interactive)

    def save_plot_difference(self, filename, interactive=False, **kwargs):
        self._plot('difference', show=False)
        self._save(filename, interactive, **kwargs)

    @value_accept(param=('actual', 'expected', 'difference'))
    def _plot(self, param='', interactive=False, show=True):
        """Plot the parameter: actual, expected, or difference."""
        plt.plot(getattr(self, param))
        plt.autoscale(axis='x', tight=True)
        if show:
            if interactive:
                mpld3 = import_mpld3()
                mpld3.show()
            else:
                plt.show()

    def _save(self, filename, interactive, **kwargs):
        """Save the figure to a file, either .png or .html."""
        if interactive:
            mpld3 = import_mpld3()
            mpld3.save_html(plt.gcf(), filename)
        else:
            plt.savefig(filename, **kwargs)


class AxisMovedMixin:
    """Mixin class for Axis."""
    @property
    @lru_cache()
    def moved(self):
        """Return whether the axis moved during treatment."""
        threshold = 0.003
        if np.std(self.actual) > threshold:
            return True
        else:
            return False


class LeafAxis(Axis, AxisMovedMixin):
    """Axis holding leaf information."""
    def __init__(self, actual, expected):
        # force expected argument to be supplied
        super().__init__(actual, expected)


class GantryAxis(Axis, AxisMovedMixin):
    """Axis holding gantry information."""
    pass


class HeadAxis(Axis, AxisMovedMixin):
    """Axis holding head information (e.g. jaw positions, collimator)."""
    pass


class CouchAxis(Axis, AxisMovedMixin):
    """Axis holding couch information."""
    pass


class BeamAxis(Axis):
    """Axis holding beam information (e.g. MU, beam hold status)."""
    pass


class Fluence(metaclass=ABCMeta):
    """An abstract base class to be used for the actual and expected fluences.

    Attributes
    ----------
    pixel_map : numpy.ndarray
        An array representing the fluence map; will be num_mlc_pairs-x-400/resolution.
        E.g., assuming a Millennium 120 MLC model and a fluence resolution of 0.1mm, the resulting
        matrix will be 60-x-4000.
    resolution : int, float
        The resolution of the fluence calculation; -1 means calculation has not been done yet.
    """
    pixel_map = object
    resolution = -1
    _fluence_type = ''  # must be specified by subclass

    def __init__(self, mlc_struct=None, mu_axis=None, jaw_struct=None):
        """
        Parameters
        ----------
        mlc_struct : MLC_Struct
        mu_axis : BeamAxis
        jaw_struct : Jaw_Struct
        """
        self._mlc = mlc_struct
        self._mu = mu_axis
        self._jaws = jaw_struct

    @property
    def map_calced(self):
        """Return a boolean specifying whether the fluence has been calculated."""
        if hasattr(self.pixel_map, 'size'):
            return True
        else:
            return False

    @lru_cache()
    def calc_map(self, resolution=0.1):
        """Calculate a fluence pixel map.

        Fluence calculation is done by adding fluence snapshot by snapshot, and leaf pair by leaf pair.
        Each leaf pair is analyzed separately. First, to optimize, it checks if the leaf is under the y-jaw.
        If so, the fluence is left at zero; if not, the leaf (or jaw) ends are determined and the MU fraction of that
        snapshot is added to the total fluence. All snapshots are iterated over for each leaf pair until the total fluence
        matrix is built.

        Parameters
        ----------
        resolution : int, float
            The resolution in mm of the fluence calculation in the leaf-moving direction.

         Returns
         -------
         numpy.ndarray
             A numpy array reconstructing the actual fluence of the log. The size will
             be the number of MLC pairs by 400 / resolution since the MLCs can move anywhere within the
             40cm-wide linac head opening.
         """
        # return if map has already been calculated under the same conditions
        # if self.map_calced and self._same_conditions(resolution):
        #     return self.pixel_map
        # preallocate arrays for expected and actual fluence of number of leaf pairs-x-4000 (40cm = 4000um, etc)
        fluence = np.zeros((self._mlc.num_pairs, int(400 / resolution)), dtype=np.float32)

        # calculate the MU delivered in each snapshot. For Tlogs this is absolute; for dynalogs it's normalized.
        mu_matrix = getattr(self._mu, self._fluence_type)
        MU_differential = np.zeros(len(mu_matrix))
        MU_differential[0] = mu_matrix[0]
        MU_differential[1:] = np.diff(mu_matrix)
        MU_differential = MU_differential / mu_matrix[-1]
        MU_cumulative = 1

        # calculate each "line" of fluence (the fluence of an MLC leaf pair, e.g. 1 & 61, 2 & 62, etc),
        # and add each "line" to the total fluence matrix
        fluence_line = np.zeros(int(400 / resolution), dtype=np.float32)
        leaf_offset = self._mlc.num_pairs
        pos_offset = int(np.round(200 / resolution))
        for pair in range(1, self._mlc.num_pairs + 1):
            if not self._mlc.leaf_under_y_jaw(pair):
                fluence_line[:] = 0  # emtpy the line values on each new leaf pair
                left_leaf_data = getattr(self._mlc.leaf_axes[pair], self._fluence_type)
                left_leaf_data = -np.round(left_leaf_data * 10 / resolution) + pos_offset
                right_leaf_data = getattr(self._mlc.leaf_axes[pair + leaf_offset], self._fluence_type)
                right_leaf_data = np.round(right_leaf_data * 10 / resolution) + pos_offset
                left_jaw_data = np.round((200 / resolution) - (self._jaws.x1.actual * 10 / resolution))
                right_jaw_data = np.round((self._jaws.x2.actual * 10 / resolution) + (200 / resolution))
                if self._mlc.pair_moved(pair):
                    for snapshot in self._mlc.snapshot_idx:
                        lt_mlc_pos = left_leaf_data[snapshot]
                        rt_mlc_pos = right_leaf_data[snapshot]
                        lt_jaw_pos = left_jaw_data[snapshot]
                        rt_jaw_pos = right_jaw_data[snapshot]
                        left_edge = int(max(lt_mlc_pos, lt_jaw_pos))
                        right_edge = int(min(rt_mlc_pos, rt_jaw_pos))
                        fluence_line[left_edge:right_edge] += MU_differential[snapshot]
                else:  # leaf didn't move; no need to calc over every snapshot
                    first_snapshot = self._mlc.snapshot_idx[0]
                    lt_mlc_pos = left_leaf_data[first_snapshot]
                    rt_mlc_pos = right_leaf_data[first_snapshot]
                    lt_jaw_pos = left_jaw_data.min()
                    rt_jaw_pos = right_jaw_data.max()
                    left_edge = max(lt_mlc_pos, lt_jaw_pos)
                    right_edge = min(rt_mlc_pos, rt_jaw_pos)
                    fluence_line[int(left_edge):int(right_edge)] = MU_cumulative
                fluence[pair - 1, :] = fluence_line

        self.pixel_map = fluence
        self.resolution = resolution
        return fluence

    def plot_map(self, show=True):
        """Plot the fluence; the fluence (pixel map) must have been calculated first."""
        if not self.map_calced:
            raise AttributeError("Map not yet calculated; use calc_map()")
        plt.clf()
        plt.imshow(self.pixel_map, aspect='auto')
        if show:
            plt.show()

    def save_map(self, filename, **kwargs):
        """Save the fluence map figure to a file."""
        self.plot_map(show=False)
        plt.savefig(filename, **kwargs)


class ActualFluence(Fluence):
    """The actual fluence object"""
    _fluence_type = 'actual'


class ExpectedFluence(Fluence):
    """The expected fluence object."""
    _fluence_type = 'expected'


class GammaFluence(Fluence):
    """Gamma object, including pixel maps of gamma, binary pass/fail pixel map, and others.

    Attributes
    ----------
    pixel_map : numpy.ndarray
        The gamma map. Only available after calling calc_map()
    passfail_map : numpy.ndarray
        The gamma pass/fail map; pixels that pass (<1.0) are set to 0, while failing pixels (>=1.0) are set to 1.
    distTA : int, float
        The distance to agreement value used in gamma calculation.
    doseTA : int, float
        The dose to agreement value used in gamma calculation.
    threshold : int, float
        The threshold percent dose value, below which gamma was not evaluated.
    pass_prcnt : float
        The percent of pixels passing gamma (<1.0).
    avg_gamma : float
        The average gamma value.
    """
    distTA = -1
    doseTA = -1
    threshold = -1
    pass_prcnt = -1
    avg_gamma = -1
    # doseTA_map = np.ndarray
    # distTA_map = np.ndarray
    passfail_map = object
    bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1]

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
        self._actual_fluence = actual_fluence
        self._expected_fluence = expected_fluence
        self._mlc = mlc_struct

    @lru_cache()
    def calc_map(self, doseTA=1, distTA=1, threshold=10, resolution=0.1, calc_individual_maps=False):
        """Calculate the gamma from the actual and expected fluences.

        The gamma calculation is based on `Bakai et al
        <http://iopscience.iop.org/0031-9155/48/21/006/>`_ eq.6,
        which is a quicker alternative to the standard Low gamma equation.

        Parameters
        ----------
        DoseTA : int, float
            Dose-to-agreement in percent; e.g. 2 is 2%.
        DistTA : int, float
            Distance-to-agreement in mm.
        threshold : int, float
            The dose threshold percentage of the maximum dose, below which is not analyzed.
        resolution : int, float
            The resolution in mm of the resulting gamma map in the leaf-movement direction.
        calc_individual_maps : bool
            Not yet implemented.
            If True, separate pixel maps for the distance-to-agreement and dose-to-agreement are created.

        Returns
        -------
        numpy.ndarray
            A num_mlc_leaves-x-400/resolution numpy array.
        """
        # calc fluences if need be
        if not self._actual_fluence.map_calced or resolution != self._actual_fluence.resolution:
            self._actual_fluence.calc_map(resolution)
        if not self._expected_fluence.map_calced or resolution != self._expected_fluence.resolution:
            self._expected_fluence.calc_map(resolution)

        actual = copy.copy(self._actual_fluence.pixel_map)
        expected = copy.copy(self._expected_fluence.pixel_map)

        # set dose values below threshold to 0 so gamma doesn't calculate over it
        actual[actual < (threshold / 100) * np.max(actual)] = 0
        expected[expected < (threshold / 100) * np.max(expected)] = 0

        # preallocation
        gamma_map = np.zeros(expected.shape)

        # image gradient in x-direction (leaf movement direction) using sobel filter
        img_x = spf.sobel(actual, 1)

        # equation: (measurement - reference) / sqrt ( doseTA^2 + distTA^2 * image_gradient^2 )
        for leaf in range(self._mlc.num_pairs):
            gamma_map[leaf] = np.abs(actual[leaf, :] - expected[leaf, :]) / np.sqrt(
                (doseTA / 100.0 ** 2) + ((distTA / resolution ** 2) * (img_x[leaf, :] ** 2)))
        # construct binary pass/fail map
        self.passfail_map = np.array(np.where(gamma_map >= 1, 1, 0)[0], dtype=bool)

        # if calc_individual_maps:
        #     # calculate DoseTA map (drops distTA calc from gamma eq)
        #     self.doseTA_map = np.zeros(gamma_map.shape)
        #     for leaf in range(self.mlc.num_pairs):
        #         self.doseTA_map[leaf] = (actual[leaf, :] - expected[leaf, :]) / np.sqrt(
        #             (DoseTA / 100.0 ** 2))
        #     # calculate DistTA map (drops DoseTA calc from gamma eq)
        #     self.distTA_map = np.zeros(gamma_map.shape)
        #     for leaf in range(self.mlc.num_pairs):
        #         self.distTA_map[leaf] = (actual[leaf, :] - expected[leaf, :]) / np.sqrt(
        #             ((DistTA / resolution ** 2) * (img_x[leaf, :] ** 2)))

        # calculate standard metrics
        self.pass_prcnt = (np.sum(gamma_map < 1) / np.sum(gamma_map >= 0)) * 100
        self.avg_gamma = np.nanmean(gamma_map)

        self.distTA = distTA
        self.doseTA = doseTA
        self.threshold = threshold
        self.resolution = resolution

        self.pixel_map = gamma_map
        return gamma_map

    def plot_map(self, show=True):
        """Plot the fluence; the fluence (pixel map) must have been calculated first."""
        if not self.map_calced:
            raise AttributeError("Map not yet calculated; use calc_map()")
        plt.imshow(self.pixel_map, aspect='auto', vmax=1)
        plt.colorbar()
        plt.show()

    def histogram(self, bins=None):
        """Return a histogram array and bin edge array of the gamma map values.

        Parameters
        ----------
        bins : sequence
            The bin edges for the gamma histogram; see numpy.histogram for more info.

        Returns
        -------
        histogram : numpy.ndarray
            A 1D histogram of the gamma values.
        bin_edges : numpy.ndarray
            A 1D array of the bin edges. If left as None, the class default will be used (self.bins).
        """
        if self.map_calced:
            if bins is None:
                bins = self.bins
            hist_arr, bin_edges = np.histogram(self.pixel_map, bins=bins)
            return hist_arr, bin_edges
        else:
            raise AttributeError("Gamma map not yet calculated")

    def plot_histogram(self, scale='log', bins=None, show=True):
        """Plot a histogram of the gamma map values.

        Parameters
        ----------
        scale : {'log', 'linear'}
            Scale of the plot y-axis.
        bins : sequence
            The bin edges for the gamma histogram; see numpy.histogram for more info.
        """
        if self.map_calced:
            if bins is None:
                bins = self.bins
            plt.clf()
            plt.hist(self.pixel_map.flatten(), bins=bins)
            plt.yscale(scale)
            if show:
                plt.show()
        else:
            raise AttributeError("Map not yet calculated; use calc_map()")

    def save_histogram(self, filename, scale='log', bins=None, **kwargs):
        """Save the histogram plot to file."""
        self.plot_histogram(scale, bins, show=False)
        plt.savefig(filename, **kwargs)

    def plot_passfail_map(self):
        """Plot the binary gamma map, only showing whether pixels passed or failed."""
        if self.map_calced:
            plt.imshow(self.passfail_map)
            plt.show()
        else:
            raise AttributeError("Map not yet calculated; use calc_map()")


class FluenceStruct:
    """Structure for data and methods having to do with fluences.

    Attributes
    ----------
    actual : :class:`~pylinac.log_analyzer.Fluence`
        The actual fluence delivered.
    expected : :class:`~pylinac.log_analyzer.Fluence`
        The expected, or planned, fluence.
    gamma : :class:`~pylinac.log_analyzer.GammaFluence`
        The gamma structure regarding the actual and expected fluences.
    """
    def __init__(self, mlc_struct=None, mu_axis=None, jaw_struct=None):
        self.actual = ActualFluence(mlc_struct, mu_axis, jaw_struct)
        self.expected = ExpectedFluence(mlc_struct, mu_axis, jaw_struct)
        self.gamma = GammaFluence(self.actual, self.expected, mlc_struct)


class MLC:
    """The MLC class holds MLC information and retrieves relevant data about the MLCs and positions."""
    def __init__(self, snapshot_idx=None, jaw_struct=None, HDMLC=False):
        """
        Parameters
        ----------
        snapshot_idx : array, list
            The snapshots to be considered for RMS and error calculations (can be all snapshots or just when beam was on).
        jaw_struct : Jaw_Struct
        HDMLC : boolean
            If False (default), indicates a regular MLC model (e.g. Millennium 120).
            If True, indicates an HD MLC model (e.g. Millennium 120 HD).

        Attributes
        ----------
        leaf_axes : dict containing :class:`~pylinac.log_analyzer.Axis`
            The dictionary is keyed by the leaf number, with the Axis as the value.

            .. warning:: Leaf numbers are 1-index based to correspond with Varian convention.
        """
        self.leaf_axes = {}
        self.snapshot_idx = snapshot_idx
        self._jaws = jaw_struct
        self.hdmlc = HDMLC

    @property
    def num_pairs(self):
        """Return the number of MLC pairs."""
        return int(self.num_leaves/2)

    @property
    def num_leaves(self):
        """Return the number of MLC leaves."""
        return len(self.leaf_axes)

    @property
    @lru_cache()
    def num_snapshots(self):
        """Return the number of snapshots used for MLC RMS & Fluence calculations.

        .. warning::
            This number may not be the same as the number of recorded snapshots in the log
            since the snapshots where the beam was off may not be included. See :meth:`MachineLog.load`
        """
        return len(self.snapshot_idx)

    @property
    @lru_cache()
    def num_moving_leaves(self):
        """Return the number of leaves that moved."""
        return len(self.moving_leaves)

    @property
    @lru_cache()
    def moving_leaves(self):
        """Return an array of the leaves that moved during treatment."""
        threshold = 0.003
        indices = ()
        for leaf_num, leafdata in self.leaf_axes.items():
            leaf_stdev = np.std(leafdata.actual[self.snapshot_idx])
            if leaf_stdev > threshold:
                indices += (leaf_num,)
        return np.array(indices)

    @type_accept(leaf_axis=LeafAxis, leaf_num=int)
    def add_leaf_axis(self, leaf_axis, leaf_num):
        """Add a leaf axis to the MLC data structure.

        Parameters
        ----------
        leaf_axis : LeafAxis
            The leaf axis to be added.
        leaf_num : int
            The leaf number.

            .. warning:: Leaf numbers are 1-index based to correspond with Varian convention.
        """
        self.leaf_axes[leaf_num] = leaf_axis

    def leaf_moved(self, leaf_num):
        """Return whether the given leaf moved during treatment.

        Parameters
        ----------
        leaf_num : int


        .. warning:: Leaf numbers are 1-index based to correspond with Varian convention.
        """
        if leaf_num in self.moving_leaves:
            return True
        else:
            return False

    def pair_moved(self, pair_num):
        """Return whether the given pair moved during treatment.

        If either leaf moved, the pair counts as moving.

        Parameters
        ----------
        pair_num : int


        .. warning:: Pair numbers are 1-index based to correspond with Varian convention.
        """
        a_leaf = pair_num
        b_leaf = pair_num + self.num_pairs
        if self.leaf_moved(a_leaf) or self.leaf_moved(b_leaf):
            return True
        else:
            return False

    @property
    @lru_cache()
    def _all_leaf_indices(self):
        """Return an array enumerated over all the leaves."""
        return np.array(range(1, len(self.leaf_axes) + 1))

    def get_RMS_avg(self, bank='both', only_moving_leaves=False):
        """Return the overall average RMS of given leaves.

        Parameters
        ----------
        bank : {'A', 'B', 'both'}
            Specifies which bank(s) is desired.
        only_moving_leaves : boolean
            If False (default), include all the leaves.
            If True, will remove the leaves that were static during treatment.

            .. warning::
                The RMS and error will nearly always be lower if all leaves are included since non-moving leaves
                have an error of 0 and will drive down the average values. Convention would include all leaves,
                but prudence would use only the moving leaves to get a more accurate assessment of error/RMS.

        Returns
        -------
        float
        """
        leaves = self.get_leaves(bank, only_moving_leaves)
        rms_array = self.create_RMS_array(leaves)
        return np.mean(rms_array)

    def get_RMS_max(self, bank='both'):
        """Return the overall maximum RMS of given leaves.

        Parameters
        ----------
        bank : {'A', 'B', 'both'}
            Specifies which bank(s) is desired.

        Returns
        -------
        float
        """
        leaves = self.get_leaves(bank)
        rms_array = self.create_RMS_array(leaves)
        return np.max(rms_array)

    def get_RMS_percentile(self, percentile=95, bank='both', only_moving_leaves=False):
        """Return the n-th percentile value of RMS for the given leaves.

        Parameters
        ----------
        percentile : int
            RMS percentile desired.
        bank : {'A', 'B', 'both'}
            Specifies which bank(s) is desired.
        only_moving_leaves : boolean
            If False (default), include all the leaves.
            If True, will remove the leaves that were static during treatment.

            .. warning::
                The RMS and error will nearly always be lower if all leaves are included since non-moving leaves
                have an error of 0 and will drive down the average values. Convention would include all leaves,
                but prudence would use only the moving leaves to get a more accurate assessment of error/RMS.
        """
        leaves = self.get_leaves(bank, only_moving_leaves)
        rms_array = self.create_RMS_array(leaves)
        return np.percentile(rms_array, percentile)

    def get_RMS(self, leaves_or_bank):
        """Return an array of leaf RMSs for the given leaves or MLC bank.

        Parameters
        ----------
        leaves_or_bank : sequence of numbers, {'a', 'b', 'both'}
            If a sequence, must be a sequence of leaf numbers desired.
            If a string, it specifies which bank (or both) is desired.

        Returns
        -------
        numpy.ndarray
            An array for the given leaves containing the RMS error.
        """
        if isinstance(leaves_or_bank, str):
            leaves_or_bank = self.get_leaves(leaves_or_bank)
        elif not is_iterable(leaves_or_bank):
            raise TypeError("Input must be iterable, or specify an MLC bank")
        return self.create_RMS_array(np.array(leaves_or_bank))

    def get_leaves(self, bank='both', only_moving_leaves=False):
        """Return a list of leaves that match the given conditions.

        Parameters
        ----------
        bank : {'A', 'B', 'both'}
            Specifies which bank(s) is desired.
        only_moving_leaves : boolean
            If False (default), include all the leaves.
            If True, will remove the leaves that were static during treatment.
        """
        # get all leaves or only the moving leaves
        if only_moving_leaves:
            leaves = copy.copy(self.moving_leaves)
        else:
            leaves = copy.copy(self._all_leaf_indices)

        # select leaves by bank if desired
        if bank is not None:
            if bank.lower() == 'a':
                leaves = leaves[leaves <= self.num_pairs]
            elif bank.lower() == 'b':
                leaves = leaves[leaves > self.num_pairs]

        return leaves

    def get_error_percentile(self, percentile=95, bank='both', only_moving_leaves=False):
        """Calculate the n-th percentile error of the leaf error.

        Parameters
        ----------
        percentile : int
            RMS percentile desired.
        bank : {'A', 'B', 'both'}
            Specifies which bank(s) is desired.
        only_moving_leaves : boolean
            If False (default), include all the leaves.
            If True, will remove the leaves that were static during treatment.

            .. warning::
                The RMS and error will nearly always be lower if all leaves are included since non-moving leaves
                have an error of 0 and will drive down the average values. Convention would include all leaves,
                but prudence would use only the moving leaves to get a more accurate assessment of error/RMS.
        """
        leaves = self.get_leaves(bank, only_moving_leaves)
        leaves -= 1
        error_array = self.create_error_array(leaves)

        abs_error = np.abs(error_array)
        return np.percentile(abs_error, percentile)

    def create_error_array(self, leaves, absolute=True):
        """Create and return an error array of only the leaves specified.

        Parameters
        ----------
        leaves : sequence
            Leaves desired.
        absolute : bool
            If True, (default) absolute error will be returned.
            If False, error signs will be retained.

        Returns
        -------
        numpy.ndarray
            An array of size leaves-x-num_snapshots
        """
        if absolute:
            error_array = self._abs_error_all_leaves
        else:
            error_array = self._error_array_all_leaves
        return error_array[leaves, :]

    def create_RMS_array(self, leaves):
        """Create an RMS array of only the leaves specified.

        Parameters
        ----------
        leaves : sequence
            Leaves desired.

        Returns
        -------
        numpy.ndarray
            An array of size leaves-x-num_snapshots
        """
        rms_array = self._RMS_array_all_leaves
        leaves -= 1
        if len(leaves) == 0:
            return np.array([0])
        return rms_array[leaves]

    @property
    @lru_cache()
    def _abs_error_all_leaves(self):
        """Absolute error of all leaves."""
        return np.abs(self._error_array_all_leaves)

    @property
    @lru_cache()
    def _error_array_all_leaves(self):
        """Error array of all leaves."""
        mlc_error = np.zeros((self.num_leaves, self.num_snapshots))
        # construct numpy array for easy array calculation
        for leaf in range(self.num_leaves):
            mlc_error[leaf, :] = self.leaf_axes[leaf+1].difference[self.snapshot_idx]
        return mlc_error

    def _snapshot_array(self, dtype='actual'):
        """Return an array of the snapshot data of all leaves."""
        arr = np.zeros((self.num_leaves, self.num_snapshots))
        # construct numpy array for easy array calculation
        for leaf in range(self.num_leaves):
            arr[leaf, :] = getattr(self.leaf_axes[leaf + 1], dtype)[self.snapshot_idx]
        return arr

    @property
    @lru_cache()
    def _RMS_array_all_leaves(self):
        """Return the RMS of all leaves."""
        rms_array = np.array([np.sqrt(np.sum(leafdata.difference[self.snapshot_idx] ** 2) / self.num_snapshots) for leafdata in self.leaf_axes.values()])
        return rms_array

    def leaf_under_y_jaw(self, leaf_num):
        """Return a boolean specifying if the given leaf is under one of the y jaws.

        Parameters
        ----------
        leaf_num : int
        """
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
            else:  # between 50 and 70
                mlc_position += outer_leaf_thickness

        y2_position = self._jaws.y2.actual.max()*10 + 200
        y1_position = 200 - self._jaws.y1.actual.max()*10
        if 10 >= leaf or leaf >= 110:
            thickness = outer_leaf_thickness
        elif 50 >= leaf or leaf >= 70:
            thickness= inner_leaf_thickness
        else:  # between 50 and 70
            thickness = outer_leaf_thickness
        if mlc_position < y1_position or mlc_position - thickness > y2_position:
            return True
        else:
            return False

    def get_snapshot_values(self, bank_or_leaf='both', dtype='actual'):
        """Retrieve the snapshot data of the given MLC bank or leaf/leaves

        Parameters
        ----------
        bank_or_leaf : str, array, list
            If a str, specifies what bank ('A', 'B', 'both').
            If an array/list, specifies what leaves (e.g. [1,2,3])
        dtype : {'actual', 'expected'}
            The type of MLC snapshot data to return.

        Returns
        -------
        ndarray
            An array of shape (number of leaves - x - number of snapshots). E.g. for an MLC bank
            and 500 snapshots, the array would be (60, 500).
        """
        if isinstance(bank_or_leaf, str):
            leaves = self.get_leaves(bank=bank_or_leaf)
            leaves -= 1
        else:
            leaves = bank_or_leaf

        arr = self._snapshot_array(dtype)
        return arr[leaves, :]

    def plot_mlc_error_hist(self, show=True):
        """Plot an MLC error histogram."""
        plt.hist(self._abs_error_all_leaves.flatten())
        if show:
            plt.show()

    def save_mlc_error_hist(self, filename, **kwargs):
        """Save the MLC error histogram to file."""
        self.plot_mlc_error_hist(show=False)
        plt.savefig(filename, **kwargs)

    def plot_rms_by_leaf(self, show=True):
        """Plot RMSs by leaf."""
        plt.clf()
        plt.bar(np.arange(len(self.get_RMS('both')))[::-1], self.get_RMS('both'), align='center')
        if show:
            plt.show()

    def save_rms_by_leaf(self, filename, **kwargs):
        """Save the RMS-leaf to file."""
        self.plot_rms_by_leaf(show=False)
        plt.savefig(filename, **kwargs)


class JawStruct:
    """Jaw Axes data structure.

    Attributes
    ----------
    x1 : :class:`~pylinac.log_analyzer.Axis`
    y1 : :class:`~pylinac.log_analyzer.Axis`
    x2 : :class:`~pylinac.log_analyzer.Axis`
    y2 : :class:`~pylinac.log_analyzer.Axis`
    """
    def __init__(self, x1, y1, x2, y2):
        if not all((
                isinstance(x1, HeadAxis),
                isinstance(y1, HeadAxis),
                isinstance(x2, HeadAxis),
                isinstance(y2, HeadAxis))):
            raise TypeError("Head_Axis not passed into Jaw structure")
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2


class CouchStruct:
    """Couch Axes data structure.

    Attributes
    ----------
    vert : :class:`~pylinac.log_analyzer.Axis`
    long : :class:`~pylinac.log_analyzer.Axis`
    latl : :class:`~pylinac.log_analyzer.Axis`
    rotn : :class:`~pylinac.log_analyzer.Axis`
    """
    def __init__(self, vertical, longitudinal, lateral, rotational, roll=None, pitch=None):
        if not all((isinstance(vertical, CouchAxis),
                    isinstance(longitudinal, CouchAxis),
                    isinstance(lateral, CouchAxis),
                    isinstance(rotational, CouchAxis))):
            raise TypeError("Couch structure must be passed Couch Axes.")
        self.vert = vertical
        self.long = longitudinal
        self.latl = lateral
        self.rotn = rotational
        if roll is not None:
            self.roll = roll
            self.pitch = pitch


class LogSection(metaclass=ABCMeta):
    @abstractproperty
    def _read(self):
        pass


class DlogSection(LogSection, metaclass=ABCMeta):
    def __init__(self, log_content):
        self._log_content = log_content


class TlogSection(LogSection, metaclass=ABCMeta):
    def __init__(self, log_content, cursor):
        self._log_content = log_content
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
            if len(output) == 1:
                output = int(output)
            self._cursor += ssize
        elif dtype == float:
            ssize = struct.calcsize('f') * num_values
            output = np.asarray(struct.unpack('f' * num_values, fc[self._cursor:self._cursor + ssize]))
            if len(output) == 1:
                output = float(output)
            self._cursor += ssize
        else:
            raise TypeError("decode_binary datatype was not valid")

        self._cursor += cursor_shift  # shift cursor if need be (e.g. if a reserved section follows)
        return output


class Subbeam(TlogSection):
    """Data structure for trajectory log "subbeams". Only applicable for auto-sequenced beams.

    Attributes
    ----------
    control_point : int
        Internally-defined marker that defines where the plan is currently executing.
    mu_delivered : float
        Dose delivered in units of MU.
    rad_time : float
        Radiation time in seconds.
    sequence_num : int
        Sequence number of the subbeam.
    beam_name : str
        Name of the subbeam.
    """
    def __init__(self, log_content, cursor, log_version):
        super().__init__(log_content, cursor)
        self._log_version = log_version

    def _read(self):
        """Read the tlog subbeam information."""
        self.control_point = self._decode_binary(self._log_content, int)
        self.mu_delivered = self._decode_binary(self._log_content, float)
        self.rad_time = self._decode_binary(self._log_content, float)
        self.sequence_num = self._decode_binary(self._log_content, int)
        # In Tlogs version 3.0 and up, beam names are 512 byte unicode strings, but in <3.0 they are 32 byte unicode strings
        if self._log_version >= 3:
            chars = 512
        else:
            chars = 32
        self.beam_name = self._decode_binary(self._log_content, str, chars, 32)

        return self, self._cursor

    @property
    def gantry_angle(self):
        """Median gantry angle of the subbeam."""
        return self._get_metadata_axis('gantry')

    @property
    def collimator_angle(self):
        """Median collimator angle of the subbeam."""
        return self._get_metadata_axis('collimator')

    @property
    def jaw_x1(self):
        """Median X1 position of the subbeam."""
        return self._get_metadata_axis('jaws', 'x1')

    @property
    def jaw_x2(self):
        """Median X2 position of the subbeam."""
        return self._get_metadata_axis('jaws', 'x2')

    @property
    def jaw_y1(self):
        """Median Y1 position of the subbeam."""
        return self._get_metadata_axis('jaws', 'y1')

    @property
    def jaw_y2(self):
        """Median Y2 position of the subbeam."""
        return self._get_metadata_axis('jaws', 'y2')

    def _get_metadata_axis(self, attr, subattr=None):
        if subattr is None:
            actual = getattr(self._axis_data, attr).actual[self._snapshots]
            expected = getattr(self._axis_data, attr).expected[self._snapshots]
        else:
            actual = getattr(getattr(self._axis_data, attr), subattr).actual[self._snapshots]
            expected = getattr(getattr(self._axis_data, attr), subattr).expected[self._snapshots]
        return Axis(np.median(actual), np.median(expected))


class SubbeamHandler:
    """One of 4 subsections of a trajectory log. Holds a list of Subbeams; only applicable for auto-sequenced beams."""
    def __init__(self, log_content, cursor, header):
        self._log_content = log_content
        self._header = header
        self._cursor = cursor
        self.subbeams = []

    def _read(self):
        """Read all the subbeams of a tlog file.

        Returns
        -------
        list
            Contains instances of Subbeam.
        cursor : int
            Internal; for tracking the cursor position in the file.
        """
        if self._header.num_subbeams > 0:
            for subbeam_num in range(self._header.num_subbeams):
                subbeam, self._cursor = Subbeam(self._log_content, self._cursor, self._header.version)._read()
                self.subbeams.append(subbeam)
        return self, self._cursor

    def post_hoc_metadata(self, axis_data):
        """From the Axis Data, perform post-hoc analysis and set metadata to the subbeams.
            Gives the subbeams more information, as not much is given directly in the logs."""
        for idx, beam in enumerate(self.subbeams):
            self._set_beamon_snapshots(axis_data, idx)

    def _set_beamon_snapshots(self, axis_data, beam_num):
        """Get the snapshot indices 1) where the beam was on and 2) between the subbeam control point values."""
        beam = self.subbeams[beam_num]
        all_snapshots = axis_data.control_point.actual
        lower_bound = beam.control_point
        try:
            upper_bound = self.subbeams[beam_num+1].control_point
        except IndexError:
            upper_bound = axis_data.control_point.actual.max()
        section_snapshots = np.logical_and(all_snapshots>=lower_bound, all_snapshots<upper_bound)
        beam_on_snapshots = axis_data.beam_hold.actual == 0
        combined_snapshot = np.logical_and(beam_on_snapshots, section_snapshots)

        beam._snapshots = combined_snapshot
        beam._axis_data = axis_data

    def __getitem__(self, item):
        return self.subbeams[item]


class TlogHeader(TlogSection):
    """A header object, one of 4 sections of a trajectory log. Holds sampling interval, version, etc.

    Attributes
    ----------
    header : str
        Header signature: 'VOSTL'.
    version : str
        Log version.
    header_size : int
        Header size; fixed at 1024.
    sampling_interval : int
        Sampling interval in milliseconds.
    num_axes : int
        Number of axes sampled.
    axis_enum : int
        Axis enumeration; see the Tlog file specification for more info.
    samples_per_axis : numpy.ndarray
        Number of samples per axis; 1 for most axes, for MLC it's # of leaves and carriages.
    num_mlc_leaves : int
        Number of MLC leaves.
    axis_scale : int
        Axis scale; 1 -> Machine scale, 2 -> Modified IEC 61217.
    num_subbeams : int
        Number of subbeams, if autosequenced.
    is_truncated : int
        Whether log was truncated due to space limitations; 0 -> not truncated, 1 -> truncated
    num_snapshots : int
        Number of snapshots, cycles, heartbeats, or whatever you'd prefer to call them.
    mlc_model : int
        The MLC model; 2 -> NDS 120 (e.g. Millennium), 3 -> NDS 120 HD (e.g. Millennium 120 HD)
    """
    def _read(self):
        """Read the header section of a tlog."""
        self.header = self._decode_binary(self._log_content, str, 16)  # for version 1.5 will be "VOSTL"
        self.version = float(self._decode_binary(self._log_content, str, 16))  # in the format of 2.x or 3.x
        self.header_size = self._decode_binary(self._log_content, int)  # fixed at 1024 in 1.5 specs
        self.sampling_interval = self._decode_binary(self._log_content, int)
        self.num_axes = self._decode_binary(self._log_content, int)
        self.axis_enum = self._decode_binary(self._log_content, int, self.num_axes)
        self.samples_per_axis = self._decode_binary(self._log_content, int, self.num_axes)
        self.num_mlc_leaves = self.samples_per_axis[-1] - 2  # subtract 2 (each carriage counts as an "axis" and must be removed)
        # self._cursor == self.num_axes * 4 # there is a reserved section after samples per axis. this moves it past it.
        self.axis_scale = self._decode_binary(self._log_content, int)
        self.num_subbeams = self._decode_binary(self._log_content, int)
        self.is_truncated = self._decode_binary(self._log_content, int)
        self.num_snapshots = self._decode_binary(self._log_content, int)
        # the section after MLC model is reserved. Cursor is moved to the end of this reserved section.
        self.mlc_model = self._decode_binary(self._log_content, int, cursor_shift=1024 - (64 + self.num_axes * 8))

        return self, self._cursor


class DlogHeader(DlogSection):
    """The Header section of a dynalog file.

    Attributes
    ----------
    version : str
        The Dynalog version letter.
    patient_name : str
        Patient information.
    plan_filename : str
        Filename if using standalone. If using Treat =<6.5 will produce PlanUID, Beam Number.
        Not yet implemented for this yet.
    tolerance : int
        Plan tolerance.
    num_mlc_leaves : int
        Number of MLC leaves.
    clinac_scale : int
        Clinac scale; 0 -> Varian scale, 1 -> IEC 60601-2-1 scale
    """
    def _read(self):
        """Read the header section of a dynalog."""
        self.version = str(next(self._log_content)[0])
        self.patient_name = next(self._log_content)
        self.plan_filename = next(self._log_content)
        self.tolerance = int(next(self._log_content)[0])
        self.num_mlc_leaves = int(
            next(self._log_content)[0]) * 2  # the # of leaves in a dynalog is actually the # of *PAIRS*, hence the *2.
        self.clinac_scale = int(next(self._log_content)[0])  # 0->Varian scale, 1->IEC scale
        return self, self._log_content


class DlogAxisData(DlogSection):
    """Axis data for dynalogs.

    Attributes
    ----------
    num_snapshots : int
        Number of snapshots recorded.
    mu : :class:`~pylinac.log_analyzer.Axis`
        Current dose fraction

        .. note:: This *can* be gantry rotation under certain conditions. See Dynalog file specs.

    previous_segment_num : :class:`~pylinac.log_analyzer.Axis`
        Previous segment *number*, starting with zero.
    beam_hold : :class:`~pylinac.log_analyzer.Axis`
        Beam hold state; 0 -> holdoff not asserted (beam on), 1 -> holdoff asserted, 2 -> carriage in transition
    beam_on : :class:`~pylinac.log_analyzer.Axis`
        Beam on state; 1 -> beam is on, 0 -> beam is off
    previous_dose_index : :class:`~pylinac.log_analyzer.Axis`
        Previous segment dose index or previous segment gantry angle.
    next_dose_index : :class:`~pylinac.log_analyzer.Axis`
        Next segment dose index.
    gantry : :class:`~pylinac.log_analyzer.Axis`
        Gantry data in degrees.
    collimator : :class:`~pylinac.log_analyzer.Axis`
        Collimator data in degrees.
    jaws : :class:`~pylinac.log_analyzer.Jaw_Struct`
        Jaw data structure. Data in cm.
    carriage_A : :class:`~pylinac.log_analyzer.Axis`
        Carriage A data. Data in cm.
    carriage_B : :class:`~pylinac.log_analyzer.Axis`
        Carriage B data. Data in cm.
    mlc : :class:`~pylinac.log_analyzer.MLC`
        MLC data structure. Data in cm.
    """
    def __init__(self, log_content, header, bfile):
        super().__init__(log_content)
        self._header = header
        self._bfile = bfile

    def _read(self, exclude_beam_off):
        """Read the dynalog axis data."""
        matrix = np.array([line for line in self._log_content], dtype=float)

        self.num_snapshots = np.size(matrix, 0)

        # assignment of snapshot values
        # There is no "expected" MU in dynalogs, but for fluence calc purposes, it is set to that of the actual
        self.mu = Axis(matrix[:, 0], matrix[:, 0])
        self.previous_segment_num = Axis(matrix[:, 1])
        self.beam_hold = Axis(matrix[:, 2])
        self.beam_on = Axis(matrix[:, 3])
        self.prior_dose_index = Axis(matrix[:, 4])  # currently not used for anything
        self.next_dose_index = Axis(matrix[:, 5])  # ditto
        self.gantry = GantryAxis(matrix[:, 6] / 10)
        self.collimator = HeadAxis(matrix[:, 7] / 10)

        # jaws are in mm; convert to cm by /10
        jaw_y1 = HeadAxis(matrix[:, 8] / 10)
        jaw_y2 = HeadAxis(matrix[:, 9] / 10)
        jaw_x1 = HeadAxis(matrix[:, 10] / 10)
        jaw_x2 = HeadAxis(matrix[:, 11] / 10)
        self.jaws = JawStruct(jaw_x1, jaw_y1, jaw_x2, jaw_y2)

        # carriages are in 100ths of mm; converted to cm.
        self.carriage_A = Axis(matrix[:, 12] / 1000)
        self.carriage_B = Axis(matrix[:, 13] / 1000)

        # remove snapshots where the beam wasn't on if flag passed
        if exclude_beam_off:
            hold_idx = np.where(self.beam_hold.actual == 0)[0]
            beamon_idx = np.where(self.beam_on.actual == 1)[0]
            snapshots = np.intersect1d(hold_idx, beamon_idx)
        else:
            snapshots = list(range(self.num_snapshots))

        self.mlc = MLC(snapshots, self.jaws)
        for leaf in range(1, (self._header.num_mlc_leaves // 2) + 1):
            axis = LeafAxis(expected=matrix[:, (leaf - 1) * 4 + 14], actual=matrix[:, (leaf - 1) * 4 + 15])
            self.mlc.add_leaf_axis(axis, leaf)

        # read in "B"-file to get bank B MLC positions. The file must be in the same folder as the "A"-file.
        # The header info is repeated but we already have that.
        with open(self._bfile) as csvf:
            dlgdata = csv.reader(csvf, delimiter=',')
            matrix = np.array([line for line in dlgdata if int(dlgdata.line_num) >= 7], dtype=float)

        # Add bank B MLC positions to mlc snapshot arrays
        for leaf in range(1, (self._header.num_mlc_leaves // 2) + 1):
            axis = LeafAxis(expected=matrix[:, (leaf - 1) * 4 + 14], actual=matrix[:, (leaf - 1) * 4 + 15])
            self.mlc.add_leaf_axis(axis, leaf+self._header.num_mlc_leaves//2)

        self._scale_dlog_mlc_pos()
        return self

    def _scale_dlog_mlc_pos(self):
        """Convert MLC leaf plane positions to isoplane positions and from 100ths of mm to cm."""
        dynalog_leaf_conversion = 1.96614  # MLC physical plane scaling factor to iso (100cm SAD) plane
        for leaf in range(1, self.mlc.num_leaves + 1):
            self.mlc.leaf_axes[leaf].actual *= dynalog_leaf_conversion / 1000
            self.mlc.leaf_axes[leaf].expected *= dynalog_leaf_conversion / 1000

    @property
    @lru_cache()
    def num_beamholds(self):
        """Return the number of times the beam was held."""
        diffmatrix = np.diff(self.beam_hold.actual)
        num_holds = int(np.sum(diffmatrix == 1))
        return num_holds


class TlogAxisData(TlogSection):
    """One of four data structures outlined in the Trajectory log file specification.
        Holds information on all Axes measured.

    Attributes
    ----------
    collimator : :class:`~pylinac.log_analyzer.Axis`
        Collimator data in degrees.
    gantry : :class:`~pylinac.log_analyzer.Axis`
        Gantry data in degrees.
    jaws : :class:`~pylinac.log_analyzer.Jaw_Struct`
        Jaw data structure. Data in cm.
    couch : :class:`~pylinac.log_analyzer.Couch_Struct`
        Couch data structure. Data in cm.
    mu : :class:`~pylinac.log_analyzer.Axis`
        MU data in MU.
    beam_hold : :class:`~pylinac.log_analyzer.Axis`
        Beam hold state. Beam *pauses* (e.g. Beam Off button pressed) are not recorded in the log.
        Data is automatic hold state.
        0 -> Normal; beam on.
        1 -> Freeze; beam on, dose servo is temporarily turned off.
        2 -> Hold; servo holding beam.
        3 -> Disabled; beam on, dose servo is disable via Service.
    control_point : :class:`~pylinac.log_analyzer.Axis`
        Current control point.
    carriage_A : :class:`~pylinac.log_analyzer.Axis`
        Carriage A data in cm.
    carriage_B : :class:`~pylinac.log_analyzer.Axis`
        Carriage B data in cm.
    mlc : :class:`~pylinac.log_analyzer.MLC`
        MLC data structure; data in cm.
    """
    def __init__(self, log_content, cursor, header):
        super().__init__(log_content, cursor)
        self._header = header

    def _read(self, exclude_beam_off):
        # step size in bytes
        step_size = sum(self._header.samples_per_axis) * 2

        # read in all snapshot data at once, then assign
        snapshot_data = self._decode_binary(self._log_content, float, step_size * self._header.num_snapshots)

        # reshape snapshot data to be a x-by-num_snapshots matrix
        snapshot_data = snapshot_data.reshape(self._header.num_snapshots, -1)

        column = snapshot_col_gen()

        self.collimator = self._get_axis(snapshot_data, next(column), HeadAxis)

        self.gantry = self._get_axis(snapshot_data, next(column), GantryAxis)

        jaw_y1 = self._get_axis(snapshot_data, next(column), HeadAxis)
        jaw_y2 = self._get_axis(snapshot_data, next(column), HeadAxis)
        jaw_x1 = self._get_axis(snapshot_data, next(column), HeadAxis)
        jaw_x2 = self._get_axis(snapshot_data, next(column), HeadAxis)
        self.jaws = JawStruct(jaw_x1, jaw_y1, jaw_x2, jaw_y2)

        vrt = self._get_axis(snapshot_data, next(column), CouchAxis)
        lng = self._get_axis(snapshot_data, next(column), CouchAxis)
        lat = self._get_axis(snapshot_data, next(column), CouchAxis)
        rtn = self._get_axis(snapshot_data, next(column), CouchAxis)
        if is_tlog_v3(self._header.version):
            rol = self._get_axis(snapshot_data, next(column), CouchAxis)
            pit = self._get_axis(snapshot_data, next(column), CouchAxis)
        else:
            rol = None
            pit = None
        self.couch = CouchStruct(vrt, lng, lat, rtn, rol, pit)

        self.mu = self._get_axis(snapshot_data, next(column), BeamAxis)

        self.beam_hold = self._get_axis(snapshot_data, next(column), BeamAxis)

        self.control_point = self._get_axis(snapshot_data, next(column), BeamAxis)

        self.carriage_A = self._get_axis(snapshot_data, next(column), HeadAxis)
        self.carriage_B = self._get_axis(snapshot_data, next(column), HeadAxis)

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
        for leaf_num in range(1, self._header.num_mlc_leaves+1):
            leaf_axis = self._get_axis(snapshot_data, next(column), LeafAxis)
            self.mlc.add_leaf_axis(leaf_axis, leaf_num)

        return self, self._cursor

    @property
    @lru_cache()
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


class CRC(TlogSection):
    """The last data section of a Trajectory log. Is a 2 byte cyclic redundancy check (CRC), specifically
        a CRC-16-CCITT. The seed is OxFFFF."""
    def __init__(self, log_content, cursor):
        super().__init__(log_content, cursor)

    def _read(self):
        # crc = self._decode_binary(self.log_content, str, 2)
        # TODO: figure this out
        pass

def is_tlog_txt_file_around(tlog_filename):
    """Boolean specifying if a Tlog *.txt file is available."""
    try:
        txt_filename = tlog_filename.replace('.bin', '.txt')
    except:
        return False
    if osp.isfile(txt_filename):
        return True
    else:
        return False


def is_log(filename):
    """Boolean specifying if filename is a valid log file."""
    if is_tlog(filename) or is_dlog(filename):
        return True
    else:
        return False

def is_tlog(filename):
    """Boolean specifying if filename is a Trajectory log file."""
    if is_valid_file(filename, raise_error=False):
        unknown_file = open_file(filename)
        header_sample = unknown_file.read(5).decode()
        if 'V' in header_sample:
            return True
        else:
            return False
    else:
        return False

def is_dlog(filename):
    """Boolean specifying if filename is a dynalog file."""
    if is_valid_file(filename, raise_error=False):
        unknown_file = open_file(filename)
        header_sample = unknown_file.read(5).decode()
        if 'B' in header_sample or 'A' in header_sample:
            return True
        else:
            return False
    else:
        return False

def is_tlog_v3(version):
    """Return whether the Tlog version is 3 or not."""
    if version >= 3:
        return True
    else:
        return False

def snapshot_col_gen():
    """Generator function for iterating through the snapshot data columns."""
    col = 0
    while True:
        yield col
        col += 2

def _return_other_dlg(dlg_filename, raise_find_error=True):
    """Return the filename of the corresponding dynalog file.

    For example, if the A*.dlg file was passed in, return the corresponding B*.dlg filename.
    Can find both A- and B-files.

    Parameters
    ----------
    dlg_filename : str
        The absolute file path of the dynalog file.

    Returns
    -------
    str
        The absolute file path to the corresponding dynalog file.
    """
    dlg_dir, dlg_file = osp.split(dlg_filename)
    if dlg_file.startswith('A'):
        file2get = dlg_file.replace("A", "B", 1)
    elif dlg_file.startswith('B'):
        file2get = dlg_file.replace("B", "A", 1)
    else:
        raise ValueError("Unable to decipher log names; ensure dynalogs start with 'A' and 'B'")
    other_filename = osp.join(dlg_dir, file2get)

    if osp.isfile(other_filename):
        return other_filename
    elif raise_find_error:
        raise FileNotFoundError("Complementary dlg file not found; ensure A and B-file are in same directory.")
    else:
        return

def write_single_value(writer, description, value, unit=None):
    writer.writerow([description, str(value), unit])

def write_array(writer, description, value, unit=None):
    # write expected
    for dtype, attr in zip((' Expected', ' Actual'), ('expected', 'actual')):
        if unit is None:
            dtype_desc = description + dtype
        else:
            dtype_desc = description + dtype + ' in units of ' + unit
        arr2write = np.insert(getattr(value, attr).astype(object), 0, dtype_desc)
        writer.writerow(arr2write)


if __name__ == '__main__':
    MachineLog().run_tlog_demo()
