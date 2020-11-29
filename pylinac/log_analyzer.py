"""
The log analyzer module reads and parses Varian linear accelerator machine logs, both Dynalogs and Trajectory logs. The module also
calculates actual and expected fluences as well as performing gamma evaluations. Data is structured to be easily accessible and
easily plottable.

Unlike most other modules of pylinac, the log analyzer module has no end goal. Data is parsed from the logs, but what is done with that
info, and which info is analyzed is up to the user.

Features:

* **Analyze Dynalogs or Trajectory logs** - Either platform is supported. Tlog versions 2.1 and 3.0 are supported.
* **Read in both .bin and .txt Trajectory log files** - Read in the machine data from both .bin and .txt files to get all the information recorded.
  See the :attr:`~pylinac.log_analyzer.MachineLog.txt` attribute.
* **Save Trajectory log data to CSV** - The Trajectory log binary data format does not allow for easy export of data. Pylinac lets you do
  that so you can use Excel or other software that you use with Dynalogs.
* **Plot or analyze any axis** - Every data axis (e.g. gantry, y1, beam holds, MLC leaves) can be accessed and plotted: the actual, expected, and even the difference.
* **Calculate fluences and gamma** - Besides reading in the MLC positions, pylinac calculates the actual and expected fluence
  as well as the gamma map; DTA and threshold values are adjustable.
* **Anonymize logs** - Both dynalogs and trajectory logs can be "anonymized" by removing the Patient ID from the filename(s)
  and file data.
"""
import collections
import concurrent.futures
import copy
import csv
from functools import lru_cache
import gc
import itertools
from io import BytesIO
import multiprocessing
import os
import os.path as osp
import shutil
from typing import Union, List, Optional, Tuple, Iterable, Sequence

import argue
import matplotlib.pyplot as plt
import numpy as np

from .settings import get_array_cmap
from .core import image
from .core import io
from .core import pdf
from .core.utilities import is_iterable, decode_binary, Structure, open_path

STATIC_IMRT = 'Static IMRT'
DYNAMIC_IMRT = 'Dynamic IMRT'
VMAT = 'VMAT'
IMAGING = 'Imaging'

MLC_FOV_WIDTH_MM = 400
MLC_FOV_HEIGHT_MM = 400
HDMLC_FOV_HEIGHT_MM = 220


class MachineLogs(list):
    """Read in machine logs from a directory. Inherits from list. Batch methods are also provided."""

    def __init__(self, folder: str, recursive: bool=True):
        """
        Parameters
        ----------
        folder : str
            The directory of interest. Will walk through and process any logs, Trajectory or dynalog, it finds.
            Non-log files will be skipped.
        recursive : bool
            Whether to walk through subfolders of passed directory. Only used if ``folder`` is a valid log directory.

        Examples
        --------
        Load a directory upon initialization::

            >>> log_folder = r'C:\path\log\directory'
            >>> logs = MachineLogs(log_folder)

        Batch methods include determining the average gamma and average gamma pass value::

            >>> logs.avg_gamma()
            >>> 0.05 # or whatever it is
            >>> logs.avg_gamma_pct()
            >>> 97.2
        """
        super().__init__()
        self.load_folder(folder, recursive)

    @classmethod
    def from_zip(cls, zfile: str):
        """Instantiate from a ZIP archive.

        Parameters
        ----------
        zfile : str
            Path to the zip archive.
        """
        with io.TemporaryZipDirectory(zfile) as tzd:
            logs = cls(tzd)
        return logs

    @property
    def num_logs(self) -> int:
        """The number of logs currently loaded."""
        return len(self)

    @property
    def num_tlogs(self) -> int:
        """The number of Trajectory logs currently loaded."""
        return sum(isinstance(log, TrajectoryLog) for log in self)

    @property
    def num_dlogs(self) -> int:
        """The number of Trajectory logs currently loaded."""
        return sum(isinstance(log, Dynalog) for log in self)

    def load_folder(self, directory: str, recursive: bool=True):
        """Load log files from a directory and append to existing list.

        Parameters
        ----------
        directory : str, None
            The directory of interest.
            If a string, will walk through and process any logs, Trajectory or dynalog, it finds.
            Non-log files will be skipped.
            If None, files must be loaded later using .load_dir() or .append().
        recursive : bool
            If True (default), will walk through subfolders of passed directory.
            If False, will only search root directory.
        """
        # get log files from directory
        log_files = _get_log_filenames(directory, recursive=recursive)
        if len(log_files) == 0:
            print("No logs found.")
            return

        # actual log loading
        print(f"{len(log_files)} logs found.")
        for idx, file in enumerate(log_files):
            self.append(file)
            print(f"Log loaded: {idx+1} of {len(log_files)}", end='\r')
        print('')

    def _check_empty(self) -> None:
        """Check if any logs have been loaded."""
        if len(self) == 0:
            raise ValueError("No logs have been loaded yet.")

    def report_basic_parameters(self) -> None:
        """Report basic parameters of the logs.

        - Number of logs
        - Average gamma value of all logs
        - Average gamma pass percent of all logs
        """
        print(f"Number of logs: {len(self)}")
        print(f"Average gamma: {self.avg_gamma():3.2f}")
        print(f"Average gamma pass percent: {self.avg_gamma_pct():3.1f}")

    def append(self, obj, recursive: bool=True):
        """Append a log. Overloads list method.

        Parameters
        ----------
        obj : str, Dynalog, TrajectoryLog
            If a string, must point to a log file.
            If a directory, must contain log files.
            If a Dynalog or Trajectory log instance, then simply appends.
        recursive : bool
            Whether to walk through subfolders of passed directory. Only applicable if obj was a directory.
        """
        if isinstance(obj, str):
            if is_log(obj):
                log = load_log(obj)
                super().append(log)
            elif osp.isdir(obj):
                files = io.retrieve_filenames(obj)
                for file in files:
                    self.append(file)
        elif isinstance(obj, (Dynalog, TrajectoryLog)):
            super().append(obj)
        else:
            raise TypeError("Can only append MachineLog or string pointing to a log or log directory.")

    def avg_gamma(self, doseTA: Union[int, float]=1, distTA: Union[int, float]=1, threshold: Union[int, float]=0.1, resolution: Union[int, float]=0.1) -> float:
        """Calculate and return the average gamma of all logs. See :meth:`~pylinac.log_analyzer.GammaFluence.calc_map()`
        for further parameter info."""
        self._check_empty()
        gamma_list = np.zeros(self.num_logs)

        for num, log in enumerate(self):
            log.fluence.gamma.calc_map(doseTA, distTA, threshold, resolution)
            gamma_list[num] = log.fluence.gamma.avg_gamma
            print(f'Calculating gammas: {num+1} of {self.num_logs}', end='\r')
        print('')
        return gamma_list.mean()

    def avg_gamma_pct(self, doseTA: Union[int, float]=1, distTA: Union[int, float]=1, threshold: Union[int, float]=0.1, resolution: Union[int, float]=0.1) -> float:
        """Calculate and return the average gamma pass percent of all logs. See :meth:`~pylinac.log_analyzer.GammaFluence.calc_map()`
        for further parameter info."""
        self._check_empty()
        gamma_list = np.zeros(self.num_logs)

        for num, log in enumerate(self):
            log.fluence.gamma.calc_map(doseTA, distTA, threshold, resolution)
            gamma_list[num] = log.fluence.gamma.pass_prcnt
            print(f"Calculating gamma pass percent: {num+1} of {self.num_logs}", end='\r')
        print('')
        return gamma_list.mean()

    def to_csv(self) -> List[str]:
        """Write trajectory logs to CSV. If there are both dynalogs and trajectory logs,
        only the trajectory logs will be written. File names will be the same as the original log file names.

        Returns
        -------
        list
            A list of all the filenames of the newly created CSV files.
        """
        tlogs_written = False
        files = []
        for log in self:
            if is_tlog(log.filename):
                file = log.to_csv()
                tlogs_written = True
                files.append(file)
        if tlogs_written:
            print('\nAll trajectory logs written to CSV files!')
        else:
            print('\nNo files written. Either no logs are loaded or all logs were dynalogs.')
        return files

    def anonymize(self, inplace: bool=False, suffix: Optional[str]=None):
        """Save anonymized versions of the logs.

        For dynalogs, this replaces the patient ID in the filename(s) and the second line of the log with 'Anonymous<suffix>`.
        This will rename both A* and B* logs if both are present in the same directory.

        For trajectory logs, the patient ID in the filename is replaced with `Anonymous<suffix>` for the .bin file. If the
        associated .txt file is in the same directory it will similarly replace the patient ID in the filename with
        `Anonymous<suffix>`. Additionally, the `Patient ID` row will be replaced with `Patient ID: Anonymous<suffix>`.

        .. note::
            Anonymization is only available for logs loaded locally (i.e. not from a URL or a data stream). To
            anonymize such a log it must be first downloaded or written to a file, then loaded in.

        .. note::
            Anonymization is done to the log *file* itself. The current instance(s) of `MachineLog` will not be anonymized.

        Parameters
        ----------
        inplace : bool
            If False (default), creates an anonymized *copy* of the log(s).
            If True, *renames and replaces* the content of the log file.
        suffix : str, optional
            An optional suffix that is added after `Anonymous` to give specificity to the log.

        Returns
        -------
        list
            A list containing the paths to the newly written files.
        """
        file_list = []
        for log in self:
            files = log.anonymize(inplace=inplace, suffix=suffix)
            file_list += files
        print("\n\nDone anonymizing!")
        return file_list


class Axis:
    """Represents an 'Axis' of a Trajectory log or dynalog file, holding actual and potentially expected and difference values.

    Attributes
    ----------
    Parameters are Attributes
    """
    def __init__(self, actual: np.ndarray, expected: Optional[np.ndarray]=None):
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
    def difference(self) -> np.ndarray:
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

    def plot_actual(self) -> None:
        """Plot the actual positions as a matplotlib figure."""
        self._plot('actual')

    def save_plot_actual(self, filename: str, **kwargs) -> None:
        self._plot('actual', show=False)
        self._save(filename, **kwargs)

    def plot_expected(self) -> None:
        """Plot the expected positions as a matplotlib figure."""
        self._plot('expected')

    def save_plot_expected(self, filename: str, **kwargs) -> None:
        self._plot('expected', show=False)
        self._save(filename, **kwargs)

    def plot_difference(self) -> None:
        """Plot the difference of positions as a matplotlib figure."""
        self._plot('difference')

    def save_plot_difference(self, filename: str, **kwargs) -> None:
        self._plot('difference', show=False)
        self._save(filename, **kwargs)

    @argue.options(param=('actual', 'expected', 'difference'))
    def _plot(self, param: str, show: bool=True):
        """Plot the parameter: actual, expected, or difference."""
        plt.plot(getattr(self, param))
        plt.grid(True)
        plt.autoscale(axis='x', tight=True)
        if show:
            plt.show()

    def _save(self, filename: str, **kwargs):
        """Save the figure to a file, either .png or .html."""
        plt.savefig(filename, **kwargs)


class AxisMovedMixin:
    """Mixin class for Axis."""
    AXIS_MOVE_THRESHOLD: float = 0.003

    @property
    @lru_cache(maxsize=1)
    def moved(self) -> bool:
        """Return whether the axis moved during treatment."""
        return np.std(self.actual) > self.AXIS_MOVE_THRESHOLD


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


class FluenceBase:
    """An abstract base class to be used for the actual and expected fluences.

    Attributes
    ----------
    array : numpy.ndarray
        An array representing the fluence map; will be num_mlc_pairs-x-400/resolution.
        E.g., assuming a Millennium 120 MLC model and a fluence resolution of 0.1mm, the resulting
        matrix will be 60-x-4000.
    resolution : int, float
        The resolution of the fluence calculation; -1 means calculation has not been done yet.
    """
    resolution = -1
    FLUENCE_TYPE = ''  # must be specified by subclass

    def __init__(self, mlc_struct=None, mu_axis: Axis=None, jaw_struct=None):
        """
        Parameters
        ----------
        mlc_struct : MLC_Struct
        mu_axis : BeamAxis
        jaw_struct : Jaw_Struct
        """
        self.array = object
        self._mlc = mlc_struct
        self._mu = mu_axis
        self._jaws = jaw_struct

    def is_map_calced(self, raise_error: bool=False) -> bool:
        """Return a boolean specifying whether the fluence has been calculated."""
        calced = hasattr(self.array, 'size')
        if (not calced) and (raise_error is True):
            raise ValueError("Map has not yet been calculated. Use .calc_map() with desired parameters first.")
        else:
            return calced

    @lru_cache(maxsize=1)
    def calc_map(self, resolution: float=0.1, equal_aspect: bool=False) -> np.ndarray:
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
        equal_aspect : bool
            If True, make the y-direction the same resolution as x. If False, the y-axis will be equal to the number of leaves.

         Returns
         -------
         numpy.ndarray
             A numpy array reconstructing the actual fluence of the log. The size will
             be the number of MLC pairs by 400 / resolution since the MLCs can move anywhere within the
             40cm-wide linac head opening.
         """
        height = MLC_FOV_HEIGHT_MM if not self._mlc.hdmlc else HDMLC_FOV_HEIGHT_MM
        if equal_aspect:
            fluence = np.zeros((int(height/resolution), int(MLC_FOV_WIDTH_MM/resolution)), dtype=np.float)
        else:
            fluence = np.zeros((self._mlc.num_pairs, int(MLC_FOV_WIDTH_MM / resolution)), dtype=np.float)

        self.array = fluence
        self.resolution = resolution

        # check if the beam was actually on at all (e.g. kV setups)
        if len(self._mlc.snapshot_idx) < 1:
            return fluence

        def create_mlc_y_positions(is_hdmlc):
            if not is_hdmlc:
                num_large_leaves = 10
                size_large_leaves = 10 / resolution
                num_small_leaves = 40
                size_small_leaves = 5 / resolution
            else:
                num_large_leaves = 14
                size_large_leaves = 5 / resolution
                num_small_leaves = 32
                size_small_leaves = 2.5 / resolution
            sizes = [size_large_leaves]*num_large_leaves + [size_small_leaves]*num_small_leaves + [size_large_leaves]*num_large_leaves
            return np.cumsum([0,] + sizes).astype(int)

        positions = create_mlc_y_positions(is_hdmlc=self._mlc.hdmlc)

        def yield_leaf_width():
            for idx in range(self._mlc.num_pairs):
                yield (positions[idx], positions[idx+1])

        # calculate the MU delivered in each snapshot. For Tlogs this is absolute; for dynalogs it's normalized.
        mu_matrix = getattr(self._mu, self.FLUENCE_TYPE)
        # if little to no MU was delivered (e.g. MV/kV setup), return
        if np.max(mu_matrix) < 0.5:
            return fluence
        MU_differential = np.array([mu_matrix[0]] + list(np.diff(mu_matrix)))
        MU_total = mu_matrix[-1]

        # calculate each "line" of fluence (the fluence of an MLC leaf pair, e.g. 1 & 61, 2 & 62, etc),
        # and add each "line" to the total fluence matrix
        leaf_offset = self._mlc.num_pairs
        fluence_line = np.zeros(int(400 / resolution), dtype=np.float32)
        pos_offset = int(np.round(200 / resolution))
        for pair, width in zip(range(1, self._mlc.num_pairs + 1), yield_leaf_width()):
            if not self._mlc.leaf_under_y_jaw(pair):
                fluence_line[:] = 0  # emtpy the line values on each new leaf pair
                right_leaf_data = getattr(self._mlc.leaf_axes[pair], self.FLUENCE_TYPE)
                right_leaf_data = np.round(right_leaf_data * 10 / resolution) + pos_offset
                left_leaf_data = getattr(self._mlc.leaf_axes[pair + leaf_offset], self.FLUENCE_TYPE)
                left_leaf_data = -np.round(left_leaf_data * 10 / resolution) + pos_offset
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
                    fluence_line[int(left_edge):int(right_edge)] = MU_total
                if equal_aspect:
                    fluence[width[0]:width[1], :] = np.tile(fluence_line, [width[1] - width[0], 1])
                else:
                    fluence[pair - 1, :] = fluence_line

        # if it's a dynalog, then normalize it because 25000 is such an arbitrary value
        if MU_total == 25000:
            fluence /= MU_total

        return fluence

    def plot_map(self, show: bool=True) -> None:
        """Plot the fluence; the fluence (pixel map) must have been calculated first."""
        self.is_map_calced(raise_error=True)
        plt.clf()
        plt.imshow(self.array, aspect='auto', cmap=get_array_cmap())
        if show:
            plt.show()

    def save_map(self, filename: str, **kwargs) -> None:
        """Save the fluence map figure to a file."""
        self.plot_map(show=False)
        plt.savefig(filename, **kwargs)


class ActualFluence(FluenceBase):
    """The actual fluence object"""
    FLUENCE_TYPE = 'actual'


class ExpectedFluence(FluenceBase):
    """The expected fluence object."""
    FLUENCE_TYPE = 'expected'


class GammaFluence(FluenceBase):
    """Gamma object, including pixel maps of gamma, binary pass/fail pixel map, and others.

    Attributes
    ----------
    array : numpy.ndarray
        The gamma map. Only available after calling calc_map()
    passfail_array : numpy.ndarray
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
    bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1]

    def __init__(self, actual_fluence: ActualFluence, expected_fluence: ExpectedFluence, mlc_struct):
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
        self.array = object
        self.passfail_array = object
        self._actual_fluence = actual_fluence
        self._expected_fluence = expected_fluence
        self._mlc = mlc_struct

    @lru_cache(maxsize=1)
    def calc_map(self, doseTA: Union[int, float] = 1, distTA: Union[int, float] = 1, threshold: Union[int, float] = 0.1,
                 resolution: Union[int, float] = 0.1, calc_individual_maps: bool = False) -> np.ndarray:
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
        if not self._actual_fluence.is_map_calced() or resolution != self._actual_fluence.resolution:
            self._actual_fluence.calc_map(resolution)
        if not self._expected_fluence.is_map_calced() or resolution != self._expected_fluence.resolution:
            self._expected_fluence.calc_map(resolution)

        actual_img = image.load(self._actual_fluence.array, dpi=25.4 / resolution)
        expected_img = image.load(self._expected_fluence.array, dpi=25.4 / resolution)
        gamma_map = actual_img.gamma(expected_img, doseTA=doseTA, distTA=distTA, threshold=threshold)

        # calculate standard metrics
        self.avg_gamma = np.nanmean(gamma_map)
        if np.isnan(self.avg_gamma):
            self.avg_gamma = 0
        pixels_passing = np.sum(gamma_map[~np.isnan(gamma_map)] < 1)
        all_calcd_pixels = np.sum(gamma_map[~np.isnan(gamma_map)] >= 0)
        self.pass_prcnt = pixels_passing / all_calcd_pixels * 100
        gamma_map = np.nan_to_num(gamma_map)
        self.passfail_array = gamma_map >= 1

        self.distTA = distTA
        self.doseTA = doseTA
        self.threshold = threshold
        self.resolution = resolution

        self.array = gamma_map
        return gamma_map

    def plot_map(self, show: bool=True):
        """Plot the fluence; the fluence (pixel map) must have been calculated first."""
        self.is_map_calced(raise_error=True)
        plt.imshow(self.array, aspect='auto', vmax=1, cmap=get_array_cmap())
        plt.colorbar()
        plt.show()

    def histogram(self, bins: Optional[list]=None) -> Tuple[np.ndarray, np.ndarray]:
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
        self.is_map_calced(raise_error=True)
        if bins is None:
            bins = self.bins
        hist_arr, bin_edges = np.histogram(self.array, bins=bins)
        return hist_arr, bin_edges

    @argue.options(scale=('log', 'linear'))
    def plot_histogram(self, scale: str='log', bins: Optional[list]=None, show: bool=True) -> None:
        """Plot a histogram of the gamma map values.

        Parameters
        ----------
        scale : {'log', 'linear'}
            Scale of the plot y-axis.
        bins : sequence
            The bin edges for the gamma histogram; see numpy.histogram for more info.
        """
        self.is_map_calced(raise_error=True)
        if bins is None:
            bins = self.bins
        plt.clf()
        plt.hist(self.array.flatten(), bins=bins)
        plt.yscale(scale)
        if show:
            plt.show()

    @argue.options(scale=('log', 'linear'))
    def save_histogram(self, filename: str, scale: str='log', bins: Optional[list]=None, **kwargs) -> None:
        """Save the histogram plot to file."""
        self.plot_histogram(scale, bins, show=False)
        plt.savefig(filename, **kwargs)

    def plot_passfail_map(self) -> None:
        """Plot the binary gamma map, only showing whether pixels passed or failed."""
        self.is_map_calced(raise_error=True)
        plt.imshow(self.passfail_array, cmap=get_array_cmap())
        plt.show()


class FluenceStruct:
    """Structure for data and methods having to do with fluences.

    Attributes
    ----------
    actual : :class:`~pylinac.log_analyzer.FluenceBase`
        The actual fluence delivered.
    expected : :class:`~pylinac.log_analyzer.FluenceBase`
        The expected, or planned, fluence.
    gamma : :class:`~pylinac.log_analyzer.GammaFluence`
        The gamma structure regarding the actual and expected fluences.
    """
    def __init__(self, mlc_struct=None, mu_axis: Axis=None, jaw_struct=None):
        self.actual = ActualFluence(mlc_struct, mu_axis, jaw_struct)
        self.expected = ExpectedFluence(mlc_struct, mu_axis, jaw_struct)
        self.gamma = GammaFluence(self.actual, self.expected, mlc_struct)


class MLC:
    """The MLC class holds MLC information and retrieves relevant data about the MLCs and positions."""
    def __init__(self, log_type, snapshot_idx: Optional[np.ndarray]=None, jaw_struct=None, hdmlc: bool=False, subbeams=None):
        """
        Parameters
        ----------
        log_type: Dynalog, TrajectoryLog
        snapshot_idx : array, list
            The snapshots to be considered for RMS and error calculations (can be all snapshots or just when beam was on).
        jaw_struct : Jaw_Struct
        hdmlc : boolean
            If False (default), indicates a regular MLC model (e.g. Millennium 120).
            If True, indicates an HD MLC model (e.g. Millennium 120 HD).

        Attributes
        ----------
        leaf_axes : dict containing :class:`~pylinac.log_analyzer.Axis`
            The dictionary is keyed by the leaf number, with the Axis as the value.

            .. warning:: Leaf numbers are 1-index based to correspond with Varian convention.
        """
        self.leaf_axes: dict = {}
        self.snapshot_idx = snapshot_idx
        self._jaws = jaw_struct
        self.hdmlc = hdmlc
        self.log_type = log_type
        self.subbeams = subbeams

    @classmethod
    def from_dlog(cls, dlog, jaws, snapshot_data: np.ndarray, snapshot_idx: np.ndarray):
        """Construct an MLC structure from a Dynalog"""
        mlc = MLC(Dynalog, snapshot_idx, jaws)
        for leaf in range(1, (dlog.header.num_mlc_leaves // 2) + 1):
            axis = LeafAxis(expected=snapshot_data[(leaf - 1) * 4 + 14], actual=snapshot_data[(leaf - 1) * 4 + 15])
            mlc.add_leaf_axis(axis, leaf)

        # read in "B"-file to get bank B MLC positions. The file must be in the same folder as the "A"-file.
        # The header info is repeated but we already have that.
        with open(dlog.b_logfile) as csvf:
            dlgdata = csv.reader(csvf, delimiter=',')
            snapshot_data = np.array([line for line in dlgdata][dlog.HEADER_LINE_LENGTH:], dtype=float).transpose()

        # Add bank B MLC positions to mlc snapshot arrays
        for leaf in range(1, (dlog.header.num_mlc_leaves // 2) + 1):
            axis = LeafAxis(expected=snapshot_data[(leaf - 1) * 4 + 14], actual=snapshot_data[(leaf - 1) * 4 + 15])
            mlc.add_leaf_axis(axis, leaf_num=leaf+dlog.header.num_mlc_leaves // 2)

        # scale dynalog leaf positions from the physical plane to the isocenter plane and from 100ths of mm to cm.
        dynalog_leaf_conversion = 1.96614  # MLC physical plane scaling factor to iso (100cm SAD) plane
        for leaf in range(1, mlc.num_leaves + 1):
            mlc.leaf_axes[leaf].actual *= dynalog_leaf_conversion / 1000
            mlc.leaf_axes[leaf].expected *= dynalog_leaf_conversion / 1000
        return mlc

    @classmethod
    def from_tlog(cls, tlog, subbeams, jaws, snapshot_data, snapshot_idx, column_iter):
        """Construct an MLC instance from a Trajectory log."""
        mlc = MLC(TrajectoryLog, snapshot_idx, jaws, tlog.is_hdmlc, subbeams=subbeams)
        for leaf_num in range(1, tlog.header.num_mlc_leaves+1):
            leaf_axis = _get_axis(snapshot_data, next(column_iter), LeafAxis)
            mlc.add_leaf_axis(leaf_axis, leaf_num)
        return mlc

    @property
    def num_pairs(self) -> int:
        """Return the number of MLC pairs."""
        return int(self.num_leaves/2)

    @property
    def num_leaves(self) -> int:
        """Return the number of MLC leaves."""
        return len(self.leaf_axes)

    @property
    def num_snapshots(self) -> int:
        """Return the number of snapshots used for MLC RMS & Fluence calculations.

        .. warning::
            This number may not be the same as the number of recorded snapshots in the log
            since the snapshots where the beam was off may not be included. See :meth:`MachineLog.load`
        """
        return len(self.snapshot_idx)

    @property
    def num_moving_leaves(self) -> int:
        """Return the number of leaves that moved."""
        return len(self.moving_leaves)

    @property
    @lru_cache(maxsize=1)
    def moving_leaves(self) -> np.ndarray:
        """Return an array of the leaves that moved during treatment."""
        threshold = 0.01
        indices = ()
        for leaf_num, leafdata in self.leaf_axes.items():
            if type(self) == TrajectoryLog:
                leaf_stdev = np.std(leafdata.actual[self.subbeams[-1]._snapshots])
            else:
                leaf_stdev = np.std(leafdata.actual[self.snapshot_idx])
            if leaf_stdev > threshold:
                indices += (leaf_num,)
        return np.array(indices)

    def add_leaf_axis(self, leaf_axis: LeafAxis, leaf_num: int) -> None:
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

    def leaf_moved(self, leaf_num: int) -> bool:
        """Return whether the given leaf moved during treatment.

        Parameters
        ----------
        leaf_num : int


        .. warning:: Leaf numbers are 1-index based to correspond with Varian convention.
        """
        return leaf_num in self.moving_leaves

    def pair_moved(self, pair_num: int) -> bool:
        """Return whether the given pair moved during treatment.

        If either leaf moved, the pair counts as moving.

        Parameters
        ----------
        pair_num : int


        .. warning:: Pair numbers are 1-index based to correspond with Varian convention.
        """
        a_leaf = pair_num
        b_leaf = pair_num + self.num_pairs
        return self.leaf_moved(a_leaf) or self.leaf_moved(b_leaf)

    @property
    def _all_leaf_indices(self) -> np.ndarray:
        """Return an array enumerated over all the leaves."""
        return np.array(range(1, len(self.leaf_axes) + 1))

    @argue.options(bank=('both', 'A', 'B'))
    def get_RMS_avg(self, bank: str='both', only_moving_leaves: bool=False):
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
        rms = np.mean(rms_array)
        if np.isnan(rms):
            return 0
        else:
            return rms

    @argue.options(bank=('both', 'A', 'B'))
    def get_RMS_max(self, bank: str='both') -> float:
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
        rms = np.max(rms_array)
        if np.isnan(rms):
            return 0
        else:
            return rms

    @argue.options(bank=('both', 'A', 'B'))
    def get_RMS_percentile(self, percentile: Union[int, float]=95, bank: str='both', only_moving_leaves: bool=False):
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

    def get_RMS(self, leaves_or_bank: Union[str, Iterable]) -> np.ndarray:
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

    @argue.options(bank=('both', 'A', 'B'))
    def get_leaves(self, bank: str='both', only_moving_leaves: bool=False) -> list:
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

    @argue.options(bank=('both', 'A', 'B'))
    def get_error_percentile(self, percentile: Union[int, float]=95, bank: str='both', only_moving_leaves: bool=False) -> float:
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

    def create_error_array(self, leaves: Sequence[int], absolute: bool=True) -> np.ndarray:
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

    def create_RMS_array(self, leaves: Sequence[int]) -> np.ndarray:
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
    def _abs_error_all_leaves(self) -> float:
        """Absolute error of all leaves."""
        return np.abs(self._error_array_all_leaves)

    @property
    @lru_cache(maxsize=1)
    def _error_array_all_leaves(self) -> np.ndarray:
        """Error array of all leaves."""
        mlc_error = np.zeros((self.num_leaves, self.num_snapshots))
        # construct numpy array for easy array calculation
        for leaf in range(self.num_leaves):
            mlc_error[leaf, :] = self.leaf_axes[leaf+1].difference[self.snapshot_idx]
        return mlc_error

    @argue.options(dtype=('actual', 'expected'))
    def _snapshot_array(self, dtype: str='actual') -> np.ndarray:
        """Return an array of the snapshot data of all leaves."""
        arr = np.zeros((self.num_leaves, self.num_snapshots))
        # construct numpy array for easy array calculation
        for leaf in range(self.num_leaves):
            arr[leaf, :] = getattr(self.leaf_axes[leaf + 1], dtype)[self.snapshot_idx]
        return arr

    @property
    @lru_cache(maxsize=1)
    def _RMS_array_all_leaves(self) -> np.ndarray:
        """Return the RMS of all leaves."""
        rms_array = np.array([np.sqrt(np.sum(leafdata.difference[self.snapshot_idx] ** 2) / self.num_snapshots) for leafdata in self.leaf_axes.values()])
        return rms_array

    def leaf_under_y_jaw(self, leaf_num: int) -> bool:
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
        return mlc_position < y1_position or mlc_position - thickness > y2_position

    @argue.options(bank_or_leaf=('A', 'B', 'both'), dtype=('actual', 'expected'))
    def get_snapshot_values(self, bank_or_leaf: str='both', dtype: str='actual') -> np.ndarray:
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

    def plot_mlc_error_hist(self, show: bool=True) -> None:
        """Plot an MLC error histogram."""
        plt.hist(self._abs_error_all_leaves.flatten())
        if show:
            plt.show()

    def save_mlc_error_hist(self, filename: str, **kwargs) -> None:
        """Save the MLC error histogram to file."""
        self.plot_mlc_error_hist(show=False)
        plt.savefig(filename, **kwargs)

    def plot_rms_by_leaf(self, show: bool=True) -> None:
        """Plot RMSs by leaf."""
        plt.clf()
        plt.bar(np.arange(len(self.get_RMS('both')))[::-1], self.get_RMS('both'), align='center')
        if show:
            plt.show()

    def save_rms_by_leaf(self, filename: str, **kwargs) -> None:
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
    def __init__(self, x1: HeadAxis, y1: HeadAxis, x2: HeadAxis, y2: HeadAxis):
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


class CouchStruct:
    """Couch Axes data structure.

    Attributes
    ----------
    vert : :class:`~pylinac.log_analyzer.Axis`
    long : :class:`~pylinac.log_analyzer.Axis`
    latl : :class:`~pylinac.log_analyzer.Axis`
    rotn : :class:`~pylinac.log_analyzer.Axis`
    """
    def __init__(self, vertical: CouchAxis, longitudinal: CouchAxis, lateral: CouchAxis, rotational: CouchAxis,
                 pitch: Optional[CouchAxis] = None, roll: Optional[CouchAxis] = None):
        if not all((isinstance(vertical, CouchAxis),
                    isinstance(longitudinal, CouchAxis),
                    isinstance(lateral, CouchAxis),
                    isinstance(rotational, CouchAxis))):
            raise TypeError("Couch structure must be passed Couch Axes.")
        self.vert = vertical
        self.long = longitudinal
        self.latl = lateral
        self.rotn = rotational
        if pitch is not None:
            self.pitch = pitch
            self.roll = roll


class Subbeam:
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
    def __init__(self, file, log_version: float):
        f = file
        self.control_point = decode_binary(f, int)
        self.mu_delivered = decode_binary(f, float)
        self.rad_time = decode_binary(f, float)
        self.sequence_num = decode_binary(f, int)
        # In Tlogs version 3.0 and up, beam names are 512 byte unicode strings, but in <3.0 they are 32 byte unicode strings
        if log_version >= 3:
            chars = 512
        else:
            chars = 32
        self.beam_name = decode_binary(f, str, chars, 32)

    @property
    def gantry_angle(self) -> float:
        """Median gantry angle of the subbeam."""
        return self._get_metadata_axis('gantry')

    @property
    def collimator_angle(self) -> float:
        """Median collimator angle of the subbeam."""
        return self._get_metadata_axis('collimator')

    @property
    def jaw_x1(self) -> float:
        """Median X1 position of the subbeam."""
        return self._get_metadata_axis('jaws', 'x1')

    @property
    def jaw_x2(self) -> float:
        """Median X2 position of the subbeam."""
        return self._get_metadata_axis('jaws', 'x2')

    @property
    def jaw_y1(self) -> float:
        """Median Y1 position of the subbeam."""
        return self._get_metadata_axis('jaws', 'y1')

    @property
    def jaw_y2(self) -> float:
        """Median Y2 position of the subbeam."""
        return self._get_metadata_axis('jaws', 'y2')

    def _get_metadata_axis(self, attr, subattr=None) -> Axis:
        if subattr is None:
            actual = getattr(self._axis_data, attr).actual[self._snapshots]
            expected = getattr(self._axis_data, attr).expected[self._snapshots]
        else:
            actual = getattr(getattr(self._axis_data, attr), subattr).actual[self._snapshots]
            expected = getattr(getattr(self._axis_data, attr), subattr).expected[self._snapshots]
        return Axis(np.median(actual), np.median(expected))


class SubbeamManager:
    """One of 4 subsections of a trajectory log. Holds a list of Subbeams; only applicable for auto-sequenced beams."""
    def __init__(self, file, header):
        self.subbeams = []
        if header.num_subbeams > 0:
            for _ in range(header.num_subbeams):
                subbeam = Subbeam(file, header.version)
                self.subbeams.append(subbeam)

    def post_hoc_metadata(self, axis_data):
        """From the Axis Data, perform post-hoc analysis and set metadata to the subbeams.
            Gives the subbeams more information, as not much is given directly in the logs."""
        for subbeam_num, subbeam in enumerate(self.subbeams):
            self._set_subbeam_snapshots(axis_data, subbeam_num)
            mlc_subsection = copy.copy(axis_data.mlc)
            mlc_subsection.snapshot_idx = subbeam._snapshots
            subbeam.fluence = FluenceStruct(mlc_subsection, axis_data.mu, axis_data.jaws)
        gc.collect()  # don't know why gc is needed; maybe something to do w/ copy?

    def _set_subbeam_snapshots(self, axis_data, beam_num: int):
        """Get the snapshot indices 1) where the beam was on and 2) between the subbeam control point values."""
        subbeam = self.subbeams[beam_num]
        cp_by_snapshot = axis_data.control_point.actual

        # find upper and lower bounds of the subbeam control points
        cp_lower_bound = subbeam.control_point
        try:
            cp_upper_bound = self.subbeams[beam_num+1].control_point
        except IndexError:
            cp_upper_bound = cp_by_snapshot[-1]

        # extract the snapshots within those control points and drop the beam holds as booleans
        snapshots_within_subbeam = np.logical_and(cp_by_snapshot>=cp_lower_bound, cp_by_snapshot<cp_upper_bound)
        beam_on_snapshots = axis_data.beam_hold.actual == 0
        combined_snaps_as_bool = np.logical_and(beam_on_snapshots, snapshots_within_subbeam)
        # convert boolean array back to snapshot indices
        combined_snapshots = [snapshot_idx for (snapshot_idx, boolean_snapshot) in enumerate(combined_snaps_as_bool) if boolean_snapshot]

        subbeam._snapshots = combined_snapshots
        subbeam._axis_data = axis_data

    def __getitem__(self, item) -> Subbeam:
        return self.subbeams[item]

    def __len__(self):
        return len(self.subbeams)


class LogBase:
    """Base class for the Dynalog and TrajectoryLog classes. Should not be called directly."""
    ANON_LINE = -1

    def __init__(self, filename: str, exclude_beam_off: bool=True):
        if is_log(filename):
            self.filename = filename
            self.exclude_beam_off = exclude_beam_off
        else:
            raise IOError(f"{filename} was not a valid log file")

    @classmethod
    def from_url(cls, url: str, exclude_beam_off: bool=True):
        """Instantiate a log from a URL."""
        filename = io.get_url(url)
        return cls(filename, exclude_beam_off)

    def plot_summary(self, show: bool=True):
        """Plot actual & expected fluence, gamma map, gamma histogram,
            MLC error histogram, and MLC RMS histogram.
        """
        self.fluence.gamma.is_map_calced(raise_error=True)

        # plot the actual fluence
        ax = plt.subplot(2, 3, 1)
        self.plot_subimage('actual', ax, show=False)

        # plot the expected fluence
        ax = plt.subplot(2, 3, 2)
        self.plot_subimage('expected', ax, show=False)

        # plot the gamma map
        ax = plt.subplot(2, 3, 3)
        self.plot_subimage('gamma', ax, show=False)

        # plot the gamma histogram
        ax = plt.subplot(2, 3, 4)
        self.plot_subgraph('gamma', ax, show=False)

        # plot the MLC error histogram
        ax = plt.subplot(2, 3, 5)
        self.plot_subgraph('leaf hist', ax, show=False)

        # plot the leaf RMSs
        ax = plt.subplot(2,3,6)
        self.plot_subgraph('rms', ax, show=False)

        if show:
            plt.show()

    def save_summary(self, filename: str, **kwargs) -> None:
        """Save the summary image to file."""
        self.plot_summary(show=False)
        plt.savefig(filename, **kwargs)
        plt.close()

    @argue.options(img=('actual', 'expected', 'gamma'))
    def plot_subimage(self, img: str, ax: plt.Axes=None, show: bool=True, fontsize: int=10):
        # img: {'actual', 'expected', 'gamma'}
        if ax is None:
            ax = plt.subplot()
        ax.tick_params(axis='both', labelsize=8)
        if img in ('actual', 'expected'):
            title = img.capitalize() + ' Fluence'
            plt.imshow(getattr(self.fluence, img).array.astype(np.float32), aspect='auto', interpolation='none',
                       cmap=get_array_cmap())
        elif img == 'gamma':
            plt.imshow(getattr(self.fluence, img).array.astype(np.float32), aspect='auto', interpolation='none', vmax=1,
                       cmap=get_array_cmap())
            plt.colorbar(ax=ax)
            title = 'Gamma Map'
        ax.autoscale(tight=True)
        ax.set_title(title, fontsize=fontsize)
        if show:
            plt.show()

    @argue.options(img=('actual', 'expected', 'gamma'))
    def save_subimage(self, filename: str, img: str, fontsize: int, **kwargs):
        self.plot_subimage(img, show=False, fontsize=fontsize)
        plt.savefig(filename, **kwargs)
        plt.close()

    def plot_subgraph(self, graph: str, ax: plt.Axes=None, show: bool=True, fontsize: int=10, labelsize: int=8):
        # graph: {'gamma hist', 'leaf hist', 'leaf rms'}
        if ax is None:
            ax = plt.subplot()
        if graph.find('gam') >= 0:
            title = 'Gamma Histogram'
            plt.hist(self.fluence.gamma.array.flatten(), bins=self.fluence.gamma.bins)
            ax.set_yscale('log')
        elif graph.find('hist') >= 0:
            title = 'Leaf Histogram'
            plt.hist(self.axis_data.mlc._abs_error_all_leaves.flatten())
        elif graph.find('rms') >= 0:
            title = 'Leaf RMS (mm)'
            ax.set_xlim([-0.5, self.axis_data.mlc.num_leaves + 0.5])  # bit of padding since bar chart alignment is center
            plt.bar(np.arange(len(self.axis_data.mlc.get_RMS('both')))[::-1], self.axis_data.mlc.get_RMS('both')*10,
                    align='center')
        ax.set_title(title, fontsize=fontsize)
        ax.tick_params(axis='both', labelsize=labelsize)
        ax.grid(True)
        if show:
            plt.show()

    def save_subgraph(self, filename: str, img: str, fontsize: int=10, labelsize: int=8, **kwargs):
        self.plot_subgraph(img, show=False, fontsize=fontsize, labelsize=labelsize)
        plt.savefig(filename, **kwargs)
        plt.close()

    def report_basic_parameters(self, printout: bool=True) -> str:
        """Print the common parameters analyzed when investigating machine logs:

        - Log type
        - Average MLC RMS
        - Maximum MLC RMS
        - 95th percentile MLC error
        - Number of beam holdoffs
        - Gamma pass percentage
        - Average gamma value
        """
        title = f"Results of file: {self.filename}\n"
        if self.treatment_type == IMAGING:
            string = title + "Log is an Imaging field; no statistics can be calculated"
        else:
            avg_rms = f"Average RMS of all leaves: {self.axis_data.mlc.get_RMS_avg(only_moving_leaves=False)*10:3.3f} mm\n"
            max_rms = f"Max RMS error of all leaves: {self.axis_data.mlc.get_RMS_max()*10:3.3f} mm\n"
            p95 = f"95th percentile error: {self.axis_data.mlc.get_error_percentile(95, only_moving_leaves=False)*10:3.3f} mm\n"
            num_holdoffs = f"Number of beam holdoffs: {self.num_beamholds:1.0f}\n"
            self.fluence.gamma.calc_map()
            gamma_pass = f"Gamma pass %: {self.fluence.gamma.pass_prcnt:2.2f}\n"
            gamma_avg = f"Gamma average: {self.fluence.gamma.avg_gamma:2.3f}\n"

            string = title + avg_rms + max_rms + p95 + num_holdoffs + gamma_pass + gamma_avg
        if printout:
            print(string)
        return string

    @property
    def treatment_type(self) -> str:
        """The treatment type of the log. Possible options:

        Returns
        -------
        str
            One of the following:
            * VMAT
            * Dynamic IMRT
            * Static IMRT
            * Imaging
        """
        if isinstance(self, TrajectoryLog):  # trajectory log
            gantry_std = max(subbeam.gantry_angle.actual.std() for subbeam in self.subbeams)
            if np.isnan(gantry_std):
                return IMAGING
        else:
            gantry_std = self.axis_data.gantry.actual.std()
        if gantry_std > 0.5:
            return VMAT
        elif self.axis_data.mu.actual.max() <= 2.1:
            return IMAGING
        elif self.axis_data.mlc.num_moving_leaves == 0 and isinstance(self, TrajectoryLog):
            return STATIC_IMRT
        else:
            return DYNAMIC_IMRT

    @property
    def _underscore_idx(self) -> int:
        base_filename = osp.basename(self.filename)
        under_index = base_filename.find('_')
        if under_index < 0:
            raise NameError(f"Filename `{base_filename}` has no underscore. "
                            "Place an underscore between the patient ID and the rest of the filename and try again.")
        return under_index

    def anonymize(self, inplace: bool=False, destination: Optional[str]=None, suffix: Optional[str]=None) -> list:
        """Save an anonymized version of the log.

        For dynalogs, this replaces the patient ID in the filename(s) and the second line of the log with 'Anonymous<suffix>`.
        This will rename both A* and B* logs if both are present in the same directory.

        For trajectory logs, the patient ID in the filename is replaced with `Anonymous<suffix>` for the .bin file. If the
        associated .txt file is in the same directory it will similarly replace the patient ID in the filename with
        `Anonymous<suffix>`. Additionally, the `Patient ID` row will be replaced with `Patient ID: Anonymous<suffix>`.

        .. note::
            Anonymization is only available for logs loaded locally (i.e. not from a URL or a data stream). To
            anonymize such a log it must be first downloaded or written to a file, then loaded in.

        .. note::
            Anonymization is done to the log *file* itself. The current instance of `MachineLog` will not be anonymized.

        Parameters
        ----------
        inplace : bool
            If False (default), creates an anonymized *copy* of the log(s).
            If True, *renames and replaces* the content of the log file.
        destination : str, optional
            A string specifying the directory where the newly anonymized logs should be placed.
            If None, will place the logs in the same directory as the originals.
        suffix : str, optional
            An optional suffix that is added after `Anonymous` to give specificity to the log.

        Returns
        -------
        list
            A list containing the paths to the newly written files.
        """
        if suffix is None:
            suffix = ''

        # determine destination directory
        if destination is None:
            dest_dir = osp.dirname(self.filename)
        else:
            if not osp.isdir(destination):
                raise NotADirectoryError(f"Specified destination `{destination}` was not a valid directory")
            dest_dir = destination

        # copy or rename the files, depending on `inplace` parameter
        anonymous_filenames = self.anon_file_renames(dest_dir, suffix)
        method = os.rename if inplace else shutil.copy
        for old_file, new_file in anonymous_filenames.items():
            method(old_file, new_file)

        # now the actual anonymization
        for file in self.anon_files(dest_dir, suffix):
            with open(file) as f:
                txtdata = f.readlines()
            txtdata[self.ANON_LINE] = 'Patient ID:\tAnonymous_' + suffix + '\n'
            with open(file, mode='w') as f:
                f.writelines(txtdata)
            print('Anonymized file written to: ', file)

        return anonymous_filenames.values()


class DynalogHeader(Structure):
    """
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
    def __init__(self, dlogdata):
        c = itertools.count()
        super().__init__(version=str(dlogdata[next(c)]),
                         patient_name=dlogdata[next(c)],
                         plan_filename=dlogdata[next(c)],
                         tolerance=int(dlogdata[next(c)][0]),
                         num_mlc_leaves=int(dlogdata[next(c)][0])*2,
                         clinac_scale=int(dlogdata[next(c)][0]))


class DynalogAxisData:
    """
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
    prior_dose_index : :class:`~pylinac.log_analyzer.Axis`
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
    def __init__(self, log, dlogdata):
        """Read the dynalog axis data."""
        snapshot_data = np.array(dlogdata[6:], dtype=np.float64).transpose()

        self.num_snapshots = np.size(snapshot_data, 1)

        c = itertools.count()
        def nx():
            return snapshot_data[next(c)]

        # assignment of snapshot values
        # There is no "expected" MU in dynalogs, but for fluence calc purposes, it is set to that of the actual
        mu = nx()

        # if treatment was vmat then MU is replaced by gantry angle (so stupid). If so, convert to normalized MU by looking at gantry movement.
        def correct_vmat_mu(mu_array):
            if mu_array[-1] == 25000:
                return mu_array
            else:
                abs_diff = list(np.abs(np.diff(mu_array)))
                # this is the cumulative gantry diff, a surrogate for MU. Normalize to 25000 to look like a "normal" dynalog
                cum_gantry_diff = np.array([0,] + list(np.cumsum(abs_diff)/np.sum(abs_diff))) * 25000
                return cum_gantry_diff

        corrected_mu = correct_vmat_mu(mu)

        self.mu = Axis(corrected_mu, corrected_mu)
        self.previous_segment_num = Axis(nx())
        self.beam_hold = Axis(nx())
        self.beam_on = Axis(nx())
        self.prior_dose_index = Axis(nx())  # currently not used for anything
        self.next_dose_index = Axis(nx())  # ditto
        self.gantry = GantryAxis(nx() / 10)
        self.collimator = HeadAxis(nx() / 10)

        # jaws are in mm; convert to cm by /10
        jaw_y1 = HeadAxis(nx() / 10)
        jaw_y2 = HeadAxis(nx() / 10)
        jaw_x1 = HeadAxis(nx() / 10)
        jaw_x2 = HeadAxis(nx() / 10)
        self.jaws = JawStruct(jaw_x1, jaw_y1, jaw_x2, jaw_y2)

        # carriages are in 100ths of mm; converted to cm.
        self.carriage_A = Axis(nx() / 1000)
        self.carriage_B = Axis(nx() / 1000)

        if log.exclude_beam_off:
            hold_idx = np.where(self.beam_hold.actual == 0)[0]
            beamon_idx = np.where(self.beam_on.actual == 1)[0]
            snapshot_idx = np.intersect1d(hold_idx, beamon_idx)
        else:
            snapshot_idx = list(range(self.num_snapshots))

        self.num_snapshots = self.num_snapshots
        self.mlc = MLC.from_dlog(log, self.jaws, snapshot_data, snapshot_idx)


class Dynalog(LogBase):
    """Class for loading, analyzing, and plotting data within a Dynalog file.

    Attributes
    ----------
    header : :class:`~pylinac.log_analyzer.DynalogHeader`
    axis_data : :class:`~pylinac.log_analyzer.DynalogAxisData`
    fluence : :class:`~pylinac.log_analyzer.FluenceStruct`
    """
    ANON_LINE = 1
    HEADER_LINE_LENGTH = 6

    def __init__(self, filename, exclude_beam_off=True):
        super().__init__(filename, exclude_beam_off)
        if not is_dlog(self.filename):
            raise NotADynalogError(f"{self.filename} was not a valid Dynalog file")
        if not self._has_other_file:
            raise DynalogMatchError("Didn't find the matching dynalog file")  # TODO: clean up

        with open(self.a_logfile) as a_log:
            dlgdata = [line for line in csv.reader(a_log, delimiter=',')]
        self.header = DynalogHeader(dlgdata)
        self.axis_data = DynalogAxisData(self, dlgdata)
        self.fluence = FluenceStruct(self.axis_data.mlc, self.axis_data.mu, self.axis_data.jaws)

    def anon_file_renames(self, destination: str, suffix: str) -> dict:
        base_a = osp.basename(self.a_logfile)
        base_b = osp.basename(self.b_logfile)
        anonymous_base_a = base_a[:self._underscore_idx] + '_Anonymous' + suffix + '.dlg'
        anonymous_base_b = base_b[:self._underscore_idx] + '_Anonymous' + suffix + '.dlg'
        anonymous_a = osp.join(destination, anonymous_base_a)
        anonymous_b = osp.join(destination, anonymous_base_b)
        filenames = collections.OrderedDict()
        filenames[self.a_logfile] = anonymous_a
        filenames[self.b_logfile] = anonymous_b
        return filenames

    def anon_files(self, destination: str, suffix: str):
        return self.anon_file_renames(destination, suffix).values()

    def snapshot_idx(self, axis_data) -> list:
        if self.exclude_beam_off:
            hold_idx = np.where(axis_data.beam_hold.actual == 0)[0]
            beamon_idx = np.where(axis_data.beam_on.actual == 1)[0]
            snapshot_idx = np.intersect1d(hold_idx, beamon_idx)
        else:
            snapshot_idx = list(range(self.num_snapshots))
        return snapshot_idx

    @property
    def _has_other_file(self) -> bool:
        """Whether the companion file (A* for B-file or vic versa)."""
        return True if self.identify_other_file(self.filename, raise_find_error=False) is not None else False

    @property
    @lru_cache(maxsize=1)
    def a_logfile(self) -> str:
        """Path of the A* dynalog file."""
        other_dlg_file = self.identify_other_file(self.filename)
        return self.filename if osp.basename(self.filename).startswith('A') else other_dlg_file

    @property
    @lru_cache(maxsize=1)
    def b_logfile(self) -> str:
        """Path of the B* dynalog file."""
        other_dlg_file = self.identify_other_file(self.filename)
        return self.filename if osp.basename(self.filename).startswith('B') else other_dlg_file

    @property
    def num_beamholds(self) -> int:
        """Return the number of times the beam was held."""
        diffmatrix = np.diff(self.axis_data.beam_hold.actual)
        num_holds = int(np.sum(diffmatrix > 0))
        return num_holds

    @classmethod
    def from_demo(cls, exclude_beam_off: bool=True):
        """Load and instantiate from the demo dynalog file included with the package."""
        demo_file = io.retrieve_demo_file(url='AQA.dlg')
        io.retrieve_demo_file(url='BQA.dlg')  # also download "B" dynalog
        return cls(demo_file, exclude_beam_off)

    @staticmethod
    def run_demo():
        """Run the Dynalog demo."""
        dlog = Dynalog.from_demo()
        dlog.report_basic_parameters()
        dlog.plot_summary()

    def publish_pdf(self, filename: str, notes: str=None, metadata: dict=None, open_file: bool=False):
        """Publish (print) a PDF containing the analysis and quantitative results.

        Parameters
        ----------
        filename : (str, file-like object}
            The file to write the results to.
        notes : str, list of strings
            Text; if str, prints single line.
            If list of strings, each list item is printed on its own line.
        open_file : bool
            Whether to open the file using the default program after creation.
        metadata : dict
            Extra data to be passed and shown in the PDF. The key and value will be shown with a colon.
            E.g. passing {'Author': 'James', 'Unit': 'TrueBeam'} would result in text in the PDF like:
            --------------
            Author: James
            Unit: TrueBeam
            --------------
        """
        self.fluence.gamma.calc_map()
        canvas = pdf.PylinacCanvas(filename, page_title="Dynalog Analysis", metadata=metadata)
        canvas.add_text(
                      text=['Dynalog results:',
                            f'Average RMS (mm): {self.axis_data.mlc.get_RMS_avg()*10:2.2f}',
                            f'Max RMS (mm): {self.axis_data.mlc.get_RMS_max()*10:2.2f}',
                            f'95th Percentile error (mm): {self.axis_data.mlc.get_error_percentile(95)*10:2.2f}',
                            f'Number of beam holdoffs: {self.num_beamholds}',
                            f'Gamma pass (%): {self.fluence.gamma.pass_prcnt:2.1f}',
                            f'Gamma average: {self.fluence.gamma.avg_gamma:2.2f}',
                            ],
                        location=(10, 25.5))
        for idx, (x, y, graph) in enumerate(zip((2, 11, 2, 11), (14, 14, 6, 6), ('actual', 'expected', 'gamma', ''))):
            data = BytesIO()
            if idx != 3:
                self.save_subimage(data, graph, fontsize=20)
            else:
                self.save_subgraph(data, 'gamma', fontsize=20, labelsize=12)
            canvas.add_image(data, location=(x, y), dimensions=(9, 9))
        if notes is not None:
            canvas.add_text(location=(1, 5.5), font_size=14, text="Notes:")
            canvas.add_text(location=(1, 5), text=notes)
        canvas.add_new_page()
        for idx, (x, y, graph) in enumerate(zip((5, 5), (13, 2), ('leaf hist', 'leaf rms'))):
            data = BytesIO()
            self.save_subgraph(data, graph, fontsize=20, labelsize=12)
            canvas.add_image(location=(x, y), dimensions=(13, 13), image_data=data)
        canvas.finish()

        if open_file:
            open_path(filename)

    @staticmethod
    def identify_other_file(first_dlg_file: str, raise_find_error: bool=True) -> str:
        """Return the filename of the corresponding dynalog file.

        For example, if the A*.dlg file was passed in, return the corresponding B*.dlg filename.
        Can find both A- and B-files.

        Parameters
        ----------
        first_dlg_file : str
            The absolute file path of the dynalog file.
        raise_find_error : bool
            Whether to raise an error if the file isn't found.

        Returns
        -------
        str
            The absolute file path to the corresponding dynalog file.
        """
        dlg_dir, dlg_file = osp.split(first_dlg_file)
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


class TrajectoryLogAxisData:
    """
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

    def __init__(self, log, file, subbeams):
        # step size in bytes
        step_size = sum(log.header.samples_per_axis) * 2

        # read in all snapshot data at once, then assign
        snapshot_data = decode_binary(file, float, step_size * log.header.num_snapshots)

        # reshape snapshot data to be a x-by-num_snapshots matrix
        snapshot_data = snapshot_data.reshape(log.header.num_snapshots, -1)

        clm_iter = itertools.count(step=2)

        self.collimator = _get_axis(snapshot_data, next(clm_iter), HeadAxis)
        self.gantry = _get_axis(snapshot_data, next(clm_iter), GantryAxis)
        jaw_y1 = _get_axis(snapshot_data, next(clm_iter), HeadAxis)
        jaw_y2 = _get_axis(snapshot_data, next(clm_iter), HeadAxis)
        jaw_x1 = _get_axis(snapshot_data, next(clm_iter), HeadAxis)
        jaw_x2 = _get_axis(snapshot_data, next(clm_iter), HeadAxis)
        self.jaws = JawStruct(jaw_x1, jaw_y1, jaw_x2, jaw_y2)

        vrt = _get_axis(snapshot_data, next(clm_iter), CouchAxis)
        lng = _get_axis(snapshot_data, next(clm_iter), CouchAxis)
        lat = _get_axis(snapshot_data, next(clm_iter), CouchAxis)
        rtn = _get_axis(snapshot_data, next(clm_iter), CouchAxis)
        if log.header.version >= 3:
            pitch = _get_axis(snapshot_data, next(clm_iter), CouchAxis)
            roll = _get_axis(snapshot_data, next(clm_iter), CouchAxis)
        else:
            pitch = None
            roll = None
        self.couch = CouchStruct(vrt, lng, lat, rtn, pitch, roll)

        self.mu = _get_axis(snapshot_data, next(clm_iter), BeamAxis)

        self.beam_hold = _get_axis(snapshot_data, next(clm_iter), BeamAxis)

        self.control_point = _get_axis(snapshot_data, next(clm_iter), BeamAxis)

        self.carriage_A = _get_axis(snapshot_data, next(clm_iter), HeadAxis)
        self.carriage_B = _get_axis(snapshot_data, next(clm_iter), HeadAxis)

        if log.exclude_beam_off:
            snapshot_idx = np.where(self.beam_hold.actual == 0)[0]
        else:
            snapshot_idx = list(range(log.header.num_snapshots))

        self.mlc = MLC.from_tlog(log, subbeams, self.jaws, snapshot_data, snapshot_idx, clm_iter)


class TrajectoryLogHeader:
    """
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

    def __init__(self, file):
        f = file
        self.header = decode_binary(f, str, 16)  # for version 1.5 will be "VOSTL"
        self.version = float(decode_binary(f, str, 16))  # in the format of 2.x or 3.x
        self.header_size = decode_binary(f, int)  # fixed at 1024 in 1.5 specs
        self.sampling_interval = decode_binary(f, int)
        self.num_axes = decode_binary(f, int)
        self.axis_enum = decode_binary(f, int, self.num_axes)
        self.samples_per_axis = decode_binary(f, int, self.num_axes)
        self.num_mlc_leaves = self.samples_per_axis[-1] - 2  # subtract 2 (each carriage counts as an "axis" and must be removed)
        self.axis_scale = decode_binary(f, int)
        self.num_subbeams = decode_binary(f, int)
        self.is_truncated = decode_binary(f, int)
        self.num_snapshots = decode_binary(f, int)
        # the section after MLC model is reserved. Cursor is moved to the end of this reserved section.
        self.mlc_model = decode_binary(f, int, cursor_shift=1024 - (64 + self.num_axes * 8))


class TrajectoryLog(LogBase):
    """A class for loading and analyzing the data of a Trajectory log.

    Attributes
    ----------
    header : `~pylinac.log_analyzer.TrajectoryLogHeader`, which has the following attributes:
    axis_data : `~pylinac.log_analyzer.TrajectoryLogAxisData`
    fluence : `~pylinac.log_analyzer.FluenceStruct`
    subbeams : `~pylinac.log_analyzer.SubbeamManager`
    """
    ANON_LINE = 0

    def __init__(self, filename: str, exclude_beam_off: bool=True):
        super().__init__(filename, exclude_beam_off)

        self._read_txt_file()

        with open(self.filename, mode='rb') as tlogfile:
            self.header = TrajectoryLogHeader(tlogfile)
            self.subbeams = SubbeamManager(tlogfile, self.header)
            self.axis_data = TrajectoryLogAxisData(self, tlogfile, self.subbeams)

        self.subbeams.post_hoc_metadata(self.axis_data)
        if not self.treatment_type == IMAGING:
            self.fluence = FluenceStruct(self.axis_data.mlc, self.axis_data.mu, self.axis_data.jaws)

    @property
    def txt_filename(self) -> str:
        """The name of the associated .txt file for the .bin file. The file may or may not be available."""
        if self.txt is not None:
            return self.filename.replace('.bin', '.txt')

    def anon_file_renames(self, destination: str, suffix: str) -> dict:
        base_filename = osp.basename(self.filename)
        anonymous_base_filename = 'Anonymous' + suffix + base_filename[self._underscore_idx:]
        anonymous_filename = osp.join(destination, anonymous_base_filename)
        filenames = collections.OrderedDict()
        filenames[self.filename] = anonymous_filename
        if self.txt_filename is not None:
            anonymous_txtfilename = anonymous_filename.replace('.bin', '.txt')
            filenames[self.txt_filename] = anonymous_txtfilename
        return filenames

    def anon_files(self, destination: str, suffix: str) -> list:
        renames = self.anon_file_renames(destination, suffix)
        if self.txt is not None:
            return [file for file in renames.values() if '.txt' in file]
        else:
            return []

    def _read_txt_file(self) -> None:
        """Read a Tlog's associated .txt file and put in under the 'txt' attribute."""
        self.txt = None
        if '.bin' in self.filename:  # files downloaded via URL may not have .bin ending
            txt_filename = self.filename.replace('.bin', '.txt')
            if osp.isfile(txt_filename):
                self.txt = {}
                with open(txt_filename) as txtfile:
                    txtdata = txtfile.readlines()
                for line in txtdata:
                    items = line.split(':')
                    if len(items) == 2:
                        self.txt[items[0].strip()] = items[1].strip()

    @classmethod
    def from_demo(cls, exclude_beam_off: bool=True):
        """Load and instantiate from the demo trajetory log file included with the package."""
        demo_file = io.retrieve_demo_file(url='Tlog.bin')
        return cls(demo_file, exclude_beam_off)

    @staticmethod
    def run_demo():
        """Run the Trajectory log demo."""
        tlog = TrajectoryLog.from_demo()
        tlog.report_basic_parameters()
        tlog.plot_summary()

    def to_csv(self, filename: Optional[str]=None) -> str:
        """Write the log to a CSV file.

        Parameters
        ----------
        filename : None, str
            If None (default), the CSV filename will be the same as the filename of the log.
            If a string, the filename will be named so.

        Returns
        -------
        str
            The full filename of the newly created CSV file.
        """
        if filename is None:
            filename = self.filename.replace('bin', 'csv')
        elif not filename.endswith('.csv'):
            filename += '.csv'

        csv_file = open(filename, mode='w')
        writer = csv.writer(csv_file, lineterminator='\n')
        # write header info
        header_titles = ('Tlog File:', 'Signature:', 'Version:', 'Header Size:', 'Sampling Inteval:',
                         'Number of Axes:', 'Axis Enumeration:', 'Samples per Axis:', 'Axis Scale:',
                         'Number of Subbeams:', 'Is Truncated?', 'Number of Snapshots:', 'MLC Model:')
        h = self.header
        header_values = (self.filename, h.header, h.version, h.header_size, h.sampling_interval,
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

        print("CSV file written to: " + filename)
        return filename

    def publish_pdf(self, filename: str, metadata: dict=None, notes: Union[str, list]=None, open_file: bool=False):
        """Publish (print) a PDF containing the analysis and quantitative results.

        Parameters
        ----------
        filename : (str, file-like object}
            The file to write the results to.
        notes : str, list of strings
            Text; if str, prints single line.
            If list of strings, each list item is printed on its own line.
        open_file : bool
            Whether to open the file using the default program after creation.
        metadata : dict
            Extra data to be passed and shown in the PDF. The key and value will be shown with a colon.
            E.g. passing {'Author': 'James', 'Unit': 'TrueBeam'} would result in text in the PDF like:
            --------------
            Author: James
            Unit: TrueBeam
            --------------
        """
        if self.treatment_type == IMAGING:
            raise ValueError("Log is of imaging type (e.g. kV setup) and does not contain relevant gamma/leaf data")
        self.fluence.gamma.calc_map()
        canvas = pdf.PylinacCanvas(filename, page_title="Trajectory Log Analysis", metadata=metadata)
        canvas.add_text(
                      text=['Trajectory Log results:',
                            f'Average RMS (mm): {self.axis_data.mlc.get_RMS_avg()*10:2.2f}',
                            f'Max RMS (mm): {self.axis_data.mlc.get_RMS_max()*10:2.2f}',
                            f'95th Percentile error (mm): {self.axis_data.mlc.get_error_percentile(95)*10:2.2f}',
                            f'Number of beam holdoffs: {self.num_beamholds}',
                            f'Gamma pass (%): {self.fluence.gamma.pass_prcnt:2.1f}',
                            f'Gamma average: {self.fluence.gamma.avg_gamma:2.2f}',
                            ],
                        location=(10, 25.5))
        for idx, (x, y, graph) in enumerate(zip((2, 11, 2, 11), (14, 14, 6, 6), ('actual', 'expected', 'gamma', ''))):
            data = BytesIO()
            if idx != 3:
                self.save_subimage(data, graph, fontsize=20)
            else:
                self.save_subgraph(data, 'gamma', fontsize=20, labelsize=12)
            canvas.add_image(data, location=(x, y), dimensions=(9, 9))
        if notes is not None:
            canvas.add_text(location=(1, 5.5), font_size=14, text="Notes:")
            canvas.add_text(location=(1, 5), text=notes)
        canvas.add_new_page()
        for idx, (x, y, graph) in enumerate(zip((5, 5), (13, 2), ('leaf hist', 'leaf rms'))):
            data = BytesIO()
            self.save_subgraph(data, graph, fontsize=20, labelsize=12)
            canvas.add_image(location=(x, y), dimensions=(13, 13), image_data=data)
        canvas.finish()

        if open_file:
            open_path(filename)

    @property
    def num_beamholds(self) -> int:
        """Return the number of times the beam was held."""
        diffmatrix = np.diff(self.axis_data.beam_hold.actual)
        num_holds = int(np.sum(diffmatrix > 0))
        return num_holds

    @property
    def is_hdmlc(self) -> bool:
        """Whether the machine has an HDMLC or not."""
        return self.header.mlc_model == 3


def anonymize(source: str, inplace: bool=False, destination: bool=None, recursive: bool=True):
    """Quickly anonymize an individual log or directory of logs.
    For directories, threaded execution is performed, making this much faster (10-20x) than loading a ``MachineLogs``
    instance of the folder and using the ``.anonymize()`` method.

    .. note::
        Because ``MachineLog`` instances are not overly memory-efficient, you *may* run into ``MemoryError`` issues.
        To avoid this, try not to anonymize more than ~3000 logs at once.

    Parameters
    ----------
    source : str
        Points to a local log file (e.g. .dlg or .bin file) or to a directory containing log files.
    inplace : bool
        Whether to edit the file itself, or created an anonymized copy and leave the original.
    destination : str, None
        Where the put the anonymized logs. Must point to an existing directory. If None, will place the logs in their original location.
    recursive : bool
        Whether to recursively enter sub-directories below the root source folder.
    """
    def _anonymize(filepath, inplace, destination):
        """Function to anonymize logs; used in the thread executor."""
        if is_tlog(filepath) or (is_dlog(filepath) and osp.basename(filepath).startswith("A")):
            log = load_log(filepath)
            log.anonymize(inplace=inplace, destination=destination)

    # if a single file, just anonymize it
    if osp.isfile(source):
        log = load_log(source)
        log.anonymize(inplace=inplace, destination=destination)
    # if a dir, start a threaded executor and walk the folder.
    elif osp.isdir(source):
        futures = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()*8) as exec:
            for pdir, sdir, files in os.walk(source):
                for file in files:
                    filepath = osp.join(pdir, file)
                    future = exec.submit(_anonymize, filepath, inplace, destination)
                    futures.append(future)
                if not recursive:
                    break
            concurrent.futures.wait(futures)
        print(f"All logs in {source} have been anonymized.")
    else:
        raise NotALogError(f"{source} is not a log file or directory.")


def load_log(file_or_dir: str, exclude_beam_off: bool = True, recursive: bool = True) -> Union[
             TrajectoryLog, Dynalog, MachineLogs]:
    """Load a log file or directory of logs, either dynalogs or Trajectory logs.

    Parameters
    ----------
    file_or_dir : str
        String pointing to a single log file or a directory that contains log files.
    exclude_beam_off : bool
        Whether to include snapshots where the beam was off.
    recursive : bool
        Whether to recursively search a directory. Irrelevant for single log files.

    Returns
    -------
    One of :class:`~pylinac.log_analyzer.Dynalog`, :class:`~pylinac.log_analyzer.TrajectoryLog`,
        :class:`~pylinac.log_analyzer.MachineLogs`.
    """
    if io.is_url(file_or_dir):
        file_or_dir = io.get_url(file_or_dir)
    if osp.isfile(file_or_dir):
        if io.is_zipfile(file_or_dir):
            return MachineLogs.from_zip(file_or_dir)
        if not is_log(file_or_dir):
            raise NotALogError("Not a valid log")
        elif is_tlog(file_or_dir):
            return TrajectoryLog(file_or_dir, exclude_beam_off)
        else:
            return Dynalog(file_or_dir, exclude_beam_off)
    elif osp.isdir(file_or_dir):
        return MachineLogs(file_or_dir, recursive)
    else:
        raise NotALogError(f"'{file_or_dir}' did not point to a valid file, directory, or ZIP archive")

 
def is_log(filename: str) -> bool:
    """Boolean specifying if filename is a valid log file."""
    return is_tlog(filename) or is_dlog(filename)


def is_tlog(filename: str) -> bool:
    """Boolean specifying if filename is a Trajectory log file."""
    return _is_log(filename, ('VOSTL',))


def is_dlog(filename: str) -> bool:
    """Boolean specifying if filename is a Dynalog file."""
    return _is_log(filename, ('B', 'A'))


def _is_log(filename: str, keys: Sequence[str]) -> bool:
    """Internal function that determines whether a file is a log.

    Parameters
    ----------
    filename : str
    keys : iterable of strings
        An iterable of strings that should be in the file. If any key
        is in the file it will return true.
    """
    if osp.isfile(filename):
        try:
            with open(filename, mode='rb') as f:
                header_sample = f.read(5).decode()
            return any(key in header_sample for key in keys)
        except:
            return False
    else:
        return False


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


def _get_log_filenames(directory: str, recursive: bool=True) -> list:
    """Extract the names of real log files from a directory."""
    tlogs = io.retrieve_filenames(directory, is_tlog, recursive=recursive)
    dlogs = io.retrieve_filenames(directory, is_dlog, recursive=recursive)
    # drop double-counted dynalogs (both A & B files; just need one of two)
    idx = 0
    while idx < len(dlogs):
        opp_file = Dynalog.identify_other_file(dlogs[idx], raise_find_error=False)
        if opp_file in dlogs:
            del dlogs[dlogs.index(opp_file)]
        else:
            del dlogs[idx]
            idx -= 1
        idx += 1
    return tlogs + dlogs


def _get_axis(snapshot_data, column, axis_type):
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


class NotALogError(IOError):
    """Machine log error. Indicates that the passed file is not a valid machine log file."""
    pass


class NotADynalogError(IOError):
    """Dynalog error. Indicates that the passed file is not a valid dynalog file."""
    pass


class DynalogMatchError(IOError):
    """Dynalog error. Indicates that the associated file of the dynalog passed in
    (A file if B passed in & vic versa) cannot be found. Ensure associated file is in the same folder
    and has the same name as the passed file, except the first letter."""
    pass
