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
from pylinac.core.utilities import is_iterable


np.seterr(invalid='ignore')  # ignore warnings for invalid numpy operations. Used for np.where() operations on partially-NaN arrays.

"""Named Constants"""
log_types = {'dlog': 'Dynalog', 'tlog': 'Trajectory log'}  # the two current log types
dlg_file_exts = ('.dlg',)  # add more if need be
tlg_file_exts = ('.bin',)  # ditto


class Axis:
    """Represents an 'Axis' of a Trajectory log or dynalog file, holding actual and possibly expected values.

    Attributes
    ----------
    See Parameters
    difference : numpy.ndarray
        An array of the difference between actual and expected values.

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


class Leaf_Axis(Axis):
    """Axis holding leaf information."""
    def __init__(self, actual, expected):
        # force expected argument to be supplied
        super().__init__(actual, expected)


class Gantry_Axis(Axis):
    """Axis holding gantry information."""
    pass


class Head_Axis(Axis):
    """Axis holding head information (e.g. jaw positions, collimator)."""
    pass


class Couch_Axis(Axis):
    """Axis holding couch information."""
    pass


class Beam_Axis(Axis):
    """Axis holding beam information (e.g. MU, beam hold status)."""
    pass


class Fluence(metaclass=ABCMeta):
    """An abstract base class to be used for the actual and expected fluences.

    Attributes
    ----------
    See Parameters
    pixel_map : numpy.ndarray
        An array representing the fluence map; will be num_mlc_pairs-x-400/resolution.
        E.g., assuming a Millenium 120 MLC model and a fluence resolution of 0.1mm, the resulting
        matrix will be 60-x-4000.
    resolution : int, float
        The resolution of the fluence calculation; -1 means calculation has not been done yet.
    """
    pixel_map = np.ndarray
    resolution = -1
    _fluence_type = ''  # must be specified by subclass

    def __init__(self, mlc_struct=None, mu_axis=None, jaw_struct=None):
        """
        Parameters
        ----------
        mlc_struct : MLC_Struct
        mu_axis : Beam_Axis
        jaw_struct : Jaw_Struct
        """
        self.mlc = mlc_struct
        self.mu = mu_axis
        self.jaws = jaw_struct

    @property
    def map_calced(self):
        """Return a boolean specifying whether the fluence has been calculated."""
        if isinstance(self.pixel_map.size, int):
            return True
        else:
            return False

    def _same_conditions(self, resolution):
        """Return whether the conditions passed are the same as prior conditions (for semi-lazy operations)."""
        if self.resolution != resolution:
            return False
        else:
            return True

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
        if self.map_calced and self._same_conditions(resolution):
            return self.pixel_map
        # preallocate arrays for expected and actual fluence of number of leaf pairs-x-4000 (40cm = 4000um, etc)
        fluence = np.zeros((self.mlc.num_pairs, 400 / resolution), dtype=float)

        # calculate the MU delivered in each snapshot. For Tlogs this is absolute; for dynalogs it's normalized.
        mu_matrix = getattr(self.mu, self._fluence_type)
        MU_differential = np.zeros(len(mu_matrix))
        MU_differential[0] = mu_matrix[0]
        MU_differential[1:] = np.diff(mu_matrix)
        MU_differential = MU_differential / mu_matrix[-1]
        MU_cumulative = 1

        # calculate each "line" of fluence (the fluence of an MLC leaf pair, e.g. 1 & 61, 2 & 62, etc),
        # and add each "line" to the total fluence matrix
        fluence_line = np.zeros((400 / resolution))
        leaf_offset = self.mlc.num_pairs
        pos_offset = int(np.round(200 / resolution))
        for pair in range(1, self.mlc.num_pairs + 1):
            if not self.mlc.leaf_under_y_jaw(pair):
                fluence_line[:] = 0  # emtpy the line values on each new leaf pair
                left_leaf_data = getattr(self.mlc.leaf_axes[pair], self._fluence_type)
                left_leaf_data = -np.round(left_leaf_data * 10 / resolution) + pos_offset
                right_leaf_data = getattr(self.mlc.leaf_axes[pair + leaf_offset], self._fluence_type)
                right_leaf_data = np.round(right_leaf_data * 10 / resolution) + pos_offset
                left_jaw_data = np.round((200 / resolution) - (self.jaws.x1.actual * 10 / resolution))
                right_jaw_data = np.round((self.jaws.x2.actual * 10 / resolution) + (200 / resolution))
                if self.mlc.pair_moved(pair):
                    for snapshot in self.mlc.snapshot_idx:
                        lt_mlc_pos = left_leaf_data[snapshot]
                        rt_mlc_pos = right_leaf_data[snapshot]
                        lt_jaw_pos = left_jaw_data[snapshot]
                        rt_jaw_pos = right_jaw_data[snapshot]
                        left_edge = max(lt_mlc_pos, lt_jaw_pos)
                        right_edge = min(rt_mlc_pos, rt_jaw_pos)
                        fluence_line[left_edge:right_edge] += MU_differential[snapshot]
                else:  # leaf didn't move; no need to calc over every snapshot
                    first_snapshot = self.mlc.snapshot_idx[0]
                    lt_mlc_pos = left_leaf_data[first_snapshot]
                    rt_mlc_pos = right_leaf_data[first_snapshot]
                    lt_jaw_pos = left_jaw_data.min()
                    rt_jaw_pos = right_jaw_data.max()
                    left_edge = max(lt_mlc_pos, lt_jaw_pos)
                    right_edge = min(rt_mlc_pos, rt_jaw_pos)
                    fluence_line[left_edge:right_edge] = MU_cumulative
                fluence[pair - 1, :] = fluence_line

        self.pixel_map = fluence
        self.resolution = resolution
        return fluence

    def plot_map(self):
        """Plot the fluence; the pixel map must have been calculated first."""
        if not self.map_calced:
            raise AttributeError("Map not yet calculated; use calc_map()")
        plt.imshow(self.pixel_map, aspect='auto')
        plt.show()

    def _zeroize_under_jaws(self, fluence):
        """Zero the fluence under the jaws.

        E.g. an MLC pair could be tucked under the jaws, and thus the fluence between tips should
        be excluded from fluence calculations since the jaw blocked it.
        """
        pass


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
    passfail_map = np.ndarray
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
        self.actual_fluence = actual_fluence
        self.expected_fluence = expected_fluence
        self.mlc = mlc_struct

    def _same_conditions(self, doseTA, distTA, threshold, resolution):
        """Determine if the passed conditions are the same as the existing ones.

        See calc_map for parameter info.
        """
        if (self.distTA == distTA and self.doseTA == doseTA and self.threshold == threshold and self.resolution == resolution):
            return True
        else:
            return False

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
            Distance-to-agreement in mm
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
        # if gamma has been calculated under same conditions, return it
        if self.map_calced and self._same_conditions(distTA, doseTA, threshold, resolution):
            return self.pixel_map

        # calc fluences if need be
        if not self.actual_fluence.map_calced or resolution != self.actual_fluence.resolution:
            self.actual_fluence.calc_map(resolution)
        if not self.expected_fluence.map_calced or resolution != self.expected_fluence.resolution:
            self.expected_fluence.calc_map(resolution)

        actual = copy.copy(self.actual_fluence.pixel_map)
        expected = copy.copy(self.expected_fluence.pixel_map)

        # set dose values below threshold to 0 so gamma doesn't calculate over it
        actual[actual < (threshold / 100) * np.max(actual)] = 0
        expected[expected < (threshold / 100) * np.max(expected)] = 0

        # preallocation
        gamma_map = np.zeros(expected.shape)

        # image gradient in x-direction (leaf movement direction) using sobel filter
        img_x = spf.sobel(actual, 1)

        # equation: (measurement - reference) / sqrt ( doseTA^2 + distTA^2 * image_gradient^2 )
        for leaf in range(self.mlc.num_pairs):
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

    @lazyproperty
    def histogram(self):
        """Return a histogram array of the gamma map values."""
        if self.map_calced:
            hist_arr, bin_edges = np.histogram(self.pixel_map, bins=self.bins)
            return hist_arr
        else:
            raise AttributeError("Gamma map not yet calculated")

    def plot_histogram(self, scale='log'):
        """Plot a histogram of the gamma map values.

        Parameters
        ----------
        scale : {'log', 'linear'}
            Scale of the plot y-axis.
        """
        ax = plt.hist(self.pixel_map.flatten(), self.bins)
        ax.set_yscale(scale)
        plt.show()

    def plot_passfail_map(self):
        """Plot the binary gamma map, only showing whether pixels passed or failed."""
        plt.imshow(self.passfail_map)
        plt.show()


class Fluence_Struct:
    """Structure for data and methods having to do with fluences."""
    def __init__(self, mlc_struct=None, mu_axis=None, jaw_struct=None):
        self.actual = ActualFluence(mlc_struct, mu_axis, jaw_struct)
        self.expected = ExpectedFluence(mlc_struct, mu_axis, jaw_struct)
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

    @type_accept(leaf_axis=Leaf_Axis, leaf_num=int)
    def add_leaf_axis(self, leaf_axis, leaf_num):
        """Add a leaf axis to the MLC data structure.

        Parameters
        ----------
        leaf_axis : Leaf_Axis
            The leaf axis to be added.
        leaf_num : int
            The leaf number.

            .. warning:: Leaf numbers are 1-index based to correspond with Varian convention.
        """
        self.leaf_axes[leaf_num] = leaf_axis

    @property
    def num_pairs(self):
        """Return the number of MLC pairs."""
        return int(self.num_leaves/2)

    @property
    def num_leaves(self):
        """Return the number of MLC leaves."""
        return len(self.leaf_axes)

    @lazyproperty
    def num_snapshots(self):
        """Return the number of snapshots used for MLC RMS & Fluence calculations.

        .. warning::
            This number may not be the same as the number of recorded snapshots in the log
            since the snapshots where the beam was off may not be included.
        """
        return len(self.snapshot_idx)

    @lazyproperty
    def moving_leaves(self):
        """Return an array of the leaves that moved during treatment."""
        threshold = 0.003
        indices = ()
        for leaf_num, leafdata in self.leaf_axes.items():
            leaf_stdev = np.std(leafdata.actual[self.snapshot_idx])
            if leaf_stdev > threshold:
                indices += (leaf_num,)
        return np.array(indices)

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

    @lazyproperty
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

    def get_RMS_max(self, bank='both', only_moving_leaves=False):
        """Return the overall maximum RMS of given leaves.

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
        return rms_array[leaves]

    @lazyproperty
    def _abs_error_all_leaves(self):
        """Absolute error of all leaves."""
        return np.abs(self._error_array_all_leaves)

    @lazyproperty
    def _error_array_all_leaves(self):
        """Error array of all leaves."""
        mlc_error = np.zeros((self.num_leaves, self.num_snapshots))
        # construct numpy array for easy array calculation
        for leaf in range(self.num_leaves):
            mlc_error[leaf, :] = self.leaf_axes[leaf+1].difference[self.snapshot_idx]
        return mlc_error

    # @lazyproperty
    # def _error_array_moving_leaves(self):
    #     """Return an error array only for leaves that moved."""
    #     moving_leaves = self.moving_leaves
    #     error_array = self._error_array_all_leaves
    #     return error_array[moving_leaves, :]

    @lazyproperty
    def _RMS_array_all_leaves(self):
        """Return the RMS of all leaves."""
        rms_array = np.array([np.sqrt(np.sum(leafdata.difference[self.snapshot_idx] ** 2) / self.num_snapshots) for leafdata in self.leaf_axes.values()])
        return rms_array

    # @lazyproperty
    # def _RMS_array_moving_leaves(self):
    #     """Return an array of RMS of only the moving leaves."""
    #     moving_leaves = self.moving_leaves
    #     rms_array = self._RMS_array_all_leaves
    #     return rms_array[moving_leaves, :]

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

        y2_position = self.jaws.y2.actual.max()*10 + 200
        y1_position = 200 - self.jaws.y1.actual.max()*10
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

    # def leaf_is_under_x_jaw(self, leaf_num):
    #     """Return a boolean specifying if the given leaf is under one of the x jaws."""
        # raise NotImplementedError()
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
                isinstance(x1, Head_Axis),
                isinstance(y1, Head_Axis),
                isinstance(x2, Head_Axis),
                isinstance(y2, Head_Axis))):
            raise TypeError("Head_Axis not passed into Jaw structure")
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2

class Couch_Struct:
    """Couch Axes data structure, holding Couch_Axis instances."""
    def __init__(self, vertical, longitudinal, lateral, rotational):
        if not all((isinstance(vertical, Couch_Axis),
                    isinstance(longitudinal, Couch_Axis),
                    isinstance(lateral, Couch_Axis),
                    isinstance(rotational, Couch_Axis))):
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
    """Data structure for trajectory log "subbeams". Only applicable for auto-sequenced beams."""
    def __init__(self, log_content, cursor, log_version):
        super().__init__(log_content, cursor)
        self.log_version = log_version

    def read(self):
        """Read the tlog subbeam information."""
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


class Subbeam_Constructor:
    """One of 4 subsections of a trajectory log. Holds a list of Subbeams; only applicable for auto-sequenced beams."""
    def __init__(self, log_content, cursor, header):
        self._log_content = log_content
        self._header = header
        self._cursor = cursor

    def read(self):
        """Read all the subbeams of a tlog file.

        Returns
        -------
        list
            Contains instances of Subbeam.
        cursor : int
            Internal; for tracking the cursor position in the file.
        """
        subbeams = []
        if self._header.num_subbeams > 0:
            for subbeam_num in range(self._header.num_subbeams):
                subbeam, self._cursor = Subbeam(self._log_content, self._cursor, self._header.version).read()
                subbeams.append(subbeam)
        return subbeams, self._cursor


class Tlog_Header(TLog_Section):
    """A header object, one of 4 sections of a trajectory log. Holds sampling interval, version, etc."""
    def read(self):
        """Read the header section of a tlog."""
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
    """The Header section of a dynalog file."""
    def read(self):
        """Read the header section of a dynalog."""
        self.version = next(self.log_content)
        self.patient_name = next(self.log_content)
        self.plan_filename = next(self.log_content)
        self.tolerance = int(next(self.log_content)[0])
        self.num_mlc_leaves = int(
            next(self.log_content)[0]) * 2  # the # of leaves in a dynalog is actually the # of *PAIRS*, hence the *2.
        self.clinac_scale = next(self.log_content)  # 0->Varian scale, 1->IEC scale
        return self, self.log_content


class Dlog_Axis_Data(DLog_Section):
    """Axis data for dynalogs."""
    def __init__(self, log_content, header, bfile):
        super().__init__(log_content)
        self.header = header
        self.bfile = bfile

    def read(self, exclude_beam_off=True):
        """Read the dynalog axis data.

        Parameters
        ----------
        exclude_beam_off : bool
            If True (default), excludes the snapshots where the beam was off and beam hold state was asserted.
            If False, no exclusion is performed.

            .. warning:: If all data is included, RMS and error calculations may be affected, since they will include
                snapshots when beam holds were asserted (e.g. moving between steps of a step-&-shoot IMRT plan).
        """
        matrix = np.array([line for line in self.log_content], dtype=float)

        self.num_snapshots = np.size(matrix, 0)

        # assignment of snapshot values
        # There is no "expected" MU in dynalogs, but for fluence calc purposes, it is set to that of the actual
        self.mu = Axis(matrix[:, 0], matrix[:, 0])
        self.DVA_segment = Axis(matrix[:, 1])
        self.beam_hold = Axis(matrix[:, 2])
        self.beam_on = Axis(matrix[:, 3])
        self.prior_dose_idx = Axis(matrix[:, 4])  # currently not used for anything
        self.next_dose_idx = Axis(matrix[:, 5])  # ditto
        self.gantry = Gantry_Axis(matrix[:, 6])
        self.collimator = Head_Axis(matrix[:, 7])

        # jaws are in mm; convert to cm by /10
        jaw_y1 = Head_Axis(matrix[:, 8] / 10)
        jaw_y2 = Head_Axis(matrix[:, 9] / 10)
        jaw_x1 = Head_Axis(matrix[:, 10] / 10)
        jaw_x2 = Head_Axis(matrix[:, 11] / 10)
        self.jaws = Jaw_Struct(jaw_x1, jaw_y1, jaw_x2, jaw_y2)

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
        for leaf in range(1, (self.header.num_mlc_leaves // 2) + 1):
            axis = Leaf_Axis(expected=matrix[:, (leaf-1)*4+14], actual=matrix[:, (leaf-1)*4+15])
            self.mlc.add_leaf_axis(axis, leaf)

        # read in "B"-file to get bank B MLC positions. The file must be in the same folder as the "A"-file.
        # The header info is repeated but we already have that.
        with open(self.bfile) as csvf:
            dlgdata = csv.reader(csvf, delimiter=',')
            matrix = np.array([line for line in dlgdata if int(dlgdata.line_num) >= 7], dtype=float)

        # Add bank B MLC positions to mlc snapshot arrays
        for leaf in range(1, (self.header.num_mlc_leaves // 2) + 1):
            axis = Leaf_Axis(expected=matrix[:, (leaf-1)*4 + 14], actual=matrix[:, (leaf-1)*4 + 15])
            self.mlc.add_leaf_axis(axis, leaf+self.header.num_mlc_leaves//2)

        self._scale_dlog_mlc_pos()
        return self

    def _scale_dlog_mlc_pos(self):
        """Convert MLC leaf plane positions to isoplane positions and from 100ths of mm to cm."""
        dynalog_leaf_conversion = 1.96614  # MLC physical plane scaling factor to iso (100cm SAD) plane
        for leaf in range(1, self.mlc.num_leaves + 1):
            self.mlc.leaf_axes[leaf].actual *= dynalog_leaf_conversion / 1000
            self.mlc.leaf_axes[leaf].expected *= dynalog_leaf_conversion / 1000

    @lazyproperty
    def num_beamholds(self):
        """Return the number of times the beam was held."""
        diffmatrix = np.diff(self.beam_hold.actual)
        num_holds = int(np.sum(diffmatrix == 1))
        return num_holds

class Tlog_Axis_Data(TLog_Section):
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
        self.collimator = self._get_axis(snapshot_data, 0, Head_Axis)

        # gantry
        self.gantry = self._get_axis(snapshot_data, 2, Gantry_Axis)

        # jaws
        jaw_y1 = self._get_axis(snapshot_data, 4, Head_Axis)
        jaw_y2 = self._get_axis(snapshot_data, 6, Head_Axis)
        jaw_x1 = self._get_axis(snapshot_data, 8, Head_Axis)
        jaw_x2 = self._get_axis(snapshot_data, 10, Head_Axis)
        self.jaws = Jaw_Struct(jaw_x1, jaw_y1, jaw_x2, jaw_y2)

        # couch
        vrt = self._get_axis(snapshot_data, 12, Couch_Axis)
        lng = self._get_axis(snapshot_data, 14, Couch_Axis)
        lat = self._get_axis(snapshot_data, 16, Couch_Axis)
        rtn = self._get_axis(snapshot_data, 18, Couch_Axis)
        self.couch = Couch_Struct(vrt, lng, lat, rtn)

        # MU
        self.mu = self._get_axis(snapshot_data, 20, Beam_Axis)

        # beam hold state
        self.beam_hold = self._get_axis(snapshot_data, 22, Beam_Axis)

        # control point
        self.control_point = self._get_axis(snapshot_data, 24, Beam_Axis)

        # carriages
        self.carriage_A = self._get_axis(snapshot_data, 26, Head_Axis)
        self.carriage_B = self._get_axis(snapshot_data, 28, Head_Axis)

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
            leaf_axis = self._get_axis(snapshot_data, 30 + 2 * (leaf_num - 1), Leaf_Axis)
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
    """Reads in and analyzes MLC log files, both dynalog and trajectory logs, from Varian linear accelerators.

    Attributes
    ----------
    header : :class:`Tlog_Header` or :class:`Dlog_header`, depending on the log type.
    axis_data : :class:`Tlog_Axis_Data` or :class:`Dlog_Axis_Data`, depending on the log type.
    subbeams : list
        Only applicable for autosequenced trajectory logs. Will contain instances of :class:`Subbeam`; will be empty if
        autosequencing was not done.
    fluence : :class:`Fluence_Struct`
        Contains actual and expected fluence data, including gamma.
    log_type : str
        The log type loaded; either 'Dynalog' or 'Trajectory log'
    log_is_loaded : bool
        Whether a log has yet been loaded.
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
        self._filename = ''
        self._cursor = 0
        self.fluence = Fluence_Struct()

        # Read file if passed in
        if filename is not '':
            self.load(filename)

    def run_tlog_demo(self):
        """Run the Trajectory log demo."""
        self.load_demo_trajectorylog()
        self.report_basic_parameters()
        self.plot_all()

    def run_dlog_demo(self):
        """Run the dynalog demo."""
        self.load_demo_dynalog()
        self.report_basic_parameters()
        self.plot_all()

    def load_demo_dynalog(self, exclude_beam_off=True):
        """Load the demo dynalog file included with the package."""
        self._filename = osp.join(osp.dirname(__file__), 'demo_files', 'log_reader', 'AQA.dlg')
        # self._filename = osp.join(osp.dirname(__file__), 'demo_files', 'log_reader', 'A20140210084210_Anonymized.dlg')
        self._read_log(exclude_beam_off)

    def load_demo_trajectorylog(self, exclude_beam_off=True):
        """Load the demo trajectory log included with the package."""
        self._filename = osp.join(osp.dirname(__file__), 'demo_files', 'log_reader', 'Tlog.bin')
        self._read_log(exclude_beam_off)

    def load_UI(self):
        """Let user load a log file with a UI dialog box. """
        filename = io.get_filepath_UI()
        if filename: # if user didn't hit cancel...
            self._filename = filename
            self._read_log()

    @type_accept(filename=str)
    def load(self, filename, exclude_beam_off=True):
        """Load the log file directly by passing the path to the file.

        Parameters
        ----------
        filename : str
            The path to the log file.
        """
        if is_valid_file(filename, raise_error=True):
            if self.is_log(filename):
                self._filename = filename
                self._read_log(exclude_beam_off)
            else:
                raise IOError("File passed is not a valid log file")

    def plot_all(self):
        """Plot actual & expected fluence, gamma map, gamma histogram,
            MLC error histogram, and MLC RMS histogram.
        """
        if self.fluence.gamma.map_calced:
            # plot the actual fluence
            ax = plt.subplot(2, 3, 1)
            ax.set_title('Actual Fluence')
            plt.imshow(self.fluence.actual.pixel_map, aspect='auto')

            # plot the expected fluence
            ax = plt.subplot(2, 3, 2)
            ax.set_title("Expected Fluence")
            plt.imshow(self.fluence.expected.pixel_map, aspect='auto')

            # plot the gamma map
            ax = plt.subplot(2, 3, 3)
            ax.set_title("Gamma Map")
            plt.imshow(self.fluence.gamma.pixel_map, aspect='auto')

            # plot the gamma histogram
            ax = plt.subplot(2, 3, 4)
            ax.set_yscale('log')
            ax.set_title("Gamma Histogram")
            plt.hist(self.fluence.gamma.pixel_map.flatten(), self.fluence.gamma.bins)

            # plot the MLC error histogram
            ax = plt.subplot(2, 3, 5)
            p95error = self.axis_data.mlc.get_error_percentile()
            ax.set_title("Leaf Error Histogram (95th Perc: {:2.3f}mm)".format(p95error))
            plt.hist(self.axis_data.mlc._abs_error_all_leaves.flatten())

            # plot the leaf RMSs
            ax = plt.subplot(2,3,6)
            ax.set_title("Leaf RMS Error")
            ax.set_xlim([-0.5, self.axis_data.mlc.num_leaves+0.5])  # bit of padding since bar chart alignment is center
            plt.bar(np.arange(len(self.axis_data.mlc.get_RMS('both')))[::-1], self.axis_data.mlc.get_RMS('both'), align='center')

            plt.show()
        else:
            raise AttributeError("Gamma map has not yet been calculated.")

    def report_basic_parameters(self):
        """Print the common parameters analyzed when investigating machine logs:

        -Log type
        -Average MLC RMS
        -Maximum MLC RMS
        -95th percentile MLC error
        -Mumber of beam holdoffs
        -Gamma pass percentage
        -Average gamma value
        """
        print("MLC log type: " + self.log_type)
        print("Average RMS of all leaves: {:3.3f} mm".format(self.axis_data.mlc.get_RMS_avg(only_moving_leaves=True)))
        print("Max RMS error of all leaves: {:3.3f} mm".format(self.axis_data.mlc.get_RMS_max(only_moving_leaves=True)))
        print("95th percentile error: {:3.3f} mm".format(self.axis_data.mlc.get_error_percentile(95, only_moving_leaves=True)))
        print("Number of beam holdoffs: {:1.0f}".format(self.axis_data.num_beamholds))
        self.fluence.gamma.calc_map()
        print("Gamma pass %: {:2.2f}".format(self.fluence.gamma.pass_prcnt))
        print("Gamma average: {:2.3f}".format(self.fluence.gamma.avg_gamma))

    def is_log(self, filename):
        """Boolean specifying if filename is a valid log file."""
        with open(filename, 'rb') as unknown_file:
            header_sample = unknown_file.read(5).decode()
            if 'B' in header_sample or 'A' in header_sample:
                return True
            elif 'V' in header_sample:
                return True
            else:
                return False

    @lazyproperty
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
        with open(self._filename, 'rb') as unknown_file:
            header_sample = unknown_file.read(5).decode()
            if 'B' in header_sample or 'A' in header_sample:
                log_type = log_types['dlog']  # dynalog
            elif 'V' in header_sample:
                log_type = log_types['tlog']  # trajectory log

        return log_type

    @property
    def log_is_loaded(self):
        """Boolean specifying if a log has been loaded in yet."""
        if self._filename == '':
            return False
        else:
            return True

    def _read_log(self, exclude_beam_off=True):
        """Read in log based on what type of log it is: Trajectory or Dynalog."""
        if not self.log_is_loaded:
            raise AttributeError('Log file has not been specified. Use load_UI() or load()')

        # read log as appropriate to type
        if self.log_type == log_types['tlog']:
            self._read_tlog(exclude_beam_off)
        elif self.log_type == log_types['dlog']:
            self._read_dlog(exclude_beam_off)

    def _read_dlog(self, exclude_beam_off=True):
        """Read in Dynalog files from .dlg files (which are renamed CSV files).
        Formatting follows from the Dynalog File Viewer Reference Guide.
        """
        # if file is B*.dlg, replace with A*.dlg
        other_dlg_file = return_other_dlg(self._filename)

        # create iterator object to read in lines
        with open(self._filename) as csvf:
            dlgdata = csv.reader(csvf, delimiter=',')
            self.header, dlgdata = Dlog_Header(dlgdata).read()
            self.axis_data = Dlog_Axis_Data(dlgdata, self.header, other_dlg_file).read(exclude_beam_off)

        self.fluence = Fluence_Struct(self.axis_data.mlc, self.axis_data.mu, self.axis_data.jaws)

    def _read_tlog(self, exclude_beam_off=True):
        """Read in Trajectory log from binary file according to TB 1.5/2.0 (i.e. Tlog v2.0/3.0) log file specifications.

        Parameters
        ----------
        exclude_beam_off : boolean
            If True (default), snapshots where the beam was not on will be removed.
            If False, all data will be included.

        .. warning::
            Including beam off data may affect fluence and gamma results. E.g. in a step-&-shoot IMRT
            delivery, leaves move between segments while the beam is held. If beam-off data is included,
            this may cause the RMS and fluence errors to rise unnecessarily.
        """
        # read in trajectory log binary data
        fcontent = open(self._filename, 'rb').read()

        # Unpack the content according to respective section and data type (see log specification file).
        self.header, self._cursor = Tlog_Header(fcontent, self._cursor).read()

        self.subbeams, self._cursor = Subbeam_Constructor(fcontent, self._cursor, self.header).read()

        self.axis_data, self._cursor = Tlog_Axis_Data(fcontent, self._cursor, self.header).read(exclude_beam_off)

        # self.crc = CRC(fcontent, self._cursor).read()

        self.fluence = Fluence_Struct(self.axis_data.mlc, self.axis_data.mu, self.axis_data.jaws)


def return_other_dlg(dlg_filename):
    """Checks that the dynalog "B-file" for a given "A-file" exists within the same directory.

    Returns the absolute locations of the B-file.
    """
    dlg_dir, dlg_file = osp.split(dlg_filename)
    if dlg_file.startswith('A'):
        file2get = dlg_file.replace("A", "B", 1)
    elif dlg_file.startswith('B'):
        file2get = dlg_file.replace("B", "A", 1)
    else:
        raise ValueError("Unable to decipher log names; ensure dynalogs start with 'A' and 'B'")
    other_filename = osp.join(dlg_dir, file2get)

    if not osp.isfile(other_filename):
        raise FileNotFoundError("Complementary dlg file not found; ensure A and B-file are in same directory.")
    else:
        return other_filename

# ------------------------
# MLC Log Viewer Example
# ------------------------
if __name__ == '__main__':
    # import cProfile
    # cProfile.run('MachineLog().run_dlog_demo()', sort=1)
    log = MachineLog()
    # log.load_demo_trajectorylog()
    log.run_tlog_demo()
    # log.run_dlog_demo()