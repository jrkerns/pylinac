"""Module for analyzing images, either film or EPID, for flatness and symmetry."""
import os.path as osp

import matplotlib.pyplot as plt
import numpy as np

from .core.decorators import type_accept
from .core import image
from .core.io import retrieve_demo_file
from .core.profile import SingleProfile
from .core.utilities import is_iterable, isnumeric
from .settings import get_dicom_cmap


class _Symmetry:
    POINT_DIFFERENCE = 'point difference'
    AREA_2 = 'area/2'
    VARIAN = 'varian'
    ELEKTA = 'elekta'
    SIEMENS = 'siemens'
    PDQ_IEC = 'pdq-iec'


class _Flatness:
    VOM80 = 'vom80'
    VARIAN = 'varian'
    SIEMENS = 'siemens'
    VOCAX = 'vocax'
    ELEKTA = 'elekta'
    IEC = 'iec'


class BeamImage:
    """Class for analyzing flatness & symmetry of a 2D beam image (perpendicular to the beam)."""
    @type_accept(filepath=str)
    def __init__(self, filepath=None):
        """
        Parameters
        ----------
        filepath : None, str
            If None, image must be loaded later.
            If a str, path to the image file.
        """
        if filepath is not None and osp.isfile(filepath):
            self.image = image.load(filepath)
        else:
            self.image = np.zeros((1,1))

    def load_demo_image(self):
        """Load the demo image."""
        demo_file = retrieve_demo_file(url='flatsym_demo.dcm')
        self.image = image.load(demo_file)

    def _plot_image(self, ax, plane, position):
        """Plot the image analyzed and a line showing where the profile was taken."""
        # position = self._convert_position(position, plane)
        ax.imshow(self.image, cmap=get_dicom_cmap())
        ax.set_title("Image")
        self._polish_plot(ax)
        if _is_crossplane(plane):
            ax.axhline(position[0], color='r')
        elif _is_inplane(plane):
            ax.axvline(position[0], color='r')
        elif _is_both_planes(plane):
            ax.axhline(position[0], color='r')  # y
            ax.axvline(position[1], color='r')  # x

        return ax

    def run_demo(self, plane='both', position='auto', method='varian'):
        """Run the BeamImage flatness & symmetry demo.

        .. seealso:: :meth:`~pylinac.flatsym.BeamImage.plot_flatsym()` for ``plane``, ``position``, and ``method`` parameter info.
        """
        self.load_demo_image()
        self.plot_flatsym(plane, position, method)

    def plot_flatsym(self, plane='both', position='auto', method='varian'):
        """Plot *both* the flatness and symmetry.

        Parameters
        ----------
        plane : str
            The plane of interest, either 'crossplane', 'inplane', or 'both'. Shortcut descriptions are also
            allowed, e.g. 'x' == 'cross' ==  'crossplane' and 'i' == 'in' == 'inplane'.
        position : str, 2-element sequence
            If the position is a str, 'auto', the position will be determined automatically,
            at the center of the FWHM.
            If the position is a tuple/list/array it must be 2 elements long, specifying the
            location desired in (y, x) to take the flatness/symmetry over. The elements may be
            either integers, specifying the actual pixel values, or floats <1.0, specifying the
            image fraction (e.g. 0.4 will result in the pixel at 40% distance along an axis).
            Combinations of int/float is also allowed. See Examples section for more.
        method : str
            The method of analysis. There are multiple methods, specified by their "real" name
            as well as the vendor that utilizes the given method. For example, Varian uses the
            variation of the mean in the 80% field width for flatness. The method could be
            specified by using 'varian', or by 'VoM80'. Another example is Elekta, who use
            the Point Difference Quotient-IEC definition for symmetry. Thus, one could use
            'elekta' or 'pdq-IEC'. Method names are case insensitive.

            For flatness, 'Varian', 'VoM80', 'Siemens' all perform the same calculation. 'Elekta' and 'IEC' both
            perform the same calculation.

            For symmetry, 'Varian', 'Point Difference' both perform the same calculation. 'Elekta' and 'PDQ-IEC' both
            perform the same calculation.

            See :ref:`analysis_definitions` for equations.

        Examples
        --------
        >>> bi = BeamImage()
        >>> bi.load_demo_image()

        Defaults:

        >>> bi.plot_flatsym()

        Specify a single plane and different method:

        >>> bi.plot_flatsym(plane='x', method='elekta')

        Specify a custom position:

        >>> bi.plot_flatsym(position=(300, 0.6))  # int/float combos allowed
        """
        if _is_both_planes(plane):
            img_ax = plt.subplot2grid((2,3), (1, -1), rowspan=2)
            f_ax1 = plt.subplot2grid((2,3), (0, 0))
            f_ax2 = plt.subplot2grid((2,3), (0, 1))
            s_ax1 = plt.subplot2grid((2,3), (1, 0))
            s_ax2 = plt.subplot2grid((2,3), (1, 1))
            f_ax = [f_ax1, f_ax2]
            s_ax = [s_ax1, s_ax2]
        else:
            fig, (f_ax, s_ax, img_ax) = plt.subplots(1, 3, figsize=(10, 5))

        self.plot_flatness(plane, position, method, ax=f_ax, show=False)
        self.plot_symmetry(plane, position, method, ax=s_ax, show=False)
        position = self._convert_position(position, plane)
        self._plot_image(img_ax, plane, position)
        plt.tight_layout()
        plt.show()

    def symmetry(self, plane='both', position='auto', method='varian'):
        """Determine and return the symmetry of the image.

        .. seealso:: :meth:`~pylinac.flatsym.BeamImage.plot_flatsym()` for ``plane``, ``position``, and ``method`` parameter info.
        """
        position = self._convert_position(position, plane)
        if _is_both_planes(plane):
            symmetry = [0, 0]
            for idx, (pl, pos) in enumerate(zip(('x', 'in'), position)):
                profile = self._get_profile(pl, pos)
                symmetry[idx], *_ = self._get_symmetry(profile, method)
        else:
            profile = self._get_profile(plane, position[0])
            symmetry, *_ = self._get_symmetry(profile, method)
        return symmetry

    def plot_symmetry(self, plane='both', position='auto', method='varian', plot_mirror=True, show=True, ax=None):
        """Plot the profile, highlighting symmetry.

        Parameters
        ----------
        show_mirror : bool
            If True (default), shows the "mirrored" profile, making visual comparison easier.
        ax : None, matplotlib.Axes, list containing matplotlib.Axes
            If None, the plot will be created on a new figure/axes, otherwise it will be plotted to the passed axes.
        show : bool
            If True (default), the plot will be drawn/shown at the end of the method call. Not showing the
            plot is useful when plotting multiple flat/sym plots.


        .. seealso:: :meth:`~pylinac.flatsym.BeamImage.plot_flatsym()` for ``plane``, ``position``, and ``method`` parameter info.
        """
        position = self._convert_position(position, plane)

        if ax is None:
            ncols = 3 if _is_both_planes(plane) else 2
            fig, (*axs, img_ax) = plt.subplots(ncols=ncols)
            self._plot_image(img_ax, plane, position)
        else:
            if not is_iterable(ax):
                axs = [ax, ]
            else:
                axs = ax

        if _is_both_planes(plane):
            planes = ('x', 'in')
        else:
            planes = (plane, )

        for axis, plane, position in zip(axs, planes, position):
            profile = self._get_profile(plane, position)
            symmetry, lt_edge, rt_edge, max_idx = self._get_symmetry(profile, method)
            # plot
            axis.plot(profile.values)

            self._plot_annotation(axis, symmetry, method, profile, 'sym')

            self._plot_title(axis, plane, 'sym')

            # Add CAX and field edge lines
            cax_idx = profile.fwxm_center()
            axis.axvline(cax_idx, color='m', linestyle='-.')
            axis.axvline(lt_edge, color='g', linestyle='-.')
            axis.axvline(rt_edge, color='g', linestyle='-.')

            # Show max variation points
            axis.plot(profile._indices[max_idx], profile.values[max_idx], 'rx')
            axis.plot(profile._indices[rt_edge - (max_idx - lt_edge)], profile.values[rt_edge - (max_idx - lt_edge)], 'rx')

            if plot_mirror:
                central_idx = int(round(profile.values.size/2))
                offset = cax_idx - central_idx
                mirror_vals = profile.values[::-1]
                axis.plot(profile._indices + 2*offset, mirror_vals)

            self._polish_plot(axis)

        if show:
            plt.tight_layout()
            plt.show()

        return axs

    def _convert_position(self, position, plane):
        """Convert position from the passed-in value to the pixel location."""
        if _is_both_planes(plane):
            if position is 'auto':
                y, x = self._determine_center(plane)
            elif len(position) is 2:
                y = self._parse_position(position[0], 'x')
                x = self._parse_position(position[1], 'y')
            else:
                raise ValueError("Position argument '{0}' must be 'auto' or 2-element sequence to do both planes".format(position))
            return [y, x]
        elif _is_crossplane(plane) or _is_inplane(plane):
            if position == 'auto':
                loc = self._determine_center(plane)
            elif isnumeric(position):
                loc = self._parse_position(position, plane)
            else:
                raise ValueError("Position argument '{0}' must be 'auto' or a number to do single planes".format(position))
            return [loc,]
        else:
            raise ValueError("Plane argument '{0}' not understood".format(plane))

    def _check_position_inbounds(self, position, plane):
        """Check that the position is within the image index bounds."""
        if _is_crossplane(plane):
            if position >= self.image.shape[1]:
                raise IndexError("Y-position {0} is out of bounds of image array".format(position))
        elif _is_inplane(plane):
            if position >= self.image.shape[0]:
                raise IndexError("X-position {0} is out of bounds of image array".format(position))

    def _parse_position(self, position, plane):
        if not _is_crossplane(plane) and not _is_inplane(plane):
            raise ValueError("Plane argument '{0}' must be either inplane or crossplane".format(plane))
        if isinstance(position, (float, np.float64)) and 0 < position < 1:
            if _is_crossplane(plane):
                arr_side = self.image.shape[0]
            elif _is_inplane(plane):
                arr_side = self.image.shape[1]
            pos = int(round(position * arr_side))
        elif isinstance(position, (int, float, np.float64)):
            pos = int(position)
        else:
            raise ValueError("Position argument '{0}' not understood.".format(position))

        self._check_position_inbounds(pos, plane)
        return pos

    def _get_symmetry(self, profile, method):
        """Get the actual symmetry of a profile using a given method"""
        if method.lower() in (_Symmetry.VARIAN, _Symmetry.POINT_DIFFERENCE):
            values = profile.field_values(field_width=0.8)
            lt, rt = profile.field_edges()
            indices = np.arange(lt, rt + 1)
            cax = profile.fwxm_center()
            dcax = profile.values[cax]
            max_val = 0
            for lt_pt, rt_pt, idx in zip(values, values[::-1], indices):
                val = abs(lt_pt - rt_pt)
                if val > max_val:
                    max_val = val
                    max_idx = idx
            symmetry = 100 * max_val / dcax
        elif method.lower() in (_Symmetry.ELEKTA, _Symmetry.PDQ_IEC):
            values = profile.field_values(field_width=0.8)
            indices = np.arange(profile.field_edges()[0], profile.field_edges()[1] + 1)
            max_val = 0
            for lt_pt, rt_pt, idx in zip(values, values[::-1], indices):
                val = max(abs(lt_pt / rt_pt), abs(rt_pt / lt_pt))
                if val > max_val:
                    max_val = val
                    max_idx = idx
            symmetry = 100 * max_val
        # elif method in (_Symmetry.SIEMENS, _Symmetry.AREA_2):
        #     lt_edge, rt_edge = profile.get_field_edges(field_width=1)
        #     lt_area, rt_area = profile.get_field_calculation(field_width=1, calculation='area')
        #     symmetry = (100 * abs(lt_area - rt_area) / (lt_area + rt_area)) / 2
        #     max_idx = 0
        else:
            raise ValueError("Method parameter '{0}' invalid".format(method))

        lt_edge, rt_edge = profile.field_edges(field_width=0.8)

        return symmetry, lt_edge, rt_edge, max_idx

    def _get_profile(self, plane, position):
        """Get a profile at the given position along the specified plane."""
        if not self._img_is_loaded:
            raise AttributeError("An image has not yet been loaded")

        position = self._convert_position(position, plane)

        # if position == 'auto':
        #     y, x = self._determine_center()
        # else:
        #     if _is_crossplane(plane):
        #         self._check_position_inbounds(position, plane)
        #         y = position
        #     elif _is_inplane(plane):
        #         self._check_position_inbounds(position, plane)
        #         x = position

        if _is_crossplane(plane):
            prof = SingleProfile(self.image[position[0], :])
        elif _is_inplane(plane):
            prof = SingleProfile(self.image[:, position[0]])
        return prof

    def flatness(self, plane='crossplane', position='auto', method='varian'):
        """Determine the flatness of the image.

        .. seealso:: :meth:`~pylinac.flatsym.BeamImage.plot_flatsym()` for ``plane``, ``position``, and ``method`` parameter info.
        """
        position = self._convert_position(position, plane)
        if _is_both_planes(plane):
            flatness = [0, 0]
            for idx, (pl, pos) in enumerate(zip(('x', 'in'), position)):
                profile = self._get_profile(pl, pos)
                flatness[idx], *_ = self._get_flatness(profile, method)
        else:
            profile = self._get_profile(plane, position[0])
            flatness, *_ = self._get_flatness(profile, method)
        return flatness

    def _get_flatness(self, profile, method):
        """Get the flatness of the profile according to the specified method."""
        if method.lower() in (_Flatness.VOM80, _Flatness.VARIAN, _Flatness.SIEMENS):
            # Variation over Mean within 80% Field Width
            dmax = profile.field_calculation(field_width=0.8, calculation='max')
            dmin = profile.field_calculation(field_width=0.8, calculation='min')
            flatness = 100 * abs(dmax - dmin) / (dmax + dmin)
        elif method.lower() in (_Flatness.ELEKTA, _Flatness.IEC):
            # IEC convention:
            # TODO: add conditionals for small/large fields
            dmax = profile.field_calculation(field_width=0.8, calculation='max')
            dmin = profile.field_calculation(field_width=0.8, calculation='min')
            flatness = 100 * (dmax / dmin)
        else:
            raise ValueError("Method parameter '{0}' invalid".format(method))

        lt_edge, rt_edge = profile.field_edges(field_width=0.8)

        return flatness, dmax, dmin, lt_edge, rt_edge

    def plot_flatness(self, plane='both', position='auto', method='varian', ax=None, show=True):
        """Plot the profile showing the min and max points.

        Parameters
        ----------
        ax : None, matplotlib.Axes, list of matplotlib.Axes
            If None, the plot will be created on a new figure/axes, otherwise it will be plotted to the passed axes.
        show : bool
            If True (default), the plot will be drawn/shown at the end of the method call. Not showing the
            plot is useful when plotting multiple flat/sym plots.


        .. seealso:: :meth:`~pylinac.flatsym.BeamImage.plot_flatsym()` for ``plane``, ``position``, and ``method`` parameter info.
        """
        position = self._convert_position(position, plane)

        if ax is None:
            ncols = 3 if _is_both_planes(plane) else 2
            fig, (*axs, img_ax) = plt.subplots(ncols=ncols)
            self._plot_image(img_ax, plane, position)
        else:
            if not is_iterable(ax):
                axs = [ax,]
            else:
                axs = ax

        if _is_both_planes(plane):
            planes = ('x', 'in')
        else:
            planes = (plane, )

        for ax, plane, pos in zip(axs, planes, position):
            profile = self._get_profile(plane, pos)
            flatness, dmax, dmin, lt_edge, rt_edge = self._get_flatness(profile, method)

            ax.plot(profile.values)

            ax.axhline(dmax, color='r')
            ax.axhline(dmin, color='r')
            ax.axvline(lt_edge, color='g', linestyle='-.')
            ax.axvline(rt_edge, color='g', linestyle='-.')

            self._plot_annotation(ax, flatness, method, profile, 'flat')

            self._plot_title(ax, plane, 'flat')

            self._polish_plot(ax)

        if show:
            plt.tight_layout()
            plt.show()

        return axs

    def _polish_plot(self, ax):
        """Polish aesthetics of the plot."""
        # tighten view around data limits
        ax.autoscale_view(tight=True, scaley=False)
        # remove ticklabels
        ax.get_xaxis().set_ticklabels([])
        ax.get_yaxis().set_ticklabels([])

    def _plot_annotation(self, ax, value, method, profile, flat_or_sym='sym'):
        """Plot the flat/sym text annotation to the axes."""
        near_top = profile.values.max() * 0.85
        near_left_edge = profile._indices.min() + (profile._indices.min() + profile._indices.max()) * 0.04
        if flat_or_sym is 'sym':
            t = 'Symmetry'
        else:
            t = 'Flatness'
        ax.text(near_left_edge, near_top,
                '{0}: {1:2.2f}%'.format(t, value) + '\nUsing ' + method.capitalize() + ' convention',
                rotation=-90)

        return ax

    def _plot_title(self, ax, plane, flat_or_sym='sym'):
        """Plot the axes title."""
        if _is_crossplane(plane):
            prefix = 'Crossplane'
        else:
            prefix = 'Inplane'
        if flat_or_sym in 'symmetry':
            suffix = ' Symmetry'
        else:
            suffix = ' Flatness'
        ax.set_title(prefix+suffix)

        return ax

    def _determine_center(self, plane):
        """Automatically find the center of the field based on FWHM."""
        if not self._img_is_loaded:
            raise AttributeError("An image has not yet been loaded")

        self.image.check_inversion()
        self.image.ground()

        col_prof = np.median(self.image, 0)
        col_prof = SingleProfile(col_prof)
        row_prof = np.median(self.image, 1)
        row_prof = SingleProfile(row_prof)

        x_cen = col_prof.fwxm_center()
        y_cen = row_prof.fwxm_center()

        if _is_crossplane(plane):
            return y_cen
        elif _is_inplane(plane):
            return x_cen
        elif _is_both_planes(plane):
            return y_cen, x_cen

    @property
    def _img_is_loaded(self):
        """Boolean specifying if an image has been loaded."""
        if self.image.size == 1:
            return False
        else:
            return True


def _is_crossplane(plane):
    if plane in 'crossplane' or plane in 'xplane':
        return True
    else:
        return False

def _is_inplane(plane):
    if plane in 'inplane' or plane in 'yplane':
        return True
    else:
        return False

def _is_both_planes(plane):
    if plane in 'both':
        return True
    else:
        return False
