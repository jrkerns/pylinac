"""Module for analyzing images, either film or EPID, for flatness and symmetry."""
import os.path as osp

import numpy as np
import matplotlib.pyplot as plt

from pylinac.core.image import Image
from pylinac.core.decorators import type_accept
from pylinac.core.io import is_valid_file
from pylinac.core.profile import SingleProfile
from pylinac.core.utilities import is_iterable


class _Symmetry:
    POINT_DIFFERENCE = 'point_difference'
    AREA_2 = 'area/2'
    VARIAN = 'varian'
    ELEKTA = 'elekta'
    SIEMENS = 'siemens'
    PDQ_IEC = 'pdq-IEC'


class _Flatness:
    VOM80 = 'VoM80'
    VARIAN = 'varian'
    SIEMENS = 'siemens'
    VOCAX = 'VoCAX'
    ELEKTA = 'elekta'
    IEC = 'IEC'


class BeamImage(Image):
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
        if filepath is not None and is_valid_file(filepath):
            self._load_file(filepath, True)
        else:
            self.array = np.zeros((1,1))

    def load_demo_image(self):
        """Load the demo image."""
        demo_file = osp.join(osp.dirname(__file__), 'demo_files', 'flatsym', 'flatsym_demo.dcm')
        self._load_file(demo_file, to_gray=True)

    def _plot_image(self, ax, plane, position):
        """Plot the image analyzed and a line showing where the profile was taken."""
        # position = self._convert_position(position, plane)
        ax.imshow(self.array, cmap=plt.cm.Greys)
        ax.set_title("Image")
        self._polish_plot(ax)
        if is_crossplane(plane):
            ax.axhline(position[0], color='r')
        elif is_inplane(plane):
            ax.axvline(position[0], color='r')
        elif is_both_planes(plane):
            ax.axhline(position[0], color='r')  # y
            ax.axvline(position[1], color='r')  # x

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
            'elekta' or 'pdq-IEC'. See the _Flatness or _Symmetry class for more method names.

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
        if is_both_planes(plane):
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
        if is_both_planes(plane):
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
            ncols = 3 if is_both_planes(plane) else 2
            fig, (*axs, img_ax) = plt.subplots(ncols=ncols)
            self._plot_image(img_ax, plane, position)
        else:
            if not is_iterable(ax):
                axs = [ax, ]
            else:
                axs = ax

        if is_both_planes(plane):
            planes = ('x', 'in')
        else:
            planes = (plane, )

        for ax, plane, position in zip(axs, planes, position):
            profile = self._get_profile(plane, position)
            symmetry, lt_edge, rt_edge, max_idx = self._get_symmetry(profile, method)
            # plot
            ax.plot(profile.x_values, profile.y_values)

            self._plot_annotation(ax, symmetry, method, profile, 'sym')

            self._plot_title(ax, plane, 'sym')

            # Add CAX and field edge lines
            cax_idx = profile.get_FWXM_center()
            ax.axvline(cax_idx, color='m', linestyle='-.')
            ax.axvline(lt_edge, color='g', linestyle='-.')
            ax.axvline(rt_edge, color='g', linestyle='-.')

            # Show max variation points
            ax.plot(profile.x_values[max_idx], profile.y_values[max_idx], 'rx')
            ax.plot(profile.x_values[rt_edge - (max_idx - lt_edge)], profile.y_values[rt_edge - (max_idx - lt_edge)], 'rx')

            if plot_mirror:
                central_idx = int(round(profile.y_values.size/2))
                offset = cax_idx - central_idx
                mirror_vals = profile.y_values[::-1]
                ax.plot(profile.x_values + 2*offset, mirror_vals)

            self._polish_plot(ax)

        if show:
            plt.tight_layout()
            plt.show()

    def _convert_position(self, position, plane):
        """Convert position from the passed-in value to the pixel location."""
        if is_both_planes(plane):
            if position is 'auto':
                y, x = self._determine_center()
            elif len(position) is 2:
                y = self._parse_position(position[0], plane)
                self._check_position_inbounds(y, plane)
                x = self._parse_position(position[1], plane)
                self._check_position_inbounds(x, plane)
            else:
                raise ValueError("Position argument '{}' not understood".format(position))
            return [y, x]
        elif is_crossplane(plane):
            if position == 'auto':
                loc, _ = self._determine_center()
            else:
                loc = self._parse_position(position, plane)
        elif is_inplane(plane):
            if position is 'auto':
                _, loc = self._determine_center()
            else:
                loc = self._parse_position(position, plane)
        else:
            raise ValueError("Position argument not understood")
        return [loc,]

    def _check_position_inbounds(self, position, plane):
        """Check that the position is within the image index bounds."""
        if is_crossplane(plane):
            if position >= self.array.shape[0]:
                raise IndexError("Y-position {} is out of bounds of image array".format(position))
        elif is_inplane(plane):
            if position >= self.array.shape[1]:
                raise IndexError("X-position {} is out of bounds of image array".format(position))

    def _parse_position(self, position, plane):
        if isinstance(position, float) and 0 < position < 1:
            if is_crossplane(plane):
                return int(round(position * self.array.shape[0]))
            elif is_inplane(plane):
                return int(round(position * self.array.shape[1]))
        elif isinstance(position, int):
            return position
        else:
            raise ValueError("Position argument '{}' not understood.".format(position))

    def _get_symmetry(self, profile, method):
        """Get the actual symmetry of a profile using a given method"""
        if method in (_Symmetry.VARIAN, _Symmetry.POINT_DIFFERENCE):
            values, x_values = profile.get_field_values(field_width=0.8)
            lt_edge, rt_edge = profile.get_field_edges(field_width=0.8)
            cax = profile.get_FWXM_center()
            dcax = profile.y_values[cax]
            max_val = 0
            for lt_pt, rt_pt, x_val in zip(values, values[::-1], x_values):
                val = abs(lt_pt - rt_pt)
                if val > max_val:
                    max_val = val
                    max_idx = x_val
            symmetry = 100 * max_val / dcax
        elif method in (_Symmetry.ELEKTA, _Symmetry.PDQ_IEC):
            values, x_values = profile.get_field_values(field_width=0.8)
            lt_edge, rt_edge = profile.get_field_edges(field_width=0.8)
            max_val = 0
            for lt_pt, rt_pt, x_val in zip(values, values[::-1], x_values):
                val = max(abs(lt_pt / rt_pt), abs(rt_pt / lt_pt))
                if val > max_val:
                    max_val = val
                    max_idx = x_val
            symmetry = 100 * max_val
        # elif method in (_Symmetry.SIEMENS, _Symmetry.AREA_2):
        #     lt_edge, rt_edge = profile.get_field_edges(field_width=1)
        #     lt_area, rt_area = profile.get_field_calculation(field_width=1, calculation='area')
        #     symmetry = (100 * abs(lt_area - rt_area) / (lt_area + rt_area)) / 2
        #     max_idx = 0
        else:
            raise ValueError("Method parameter '{}' invalid".format(method))
        return symmetry, lt_edge, rt_edge, max_idx

    def _get_profile(self, plane, position):
        """Get a profile at the given position along the specified plane."""
        if not self._img_is_loaded:
            raise AttributeError("An image has not yet been loaded")

        if position == 'auto':
            y, x = self._determine_center()
        else:
            if is_crossplane(plane):
                self._check_position_inbounds(position, plane)
                y = position
            elif is_inplane(plane):
                self._check_position_inbounds(position, plane)
                x = position

        if is_crossplane(plane):
            prof = SingleProfile(self.array[y, :])
        elif is_inplane(plane):
            prof = SingleProfile(self.array[:, x])
        else:
            raise ValueError("Plane parameter '{}' not understood".format(plane))
        return prof

    def flatness(self, plane='crossplane', position='auto', method='varian'):
        """Determine the flatness of the image.

        .. seealso:: :meth:`~pylinac.flatsym.BeamImage.plot_flatsym()` for ``plane``, ``position``, and ``method`` parameter info.
        """
        position = self._convert_position(position, plane)
        if is_both_planes(plane):
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
        if method in (_Flatness.VOM80, _Flatness.VARIAN, _Flatness.SIEMENS):
            # Variation over Mean within 80% Field Width
            dmax = profile.get_field_calculation(field_width=0.8, calculation='max')
            dmin = profile.get_field_calculation(field_width=0.8, calculation='min')
            lt_edge, rt_edge = profile.get_field_edges(field_width=0.8)
            flatness = 100 * abs(dmax - dmin) / (dmax + dmin)
        elif method in (_Flatness.ELEKTA, _Flatness.IEC):
            # IEC convention:
            # TODO: add conditionals for small/large fields
            dmax = profile.get_field_calculation(field_width=0.8, calculation='max')
            dmin = profile.get_field_calculation(field_width=0.8, calculation='min')
            lt_edge, rt_edge = profile.get_field_edges(field_width=0.8)
            flatness = 100 * (dmax / dmin)
        else:
            raise ValueError("Method parameter '{}' invalid".format(method))
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
            ncols = 3 if is_both_planes(plane) else 2
            fig, (*axs, img_ax) = plt.subplots(ncols=ncols)
            self._plot_image(img_ax, plane, position)
        else:
            if not is_iterable(ax):
                axs = [ax,]
            else:
                axs = ax

        if is_both_planes(plane):
            planes = ('x', 'in')
        else:
            planes = (plane, )

        for ax, plane, pos in zip(axs, planes, position):
            profile = self._get_profile(plane, pos)
            flatness, dmax, dmin, lt_edge, rt_edge = self._get_flatness(profile, method)

            ax.plot(profile.x_values, profile.y_values)

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

    def _polish_plot(self, ax):
        """Polish aesthetics of the plot."""
        # tighten view around data limits
        ax.autoscale_view(tight=True, scaley=False)
        # remove ticklabels
        ax.get_xaxis().set_ticklabels([])
        ax.get_yaxis().set_ticklabels([])

    def _plot_annotation(self, ax, value, method, profile, flat_or_sym='sym'):
        """Plot the flat/sym text annotation to the axes."""
        near_top = profile.y_values.max() * 0.85
        near_left_edge = profile.x_values.min() + (profile.x_values.min() + profile.x_values.max()) * 0.04
        if flat_or_sym is 'sym':
            t = 'Symmetry'
        else:
            t = 'Flatness'
        ax.text(near_left_edge, near_top,
                '{}: {:2.2f}%'.format(t, value) + '\nUsing ' + method.capitalize() + ' convention',
                rotation=-90)

    def _plot_title(self, ax, plane, flat_or_sym='sym'):
        """Plot the axes title."""
        if is_crossplane(plane):
            prefix = 'Crossplane'
        else:
            prefix = 'Inplane'
        if flat_or_sym is 'sym':
            suffix = ' Symmetry'
        else:
            suffix = ' Flatness'
        ax.set_title(prefix+suffix)

    def _determine_center(self):
        """Automatically find the center of the field based on FWHM."""
        if not self._img_is_loaded:
            raise AttributeError("An image has not yet been loaded")

        self._check_inversion()
        self.ground()

        col_prof = np.median(self.array, 0)
        col_prof = SingleProfile(col_prof)
        row_prof = np.median(self.array, 1)
        row_prof = SingleProfile(row_prof)

        x_cen = col_prof.get_FWXM_center(round=True)
        y_cen = row_prof.get_FWXM_center(round=True)

        return y_cen, x_cen

    def _check_inversion(self):
        """Check the image for inversion (pickets are valleys, not peaks) by sampling the 4 image corners.
        If the average value of the four corners is above the average pixel value, then it is very likely inverted.
        """
        outer_edge = 10
        inner_edge = 30
        TL_corner = self.array[outer_edge:inner_edge, outer_edge:inner_edge]
        BL_corner = self.array[-inner_edge:-outer_edge, -inner_edge:-outer_edge]
        TR_corner = self.array[outer_edge:inner_edge, outer_edge:inner_edge]
        BR_corner = self.array[-inner_edge:-outer_edge, -inner_edge:-outer_edge]
        corner_avg = np.mean((TL_corner, BL_corner, TR_corner, BR_corner))
        if corner_avg > np.mean(self.array.flatten()):
            self.invert()

    @property
    def _img_is_loaded(self):
        """Boolean specifying if an image has been loaded."""
        if self.array.size == 1:
            return False
        else:
            return True


def is_crossplane(plane):
    if plane in 'crossplane' or plane in 'xplane':
        return True
    else:
        return False

def is_inplane(plane):
    if plane in 'inplane':
        return True
    else:
        return False

def is_both_planes(plane):
    if plane is 'both':
        return True
    else:
        return False


if __name__ == '__main__':
    img = BeamImage()
    img.load_demo_image()
    # img.flatness()
    # print(img.symmetry(plane='x'))
    # print(img.flatness(plane='i'))
    # img.plot_flatsym(plane='both')
    # img.plot_symmetry(plane='both')
    img.plot_flatness()
    # img.run_demo(plane='x')
