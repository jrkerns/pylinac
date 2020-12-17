"""Module for analyzing images or 2D arrays for parameters such as flatness and symmetry.
   Adapted from FlatSym by Alan Chamberlain"""

import io
import os.path as osp
from typing import Union, Optional

import matplotlib.pyplot as plt
import numpy as np

from pylinac.core.utilities import open_path
from .core.exceptions import NotAnalyzed
from .core.io import retrieve_demo_file
from .core import image
from .core.profile import SingleProfile
from .core import pdf
from .settings import get_dicom_cmap

interpolate: bool = True
norm: str = 'max grounded'               # one of 'cax', 'max', 'cax grounded', 'max grounded'

# ----------------------------------------------------------------------------------------------------------------------
# Parameter functions. Each parameter defined in a protocol list must have an associated parameter function that
# returns the value of the parameter.
# ----------------------------------------------------------------------------------------------------------------------


def left_edge_50(profile: SingleProfile, *args) -> float:
    """Return the position of the 50% of max dose value on the left of the profile"""
    left_edge = abs(profile.distance_to_dose(50, norm, interpolate)[0] - profile.profile_center)/profile.dpmm
    return left_edge


def right_edge_50(profile: SingleProfile, *args):
    """Return the position of the 50% of max dose value on the right of the profile"""
    right_edge = abs(profile.distance_to_dose(50, norm, interpolate)[1] - profile.profile_center)/profile.dpmm
    return right_edge


def field_size_50(profile: SingleProfile, *args):
    """Return the field size at 50% of max dose. Not affected by the normalisation mode.
    Included for testing purposes"""
    return profile.fwxm(50)/profile.dpmm


def field_size_edge_50(profile: SingleProfile, *args):
    """Return the field size at 50% of max dose"""
    return right_edge_50(profile) + left_edge_50(profile)


def field_center_fwhm(profile: SingleProfile, *args):
    """Field center as given by the center of the profile FWHM. Not affected by the normalisation mode.
    Included for testing purposes"""
    field_center = (profile.fwxm_center(50, interpolate)[0] - profile.profile_center)/profile.dpmm
    return field_center


def field_center_edge_50(profile: SingleProfile, *args):
    """Calculates the field center from the 50 dose max field edges. May be different from the field_center_fwxm"""
    return (right_edge_50(profile) - left_edge_50(profile))/2


def penumbra_left_80_20(profile: SingleProfile, *args):
    """Return the distance between the 80% and 20% max dose values on the left side of the profile"""
    left_penum = abs(profile.distance_to_dose(80, norm, interpolate)[0]
                     - profile.distance_to_dose(20, norm, interpolate)[0])/profile.dpmm
    return left_penum


def penumbra_right_80_20(profile: SingleProfile, *args):
    """Return the distance between the 80% and 20% max dose values on the right side of the profile"""
    right_penum = abs(profile.distance_to_dose(80, norm, interpolate)[1]
                      - profile.distance_to_dose(20, norm, interpolate)[1])/profile.dpmm
    return right_penum


def flatness_dose_difference(profile: SingleProfile, ifa: float=0.8):
    """The Varian specification for calculating flatness"""
    try:
        dmax = profile.field_calculation(field_width=ifa, calculation='max')
        dmin = profile.field_calculation(field_width=ifa, calculation='min')
    except ValueError:
        raise ValueError("An error was encountered in the flatness calculation. The image is likely inverted. Try inverting the image before analysis with <instance>.image.invert().")
    flatness = 100 * abs(dmax - dmin)/(dmax + dmin)
    return flatness


def flatness_dose_ratio(profile: SingleProfile, ifa: float=0.8):
    """The Elekta specification for calculating flatness"""
    try:
        dmax = profile.field_calculation(field_width=ifa, calculation='max')
        dmin = profile.field_calculation(field_width=ifa, calculation='min')
    except ValueError:
        raise ValueError("An error was encountered in the flatness calculation. The image is likely inverted. Try inverting the image before analysis with <instance>.image.invert().")
    flatness = 100 * (dmax/dmin)
    return flatness


def symmetry_point_difference(profile: SingleProfile, ifa: float=0.8):
    """Calculation of symmetry by way of point difference equidistant from the CAX. Field calculation is
    automatically centred."""
    values = profile.field_values(field_width=ifa)
    _, cax_val = profile.fwxm_center()
    sym_array = []
    for lt_pt, rt_pt in zip(values, values[::-1]):
        val = 100 * abs(lt_pt - rt_pt) / cax_val
        sym_array.append(val)
    symmetry = max(sym_array)
    return symmetry


def symmetry_pdq_iec(profile: SingleProfile, ifa: float = 0.8):
    """Symmetry calculation by way of PDQ IEC. Field calculation is automatically centred"""
    values = profile.field_values(field_width=ifa)
    max_val = 0
    for lt_pt, rt_pt in zip(values, values[::-1]):
        val = max(abs(lt_pt / rt_pt), abs(rt_pt / lt_pt))
        if val > max_val:
            max_val = val
    symmetry = 100 * max_val
    return symmetry


# ----------------------------------------------------------------------------------------------------------------------
# Predefined Protocols - Do not change these. Instead copy a protocol, give it a new name, put it after these protocols
# and add the protocol name to the dictionary PROTOCOLS.
# ----------------------------------------------------------------------------------------------------------------------

ALL = {
    'Left edge (50%)': left_edge_50,
    'Right edge (50%)': right_edge_50,
    'Field size FWHM': field_size_50,
    'Field size edge': field_size_edge_50,
    'Field center edge': field_center_edge_50,
    'Field center FWHM': field_center_fwhm,
    'Penumbra 80-20% left': penumbra_left_80_20,
    'Penumbra 80-20% right': penumbra_right_80_20,
    'Flatness diff': flatness_dose_difference,
    'Flatness ratio': flatness_dose_ratio,
    'Symmetry diff': symmetry_point_difference,
    'Symmetry ratio': symmetry_pdq_iec
}

VARIAN = {
    'Left edge (50%)': left_edge_50,
    'Right edge (50%)': right_edge_50,
    'Field size': field_size_edge_50,
    'Field center': field_center_edge_50,
    'Penumbra 80-20% left': penumbra_left_80_20,
    'Penumbra 80-20% right': penumbra_right_80_20,
    'Flatness': flatness_dose_difference,
    'Symmetry': symmetry_point_difference,
}

ORIGINAL_VARIAN = {
    'Left edge (50%)': left_edge_50,
    'Right edge (50%)': right_edge_50,
    'Field size': field_size_50,
    'Field center': field_center_fwhm,
    'Penumbra 80-20% left': penumbra_left_80_20,
    'Penumbra 80-20% right': penumbra_right_80_20,
    'Flatness': flatness_dose_difference,
    'Symmetry': symmetry_point_difference,
}

ELEKTA = {
    'Left edge (50%)': left_edge_50,
    'Right edge (50%)': right_edge_50,
    'Field size': field_size_50,
    'Field center': field_center_fwhm,
    'Penumbra 80-20% left': penumbra_left_80_20,
    'Penumbra 80-20% right': penumbra_right_80_20,
    'Flatness': flatness_dose_ratio,
    'Symmetry': symmetry_pdq_iec,
}

SIEMENS = {
    'Left edge (50%)': left_edge_50,
    'Right edge (50%)': right_edge_50,
    'Field size': field_size_50,
    'Field center': field_center_fwhm,
    'Penumbra 80-20% left': penumbra_left_80_20,
    'Penumbra 80-20% right': penumbra_right_80_20,
    'Flatness': flatness_dose_difference
}

VOM80 = {
    'Left edge (50%)': left_edge_50,
    'Right edge (50%)': right_edge_50,
    'Field size': field_size_50,
    'Field center': field_center_fwhm,
    'Penumbra 80-20% left': penumbra_left_80_20,
    'Penumbra 80-20% right': penumbra_right_80_20,
    'Flatness': flatness_dose_difference
}

IEC9076 = {
    'Left edge (50%)': left_edge_50,
    'Right edge (50%)': right_edge_50,
    'Field size': field_size_50,
    'Field center': field_center_fwhm,
    'Penumbra 80-20% left': penumbra_left_80_20,
    'Penumbra 80-20% right': penumbra_right_80_20,
    'Symmetry': symmetry_pdq_iec,
}

# ----------------------------------------------------------------------------------------------------------------------
# End of predefined protocols - Do not change these. Instead copy a protocol, give it a new name, put it after these
# protocols and add the protocol name to the dictionary PROTOCOLS.
# ----------------------------------------------------------------------------------------------------------------------

PROTOCOLS = {
    'all': ALL,
    'varian': VARIAN,
    'elekta': ELEKTA,
    'siemens': SIEMENS,
    'vom80': VOM80,
    'iec9076': IEC9076
}


class FieldParams:
    """Class for analyzing the various parameters of a radiation image, most commonly an open image from a linac.

    Attributes
    ----------
    parameters : dict
        Contains the calculated parameters.
    positions : dict
        The position ratio used for analysis for vertical and horizontal.
    """

    def __init__(self, path: str, filter: Optional[int]=None):
        """
        Parameters
        ----------
        path : str
            The path to the image.
        filter : None or int
            If None, no filter is applied. If an int, a median filter of size n pixels is applied. Generally, a good idea.
            Default is None for backwards compatibility.
        """
        self.image = image.load(path)
        if filter:
            self.image.filter(size=filter)
        self.vert_profile = SingleProfile(np.empty(0))
        self.horiz_profile = SingleProfile(np.empty(0))
        self.infield_area: float = 0.8
        self.parameters: dict = {}
        self.positions: dict = {}
        self._is_analyzed: bool = False
        self.image.check_inversion_by_histogram()

    @classmethod
    def from_demo_image(cls):
        """Load the demo image into an instance"""
        demo_file = retrieve_demo_file(url='flatsym_demo.dcm')
        return cls(demo_file)

    @staticmethod
    def run_demo():
        """Run the Flat/Sym demo by loading the demo image, print results, and plot the profiles."""
        fs = FieldParams.from_demo_image()
        fs.analyze(protocol='varian')
        print(fs.results())
        fs.plot_analyzed_image()

    def analyze(self, protocol: str, vert_position: float=0.5, horiz_position: float=0.5,
                vert_width: float=0, horiz_width: float=0, invert=False):
        """Analyze the image to determine parameters such as flatness and symmetry.

        Parameters
        ----------
        protocol : {'varian', 'elekta', 'vom80', 'siemens', 'iec'}
            The analysis protocol. See :ref:`analysis_definitions` for equations.
        vert_position : float (0.0-1.0)
            The distance ratio of the image to sample. E.g. at the default of 0.5 the profile is extracted
            in the middle of the image. 0.0 is at the left edge of the image and 1.0 is at the right edge of the image.
        horiz_position : float (0.0-1.0)
            The distance ratio of the image to sample. E.g. at the default of 0.5 the profile is extracted
            in the middle of the image. 0.0 is at the top edge of the image and 1.0 is at the bottom edge of the image.
        vert_width : float (0.0-1.0)
            The width ratio of the image to sample. E.g. at the default of 0.0 a 1 pixel wide profile is extracted.
            0.0 would be 1 pixel wide and 1.0 would be the vertical image width.
        horiz_width : float (0.0-1.0)
            The width ratio of the image to sample. E.g. at the default of 0.0 a 1 pixel wide profile is extracted.
            0.0 would be 1 pixel wide and 1.0 would be the horizontal image width.
        invert : bool
            Whether to invert the image. Setting this to True will override the default inversion. This is useful if
            pylinac's automatic inversion is incorrect.
        """
        if invert:
            self.image.invert()
        protocol_params = PROTOCOLS[protocol]
        self.parameters['Method'] = protocol
        self.infield_area = self._get_infield_area(protocol)

        # get vertical (inline) parameters
        self.vert_profile = self._get_vert_profile(vert_position, vert_width)
        self.parameters['vertical'] = {}
        for param in protocol_params:
            calc = protocol_params[param]
            vresult = calc(self.vert_profile, self.infield_area)
            self.parameters['vertical'][param] = vresult

        # get horizontal (crossline) parameters
        self.horiz_profile = self._get_horiz_profile(horiz_position, horiz_width)
        self.parameters['horizontal'] = {}
        for param in protocol_params:
            calc = protocol_params[param]
            hresult = calc(self.horiz_profile, self.infield_area)
            self.parameters['horizontal'][param] = hresult

        self._is_analyzed = True

    def results(self, as_str=True) -> Union[str, list]:
        """Get the results of the analysis.

        Parameters
        ----------
        as_str : bool
            If True, return a single string.
            If False, return a list. Useful for PDF publishing.

        Return
        ------
        results : str or list
        """
        if not self._is_analyzed:
            raise NotAnalyzed("Image is not analyzed yet. Use analyze() first.")
        # return parameters
        results = [
            'Field Parameters',
            '',
            f"Analysis protocol: {self.parameters['Method']}",
            f"Normalisation method: {norm}",
            f"Interpolation is {'on' if interpolate else 'off'}",
            f'File: {self.image.truncated_path}',
            '']
        vert_params = self.parameters['vertical']
        horiz_params = self.parameters['horizontal']
        for param in vert_params:
            results.append(f"{param} vertical: {vert_params[param]:.2f}")
            results.append(f"{param} horizontal: {horiz_params[param]:.2f}")
            results.append('')

        if as_str:
            results = '\n'.join(result for result in results)
        return results

    def _get_infield_area(self, protocol):
        """Return the in field area as a proportion 0.0-1.0 of the field size. The in field area depends on the
        protocol and field size, but for now define it as 80%"""
        return 0.8

    def _get_vert_profile(self, vert_position: float, vert_width: float):
        left_edge = int(round(self.image.array.shape[1]*vert_position - self.image.array.shape[1]*vert_width/2))
        left_edge = max(left_edge, 0)  # clip to 0
        right_edge = int(round(self.image.array.shape[1]*vert_position + self.image.array.shape[1]*vert_width/2) + 1)
        right_edge = min(right_edge, self.image.array.shape[1])  # clip to image limit
        self.positions['vertical left'] = left_edge
        self.positions['vertical right'] = right_edge
        vert_profile = SingleProfile(np.sum(self.image.array[:, left_edge:right_edge], 1))
        vert_profile.dpmm = self.image.dpmm
        return vert_profile

    def _get_horiz_profile(self, horiz_position: float, horiz_width: float):
        bottom_edge = int(round(self.image.array.shape[0] * horiz_position - self.image.array.shape[0] * horiz_width / 2))
        bottom_edge = max(bottom_edge, 0)
        top_edge = int(round(self.image.array.shape[0] * horiz_position + self.image.array.shape[0] * horiz_width / 2) + 1)
        top_edge = min(top_edge, self.image.array.shape[0])
        self.positions['horizontal bottom'] = bottom_edge
        self.positions['horizontal top'] = top_edge
        horiz_profile = SingleProfile(np.sum(self.image.array[bottom_edge:top_edge, :], 0))
        horiz_profile.dpmm = self.image.dpmm
        return horiz_profile

    def publish_pdf(self, filename: str, notes: Union[str, list]=None, open_file: bool=False, metadata: dict=None):
        """Publish (print) a PDF containing the analysis, images, and quantitative results.

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
        if not self._is_analyzed:
            raise NotAnalyzed("Image is not analyzed yet. Use analyze() first.")
        canvas = pdf.PylinacCanvas(filename, page_title="Beam Parameter Analysis",
                                   metadata=metadata, metadata_location=(2, 5))
        # draw result text
        text = self.results(as_str=False)
        number_of_lines = len(text)
        i = 0
        while i < number_of_lines:
            if i > number_of_lines - 1:
                i = number_of_lines - 1
            canvas.add_text(text=text[i:i+39], location=(2, 25.5), font_size=14)
            canvas.add_new_page()
            i = i + 40

        # draw vertical profile
        data = io.BytesIO()
        self._save_plot(self._plot_vert, data)
        canvas.add_image(data, location=(-4, 12.5), dimensions=(28, 12))

        # draw horizontal profile
        data = io.BytesIO()
        self._save_plot(self._plot_horiz, data)
        canvas.add_image(data, location=(-4, 1), dimensions=(28, 12))

        # draw image on last page
        canvas.add_new_page()
        data = io.BytesIO()
        self._save_plot(self._plot_image, data, title="Image")
        canvas.add_image(data, location=(1, 2), dimensions=(18, 20))
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 5.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 5))
        canvas.finish()

        if open_file:
            open_path(filename)

    def plot_analyzed_image(self, show: bool=True):
        """Plot the analyzed image. Shows parameters such as flatness & symmetry.

        Parameters
        ----------
        show : bool
            Whether to show the plot when called.
        """
        if not self._is_analyzed:
            raise NotAnalyzed("Image is not analyzed yet. Use analyze() first.")
        # set up axes
        plt.ioff()
        image_ax = plt.subplot2grid((2, 2), (0, 1))
        vert_ax = plt.subplot2grid((2, 2), (1, 1))
        horiz_ax = plt.subplot2grid((2, 2), (0, 0))

        # plot parameters

        # plot image and profile lines
        self._plot_image(image_ax, title=osp.basename(self.image.path))
        self._plot_horiz(horiz_ax)
        self._plot_vert(vert_ax)

        plt.suptitle("Beam Parameters")
        if show:
            plt.show()

    def _plot_image(self, axis: plt.Axes=None, title: str=''):
        plt.ioff()
        if axis is None:
            fig, axis = plt.subplots()
        axis.imshow(self.image.array, cmap=get_dicom_cmap())
        # show vertical/axial profiles
        left_profile = self.positions['vertical left']
        right_profile = self.positions['vertical right']
        axis.axvline(left_profile, color='b')
        axis.axvline(right_profile, color='b')
        # show horizontal/transverse profiles
        bottom_profile = self.positions['horizontal bottom']
        top_profile = self.positions['horizontal top']
        axis.axhline(bottom_profile, color='b')
        axis.axhline(top_profile, color='b')
        _remove_ticklabels(axis)
        axis.set_title(title)

    def _plot_vert(self, axis: plt.Axes=None):
        plt.ioff()
        if axis is None:
            fig, axis = plt.subplots()
        axis.set_title("Vertical Profile")
        axis.plot(self.vert_profile.values)
        left_ifa, right_ifa = self.vert_profile.field_edges(self.infield_area)
        axis.axvline(left_ifa, color='g', linestyle='-.')
        axis.axvline(right_ifa, color='g', linestyle='-.')
        _remove_ticklabels(axis)

    def _plot_horiz(self, axis: plt.Axes=None):
        plt.ioff()
        if axis is None:
            fig, axis = plt.subplots()
        axis.set_title("Horizontal Profile")
        axis.plot(self.horiz_profile.values)
        left_ifa, right_ifa = self.horiz_profile.field_edges(self.infield_area)
        axis.axvline(left_ifa, color='g', linestyle='-.')
        axis.axvline(right_ifa, color='g', linestyle='-.')
        _remove_ticklabels(axis)

    def _plot_symmetry(self, direction: str, axis: plt.Axes=None):
        plt.ioff()
        if axis is None:
            fig, axis = plt.subplots()
        data = self.symmetry[direction.lower()]
        axis.set_title(direction.capitalize() + " Symmetry")
        axis.plot(data['profile'].values)
        # plot lines
        cax_idx = data['profile'].fwxm_center()
        axis.axvline(data['profile left'], color='g', linestyle='-.')
        axis.axvline(data['profile right'], color='g', linestyle='-.')
        axis.axvline(cax_idx, color='m', linestyle='-.')
        # plot symmetry array
        if not data['array'] == 0:
            twin_axis = axis.twinx()
            twin_axis.plot(range(cax_idx, data['profile right']), data['array'][int(round(len(data['array'])/2)):])
            twin_axis.set_ylabel("Symmetry (%)")
        _remove_ticklabels(axis)
        # plot profile mirror
        central_idx = int(round(data['profile'].values.size / 2))
        offset = cax_idx - central_idx
        mirror_vals = data['profile'].values[::-1]
        axis.plot(data['profile']._indices + 2 * offset, mirror_vals)

    @staticmethod
    def _save_plot(func, filename: str, **kwargs):
        func(**kwargs)
        plt.savefig(filename)


def _remove_ticklabels(axis: plt.Axes):
    axis.get_yaxis().set_ticklabels([])
    axis.get_xaxis().set_ticklabels([])
