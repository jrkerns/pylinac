"""Module for analyzing images or 2D arrays for parameters such as flatness and symmetry.

   Adapted from FlatSym by Alan Chamberlain
   """

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
from .core.hillreg import hill_reg, hill_func, inv_hill_func, deriv_hill_func

interpolate: bool = True
norm: str = 'max grounded'               # one of 'cax', 'max', 'cax grounded', 'max grounded'
pen_width: float = 20                    # penumbra width for sigmoid analysis

# ----------------------------------------------------------------------------------------------------------------------
# Parameter functions. Each parameter defined in a protocol list must have an associated parameter function that
# returns the value of the parameter.
# ----------------------------------------------------------------------------------------------------------------------


# Field edge parameters ------------------------------------------------------------------------------------------------
def left_edge_50(profile: SingleProfile, *args) -> float:
    """Return the position of the 50% of max dose value on the left of the profile."""
    left_edge = abs(profile.field_edges(1.0, 50, norm, interpolate)[0] - profile.center()[0])/profile.dpmm
    return left_edge


def right_edge_50(profile: SingleProfile, *args):
    """Return the position of the 50% of max dose value on the right of the profile."""
    right_edge = abs(profile.field_edges(1.0, 50, norm, interpolate)[1] - profile.center()[0])/profile.dpmm
    return right_edge


def left_edge_infl(profile: SingleProfile, *args):
    """Return the position of the inflection point on the left of the profile."""
    left_edge_idx = profile.infl_points(pen_width, 'left')[0]
    left_edge = abs(left_edge_idx - profile.center()[0])/profile.dpmm
    return left_edge


def right_edge_infl(profile: SingleProfile, *args):
    """Return the position of the inflection point on the right of the profile."""
    right_edge_idx = profile.infl_points(pen_width, 'right')[0]
    right_edge = abs(right_edge_idx - profile.center()[0])/profile.dpmm
    return right_edge


# Field size parameters ------------------------------------------------------------------------------------------------
def field_size_50(profile: SingleProfile, *args):
    """Return the field size at 50% of max dose. Not affected by the normalisation mode.

       Included for testing purposes"""
    return profile.fwxm(50)/profile.dpmm


def field_size_edge_50(profile: SingleProfile, *args):
    """Return the field size at 50% of max dose."""
    return right_edge_50(profile) + left_edge_50(profile)


def field_size_edge_infl(profile: SingleProfile, *args):
    """Return the field size at 50% of max dose."""
    return right_edge_infl(profile) + left_edge_infl(profile)


# Field centre parameters ----------------------------------------------------------------------------------------------
def field_center_fwhm(profile: SingleProfile, *args):
    """Field center as given by the center of the profile FWHM. Not affected by the normalisation mode.

       Included for testing purposes."""
    field_center = (profile.fwxm_center(50, interpolate)[0] - profile.center()[0])/profile.dpmm
    return field_center


def field_center_edge_50(profile: SingleProfile, *args):
    """Calculates the field center from the 50 dose max field edges. May be different from the field_center_fwxm."""
    return (right_edge_50(profile) - left_edge_50(profile))/2


def field_center_edge_infl(profile: SingleProfile, *args):
    """Calculates the field center from the inflection point field edges. May be different from the field_center_fwxm."""
    return (right_edge_infl(profile) - left_edge_infl(profile))/2


# Field penumbra parameters --------------------------------------------------------------------------------------------
def penumbra_left_80_20(profile: SingleProfile, *args):
    """Return the distance between the 80% and 20% max dose values on the left side of the profile."""
    left_penum = abs(profile.field_edges(1.0, 80, norm, interpolate)[0]
                     - profile.field_edges(1.0, 20, norm, interpolate)[0])/profile.dpmm
    return left_penum


def penumbra_right_80_20(profile: SingleProfile, *args):
    """Return the distance between the 80% and 20% max dose values on the right side of the profile."""
    right_penum = abs(profile.field_edges(1.0, 80, norm, interpolate)[1]
                      - profile.field_edges(1.0, 20, norm, interpolate)[1])/profile.dpmm
    return right_penum


def penumbra_left_infl(profile: SingleProfile, *args):
    """Left penumbra value.

    Returns the distance between the locations where the dose equals 0.4 times the dose at the inflection point
    and 1.6 times that dose on the left side of the profile."""
    infl_idx, fit_params = profile.infl_points(pen_width, 'left')
    infl_val = hill_func(infl_idx, fit_params[0], fit_params[1], fit_params[2], fit_params[3])
    upper_idx = inv_hill_func(infl_val*1.6, fit_params)
    lower_idx = inv_hill_func(infl_val*0.4, fit_params)
    left_penum = abs(upper_idx - lower_idx)/profile.dpmm
    return left_penum


def penumbra_right_infl(profile: SingleProfile, *args):
    """Right penumbra value.

    Returns  the distance between the locations where the dose equals 0.4 times the dose at the inflection point
    and 1.6 times that dose on the right side of the profile."""
    infl_idx, fit_params = profile.infl_points(pen_width, 'right')
    infl_val = hill_func(infl_idx, fit_params[0], fit_params[1], fit_params[2], fit_params[3])
    upper_idx = inv_hill_func(infl_val*1.6, fit_params)
    lower_idx = inv_hill_func(infl_val*0.4, fit_params)
    right_penum = abs(upper_idx - lower_idx)/profile.dpmm
    return right_penum


def penumbra_slope_left_infl(profile: SingleProfile, *args):
    """Slope at the inflection point on the left penumbra"""
    if norm in ['max', 'max grounded']:
        cax_val = profile.values.max()
    else:
        _, cax_val = profile.fwxm_center()
    left_edge_idx, fit_params = profile.infl_points(pen_width, 'left')
    inf_slope = 100*deriv_hill_func(left_edge_idx, fit_params)*profile.dpmm/cax_val
    return inf_slope


def penumbra_slope_right_infl(profile: SingleProfile, *args):
    """Slope at the inflection point on the right penumbra"""
    if norm in ['max', 'max grounded']:
        cax_val = profile.values.max()
    else:
        _, cax_val = profile.fwxm_center()
    right_edge_idx, fit_params = profile.infl_points(pen_width, 'right')
    inf_slope = 100*deriv_hill_func(right_edge_idx, fit_params)*profile.dpmm/cax_val
    return inf_slope


# Dose point values ----------------------------------------------------------------------------------------------------
def dose_point_left_20(profile: SingleProfile, *args):
    """Dose value relative to CAX at 20% of field size from CAX"""
    return profile.dose_point(rel_dist=0.2, pen_width=20, side='left', norm=norm, interpolate=interpolate)


def dose_point_right_20(profile: SingleProfile, *args):
    """Dose value relative to CAX at 20% of field size from CAX"""
    return profile.dose_point(rel_dist=0.2, pen_width=20, side='right', norm=norm, interpolate=interpolate)


def dose_point_left_50(profile: SingleProfile, *args):
    """Dose value relative to CAX at 20% of field size from CAX"""
    return profile.dose_point(rel_dist=0.5, pen_width=20, side='left', norm=norm, interpolate=interpolate)


def dose_point_right_50(profile: SingleProfile, *args):
    """Dose value relative to CAX at 20% of field size from CAX"""
    return profile.dose_point(rel_dist=0.5, pen_width=20, side='right', norm=norm, interpolate=interpolate)


def dose_point_left_80(profile: SingleProfile, *args):
    """Dose value relative to CAX at 20% of field size from CAX"""
    return profile.dose_point(rel_dist=0.8, pen_width=20, side='left', norm=norm, interpolate=interpolate)


def dose_point_right_80(profile: SingleProfile, *args):
    """Dose value relative to CAX at 20% of field size from CAX"""
    return profile.dose_point(rel_dist=0.8, pen_width=20, side='right', norm=norm, interpolate=interpolate)


# Field flatness parameters --------------------------------------------------------------------------------------------
def flatness_dose_difference(profile: SingleProfile, ifa: float=0.8):
    """The Varian specification for calculating flatness."""
    try:
        dmax = profile.field_calculation(field_width=ifa, calculation='max')
        dmin = profile.field_calculation(field_width=ifa, calculation='min')
    except ValueError:
        raise ValueError("An error was encountered in the flatness calculation. The image is likely inverted. Try inverting the image before analysis with <instance>.image.invert().")
    flatness = 100 * abs(dmax - dmin)/(dmax + dmin)
    return flatness


def flatness_dose_ratio(profile: SingleProfile, ifa: float=0.8):
    """The Elekta specification for calculating flatness."""
    try:
        dmax = profile.field_calculation(field_width=ifa, calculation='max')
        dmin = profile.field_calculation(field_width=ifa, calculation='min')
    except ValueError:
        raise ValueError("An error was encountered in the flatness calculation. The image is likely inverted. Try inverting the image before analysis with <instance>.image.invert().")
    flatness = 100 * (dmax/dmin)
    return flatness


# Field symmetry parameters --------------------------------------------------------------------------------------------
def symmetry_point_difference(profile: SingleProfile, ifa: float=0.8):
    """Calculation of symmetry by way of point difference equidistant from the CAX.

    Field calculation is automatically centered."""
    values = profile.field_values(field_width=ifa)
    if norm in ['max', 'max grounded']:
        _, cax_val = profile.fwxm_center()
    else:
        _, cax_val = profile.center()
    sym_array = []
    for lt_pt, rt_pt in zip(values, values[::-1]):
        val = 100 * abs(lt_pt - rt_pt) / cax_val
        sym_array.append(val)
    symmetry = max(sym_array)
    return symmetry


def symmetry_pdq_iec(profile: SingleProfile, ifa: float = 0.8):
    """Symmetry calculation by way of PDQ IEC. Field calculation is automatically centered."""
    values = profile.field_values(field_width=ifa)
    max_val = 0
    for lt_pt, rt_pt in zip(values, values[::-1]):
        val = max(abs(lt_pt / rt_pt), abs(rt_pt / lt_pt))
        if val > max_val:
            max_val = val
    symmetry = 100 * max_val
    return symmetry


def symmetry_area(profile: SingleProfile, ifa: float = 0.8):
    """Ratio of the area under the left and right profile segments. Field is automatically centered."""
    values = profile.field_values(field_width=ifa)
    plen = len(values)
    cax_idx = round((plen - 1)/2)
    if plen % 2 == 0:                         # even number of values, CAX is straddled by centre values.
        area_left = np.sum(values[:cax_idx])
        area_right = np.sum(values[cax_idx:])
    else:                                     # include centre value on CAX
        area_left = np.sum(values[:cax_idx + 1])
        area_right = np.sum(values[cax_idx:])
    symmetry = 100*abs(area_left - area_right)/(area_left + area_right)
    return symmetry


# Field deviation parameters -------------------------------------------------------------------------------------------
def deviation_diff(profile: SingleProfile, ifa: float = 0.8):
    """Difference between the minimum and maximum."""
    if norm in ['max', 'max grounded']:
        _, cax_val = profile.fwxm_center()
    else:
        _, cax_val = profile.center()
    try:
        dmax = profile.field_calculation(field_width=ifa, calculation='max')
        dmin = profile.field_calculation(field_width=ifa, calculation='min')
    except ValueError:
        raise ValueError("An error was encountered in the deviation calculation. The image is likely inverted. Try inverting the image before analysis with <instance>.image.invert().")
    deviation = 100*(dmax - dmin)/cax_val
    return deviation


def deviation_max(profile: SingleProfile, ifa: float = 0.8):
    """Maximum deviation."""
    if norm in ['max', 'max grounded']:
        _, cax_val = profile.fwxm_center()
    else:
        _, cax_val = profile.center()
    try:
        dmax = profile.field_calculation(field_width=ifa, calculation='max')
    except ValueError:
        raise ValueError("An error was encountered in the deviation calculation. The image is likely inverted. Try inverting the image before analysis with <instance>.image.invert().")
    deviation = 100*dmax/cax_val
    return deviation


# ----------------------------------------------------------------------------------------------------------------------
# Predefined Protocols - Do not change these. Instead copy a protocol, give it a new name, put it after these protocols
# and add the protocol name to the dictionary PROTOCOLS. Parameter labels should include decimal places and units.
# ----------------------------------------------------------------------------------------------------------------------

ALL = {
    'left edge 50%: {:.1f} mm': left_edge_50,
    'right edge 50%: {:.1f} mm': right_edge_50,
    'left edge Inf: {:.1f} mm': left_edge_infl,
    'right edge Inf: {:.1f} mm': right_edge_infl,
    'field size FWHM: {:.1f} mm': field_size_50,
    'field size edge: {:.1f} mm': field_size_edge_50,
    'field size infl: {:.1f} mm': field_size_edge_infl,
    'field center edge: {:.1f} mm': field_center_edge_50,
    'field center FWHM: {:.1f} mm': field_center_fwhm,
    'field center infl: {:.1f} mm': field_center_edge_infl,
    'left penumbra 80-20%: {:.1f} mm': penumbra_left_80_20,
    'right penumbra 80-20%: {:.1f} mm': penumbra_right_80_20,
    'left penumbra infl: {:.1f} mm': penumbra_left_infl,
    'right penumbra infl: {:.1f} mm': penumbra_right_infl,
    'left penumbra slope {:.1f} %/mm': penumbra_slope_left_infl,
    'right penumbra slope {:.1f} %/mm': penumbra_slope_right_infl,
    'flatness diff: {:.2f} %': flatness_dose_difference,
    'flatness ratio: {:.2f} %': flatness_dose_ratio,
    'symmetry diff: {:.2f} %': symmetry_point_difference,
    'symmetry ratio: {:.2f} %': symmetry_pdq_iec,
    'symmetry area: {:.2f} %': symmetry_area,
    'deviation max: {:.2f} %': deviation_max,
    'deviation diff: {:.2f} %': deviation_diff,
    'left dose point 20%: {:.1f} %': dose_point_left_20,
    'right dose point 20%: {:.1f} %': dose_point_right_20,
    'left dose point 50%: {:.1f} %': dose_point_left_50,
    'right dose point 50%: {:.1f} %': dose_point_right_50,
    'left dose point 80%: {:.1f} %': dose_point_left_80,
    'right dose point 80%: {:.1f} %': dose_point_right_80
}

VARIAN = {
    'left edge 50%: {:.1f} mm': left_edge_50,
    'right edge 50%: {:.1f} mm': right_edge_50,
    'field size: {:.1f} mm': field_size_edge_50,
    'field center: {:.1f} mm': field_center_edge_50,
    'penumbra 80-20% left: {:.1f} mm': penumbra_left_80_20,
    'penumbra 80-20% right: {:.1f} mm': penumbra_right_80_20,
    'flatness: {:.2f} %': flatness_dose_difference,
    'symmetry: {:.2f} %': symmetry_point_difference,
}

FLATSYM_VARIAN = {
    'left edge 50%: {:.1f} mm': left_edge_50,
    'right edge 50%: {:.1f} mm': right_edge_50,
    'field size: {:.1f} mm': field_size_50,
    'field center: {:.1f} mm': field_center_fwhm,
    'penumbra 80-20% left: {:.1f} mm': penumbra_left_80_20,
    'penumbra 80-20% right: {:.1f} mm': penumbra_right_80_20,
    'flatness: {:.2f} %': flatness_dose_difference,
    'symmetry: {:.2f} %': symmetry_point_difference,
}

ELEKTA = {
    'left edge 50%: {:.1f} mm': left_edge_50,
    'right edge 50%: {:.1f} mm': right_edge_50,
    'field size: {:.1f} mm': field_size_50,
    'field center: {:.1f} mm': field_center_fwhm,
    'penumbra 80-20% left: {:.1f} mm': penumbra_left_80_20,
    'penumbra 80-20% right: {:.1f} mm': penumbra_right_80_20,
    'flatness: {:.2f} %': flatness_dose_ratio,
    'symmetry: {:.2f} %': symmetry_pdq_iec,
    'deviation diff: {:.2f} %': deviation_diff
}

SIEMENS = {
    'left edge 50%: {:.1f} mm': left_edge_50,
    'right edge 50%: {:.1f} mm': right_edge_50,
    'field size: {:.1f} mm': field_size_50,
    'field center: {:.1f} mm': field_center_fwhm,
    'penumbra 80-20% left: {:.1f} mm': penumbra_left_80_20,
    'penumbra 80-20% right: {:.1f} mm': penumbra_right_80_20,
    'flatness: {:.2f} %': flatness_dose_difference,
    'symmetry: {:.2f} %': symmetry_area,
    'deviation max: {:.2f} %': deviation_max
}

VOM80 = {
    'left edge 50%: {:.1f} mm': left_edge_50,
    'right edge 50%: {:.1f} mm': right_edge_50,
    'field size: {:.1f} mm': field_size_50,
    'field center: {:.1f} mm': field_center_fwhm,
    'penumbra 80-20% left: {:.1f} mm': penumbra_left_80_20,
    'penumbra 80-20% right: {:.1f} mm': penumbra_right_80_20,
    'flatness: {:.2f} %': flatness_dose_difference
}

IEC9076 = {
    'left edge 50%: {:.1f} mm': left_edge_50,
    'right edge 50%: {:.1f} mm': right_edge_50,
    'field size: {:.1f} mm': field_size_50,
    'field center: {:.1f} mm': field_center_fwhm,
    'penumbra 80-20% left: {:.1f} mm': penumbra_left_80_20,
    'penumbra 80-20% right: {:.1f} mm': penumbra_right_80_20,
    'flatness: {:.2f} %': flatness_dose_ratio,
    'symmetry: {:.2f} %': symmetry_pdq_iec,
}

AFSSAPS_JORF = {
    'left edge 50%: {:.1f} mm': left_edge_50,
    'right edge 50%: {:.1f} mm': right_edge_50,
    'field size: {:.1f} mm': field_size_edge_50,
    'field center: {:.1f} mm': field_center_edge_50,
    'penumbra 80-20% left: {:.1f} mm': penumbra_left_80_20,
    'penumbra 80-20% right: {:.1f} mm': penumbra_right_80_20,
    'flatness: {:.2f} %': flatness_dose_difference,
    'symmetry: {:.2f} %': symmetry_pdq_iec,
    'deviation max: {:.2f} %': deviation_max
}

DIN = {
    'left edge 50%: {:.1f} mm': left_edge_50,
    'right edge 50%: {:.1f} mm': right_edge_50,
    'field size: {:.1f} mm': field_size_50,
    'field center: {:.1f} mm': field_center_fwhm,
    'penumbra 80-20% left: {:.1f} mm': penumbra_left_80_20,
    'penumbra 80-20% right: {:.1f} mm': penumbra_right_80_20,
    'flatness: {:.2f} %': flatness_dose_ratio,
    'symmetry: {:.2f} %': symmetry_pdq_iec,
}


FFF = {
    'left edge inf: {:.1f} mm': left_edge_infl,
    'right edge inf: {:.1f} mm': right_edge_infl,
    'field center infl: {:.1f} mm': field_center_edge_infl,
    'field size infl: {:.1f} mm': field_size_edge_infl,
    'left penumbra: {:.1f} mm': penumbra_left_infl,
    'right penumbra: {:.1f} mm': penumbra_right_infl,
    'left penumbra slope {:.1f} %/mm': penumbra_slope_left_infl,
    'right penumbra slope {:.1f} %/mm': penumbra_slope_right_infl,
    'left dose point 20%: {:.1f} %': dose_point_left_20,
    'right dose point 20%: {:.1f} %': dose_point_right_20,
    'left dose point 50%: {:.1f} %': dose_point_left_50,
    'right dose point 50%: {:.1f} %': dose_point_right_50,
    'left dose point 80%: {:.1f} %': dose_point_left_80,
    'right dose point 80%: {:.1f} %': dose_point_right_80,
}
# ----------------------------------------------------------------------------------------------------------------------
# End of predefined protocols - Do not change these. Instead copy a protocol, give it a new name, put it after these
# protocols and add the protocol name to the dictionary PROTOCOLS.
# ----------------------------------------------------------------------------------------------------------------------

PROTOCOLS = {
    'all': ALL,
    'default': VARIAN,
    'varian': VARIAN,
    'flatsym_varian': FLATSYM_VARIAN,
    'elekta': ELEKTA,
    'siemens': SIEMENS,
    'vom80': VOM80,
    'iec9076': IEC9076,
    'afssaps-jorf': AFSSAPS_JORF,
    'din': DIN,
    'fff': FFF,
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
        """Load the demo image into an instance."""
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
            s_vert = param.format(vert_params[param])
            results.append(f"Vertical {s_vert}")
            s_horiz = param.format(horiz_params[param])
            results.append(f"Horizontal {s_horiz}")
            results.append('')

        if as_str:
            results = '\n'.join(result for result in results)
        return results

    def _get_infield_area(self, protocol):
        """Return the in field area as a proportion 0.0-1.0 of the field size. The in field area depends on the
        protocol and field size, but for now define it as 80%."""
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
        """Plot vertical profile"""
        plt.ioff()
        if axis is None:
            fig, axis = plt.subplots()
        axis.set_title("Vertical Profile")
        axis.plot(self.vert_profile.values)

        # plot parameters on profile
        protocol = self.parameters['Method']
        protocol_params = PROTOCOLS[protocol]
        if bool([val for key, val in protocol_params.items() if 'flatness' in key]):
            plot_flatness(self.vert_profile, self.infield_area, axis)
        if bool([val for key, val in protocol_params.items() if 'inf' in key]):
            plot_infl(self.vert_profile, axis)
        # _remove_ticklabels(axis)

    def _plot_horiz(self, axis: plt.Axes=None):
        """Plot horizontal profile"""
        plt.ioff()
        if axis is None:
            fig, axis = plt.subplots()
        axis.set_title("Horizontal Profile")
        axis.plot(self.horiz_profile.values)

        # plot parameters on profile
        protocol = self.parameters['Method']
        protocol_params = PROTOCOLS[protocol]
        if bool([val for key, val in protocol_params.items() if 'flatness' in key]):
            plot_flatness(self.horiz_profile, self.infield_area, axis)
        if bool([val for key, val in protocol_params.items() if 'infl' in key]):
            plot_infl(self.horiz_profile, axis)
        # _remove_ticklabels(axis)

    @staticmethod
    def _save_plot(func, filename: str, **kwargs):
        func(**kwargs)
        plt.savefig(filename)


# ----------------------------------------------------------------------------------------------------------------------
# Plot functions. These display various parameters graphically and are called depending on the protocol.
# ----------------------------------------------------------------------------------------------------------------------

def plot_flatness(profile: SingleProfile, ifa: float=0.8, axis: plt.Axes = None):
    """Plot flatness parameters"""
    left_ifa, right_ifa = profile.field_edges(ifa)
    axis.axvline(left_ifa, color='g', linestyle='-.')
    axis.axvline(right_ifa, color='g', linestyle='-.')


def plot_infl(profile: SingleProfile, axis: plt.Axes = None):
    """Plot the non-linear regression fit against the profile"""

    # plot left penumbra
    indices, values = profile.penumbra_values('left', pen_width)
    fit_params = hill_reg(indices, values)
    xModel = np.linspace(min(indices), max(indices))
    yModel = hill_func(xModel, *fit_params)
    axis.plot(xModel, yModel, color='r')

    # plot right penumbra
    indices, values = profile.penumbra_values('right', pen_width)
    fit_params = hill_reg(indices, values)
    xModel = np.linspace(min(indices), max(indices))
    yModel = hill_func(xModel, *fit_params)
    axis.plot(xModel, yModel, color='r')


def _plot_symmetry(profile: SingleProfile, axis: plt.Axes = None):
    plt.ioff()
    if axis is None:
        fig, axis = plt.subplots()
    data = profile.values
    axis.plot(data['profile'].values)
    # plot lines
    cax_idx = data['profile'].fwxm_center()
    axis.axvline(data['profile left'], color='g', linestyle='-.')
    axis.axvline(data['profile right'], color='g', linestyle='-.')
    axis.axvline(cax_idx, color='m', linestyle='-.')
    # plot symmetry array
    if not data['array'] == 0:
        twin_axis = axis.twinx()
        twin_axis.plot(range(cax_idx, data['profile right']), data['array'][int(round(len(data['array']) / 2)):])
        twin_axis.set_ylabel("Symmetry (%)")
    _remove_ticklabels(axis)
    # plot profile mirror
    central_idx = int(round(data['profile'].values.size / 2))
    offset = cax_idx - central_idx
    mirror_vals = data['profile'].values[::-1]
    axis.plot(data['profile']._indices + 2 * offset, mirror_vals)


def _remove_ticklabels(axis: plt.Axes):
    axis.get_yaxis().set_ticklabels([])
    axis.get_xaxis().set_ticklabels([])
