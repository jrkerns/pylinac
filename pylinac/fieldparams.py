"""Module for analyzing images or 2D arrays for parameters such as flatness and symmetry.

   Adapted from FlatSym by Alan Chamberlain
   """
import enum
import io
import os.path as osp
from math import floor, ceil
from typing import Union, Optional

import matplotlib.pyplot as plt
import numpy as np

from pylinac.core.utilities import open_path
from . import __version__
from .core import image, pdf
from .core.exceptions import NotAnalyzed
from .core.geometry import Rectangle
from .core.io import retrieve_demo_file
from .core.profile import SingleProfile, Edge, Interpolation, Normalization
from .settings import get_dicom_cmap


# Field flatness parameters --------------------------------------------------------------------------------------------
def flatness_dose_difference(profile: SingleProfile, ifa: float = 0.8, **kwargs):
    """The Varian specification for calculating flatness."""
    try:
        dmax = profile.field_calculation(in_field_ratio=ifa)
        dmin = profile.field_calculation(in_field_ratio=ifa)
    except IOError:
        raise ValueError(
            "An error was encountered in the flatness calculation. The image is likely inverted. Try inverting the image before analysis with <instance>.image.invert().")
    flatness = 100 * abs(dmax - dmin) / (dmax + dmin)
    return flatness


def flatness_dose_ratio(profile: SingleProfile, ifa: float = 0.8, **kwargs):
    """The Elekta specification for calculating flatness."""
    try:
        dmax = profile.field_calculation(in_field_ratio=ifa, calculation='max')
        dmin = profile.field_calculation(in_field_ratio=ifa, calculation='min')
    except ValueError:
        raise ValueError(
            "An error was encountered in the flatness calculation. The image is likely inverted. Try inverting the image before analysis with <instance>.image.invert().")
    flatness = 100 * (dmax / dmin)
    return flatness


def plot_flatness(instance, profile: SingleProfile, axis: plt.Axes):
    """Plot flatness parameters"""
    data = profile.field_data(in_field_ratio=instance._in_field_ratio)
    axis.axhline(np.max(data['field values']), color='g', linestyle='-.', label='Flatness region')
    axis.axhline(np.min(data['field values']), color='g', linestyle='-.')


# Field symmetry parameters --------------------------------------------------------------------------------------------
def symmetry_point_difference(profile: SingleProfile, in_field_ratio: float, **kwargs) -> float:
    """Calculation of symmetry by way of point difference equidistant from the CAX."""
    field = profile.field_data(in_field_ratio=in_field_ratio)
    field_values = field['field values']
    cax_value = field['beam center value (@rounded)']

    def calc_sym(lt, rt, cax) -> float:
        return 100 * abs(lt - rt) / cax

    return max(calc_sym(lt, rt, cax_value) for lt, rt in zip(field_values, field_values[::-1]))


def plot_symmetry_point_difference(instance, profile: SingleProfile, axis: plt.Axes) -> None:
    """Calculation of symmetry by way of point difference equidistant from the CAX."""

    def calc_sym(lt, rt, cax) -> float:
        return 100 * abs(lt - rt) / cax

    _plot_sym_common(instance, calc_sym, profile, axis, label='Symmetry (%)', padding=(5, 0.5))


def _plot_sym_common(instance, calc_func, profile, axis, label: str, padding: tuple):
    field = profile.field_data(in_field_ratio=instance._in_field_ratio)
    field_values = field['field values']
    left_idx = field['left index (rounded)']
    right_idx = field['right index (rounded)']
    cax_value = field['beam center value (@rounded)']

    # same calc as PDQ and point difference, except we find the INDEX where the symmetry is maximum
    sym_values = [calc_func(lt, rt, cax_value) for lt, rt, _ in zip(field_values, field_values[::-1], range(int(round(len(field_values)/2))))]

    idx = np.argmax(sym_values)
    axis.plot(field['left index (rounded)']+idx, profile.values[field['left index (rounded)']+idx], '*', color='red', label='Symmetry max')
    axis.plot(field['right index (rounded)']-idx, profile.values[field['right index (rounded)']-idx], '*', color='red')
    sec_ax = axis.twinx()
    sec_ax.set_ylabel(label)

    # squish the secondary graph so it's not so large-looking
    ylim_top = max(sym_values) + padding[0]
    ylim_bottom = min(sym_values) - padding[1]
    sec_ax.set_ylim(ylim_bottom, ylim_top)
    left_end = int(round(left_idx+(right_idx-left_idx)/2))
    sec_ax.plot(range(left_end, left_end + len(sym_values)), sym_values[::-1])


def plot_symmetry_pdq(instance, profile: SingleProfile, axis: plt.Axes) -> None:
    """Calculation of symmetry by way of point difference equidistant from the CAX."""

    def calc_sym(lt, rt, _) -> float:
        return max(abs(lt / rt), abs(rt / lt))

    _plot_sym_common(instance, calc_sym, profile, axis, label='Symmetry (AU)', padding=(0.05, 0.01))


def symmetry_pdq_iec(profile: SingleProfile, in_field_ratio: float, **kwargs):
    """Symmetry calculation by way of PDQ IEC. Field calculation is automatically centered."""
    field = profile.field_data(in_field_ratio=in_field_ratio)
    field_values = field['field values']

    def calc_sym(lt, rt) -> float:
        return max(abs(lt / rt), abs(rt / lt))

    return max(calc_sym(lt, rt) for lt, rt in zip(field_values, field_values[::-1]))


def symmetry_area(profile: SingleProfile, in_field_ratio: float, **kwargs):
    """Ratio of the area under the left and right profile segments. Field is automatically centered."""
    data = profile.field_data(in_field_ratio=in_field_ratio)
    cax_idx = data['beam center index (exact)'] - data['left index (exact)']
    area_left = np.sum(data['field values'][:floor(cax_idx)])
    area_right = np.sum(data['field values'][ceil(cax_idx):])
    symmetry = 100 * abs(area_left - area_right) / (area_left + area_right)
    return symmetry


def plot_symmetry_area(instance, profile, axis):
    data = profile.field_data(in_field_ratio=instance._in_field_ratio)
    cax_idx = data['beam center index (exact)']
    left_idx = data['left index (rounded)']
    right_idx = data['right index (rounded)']

    axis.fill_between(range(left_idx, floor(cax_idx)), data['field values'][:floor(cax_idx)-left_idx], color='green', alpha=0.1, label='Left Area')
    axis.fill_between(range(ceil(cax_idx), right_idx), data['field values'][ceil(cax_idx)-left_idx:], color='slateblue', alpha=0.1, label='Right Area')


# VARIAN = {
    # 'left edge 50%: {:.1f} mm': left_edge_50,
    # 'right edge 50%: {:.1f} mm': right_edge_50,
    # 'field size: {:.1f} mm': field_size_edge_50,
    # 'field center: {:.1f} mm': field_center_edge_50,
    # 'penumbra 80-20% left: {:.1f} mm': penumbra_left_80_20,
    # 'penumbra 80-20% right: {:.1f} mm': penumbra_right_80_20,
#     'flatness: {:.2f} %': flatness_dose_difference,
#     'symmetry: {:.2f} %': symmetry_point_difference,
# }

# FLATSYM_VARIAN = {
#     'left edge 50%: {:.1f} mm': left_edge_50,
#     'right edge 50%: {:.1f} mm': right_edge_50,
#     'field size: {:.1f} mm': field_size_50,
#     'field center: {:.1f} mm': field_center_fwhm,
#     'penumbra 80-20% left: {:.1f} mm': penumbra_left_80_20,
#     'penumbra 80-20% right: {:.1f} mm': penumbra_right_80_20,
#     'flatness: {:.2f} %': flatness_dose_difference,
#     'symmetry: {:.2f} %': symmetry_point_difference,
# }
#
# ELEKTA = {
#     'left edge 50%: {:.1f} mm': left_edge_50,
#     'right edge 50%: {:.1f} mm': right_edge_50,
#     'field size: {:.1f} mm': field_size_50,
#     'field center: {:.1f} mm': field_center_fwhm,
#     'penumbra 80-20% left: {:.1f} mm': penumbra_left_80_20,
#     'penumbra 80-20% right: {:.1f} mm': penumbra_right_80_20,
#     'flatness: {:.2f} %': flatness_dose_ratio,
#     'symmetry: {:.2f} %': symmetry_pdq_iec,
#     # 'deviation diff: {:.2f} %': deviation_diff
# }
#
# SIEMENS = {
#     'left edge 50%: {:.1f} mm': left_edge_50,
#     'right edge 50%: {:.1f} mm': right_edge_50,
#     'field size: {:.1f} mm': field_size_50,
#     'field center: {:.1f} mm': field_center_fwhm,
#     'penumbra 80-20% left: {:.1f} mm': penumbra_left_80_20,
#     'penumbra 80-20% right: {:.1f} mm': penumbra_right_80_20,
#     'flatness: {:.2f} %': flatness_dose_difference,
#     'symmetry: {:.2f} %': symmetry_area,
#     # 'deviation max: {:.2f} %': deviation_max
# }
#
# VOM80 = {
#     'left edge 50%: {:.1f} mm': left_edge_50,
#     'right edge 50%: {:.1f} mm': right_edge_50,
#     'field size: {:.1f} mm': field_size_50,
#     'field center: {:.1f} mm': field_center_fwhm,
#     'penumbra 80-20% left: {:.1f} mm': penumbra_left_80_20,
#     'penumbra 80-20% right: {:.1f} mm': penumbra_right_80_20,
#     'flatness: {:.2f} %': flatness_dose_difference
# }
#
# IEC9076 = {
#     'left edge 50%: {:.1f} mm': left_edge_50,
#     'right edge 50%: {:.1f} mm': right_edge_50,
#     'field size: {:.1f} mm': field_size_50,
#     'field center: {:.1f} mm': field_center_fwhm,
#     'penumbra 80-20% left: {:.1f} mm': penumbra_left_80_20,
#     'penumbra 80-20% right: {:.1f} mm': penumbra_right_80_20,
#     'flatness: {:.2f} %': flatness_dose_ratio,
#     'symmetry: {:.2f} %': symmetry_pdq_iec,
# }
#
# AFSSAPS_JORF = {
#     'left edge 50%: {:.1f} mm': left_edge_50,
#     'right edge 50%: {:.1f} mm': right_edge_50,
#     'field size: {:.1f} mm': field_size_edge_50,
#     'field center: {:.1f} mm': field_center_edge_50,
#     'penumbra 80-20% left: {:.1f} mm': penumbra_left_80_20,
#     'penumbra 80-20% right: {:.1f} mm': penumbra_right_80_20,
#     'flatness: {:.2f} %': flatness_dose_difference,
#     'symmetry: {:.2f} %': symmetry_pdq_iec,
#     # 'deviation max: {:.2f} %': deviation_max
# }
#
# DIN = {
#     'left edge 50%: {:.1f} mm': left_edge_50,
#     'right edge 50%: {:.1f} mm': right_edge_50,
#     'field size: {:.1f} mm': field_size_50,
#     'field center: {:.1f} mm': field_center_fwhm,
#     'penumbra 80-20% left: {:.1f} mm': penumbra_left_80_20,
#     'penumbra 80-20% right: {:.1f} mm': penumbra_right_80_20,
#     'flatness: {:.2f} %': flatness_dose_ratio,
#     'symmetry: {:.2f} %': symmetry_pdq_iec,
# }
#
# FFF = {
#     'left edge inf: {:.1f} mm': left_edge_infl,
#     'right edge inf: {:.1f} mm': right_edge_infl,
#     'field center infl: {:.1f} mm': field_center_edge_infl,
#     'field center in-field slope: {:.1f} mm': field_center_infield_slope,
#     'max position: {:.1f} mm': field_top,
#     'field size infl: {:.1f} mm': field_size_edge_infl,
#     'left penumbra: {:.1f} mm': penumbra_left_infl,
#     'right penumbra: {:.1f} mm': penumbra_right_infl,
#     'left penumbra slope {:.1f} %/mm': penumbra_slope_left_infl,
#     'right penumbra slope {:.1f} %/mm': penumbra_slope_right_infl,
#     'left dose point 20%: {:.1f} %': dose_point_left_20,
#     'right dose point 20%: {:.1f} %': dose_point_right_20,
#     'left dose point 50%: {:.1f} %': dose_point_left_50,
#     'right dose point 50%: {:.1f} %': dose_point_right_50,
#     'left dose point 80%: {:.1f} %': dose_point_left_80,
#     'right dose point 80%: {:.1f} %': dose_point_right_80,
#     'in-field slope left: {:.2f} %/mm': infield_slope_left,
#     'in-field slope right: {:.2f} %/mm': infield_slope_right
# }

varian_protocol = {
    'symmetry': {'calc': symmetry_point_difference, 'unit': '%', 'plot': plot_symmetry_point_difference},
    'flatness': {'calc': flatness_dose_difference, 'unit': '%', 'plot': plot_flatness},
}
elekta_protocol = {
    'symmetry': {'calc': symmetry_pdq_iec, 'unit': '', 'plot': plot_symmetry_pdq},
    'flatness': {'calc': flatness_dose_ratio, 'unit': '', 'plot': plot_flatness},
}
# TODO: check below
siemens_protocol = {
    'symmetry': {'calc': symmetry_area, 'unit': '', 'plot': plot_symmetry_area},
    'flatness': {'calc': flatness_dose_difference, 'unit': '', 'plot': plot_flatness},
}


class Protocol(enum.Enum):
    NONE = {}
    VARIAN = varian_protocol
    SIEMENS = siemens_protocol
    ELEKTA = elekta_protocol


class FieldParams:
    """Class for analyzing the various parameters of a radiation image, most commonly an open image from a linac.
    """

    def __init__(self, path: str, filter: Optional[int] = None):
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
        self.vert_profile: SingleProfile
        self.horiz_profile: SingleProfile
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
        # fs.analyze(protocol='varian')
        print(fs.results())
        fs.plot_analyzed_image()

    def analyze(self, protocol: Protocol,
                vert_position: float = 0.5, horiz_position: float = 0.5,
                vert_width: float = 0, horiz_width: float = 0,
                in_field_ratio=0.8,
                is_FFF=False,
                invert=False,
                penumbra=(20, 80), interpolation=Interpolation.LINEAR, interpolation_resolution=0.1,
                ground=True, normalization_method=Normalization.BEAM_CENTER, edge_detection_method=Edge.FWHM,
                edge_smoothing_ratio=0.003,
                slope_exclusion_ratio=0.2, **kwargs):
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
        self._protocol = protocol
        self._penumbra = penumbra
        self._is_FFF: bool = is_FFF
        self._edge_detection = edge_detection_method
        self._in_field_ratio = in_field_ratio
        self._slope_exclusion_ratio = slope_exclusion_ratio

        # calculate the profiles
        horiz_values, upper_h_idx, lower_h_idx = self._get_horiz_values(horiz_position, horiz_width)
        self._upper_h_index = upper_h_idx
        self._lower_h_index = lower_h_idx
        self.horiz_profile = SingleProfile(horiz_values, dpmm=self.image.dpmm, interpolation=interpolation,
                                           interpolation_resolution_mm=interpolation_resolution, ground=ground,
                                           edge_detection_method=edge_detection_method,
                                           normalization_method=normalization_method,
                                           edge_smoothing_ratio=edge_smoothing_ratio)

        vert_values, left_v_idx, right_v_idx = self._get_vert_values(vert_position, vert_width)
        self._left_v_index = left_v_idx
        self._right_v_index = right_v_idx
        self.vert_profile = SingleProfile(vert_values, dpmm=self.image.dpmm, interpolation=interpolation,
                                          interpolation_resolution_mm=interpolation_resolution, ground=ground,
                                          edge_detection_method=edge_detection_method,
                                          normalization_method=normalization_method,
                                          edge_smoothing_ratio=edge_smoothing_ratio)

        self._results = {}
        # calculate common field info
        v_pen = self.vert_profile.penumbra(penumbra[0], penumbra[1])
        h_pen = self.horiz_profile.penumbra(penumbra[0], penumbra[1])
        self._results['top penumbra (mm)'] = v_pen['left penumbra width (exact) mm']
        self._results['bottom penumbra (mm)'] = v_pen['right penumbra width (exact) mm']
        self._results['left penumbra (mm)'] = h_pen['left penumbra width (exact) mm']
        self._results['right penumbra (mm)'] = h_pen['right penumbra width (exact) mm']
        if edge_detection_method == Edge.INFLECTION_HILL:
            self._results['top penumbra (%/mm)'] = abs(v_pen['left gradient (exact) %/mm'])
            self._results['bottom penumbra (%/mm)'] = abs(v_pen['right gradient (exact) %/mm'])
            self._results['left penumbra (%/mm)'] = abs(h_pen['left gradient (exact) %/mm'])
            self._results['right penumbra (%/mm)'] = abs(h_pen['right gradient (exact) %/mm'])

        self._results['geometric center index (x, y)'] = (self.horiz_profile.geometric_center()['index (exact)'],
                                                          self.vert_profile.geometric_center()['index (exact)'])
        self._results['beam center index (x, y)'] = (self.horiz_profile.beam_center()['index (exact)'],
                                                     self.vert_profile.beam_center()['index (exact)'])
        self._results['field size vertical (mm)'] = self.vert_profile.field_data(in_field_ratio=1.0)['width (exact) mm']
        self._results['field size horizontal (mm)'] = self.horiz_profile.field_data(in_field_ratio=1.0)[
            'width (exact) mm']
        self._results['beam center->Top (mm)'] = self.vert_profile.field_data(in_field_ratio=1.0)[
            'left distance->beam center (exact) mm']
        self._results['beam center->Bottom (mm)'] = self.vert_profile.field_data(in_field_ratio=1.0)[
            'right distance->beam center (exact) mm']
        self._results['beam center->Left (mm)'] = self.horiz_profile.field_data(in_field_ratio=1.0)[
            'left distance->beam center (exact) mm']
        self._results['beam center->Right (mm)'] = self.horiz_profile.field_data(in_field_ratio=1.0)[
            'right distance->beam center (exact) mm']
        self._results['CAX->Top (mm)'] = self.vert_profile.field_data(in_field_ratio=1.0)[
            'left distance->CAX (exact) mm']
        self._results['CAX->Bottom (mm)'] = self.vert_profile.field_data(in_field_ratio=1.0)[
            'right distance->CAX (exact) mm']
        self._results['CAX->Left (mm)'] = self.horiz_profile.field_data(in_field_ratio=1.0)[
            'left distance->CAX (exact) mm']
        self._results['CAX->Right (mm)'] = self.horiz_profile.field_data(in_field_ratio=1.0)[
            'right distance->CAX (exact) mm']

        # if is_FFF:
        h_field_data = self.horiz_profile.field_data(in_field_ratio=in_field_ratio,
                                                     slope_exclusion_ratio=slope_exclusion_ratio)
        v_field_data = self.vert_profile.field_data(in_field_ratio=in_field_ratio,
                                                    slope_exclusion_ratio=slope_exclusion_ratio)
        self._results['"top" position index (x, y)'] = (
            h_field_data['"top" index (exact)'], v_field_data['"top" index (exact)'])
        self._results['"top" horizontal distance from CAX (mm)'] = h_field_data['"top"->CAX (exact) mm']
        self._results['"top" vertical distance from CAX (mm)'] = v_field_data['"top"->CAX (exact) mm']
        self._results['"top" horizontal distance from beam center (mm)'] = h_field_data['"top"->beam center (exact) mm']
        self._results['"top" vertical distance from beam center (mm)'] = v_field_data['"top"->beam center (exact) mm']
        self._results['left slope'] = h_field_data['left slope (%/mm)']
        self._results['right slope'] = h_field_data['right slope (%/mm)']
        self._results['top slope'] = v_field_data['left slope (%/mm)']
        self._results['bottom slope'] = v_field_data['right slope (%/mm)']

        # calculate protocol info
        for name, item in protocol.value.items():
            self._results[name + ' horizontal'] = item['calc'](self.horiz_profile, in_field_ratio, **kwargs)
            self._results[name + ' vertical'] = item['calc'](self.vert_profile, in_field_ratio, **kwargs)

        self._is_analyzed = True

    def results(self, as_str=True) -> str:
        """Get the results of the analysis."""
        if not self._is_analyzed:
            raise NotAnalyzed("Image is not analyzed yet. Use analyze() first.")

        results = [
            'Field Analysis Results',
            '----------------------',
            f'File: {self.image.truncated_path}',
            f"Protocol: {self._protocol.name}",
            f"Normalization method: {self.horiz_profile._norm_method.value}",
            f"Interpolation: {self.horiz_profile._interp_method.value}",
            f"Edge detection method: {self.horiz_profile._edge_method.value}",
            "",
            f"Penumbra width ({self._penumbra[0]}/{self._penumbra[1]}):",
            f"Left: {self._results['left penumbra (mm)']:3.1f}mm",
            f"Right: {self._results['right penumbra (mm)']:3.1f}mm",
            f"Top: {self._results['top penumbra (mm)']:3.1f}mm",
            f"Bottom: {self._results['bottom penumbra (mm)']:3.1f}mm",
            "",
            ]
        if self._edge_detection == Edge.INFLECTION_HILL:
            results += [
                "Penumbra gradients:",
                f"Left gradient: {self._results['left penumbra (%/mm)']:3.4f}%/mm",
                f"Right gradient: {self._results['right penumbra (%/mm)']:3.4f}%/mm",
                f"Top gradient: {self._results['top penumbra (%/mm)']:3.4f}%/mm",
                f"Bottom gradient: {self._results['bottom penumbra (%/mm)']:3.4f}%/mm",
            ]
        results += [
            "",
            "Field Size:",
            f"Horizontal: {self._results['field size horizontal (mm)']:3.1f}mm",
            f"Vertical: {self._results['field size vertical (mm)']:3.1f}mm",
            "",
            "CAX to edge distances:",
            f"CAX -> Top edge: {self._results['CAX->Top (mm)']:3.1f}mm",
            f"CAX -> Bottom edge: {self._results['CAX->Bottom (mm)']:3.1f}mm",
            f"CAX -> Left edge: {self._results['CAX->Left (mm)']:3.1f}mm",
            f"CAX -> Right edge: {self._results['CAX->Right (mm)']:3.1f}mm",
            ""]
        if self._is_FFF:
            results += [
                "'Top' vertical distance from CAX: {:3.1f}mm".format(
                        self._results['"top" vertical distance from CAX (mm)']),
                "'Top' horizontal distance from CAX: {:3.1f}mm".format(
                        self._results['"top" horizontal distance from CAX (mm)']),
                "'Top' vertical distance from beam center: {:3.1f}mm".format(
                        self._results['"top" vertical distance from beam center (mm)']),
                "'Top' horizontal distance from beam center: {:3.1f}mm".format(
                        self._results['"top" horizontal distance from beam center (mm)']),
                "", ]
        results += [
            f"Top slope: {self._results['top slope']:3.3f}%/mm",
            f"Bottom slope: {self._results['bottom slope']:3.4f}%/mm",
            f"Left slope: {self._results['left slope']:3.5f}%/mm",
            f"Right slope: {self._results['right slope']:3.6f}%/mm",
            "",
            "Protocol data:",
            "--------------",
        ]

        for name, item in self._protocol.value.items():
            results.append(f"Vertical {name}: {self._results[name + ' vertical']:3.3f}{item['unit']}")
            results.append(
                f"Horizontal {name}: {self._results[name + ' horizontal']:3.3f}{item['unit']}")
            results.append('')

        if as_str:
            results = '\n'.join(result for result in results)
        return results

    def results_data(self):
        """Present the results data and metadata as a dict."""
        data = dict()
        data['pylinac version'] = __version__
        data['protocol'] = self._protocol.name
        data['normalization method'] = self.horiz_profile._norm_method.value
        data['interpolation'] = self.horiz_profile._interp_method.value
        data['edge detection method'] = self.horiz_profile._edge_method.value
        data.update(self._results)
        return data

    def _get_vert_values(self, vert_position: float, vert_width: float) -> (np.ndarray, float, float):
        left_edge = int(round(self.image.array.shape[1] * vert_position - self.image.array.shape[1] * vert_width / 2))
        left_edge = max(left_edge, 0)  # clip to 0
        right_edge = int(
                round(self.image.array.shape[1] * vert_position + self.image.array.shape[1] * vert_width / 2) + 1)
        right_edge = min(right_edge, self.image.array.shape[1])  # clip to image limit
        return np.mean(self.image.array[:, left_edge:right_edge], 1), left_edge, right_edge

    def _get_horiz_values(self, horiz_position: float, horiz_width: float) -> (np.ndarray, float, float):
        bottom_edge = int(round(self.image.array.shape[0] * horiz_position - self.image.array.shape[0] * horiz_width / 2))
        bottom_edge = max(bottom_edge, 0)
        top_edge = int(
                round(self.image.array.shape[0] * horiz_position + self.image.array.shape[0] * horiz_width / 2) + 1)
        top_edge = min(top_edge, self.image.array.shape[0])
        return np.mean(self.image.array[bottom_edge:top_edge, :], 0), bottom_edge, top_edge

    def publish_pdf(self, filename: str, notes: Union[str, list] = None, open_file: bool = False,
                    metadata: dict = None):
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
            canvas.add_text(text=text[i:i + 39], location=(2, 25.5), font_size=12)
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

    def plot_analyzed_image(self, show: bool = True, grid: bool = True):
        """Plot the analyzed image. Shows parameters such as flatness & symmetry.

        Parameters
        ----------
        show
            Whether to show the plot when called.
        grid
            Whether to show a grid on the profile plots
        """
        if not self._is_analyzed:
            raise NotAnalyzed("Image is not analyzed yet. Use analyze() first.")
        # set up axes
        plt.ioff()
        image_ax = plt.subplot2grid((2, 2), (0, 1))
        vert_ax = plt.subplot2grid((2, 2), (1, 1))
        horiz_ax = plt.subplot2grid((2, 2), (0, 0))
        if grid:
            vert_ax.grid(True)
            horiz_ax.grid(True)

        # plot image and profile lines
        self._plot_image(image_ax, title=osp.basename(self.image.path))
        self._plot_vert(vert_ax)
        self._plot_horiz(horiz_ax)

        # plot legend
        lines = []
        labels = []
        v_lines, v_labels = vert_ax.get_legend_handles_labels()
        h_lines, h_labels = horiz_ax.get_legend_handles_labels()
        for line, label in zip((v_lines + h_lines), (v_labels + h_labels)):
            if line not in lines: lines.append(line)
            if label not in labels: labels.append(label)
        legend_ax = plt.subplot2grid((2, 2), (1, 0))
        legend_ax.legend(lines, labels, loc="center")
        legend_ax.axis('off')

        _remove_ticklabels(legend_ax)
        plt.suptitle("Field Profile Analysis")
        if show:
            plt.show()

    def _plot_image(self, axis: plt.Axes = None, title: str = ''):
        plt.ioff()
        if axis is None:
            fig, axis = plt.subplots()
        axis.imshow(self.image.array, cmap=get_dicom_cmap())

        # vertical line/rect
        width = abs(self._left_v_index-self._right_v_index)
        center = (width / 2 + self._left_v_index, self.image.shape[0]/2)
        r = Rectangle(width=width,
                      height=self.image.shape[0],
                      center=center)
        r.plot2axes(axis, edgecolor='b', fill=True, alpha=0.2, facecolor='b', label="Profile Extraction Area")

        # show horizontal line/rect
        width_h = abs(self._upper_h_index-self._lower_h_index)
        center_h = (self.image.shape[1]/2, width / 2 + self._upper_h_index)
        r = Rectangle(width=self.image.shape[1],
                      height=width_h,
                      center=center_h)
        r.plot2axes(axis, edgecolor='b', fill=True, alpha=0.2, facecolor='b')
        _remove_ticklabels(axis)
        axis.set_title(title)
        axis.legend()

    def _plot_vert(self, axis: plt.Axes = None):
        """Plot vertical profile"""
        if axis is None:
            fig, axis = plt.subplots()
        axis.set_title("Vertical Profile")
        axis.plot(self.vert_profile.values, label='Profile')
        axis.set_xlabel("pixels")
        axis.set_ylabel("Normalized Response")

        # plot second axis w/ physical distance
        sec_y = axis.twiny()
        physical_distance = np.array(range(int(len(self.vert_profile.values)))) / self.vert_profile.dpmm
        sec_y.plot(physical_distance, self.vert_profile.values)
        sec_y.set_xlabel("mm")

        # plot basic parameters on profile
        self._plot_penumbra(self.vert_profile, axis)
        self._plot_field_edges(self.vert_profile, axis)
        if self._is_FFF:
            self._plot_top(self.vert_profile, axis)
            self._plot_infield_slope(self.vert_profile, axis)

        for name, item in self._protocol.value.items():
            if item.get("plot"):
                item['plot'](self, self.vert_profile, axis)

    def _plot_horiz(self, axis: plt.Axes = None):
        """Plot horizontal profile"""
        if axis is None:
            fig, axis = plt.subplots()
        axis.set_title("Horizontal Profile")
        axis.plot(self.horiz_profile.values, label='Profile')
        axis.set_xlabel("pixels")
        axis.set_ylabel("Normalized Response")

        # plot second axis w/ physical distance
        sec_y = axis.twiny()
        physical_distance = np.array(range(int(len(self.horiz_profile.values)))) / self.horiz_profile.dpmm
        sec_y.plot(physical_distance, self.horiz_profile.values)
        sec_y.set_xlabel("mm")

        # plot basic parameters on profile
        self._plot_penumbra(self.horiz_profile, axis)
        self._plot_field_edges(self.horiz_profile, axis)
        if self._is_FFF:
            self._plot_top(self.horiz_profile, axis)
            self._plot_infield_slope(self.horiz_profile, axis)

        for name, item in self._protocol.value.items():
            if item.get("plot"):
                item['plot'](self, self.horiz_profile, axis)

    @staticmethod
    def _save_plot(func, filename: str, **kwargs):
        func(**kwargs)
        plt.savefig(filename)

    def _plot_penumbra(self, profile: SingleProfile, axis: plt.Axes = None):
        """Plot the non-linear regression fit against the profile"""
        data = profile.penumbra(self._penumbra[0], self._penumbra[1])
        axis.axvline(x=data[f'left {self._penumbra[0]}% index (exact)'], color='pink')
        axis.axvline(x=data[f'left {self._penumbra[1]}% index (exact)'], color='pink', label='Penumbra region')
        axis.axvline(x=data[f'right {self._penumbra[0]}% index (exact)'], color='pink')
        axis.axvline(x=data[f'right {self._penumbra[1]}% index (exact)'], color='pink')

    def _plot_field_edges(self, profile: SingleProfile, axis: plt.Axes):
        data = profile.field_data(in_field_ratio=1.0, slope_exclusion_ratio=self._slope_exclusion_ratio)
        axis.plot(data['left index (rounded)'], data['left value (@rounded)'], 'x', color='green', label='Field edge')
        axis.plot(data['right index (rounded)'], data['right value (@rounded)'], 'x', color='green')

    def _plot_infield_slope(self, profile: SingleProfile, axis: plt.Axes):
        data = profile.field_data(self._in_field_ratio, self._slope_exclusion_ratio)
        # left slope
        left_x_values = range(data['left index (rounded)'], data['left inner index (rounded)'])
        left_y_values = data['left slope'] * left_x_values + data['left intercept']
        axis.plot(left_x_values, left_y_values, color='yellow', label='in-field slope')
        # right slope
        right_x_values = range(data['right inner index (rounded)'], data['right index (rounded)'])
        right_y_values = data['right slope'] * right_x_values + data['right intercept']
        axis.plot(right_x_values, right_y_values, color='yellow')

    def _plot_top(self, profile: SingleProfile, axis: plt.Axes = None):
        """Plot a second order polynomial to the peak of the FFF field"""
        data = profile.field_data(self._in_field_ratio, self._slope_exclusion_ratio)
        x_model = np.linspace(data['left inner index (rounded)'], data['right inner index (rounded)'])
        y_model = data['top params'][0] * x_model ** 2 + data['top params'][1] * x_model + data['top params'][2]
        axis.plot(x_model, y_model, color='magenta', label='Polynomial fit')
        axis.plot(data['"top" index (rounded)'], data['"top" value (@rounded)'], 'x', color='magenta', label='Max position')


def _remove_ticklabels(axis: plt.Axes):
    axis.get_yaxis().set_ticklabels([])
    axis.get_xaxis().set_ticklabels([])
