"""Module for performing analysis of images or 2D arrays for parameters such as flatness and symmetry."""
import io
import os.path as osp
import warnings
from enum import Enum
from math import floor, ceil
from typing import Union, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from pylinac.core.utilities import open_path
from . import __version__
from .core import image, pdf
from .core.exceptions import NotAnalyzed
from .core.geometry import Rectangle
from .core.hill import Hill
from .core.io import retrieve_demo_file, SNCProfiler
from .core.profile import SingleProfile, Edge, Interpolation, Normalization
from .settings import get_dicom_cmap


def flatness_dose_difference(profile: SingleProfile, in_field_ratio: float = 0.8, **kwargs) -> float:
    """The Varian specification for calculating flatness. See :ref:`varian_protocol`. """
    try:
        dmax = profile.field_calculation(in_field_ratio=in_field_ratio, calculation='max')
        dmin = profile.field_calculation(in_field_ratio=in_field_ratio, calculation='min')
    except IOError:
        raise ValueError(
            "An error was encountered in the flatness calculation. The image is likely inverted. Try inverting the image before analysis with <instance>.image.invert().")
    flatness = 100 * abs(dmax - dmin) / (dmax + dmin)
    return flatness


def flatness_dose_ratio(profile: SingleProfile, in_field_ratio: float = 0.8, **kwargs) -> float:
    """The Elekta specification for calculating flatness. See :ref:`elekta_protocol`. """
    try:
        dmax = profile.field_calculation(in_field_ratio=in_field_ratio, calculation='max')
        dmin = profile.field_calculation(in_field_ratio=in_field_ratio, calculation='min')
    except ValueError:
        raise ValueError(
            "An error was encountered in the flatness calculation. The image is likely inverted. Try inverting the image before analysis with <instance>.image.invert().")
    flatness = 100 * (dmax / dmin)
    return flatness


def plot_flatness(instance, profile: SingleProfile, axis: plt.Axes) -> None:
    """Plot flatness parameters. Applies to both flatness dose ratio and dose difference."""
    data = profile.field_data(in_field_ratio=instance._in_field_ratio)
    axis.axhline(np.max(data['field values']), color='g', linestyle='-.', label='Flatness region')
    axis.axhline(np.min(data['field values']), color='g', linestyle='-.')


def symmetry_point_difference(profile: SingleProfile, in_field_ratio: float, **kwargs) -> float:
    """Calculation of symmetry by way of point difference equidistant from the CAX. See :ref:`varian_protocol`.

    A negative value means the right side is higher. A positive value means the left side is higher.
    """
    field = profile.field_data(in_field_ratio=in_field_ratio)
    field_values = field['field values']
    cax_value = field['beam center value (@rounded)']

    def calc_sym(lt, rt, cax) -> float:
        return 100 * (lt - rt) / cax

    # get value w/ max magnitude
    sym_vals = [calc_sym(lt, rt, cax_value) for lt, rt in zip(field_values, field_values[::-1])]
    max_sym_idx = np.argmax(np.abs(sym_vals))
    return sym_vals[max_sym_idx]


def plot_symmetry_point_difference(instance, profile: SingleProfile, axis: plt.Axes) -> None:
    """Plotting of the symmetry point difference."""

    def calc_sym(lt, rt, cax) -> float:
        return 100 * abs(lt - rt) / cax

    _plot_sym_common(instance, calc_sym, profile, axis, label='Symmetry (%)', padding=(5, 0.5))


def _plot_sym_common(instance, calc_func: callable, profile: SingleProfile, axis: plt.Axes, label: str, padding: tuple) -> None:
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
    """Plotting of the symmetry point difference quotient."""

    def calc_sym(lt, rt, _) -> float:
        return max(abs(lt / rt), abs(rt / lt))

    _plot_sym_common(instance, calc_sym, profile, axis, label='Symmetry (AU)', padding=(0.05, 0.01))


def symmetry_pdq_iec(profile: SingleProfile, in_field_ratio: float, **kwargs) -> float:
    """Symmetry calculation by way of PDQ IEC. See :ref:`elekta_protocol`.

    A negative value means the right side is higher. A positive value means the left side is higher.
    """
    field = profile.field_data(in_field_ratio=in_field_ratio)
    field_values = field['field values']

    def calc_sym(lt, rt) -> float:
        sym1 = lt/rt
        sym2 = rt/lt
        if abs(sym1) > abs(sym2):
            sign = np.sign(sym1)
        else:
            sign = np.sign(sym2)
        return max(abs(lt / rt), abs(rt / lt)) * sign

    sym_values = [calc_sym(lt, rt) for lt, rt in zip(field_values, field_values[::-1])]
    sym_idx = np.argmax(np.abs(sym_values))

    return sym_values[sym_idx]


def symmetry_area(profile: SingleProfile, in_field_ratio: float, **kwargs) -> float:
    """Ratio of the area under the left and right profile segments.  See :ref:`siemens_protocol`.

    A negative value indicates the right side is higher; a positive value indicates the left side is higher.
    """
    data = profile.field_data(in_field_ratio=in_field_ratio)
    cax_idx = data['beam center index (exact)'] - data['left index (exact)']
    area_left = np.sum(data['field values'][:floor(cax_idx)])
    area_right = np.sum(data['field values'][ceil(cax_idx):])
    symmetry = 100 * (area_left - area_right) / (area_left + area_right)
    return symmetry


def plot_symmetry_area(instance, profile: SingleProfile, axis: plt.Axes) -> None:
    """PLot the symmetry area."""
    data = profile.field_data(in_field_ratio=instance._in_field_ratio)
    cax_idx = data['beam center index (exact)']
    left_idx = data['left index (rounded)']
    right_idx = data['right index (rounded)']

    axis.fill_between(range(left_idx, floor(cax_idx)), data['field values'][:floor(cax_idx)-left_idx], color='green', alpha=0.1, label='Left Area')
    axis.fill_between(range(ceil(cax_idx), right_idx), data['field values'][ceil(cax_idx)-left_idx:], color='slateblue', alpha=0.1, label='Right Area')


varian_protocol = {
    'symmetry': {'calc': symmetry_point_difference, 'unit': '%', 'plot': plot_symmetry_point_difference},
    'flatness': {'calc': flatness_dose_difference, 'unit': '%', 'plot': plot_flatness},
}
elekta_protocol = {
    'symmetry': {'calc': symmetry_pdq_iec, 'unit': '', 'plot': plot_symmetry_pdq},
    'flatness': {'calc': flatness_dose_ratio, 'unit': '', 'plot': plot_flatness},
}
siemens_protocol = {
    'symmetry': {'calc': symmetry_area, 'unit': '', 'plot': plot_symmetry_area},
    'flatness': {'calc': flatness_dose_difference, 'unit': '', 'plot': plot_flatness},
}


class Protocol(Enum):
    """Protocols to analyze additional metrics of the field. See :ref:`analysis_definitions`"""
    NONE = {}
    VARIAN = varian_protocol
    SIEMENS = siemens_protocol
    ELEKTA = elekta_protocol


class Centering(Enum):
    """See :ref:`centering`"""
    MANUAL = 'manual'
    BEAM_CENTER = 'Beam center'
    GEOMETRIC_CENTER = 'Geometric center'


class FieldAnalysis:
    """Class for analyzing the various parameters of a radiation image, most commonly an open image from a linac.


    """

    def __init__(self, path: str, filter: Optional[int] = None):
        """


        Parameters
        ----------
        path
            The path to the image.
        filter
            If None, no filter is applied. If an int, a median filter of size n pixels is applied. Generally, a good idea.
            Default is None for backwards compatibility.


        Attributes
        ----------
        vert_profile : :class:`pylinac.core.profile.SingleProfile`
        horiz_profile : :class:`pylinac.core.profile.SingleProfile`
        image : :class:`pylinac.core.image.DicomImage`
        """
        self._path: str = path
        self.image: image.ImageLike = image.load(path)
        if filter:
            self.image.filter(size=filter)
        self.vert_profile: SingleProfile
        self.horiz_profile: SingleProfile
        self._is_analyzed: bool = False
        self._from_device: bool = False
        self.image.check_inversion_by_histogram()

    @classmethod
    def from_demo_image(cls):
        """Load the demo image into an instance."""
        demo_file = retrieve_demo_file(url='flatsym_demo.dcm')
        return cls(demo_file)

    @staticmethod
    def run_demo() -> None:
        """Run the Flat/Sym demo by loading the demo image, print results, and plot the profiles."""
        fs = FieldAnalysis.from_demo_image()
        fs.analyze(protocol=Protocol.VARIAN)
        print(fs.results())
        fs.plot_analyzed_image()

    def _determine_center(self, centering: Centering) -> Tuple[float, float]:
        """Determine the position ratio using a centering technique."""
        vert_sum = np.sum(self.image.array, axis=1)
        horiz_sum = np.sum(self.image.array, axis=0)
        v_prof = SingleProfile(vert_sum)
        h_prof = SingleProfile(horiz_sum)
        if centering == Centering.GEOMETRIC_CENTER:
            # horiz and vert appear switched, but it's because the center of the vert profile
            # is where to take the horizontal profile and vic versa
            horiz_ratio = v_prof.geometric_center()['index (exact)'] / len(v_prof.values)
            vert_ratio = h_prof.geometric_center()['index (exact)'] / len(h_prof.values)
        elif centering == Centering.BEAM_CENTER:
            horiz_ratio = v_prof.beam_center()['index (exact)'] / len(v_prof.values)
            vert_ratio = h_prof.beam_center()['index (exact)'] / len(h_prof.values)
        return vert_ratio, horiz_ratio

    def _extract_profiles(self, horiz_position, horiz_width, interpolation_resolution_mm, vert_position, vert_width,
                          edge_detection_method, edge_smoothing_ratio, ground, interpolation,
                          interpolation_resolution, normalization_method, centering, hill_window_ratio) -> None:
        """Figures out 1) where to extract the profiles from the image and 2) sets the profiles to instance attrs"""

        # calculate the horiz/vert extraction positions if necessary
        if centering in (Centering.BEAM_CENTER, Centering.GEOMETRIC_CENTER):
            vert_position, horiz_position = self._determine_center(centering)

        # calculate the profiles
        horiz_values, upper_h_idx, lower_h_idx = self._get_horiz_values(horiz_position, horiz_width)
        self._upper_h_index = upper_h_idx
        self._lower_h_index = lower_h_idx
        self.horiz_profile = SingleProfile(horiz_values, dpmm=self.image.dpmm, interpolation=interpolation,
                                           interpolation_resolution_mm=interpolation_resolution_mm, ground=ground,
                                           edge_detection_method=edge_detection_method,
                                           normalization_method=normalization_method,
                                           edge_smoothing_ratio=edge_smoothing_ratio,
                                           hill_window_ratio=hill_window_ratio)

        vert_values, left_v_idx, right_v_idx = self._get_vert_values(vert_position, vert_width)
        self._left_v_index = left_v_idx
        self._right_v_index = right_v_idx
        self.vert_profile = SingleProfile(vert_values, dpmm=self.image.dpmm, interpolation=interpolation,
                                          interpolation_resolution_mm=interpolation_resolution_mm, ground=ground,
                                          edge_detection_method=edge_detection_method,
                                          normalization_method=normalization_method,
                                          edge_smoothing_ratio=edge_smoothing_ratio,
                                          hill_window_ratio=hill_window_ratio)

    def analyze(self, protocol: Enum = Protocol.VARIAN,
                centering: Centering = Centering.BEAM_CENTER,
                vert_position: float = 0.5, horiz_position: float = 0.5,
                vert_width: float = 0, horiz_width: float = 0,
                in_field_ratio: float = 0.8,
                slope_exclusion_ratio: float = 0.2,
                invert: bool = False,
                is_FFF: bool = False,
                penumbra: Tuple[float, float] = (20, 80),
                interpolation: Interpolation = Interpolation.LINEAR, interpolation_resolution_mm: float = 0.1,
                ground: bool = True,
                normalization_method: Normalization = Normalization.BEAM_CENTER,
                edge_detection_method: Edge = Edge.INFLECTION_DERIVATIVE,
                edge_smoothing_ratio: float = 0.003,
                hill_window_ratio: float = 0.15,
                **kwargs) -> None:
        """Analyze the image to determine parameters such as field edges, penumbra, and/or flatness & symmetry.

        Parameters
        ----------
        protocol : :class:`~pylinac.field_analysis.Protocol`
            The analysis protocol. See :ref:`analysis_definitions` for equations.
        centering : :class:`~pylinac.field_analysis.Centering`
            The profile extraction position technique. Beam center will determine the beam center and take profiles through the middle.
            Geometric center will simply take profiles centered about the image in both axes.
            Manual will use the values of `vert_position` and `horiz_position` as the position.
            See :ref:`centering`.
        vert_position
            The distance ratio of the image to sample. E.g. at the default of 0.5 the profile is extracted
            in the middle of the image. 0.0 is at the left edge of the image and 1.0 is at the right edge of the image.

            .. note::

                This value only applies when centering is MANUAL.

        horiz_position
            The distance ratio of the image to sample. E.g. at the default of 0.5 the profile is extracted
            in the middle of the image. 0.0 is at the top edge of the image and 1.0 is at the bottom edge of the image.

            .. note::

                This value only applies when centering is MANUAL.

        vert_width
            The width ratio of the image to sample. E.g. at the default of 0.0 a 1 pixel wide profile is extracted.
            0.0 would be 1 pixel wide and 1.0 would be the vertical image width.
        horiz_width
            The width ratio of the image to sample. E.g. at the default of 0.0 a 1 pixel wide profile is extracted.
            0.0 would be 1 pixel wide and 1.0 would be the horizontal image width.
        in_field_ratio
            The ratio of the field width to use for protocol values. E.g. 0.8 means use the 80% field width.
        slope_exclusion_ratio
            This is the ratio of the field to use to 1) calculate the "top" of an FFF field as well as 2) exclude from the
            "slope" calculation of each side of the field. Alternatively, this also defines the area to use for the
            slope calculation. E.g. an `in_field_ratio` of 0.8 and `slope_exclusion_ratio` of 0.2 means the central 20% of the
            field is used to fit and calculate the "top", while the region on either side of the central 20% between the central
            80% is used to calculate a slope on either side using linear regression.

            .. note::

                While the "top" is always calculated, it will not be displayed in plots if the `is_FFF` parameter is false.

        invert
            Whether to invert the image. Setting this to True will override the default inversion. This is useful if
            pylinac's automatic inversion is incorrect.
        is_FFF
            This is a flag to display the "top" calculation and slopes on either side of the field.
        penumbra
            A tuple of (lower, higher) % of the penumbra to calculate. E.g. (20, 80) will calculate the penumbra width at 20% and 80%.

            .. note::

                The exact height of the penumbra depends on the edge detection method. E.g. FWHM will result in
                calculating penumbra at 20/80% of the field max, but if something like inflection is used, the penumbra
                height will be 20/50*100*inflection height and 80/50*100*inflection height.

        ground
            Whether to ground the profile (set min value to 0). Helpful most of the time.
        interpolation
            Interpolation technique to use. See :ref:`Interpolation`.
        interpolation_resolution_mm
            The resolution that the interpolation will scale to.
            E.g. if the native dpmm is 2 and the resolution is set to 0.1mm the data will be interpolated to have a new dpmm of 10 (1/0.1).
        normalization_method
            How to pick the point to normalize the data to. See :ref:`Normalization`.
        edge_detection_method
            The method by which to detect the field edge. FWHM is reasonable most of the time except for FFF beams.
            Inflection-derivative will use the max gradient to determine the field edge. Note that this may not be the
            50% height. In fact, for FFF beams it shouldn't be. Inflection methods are better for FFF and other unusual
            beam shapes. See :ref:`edge`.
        edge_smoothing_ratio
            The ratio of the length of the values to use as the sigma for a Gaussian filter applied before searching for
            the inflection. E.g. 0.005 with a profile of 1000 points will result in a sigma of 5.
            This helps make the inflection point detection more robust to noise. Increase for noisy data.
        hill_window_ratio
            The ratio of the field size to use as the window to fit the Hill function. E.g. 0.2 will using a window
            centered about each edge with a width of 20% the size of the field width. Only applies when the edge
            detection is ``INFLECTION_HILL``.
        kwargs
            Use these to pass parameters to custom protocol functions. See :ref:`custom_protocols`.
        """
        if is_FFF and edge_detection_method == Edge.FWHM:
            warnings.warn("Using FWHM for an FFF beam is not advised. Consider using INFLECTION_DERIVATIVE or INFLECTION_HILL")
        if invert:
            self.image.invert()

        self._analyze(edge_detection_method, edge_smoothing_ratio, ground, horiz_position, horiz_width, in_field_ratio,
                      interpolation, interpolation_resolution_mm, is_FFF, kwargs, normalization_method, penumbra,
                      protocol, slope_exclusion_ratio, vert_position, vert_width, centering, hill_window_ratio)

    def _analyze(self, edge_detection_method, edge_smoothing_ratio, ground, horiz_position, horiz_width, in_field_ratio,
                 interpolation, interpolation_resolution_mm, is_FFF, kwargs, normalization_method, penumbra, protocol,
                 slope_exclusion_ratio, vert_position, vert_width, centering, hill_window_ratio):
        self._protocol = protocol
        self._penumbra = penumbra
        self._centering = centering
        self._is_FFF: bool = is_FFF
        self._edge_detection = edge_detection_method
        self._in_field_ratio = in_field_ratio
        self._slope_exclusion_ratio = slope_exclusion_ratio
        self._hill_window_ratio = hill_window_ratio
        self._interpolation_method = interpolation
        self._extract_profiles(horiz_position, horiz_width, interpolation_resolution_mm, vert_position, vert_width,
                               edge_detection_method, edge_smoothing_ratio, ground, interpolation,
                               interpolation_resolution_mm,
                               normalization_method, centering, hill_window_ratio)
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
        """Get the results of the analysis.

        Parameters
        ----------
        as_str
            If True, return a simple string. If False, return a list of each line of text.
        """
        if not self._is_analyzed:
            raise NotAnalyzed("Image is not analyzed yet. Use analyze() first.")

        results = [
            'Field Analysis Results',
            '----------------------',
            f'File: {self._path}',
            f"Protocol: {self._protocol.name}",
            ]
        if not self._from_device:
            results += [f"Centering method: {self._centering.value}",]
        results += [
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
                f"Left gradient: {self._results['left penumbra (%/mm)']:3.2f}%/mm",
                f"Right gradient: {self._results['right penumbra (%/mm)']:3.2f}%/mm",
                f"Top gradient: {self._results['top penumbra (%/mm)']:3.2f}%/mm",
                f"Bottom gradient: {self._results['bottom penumbra (%/mm)']:3.2f}%/mm",
                "",
            ]
        results += [
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
            f"Bottom slope: {self._results['bottom slope']:3.3f}%/mm",
            f"Left slope: {self._results['left slope']:3.3f}%/mm",
            f"Right slope: {self._results['right slope']:3.3f}%/mm",
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

    def results_data(self) -> dict:
        """Present the results data and metadata as a dict."""
        data = dict()
        data['pylinac version'] = __version__
        data['protocol'] = self._protocol.name
        data['centering method'] = getattr(self._centering, 'value', None)
        data['normalization method'] = self.horiz_profile._norm_method.value
        data['interpolation'] = self.horiz_profile._interp_method.value
        data['edge detection method'] = self.horiz_profile._edge_method.value
        data.update(self._results)
        return data

    def _get_vert_values(self, vert_position: float, vert_width: float) -> (np.ndarray, float, float):
        """Get the raw values of the profile to pass to SingleProfile"""
        left_edge = int(round(self.image.array.shape[1] * vert_position - self.image.array.shape[1] * vert_width / 2))
        left_edge = max(left_edge, 0)  # clip to 0
        right_edge = int(
                round(self.image.array.shape[1] * vert_position + self.image.array.shape[1] * vert_width / 2) + 1)
        right_edge = min(right_edge, self.image.array.shape[1])  # clip to image limit
        return np.mean(self.image.array[:, left_edge:right_edge], 1), left_edge, right_edge

    def _get_horiz_values(self, horiz_position: float, horiz_width: float) -> (np.ndarray, float, float):
        """Get the raw values of the profile to pass to SingleProfile"""
        bottom_edge = int(round(self.image.array.shape[0] * horiz_position - self.image.array.shape[0] * horiz_width / 2))
        bottom_edge = max(bottom_edge, 0)
        top_edge = int(
                round(self.image.array.shape[0] * horiz_position + self.image.array.shape[0] * horiz_width / 2) + 1)
        top_edge = min(top_edge, self.image.array.shape[0])
        return np.mean(self.image.array[bottom_edge:top_edge, :], 0), bottom_edge, top_edge

    def publish_pdf(self, filename: str, notes: Union[str, list] = None, open_file: bool = False,
                    metadata: dict = None) -> None:
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
            Extra stream to be passed and shown in the PDF. The key and value will be shown with a colon.
            E.g. passing {'Author': 'James', 'Unit': 'TrueBeam'} would result in text in the PDF like:
            --------------
            Author: James
            Unit: TrueBeam
            --------------
        """
        plt.ioff()
        if not self._is_analyzed:
            raise NotAnalyzed("Image is not analyzed yet. Use analyze() first.")
        canvas = pdf.PylinacCanvas(filename, page_title="Field Analysis",
                                   metadata=metadata, metadata_location=(2, 5))
        # draw result text
        text = self.results(as_str=False)
        number_of_lines = len(text)
        i = 0
        while i < number_of_lines:
            if i > number_of_lines - 1:
                i = number_of_lines - 1
            canvas.add_text(text=text[i:i + 50], location=(2, 25.5), font_size=10)
            canvas.add_new_page()
            i = i + 50

        # draw vertical profile
        stream = io.BytesIO()
        self._save_plot(self._plot_vert, stream)
        canvas.add_image(stream, location=(-4, 13), dimensions=(28, 12))

        # draw horizontal profile
        stream = io.BytesIO()
        self._save_plot(self._plot_horiz, stream)
        canvas.add_image(stream, location=(-4, 1), dimensions=(28, 12))

        # draw image on last page
        canvas.add_new_page()
        stream = io.BytesIO()
        self._save_plot(self._plot_image, stream, title="Image")
        canvas.add_image(stream, location=(1, 2), dimensions=(18, 20))
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
        if not self._from_device:
            image_ax = plt.subplot2grid((2, 2), (0, 1))
            vert_ax = plt.subplot2grid((2, 2), (1, 1))
            horiz_ax = plt.subplot2grid((2, 2), (0, 0))
        else:
            vert_ax = plt.subplot2grid((1, 2), (0, 1))
            horiz_ax = plt.subplot2grid((1, 2), (0, 0))

        # plot image and profile lines
        if not self._from_device:
            self._plot_image(image_ax, title=osp.basename(self._path))
        self._plot_vert(vert_ax, grid)
        self._plot_horiz(horiz_ax, grid)

        # plot legend
        lines = []
        labels = []
        v_lines, v_labels = vert_ax.get_legend_handles_labels()
        h_lines, h_labels = horiz_ax.get_legend_handles_labels()
        for label, line in zip(v_labels, v_lines):
            if label not in labels:
                lines.append(line)
                labels.append(label)
        for label, line in zip(h_labels, h_lines):
            if label not in labels:
                lines.append(line)
                labels.append(label)
        if not self._from_device:
            legend_ax = plt.subplot2grid((2, 2), (1, 0))
            legend_ax.legend(lines, labels, loc="center")
            legend_ax.axis('off')

            _remove_ticklabels(legend_ax)

        else:
            vert_ax.legend(lines, labels, loc='best',)
        plt.suptitle("Field Profile Analysis")
        if show:
            plt.show()

    def _plot_image(self, axis: plt.Axes = None, title: str = '') -> None:
        """Plot the image and profile extraction overlay"""
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

        # horizontal line/rect
        width_h = abs(self._upper_h_index-self._lower_h_index)
        center_h = (self.image.shape[1]/2, width / 2 + self._upper_h_index)
        r = Rectangle(width=self.image.shape[1],
                      height=width_h,
                      center=center_h)
        r.plot2axes(axis, edgecolor='b', fill=True, alpha=0.2, facecolor='b')

        # cleanup
        _remove_ticklabels(axis)
        axis.set_title(title)
        axis.legend()

    def _plot_vert(self, axis: plt.Axes = None, grid: bool = True) -> None:
        """Plot vertical profile"""
        if axis is None:
            fig, axis = plt.subplots()
        axis.grid(grid)
        axis.set_title("Vertical Profile")
        if self._from_device:
            axis.set_xlabel("detector")
            if self._interpolation_method == Interpolation.NONE:
                markers = "b+"
            else:
                markers = "b"
        else:
            axis.set_xlabel("pixels")
            markers = "b"
        axis.plot(self.vert_profile.values, markers, label='Profile')
        axis.set_ylabel("Normalized Response")

        # plot second axis w/ physical distance
        sec_y = axis.twiny()
        physical_distance = np.array(range(int(len(self.vert_profile.values)))) / self.vert_profile.dpmm
        sec_y.plot(physical_distance, self.vert_profile.values, markers)
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

    def _plot_horiz(self, axis: plt.Axes = None, grid: bool = True) -> None:
        """Plot horizontal profile"""
        if axis is None:
            fig, axis = plt.subplots()
        axis.grid(grid)
        axis.set_title("Horizontal Profile")
        if self._from_device:
            axis.set_xlabel("detector")
            if self._interpolation_method == Interpolation.NONE:
                markers = "b+"
            else:
                markers = "b"
        else:
            axis.set_xlabel("pixels")
            markers = "b"
        axis.plot(self.horiz_profile.values, markers, label='Profile')
        axis.set_ylabel("Normalized Response")

        # plot second axis w/ physical distance
        sec_y = axis.twiny()
        physical_distance = np.array(range(int(len(self.horiz_profile.values)))) / self.horiz_profile.dpmm
        sec_y.plot(physical_distance, self.horiz_profile.values, markers)
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
    def _save_plot(func, filename: Union[str, io.BytesIO], **kwargs) -> None:
        func(**kwargs)
        # plt.tight_layout(1.2)
        plt.savefig(filename)

    def _plot_penumbra(self, profile: SingleProfile, axis: plt.Axes = None) -> None:
        """Plot the non-linear regression fit against the profile"""
        data = profile.penumbra(self._penumbra[0], self._penumbra[1])
        axis.axvline(x=data[f'left {self._penumbra[0]}% index (exact)'], color='pink')
        axis.axvline(x=data[f'left {self._penumbra[1]}% index (exact)'], color='pink', label='Penumbra region')
        axis.axvline(x=data[f'right {self._penumbra[0]}% index (exact)'], color='pink')
        axis.axvline(x=data[f'right {self._penumbra[1]}% index (exact)'], color='pink')
        if self._edge_detection == Edge.INFLECTION_HILL:
            # plot left side Hill fit
            fw = profile.field_data(in_field_ratio=1.0)['width (exact)'] * self._hill_window_ratio / 2
            left_hill_idx = int(round(data[f'left {self._penumbra[0]}% index (exact)'] - fw))
            right_hill_idx = int(round(data[f'left {self._penumbra[1]}% index (exact)'] + fw))
            infl_data = profile.inflection_data()
            hill_fit = Hill.from_params(infl_data['left Hill params'])
            l_x_data = np.linspace(left_hill_idx, right_hill_idx, 200)
            axis.plot(l_x_data, hill_fit.y(l_x_data), color='black', label="Hill fit")

            # plot right side Hill fit
            left_hill_idx = int(round(data[f'right {self._penumbra[1]}% index (exact)'] - fw))
            right_hill_idx = int(round(data[f'right {self._penumbra[0]}% index (exact)'] + fw))
            hill_fit = Hill.from_params(infl_data['right Hill params'])
            r_x_data = np.linspace(left_hill_idx, right_hill_idx, 200)
            axis.plot(r_x_data, hill_fit.y(r_x_data), color='black', label="Hill fit")

    def _plot_field_edges(self, profile: SingleProfile, axis: plt.Axes) -> None:
        data = profile.field_data(in_field_ratio=1.0, slope_exclusion_ratio=self._slope_exclusion_ratio)
        axis.plot(data['left index (rounded)'], data['left value (@rounded)'], 'x', color='green', label='Field edge')
        axis.plot(data['right index (rounded)'], data['right value (@rounded)'], 'x', color='green')

    def _plot_infield_slope(self, profile: SingleProfile, axis: plt.Axes) -> None:
        data = profile.field_data(self._in_field_ratio, self._slope_exclusion_ratio)
        # left slope
        left_x_values = range(data['left index (rounded)'], data['left inner index (rounded)'])
        left_y_values = data['left slope'] * left_x_values + data['left intercept']
        axis.plot(left_x_values, left_y_values, color='tomato', label='in-field slope')
        # right slope
        right_x_values = range(data['right inner index (rounded)'], data['right index (rounded)'])
        right_y_values = data['right slope'] * right_x_values + data['right intercept']
        axis.plot(right_x_values, right_y_values, color='tomato')

    def _plot_top(self, profile: SingleProfile, axis: plt.Axes = None) -> None:
        """Plot a second order polynomial to the peak of the FFF field"""
        data = profile.field_data(self._in_field_ratio, self._slope_exclusion_ratio)
        x_model = np.linspace(data['left inner index (rounded)'], data['right inner index (rounded)'], 1000)
        y_model = data['top params'][0] * x_model ** 2 + data['top params'][1] * x_model + data['top params'][2]
        axis.plot(x_model, y_model, color='magenta', label='"top" polynomial fit')
        axis.plot(data['"top" index (exact)'], data['"top" value (@exact)'], 'x', color='magenta', label='"top" position')


class Device(Enum):
    """2D array device Enum.

    Attributes
    ----------
    PROFILER
    """
    PROFILER = {'device': SNCProfiler, 'detector spacing (mm)': 5}


class DeviceFieldAnalysis(FieldAnalysis):
    """Field analysis using a device array."""

    def __init__(self, path: str, device: Device):
        """
        Parameters
        ----------
        path
            Path to the file of the device output
        device
            The array device. Currently, the Profiler is supported. See :ref:`loading_device_data`.


        Attributes
        ----------
        vert_profile : :class:`pylinac.core.profile.SingleProfile`
        horiz_profile : :class:`pylinac.core.profile.SingleProfile`
        device : :class:`pylinac.field_analysis.Device`
        """
        self.device = device.value['device'](path=path)
        self._path = path
        self._from_device = True
        self._dpmm = 1/device.value['detector spacing (mm)']

    @classmethod
    def from_demo_image(cls) -> None:
        """Load the demo image into an instance."""
        demo_file = retrieve_demo_file(url='6fff.prm')
        return cls(demo_file, device=Device.PROFILER)

    @staticmethod
    def run_demo() -> None:
        """Run the Field analysis demo by loading the demo device dataset, print results, and plot the profiles."""
        fs = DeviceFieldAnalysis.from_demo_image()
        fs.analyze(protocol=Protocol.VARIAN, is_FFF=True)
        print(fs.results())
        fs.plot_analyzed_image()

    def analyze(self, protocol: Protocol = Protocol.VARIAN,
                in_field_ratio: float = 0.8,
                slope_exclusion_ratio: float = 0.3,
                is_FFF: bool = False,
                penumbra: tuple = (20, 80),
                interpolation: Interpolation = Interpolation.NONE, interpolation_resolution_mm: float = 0.1,
                ground: bool = True,
                normalization_method: Normalization = Normalization.GEOMETRIC_CENTER,
                edge_detection_method: Edge = Edge.INFLECTION_HILL,
                edge_smoothing_ratio: float = 0.003,
                hill_window_ratio: float = 0.15,
                **kwargs) -> None:
        """Analyze the device profiles to determine parameters such as field edges, penumbra, and/or flatness & symmetry.

        Parameters
        ----------
        protocol
            The analysis protocol. See :ref:`analysis_definitions` for equations and options.
        in_field_ratio
            The ratio of the field width to use for protocol values. E.g. 0.8 means use the 80% field width.
        slope_exclusion_ratio
            This is the ratio of the field to use to 1) calculate the "top" of an FFF field as well as 2) exclude from the
            "slope" calculation of each side of the field. Alternatively, this also defines the area to use for the
            slope calculation. E.g. an `in_field_ratio` of 0.8 and `slope_exclusion_ratio` of 0.2 means the central 20% of the
            field is used to fit and calculate the "top", while the region on either side of the central 20% between the central
            80% is used to calculate a slope on either side using linear regression.

            .. note::

                While the "top" is always calculated, it will not be displayed in plots if the `is_FFF` parameter is false.

        is_FFF
            This is a flag to display the "top" calculation and slopes on either side of the field.
        penumbra
            A tuple of (lower, higher) % of the penumbra to calculate. E.g. (20, 80) will calculate the penumbra width at 20% and 80%.

            .. note::

                The exact height of the penumbra depends on the edge detection method. E.g. FWHM will result in
                calculating penumbra at 20/80% of the field max, but if something like inflection is used, the penumbra
                height will be 20/50*100*inflection height and 80/50*100*inflection height.

        interpolation
            Interpolation technique to use. Must be one of the enum options of ``Interpolation``.
        ground
            Whether to ground the profile (set min value to 0). Helpful most of the time.
        interpolation_resolution_mm
            The resolution that the interpolation will scale to.
            E.g. if the native dpmm is 2 and the resolution is set to 0.1mm the data will be interpolated to have a new dpmm of 10 (1/0.1).
        normalization_method
            How to pick the point to normalize the data to.
        edge_detection_method
            The method by which to detect the field edge. FWHM is reasonable most of the time except for FFF beams.
            Inflection-derivative will use the max gradient to determine the field edge. Note that this may not be the
            50% height. In fact, for FFF beams it shouldn't be. Inflection methods are better for FFF and other unusual
            beam shapes.
        edge_smoothing_ratio
            The ratio of the length of the values to use as the sigma for a Gaussian filter applied before searching for
            the inflection. E.g. 0.005 with a profile of 1000 points will result in a sigma of 5.
            This helps make the inflection point detection more robust to noise. Increase for noisy data.
        hill_window_ratio
            The ratio of the field size to use as the window to fit the Hill function. E.g. 0.2 will using a window
            centered about each edge with a width of 20% the size of the field width. Only applies when the edge
            detection is ``INFLECTION_HILL``.
        kwargs
            Use these to pass parameters to custom protocol functions. See :ref:`custom_protocols`.
        """
        self._analyze(edge_detection_method, edge_smoothing_ratio, ground, None, None, in_field_ratio, interpolation,
                      interpolation_resolution_mm, is_FFF, kwargs, normalization_method, penumbra, protocol,
                      slope_exclusion_ratio, None, None, None, hill_window_ratio)

    def _extract_profiles(self, horiz_position, horiz_width, interpolation_resolution_mm, vert_position, vert_width,
                          edge_detection_method, edge_smoothing_ratio, ground, interpolation,
                          interpolation_resolution, normalization_method, centering, hill_window_ratio):
        # calculate the profiles from the device. Since it's not an image, there's no position values
        x_prof, y_prof, _, _ = self.device.to_profiles(dpmm=self._dpmm, interpolation=interpolation,
                                                       interpolation_resolution_mm=interpolation_resolution,
                                                       ground=ground, edge_detection_method=edge_detection_method,
                                                       normalization_method=normalization_method,
                                                       edge_smoothing_ratio=edge_smoothing_ratio,
                                                       hill_window_ratio=hill_window_ratio)
        self.vert_profile = y_prof
        self.horiz_profile = x_prof


def _remove_ticklabels(axis: plt.Axes):
    axis.get_yaxis().set_ticklabels([])
    axis.get_xaxis().set_ticklabels([])
