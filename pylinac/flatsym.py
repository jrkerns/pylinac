"""Module for analyzing images, either film or EPID, for flatness and symmetry."""
import io
import os.path as osp
from typing import Tuple, Union

import matplotlib.pyplot as plt

from pylinac.core.utilities import open_path
from .core.exceptions import NotAnalyzed
from .core.io import retrieve_demo_file
from .core.image import DicomImage
from .core.profile import SingleProfile
from .core import pdf
from .settings import get_dicom_cmap


def flatness_varian(profile: SingleProfile):
    """The Varian specification for calculating flatness"""
    dmax = profile.field_calculation(field_width=0.8, calculation='max')
    dmin = profile.field_calculation(field_width=0.8, calculation='min')
    flatness = 100 * abs(dmax - dmin) / (dmax + dmin)
    lt_edge, rt_edge = profile.field_edges()
    return flatness, dmax, dmin, lt_edge, rt_edge


def flatness_elekta(profile: SingleProfile):
    """The Elekta specification for calculating flatness"""
    dmax = profile.field_calculation(field_width=0.8, calculation='max')
    dmin = profile.field_calculation(field_width=0.8, calculation='min')
    flatness = 100 * (dmax / dmin)
    lt_edge, rt_edge = profile.field_edges()
    return flatness, dmax, dmin, lt_edge, rt_edge


def symmetry_point_difference(profile: SingleProfile):
    """Calculation of symmetry by way of point difference equidistant from the CAX"""
    values = profile.field_values(field_width=0.8)
    lt_edge, rt_edge = profile.field_edges(field_width=0.8)
    cax = profile.fwxm_center()
    dcax = profile.values[cax]
    max_val = 0
    sym_array = []
    for lt_pt, rt_pt in zip(values, values[::-1]):
        val = 100 * abs(lt_pt - rt_pt) / dcax
        sym_array.append(val)
        if val > max_val:
            max_val = val
    symmetry = max_val
    return symmetry, sym_array, lt_edge, rt_edge


def symmetry_pdq_iec(profile: SingleProfile):
    """Symmetry calculation by way of PDQ IEC"""
    values = profile.field_values(field_width=0.8)
    lt_edge, rt_edge = profile.field_edges(field_width=0.8)
    max_val = 0
    sym_array = []
    for lt_pt, rt_pt in zip(values, values[::-1]):
        val = max(abs(lt_pt / rt_pt), abs(rt_pt / lt_pt))
        sym_array.append(val)
        if val > max_val:
            max_val = val
    symmetry = 100 * max_val
    return symmetry, sym_array, lt_edge, rt_edge


SYMMETRY_EQUATIONS = {
    'varian': symmetry_point_difference,
    'point difference': symmetry_point_difference,
    'elekta': symmetry_pdq_iec,
    'pdq iec': symmetry_pdq_iec,
}
FLATNESS_EQUATIONS = {
    'varian': flatness_varian,
    'elekta': flatness_elekta,
    'vom80': flatness_varian,
    'siemens': flatness_varian,
    'iec': flatness_elekta,
}


class FlatSym(DicomImage):
    """Class for analyzing the flatness and symmetry of a radiation image, most commonly an open image from a linac.

    Attributes
    ----------
    symmetry : dict
        Contains the method of calculation and the vertical and horizontal symmetry data including the value.
    flatness : dict
        Contains the method of calculation and the vertical and horizontal flatness data including the value.
    positions : dict
        The position ratio used for analysis for vertical and horizontal.
    """

    def __init__(self, path: str):
        """
        Parameters
        ----------
        path : str
            The path to the image.
        """
        super().__init__(path)
        self.symmetry: dict = {}
        self.flatness: dict = {}
        self.positions: dict = {}
        self._is_analyzed: bool = False

    @classmethod
    def from_demo_image(cls):
        """Load the demo image into an instance"""
        demo_file = retrieve_demo_file(url='flatsym_demo.dcm')
        return cls(demo_file)

    @staticmethod
    def run_demo():
        """Run the Flat/Sym demo by loading the demo image, print results, and plot the profiles."""
        fs = FlatSym.from_demo_image()
        fs.analyze(flatness_method='varian', symmetry_method='varian')
        print(fs.results())
        fs.plot()

    def analyze(self, flatness_method: str, symmetry_method: str, vert_position: float=0.5, horiz_position: float=0.5):
        """Analyze the image to determine flatness & symmetry.

        Parameters
        ----------
        flatness_method : {'varian', 'elekta', 'vom80', 'siemens', 'iec'}
            The flatness algorithm. See :ref:`analysis_definitions` for equations.
        symmetry_method : {'varian', 'elekta', 'point difference', 'pdq iec'}
            The symmetry algorithm. See :ref:`analysis_definitions` for equations.
        vert_position : float (0.0-1.0)
            The distance ratio of the image to sample. E.g. at the default of 0.5 the profile is extracted
            in the middle of the image. 0.0 is at the left edge of the image and 1.0 is at the right edge of the image.
        horiz_position : float (0.0-1.0)
            The distance ratio of the image to sample. E.g. at the default of 0.5 the profile is extracted
            in the middle of the image. 0.0 is at the top edge of the image and 1.0 is at the bottom edge of the image.
        """
        self.symmetry = self._calc_symmetry(symmetry_method, vert_position, horiz_position)
        self.flatness = self._calc_flatness(flatness_method, vert_position, horiz_position)
        self.positions = {'vertical': vert_position, 'horizontal': horiz_position}
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
        # do some calculations
        horiz_penum = self.symmetry['horizontal']['profile'].penumbra_width() / self.dpmm
        vert_penum = self.symmetry['vertical']['profile'].penumbra_width() / self.dpmm
        horiz_width = self.symmetry['horizontal']['profile'].fwxm() / self.dpmm
        vert_width = self.symmetry['vertical']['profile'].fwxm() / self.dpmm
        upper_dist = abs(self.symmetry['vertical']['profile']._penumbra_point('left') - self.center.y) / self.dpmm
        lower_dist = abs(self.symmetry['vertical']['profile']._penumbra_point('right') - self.center.y) / self.dpmm
        left_dist = abs(self.symmetry['horizontal']['profile']._penumbra_point('left') - self.center.x) / self.dpmm
        right_dist = abs(self.symmetry['horizontal']['profile']._penumbra_point('right') - self.center.x) / self.dpmm
        results = [f'Flatness & Symmetry',
                   f'File: {self.truncated_path}',
                   "",
                   f'Flatness method: {self.flatness["method"].capitalize()}',
                   f"Vertical flatness: {self.flatness['vertical']['value']:3.3}%",
                   f"Horizontal flatness: {self.flatness['horizontal']['value']:3.3}%",
                   f'Symmetry method: {self.symmetry["method"].capitalize()}',
                   f"Vertical symmetry: {self.symmetry['vertical']['value']:3.3}%",
                   f"Horizontal symmetry: {self.symmetry['horizontal']['value']:3.3}%",
                   "",
                   "Penumbra (80/20):",
                   f"Horizontal: {horiz_penum:3.1f}mm",
                   f"Vertical: {vert_penum:3.1f}mm",
                   "",
                   "Field Size:",
                   f'Horizontal: {horiz_width:3.1f}mm',
                   f"Vertical: {vert_width:3.1f}mm",
                   "",
                   "CAX to edge distances:",
                   f"CAX -> Upper edge: {upper_dist:3.1f}mm",
                   f"CAX -> Lower edge: {lower_dist:3.1f}mm",
                   f"CAX -> Left edge: {left_dist:3.1f}mm",
                   f"CAX -> Right edge: {right_dist:3.1f}mm",
                   ]
        if as_str:
            results = '\n'.join(result for result in results)
        return results

    def _calc_symmetry(self, method: str, vert_position: float, horiz_position: float):
        vert_profile = SingleProfile(self.array[:, int(round(self.array.shape[1]*vert_position))])
        horiz_profile = SingleProfile(self.array[int(round(self.array.shape[0]*horiz_position)), :])
        # calc sym from profile
        symmetry_calculation = SYMMETRY_EQUATIONS[method.lower()]
        vert_sym, vert_sym_array, vert_lt, vert_rt = symmetry_calculation(vert_profile)
        horiz_sym , horiz_sym_array, horiz_lt, horiz_rt = symmetry_calculation(horiz_profile)
        return {
            'method': method,
            'horizontal': {
                'profile': horiz_profile, 'value': horiz_sym, 'array': horiz_sym_array, 'profile left': horiz_lt, 'profile right': horiz_rt,
                },
            'vertical': {
                'profile': vert_profile, 'value': vert_sym, 'array': vert_sym_array, 'profile left': vert_lt, 'profile right': vert_rt,
                },
        }

    def _calc_flatness(self, method: str, vert_position: float, horiz_position: float):
        vert_profile = SingleProfile(self.array[:, int(round(self.array.shape[1] * vert_position))])
        horiz_profile = SingleProfile(self.array[int(round(self.array.shape[0] * horiz_position)), :])
        # calc flatness from profile
        flatness_calculation = FLATNESS_EQUATIONS[method.lower()]
        vert_flatness, vert_max, vert_min, vert_lt, vert_rt = flatness_calculation(vert_profile)
        horiz_flatness, horiz_max, horiz_min, horiz_lt, horiz_rt = flatness_calculation(horiz_profile)
        return {
            'method': method,
            'horizontal': {
                'value': horiz_flatness, 'profile': horiz_profile, 'profile max': horiz_max, 'profile min': horiz_min, 'profile left': horiz_lt, 'profile right': horiz_rt,
            },
            'vertical': {
                'value': vert_flatness, 'profile': vert_profile, 'profile max': vert_max, 'profile min': vert_min, 'profile left': vert_lt, 'profile right': vert_rt,
            },
        }

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
        canvas = pdf.PylinacCanvas(filename, page_title="Flatness & Symmetry Analysis",
                                   metadata=metadata, metadata_location=(2, 5))
        # draw result text
        text = self.results(as_str=False)
        canvas.add_text(text=text, location=(2, 25.5), font_size=14)
        canvas.add_new_page()
        # draw flatness & symmetry on two pages
        for method in (self._plot_symmetry, self._plot_flatness):
            for height, direction in zip((1, 12.5), ('vertical', 'horizontal')):
                data = io.BytesIO()
                self._save_plot(method, data, direction=direction)
                canvas.add_image(data, location=(-4, height), dimensions=(28, 12))
            canvas.add_new_page()
        # draw image on last page
        data = io.BytesIO()
        self._save_plot(self._plot_image, data, title="Image")
        canvas.add_image(data, location=(1, 2), dimensions=(18, 20))
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 5.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 5))
        canvas.finish()

        if open_file:
            open_path(filename)

    def plot(self, show: bool=True):
        """Plot the analyzed image. Shows both flatness & symmetry.

        Parameters
        ----------
        show : bool
            Whether to show the plot when called.
        """
        if not self._is_analyzed:
            raise NotAnalyzed("Image is not analyzed yet. Use analyze() first.")
        # set up axes
        plt.ioff()
        vert_flat_ax = plt.subplot2grid((3, 3), (0, 0))
        vert_sym_ax = plt.subplot2grid((3, 3), (1, 0))
        image_ax = plt.subplot2grid((3, 3), (0, 1), colspan=2, rowspan=2)
        horiz_flat_ax = plt.subplot2grid((3, 3), (2, 1))
        horiz_sym_ax = plt.subplot2grid((3, 3), (2, 2))

        # plot flat/sym axis
        self._plot_flatness('vertical', vert_flat_ax)
        self._plot_flatness('horizontal', horiz_flat_ax)
        self._plot_symmetry('vertical', vert_sym_ax)
        self._plot_symmetry('horizontal', horiz_sym_ax)

        # plot image and profile lines
        self._plot_image(image_ax, title=osp.basename(self.path))

        plt.suptitle("Flatness & Symmetry")
        if show:
            plt.show()

    def _plot_image(self, axis: plt.Axes=None, title: str=''):
        plt.ioff()
        if axis is None:
            fig, axis = plt.subplots()
        axis.imshow(self.array, cmap=get_dicom_cmap())
        axis.axhline(self.positions['horizontal']*self.array.shape[0], color='r')  # y
        axis.axvline(self.positions['vertical']*self.array.shape[1], color='r')  # x
        _remove_ticklabels(axis)
        axis.set_title(title)

    def _plot_flatness(self, direction: str, axis: plt.Axes=None):
        plt.ioff()
        if axis is None:
            fig, axis = plt.subplots()
        data = self.flatness[direction.lower()]
        axis.set_title(direction.capitalize() + " Flatness")
        axis.plot(data['profile'].values)
        _remove_ticklabels(axis)
        axis.axhline(data['profile max'], color='r')
        axis.axhline(data['profile min'], color='r')
        axis.axvline(data['profile left'], color='g', linestyle='-.')
        axis.axvline(data['profile right'], color='g', linestyle='-.')

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
