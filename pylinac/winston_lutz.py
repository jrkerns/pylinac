"""The Winston-Lutz module loads and processes EPID images that have acquired Winston-Lutz type images.

Features:

* **Automatic field & BB positioning** - When an image or directory is loaded, the field CAX and the BB
  are automatically found, along with the vector and scalar distance between them.
* **Isocenter size determination** - Using backprojections of the EPID images, the 3D gantry isocenter size
  and position can be determined *independent of the BB position*. Additionally, the 2D planar isocenter size
  of the collimator and couch can also be determined.
* **BB shift instructions** - Direct shift instructions can be printed for iterative BB placement.
  The current couch position can even be input to get the new couch values.
* **Axis deviation plots** - Plot the variation of the gantry, collimator, couch, and EPID in each plane
  as well as RMS variation.
* **File name interpretation** - Rename DICOM filenames to include axis information for linacs that don't include
  such information in the DICOM tags. E.g. "myWL_gantry45_coll0_couch315.dcm".
"""
from functools import lru_cache
from itertools import zip_longest
import io
import math
import os.path as osp
from typing import Union, List, Tuple
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage, optimize

from .core import image
from .core.decorators import value_accept
from .core.geometry import Point, Line, Circle, Vector, cos, sin
from .core.io import TemporaryZipDirectory, get_url, retrieve_demo_file, is_dicom_image
from .core.mask import filled_area_ratio, bounding_box
from .core import pdf
from .core.profile import SingleProfile
from .core.utilities import is_close, open_path

GANTRY = 'Gantry'
COLLIMATOR = 'Collimator'
COUCH = 'Couch'
COMBO = 'Combo'
EPID = 'Epid'
REFERENCE = 'Reference'
ALL = 'All'


class ImageManager(list):
    """Manages the images of a Winston-Lutz test."""
    def __init__(self, directory: str, use_filenames: bool):
        """
        Parameters
        ----------
        directory : str
            The path to the images.
        use_filenames: bool
            Whether to try to use the file name to determine axis values.
            Useful for Elekta machines that do not include that info in the DICOM data.
        """
        super().__init__()
        if isinstance(directory, list):
            for file in directory:
                if is_dicom_image(file):
                    img = WLImage(file, use_filenames)
                    self.append(img)
        elif not osp.isdir(directory):
            raise ValueError("Invalid directory passed. Check the correct method and file was used.")
        else:
            image_files = image.retrieve_image_files(directory)
            for file in image_files:
                img = WLImage(file, use_filenames)
                self.append(img)
        if len(self) < 2:
            raise ValueError("<2 valid WL images were found in the folder/file. Ensure you chose the correct folder/file for analysis")
        # reorder list based on increasing gantry angle
        self.sort(key=lambda i: (i.gantry_angle, i.collimator_angle, i.couch_angle))


class WinstonLutz:
    """Class for performing a Winston-Lutz test of the radiation isocenter."""
    images: ImageManager

    def __init__(self, directory: str, use_filenames: bool = False):
        """
        Parameters
        ----------
        directory : str
            Path to the directory of the Winston-Lutz EPID images.
        use_filenames: bool
            Whether to try to use the file name to determine axis values.
            Useful for Elekta machines that do not include that info in the DICOM data.

        Examples
        --------
        Run the demo:

            >>> WinstonLutz.run_demo()

        Load a directory with Winston-Lutz EPID images::

            >>> wl = WinstonLutz('path/to/directory')

        Load from a zip file::

            >>> wl = WinstonLutz.from_zip('path/to/images.zip')

        Or use the demo images provided::

            >>> wl = WinstonLutz.from_demo_images()

        Attributes
        ----------
        images : :class:`~pylinac.winston_lutz.ImageManager` instance
        """
        self.images = ImageManager(directory, use_filenames)

    @classmethod
    def from_demo_images(cls):
        """Instantiate using the demo images."""
        demo_file = retrieve_demo_file(url='winston_lutz.zip')
        return cls.from_zip(demo_file)

    @classmethod
    def from_zip(cls, zfile: str, use_filenames: bool=False):
        """Instantiate from a zip file rather than a directory.

        Parameters
        ----------
        zfile : str
            Path to the archive file.
        use_filenames : bool
            Whether to interpret axis angles using the filenames.
            Set to true for Elekta machines where the gantry/coll/couch data is not in the DICOM metadata.
        """
        with TemporaryZipDirectory(zfile) as tmpz:
            obj = cls(tmpz, use_filenames=use_filenames)
        return obj

    @classmethod
    def from_url(cls, url: str, use_filenames: bool = False):
        """Instantiate from a URL.

        Parameters
        ----------
        url : str
            URL that points to a zip archive of the DICOM images.
        use_filenames : bool
            Whether to interpret axis angles using the filenames.
            Set to true for Elekta machines where the gantry/coll/couch data is not in the DICOM metadata.
        """
        zfile = get_url(url)
        return cls.from_zip(zfile, use_filenames=use_filenames)

    @staticmethod
    def run_demo():
        """Run the Winston-Lutz demo, which loads the demo files, prints results, and plots a summary image."""
        wl = WinstonLutz.from_demo_images()
        print(wl.results())
        wl.plot_summary()

    @lru_cache()
    def _minimize_axis(self, axis=GANTRY):
        """Return the minimization result of the given axis."""
        def distance(p, things):
            """Calculate the maximum distance to any line from the given point."""
            other_thing = Point(p[0], p[1], p[2])
            if axis == COUCH:
                other_thing = Circle(other_thing, radius=p[3])
            return max(thing.distance_to(other_thing) for thing in things)

        if axis == GANTRY:
            attr = 'cax_line_projection'
        else:
            attr = 'cax2bb_vector'

        things = [getattr(image, attr) for image in self.images if image.variable_axis in (axis, REFERENCE)]
        if len(things) <= 1:
            raise ValueError("Not enough images of the given type to identify the axis isocenter")
        bounds = [(-30, 30), (-30, 30), (-30, 30), (0, 28)]  # search bounds for the optimization
        initial_guess = np.array([0, 0, 0, 0])
        result = optimize.minimize(distance, initial_guess, args=things, bounds=bounds)
        return result

    @property
    def gantry_iso_size(self) -> float:
        """The diameter of the 3D gantry isocenter size in mm. Only images where the collimator
        and couch were at 0 are used to determine this value."""
        num_gantry_like_images = self._get_images((GANTRY, REFERENCE))[0]
        if num_gantry_like_images > 1:
            return self._minimize_axis(GANTRY).fun * 2
        else:
            return 0

    @property
    def collimator_iso_size(self) -> float:
        """The 2D collimator isocenter size (diameter) in mm. The iso size is in the plane
        normal to the gantry."""
        num_collimator_like_images = self._get_images((COLLIMATOR, REFERENCE))[0]
        if num_collimator_like_images > 1:
            return self._minimize_axis(COLLIMATOR).fun * 2
        else:
            return 0

    @property
    def couch_iso_size(self) -> float:
        """The diameter of the 2D couch isocenter size in mm. Only images where
        the gantry and collimator were at zero are used to determine this value."""
        num_couch_like_images = self._get_images((COUCH, REFERENCE))[0]
        if num_couch_like_images > 1:
            return self._minimize_axis(COUCH).x[3] * 2
        else:
            return 0

    @property
    def bb_shift_vector(self) -> Vector:
        """The shift necessary to place the BB at the radiation isocenter"""
        vs = [img.cax2bb_vector3d for img in self.images]
        # only include the values that are clinically significant; rules out values that are along the beam axis
        xs = np.mean([v.x for v in vs if (0.02 < v.x or v.x < -0.02)])
        ys = np.mean([v.y for v in vs if (0.02 < v.y or v.y < -0.02)])
        zs = np.mean([v.z for v in vs if (0.02 < v.z or v.z < -0.02)])
        return Vector(-xs, -ys, -zs)

    def bb_shift_instructions(self, couch_vrt: float=None, couch_lng: float=None, couch_lat: float=None) -> str:
        """A string describing how to shift the BB to the radiation isocenter"""
        sv = self.bb_shift_vector
        x_dir = 'LEFT' if sv.x < 0 else 'RIGHT'
        y_dir = 'UP' if sv.y > 0 else 'DOWN'
        z_dir = 'IN' if sv.z < 0 else 'OUT'
        move = "{} {:2.2f}mm; {} {:2.2f}mm; {} {:2.2f}mm".format(x_dir, abs(sv.x), y_dir, abs(sv.y), z_dir, abs(sv.z))
        if all(val is not None for val in [couch_vrt, couch_lat, couch_lng]):
            new_lat = round(couch_lat + sv.x/10, 2)
            new_vrt = round(couch_vrt + sv.y/10, 2)
            new_lng = round(couch_lng - sv.z/10, 2)
            move += "\nNew couch coordinates (mm): VRT: {:3.2f}; LNG: {:3.2f}; LAT: {:3.2f}".format(new_vrt, new_lng, new_lat)
        return move

    @value_accept(axis=(GANTRY, COLLIMATOR, COUCH, EPID, COMBO), value=('all', 'range'))
    def axis_rms_deviation(self, axis: str=GANTRY, value: str='all'):
        """The RMS deviations of a given axis/axes.

        Parameters
        ----------
        axis : ('Gantry', 'Collimator', 'Couch', 'Epid', 'Combo'}
            The axis desired.
        value : {'all', 'range'}
            Whether to return all the RMS values from all images for that axis, or only return the maximum range of
            values, i.e. the 'sag'.
        """
        if axis != EPID:
            attr = 'bb_'
        else:
            attr = 'epid_'
        imgs = self._get_images(axis=axis)[1]
        if len(imgs) <= 1:
            return (0, )
        rms = []
        for img in imgs:
            rms.append(np.sqrt(sum(getattr(img, attr + ax + '_offset') ** 2 for ax in ('x', 'y', 'z'))))
        if value == 'range':
            rms = max(rms) - min(rms)
        return rms

    def cax2bb_distance(self, metric: str='max'):
        """The distance in mm between the CAX and BB for all images according to the given metric.

        Parameters
        ----------
        metric : {'max', 'median'}
            The metric of distance to use.
        """
        if metric == 'max':
            return max(image.cax2bb_distance for image in self.images)
        elif metric == 'median':
            return np.median([image.cax2bb_distance for image in self.images])

    def cax2epid_distance(self, metric: str='max'):
        """The distance in mm between the CAX and EPID center pixel for all images according to the given metric.

        Parameters
        ----------
        metric : {'max', 'median'}
            The metric of distance to use.
        """
        if metric == 'max':
            return max(image.cax2epid_distance for image in self.images)
        elif metric == 'median':
            return np.median([image.cax2epid_distance for image in self.images])

    @value_accept(item=(GANTRY, EPID, COLLIMATOR, COUCH))
    def _plot_deviation(self, item: str, ax: plt.Axes=None, show: bool=True):
        """Helper function: Plot the sag in Cartesian coordinates.

        Parameters
        ----------
        item : {'gantry', 'epid', 'collimator', 'couch'}
            The axis to plot.
        ax : None, matplotlib.Axes
            The axis to plot to. If None, creates a new plot.
        show : bool
            Whether to show the image.
        """
        title = 'Relative {} displacement'.format(item)
        if item == EPID:
            attr = 'epid'
            item = GANTRY
        else:
            attr = 'bb'
        # get axis images, angles, and shifts
        imgs = [image for image in self.images if image.variable_axis in (item, REFERENCE)]
        angles = [getattr(image, '{}_angle'.format(item.lower())) for image in imgs]
        z_sag = np.array([getattr(image, attr + '_z_offset') for image in imgs])
        y_sag = np.array([getattr(image, attr + '_y_offset') for image in imgs])
        x_sag = np.array([getattr(image, attr + '_x_offset') for image in imgs])
        rms = np.sqrt(x_sag**2+y_sag**2+z_sag**2)

        # plot the axis deviation
        if ax is None:
            ax = plt.subplot(111)
        ax.plot(angles, z_sag, 'bo', label='In/Out', ls='-.')
        ax.plot(angles, x_sag, 'm^', label='Left/Right', ls='-.')
        if item not in (COUCH, COLLIMATOR):
            ax.plot(angles, y_sag, 'r*', label='Up/Down', ls='-.')
        ax.plot(angles, rms, 'g+', label='RMS', ls='-')
        ax.set_title(title)
        ax.set_ylabel('mm')
        ax.set_xlabel("{} angle".format(item))
        ax.set_xticks(np.arange(0, 361, 45))
        ax.set_xlim([-15, 375])
        ax.grid(True)
        ax.legend(numpoints=1)
        if show:
            plt.show()

    def _get_images(self, axis: tuple=(GANTRY,)) -> Tuple[float, list]:
        if isinstance(axis, str):
            axis = (axis,)
        images = [image for image in self.images if image.variable_axis in axis]
        return len(images), images

    @value_accept(axis=(GANTRY, COLLIMATOR, COUCH, COMBO))
    def plot_axis_images(self, axis: str=GANTRY, show: bool=True, ax: plt.Axes=None):
        """Plot all CAX/BB/EPID positions for the images of a given axis.

        For example, axis='Couch' plots a reference image, and all the BB points of the other
        images where the couch was moving.

        Parameters
        ----------
        axis : {'Gantry', 'Collimator', 'Couch', 'Combo'}
            The images/markers from which accelerator axis to plot.
        show : bool
            Whether to actually show the images.
        ax : None, matplotlib.Axes
            The axis to plot to. If None, creates a new plot.
        """
        images = [image for image in self.images if image.variable_axis in (axis, REFERENCE)]
        ax = images[0].plot(show=False, ax=ax)  # plots the first marker; plot the rest of the markers below
        if axis != COUCH:
            # plot EPID
            epid_xs = [img.epid.x for img in images[1:]]
            epid_ys = [img.epid.y for img in images[1:]]
            ax.plot(epid_xs, epid_ys, 'b+', ms=8)
            # get CAX positions
            xs = [img.field_cax.x for img in images[1:]]
            ys = [img.field_cax.y for img in images[1:]]
            marker = 'gs'
        else:
            # get BB positions
            xs = [img.bb.x for img in images[1:]]
            ys = [img.bb.y for img in images[1:]]
            marker = 'ro'
        ax.plot(xs, ys, marker, ms=8)
        # set labels
        ax.set_title(axis + ' wobble')
        ax.set_xlabel(axis + ' positions superimposed')
        ax.set_ylabel(axis + " iso size: {0:3.2f}mm".format(getattr(self, axis.lower() + '_iso_size')))
        if show:
            plt.show()

    @value_accept(axis=(GANTRY, COLLIMATOR, COUCH, COMBO, ALL))
    def plot_images(self, axis: str=ALL, show: bool=True):
        """Plot a grid of all the images acquired.

        Four columns are plotted with the titles showing which axis that column represents.

        Parameters
        ----------
        axis : {'Gantry', 'Collimator', 'Couch', 'Combo', 'All'}
        show : bool
            Whether to show the image.
        """

        def plot_image(image, axis):
            """Helper function to plot a WLImage to an axis."""
            if image is None:
                axis.set_frame_on(False)
                axis.axis('off')
            else:
                image.plot(ax=axis, show=False)

        # get axis images
        if axis == GANTRY:
            images = [image for image in self.images if image.variable_axis in (GANTRY, REFERENCE)]
        elif axis == COLLIMATOR:
            images = [image for image in self.images if image.variable_axis in (COLLIMATOR, REFERENCE)]
        elif axis == COUCH:
            images = [image for image in self.images if image.variable_axis in (COUCH, REFERENCE)]
        elif axis == COMBO:
            images = [image for image in self.images if image.variable_axis in (COMBO,)]
        elif axis == ALL:
            images = self.images

        # create plots
        max_num_images = math.ceil(len(images)/4)
        dpi = 72
        width_px = 1080
        width_in = width_px/dpi
        height_in = (width_in / 4) * max_num_images
        fig, axes = plt.subplots(nrows=max_num_images, ncols=4, figsize=(width_in, height_in))
        for mpl_axis, wl_image in zip_longest(axes.flatten(), images):
            plot_image(wl_image, mpl_axis)

        # set titles
        fig.suptitle("{} images".format(axis), fontsize=14, y=1)
        plt.tight_layout()
        if show:
            plt.show()

    @value_accept(axis=(GANTRY, COLLIMATOR, COUCH, COMBO, ALL))
    def save_images(self, filename: str, axis: str=ALL, **kwargs):
        """Save the figure of `plot_images()` to file. Keyword arguments are passed to `matplotlib.pyplot.savefig()`.

        Parameters
        ----------
        filename : str
            The name of the file to save to.
        """
        self.plot_images(axis=axis, show=False)
        plt.savefig(filename, **kwargs)

    def plot_summary(self, show: bool=True):
        """Plot a summary figure showing the gantry sag and wobble plots of the three axes."""
        plt.figure(figsize=(11, 9))
        grid = (3, 6)
        gantry_sag_ax = plt.subplot2grid(grid, (0, 0), colspan=3)
        self._plot_deviation(GANTRY, gantry_sag_ax, show=False)
        epid_sag_ax = plt.subplot2grid(grid, (0, 3), colspan=3)
        self._plot_deviation(EPID, epid_sag_ax, show=False)
        if self._get_images((COLLIMATOR, REFERENCE))[0] > 1:
            coll_sag_ax = plt.subplot2grid(grid, (1, 0), colspan=3)
            self._plot_deviation(COLLIMATOR, coll_sag_ax, show=False)
        if self._get_images((COUCH, REFERENCE))[0] > 1:
            couch_sag_ax = plt.subplot2grid(grid, (1, 3), colspan=3)
            self._plot_deviation(COUCH, couch_sag_ax, show=False)

        for axis, axnum in zip((GANTRY, COLLIMATOR, COUCH), (0, 2, 4)):
            if self._get_images((axis, REFERENCE))[0] > 1:
                ax = plt.subplot2grid(grid, (2, axnum), colspan=2)
                self.plot_axis_images(axis=axis, ax=ax, show=False)
        if show:
            plt.tight_layout()
            plt.show()

    def save_summary(self, filename: str, **kwargs):
        """Save the summary image."""
        self.plot_summary(show=False)
        plt.tight_layout()
        plt.savefig(filename, **kwargs)

    def results(self) -> str:
        """Return the analysis results summary."""
        result = "\nWinston-Lutz Analysis\n" \
                 "=================================\n" \
                 "Number of images: {}\n" \
                 "Maximum 2D CAX->BB distance: {:.2f}mm\n" \
                 "Median 2D CAX->BB distance: {:.2f}mm\n" \
                 "Shift BB to iso, facing gantry: {}\n" \
                 "Gantry 3D isocenter diameter: {:.2f}mm\n" \
                 "Maximum Gantry RMS deviation (mm): {:.2f}mm\n" \
                 "Maximum EPID RMS deviation (mm): {:.2f}mm\n" \
                 "Collimator 2D isocenter diameter: {:.2f}mm\n" \
                 "Maximum Collimator RMS deviation (mm): {:.2f}\n" \
                 "Couch 2D isocenter diameter: {:.2f}mm\n" \
                 "Maximum Couch RMS deviation (mm): {:.2f}".format(
                    len(self.images), self.cax2bb_distance('max'), self.cax2bb_distance('median'),
                    self.bb_shift_instructions(), self.gantry_iso_size, max(self.axis_rms_deviation(GANTRY)),
                    max(self.axis_rms_deviation(EPID)),
                    self.collimator_iso_size, max(self.axis_rms_deviation(COLLIMATOR)),
                    self.couch_iso_size,
                    max(self.axis_rms_deviation(COUCH)),
                 )

        return result

    def publish_pdf(self, filename: str, notes: Union[str, List[str]]=None, open_file: bool=False, metadata: dict=None):
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
        plt.ioff()
        title = "Winston-Lutz Analysis"
        canvas = pdf.PylinacCanvas(filename, page_title=title, metadata=metadata)
        avg_sid = np.mean([image.metadata.RTImageSID for image in self.images])
        text = ['Winston-Lutz results:',
                'Average SID (mm): {:2.0f}'.format(avg_sid),
                'Number of images: {}'.format(len(self.images)),
                'Maximum distance to BB (mm): {:2.2f}'.format(self.cax2bb_distance('max')),
                'Median distances to BB (mm): {:2.2f}'.format(self.cax2bb_distance('median')),
                'Gantry 3D isocenter diameter (mm): {:2.2f}'.format(self.gantry_iso_size),
                ]
        if self._contains_axis_images(COLLIMATOR):
            text.append('Collimator 2D isocenter diameter (mm): {:2.2f}'.format(self.collimator_iso_size),)
        if self._contains_axis_images(COUCH):
            text.append('Couch 2D isocenter diameter (mm): {:2.2f}'.format(self.couch_iso_size), )
        canvas.add_text(text=text, location=(10, 25.5))
        # draw summary image on 1st page
        data = io.BytesIO()
        self.save_summary(data, figsize=(10, 10))
        canvas.add_image(image_data=data, location=(2, 3), dimensions=(18, 18))
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 4.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 4))
        # add more pages showing individual axis images
        for ax in (GANTRY, COLLIMATOR, COUCH, COMBO):
            if self._contains_axis_images(ax):
                canvas.add_new_page(metadata=metadata)
                data = io.BytesIO()
                self.save_images(data, axis=ax, figsize=(10, 10))
                canvas.add_image(data, location=(2, 7), dimensions=(18, 18))

        canvas.finish()

        if open_file:
            open_path(filename)

    def _contains_axis_images(self, axis: str=GANTRY) -> bool:
        """Return whether or not the set of WL images contains images pertaining to a given axis"""
        return any(True for image in self.images if image.variable_axis in (axis,))


class WLImage(image.LinacDicomImage):
    """Holds individual Winston-Lutz EPID images, image properties, and automatically finds the field CAX and BB."""

    def __init__(self, file: str, use_filenames: bool):
        """
        Parameters
        ----------
        file : str
            Path to the image file.
        use_filenames: bool
            Whether to try to use the file name to determine axis values.
            Useful for Elekta machines that do not include that info in the DICOM data.
        """
        super().__init__(file, use_filenames=use_filenames)
        self.file = osp.basename(file)
        self.check_inversion()
        self.flipud()
        self._clean_edges()
        self.field_cax, self.bounding_box = self._find_field_centroid()
        self.bb = self._find_bb()

    def __repr__(self):
        return "WLImage(G={0:.1f}, B={1:.1f}, P={2:.1f})".format(self.gantry_angle, self.collimator_angle, self.couch_angle)

    def _clean_edges(self, window_size: int=2):
        """Clean the edges of the image to be near the background level."""
        def has_noise(self, window_size):
            """Helper method to determine if there is spurious signal at any of the image edges.

            Determines if the min or max of an edge is within 10% of the baseline value and trims if not.
            """
            near_min, near_max = np.percentile(self.array, [5, 99.5])
            img_range = near_max - near_min
            top = self[:window_size, :]
            left = self[:, :window_size]
            bottom = self[-window_size:, :]
            right = self[:, -window_size:]
            edge_array = np.concatenate((top.flatten(), left.flatten(), bottom.flatten(), right.flatten()))
            edge_too_low = edge_array.min() < (near_min - img_range / 10)
            edge_too_high = edge_array.max() > (near_max + img_range / 10)
            return edge_too_low or edge_too_high

        safety_stop = np.min(self.shape)/10
        while has_noise(self, window_size) and safety_stop > 0:
            self.remove_edges(window_size)
            safety_stop -= 1

    def _find_field_centroid(self) -> Tuple[Point, List]:
        """Find the centroid of the radiation field based on a 50% height threshold.

        Returns
        -------
        p
            The CAX point location.
        edges
            The bounding box of the field, plus a small margin.
        """
        min, max = np.percentile(self.array, [5, 99.9])
        threshold_img = self.as_binary((max - min)/2 + min)
        # clean single-pixel noise from outside field
        cleaned_img = ndimage.binary_erosion(threshold_img)
        [*edges] = bounding_box(cleaned_img)
        edges[0] -= 10
        edges[1] += 10
        edges[2] -= 10
        edges[3] += 10
        coords = ndimage.measurements.center_of_mass(threshold_img)
        p = Point(x=coords[-1], y=coords[0])
        return p, edges

    def _find_bb(self) -> Point:
        """Find the BB within the radiation field. Iteratively searches for a circle-like object
        by lowering a low-pass threshold value until found.

        Returns
        -------
        Point
            The weighted-pixel value location of the BB.
        """
        # get initial starting conditions
        hmin, hmax = np.percentile(self.array, [5, 99.9])
        spread = hmax - hmin
        max_thresh = hmax
        lower_thresh = hmax - spread / 1.5
        # search for the BB by iteratively lowering the low-pass threshold value until the BB is found.
        found = False
        while not found:
            try:
                binary_arr = np.logical_and((max_thresh > self), (self >= lower_thresh))
                labeled_arr, num_roi = ndimage.measurements.label(binary_arr)
                roi_sizes, bin_edges = np.histogram(labeled_arr, bins=num_roi + 1)
                bw_bb_img = np.where(labeled_arr == np.argsort(roi_sizes)[-3], 1, 0)

                if not is_round(bw_bb_img):
                    raise ValueError
                if not is_modest_size(bw_bb_img, self.bounding_box):
                    raise ValueError
                if not is_symmetric(bw_bb_img):
                    raise ValueError
            except (IndexError, ValueError):
                max_thresh -= 0.05 * spread
                if max_thresh < hmin:
                    raise ValueError("Unable to locate the BB. Make sure the field edges do not obscure the BB and that there is no artifacts in the images.")
            else:
                found = True

        # determine the center of mass of the BB
        inv_img = image.load(self.array)
        inv_img.invert()
        x_arr = np.abs(np.average(bw_bb_img, weights=inv_img, axis=0))
        x_com = SingleProfile(x_arr).fwxm_center(interpolate=True)
        y_arr = np.abs(np.average(bw_bb_img, weights=inv_img, axis=1))
        y_com = SingleProfile(y_arr).fwxm_center(interpolate=True)
        return Point(x_com, y_com)

    @property
    def epid(self) -> Point:
        """Center of the EPID panel"""
        return self.center

    @property
    def epid_y_offset(self) -> float:
        """The offset or distance between the field CAX and EPID in the y-direction (AP)."""
        return -sin(self.gantry_angle) * self.cax2epid_vector.x

    @property
    def bb_y_offset(self) -> float:
        """The offset or distance between the field CAX and BB in the y-direction (AP)."""
        return -sin(self.gantry_angle) * self.cax2bb_vector.x

    @property
    def epid_x_offset(self) -> float:
        """The offset or distance between the field CAX and EPID in the x-direction (LR)."""
        return cos(self.gantry_angle) * self.cax2epid_vector.x

    @property
    def bb_x_offset(self) -> float:
        """The offset or distance between the field CAX and BB in the x-direction (LR)."""
        return cos(self.gantry_angle) * self.cax2bb_vector.x

    @property
    def epid_z_offset(self) -> float:
        """The offset or distance between the field CAX and EPID in z-direction (SI)."""
        if is_close(self.couch_angle, [0, 360], delta=2):
            return -self.cax2epid_vector.y

    @property
    def bb_z_offset(self) -> float:
        """The offset or distance between the field CAX and BB in z-direction (SI)."""
        return -self.cax2bb_vector.y

    @property
    def cax_line_projection(self) -> Line:
        """The projection of the field CAX through space around the area of the BB.
        Used for determining gantry isocenter size.

        Returns
        -------
        Line
            The virtual line in space made by the beam CAX.
        """
        p1 = Point()
        p2 = Point()
        # point 1 - ray origin
        p1.x = self.bb_x_offset + 20 * sin(self.gantry_angle)
        p1.y = self.bb_y_offset + 20 * cos(self.gantry_angle)
        p1.z = self.bb_z_offset
        # point 2 - ray destination
        p2.x = self.bb_x_offset - 20 * sin(self.gantry_angle)
        p2.y = self.bb_y_offset - 20 * cos(self.gantry_angle)
        p2.z = self.bb_z_offset
        l = Line(p1, p2)
        return l

    @property
    def cax2bb_vector(self) -> Vector:
        """The vector in mm from the CAX to the BB."""
        dist = (self.bb - self.field_cax) / self.dpmm
        # translate vector by couch angle
        c_ang = self.couch_angle
        new_x = dist.x * cos(c_ang) - dist.y * sin(c_ang)
        new_y = dist.x * sin(c_ang) + dist.y * cos(c_ang)
        return Vector(new_x, new_y, dist.z)

    @property
    def cax2bb_vector3d(self) -> Vector:
        """The vector in mm from the CAX to the BB."""
        return Vector(self.bb_x_offset, self.bb_y_offset, self.bb_z_offset)

    @property
    def cax2epid_vector(self) -> Vector:
        """The vector in mm from the CAX to the EPID center pixel"""
        dist = (self.epid - self.field_cax) / self.dpmm
        return Vector(dist.x, dist.y, dist.z)

    @property
    def cax2bb_distance(self) -> float:
        """The scalar distance in mm from the CAX to the BB."""
        dist = self.field_cax.distance_to(self.bb)
        return dist / self.dpmm

    @property
    def cax2epid_distance(self) -> float:
        """The scalar distance in mm from the CAX to the EPID center pixel"""
        return self.field_cax.distance_to(self.epid) / self.dpmm

    def plot(self, ax=None, show=True, clear_fig=False):
        """Plot the image, zoomed-in on the radiation field, along with the detected
        BB location and field CAX location.

        Parameters
        ----------
        ax : None, matplotlib Axes instance
            The axis to plot to. If None, will create a new figure.
        show : bool
            Whether to actually show the image.
        clear_fig : bool
            Whether to clear the figure first before drawing.
        """
        ax = super().plot(ax=ax, show=False, clear_fig=clear_fig)
        ax.plot(self.field_cax.x, self.field_cax.y, 'gs', ms=8)
        ax.plot(self.bb.x, self.bb.y, 'ro', ms=8)
        ax.plot(self.epid.x, self.epid.y, 'b+', ms=8)
        ax.set_ylim([self.bounding_box[0], self.bounding_box[1]])
        ax.set_xlim([self.bounding_box[2], self.bounding_box[3]])
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_title('\n'.join(wrap(self.file, 30)), fontsize=10)
        ax.set_xlabel("G={0:.0f}, B={1:.0f}, P={2:.0f}".format(self.gantry_angle, self.collimator_angle, self.couch_angle))
        ax.set_ylabel("CAX to BB: {0:3.2f}mm".format(self.cax2bb_distance))
        if show:
            plt.show()
        return ax

    def save_plot(self, filename: str, **kwargs):
        """Save the image plot to file."""
        self.plot(show=False)
        plt.savefig(filename, **kwargs)

    @property
    def variable_axis(self) -> str:
        """The axis that is varying.

        There are five types of images:

        * Reference : All axes are at 0.
        * Gantry: All axes but gantry at 0.
        * Collimator : All axes but collimator at 0.
        * Couch : All axes but couch at 0.
        * Combo : More than one axis is not at 0.
        """
        G0 = is_close(self.gantry_angle, [0, 360])
        B0 = is_close(self.collimator_angle, [0, 360])
        P0 = is_close(self.couch_angle, [0, 360])
        if G0 and B0 and not P0:
            return COUCH
        elif G0 and P0 and not B0:
            return COLLIMATOR
        elif P0 and B0 and not G0:
            return GANTRY
        elif P0 and B0 and G0:
            return REFERENCE
        else:
            return COMBO


def is_symmetric(logical_array: np.ndarray) -> bool:
    """Whether the binary object's dimensions are symmetric, i.e. a perfect circle. Used to find the BB."""
    ymin, ymax, xmin, xmax = bounding_box(logical_array)
    y = abs(ymax - ymin)
    x = abs(xmax - xmin)
    if x > max(y * 1.05, y + 3) or x < min(y * 0.95, y - 3):
        return False
    return True


def is_modest_size(logical_array: np.ndarray, field_bounding_box) -> bool:
    """Decide whether the ROI is roughly the size of a BB; not noise and not an artifact. Used to find the BB."""
    bbox = field_bounding_box
    rad_field_area = (bbox[1] - bbox[0]) * (bbox[3] - bbox[2])
    return rad_field_area * 0.003 < np.sum(logical_array) < rad_field_area * 0.3


def is_round(logical_array: np.ndarray) -> bool:
    """Decide if the ROI is circular in nature by testing the filled area vs bounding box. Used to find the BB."""
    expected_fill_ratio = np.pi / 4
    actual_fill_ratio = filled_area_ratio(logical_array)
    return expected_fill_ratio * 1.2 > actual_fill_ratio > expected_fill_ratio * 0.8
