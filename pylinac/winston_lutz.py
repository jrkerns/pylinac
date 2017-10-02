"""The Winston-Lutz module loads and processes EPID images that have acquired Winston-Lutz type images.

Features:

* **Automatic field & BB positioning** - When an image or directory is loaded, the field CAX and the BB
  are automatically found, along with the vector and scalar distance between them.
* **Isocenter size determination** - Using backprojections of the EPID images, the 3D gantry isocenter size
  and position can be determined *independent of the BB position*. Additionally, the 2D planar isocenter size
  of the collimator and couch can also be determined.
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
import re

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage, optimize

from .core import image
from .core.decorators import value_accept
from .core.geometry import Point, Line, Circle, Vector, cos, sin
from .core.io import TemporaryZipDirectory, get_url, retrieve_demo_file
from .core.mask import filled_area_ratio, bounding_box
from .core import pdf
from .core.profile import SingleProfile
from .core.utilities import is_close, is_dicom_image

GANTRY = 'Gantry'
COLLIMATOR = 'Collimator'
COUCH = 'Couch'
COMBO = 'Combo'
EPID = 'Epid'
REFERENCE = 'Reference'
ALL = 'All'


class WinstonLutz:
    """Class for performing a Winston-Lutz test of the radiation isocenter."""

    def __init__(self, directory):
        """
        Parameters
        ----------
        directory : str
            Path to the directory of the Winston-Lutz EPID images.

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
        self.images = ImageManager(directory)

    @classmethod
    def from_demo_images(cls):
        """Instantiate using the demo images."""
        demo_file = retrieve_demo_file(url='winston_lutz.zip')
        return cls.from_zip(demo_file)

    @classmethod
    def from_zip(cls, zfile):
        """Instantiate from a zip file rather than a directory.

        Parameters
        ----------
        zfile : str
            Path to the archive file.
        """
        with TemporaryZipDirectory(zfile) as tmpz:
            obj = cls(tmpz)
        return obj

    @classmethod
    def from_url(cls, url):
        """Instantiate from a URL.

        Parameters
        ----------
        url : str
            URL that points to a zip archive of the DICOM images.
        """
        zfile = get_url(url)
        return cls.from_zip(zfile)

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
    def gantry_iso_size(self):
        """The diameter of the 3D gantry isocenter size in mm. Only images where the collimator
        and couch were at 0 are used to determine this value."""
        return self._minimize_axis(GANTRY).fun * 2

    @property
    def gantry_iso2bb_vector(self):
        """The 3D vector from the isocenter to the BB (located at the origin)."""
        min_fun = self._minimize_axis(GANTRY)
        # optimization result goes from origin TO the iso, thus to go from iso to bb (origin), invert sign
        return Vector(min_fun.x[0], min_fun.x[1], min_fun.x[2])

    @property
    def collimator_iso_size(self):
        """The 2D collimator isocenter size (diameter) in mm. The iso size is in the plane
        normal to the gantry."""
        if self._get_images((COLLIMATOR, REFERENCE))[0] > 1:
            return self._minimize_axis(COLLIMATOR).fun * 2
        else:
            return 0

    @property
    def collimator_iso2bb_vector(self):
        """The 2D vector from the collimator isocenter to the BB (located at the origin)."""
        if self._get_images((COLLIMATOR, REFERENCE))[0] > 1:
            min_col = self._minimize_axis(COLLIMATOR)
            return Vector(min_col.x[0], min_col.x[1])
        else:
            return Vector()

    @property
    def couch_iso_size(self):
        """The diameter of the 2D couch isocenter size in mm. Only images where
        the gantry and collimator were at zero are used to determine this value."""
        if self._get_images((COUCH, REFERENCE))[0] > 1:
            return self._minimize_axis(COUCH).x[3] * 2
        else:
            return 0

    @property
    def couch_iso2bb_vector(self):
        """The 2D vector from the couch isocenter to the BB (located at the origin)."""
        if self._get_images((COUCH, REFERENCE))[0] > 1:
            min_col = self._minimize_axis(COUCH)
            return Vector(min_col.x[0], min_col.x[1])
        else:
            return Vector()

    @value_accept(axis=(GANTRY, COLLIMATOR, COUCH, EPID, COMBO), value=('all', 'range'))
    def axis_rms_deviation(self, axis=GANTRY, value='all'):
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

    def cax2bb_distance(self, metric='max'):
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

    def cax2epid_distance(self, metric='max'):
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
    def _plot_deviation(self, item, ax=None, show=True):
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
        ax.grid('on')
        ax.legend(numpoints=1)
        if show:
            plt.show()

    def _get_images(self, axis=(GANTRY,)):
        if isinstance(axis, str):
            axis = (axis,)
        images = [image for image in self.images if image.variable_axis in axis]
        return len(images), images

    @value_accept(axis=(GANTRY, COLLIMATOR, COUCH, COMBO))
    def plot_axis_images(self, axis=GANTRY, show=True, ax=None):
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
    def plot_images(self, axis=ALL, show=True):
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
        fig, axes = plt.subplots(nrows=max_num_images, ncols=4)
        for mpl_axis, wl_image in zip_longest(axes.flatten(), images):
            plot_image(wl_image, mpl_axis)

        # set titles
        fig.suptitle("{} images".format(axis), fontsize=14, y=1)
        plt.tight_layout()
        if show:
            plt.show()

    @value_accept(axis=(GANTRY, COLLIMATOR, COUCH, COMBO, ALL))
    def save_images(self, filename, axis=ALL, **kwargs):
        """Save the figure of `plot_images()` to file. Keyword arguments are passed to `matplotlib.pyplot.savefig()`.

        Parameters
        ----------
        filename : str
            The name of the file to save to.
        """
        self.plot_images(axis=axis, show=False)
        plt.savefig(filename, **kwargs)

    def plot_summary(self, show=True):
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

    def save_summary(self, filename, **kwargs):
        """Save the summary image."""
        self.plot_summary(show=False)
        plt.tight_layout()
        plt.savefig(filename, **kwargs)

    def results(self):
        """Return the analysis results summary."""
        result = "\nWinston-Lutz Analysis\n\n" \
                 "Number of images: {}\n" \
                 "Maximum 2D CAX->BB distance: {:.2f}mm\n" \
                 "Median 2D CAX->BB distance: {:.2f}mm\n" \
                 "Gantry 3D isocenter diameter: {:.2f}mm\n" \
                 "Gantry iso->BB vector: {}\n" \
                 "Maximum Gantry RMS deviation (mm): {:.2f}mm\n" \
                 "Maximum EPID RMS deviation (mm): {:.2f}mm\n" \
                 "Collimator 2D isocenter diameter: {:.2f}mm\n" \
                 "Collimator 2D iso->BB vector: {}\n" \
                 "Maximum Collimator RMS deviation (mm): {:.2f}\n" \
                 "Couch 2D isocenter diameter: {:.2f}mm\n" \
                 "Couch 2D iso->BB vector: {}\n" \
                 "Maximum Couch RMS deviation (mm): {:.2f}".format(
                    len(self.images), self.cax2bb_distance('max'), self.cax2bb_distance('median'),
                    self.gantry_iso_size, self.gantry_iso2bb_vector, max(self.axis_rms_deviation(GANTRY)),
                    max(self.axis_rms_deviation(EPID)),
                    self.collimator_iso_size, self.collimator_iso2bb_vector, max(self.axis_rms_deviation(COLLIMATOR)),
                    self.couch_iso_size, self.couch_iso2bb_vector,
                    max(self.axis_rms_deviation(COUCH)),
                 )

        return result

    def publish_pdf(self, filename, unit=None, notes=None, open_file=False):
        """Publish (print) a PDF containing the analysis and quantitative results.

        Parameters
        ----------
        filename : (str, file-like object}
            The file to write the results to.
        """
        from reportlab.lib.units import cm
        title = "Winston-Lutz Analysis"
        canvas = pdf.create_pylinac_page_template(filename, analysis_title=title, unit=unit)
        avg_sid = np.mean([image.metadata.RTImageSID for image in self.images])
        text = ['Winston-Lutz results:',
                'Key, looking from foot of table:',
                '+x: right, +y: up, +z:out',
                'Average SID (mm): {:2.0f}'.format(avg_sid),
                'Number of images: {}'.format(len(self.images)),
                'Maximum 2D CAX->BB distance (mm): {:2.2f}'.format(self.cax2bb_distance('max')),
                'Median 2D CAX->BB distance (mm): {:2.2f}'.format(self.cax2bb_distance('median')),
                'Gantry iso->BB vector (mm): x={:2.2f}, y={:2.2f}, z={:2.2f}'.format(self.gantry_iso2bb_vector.x,
                                                                 self.gantry_iso2bb_vector.y,
                                                                 self.gantry_iso2bb_vector.z),
                'Gantry 3D isocenter diameter (mm): {:2.2f}'.format(self.gantry_iso_size),
                ]
        if self._contains_axis_images(COLLIMATOR):
            text.append('Collimator 2D isocenter diameter (mm): {:2.2f}'.format(self.collimator_iso_size),)
        if self._contains_axis_images(COUCH):
            text.append('Couch 2D isocenter diameter (mm): {:2.2f}'.format(self.couch_iso_size), )
        pdf.draw_text(canvas, x=10*cm, y=25.5*cm, text=text)
        # draw summary image on 1st page
        data = io.BytesIO()
        self.save_summary(data, figsize=(10, 10))
        img = pdf.create_stream_image(data)
        canvas.drawImage(img, 2 * cm, 3 * cm, width=18 * cm, height=18 * cm, preserveAspectRatio=True)
        if notes is not None:
            pdf.draw_text(canvas, x=1*cm, y=4.5*cm, fontsize=14, text="Notes:")
            pdf.draw_text(canvas, x=1*cm, y=4*cm, text=notes)
        canvas.showPage()
        # add more pages showing individual axis images
        for ax in (GANTRY, COLLIMATOR, COUCH, COMBO):
            if self._contains_axis_images(ax):
                pdf.add_pylinac_page_template(canvas, analysis_title=title)
                data = io.BytesIO()
                self.save_images(data, axis=ax, figsize=(10, 10))
                img = pdf.create_stream_image(data)
                canvas.drawImage(img, 2*cm, 7*cm, width=18*cm, height=18*cm, preserveAspectRatio=True)
                canvas.showPage()
        pdf.finish(canvas, open_file=open_file, filename=filename)

    def _contains_axis_images(self, axis=GANTRY):
        """Return whether or not the set of WL images contains images pertaining to a given axis"""
        return any(True for image in self.images if image.variable_axis in (axis,))


class ImageManager(list):
    """Manages the images of a Winston-Lutz test."""
    def __init__(self, directory):
        """
        Parameters
        ----------
        directory : str
            The path to the images.
        """
        super().__init__()
        if isinstance(directory, list):
            for file in directory:
                if is_dicom_image(file):
                    img = WLImage(file)
                    self.append(img)
        elif not osp.isdir(directory):
            raise ValueError("Invalid directory passed. Check the correct method and file was used.")
        else:
            image_files = image.retrieve_image_files(directory)
            for file in image_files:
                img = WLImage(file)
                self.append(img)
        # reorder list based on increasing gantry angle
        self.sort(key=lambda i: (i.gantry_angle, i.collimator_angle, i.couch_angle))


class WLImage(image.DicomImage):
    """Holds individual Winston-Lutz EPID images, image properties, and automatically finds the field CAX and BB."""

    def __init__(self, file):
        """
        Parameters
        ----------
        file : str
            Path to the image file.
        """
        super().__init__(file)
        self.check_inversion()
        self.flipud()
        self._clean_edges()
        self.field_cax, self.bounding_box = self._find_field_centroid()
        self.bb = self._find_bb()

    def __repr__(self):
        return "WLImage(G={0:.1f}, B={1:.1f}, P={2:.1f})".format(self.gantry_angle, self.collimator_angle, self.couch_angle)

    def _clean_edges(self, window_size=2):
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

    def _find_field_centroid(self):
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

    def _find_bb(self):
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
    def epid(self):
        """Center of the EPID panel"""
        return self.center

    @property
    @lru_cache(1)
    def gantry_angle(self):
        """Gantry angle of the irradiation."""
        return self.get_axis(GANTRY.lower(), 'GantryAngle')

    @property
    @lru_cache(1)
    def collimator_angle(self):
        """Collimator angle of the irradiation."""
        return self.get_axis('coll', 'BeamLimitingDeviceAngle')

    @property
    @lru_cache(1)
    def couch_angle(self):
        """Couch angle of the irradiation."""
        return self.get_axis(COUCH.lower(), 'PatientSupportAngle')

    def get_axis(self, axis_str, axis_dcm_attr):
        """Retrieve the value of the axis. This will first look in the file name for the value. 
        If not in the filename then it will look in the DICOM metadata. If the value can be found in neither 
        then a value of 0 is assumed.

        Parameters
        ----------
        axis_str : str
            The string to look for in the filename.
        axis_dcm_attr : str
            The DICOM attribute that should contain the axis value.

        Returns
        -------
        float
        """
        filename = osp.basename(self.path)
        if axis_str in filename:
            try:
                match = re.search('(?<={})\d+'.format(axis_str), filename)
                axis = float(match.group())
            except:
                raise ValueError(
                    "The filename contains '{}' but could not read a number following it. Use the format '{}45' e.g.".format(
                        axis_str, axis_str))
        else:
            try:
                axis = float(getattr(self.metadata, axis_dcm_attr))
            except AttributeError:
                axis = 0
        if is_close(axis, [0, 360], delta=1):
            return 0
        else:
            return axis

    @property
    def epid_y_offset(self):
        """The offset or distance between the field CAX and EPID in the y-direction (AP)."""
        return -sin(self.gantry_angle) * self.cax2epid_vector.x

    @property
    def bb_y_offset(self):
        """The offset or distance between the field CAX and BB in the y-direction (AP)."""
        return -sin(self.gantry_angle) * self.cax2bb_vector.x

    @property
    def epid_x_offset(self):
        """The offset or distance between the field CAX and EPID in the x-direction (LR)."""
        return cos(self.gantry_angle) * self.cax2epid_vector.x

    @property
    def bb_x_offset(self):
        """The offset or distance between the field CAX and BB in the x-direction (LR)."""
        return cos(self.gantry_angle) * self.cax2bb_vector.x

    @property
    def epid_z_offset(self):
        """The offset or distance between the field CAX and EPID in z-direction (SI)."""
        if is_close(self.couch_angle, [0, 360], delta=2):
            return -self.cax2epid_vector.y

    @property
    def bb_z_offset(self):
        """The offset or distance between the field CAX and BB in z-direction (SI)."""
        return -self.cax2bb_vector.y

    @property
    def cax_line_projection(self):
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
    def cax2bb_vector(self):
        """The vector in mm from the CAX to the BB."""
        dist = (self.bb - self.field_cax) / self.dpmm
        return Vector(dist.x, dist.y, dist.z)

    @property
    def cax2epid_vector(self):
        """The vector in mm from the CAX to the EPID center pixel"""
        dist = (self.epid - self.field_cax) / self.dpmm
        return Vector(dist.x, dist.y, dist.z)

    @property
    def cax2bb_distance(self):
        """The scalar distance in mm from the CAX to the BB."""
        dist = self.field_cax.distance_to(self.bb)
        return dist / self.dpmm

    @property
    def cax2epid_distance(self):
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
        ax.set_xlabel("G={0:.0f}, B={1:.0f}, P={2:.0f}".format(self.gantry_angle, self.collimator_angle, self.couch_angle))
        ax.set_ylabel("CAX to BB: {0:3.2f}mm".format(self.cax2bb_distance))
        if show:
            plt.show()
        return ax

    def save_plot(self, filename, **kwargs):
        """Save the image plot to file."""
        self.plot(show=False)
        plt.savefig(filename, **kwargs)

    @property
    def variable_axis(self):
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


def is_symmetric(logical_array):
    """Whether the binary object's dimensions are symmetric, i.e. a perfect circle. Used to find the BB."""
    ymin, ymax, xmin, xmax = bounding_box(logical_array)
    y = abs(ymax - ymin)
    x = abs(xmax - xmin)
    if x > max(y * 1.05, y + 3) or x < min(y * 0.95, y - 3):
        return False
    return True


def is_modest_size(logical_array, field_bounding_box):
    """Decide whether the ROI is roughly the size of a BB; not noise and not an artifact. Used to find the BB."""
    bbox = field_bounding_box
    rad_field_area = (bbox[1] - bbox[0]) * (bbox[3] - bbox[2])
    return rad_field_area * 0.003 < np.sum(logical_array) < rad_field_area * 0.3


def is_round(logical_array):
    """Decide if the ROI is circular in nature by testing the filled area vs bounding box. Used to find the BB."""
    expected_fill_ratio = np.pi / 4
    actual_fill_ratio = filled_area_ratio(logical_array)
    return expected_fill_ratio * 1.2 > actual_fill_ratio > expected_fill_ratio * 0.8
