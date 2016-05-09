"""The Winston-Lutz module loads and processes EPID images that have acquired Winston-Lutz type images.

Features:

* **Automatic field & BB positioning** - When an image or directory is loaded, the field CAX and the BB
  are automatically found, along with the vector and scalar distance between them.
* **Isocenter size determination** - Using backprojections of the EPID images, the 3D gantry isocenter size
  and position can be determined *independent of the BB position*. Additionally, the 2D planar isocenter size
  of the collimator and couch can also be determined.
* **Image plotting** - WL images can be plotted separately or together, each of which shows the field CAX, BB and
  scalar distance from BB to CAX.
* **Gantry sag** - The sag of the gantry is also quantified and can be plotted.
"""
from functools import lru_cache
from itertools import zip_longest
import os.path as osp

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage, optimize

from .core import image
from .core.geometry import Point, Line, Circle, Vector, cos, sin
from .core.io import TemporaryZipDirectory, get_url
from .core.mask import filled_area_ratio, bounding_box
from .core.profile import SingleProfile
from .core.utilities import is_close, retrieve_demo_file

GANTRY = 'Gantry'
COLLIMATOR = 'Collimator'
COUCH = 'Couch'
COMBO = 'Combo'
REFERENCE = 'Reference'


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
    def _minimize_axis(self, axis='Gantry'):
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
        return self._minimize_axis(COLLIMATOR).fun * 2

    @property
    def collimator_iso2bb_vector(self):
        """The 2D vector from the collimator isocenter to the BB (located at the origin)."""
        min_col = self._minimize_axis(COLLIMATOR)
        return Vector(min_col.x[0], min_col.x[1])

    @property
    def couch_iso_size(self):
        """The diameter of the 2D couch isocenter size in mm. Only images where
        the gantry and collimator were at zero are used to determine this value."""
        return self._minimize_axis(COUCH).x[3] * 2

    @property
    def couch_iso2bb_vector(self):
        """The 2D vector from the couch isocenter to the BB (located at the origin)."""
        min_col = self._minimize_axis(COUCH)
        return Vector(min_col.x[0], min_col.x[1])

    def gantry_sag(self, axis='z'):
        """The range of gantry sag along the given axis."""
        attr = axis + '_offset'
        sag_array = [getattr(image, attr) for image in self.images if image.variable_axis in (GANTRY, REFERENCE)]
        return max(sag_array) - min(sag_array)

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

    def plot_gantry_sag(self, ax=None, show=True):
        """Plot the gantry sag in Cartesian coordinates.

        Parameters
        ----------
        ax : None, matplotlib.Axes
            The axis to plot to. If None, creates a new plot.
        show : bool
            Whether to show the image.
        """
        # get gantry images, angles, and shifts
        gantry_imgs = [image for image in self.images if image.variable_axis in (GANTRY, REFERENCE)]
        gantry_angles = [image.gantry_angle for image in gantry_imgs]
        gantry_zsag = np.array([image.z_offset for image in gantry_imgs])
        gantry_ysag = np.array([image.y_offset for image in gantry_imgs])
        gantry_xsag = np.array([image.x_offset for image in gantry_imgs])

        # plot the gantry sag
        if ax is None:
            ax = plt.subplot(111)
        ax.plot(gantry_angles, gantry_zsag, 'bo', label='In/Out')
        ax.plot(gantry_angles, gantry_ysag, 'r*', label='Up/Down')
        ax.plot(gantry_angles, gantry_xsag, 'm^', label='Left/Right')
        ax.set_title("Relative gantry displacement")
        ax.set_ylabel('mm')
        ax.set_xticks(np.arange(0, 361, 45))
        ax.set_xlim([-15, 375])
        ax.grid('on')
        ax.legend(numpoints=1)
        if show:
            plt.show()

    def plot_axis_images(self, axis='Gantry', show=True, ax=None):
        """Plot all CAX/BB positions for the images of a given axis.

        For example, axis='Couch' plots a reference image, and all the BB points of the other
        images where the couch was moving.

        Parameters
        ----------
        axis : {'Gantry', 'Collimator', 'Couch'}
            The images/markers from which accelerator axis to plot.
        show : bool
            Whether to actually show the images.
        ax : None, matplotlib.Axes
            The axis to plot to. If None, creates a new plot.
        """
        images = [image for image in self.images if image.variable_axis in (axis, REFERENCE)]
        ax = images[0].plot(show=False, ax=ax)
        if axis != COUCH:
            xs = [img.field_cax.x for img in images[1:]]
            ys = [img.field_cax.y for img in images[1:]]
            marker = 'gs'
        else:
            xs = [img.bb.x for img in images[1:]]
            ys = [img.bb.y for img in images[1:]]
            marker = 'ro'
        ax.plot(xs, ys, marker, ms=8)
        # set labels
        ax.set_title(axis + ' wobble')
        ax.set_xlabel(axis + ' axis images superimposed')
        ax.set_ylabel(axis + " iso size: {0:3.2f}mm".format(getattr(self, axis.lower() + '_iso_size')))
        if show:
            plt.show()

    def plot_images(self, show=True):
        """Plot a grid of all the images acquired.

        Four columns are plotted with the titles showing which axis that column represents.

        Parameters
        ----------
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
        gantry_images = [image for image in self.images if image.variable_axis in (GANTRY, REFERENCE)]
        bld_images = [image for image in self.images if image.variable_axis in (COLLIMATOR, REFERENCE)]
        psd_images = [image for image in self.images if image.variable_axis in (COUCH, REFERENCE)]
        combo_images = [image for image in self.images if image.variable_axis in (COMBO,)]

        # create plots
        max_num_images = max([len(bld_images), len(gantry_images), len(psd_images), len(combo_images)])
        fig, axes = plt.subplots(nrows=max_num_images, ncols=4)
        for axis, column in zip([gantry_images, bld_images, psd_images, combo_images], range(4)):
            for image, axis in zip_longest(axis, axes[:, column]):
                plot_image(image, axis)

        # set titles
        axes[0, 0].set_title("Gantry shots")
        axes[0, 1].set_title("Collimator shots")
        axes[0, 2].set_title("Couch shots")
        axes[0, 3].set_title("Combination shots")
        # plt.tight_layout()
        if show:
            plt.show()

    def save_images(self, filename, **kwargs):
        """Save the figure of `plot_images()` to file. Keyword arguments are passed to `matplotlib.pyplot.savefig()`.

        Parameters
        ----------
        filename : str
            The name of the file to save to.
        """
        self.plot_images(show=False)
        plt.savefig(filename, **kwargs)

    def plot_summary(self, show=True):
        """Plot a summary figure showing the gantry sag and wobble plots of the three axes."""
        fig, ax = plt.subplots(2, 2, figsize=(9, 7))
        self.plot_gantry_sag(ax.flatten()[0], show=False)
        for axis, ax in zip((GANTRY, COLLIMATOR, COUCH), ax.flatten()[1:]):
            self.plot_axis_images(axis=axis, ax=ax, show=False)
        if show:
            plt.show()

    def save_summary(self, filename, **kwargs):
        """Save the summary image."""
        self.plot_summary(show=False)
        plt.savefig(filename, **kwargs)

    def results(self):
        """Return the analysis results summary."""
        result = "\nWinston-Lutz Analysis\n\n" \
                 "Number of images: {}\n" \
                 "Maximum 2D CAX->BB distance: {:.2f}mm\n" \
                 "Median 2D CAX->BB distance: {:.2f}mm\n" \
                 "Gantry 3D isocenter size: {:.2f}mm\n" \
                 "Gantry iso->BB vector: {}\n" \
                 "Gantry sag in the z-direction: {:.2f}mm\n" \
                 "Collimator 2D isocenter size: {:.2f}mm\n" \
                 "Collimator 2D iso->BB vector: {}\n" \
                 "Couch 2D isocenter size: {:.2f}mm\n" \
                 "Couch 2D iso->BB vector: {}".format(
                    len(self.images), self.cax2bb_distance('max'), self.cax2bb_distance('median'),
                    self.gantry_iso_size, self.gantry_iso2bb_vector, self.gantry_sag(),
                    self.collimator_iso_size, self.collimator_iso2bb_vector,
                    self.couch_iso_size, self.couch_iso2bb_vector
                 )
        return result


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
        if not osp.isdir(directory):
            raise ValueError("Invalid directory passed. Check the correct method and file was used.")
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
    def gantry_angle(self):
        """Gantry angle of the irradiation."""
        if is_close(self.metadata.GantryAngle, [0, 360], delta=1):
            return 0
        else:
            return self.metadata.GantryAngle

    @property
    def collimator_angle(self):
        """Collimator angle of the irradiation."""
        if is_close(self.metadata.BeamLimitingDeviceAngle, [0, 360], delta=1):
            return 0
        else:
            return self.metadata.BeamLimitingDeviceAngle

    @property
    def couch_angle(self):
        """Couch angle of the irradiation."""
        if is_close(self.metadata.PatientSupportAngle, [0, 360], delta=1):
            return 0
        else:
            return self.metadata.PatientSupportAngle

    @property
    def y_offset(self):
        """The offset or distance between the field CAX and BB in the y-direction (AP)."""
        return -sin(self.gantry_angle) * self.cax2bb_vector.x

    @property
    def x_offset(self):
        """The offset or distance between the field CAX and BB in the x-direction (LR)."""
        return cos(self.gantry_angle) * self.cax2bb_vector.x

    @property
    def z_offset(self):
        """The offset or distance between the field CAX and BB in z-direction (SI)."""
        if is_close(self.couch_angle, [0, 360], delta=2):
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
        p1.x = self.x_offset + 20 * sin(self.gantry_angle)
        p1.y = self.y_offset + 20 * cos(self.gantry_angle)
        p1.z = self.z_offset
        # point 2 - ray destination
        p2.x = self.x_offset - 20 * sin(self.gantry_angle)
        p2.y = self.y_offset - 20 * cos(self.gantry_angle)
        p2.z = self.z_offset
        l = Line(p1, p2)
        return l

    @property
    def cax2bb_vector(self):
        """The vector in mm from the CAX to the BB."""
        dist = (self.bb - self.field_cax) / self.dpmm
        return Vector(dist.x, dist.y, dist.z)

    @property
    def cax2bb_distance(self):
        """The scalar distance in mm from the CAX to the BB."""
        dist = self.field_cax.distance_to(self.bb)
        return dist / self.dpmm

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
