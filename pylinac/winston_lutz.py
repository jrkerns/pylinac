"""The Winston-Lutz module loads and processes EPID images that have acquired Winston-Lutz type images.

Features:

* **Couch shift instructions** - After running a WL test, get immediate feedback on how to shift the couch.
  Couch values can also be passed in and the new couch values will be presented so you don't have to do that pesky conversion.
  "Do I subtract that number or add it?"
* **Automatic field & BB positioning** - When an image or directory is loaded, the field CAX and the BB
  are automatically found, along with the vector and scalar distance between them.
* **Isocenter size determination** - Using backprojections of the EPID images, the 3D gantry isocenter size
  and position can be determined *independent of the BB position*. Additionally, the 2D planar isocenter size
  of the collimator and couch can also be determined.
* **Image plotting** - WL images can be plotted separately or together, each of which shows the field CAX, BB and
  scalar distance from BB to CAX.
* **Axis deviation plots** - Plot the variation of the gantry, collimator, couch, and EPID in each plane
  as well as RMS variation.
* **File name interpretation** - Rename DICOM filenames to include axis information for linacs that don't include
  such information in the DICOM tags. E.g. "myWL_gantry45_coll0_couch315.dcm".
"""
import dataclasses
import enum
import io
import math
import os.path as osp
from dataclasses import dataclass
from itertools import zip_longest
from textwrap import wrap
from typing import Union, List, Tuple, Optional, BinaryIO, Iterable

import argue
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage, optimize, linalg
from skimage import measure
from skimage.measure import label
from skimage.measure._regionprops import RegionProperties
from skimage.morphology import remove_small_objects, remove_small_holes

from .core import image, pdf
from .core.decorators import lru_cache
from .core.geometry import Point, Line, Vector, cos, sin
from .core.image import ArrayImage
from .core.io import TemporaryZipDirectory, get_url, retrieve_demo_file, is_dicom_image
from .core.mask import bounding_box
from .core.utilities import is_close, open_path, ResultBase, convert_to_enum


class Axis(enum.Enum):
    GANTRY = 'Gantry'  #:
    COLLIMATOR = 'Collimator'  #:
    COUCH = 'Couch'  #:
    GB_COMBO = 'GB Combo'  #:
    GBP_COMBO = 'GBP Combo'  #:
    EPID = 'Epid'  #:
    REFERENCE = 'Reference'  #:


@dataclass
class WinstonLutz2DResult(ResultBase):
    variable_axis: str  #:
    cax2epid_vector: Vector  #:
    cax2epid_distance: float  #:
    cax2bb_distance: float  #:
    cax2bb_vector: Vector  #:
    bb_location: Point  #:
    field_cax: Point  #:


@dataclass
class WinstonLutzResult(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""
    num_gantry_images: int  #:
    num_gantry_coll_images: int  #:
    num_coll_images: int  #:
    num_couch_images: int  #:
    num_total_images: int  #:
    max_2d_cax_to_bb_mm: float  #:
    median_2d_cax_to_bb_mm: float  #:
    max_2d_cax_to_epid_mm: float  #:
    median_2d_cax_to_epid_mm: float  #:
    gantry_3d_iso_diameter_mm: float  #:
    max_gantry_rms_deviation_mm: float  #:
    max_epid_rms_deviation_mm: float  #:
    gantry_coll_3d_iso_diameter_mm: float  #:
    coll_2d_iso_diameter_mm: float  #:
    max_coll_rms_deviation_mm: float  #:
    couch_2d_iso_diameter_mm: float  #:
    max_couch_rms_deviation_mm: float  #:


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
                    img = WinstonLutz2D(file, use_filenames)
                    self.append(img)
        elif not osp.isdir(directory):
            raise ValueError("Invalid directory passed. Check the correct method and file was used.")
        else:
            image_files = image.retrieve_image_files(directory)
            for file in image_files:
                img = WinstonLutz2D(file, use_filenames)
                self.append(img)
        if len(self) < 2:
            raise ValueError("<2 valid WL images were found in the folder/file. Ensure you chose the correct folder/file for analysis")
        # reorder list based on increasing gantry angle, collimator angle, then couch angle
        self.sort(key=lambda i: (i.gantry_angle, i.collimator_angle, i.couch_angle))


class WinstonLutz:
    """Class for performing a Winston-Lutz test of the radiation isocenter."""
    images: ImageManager  #:
    _find_field: bool

    def __init__(self, directory: str, use_filenames: bool = False):
        """
        Parameters
        ----------
        directory : str
            Path to the directory of the Winston-Lutz EPID images.
        use_filenames: bool
            Whether to try to use the file name to determine axis values.
            Useful for Elekta machines that do not include that info in the DICOM data.
        """
        self.images = ImageManager(directory, use_filenames)
        self._is_analyzed = False

    @classmethod
    def from_demo_images(cls):
        """Instantiate using the demo images."""
        demo_file = retrieve_demo_file(url='winston_lutz.zip')
        return cls.from_zip(demo_file)

    @classmethod
    def from_zip(cls, zfile: Union[str, BinaryIO], use_filenames: bool = False):
        """Instantiate from a zip file rather than a directory.

        Parameters
        ----------
        zfile
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
        wl.analyze()
        print(wl.results())
        wl.plot_summary()

    def analyze(self, bb_size_mm: float = 5, find_field: bool = True):
        """Analyze the WL images.

        Parameters
        ----------
        bb_size_mm
            The expected size of the BB in mm. The actual size of the BB can be +/-2mm from the passed value.
        find_field
            Whether to find the radiation field. Set this to False if performing a fully-open field such as a kV WL.
        """
        self._find_field = find_field
        for img in self.images:
            img.analyze(bb_size_mm, find_field=find_field)
        self._is_analyzed = True

    @lru_cache()
    def _minimize_axis(self, axes: Union[Axis, Tuple[Axis, ...]] = (Axis.GANTRY,)):
        """Return the minimization result of the given axis."""
        if isinstance(axes, Axis):
            axes = (axes,)

        def max_distance_to_lines(p, lines) -> float:
            """Calculate the maximum distance to any line from the given point."""
            point = Point(p[0], p[1], p[2])
            return max(line.distance_to(point) for line in lines)

        things = [image.cax_line_projection for image in self.images if image.variable_axis in (axes + (Axis.REFERENCE,))]
        if len(things) <= 1:
            raise ValueError("Not enough images of the given type to identify the axis isocenter")
        initial_guess = np.array([0, 0, 0])
        bounds = [(-20, 20), (-20, 20), (-20, 20)]
        result = optimize.minimize(max_distance_to_lines, initial_guess, args=things, bounds=bounds)
        return result

    @property
    def gantry_iso_size(self) -> float:
        """The diameter of the 3D gantry isocenter size in mm. Only images where the collimator
        and couch were at 0 are used to determine this value."""
        num_gantry_like_images = self._get_images((Axis.GANTRY, Axis.REFERENCE))[0]
        if num_gantry_like_images > 1:
            return self._minimize_axis(Axis.GANTRY).fun * 2
        else:
            return 0

    @property
    def gantry_coll_iso_size(self) -> float:
        """The diameter of the 3D gantry isocenter size in mm *including collimator and gantry/coll combo images*.
        Images where the couch!=0 are excluded."""
        num_gantry_like_images = self._get_images((Axis.GANTRY, Axis.COLLIMATOR, Axis.GB_COMBO, Axis.REFERENCE))[0]
        if num_gantry_like_images > 1:
            return self._minimize_axis((Axis.GANTRY, Axis.COLLIMATOR, Axis.GB_COMBO)).fun * 2
        else:
            return 0

    @staticmethod
    def _find_max_distance_between_points(images) -> float:
        """Find the maximum distance between a set of points. Used for 2D images like collimator and couch."""
        points = [Point(image.cax2bb_vector.x, image.cax2bb_vector.y) for image in images]
        dists = []
        for point1 in points:
            for point2 in points:
                p = point1.distance_to(point2)
                dists.append(p)
        return max(dists)

    @property
    def collimator_iso_size(self) -> float:
        """The 2D collimator isocenter size (diameter) in mm. The iso size is in the plane
        normal to the gantry."""
        num_collimator_like_images, images = self._get_images((Axis.COLLIMATOR, Axis.REFERENCE))
        if num_collimator_like_images > 1:
            return self._find_max_distance_between_points(images)
        else:
            return 0

    @property
    def couch_iso_size(self) -> float:
        """The diameter of the 2D couch isocenter size in mm. Only images where
        the gantry and collimator were at zero are used to determine this value."""
        num_couch_like_images, images = self._get_images((Axis.COUCH, Axis.REFERENCE))
        if num_couch_like_images > 1:
            return self._find_max_distance_between_points(images)
        else:
            return 0

    @property
    def bb_shift_vector(self) -> Vector:
        """The shift necessary to place the BB at the radiation isocenter.
        The values are in the coordinates defined in the documentation.

        The shift is based on the paper by Low et al. See online documentation for more.
        """
        A = np.empty([2 * len(self.images), 3])
        epsilon = np.empty([2 * len(self.images), 1])
        for idx, img in enumerate(self.images):
            g = img.gantry_angle
            c = img.couch_angle_varian_scale
            A[2 * idx:2 * idx + 2, :] = np.array([[-cos(c), -sin(c), 0],
                                                  [-cos(g) * sin(c), cos(g) * cos(c), -sin(g)],
                                                  ])  # equation 6 (minus delta)
            epsilon[2 * idx:2 * idx + 2] = np.array([[img.cax2bb_vector.y], [img.cax2bb_vector.x]])  # equation 7

        B = linalg.pinv(A)
        delta = B.dot(epsilon)  # equation 9
        return Vector(x=delta[1][0], y=-delta[0][0], z=-delta[2][0])

    def bb_shift_instructions(self, couch_vrt: Optional[float] = None, couch_lng: Optional[float] = None,
                              couch_lat: Optional[float] = None) -> str:
        """Returns a string describing how to shift the BB to the radiation isocenter looking from the foot of the couch.
        Optionally, the current couch values can be passed in to get the new couch values. If passing the current
        couch position all values must be passed.

        Parameters
        ----------
        couch_vrt : float
            The current couch vertical position in cm.
        couch_lng : float
            The current couch longitudinal position in cm.
        couch_lat : float
            The current couch lateral position in cm.
        """
        sv = self.bb_shift_vector
        x_dir = 'LEFT' if sv.x < 0 else 'RIGHT'
        y_dir = 'IN' if sv.y > 0 else 'OUT'
        z_dir = 'UP' if sv.z > 0 else 'DOWN'
        move = f"{x_dir} {abs(sv.x):2.2f}mm; {y_dir} {abs(sv.y):2.2f}mm; {z_dir} {abs(sv.z):2.2f}mm"
        if all(val is not None for val in [couch_vrt, couch_lat, couch_lng]):
            new_lat = round(couch_lat + sv.x/10, 2)
            new_vrt = round(couch_vrt + sv.z/10, 2)
            new_lng = round(couch_lng + sv.y/10, 2)
            move += f"\nNew couch coordinates (mm): VRT: {new_vrt:3.2f}; LNG: {new_lng:3.2f}; LAT: {new_lat:3.2f}"
        return move

    @argue.options(value=('all', 'range'))
    def axis_rms_deviation(self, axis: Axis = Axis.GANTRY, value: str = 'all') -> Union[Iterable, float]:
        """The RMS deviations of a given axis/axes.

        Parameters
        ----------
        axis : ('Gantry', 'Collimator', 'Couch', 'Epid', 'GB Combo',  'GBP Combo'}
            The axis desired.
        value : {'all', 'range'}
            Whether to return all the RMS values from all images for that axis, or only return the maximum range of
            values, i.e. the 'sag'.
        """
        axis = convert_to_enum(axis, Axis)
        if axis != Axis.EPID:
            attr = 'cax2bb_vector'
        else:
            attr = 'cax2epid_vector'
            axis = (Axis.GANTRY, Axis.COLLIMATOR, Axis.REFERENCE)
        imgs = self._get_images(axis=axis)[1]
        if len(imgs) <= 1:
            return (0, )
        rms = [getattr(img, attr).as_scalar() for img in imgs]
        if value == 'range':
            rms = max(rms) - min(rms)
        return rms

    @argue.options(metric=('max', 'median'))
    def cax2bb_distance(self, metric: str = 'max') -> float:
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

    @argue.options(metric=('max', 'median'))
    def cax2epid_distance(self, metric: str = 'max') -> float:
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

    def _plot_deviation(self, axis: Axis, ax: Optional[plt.Axes] = None, show: bool=True):
        """Helper function: Plot the sag in Cartesian coordinates.

        Parameters
        ----------
        axis : {'gantry', 'epid', 'collimator', 'couch'}
            The axis to plot.
        ax : None, matplotlib.Axes
            The axis to plot to. If None, creates a new plot.
        show : bool
            Whether to show the image.
        """
        title = f'In-plane {axis.value} displacement'

        if axis == Axis.EPID:
            attr = 'cax2epid_vector'
            axis = Axis.GANTRY
        else:
            attr = 'cax2bb_vector'
        if not self._find_field:
            attr = 'epid2bb_vector'
        # get axis images, angles, and shifts
        imgs = [image for image in self.images if image.variable_axis in (axis, Axis.REFERENCE)]
        angles = [getattr(image, '{}_angle'.format(axis.value.lower())) for image in imgs]
        xz_sag = np.array([getattr(img, attr).x for img in imgs])
        y_sag = np.array([getattr(img, attr).y for img in imgs])
        rms = np.sqrt(xz_sag**2+y_sag**2)

        # plot the axis deviation
        if ax is None:
            ax = plt.subplot(111)
        ax.plot(angles, y_sag, 'bo', label='Y-axis', ls='-.')
        ax.plot(angles, xz_sag, 'm^', label='X/Z-axis', ls='-.')
        ax.plot(angles, rms, 'g+', label='RMS', ls='-')
        ax.set_title(title)
        ax.set_ylabel('mm')
        ax.set_xlabel(f"{axis.value} angle")
        ax.set_xticks(np.arange(0, 361, 45))
        ax.set_xlim(-15, 375)
        ax.grid(True)
        ax.legend(numpoints=1)
        if show:
            plt.show()

    def _get_images(self, axis: Union[Axis, Tuple[Axis, ...]] = (Axis.GANTRY,)) -> Tuple[float, list]:
        if isinstance(axis, Axis):
            axis = (axis,)
        images = [image for image in self.images if image.variable_axis in axis]
        return len(images), images

    def plot_axis_images(self, axis: Axis = Axis.GANTRY, show: bool = True, ax: Optional[plt.Axes] = None):
        """Plot all CAX/BB/EPID positions for the images of a given axis.

        For example, axis='Couch' plots a reference image, and all the BB points of the other
        images where the couch was moving.

        Parameters
        ----------
        axis : {'Gantry', 'Collimator', 'Couch', 'GB Combo',  'GBP Combo'}
            The images/markers from which accelerator axis to plot.
        show : bool
            Whether to actually show the images.
        ax : None, matplotlib.Axes
            The axis to plot to. If None, creates a new plot.
        """
        axis = convert_to_enum(axis, Axis)
        images = [image for image in self.images if image.variable_axis in (axis, Axis.REFERENCE)]
        ax = images[0].plot(show=False, ax=ax)  # plots the first marker; plot the rest of the markers below
        if axis != Axis.COUCH:
            # plot EPID
            epid_xs = [img.epid.x for img in images[1:]]
            epid_ys = [img.epid.y for img in images[1:]]
            ax.plot(epid_xs, epid_ys, 'b+', ms=8)
            # get CAX positions
            if self._find_field:
                xs = [img.field_cax.x for img in images[1:]]
                ys = [img.field_cax.y for img in images[1:]]
                marker = 'gs'
        else:
            # get BB positions
            xs = [img.bb.x for img in images[1:]]
            ys = [img.bb.y for img in images[1:]]
            marker = 'ro'
        if self._find_field:
            ax.plot(xs, ys, marker, ms=8)
        # set labels
        ax.set_title(axis.value + ' wobble')
        ax.set_xlabel(axis.value + ' positions superimposed')
        if self._find_field:
            ax.set_ylabel(axis.value + f" iso size: {getattr(self, axis.value.lower() + '_iso_size'):3.2f}mm")
        if show:
            plt.show()

    def plot_images(self, axis: Axis = Axis.GANTRY, show: bool = True):
        """Plot a grid of all the images acquired.

        Four columns are plotted with the titles showing which axis that column represents.

        Parameters
        ----------
        axis : {'Gantry', 'Collimator', 'Couch', 'GB Combo', 'GBP Combo', 'All'}
        show : bool
            Whether to show the image.
        """
        axis = convert_to_enum(axis, Axis)
        if not self._is_analyzed:
            raise ValueError("The set is not analyzed. Use .analyze() first.")

        def plot_image(image, axis):
            """Helper function to plot a WLImage to an axis."""
            if image is None:
                axis.set_frame_on(False)
                axis.axis('off')
            else:
                image.plot(ax=axis, show=False)

        # get axis images
        if axis == Axis.GANTRY:
            images = [image for image in self.images if image.variable_axis in (Axis.GANTRY, Axis.REFERENCE)]
        elif axis == Axis.COLLIMATOR:
            images = [image for image in self.images if image.variable_axis in (Axis.COLLIMATOR, Axis.REFERENCE)]
        elif axis == Axis.COUCH:
            images = [image for image in self.images if image.variable_axis in (Axis.COUCH, Axis.REFERENCE)]
        elif axis == Axis.GB_COMBO:
            images = [image for image in self.images if image.variable_axis in (Axis.GB_COMBO, Axis.GANTRY, Axis.COLLIMATOR, Axis.REFERENCE)]
        elif axis == Axis.GBP_COMBO:
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
        fig.suptitle(f"{axis.value} images", fontsize=14, y=1)
        plt.tight_layout()
        if show:
            plt.show()

    def save_images(self, filename: Union[str, BinaryIO], axis: Axis = Axis.GANTRY, **kwargs):
        """Save the figure of `plot_images()` to file. Keyword arguments are passed to `matplotlib.pyplot.savefig()`.

        Parameters
        ----------
        filename : str
            The name of the file to save to.
        axis
            The axis to save.
        """
        self.plot_images(axis=axis, show=False)
        plt.savefig(filename, **kwargs)

    def plot_summary(self, show: bool = True):
        """Plot a summary figure showing the gantry sag and wobble plots of the three axes."""
        if not self._is_analyzed:
            raise ValueError("The set is not analyzed. Use .analyze() first.")
        plt.figure(figsize=(11, 9))
        grid = (3, 6)
        if self._find_field:
            gantry_sag_ax = plt.subplot2grid(grid, (0, 0), colspan=3)
            self._plot_deviation(Axis.GANTRY, gantry_sag_ax, show=False)
        epid_sag_ax = plt.subplot2grid(grid, (0, 3), colspan=3)
        self._plot_deviation(Axis.EPID, epid_sag_ax, show=False)
        if self._get_images((Axis.COLLIMATOR, Axis.REFERENCE))[0] > 1:
            coll_sag_ax = plt.subplot2grid(grid, (1, 0), colspan=3)
            self._plot_deviation(Axis.COLLIMATOR, coll_sag_ax, show=False)
        if self._get_images((Axis.COUCH, Axis.REFERENCE))[0] > 1:
            couch_sag_ax = plt.subplot2grid(grid, (1, 3), colspan=3)
            self._plot_deviation(Axis.COUCH, couch_sag_ax, show=False)

        for axis, axnum in zip((Axis.GANTRY, Axis.COLLIMATOR, Axis.COUCH), (0, 2, 4)):
            if self._get_images((axis, Axis.REFERENCE))[0] > 1:
                ax = plt.subplot2grid(grid, (2, axnum), colspan=2)
                self.plot_axis_images(axis=axis, ax=ax, show=False)
        if show:
            plt.tight_layout()
            plt.show()

    def save_summary(self, filename: Union[str, BinaryIO], **kwargs):
        """Save the summary image."""
        self.plot_summary(show=False)
        plt.tight_layout()
        plt.savefig(filename, **kwargs)

    def results(self, as_list: bool = False) -> str:
        """Return the analysis results summary.

        Parameters
        ----------
        as_list : bool
            Whether to return as a list of strings vs single string. Pretty much for internal usage.
        """
        if not self._is_analyzed:
            raise ValueError("The set is not analyzed. Use .analyze() first.")
        num_gantry_imgs = self._get_images(axis=(Axis.GANTRY, Axis.REFERENCE))[0]
        num_gantry_coll_imgs = self._get_images(axis=(Axis.GANTRY, Axis.COLLIMATOR, Axis.GB_COMBO, Axis.REFERENCE))[0]
        num_coll_imgs = self._get_images(axis=(Axis.COLLIMATOR, Axis.REFERENCE))[0]
        num_couch_imgs = self._get_images(axis=(Axis.COUCH, Axis.REFERENCE))[0]
        num_imgs = len(self.images)
        result = ["Winston-Lutz Analysis",
                  "=================================",
                  f"Number of images: {num_imgs}",
                  f"Maximum 2D CAX->BB distance: {self.cax2bb_distance('max'):.2f}mm",
                  f"Median 2D CAX->BB distance: {self.cax2bb_distance('median'):.2f}mm",
                  f"Shift to iso: facing gantry, move BB: {self.bb_shift_instructions()}",
                  f"Gantry 3D isocenter diameter: {self.gantry_iso_size:.2f}mm ({num_gantry_imgs}/{num_imgs} images considered)",
                  f"Maximum Gantry RMS deviation (mm): {max(self.axis_rms_deviation(Axis.GANTRY)):.2f}mm",
                  f"Maximum EPID RMS deviation (mm): {max(self.axis_rms_deviation(Axis.EPID)):.2f}mm",
                  f"Gantry+Collimator 3D isocenter diameter: {self.gantry_coll_iso_size:.2f}mm ({num_gantry_coll_imgs}/{num_imgs} images considered)",
                  f"Collimator 2D isocenter diameter: {self.collimator_iso_size:.2f}mm ({num_coll_imgs}/{num_imgs} images considered)",
                  f"Maximum Collimator RMS deviation (mm): {max(self.axis_rms_deviation(Axis.COLLIMATOR)):.2f}",
                  f"Couch 2D isocenter diameter: {self.couch_iso_size:.2f}mm ({num_couch_imgs}/{num_imgs} images considered)",
                  f"Maximum Couch RMS deviation (mm): {max(self.axis_rms_deviation(Axis.COUCH)):.2f}"
        ]
        if not as_list:
            result = '\n'.join(result)
        return result

    def results_data(self, as_dict=False) -> Union[WinstonLutzResult, dict]:
        """Present the results data and metadata as a dataclass or dict.
        The default return type is a dataclass."""
        if not self._is_analyzed:
            raise ValueError("The set is not analyzed. Use .analyze() first.")
        num_gantry_imgs = self._get_images(axis=(Axis.GANTRY, Axis.REFERENCE))[0]
        num_gantry_coll_imgs = self._get_images(axis=(Axis.GANTRY, Axis.COLLIMATOR, Axis.GB_COMBO, Axis.REFERENCE))[0]
        num_coll_imgs = self._get_images(axis=(Axis.COLLIMATOR, Axis.REFERENCE))[0]
        num_couch_imgs = self._get_images(axis=(Axis.COUCH, Axis.REFERENCE))[0]

        data = WinstonLutzResult(
                num_total_images=len(self.images),
                num_gantry_images=num_gantry_imgs,
                num_coll_images=num_coll_imgs,
                num_gantry_coll_images=num_gantry_coll_imgs,
                num_couch_images=num_couch_imgs,
                max_2d_cax_to_bb_mm=self.cax2bb_distance('max'),
                median_2d_cax_to_bb_mm=self.cax2bb_distance('median'),
                max_2d_cax_to_epid_mm=self.cax2epid_distance('max'),
                median_2d_cax_to_epid_mm=self.cax2epid_distance('median'),
                coll_2d_iso_diameter_mm=self.collimator_iso_size,
                couch_2d_iso_diameter_mm=self.couch_iso_size,
                gantry_3d_iso_diameter_mm=self.gantry_iso_size,
                gantry_coll_3d_iso_diameter_mm=self.gantry_coll_iso_size,
                max_gantry_rms_deviation_mm=max(self.axis_rms_deviation(axis=Axis.GANTRY)),
                max_coll_rms_deviation_mm=max(self.axis_rms_deviation(axis=Axis.COLLIMATOR)),
                max_couch_rms_deviation_mm=max(self.axis_rms_deviation(axis=Axis.COUCH)),
                max_epid_rms_deviation_mm=max(self.axis_rms_deviation(axis=Axis.EPID)),
        )
        if as_dict:
            return dataclasses.asdict(data)
        return data

    def publish_pdf(self, filename: str, notes: Optional[Union[str, List[str]]] = None, open_file: bool = False,
                    metadata: Optional[dict] = None):
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
            raise ValueError("The set is not analyzed. Use .analyze() first.")
        plt.ioff()
        title = "Winston-Lutz Analysis"
        canvas = pdf.PylinacCanvas(filename, page_title=title, metadata=metadata)
        text = self.results(as_list=True)
        canvas.add_text(text=text, location=(7, 25.5))
        # draw summary image on 1st page
        data = io.BytesIO()
        self.save_summary(data, figsize=(8, 8))
        canvas.add_image(image_data=data, location=(2, 3), dimensions=(16, 16))
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 4.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 4))
        # add more pages showing individual axis images
        for ax in (Axis.GANTRY, Axis.COLLIMATOR, Axis.COUCH, Axis.GB_COMBO, Axis.GBP_COMBO):
            if self._contains_axis_images(ax):
                canvas.add_new_page()
                data = io.BytesIO()
                self.save_images(data, axis=ax, figsize=(10, 10))
                canvas.add_image(data, location=(2, 7), dimensions=(18, 18))

        canvas.finish()

        if open_file:
            open_path(filename)

    def _contains_axis_images(self, axis: Axis = Axis.GANTRY) -> bool:
        """Return whether or not the set of WL images contains images pertaining to a given axis"""
        return any(True for image in self.images if image.variable_axis in (axis,))


class WinstonLutz2D(image.LinacDicomImage):
    """Holds individual Winston-Lutz EPID images, image properties, and automatically finds the field CAX and BB."""
    bb: Point
    field_cax: Point
    _find_field: bool
    _rad_field_bounding_box: list

    def __init__(self, file: str, use_filenames: bool = False):
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
        self._is_analyzed = False

    def analyze(self, bb_size_mm: float = 5, find_field=True) -> None:
        """Analyze the image."""
        self.flipud()
        self._clean_edges()
        self.ground()
        self.normalize()
        self._find_field = find_field
        if find_field:
            self.field_cax, self._rad_field_bounding_box = self._find_field_centroid()
        self.bb = self._find_bb(bb_size_mm)
        self._is_analyzed = True

    def __repr__(self):
        return f"WLImage(G={self.gantry_angle:.1f}, B={self.collimator_angle:.1f}, P={self.couch_angle:.1f})"

    def _clean_edges(self, window_size: int = 2) -> None:
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
        filled_img = ndimage.binary_fill_holes(threshold_img)
        # clean single-pixel noise from outside field
        cleaned_img = ndimage.binary_erosion(threshold_img)
        [*edges] = bounding_box(cleaned_img)
        edges[0] -= 10
        edges[1] += 10
        edges[2] -= 10
        edges[3] += 10
        coords = ndimage.measurements.center_of_mass(filled_img)
        p = Point(x=coords[-1], y=coords[0])
        return p, edges

    def _find_bb(self, bb_size: float) -> Point:
        """Find the BB within the radiation field. Crops down the image to the central 30x30mm square.
        Iteratively searches for a circle-like object
        by lowering a low-pass threshold value until found.

        Returns
        -------
        Point
            The weighted-pixel value location of the BB.
        """
        # crop to the central region to avoid edge values influencing histogram values
        centered_box_side_pixels = 15 * self.dpmm
        left, right = int(self.shape[1]/2 - centered_box_side_pixels), int(self.shape[1]/2 + centered_box_side_pixels)
        bottom, top = int(self.shape[0] / 2 - centered_box_side_pixels), int(self.shape[0] / 2 + centered_box_side_pixels)
        central_array = self.array[bottom:top, left:right]

        # create an intensity image; this is used later to weight the BB center based on pixel value
        intensity_image = ArrayImage(central_array)
        intensity_image.invert()

        # search for the BB by iteratively lowering the low-pass threshold value until the BB is found.
        spread = central_array.max() - central_array.min()
        threshold = central_array.max()
        found = False
        while not found:
            try:
                # this is an inversion threshold; values less than the threshold are set to positive values; this makes the BB show up as positive
                binary_arr = (central_array < threshold)
                min_size_1mm_circle = np.pi * self.dpmm
                array = remove_small_objects(binary_arr, min_size=min_size_1mm_circle)
                max_size_bb = np.pi * bb_size * self.dpmm
                array = remove_small_holes(array, area_threshold=max_size_bb)
                # use below for debugging
                # plt.imshow(array)
                # plt.show()
                labeled_arr = label(array)
                bb_regionprops = measure.regionprops(labeled_arr, intensity_image=intensity_image)
                # look for the BB in the regionprops
                for region in bb_regionprops:
                    if is_round(region) and is_modest_size(region, self.dpmm, bb_size) and is_symmetric(region):
                        found = True
                        break
                if not found:
                    raise ValueError
            except (IndexError, ValueError):
                threshold -= 0.03 * spread
                if threshold < central_array.min():
                    raise ValueError("Unable to locate the BB. Make sure the field edges do not obscure the BB and that there is no artifacts in the images.")
            else:
                found = True

        # determine the center of mass of the BB
        # restore the centroid position to the full image coordinates
        x, y = region.weighted_centroid[1] + left, region.weighted_centroid[0] + bottom
        return Point(x, y)

    @property
    def epid(self) -> Point:
        """Center of the EPID panel"""
        return self.center

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
        p1.x = self.cax2bb_vector.x*cos(self.gantry_angle) + 20 * sin(self.gantry_angle)
        p1.y = self.cax2bb_vector.x*-sin(self.gantry_angle) + 20 * cos(self.gantry_angle)
        p1.z = self.cax2bb_vector.y
        # point 2 - ray destination
        p2.x = self.cax2bb_vector.x*cos(self.gantry_angle) - 20 * sin(self.gantry_angle)
        p2.y = self.cax2bb_vector.x*-sin(self.gantry_angle) - 20 * cos(self.gantry_angle)
        p2.z = self.cax2bb_vector.y
        l = Line(p1, p2)
        return l

    @property
    def couch_angle_varian_scale(self) -> float:
        """The couch angle converted from IEC 61217 scale to "Varian" scale. Note that any new Varian machine uses 61217."""
        #  convert to Varian scale per Low paper scale
        if super().couch_angle > 250:
            return 2*270-super().couch_angle
        else:
            return 180 - super().couch_angle

    @property
    def epid2bb_distance(self) -> Vector:
        """The distance from the BB to the EPID center in mm."""
        return self.epid.distance_to(self.bb)

    @property
    def epid2bb_vector(self) -> Vector:
        """The distance from the BB to the EPID center in mm."""
        dist = (self.bb - self.epid) / self.dpmm
        return Vector(dist.x, dist.y, dist.z)

    @property
    def cax2bb_vector(self) -> Vector:
        """The vector in mm from the CAX to the BB."""
        if not self._find_field:
            raise AnalysisError("Field-only was set to True; any property or method involving the field CAX is forbidden.")
        dist = (self.bb - self.field_cax) / self.dpmm
        return Vector(dist.x, dist.y, dist.z)

    @property
    def cax2bb_distance(self) -> float:
        """The scalar distance in mm from the CAX to the BB."""
        if not self._find_field:
            raise AnalysisError("Field-only was set to True; any property or method involving the field CAX is forbidden.")
        dist = self.field_cax.distance_to(self.bb)
        return dist / self.dpmm

    @property
    def cax2epid_vector(self) -> Vector:
        """The vector in mm from the CAX to the EPID center pixel"""
        if not self._find_field:
            raise AnalysisError("Field-only was set to True; any property or method involving the field CAX is forbidden.")
        dist = (self.epid - self.field_cax) / self.dpmm
        return Vector(dist.x, dist.y, dist.z)

    @property
    def cax2epid_distance(self) -> float:
        """The scalar distance in mm from the CAX to the EPID center pixel"""
        if not self._find_field:
            raise AnalysisError("Field-only was set to True; any property or method involving the field CAX is forbidden.")
        return self.field_cax.distance_to(self.epid) / self.dpmm

    def plot(self, ax: Optional[plt.Axes] = None, show: bool = True, clear_fig: bool = False):
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
        if self._find_field:
            ax.plot(self.field_cax.x, self.field_cax.y, 'gs', ms=8)
        ax.plot(self.bb.x, self.bb.y, 'ro', ms=8)
        ax.axvline(self.epid.x, color='b', linewidth=0.5)
        ax.axhline(self.epid.y, color='b', linewidth=0.5)
        if self._find_field:
            ax.set_ylim([self._rad_field_bounding_box[0], self._rad_field_bounding_box[1]])
            ax.set_xlim([self._rad_field_bounding_box[2], self._rad_field_bounding_box[3]])
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_title('\n'.join(wrap(self.file, 30)), fontsize=10)
        ax.set_xlabel(f"G={self.gantry_angle:.0f}, B={self.collimator_angle:.0f}, P={self.couch_angle:.0f}")
        if self._find_field:
            ax.set_ylabel(f"CAX to BB: {self.cax2bb_distance:3.2f}mm")
        else:
            ax.set_ylabel(f"EPID to BB: {self.epid2bb_distance:3.2f}mm")
        if show:
            plt.show()
        return ax

    def save_plot(self, filename: str, **kwargs):
        """Save the image plot to file."""
        self.plot(show=False)
        plt.savefig(filename, **kwargs)

    @property
    def variable_axis(self) -> Axis:
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
            return Axis.COUCH
        elif G0 and P0 and not B0:
            return Axis.COLLIMATOR
        elif P0 and B0 and not G0:
            return Axis.GANTRY
        elif P0 and B0 and G0:
            return Axis.REFERENCE
        elif P0:
            return Axis.GB_COMBO
        else:
            return Axis.GBP_COMBO

    def results_data(self, as_dict: bool = False) -> Union[WinstonLutz2DResult, dict]:
        """Present the results data and metadata as a dataclass or dict.
        The default return type is a dataclass."""
        if not self._is_analyzed:
            raise ValueError("The image is not analyzed. Use .analyze() first.")

        data = WinstonLutz2DResult(
                variable_axis=self.variable_axis.value,
                cax2bb_vector=self.cax2bb_vector,
                cax2epid_vector=self.cax2epid_vector,
                cax2bb_distance=self.cax2bb_distance,
                cax2epid_distance=self.cax2epid_distance,
                bb_location=self.bb,
                field_cax=self.field_cax,
        )
        if as_dict:
            return dataclasses.asdict(data)
        return data


def is_symmetric(region: RegionProperties) -> bool:
    """Whether the binary object's dimensions are symmetric, i.e. a perfect circle. Used to find the BB."""
    ymin, xmin, ymax, xmax = region.bbox
    y = abs(ymax - ymin)
    x = abs(xmax - xmin)
    if x > max(y * 1.05, y + 3) or x < min(y * 0.95, y - 3):
        return False
    return True


def is_near_center(region: RegionProperties, dpmm: float, shape: Tuple[float, float]) -> bool:
    """Whether the bb is <2cm from the center of the field"""
    extent_limit_mm = 20
    bottom, left, top, right = region.bbox
    bb_center_x = left + (right - left) / 2
    bb_center_y = bottom + (top - bottom) / 2
    x_lo_limit = shape[1]/2 - dpmm * extent_limit_mm
    x_hi_limit = shape[1]/2 + dpmm * extent_limit_mm
    is_bb_x_centered = x_lo_limit < bb_center_x < x_hi_limit
    y_lo_limit = shape[0]/2 - dpmm * extent_limit_mm
    y_hi_limit = shape[0]/2 + dpmm * extent_limit_mm
    is_bb_y_centered = y_lo_limit < bb_center_y < y_hi_limit
    return is_bb_x_centered and is_bb_y_centered


def is_modest_size(region: RegionProperties, dpmm: float, bb_size: float) -> bool:
    """Decide whether the ROI is roughly the size of a BB; not noise and not an artifact. Used to find the BB."""
    bb_area = region.filled_area / (dpmm ** 2)
    bb_size = max((bb_size, 2.1))
    larger_bb_area = np.pi * ((bb_size+2)/2)**2
    smaller_bb_area = np.pi * ((bb_size-2)/2)**2
    return smaller_bb_area < bb_area < larger_bb_area


def is_round(region: RegionProperties) -> bool:
    """Decide if the ROI is circular in nature by testing the filled area vs bounding box. Used to find the BB."""
    expected_fill_ratio = np.pi / 4  # area of a circle inside a square
    actual_fill_ratio = region.filled_area / region.bbox_area
    return expected_fill_ratio * 1.2 > actual_fill_ratio > expected_fill_ratio * 0.8


class AnalysisError(Exception):
    pass
