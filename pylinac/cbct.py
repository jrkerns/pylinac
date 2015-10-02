
"""The CBCT module automatically analyzes DICOM images of a CatPhan 504 acquired when doing CBCT or regular CT quality assurance.
It can load a folder or zip file that the images are in and automatically correct for phantom setup in 6 degrees.
It can analyze the HU regions and image scaling (CTP404), the high-contrast line pairs (CTP528) to calculate the modulation transfer function (MTF),
the HU uniformity (CTP486), and Low Contrast (CTP515) on the corresponding slices.

Currently only the CatPhan 504 (Varian's default phantom) is supported, but CatPhan 503 (Elekta) support is being worked on.

Features:

* **Automatic phantom registration** - Your phantom can be tilted, rotated, or translated--pylinac will register the phantom.
* **Automatic testing of 4 major modules** - Major modules are automatically registered and analyzed.
* **Any scan protocol** - Scan your CatPhan 504 with any Varian protocol; or even scan it in a regular CT scanner.
  Any field size or field extent is allowed.
"""
import os.path as osp
import zipfile
from io import BytesIO
from functools import lru_cache
from collections import OrderedDict
from abc import abstractmethod

import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt

from pylinac.core.decorators import value_accept
from pylinac.core.image import Image, DICOMStack
from pylinac.core.geometry import Point, Circle, sector_mask, Line, Rectangle
from pylinac.core.profile import MultiProfile, CollapsedCircleProfile, SingleProfile
from pylinac.core.io import get_folder_UI, get_filepath_UI
from pylinac.core.utilities import simple_round, import_mpld3, get_url

np.seterr(invalid='ignore')  # ignore warnings for invalid numpy operations. Used for np.where() operations on partially-NaN arrays.


class CBCT:
    """A class for loading and analyzing CT DICOM files of a CatPhan 504. Can be from a CBCT or CT scanner
    Analyzes: Uniformity (CTP486), High-Contrast Spatial Resolution (CTP528), Image Scaling & HU Linearity (CTP404).

    Attributes
    ----------
    dicom_stack : :class:`~pylinac.core.image.DICOMStack`
    settings : :class:`~pylinac.cbct.Settings`
    hu : :class:`~pylinac.cbct.HUSlice`
    uniformity : :class:`~pylinac.cbct.UniformitySlice`
    geometry : :class:`~pylinac.cbct.GeometrySlice`
    lowcontrast: :class:`~pylinac.cbct.LowContrastSlice`
    spatialres : :class:`~pylinac.cbct.SpatialResolutionSlice`
    thickness : :class:`~pylinac.cbct.ThicknessSlice`

    Examples
    --------
    Run the demo:
        >>> CBCT().run_demo()

    Typical session:
        >>> cbct_folder = r"C:/QA/CBCT/June"
        >>> mycbct = CBCT(cbct_folder)
        >>> mycbct.analyze()
        >>> print(mycbct.return_results())
        >>> mycbct.plot_analyzed_image()
    """
    def __init__(self, folderpath=None):
        """
        Parameters
        ----------
        folderpath : str, None
            If None, the images must be loaded later.
            If a string, must point to the CBCT image folder location.
        """
        self.settings = None
        self.hu = None
        self.uniformity = None
        self.geometry = None
        self.lowcontrast = None
        self.spatialres = None
        if folderpath is not None:
            self.load_folder(folderpath)

    @classmethod
    def from_demo_images(cls):
        """Construct a CBCT object and load the demo images.

        .. versionadded:: 0.6
        """
        obj = cls()
        obj.load_demo_images()
        return obj

    def load_demo_images(self):
        """Load the CBCT demo images."""
        cbct_demo_zip = osp.join(osp.dirname(osp.abspath(__file__)), 'demo_files', 'cbct', 'High quality head.zip')
        self.load_zip_file(cbct_demo_zip)

    @classmethod
    def from_url(cls, url):
        """Instantiate from a URL.

        .. versionadded:: 0.7.1
        .. note:: ``requests`` must be installed to use this feature.
        """
        obj = cls()
        obj.load_url(url)
        return obj

    def load_url(self, url):
        """Load from a URL.

        .. versionadded:: 0.7.1
        .. note:: ``requests`` must be installed to use this feature.
        """
        response = get_url(url)
        stream = zipfile.ZipFile(BytesIO(response.content))
        self.load_zip_file(stream)

    @classmethod
    def from_folder_UI(cls):
        """Construct a CBCT object an get the files using a UI dialog box.

        .. versionadded:: 0.6
        """
        obj = cls()
        obj.load_folder_UI()
        return obj

    def load_folder_UI(self):
        """Load the CT DICOM files from a folder using a UI dialog box."""
        folder = get_folder_UI()
        if folder:
            self.load_folder(folder)

    def load_folder(self, folder):
        """Load the CT DICOM files string input.

        Parameters
        ----------
        folder : str
            Path to the folder.

        Raises
        ------
        NotADirectoryError : If folder str passed is not a valid directory.
        FileNotFoundError : If no CT images are found in the folder
        """
        if not osp.isdir(folder):
            raise NotADirectoryError("Path given was not a Directory/Folder")
        self.dicom_stack = DICOMStack(folder)
        self.settings = Settings(self.dicom_stack)

    @classmethod
    def from_zip_file_UI(cls):
        """Construct a CBCT object and pass the zip file.

        .. versionadded:: 0.6
        """
        obj = cls()
        obj.load_zip_file_UI()
        return obj

    def load_zip_file_UI(self):
        """Load a zip file using a UI dialog box.

        .. versionadded:: 0.6
        """
        zfile = get_filepath_UI()
        self.load_zip_file(zfile)

    @classmethod
    def from_zip_file(cls, zip_file):
        """Construct a CBCT object and pass the zip file.

        .. versionadded:: 0.6
        """
        obj = cls()
        obj.load_zip_file(zip_file)
        return obj

    def load_zip_file(self, zip_file):
        """Load a CBCT dataset from a zip file.

        Parameters
        ----------
        zip_file : str, ZipFile
            Path to the zip file or a ZipFile object.

        Raises
        ------
        FileExistsError : If zip_file passed was not a legitimate zip file.
        FileNotFoundError : If no CT images are found in the folder
        """
        self.dicom_stack = DICOMStack.from_zip(zip_file)
        self.settings = Settings(self.dicom_stack)

    def plot_analyzed_image(self, hu=True, geometry=True, uniformity=True, spatial_res=True, thickness=True, low_contrast=True, show=True):
        """Plot the images used in the calculate and summary data.

        Parameters
        ----------
        hu : bool
            Whether to show the HU linearity circles.
        geometry : bool
            Whether to show the geometric lines found.
        uniformity : bool
            Whether to show the uniformity ROIs.
        spatial_res : bool
            Whether to show the spatial resolution circle outline.
        thickness : bool
            Whether to show the wire ramp boxes.
        low_contrast : bool
            Whether to show the low contrast bubble circles.
        show : bool
            Whether to plot the image or not.
        """
        # set up grid and axes
        grid_size = (2,4)
        unif_ax = plt.subplot2grid(grid_size, (0, 0))
        hu_ax = plt.subplot2grid(grid_size, (0, 1))
        sr_ax = plt.subplot2grid(grid_size, (1, 0))
        lowcon_ax = plt.subplot2grid(grid_size, (1, 1))
        hu_lin_ax = plt.subplot2grid(grid_size, (0, 2))
        mtf_ax = plt.subplot2grid(grid_size, (0, 3))
        unif_prof_ax = plt.subplot2grid(grid_size, (1, 2), colspan=2)

        # plot the images
        show_section = [uniformity, hu, spatial_res, low_contrast]
        axes = (unif_ax, hu_ax, sr_ax, lowcon_ax)
        items = ('uniformity', 'hu', 'spatialres', 'lowcontrast')
        titles = ('Uniformity', 'HU & Geometry', 'Spatial Resolution', 'Low Contrast')
        for show_unit, axis, title, item in zip(show_section, axes, titles, items):
            if show_unit:
                klass = getattr(self, item)
                axis.imshow(klass.image.array, cmap=plt.cm.Greys)
                klass.plot_rois(axis)
                axis.autoscale(tight=True)
                axis.set_title(title)
                axis.axis('off')
        # geometry & slice thickness don't have their own axes
        if geometry:
            self.geometry.plot_lines(hu_ax)
        if thickness:
            self.thickness.plot_rois(hu_ax)

        # plot the other sections
        self.hu.plot_linearity(hu_lin_ax)
        hu_lin_ax.set_title("HU linearity")
        self.uniformity.plot_profiles(unif_prof_ax)
        unif_prof_ax.set_title("Uniformity Profiles")
        self.spatialres.plot_mtf(mtf_ax)
        mtf_ax.set_title('RMTF')

        # finish up
        plt.tight_layout()
        if show:
            plt.show()

    def save_analyzed_image(self, filename, **kwargs):
        """Save the analyzed summary plot.

        Parameters
        ----------
        filename : str, file object
            The name of the file to save the image to.
        kwargs :
            Any valid matplotlib kwargs.
        """
        self.plot_analyzed_image(show=False)
        plt.savefig(filename, **kwargs)

    def plot_analyzed_subimage(self, subimage='hu', delta=True, interactive=False, show=True):
        """Plot a specific component of the CBCT analysis.

        Parameters
        ----------
        subimage : {'hu', 'un', 'sp', 'lc', 'mtf', 'lin', 'prof'}
            The subcomponent to plot. Values must contain one of the following letter combinations.
            E.g. ``linearity``, ``linear``, and ``lin`` will all draw the HU linearity values.

            * ``hu`` draws the HU linearity image.
            * ``un`` draws the HU uniformity image.
            * ``sp`` draws the Spatial Resolution image.
            * ``lc`` draws the Low Contrast image.
            * ``mtf`` draws the RMTF plot.
            * ``lin`` draws the HU linearity values. Used with ``delta``.
            * ``prof`` draws the HU uniformity profiles.
        delta : bool
            Only for use with ``lin``. Whether to plot the HU delta or actual values.
        interactive : bool
            If False (default), the figure is plotted statically as a matplotlib figure.
            If True, the figure is plotted as an HTML page in the default browser that can have tooltips.
            This setting is only applicable for ``mtf``, ``lin``, and ``prof``, i.e. the "regular" plots. Image
            plotting does not work with mpld3.

            .. note:: mpld3 must be installed to use this feature.

        show : bool
            Whether to actually show the plot.
        """
        subimage = subimage.lower()
        plt.clf()
        plt.axis('off')
        if interactive:
            mpld3 = import_mpld3()

        if 'hu' in subimage:  # HU, GEO & thickness objects
            plt.imshow(self.hu.image.array, cmap=plt.cm.Greys)
            self.hu.plot_rois(plt.gca())
            self.geometry.plot_lines(plt.gca())
            self.thickness.plot_rois(plt.gca())
            plt.autoscale(tight=True)
        elif 'un' in subimage:  # uniformity
            plt.imshow(self.uniformity.image.array, cmap=plt.cm.Greys)
            self.uniformity.plot_rois(plt.gca())
            plt.autoscale(tight=True)
        elif 'sp' in subimage:  # SR objects
            plt.imshow(self.spatialres.image.array, cmap=plt.cm.Greys)
            self.spatialres.plot_rois(plt.gca())
            plt.autoscale(tight=True)
        elif 'lc' in subimage:  # low contrast objects
            plt.imshow(self.lowcontrast.image.array, cmap=plt.cm.Greys)
            self.lowcontrast.plot_rois(plt.gca())
            plt.autoscale(tight=True)
        elif 'mtf' in subimage:
            subimage = 'mtf'
            plt.axis('on')
            points = self.spatialres.plot_mtf(plt.gca())
            if interactive:
                labels = ['MTF: {:3.3f}'.format(i) for i in self.spatialres.line_pair_mtfs]
                tooltip = mpld3.plugins.PointLabelTooltip(points[0], labels, location='top right')
                mpld3.plugins.connect(plt.gcf(), tooltip)
        elif 'lin' in subimage:
            subimage = 'lin'
            plt.axis('on')
            points = self.hu.plot_linearity(plt.gca(), delta)
            if interactive:
                delta_values = [roi.value_diff for roi in self.hu.rois.values()]
                actual_values = [roi.pixel_value for roi in self.hu.rois.values()]
                names = [roi for roi in self.hu.rois.keys()]
                labels = ['{} -- Actual: {:3.1f}; Difference: {:3.1f}'.format(name, actual, delta) for name, actual, delta in zip(names, actual_values, delta_values)]
                tooltip = mpld3.plugins.PointLabelTooltip(points[0], labels, location='top right')
                mpld3.plugins.connect(plt.gcf(), tooltip)
        elif 'prof' in subimage:
            subimage = 'prof'
            plt.axis('on')
            self.uniformity.plot_profiles(plt.gca())
        else:
            raise ValueError("Subimage parameter {} not understood".format(subimage))

        if show:
            if interactive and (subimage in ('mtf', 'lin', 'prof')):
                mpld3.show()
            else:
                plt.show()

    def save_analyzed_subimage(self, filename, subimage='hu', interactive=False, **kwargs):
        """Save a component image to file.

        Parameters
        ----------
        filename : str, file object
            The file to write the image to.
        subimage : str
            See :meth:`~pylinac.cbct.CBCT.plot_analyzed_subimage` for parameter info.
        interactive : bool
            If False (default), saves a matplotlib figure as a .png image.
            If True, saves an html file, which can be opened in a browser, etc.

            .. note:: mpld3 must be installed to use this feature.
        """
        self.plot_analyzed_subimage(subimage, interactive=interactive, show=False)
        if interactive and (subimage in ('mtf', 'lin', 'prof')):
            mpld3 = import_mpld3()
            mpld3.save_html(plt.gcf(), filename)
        else:
            plt.savefig(filename, **kwargs)
        print("CBCT subimage figure saved to {}".format(osp.abspath(filename)))

    def return_results(self):
        """Return the results of the analysis as a string. Use with print()."""
        #TODO: make prettier
        string = ('\n - CBCT QA Test - \n'
                  'HU Regions: {}\n'
                  'HU Passed?: {}\n'
                  'Uniformity: {}\n'
                  'Uniformity Passed?: {}\n'
                  'MTF 80% (lp/mm): {}\n'
                  'Geometric Line Average (mm): {}\n'
                  'Geometry Passed?: {}\n'
                  'Low Contrast ROIs visible: {}\n'
                  'Low Contrast Passed? {}\n'
                  'Slice Thickness (mm): {}\n'
                  'Slice Thickeness Passed? {}\n').format(self.hu.get_ROI_vals(), self.hu.overall_passed,
                                                   self.uniformity.get_ROI_vals(), self.uniformity.overall_passed, self.spatialres.mtf(80),
                                                   self.geometry.avg_line_length, self.geometry.overall_passed, self.lowcontrast.rois_visible,
                                                   self.lowcontrast.overall_passed, self.thickness.avg_slice_thickness,
                                                          self.thickness.passed)
        return string

    def _return_results(self):
        """Helper function to spit out values that will be tested."""
        print(self.return_results())
        print("Phantom roll: {}".format(self.settings.phantom_roll))
        mtfs = {}
        for mtf in (60, 70, 80, 90, 95):
            mtfval = self.spatialres.mtf(mtf)
            mtfs[mtf] = mtfval
        print(mtfs)
        for slice in ('hu_slice_num', 'un_slice_num', 'sr_slice_num', 'lc_slice_num'):
            slicenum = getattr(self.settings, slice)
            print(slice, slicenum)

    def analyze(self, hu_tolerance=40, scaling_tolerance=1, thickness_tolerance=0.2, low_contrast_tolerance=1, contrast_threshold=10):
        """Single-method full analysis of CBCT DICOM files.

        Parameters
        ----------
        hu_tolerance : int
            The HU tolerance value for both HU uniformity and linearity.
        scaling_tolerance : float, int
            The scaling tolerance in mm of the geometric nodes on the HU linearity slice (CTP404 module).
        thickness_tolerance : float, int
            The tolerance of the thickness calculation in mm, based on the wire ramps in the CTP404 module.

            .. warning:: Thickness accuracy degrades with image noise; i.e. low mAs images are less accurate.
        low_contrast_tolerance : int
            The number of low-contrast bubbles needed to be "seen" to pass.
        contrast_threshold : float, int
            The threshold for "detecting" low-contrast image. See RTD for calculation info.
        """
        if not self.images_loaded:
            raise AttributeError("Images not yet loaded")

        # set various setting values
        self.settings.hu_tolerance = hu_tolerance
        self.settings.scaling_tolerance = scaling_tolerance
        self.settings.thickness_tolerance = thickness_tolerance
        self.settings.low_contrast_tolerance = low_contrast_tolerance
        self.settings.contrast_threshold = contrast_threshold

        # Analysis
        self.hu = HUSlice(self.dicom_stack, self.settings)
        self.uniformity = UniformitySlice(self.dicom_stack, self.settings)
        self.spatialres = SpatialResolutionSlice(self.dicom_stack, self.settings)
        self.lowcontrast = LowContrastSlice(self.dicom_stack, self.settings)
        self.thickness = ThicknessSlice(self.dicom_stack, self.settings)
        self.geometry = GeometrySlice(self.dicom_stack, self.settings)

    def run_demo(self, show=True):
        """Run the CBCT demo using high-quality head protocol images."""
        self.load_demo_images()
        self.analyze()
        print(self.return_results())
        self.plot_analyzed_image(show)

    @property
    def images_loaded(self):
        """Boolean property specifying if the images have been loaded."""
        return False if self.settings is None else True


class Settings:
    """Data structure for retaining certain settings and information regarding the CBCT algorithm and image data.
    This class is initialized during the CBCT image loading.

    Attributes
    ----------
    threshold : int
        The threshold for converting the image to binary (for things like phantom position locating). Default is -800.
    air_bubble_radius_mm : int, float
        The size of the "Air" HU ROIs in mm; for finding the phantom roll.
    """
    threshold = -800
    hu_tolerance = 40
    scaling_tolerance = 1
    thickness_tolerance = 0.2
    low_contrast_tolerance = 1
    contrast_threshold = 10
    air_bubble_radius_mm = 6

    def __init__(self, dicom_stack):
        self.dicom_stack = dicom_stack

    @property
    def air_bubble_size(self):
        """The size of the "Air" HU ROIs in pixels; for finding the phantom roll."""
        return (np.pi*self.air_bubble_radius_mm**2)/self.mm_per_pixel**2

    @property
    def mm_per_pixel(self):
        """The millimeters per pixel of the DICOM images."""
        return self.dicom_stack.metadata.PixelSpacing[0]

    @property
    @lru_cache()
    def hu_slice_num(self):
        """Using a brute force search of the images, find the median HU linearity slice.

        This method walks through all the images and takes a collapsed circle profile where the HU
        linearity ROIs are. If the profile contains both low (<800) and high (>800) HU values and most values are the same
        (i.e. it's not an artifact), then
        it can be assumed it is an HU linearity slice. The median of all applicable slices is the
        center of the HU slice.

        Returns
        -------
        int
            The middle slice of the HU linearity module.
        """
        hu_slices = []
        for image_number in range(self.num_images):
            slice = Slice(self.dicom_stack, self, image_number, combine=False)
            try:
                center = slice.phan_center
            except ValueError:  # a slice without the phantom in view
                pass
            else:
                circle_prof = CollapsedCircleProfile(center, radius=59/self.mm_per_pixel, image_array=slice.image, width_ratio=0.05, num_profiles=5)
                prof = circle_prof.values
                # determine if the profile contains both low and high values and that most values are the same
                if (np.percentile(prof, 2) < 800) and (np.percentile(prof, 98) > 800) and (
                        np.percentile(prof, 80) - np.percentile(prof, 30) < 40):
                    hu_slices.append(image_number)

        center_hu_slice = int(np.median(hu_slices))
        if self._is_within_image_extent(center_hu_slice):
            return center_hu_slice

    @property
    def un_slice_num(self):
        """The HU uniformity slice number."""
        slice = int(self.hu_slice_num - round(73/self.dicom_stack.metadata.SliceThickness))
        if self._is_within_image_extent(slice):
            return slice

    @property
    def sr_slice_num(self):
        """The Spatial Resolution slice number."""
        slice = int(self.hu_slice_num + round(30/self.dicom_stack.metadata.SliceThickness))
        if self._is_within_image_extent(slice):
            return slice

    @property
    def lc_slice_num(self):
        """The low contrast slice number."""
        slice = int(self.hu_slice_num - round(30/self.dicom_stack.metadata.SliceThickness))
        if self._is_within_image_extent(slice):
            return slice

    @property
    @lru_cache()
    def phantom_roll(self):
        """Determine the "roll" of the phantom.

         This algorithm uses the two air bubbles in the HU slice and the resulting angle between them.

        Returns
        -------
        float : the angle of the phantom in **degrees**.
        """
        slice_of_interest = Image.from_array(self.dicom_stack.slice(self.hu_slice_num)).threshold(self.threshold)
        slice_of_interest.invert()
        labels, no_roi = ndimage.measurements.label(slice_of_interest.array)
        # calculate ROI sizes of each label TODO: simplify the air bubble-finding
        roi_sizes = [ndimage.measurements.sum(slice_of_interest.array, labels, index=item) for item in range(1, no_roi + 1)]
        # extract air bubble ROIs (based on size threshold)
        bubble_thresh = self.air_bubble_size
        air_bubbles = [idx + 1 for idx, item in enumerate(roi_sizes) if
                       item < bubble_thresh * 1.5 and item > bubble_thresh / 1.5]
        # if the algo has worked correctly, it has found 2 and only 2 ROIs (the air bubbles)
        if len(air_bubbles) == 2:
            air_bubble_CofM = ndimage.measurements.center_of_mass(slice_of_interest.array, labels, air_bubbles)
            y_dist = air_bubble_CofM[0][0] - air_bubble_CofM[1][0]
            x_dist = air_bubble_CofM[0][1] - air_bubble_CofM[1][1]
            angle = np.arctan2(y_dist, x_dist)
            if angle < 0:
                roll = abs(angle) - np.pi / 2
            else:
                roll = angle - np.pi / 2
            phan_roll = roll
        else:
            phan_roll = 0
            print("Warning: CBCT phantom roll unable to be determined; assuming 0")
        return np.rad2deg(phan_roll)

    @property
    def expected_phantom_size(self):
        """The expected size of the phantom in pixels, based on a 20cm wide phantom."""
        phan_area = np.pi*(101**2)  # Area = pi*r^2; slightly larger value used based on actual values acquired
        return phan_area/(self.mm_per_pixel**2)

    @property
    def num_images(self):
        """Return the number of images loaded."""
        return self.dicom_stack.shape[-1]

    def _is_within_image_extent(self, image_num):
        """Determine if the image number is beyond the edges of the images (negative or past last image)."""
        if self.num_images - 1 > image_num > 1:
            return True
        else:
            raise ValueError("The determined image number is beyond the image extent. Either the entire dataset "
                             "wasn't loaded or the entire phantom wasn't scanned.")


class RectangleROI(Rectangle):
    """Class that represents a rectangular ROI."""
    def __init__(self, array, width, height, angle, dist_from_center, phantom_center):
        y_shift = -np.sin(np.deg2rad(angle)) * dist_from_center
        x_shift = np.cos(np.deg2rad(angle)) * dist_from_center
        center = Point(phantom_center.x + x_shift, phantom_center.y + y_shift)
        super().__init__(width, height, center, as_int=True)
        self._array = array

    @property
    def pixel_array(self):
        """The pixel array within the ROI."""
        return self._array[self.bl_corner.x:self.tr_corner.x, self.bl_corner.y:self.tr_corner.y]


class DiskROI(Circle):
    """An class representing a disk-shaped Region of Interest on a CBCT slice."""
    def __init__(self, array, angle, roi_radius, dist_from_center, phantom_center):
        """
        Parameters
        ----------
        array : ndarray
            The 2D array representing the image the disk is on.
        angle : int, float
            The angle of the ROI in degrees from the phantom center.
        roi_radius : int, float
            The radius of the ROI from the center of the phantom.
        dist_from_center : int, float
            The distance of the ROI from the phantom center.
        phantom_center : tuple
            The location of the phantom center.
        """
        center = self._get_shifted_center(angle, dist_from_center, phantom_center)
        super().__init__(center_point=center, radius=roi_radius)
        self._array = array

    def _get_shifted_center(self, angle, dist_from_center, phantom_center):
        """The center of the ROI; corrects for phantom dislocation and roll."""
        y_shift = -np.sin(np.deg2rad(angle)) * dist_from_center
        x_shift = np.cos(np.deg2rad(angle)) * dist_from_center
        return Point(phantom_center.x + x_shift, phantom_center.y + y_shift)

    @property
    def pixel_value(self):
        """The median pixel value of the ROI."""
        masked_img = self._get_roi_mask()
        return np.nanmedian(masked_img)

    @property
    def std(self):
        """The standard deviation of the pixel values."""
        masked_img = self._get_roi_mask()
        return np.nanstd(masked_img)

    def _get_roi_mask(self):
        """Return a masked array of the ROI."""
        mask = sector_mask(self._array.shape, self.center, self.radius)
        masked_img = np.where(mask == True, self._array, np.NaN)
        return masked_img


class HUDiskROI(DiskROI):
    """An HU ROI object. Represents a circular area measuring either HU sample (Air, Poly, ...)
    or HU uniformity (bottom, left, ...).
    """
    def __init__(self, array, angle, roi_radius, dist_from_center, phantom_center, nominal_value=None, tolerance=None):
        """
        Parameters
        ----------
        nominal_value : int
            The nominal pixel value of the HU ROI.
        tolerance : int
            The roi pixel value tolerance.
        """
        super().__init__(array, angle, roi_radius, dist_from_center, phantom_center)
        self.nominal_val = nominal_value
        self.tolerance = tolerance

    @property
    def value_diff(self):
        """The difference in HU between measured and nominal."""
        return abs(self.pixel_value - self.nominal_val)

    @property
    def passed(self):
        """Boolean specifying if ROI pixel value was within tolerance of the nominal value."""
        return self.value_diff <= self.tolerance

    @property
    def plot_color(self):
        """Return one of two colors depending on if ROI passed."""
        return 'blue' if self.passed else 'red'


class LowContrastDiskROI(DiskROI):
    """A class for analyzing the low-contrast disks."""
    def __init__(self, array, angle, roi_radius, dist_from_center, phantom_center, contrast_threshold):
        """
        Parameters
        ----------
        contrast_threshold : float, int
            The threshold for considering a bubble to be "seen".
        """
        super().__init__(array, angle, roi_radius, dist_from_center, phantom_center)
        self.contrast_threshold = contrast_threshold

    @property
    def contrast(self):
        """The true contrast of the bubble: (Signal - Background)/Stdev."""
        return (self.pixel_value - self.background)/self.std

    @property
    def contrast_constant(self):
        """The contrast value times the bubble diameter."""
        return self.contrast*self.diameter

    @property
    def passed(self):
        """Boolean specifying if ROI pixel value was within tolerance of the nominal value."""
        return self.contrast_constant > self.contrast_threshold

    @property
    def plot_color(self):
        """Return one of two colors depending on if ROI passed."""
        return 'blue' if self.passed else 'red'


class ThicknessROI(RectangleROI):
    """A rectangular ROI that measures the angled wire rod in the HU linearity slice which determines slice thickness."""

    @property
    def long_profile(self):
        """The profile along the axis perpendicular to ramped wire."""
        img = Image.from_array(self.pixel_array)
        img.median_filter()
        prof = SingleProfile(img.array.max(axis=np.argmin(img.shape)))
        prof.filter(0.05)
        return prof

    @property
    @lru_cache()
    def wire_fwhm(self):
        """The FWHM of the wire in pixels."""
        prof = self.long_profile
        return prof.fwxm(interpolate=True)

    @property
    def plot_color(self):
        """The plot color."""
        return 'blue'


class ROIManagerMixin:
    """Class for handling multiple ROIs. Used for the HU linearity, Uniformity, Geometry, Low-contrast, and Thickness slices.

    Attributes
    ----------
    dist2rois_mm : int, float
        The distance from the phantom center to the ROIs, in mm.
    roi_radius_mm : int, float
        The radius of the ROIs, in mm.
    roi_names : list
        The names of the ROIs.
    roi_nominal_angles : list
        The nominal angles of the ROIs; must be same order as ``roi_names``.
    """
    dist2rois_mm = 0
    roi_radius_mm = 0
    roi_names = []
    roi_nominal_angles = []

    @abstractmethod
    def _setup_rois(self):
        pass

    @property
    def roi_angles(self):
        """The ROI angles, corrected for phantom roll."""
        return np.array(self.roi_nominal_angles) + self.settings.phantom_roll

    @property
    def dist2rois(self):
        """Distance from the phantom center to the ROIs, corrected for pixel spacing."""
        return self.dist2rois_mm / self.settings.mm_per_pixel

    @property
    def roi_radius(self):
        """ROI radius, corrected for pixel spacing."""
        return self.roi_radius_mm / self.settings.mm_per_pixel

    def get_ROI_vals(self):
        """Return a dict of the HU values of the HU ROIs."""
        return {key: val.pixel_value for key, val in self.rois.items()}

    def plot_rois(self, axis):
        """Plot the ROIs to the axis."""
        for roi in self.rois.values():
            roi.add_to_axes(axis, edgecolor=roi.plot_color)


class Slice:
    """Base class for analyzing specific slices of a CBCT dicom set."""
    def __init__(self, dicom_stack, settings, slice_num=None, combine=True, combine_method='mean', num_slices=1):
        """
        Parameters
        ----------
        dicom_stack : :class:`~pylinac.core.image.DICOMStack`
        settings : :class:`~pylinac.cbct.Settings`
        slice_num : int
            The slice number of the DICOM array desired. If None, will use the ``slice_num`` property of subclass.
        combine : bool
            If True, combines the slices +/- ``num_slices`` around the slice of interest to improve signal/noise.
        combine_method : {'mean', 'max'}
            How to combine the slices if ``combine`` is True.
        num_slices : int
            The number of slices on either side of the nominal slice to combine to improve signal/noise; only
            applicable if ``combine`` is True.
        """
        self.settings = settings
        if slice_num is not None:
            self.slice_num = slice_num
        if combine:
            array = combine_surrounding_slices(dicom_stack.array, self.slice_num, mode=combine_method, slices_plusminus=num_slices)
        else:
            array = dicom_stack.slice(self.slice_num)
        self.image = Image.from_array(array)

    @property
    def __getitem__(self, item):
        return self.image.array[item]

    @property
    @lru_cache()
    def phan_center(self):
        """Determine the location of the center of the phantom.

        The image is analyzed to see if:
        1) the catphan is even in the image (if there were any ROIs detected)
        2) an ROI is within the size criteria of the catphan
        3) the ROI area that is filled compared to the bounding box area is close to that of a circle

        Raises
        ------
        ValueError
            If any of the above conditions are not met.
        """
        SOI_bw = self.image.threshold(self.settings.threshold)  # convert slice to binary based on threshold
        SOI_labeled, num_roi = ndimage.label(SOI_bw)  # identify the ROIs
        if num_roi < 1 or num_roi is None:
            raise ValueError("Unable to locate the CatPhan")
        roi_sizes, bin_edges = np.histogram(SOI_labeled, bins=num_roi+1)  # hist will give the size of each label
        rois_in_size_criteria = [self.settings.expected_phantom_size * 0.98 < roi_size < self.settings.expected_phantom_size * 1.02 for roi_size in roi_sizes]
        if not any(rois_in_size_criteria):
            raise ValueError("Unable to locate the CatPhan")
        else:
            catphan_idx = np.abs(roi_sizes - self.settings.expected_phantom_size).argmin()
        SOI_bw_clean = np.where(SOI_labeled == catphan_idx, 1, 0)  # remove all ROIs except CatPhan
        expected_fill_ratio = np.pi / 4  # the area of a circle divided by its bounding box
        actual_fill_ratio = get_filled_area_ratio(SOI_bw_clean)
        if (expected_fill_ratio * 1.02 < actual_fill_ratio) or (actual_fill_ratio < expected_fill_ratio * 0.98):
            raise ValueError("Unable to locate the CatPhan")
        center_pixel = get_center_of_mass(SOI_bw_clean)
        return Point(center_pixel[1], center_pixel[0])


class HUSlice(Slice, ROIManagerMixin):
    """Class for analysis of the HU linearity slice of the CBCT dicom data set. Although the same CBCT slice is used for other
    analyses (Geometry, Thickness), this class is only interested in the linearity ROIs. The class analyzes 7 of the
    ROIs (skipping the 2nd air ROI); these values should be linear.

    Attributes
    ----------
    roi_nominal_values : list
        A list of the nominal HU values for the ROIs; in same order as ``roi_names``.
    """
    dist2rois_mm = 58.7
    roi_radius_mm = 5
    roi_names = ['Air', 'PMP', 'LDPE', 'Poly', 'Acrylic', 'Delrin', 'Teflon']
    roi_nominal_values = [-1000, -200, -100, -35, 120, 340, 990]
    roi_nominal_angles = [90, 120, 180, -120, -60, 0, 60]

    def __init__(self, dicom_stack, settings):
        super().__init__(dicom_stack, settings)
        self._setup_rois()

    def _setup_rois(self):
        self.rois = OrderedDict()
        for name, angle, nominal_value in zip(self.roi_names, self.roi_angles, self.roi_nominal_values):
            self.rois[name] = HUDiskROI(self.image, angle, self.roi_radius, self.dist2rois,
                                        self.phan_center, nominal_value, self.tolerance)

    @property
    def slice_num(self):
        """The CBCT slice number for the HU linearity ROIs."""
        return self.settings.hu_slice_num

    @property
    def tolerance(self):
        """The tolerance of the HU linearity ROIs."""
        return self.settings.hu_tolerance

    def plot_linearity(self, axis=None, plot_delta=True):
        """Plot the HU linearity values to an axis.

        Parameters
        ----------
        axis : None, matplotlib.Axes
            The axis to plot the values on. If None, will create a new figure.
        plot_delta : bool
            Whether to plot the actual measured HU values (False), or the difference from nominal (True).
        """
        if axis is None:
            fig, axis = plt.subplots()
        if plot_delta:
            values = [roi.value_diff for roi in self.rois.values()]
            nominal_values = [0]*len(values)
            ylabel = 'HU Delta'
        else:
            values = [roi.pixel_value for roi in self.rois.values()]
            nominal_values = self.roi_nominal_values
            ylabel = 'Measured Values'
        points = axis.plot(self.roi_nominal_values, values, 'g+', markersize=15, mew=2)
        axis.plot(self.roi_nominal_values, nominal_values)
        axis.plot(self.roi_nominal_values, np.array(nominal_values) + self.tolerance, 'r--')
        axis.plot(self.roi_nominal_values, np.array(nominal_values) - self.tolerance, 'r--')
        axis.margins(0.05)
        axis.grid('on')
        axis.set_xlabel("Nominal Values")
        axis.set_ylabel(ylabel)
        return points

    @property
    def overall_passed(self):
        """Boolean specifying whether all the ROIs passed within tolerance."""
        return all(roi.passed for roi in self.rois.values())


class UniformitySlice(Slice, ROIManagerMixin):
    """Class for analysis of the Uniformity slice of the CBCT dicom data set. Measures 5 ROIs around the slice that
    should all be close to the same value.
    """
    dist2rois_mm = 53
    roi_radius_mm = 10
    roi_names = ['Top', 'Right', 'Bottom', 'Left', 'Center']
    roi_nominal_values = [0, 0, 0, 0, 0]
    roi_nominal_angles = [-90, 0, 90, 180, 0]

    def __init__(self, dicom_stack, settings):
        super().__init__(dicom_stack, settings)
        self._setup_rois()

    def _setup_rois(self):
        self.rois = OrderedDict()
        for name, angle, nominal_value in zip(self.roi_names, self.roi_angles, self.roi_nominal_values):
            distance = self.dist2rois if name != 'Center' else 0
            self.rois[name] = HUDiskROI(self.image, angle, self.roi_radius, distance,
                                        self.phan_center, nominal_value, self.tolerance)

    @property
    def slice_num(self):
        """The CBCT slice number for the HU uniformity ROIs."""
        return self.settings.un_slice_num

    @property
    def tolerance(self):
        """The tolerance of the HU uniformity ROIs."""
        return self.settings.hu_tolerance

    def plot_profiles(self, axis=None):
        """Plot the horizontal and vertical profiles of the Uniformity slice.

        Parameters
        ----------
        axis : None, matplotlib.Axes
            The axis to plot on; if None, will create a new figure.
        """
        if axis is None:
            fig, axis = plt.subplots()
        horiz_data = self.image[int(self.phan_center.y), :]
        vert_data = self.image[:, int(self.phan_center.x)]
        axis.plot(horiz_data, 'g', label='Horizontal')
        axis.plot(vert_data, 'b', label='Vertical')
        axis.autoscale(tight=True)
        axis.set_ylim([-100, 100])
        # TODO: replace .plot() calls with .axhline() calls when mpld3 fixes functionality
        axis.plot([i for i in range(len(horiz_data))], [self.tolerance] * len(horiz_data), 'r-', linewidth=3)
        axis.plot([i for i in range(len(horiz_data))], [-self.tolerance] * len(horiz_data), 'r-', linewidth=3)
        axis.grid('on')
        axis.set_ylabel("HU")
        axis.legend(loc=8, fontsize='small', title="")

    @property
    def overall_passed(self):
        """Boolean specifying whether all the ROIs passed within tolerance."""
        return all(roi.passed for roi in self.rois.values())


class LowContrastSlice(Slice, ROIManagerMixin):
    """Class for analysis of the low contrast slice of the CBCT dicom data set. Low contrast is measured by obtaining
    the average pixel value of the contrast ROIs and comparing that value to the average background value. To obtain
    a more "human" detection level, the contrast (which is largely the same across different-sized ROIs) is multiplied
    by the diameter. This value is compared to the contrast threshold to decide if it can be "seen".

    Attributes
    ----------
    roi_background_names : list
        A list identifying which of the ``roi_names`` are the background ROIs.

    """
    dist2rois_mm = 50
    roi_radius_mm = [7, 7, 4.5, 4, 3.5, 3, 2.5, 7]
    roi_names = ['bg', '15', '9', '8', '7', '6', '5', 'bg2']
    roi_background_names = ['bg', 'bg2']
    roi_nominal_angles = [110, 89, 68.6, 52.4, 37.6, 26.6, 12.9, -10]

    def __init__(self, dicom_stack, settings):
        super().__init__(dicom_stack, settings)
        self._setup_rois()

    def _setup_rois(self):
        self.rois = OrderedDict()
        self.bg_rois = {}
        for name, angle, radius in zip(self.roi_names, self.roi_angles, self.roi_radius):
            if name not in self.roi_background_names:
                self.rois[name] = LowContrastDiskROI(self.image, angle, radius, self.dist2rois,
                                                     self.phan_center, self.settings.contrast_threshold)
            else:
                self.bg_rois[name] = LowContrastDiskROI(self.image, angle, radius, self.dist2rois,
                                                        self.phan_center, self.settings.contrast_threshold)
        for roi in self.rois.values():
            roi.background = np.mean([roi.pixel_value for roi in self.bg_rois.values()])

    @property
    def tolerance(self):
        return self.settings.low_contrast_tolerance

    @property
    def slice_num(self):
        return self.settings.lc_slice_num

    @property
    def roi_radius(self):
        """A list of the ROI radii, scaled to pixels."""
        return [radius / self.settings.mm_per_pixel for radius in self.roi_radius_mm]

    @property
    def rois_visible(self):
        """The number of ROIs "visible"."""
        return sum(roi.passed for roi in self.rois.values())

    def plot_rois(self, axis):
        """Plot the ROIs to an axis."""
        super().plot_rois(axis)
        for roi in self.bg_rois.values():
            roi.add_to_axes(axis, 'blue')

    @property
    def overall_passed(self):
        """Whether there were enough low contrast ROIs "seen"."""
        return sum(roi.passed for roi in self.rois.values()) >= self.tolerance


class SpatialResolutionSlice(Slice):
    """Class for analysis of the Spatial Resolution slice of the CBCT dicom data set.

    A collapsed circle profile is taken of the line-pair region. This profile is search for
    peaks and valleys. The MTF is calculated from those peaks & valleys.

    Attributes
    ----------
    line_pair_frequency : tuple
        The frequency of line pairs for the first 6 regions.
    num_peaks : array
        The index of the peaks for the various regions.
    num_valleys : array
        The index of the valleys for the various regions.
    radius2linepairs_mm : float
        The radius in mm to the line pairs.
    line_pair_cutoff : int
        The approximate index to stop analyzing the line profile (higher indices are too hard to analyze).
    """
    line_pair_frequency = (0.2, 0.4, 0.6, 0.8, 1, 1.2)
    # num_peaks = np.array((0, 2, 3, 3, 4, 4, 4)).cumsum()
    # num_valleys = np.array((0, 1, 2, 2, 3, 3, 3)).cumsum()
    radius2linepairs_mm = 47.5
    line_pair_cutoff = 0.34

    def __init__(self, *args, **kwargs):
        super().__init__(combine_method='max', *args, **kwargs)

    @property
    def slice_num(self):
        """The slice number of the spatial resolution module."""
        return self.settings.sr_slice_num

    @property
    def num_peaks(self):
        return np.array((0, 2, 3, 3, 4, 4, 4)).cumsum()

    @property
    def num_valleys(self):
        return np.array((0, 1, 2, 2, 3, 3, 3)).cumsum()


    @property
    def radius2linepairs(self):
        """Radius from the phantom center to the line-pair region, corrected for pixel spacing."""
        return self.radius2linepairs_mm / self.settings.mm_per_pixel

    def plot_rois(self, axis):
        """Plot the circles where the profile was taken within."""
        self.circle_profile.add_to_axes(axis, edgecolor='blue')

    @property
    @lru_cache()
    def circle_profile(self):
        """Calculate the median profile of the Line Pair region.

        Returns
        -------
        `pylinac.core.profile.CollapsedCircleProfile` : A 1D profile of the Line Pair region.
        """
        circle_profile = CollapsedCircleProfile(self.phan_center, self.radius2linepairs, image_array=self.image,
                                                start_angle=np.pi + np.deg2rad(self.settings.phantom_roll),
                                                width_ratio=0.04, sampling_ratio=2)
        circle_profile.filter(0.001, kind='gaussian')
        return circle_profile

    @property
    @lru_cache()
    def spaced_circle_profile(self):
        """The median profile of the line pair region with equalized spacing applied.

        The line pairs change in spacing, and this function stretches out the later sections to be similar
        in peak spacing as those of the larger areas.
        """
        line_cutoff = int(self.line_pair_cutoff * len(self.circle_profile.values))
        spacing_array = np.linspace(1, 12, num=line_cutoff, dtype=int)
        spaced_array = np.repeat(self.circle_profile.values[:line_cutoff], spacing_array)
        profile = MultiProfile(spaced_array)
        profile.ground()
        return profile

    @property
    @lru_cache()
    def profile_peaks_idxs(self):
        """Return the peak values and indices of the spaced-out profile."""
        max_idxs = self.spaced_circle_profile.find_peaks(min_distance=0.025, max_number=17, threshold=0.05)
        if len(max_idxs) != 17:
            # TODO: add some robustness here
            raise ArithmeticError("Did not find the correct number of line pairs")
        return max_idxs

    @property
    @lru_cache()
    def profile_peaks_vals(self):
        """Return the peak values and indices of the spaced-out profile."""
        max_vals = self.spaced_circle_profile.find_peaks(min_distance=0.025, max_number=17, threshold=0.05, kind='value')
        if len(max_vals) != 17:
            raise ArithmeticError("Did not find the correct number of line pairs")
        return max_vals

    @property
    @lru_cache()
    def profile_valleys_idxs(self):
        """Return the valley values and indices of the spaced-out profile,
            with valleys in-between line-pair regions removed."""
        idx2del = np.array((1, 4, 8, 12))
        min_idxs = np.zeros(16)
        max_idxs = sorted(self.profile_peaks_idxs)
        for idx in range(len(max_idxs) - 1):
            min_idx = self.spaced_circle_profile.find_valleys(search_region=(int(max_idxs[idx]),
                                                                             int(max_idxs[idx + 1])),
                                                              max_number=1)
            if len(min_idx) > 0:
                min_idxs[idx] = min_idx[0]
        min_idxs = np.delete(min_idxs, idx2del)
        return min_idxs

    @property
    @lru_cache()
    def profile_valleys_vals(self):
        """Return the valley values and indices of the spaced-out profile,
            with valleys in-between line-pair regions removed."""
        idx2del = np.array((1, 4, 8, 12))
        min_vals = np.zeros(16)
        max_idxs = sorted(self.profile_peaks_idxs)
        for idx in range(len(max_idxs) - 1):
            min_idx = self.spaced_circle_profile.find_valleys(search_region=(int(max_idxs[idx]),
                                                                             int(max_idxs[idx + 1])),
                                                              max_number=1, kind='value')
            if len(min_idx) > 0:
                min_vals[idx] = min_idx[0]
        min_vals = np.delete(min_vals, idx2del)
        return min_vals

    def mtf(self, percent=None, region=None):
        """Return the MTF value of the spatial resolution. Only one of the two parameters may be used.

        Parameters
        ----------
        percent : int, float
            The percent relative MTF; i.e. 0-100.
        region : int
            The line-pair region desired (1-6).

        Returns
        -------
        float : the line-pair resolution at the given MTF percent or region.
        """
        if (region is None and percent is None) or (region is not None and percent is not None):
            raise ValueError("Must pass in either region or percent")
        if percent is not None:
            x_vals_intrp = np.arange(self.line_pair_frequency[0], self.line_pair_frequency[-1], 0.01)
            x_vals = np.array(self.line_pair_frequency)
            y_vals = self.line_pair_mtfs
            y_vals_intrp = np.interp(x_vals_intrp, x_vals, y_vals)
            # TODO: warn user if value at MTF edge; may not be true MTF
            mtf_percent = x_vals_intrp[np.argmin(np.abs(y_vals_intrp - (percent / 100)))]
            return simple_round(mtf_percent, 2)
        elif region is not None:
            return self.line_pair_mtfs[region]

    @property
    @lru_cache()
    def line_pair_mtfs(self):
        """The discrete MTF values at the line-pair regions."""
        num_peaks = self.num_peaks
        num_valleys = self.num_valleys
        mtfs = []
        for frequency, region in zip(self.line_pair_frequency, range(len(num_peaks) - 1)):
            region_max = self.profile_peaks_vals[num_peaks[region]:num_peaks[region + 1]].mean()
            region_min = self.profile_valleys_vals[num_valleys[region]:num_valleys[region + 1]].mean()
            mtfs.append((region_max - region_min) / (region_max + region_min))
        # normalize the values by the first LP
        mtfs = np.array(mtfs) / max(mtfs)
        return mtfs

    def plot_mtf(self, axis=None):
        """Plot the Relative MTF.

        Parameters
        ----------
        axis : None, matplotlib.Axes
            The axis to plot the MTF on. If None, will create a new figure.
        """
        if axis is None:
            fig, axis = plt.subplots()
        line_pairs = self.line_pair_frequency
        mtf_vals = self.line_pair_mtfs
        points = axis.plot(line_pairs, mtf_vals, marker='o')
        axis.margins(0.05)
        axis.grid('on')
        axis.set_xlabel('Line pairs / mm')
        axis.set_ylabel("Relative MTF")
        return points


class ThicknessSlice(Slice, ROIManagerMixin):
    """This class analyzes the angled wire on the HU slice to determine the slice thickness.

    Attributes
    ----------
    roi_widths_mm : list
        The widths of the rectangular ROIs in mm. Follows the order of ``roi_names``.
    roi_heights_mm : list
        The heights of the rectangular ROIs in mm. Follows the order of ``roi_names``.
    """
    roi_names = ['Left', 'Top', 'Right', 'Bottom']
    roi_nominal_angles = [180, 90, 0, -90]
    roi_widths_mm = [8, 40, 8, 40]
    roi_heights_mm = [40, 8, 40, 8]
    dist2rois_mm = 38

    def __init__(self, dicom_stack, settings):
        super().__init__(dicom_stack, settings)
        self._setup_rois()

    def _setup_rois(self):
        self.rois = OrderedDict()
        for name, angle, height, width in zip(self.roi_names, self.roi_angles, self.roi_heights, self.roi_widths):
            self.rois[name] = ThicknessROI(self.image, width, height, angle, self.dist2rois, self.phan_center)

    @property
    def tolerance(self):
        """Tolerance of the slice thickness."""
        return self.settings.thickness_tolerance

    @property
    def slice_num(self):
        """Slice number of the thickness slice."""
        return self.settings.hu_slice_num

    @property
    def avg_slice_thickness(self):
        """The average slice thickness for the 4 wire measurements in mm."""
        return np.mean(sorted(roi.wire_fwhm*self.settings.mm_per_pixel*0.42 for roi in self.rois.values())[-2:])/3

    @property
    def nominal_slice_thickness(self):
        """The nominal slice thickness from the DICOM data."""
        return self.settings.dicom_stack.metadata.SliceThickness

    @property
    def passed(self):
        """Whether the slice thickness was within tolerance from nominal."""
        return self.nominal_slice_thickness - self.tolerance < self.avg_slice_thickness < self.nominal_slice_thickness + self.tolerance

    @property
    def roi_heights(self):
        """The ROI heights, corrected for pixel spacing."""
        return [height / self.settings.mm_per_pixel for height in self.roi_heights_mm]

    @property
    def roi_widths(self):
        """The ROI widths, corrected for pixel spacing."""
        return [width / self.settings.mm_per_pixel for width in self.roi_widths_mm]


class GeoDiskROI(DiskROI):
    """A circular ROI, much like the HU ROI, but with methods to find the center of the geometric "node"."""

    def _threshold_node(self):
        """Threshold the ROI to find node.

        Three of the four nodes have a positive value, while one node is air and
        thus has a low value. The algorithm thus thresholds for extreme values relative
        to the median value (which is the node).

        Returns
        -------
        bw_node : numpy.array
            A masked 2D array the size of the Slice image, where only the node pixels have a value.
        """
        masked_img = self._get_roi_mask()
        # threshold image
        upper_band_pass = np.where(masked_img > np.nanmedian(masked_img) * 1.4, 1, 0)
        lower_band_pass = np.where(masked_img < np.nanmedian(masked_img) * 0.6, 1, 0)
        bw_node = upper_band_pass + lower_band_pass
        return bw_node

    @property
    @lru_cache()
    def node_center(self):
        """Find the center of the geometric node within the ROI."""
        bw_node = self._threshold_node()
        # label ROIs found
        labeled_arr, num_roi = ndimage.measurements.label(bw_node)
        roi_sizes, bin_edges = np.histogram(labeled_arr, bins=num_roi+1)  # hist will give the size of each label
        bw_node_cleaned = np.where(labeled_arr == np.argsort(roi_sizes)[-2], 1, 0)  # remove all ROIs except the second largest one (largest one is the air itself)
        labeled_arr, num_roi = ndimage.measurements.label(bw_node_cleaned)
        if num_roi != 1:
            raise ValueError("Did not find the geometric node.")
        # determine the center of mass of the geometric node
        x_arr = np.abs(np.average(bw_node_cleaned, weights=self._array, axis=0))
        x_com = SingleProfile(x_arr).fwxm_center(interpolate=True)
        y_arr = np.abs(np.average(bw_node_cleaned, weights=self._array, axis=1))
        y_com = SingleProfile(y_arr).fwxm_center(interpolate=True)
        return Point(x_com, y_com)


class GeometricLine(Line):
    """Represents a line connecting two nodes/ROIs on the Geometry Slice.

    Attributes
    ----------
    nominal_length_mm : int, float
        The nominal distance between the geometric nodes, in mm.
    """
    nominal_length_mm = 50

    def __init__(self, geo_roi1, geo_roi2, mm_per_pixel, tolerance):
        """
        Parameters
        ----------
        geo_roi1 : GEO_ROI
            One of two ROIs representing one end of the line.
        geo_roi2 : GEO_ROI
            The other ROI which is the other end of the line.
        mm_per_pixel : float
            The mm/pixel value.
        tolerance : int, float
            The tolerance of the geometric line, in mm.
        """
        super().__init__(geo_roi1.node_center, geo_roi2.node_center)
        self.mm_per_pixel = mm_per_pixel
        self.tolerance = tolerance

    @property
    def passed(self):
        """Whether the line passed tolerance."""
        return self.nominal_length_mm - self.tolerance < self.length_mm < self.nominal_length_mm + self.tolerance

    @property
    def pass_fail_color(self):
        """Plot color for the line, based on pass/fail status."""
        return 'blue' if self.passed else 'red'

    @property
    def length_mm(self):
        """Return the length of the line in mm."""
        return self.length*self.mm_per_pixel


class GeometrySlice(Slice, ROIManagerMixin):
    """Class for analysis of the Geometry slice of the CBCT dicom data set.

    Four ROIs are set, which correspond to the locations of the 1 air and 3 acrylic "nodes".
    Within these ROIs the center of the nodes must be found.

    Once the nodes centers are found four lines are constructed by linking the node centers.

    Attributes
    ----------
    line_assignments : dict
        The assignment of line -> nodes.
    lines : dict
        Holds `~pylinac.cbct.GeometricLine` instances.
    """
    line_assignments = {'Top-Horizontal': ('Top-Left', 'Top-Right'),
                        'Bottom-Horizontal': ('Bottom-Left', 'Bottom-Right'),
                        'Left-Vertical': ('Top-Left', 'Bottom-Left'),
                        'Right-Vertical': ('Top-Right', 'Bottom-Right')}
    lines = {}
    dist2rois_mm = 35
    roi_radius_mm = 6
    roi_names = ['Top-Left', 'Top-Right', 'Bottom-Right', 'Bottom-Left']
    roi_nominal_angles = [-135, -45, 45, 135]

    def __init__(self, dicom_stack, settings):
        super().__init__(dicom_stack, settings)
        # apply a filter to reduce salt & pepper noise around nodes
        self.image.median_filter(size=3)
        self._setup_rois()

    def _setup_rois(self):
        self.rois = OrderedDict()
        for name, angle in zip(self.roi_names, self.roi_angles):
            self.rois[name] = GeoDiskROI(self.image, angle, self.roi_radius, self.dist2rois, self.phan_center)
        # setup the geometric lines
        for name, (node1, node2) in self.line_assignments.items():
            self.lines[name] = GeometricLine(self.rois[node1], self.rois[node2], self.settings.mm_per_pixel,
                                             self.settings.scaling_tolerance)

    @property
    def tolerance(self):
        """Tolerance of the geometric lines."""
        return self.settings.scaling_tolerance

    @property
    def slice_num(self):
        """Geometric slice number."""
        return self.settings.hu_slice_num

    def plot_lines(self, axis):
        """Plot the geometric node connection lines."""
        for line in self.lines.values():
            line.add_to_axes(axis, color=line.pass_fail_color)

    @property
    def avg_line_length(self):
        return np.mean([line.length_mm for line in self.lines.values()])

    @property
    def overall_passed(self):
        """Teturns whether all the line lengths were within tolerance."""
        return all(line.passed for line in self.lines.values())


def get_bounding_box(array):
    """Get the bounding box values of an ROI in a 2D array."""
    B = np.argwhere(array)
    (ymin, xmin), (ymax, xmax) = B.min(0), B.max(0) + 1
    return ymin, ymax, xmin, xmax

def get_center_of_mass(array):
    """Return the center of mass of the ROI."""
    return ndimage.measurements.center_of_mass(array)  # returns (y,x)

def get_filled_area_ratio(array):
    """Return the ratio of the bounding box area and the area actually filled."""
    ymin, ymax, xmin, xmax = get_bounding_box(array)
    box_area = (ymax - ymin) * (xmax - xmin)
    filled_area = np.sum(array)
    return filled_area/box_area

@value_accept(mode=('mean','median','max'))
def combine_surrounding_slices(slice_array, nominal_slice_num, slices_plusminus=1, mode='mean'):
    """Return an array that is the combination of a given slice and a number of slices surrounding it.

    Parameters
    ----------
    im_array : numpy.array
        The original 3D numpy array of images.
    nominal_slice_num : int
        The slice of interest (along 3rd dim).
    slices_plusminus: int
        How many slices plus and minus to combine (also along 3rd dim).
    mode : {'mean', 'median', 'max}
        Specifies the method of combination.

    Returns
    -------
    comb_slice : numpy.array
        An array the same size in the first two dimensions of im_array, combined.
    """
    slices = slice_array[:,:,nominal_slice_num-slices_plusminus:nominal_slice_num+slices_plusminus+1]
    if mode == 'mean':
        comb_slice = np.mean(slices, 2)
    elif mode == 'median':
        comb_slice = np.median(slices, 2)
    else:
        comb_slice = np.max(slices, 2)
    return comb_slice


# ----------------------------------------
# CBCT Demo
# ----------------------------------------
if __name__ == '__main__':
    CBCT().run_demo()

