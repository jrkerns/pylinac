# -*- coding: utf-8 -*-
"""The VMAT module consists of the class VMAT, which is capable of loading an EPID DICOM Open field image and MLC field image and analyzing the
images according to the Varian RapidArc QA tests and procedures, specifically the Dose-Rate & Gantry-Speed (DRGS)
and Dose-Rate & MLC speed (DRMLC) tests.

Features:

* **Do both tests** - Pylinac can handle either DRGS or DRMLC tests.
* **Adjust for offsets** - Older VMAT patterns were off-center. Easily account for the offset by passing it in.
* **Automatic identification using file names** - If your file names are clear, the image type and test type don't even
  have to be specified; just load and analyze.
"""
from io import BytesIO
import os.path as osp
import zipfile

import matplotlib.pyplot as plt
import numpy as np

from pylinac.core.decorators import value_accept, type_accept
from pylinac.core.geometry import Point, Rectangle
from pylinac.core.image import Image
from pylinac.core.io import get_filepath_UI, get_filenames_UI
from pylinac.core.profile import SingleProfile
from pylinac.core.utilities import typed_property, get_url, import_mpld3, load_zipfile

# test types
DRGS = 'drgs'
DRMLC = 'drmlc'
# image types
DMLC = 'dmlc'
OPEN = 'open'
PROFILE = 'profile'

# ROI settings according to Varian
DRGS_SETTINGS = {'X-plane offsets (mm)': (-60, -40, -20, 0, 20, 40, 60)}
MLCS_SETTINGS = {'X-plane offsets (mm)': (-45, -15, 15, 45)}


class VMAT:
    """The VMAT class analyzes two DICOM images acquired via a linac's EPID and analyzes
        regions of interest (segments) based on the Varian RapidArc QA specifications,
        specifically, the Dose Rate & Gantry Speed (DRGS) and Dose Rate & MLC speed (DRMLC) tests.

        Attributes
        ----------
        image_open : :class:`~pylinac.core.image.Image`
        image_dmlc : :class:`~pylinac.core.image.Image`
        segments : :class:`~pylinac.vmat.SegmentManager`
        settings : :class:`~pylinac.vmat.Settings`

        Examples
        --------
        Run the DRGS demo:
            >>> VMAT().run_demo_drgs()

        Run the DRMLC demo:
            >>> VMAT().run_demo_drmlc()

        A typical use case:
            >>> open_img = "C:/QA Folder/VMAT/open_field.dcm"
            >>> dmlc_img = "C:/QA Folder/VMAT/dmlc_field.dcm"
            >>> myvmat = VMAT((open_img, dmlc_img))
            >>> myvmat.analyze(test='drmlc', tolerance=1.5)
            >>> print(myvmat.return_results())
            >>> myvmat.plot_analyzed_image()
    """
    def __init__(self, images=None):
        self.settings = Settings('', 1.5)
        if images is not None:
            self.load_images(images)

    @classmethod
    def from_images_UI(cls):
        """Construct a VMAT instance and select the files via a UI dialog box.

        .. versionadded:: 0.6
        """
        obj = cls()
        obj.load_images_UI()
        return obj

    def load_images_UI(self):
        """Load images via a UI dialog box. The open field must have 'open' in the name."""
        fs = get_filenames_UI()
        self.load_images(fs)

    @value_accept(im_type=(OPEN, DMLC))
    def load_image_UI(self, im_type='open'):
        """Open a UI file browser to load dicom image.

        Parameters
        ----------
        im_type : {'open', 'mlcs'}
            Specifies which file/image type is being loaded in.
        """
        if _is_open_type(im_type):
            caption = "Select Open Field EPID Image..."
        elif _is_dmlc_type(im_type):
            caption = "Select MLC Field EPID Image..."

        fs = get_filepath_UI(self, caption=caption)
        if fs:  # if user didn't hit cancel
            self.load_image(fs, im_type=im_type)

    @value_accept(im_type=(OPEN, DMLC))
    def load_image(self, file_path, im_type=None):
        """Load the image directly by the file path.

        Parameters
        ----------
        file_path : str, file-like object
            The path to the DICOM image or I/O stream.
        im_type : {'open', 'mlcs', None}
            Specifies what image type is being loaded in. If None, will try to determine the type from the name.
            The name must have 'open' or 'dmlc' in the name.
        """
        img = Image.load(file_path)
        if im_type is not None:
            if _is_open_type(im_type):
                self.image_open = img
            elif _is_dmlc_type(im_type):
                self.image_dmlc = img
        else:
            # try to guess type by the name
            imtype = self._try_to_guess_image_type(osp.basename(file_path))
            if imtype is not None:
                self.load_image(file_path, imtype)
            else:
                raise ValueError("Image type was not given nor could it be determined from the path name. Please enter and image type.")

    @staticmethod
    def _try_to_guess_image_type(name):
        if 'open' in name.lower():
            return 'open'
        elif 'dmlc' in name.lower():
            return 'dmlc'

    def _try_to_guess_test(self, name):
        # try to determine test type
        if DRMLC in name.lower():
            self.settings.test_type = DRMLC
        elif DRGS in name.lower():
            self.settings.test_type = DRGS

    def load_images(self, images, names=None):
        """Load both images simultaneously, assuming a clear name convention:
        1) The open image must have 'open' in the filename.
        2) Optionally, the test type can be in the image names as well. If so, the test type is automatically determined.

        .. versionadded:: 0.6

        Parameters
        ----------
        images : str
            File paths to the images. Order does not matter. The open image must have ``open`` somewhere in the name.
        names : None, str
            Internal keyword. If the passed images are URLs, then the names must also be passed. If loading URLs, use .from_urls() instead.
        """
        if (names is not None) and (len(names) != 2):
            raise ValueError("Exactly 2 images must be passed")
        elif (names is None) and (len(images) != 2):
            raise ValueError("Exactly 2 images must be passed")

        for idx, image in enumerate(images):
            if names is not None:
                name = names[idx]
            else:
                name = image
            if OPEN in name.lower():
                im_type = OPEN
            elif DMLC in name.lower():
                im_type = DMLC
            else:
                raise ValueError("Could not determine which image was open and DMLC; image names must have 'open' or 'dmlc' in the names when using `load_images()`")
            self.load_image(image, im_type=im_type)
        if names is not None:
            self._try_to_guess_test(name)
        else:
            self._try_to_guess_test(image)

    @classmethod
    def from_demo_images(cls, type='drgs'):
        """Construct a VMAT instance using the demo images.

        .. versionadded:: 0.6
        """
        obj = cls()
        obj.load_demo_image(type=type)
        return obj

    @value_accept(type=(DRGS, DRMLC))
    def load_demo_image(self, type='drgs'):
        """Load the demo DICOM images from demo files folder.

        Parameters
        ----------
        type : {'drgs', 'drmlc'}
            Test type of images to load.
        """
        demo_folder = osp.join(osp.dirname(osp.abspath(__file__)), "demo_files", 'vmat')
        if _test_is_mlcs(type):
            im_open_path = osp.join(demo_folder, "DRMLC_open.dcm")
            im_dmlc_path = osp.join(demo_folder, 'DRMLC_dmlc.dcm')
        else:
            im_open_path = osp.join(demo_folder, "DRGS_open.dcm")
            im_dmlc_path = osp.join(demo_folder, 'DRGS_dmlc.dcm')

        self.load_images((im_open_path, im_dmlc_path))

    @classmethod
    def from_zip(cls, path):
        """Load VMAT images from a ZIP file that contains both images.

        .. versionadded:: 0.8
        """
        obj = cls()
        obj.load_zip(path)
        return obj

    def load_zip(self, zip_file):
        """Load VMAT images from a ZIP file that contains both images.

        .. versionadded:: 0.8
        """
        zfiles = load_zipfile(zip_file)
        images = [BytesIO(zfiles.read(name)) for name in zfiles.namelist()]
        self.load_images(images, names=zfiles.namelist())

    @classmethod
    def from_urls(cls, urls):
        """Load from URLs.

        .. versionadded:: 0.7.1
        """
        obj = cls()
        obj.load_urls(urls)
        return obj

    def load_urls(self, urls):
        """Load from URLs.

        .. versionadded:: 0.7.1
        """
        url_dict = {}
        for idx, url in enumerate(urls):
            response = get_url(url)
            url_dict[idx] = BytesIO(response.content)
        imgs = [item for item in url_dict.values()]
        self.load_images(imgs, names=urls)

    @type_accept(tolerance=(int, float))
    def run_demo_drgs(self, tolerance=1.5):
        """Run the VMAT demo for the Dose Rate & Gantry Speed test."""
        self.load_demo_image('drgs')
        self.analyze(test='drgs', tolerance=tolerance, x_offset=20)  # old images (rev1, not new rev2's), which are offset
        print(self.return_results())
        self.plot_analyzed_image()

    @type_accept(tolerance=(int, float))
    def run_demo_drmlc(self, tolerance=1.5):
        """Run the VMAT demo for the MLC leaf speed test."""
        self.load_demo_image(DRMLC)
        self.analyze(test=DRMLC, tolerance=tolerance)
        print(self.return_results())
        self.plot_analyzed_image()

    @type_accept(test=str, x_offset=int)
    @value_accept(test=(DRMLC, DRGS), tolerance=(0, 8))
    def analyze(self, test=None, tolerance=1.5, x_offset=0):
        """Analyze the open and DMLC field VMAT images, according to 1 of 2 possible tests.

        Parameters
        ----------
        test : {'drgs', 'mlcs', None}
            The test to perform, either Dose Rate Gantry Speed ('drgs') or MLC Speed ('mlcs').
            The default is None, which assumes a clear naming convention; see load_image() for info.
        tolerance : float, int, optional
            The tolerance of the sample deviations in percent. Default is 1.5.
            Must be between 0 and 8.
        x_offset : int

            .. versionadded:: 0.8

            If the VMAT segments aren't aligned to the CAX (older VMAT images), the segments need a shift.
            Positive moves the segments right, negative to the left.
        """
        if not self.open_img_is_loaded or not self.dmlc_img_is_loaded:
            raise AttributeError("Open or MLC Image not loaded yet. Use .load_image() or .load_image_UI()")

        self.settings.tolerance = tolerance/100
        if test is None and self.settings.test_type is '':
            raise ValueError("Test parameter must be given, either 'drgs' or 'drmlc'")
        elif test is not None:
            self.settings.test_type = test
        self.settings.x_offset = x_offset

        """Pre-Analysis"""
        self._check_img_inversion()

        """Analysis"""
        self.segments = SegmentManager(self.image_open, self.image_dmlc, self.settings)

    def _check_img_inversion(self):
        """Check that the images are correctly inverted."""
        for image in [self.image_open, self.image_dmlc]:
            image.check_inversion()

    @property
    def avg_abs_r_deviation(self):
        """Return the average of the absolute R_deviation values."""
        return self.segments.avg_abs_r_deviation

    @property
    def avg_r_deviation(self):
        """Return the average of the R_deviation values, including the sign."""
        return self.segments.avg_r_deviation

    @property
    def max_r_deviation(self):
        """Return the value of the maximum R_deviation segment."""
        return self.segments.max_r_deviation

    @property
    def open_img_is_loaded(self):
        """Status of open image."""
        return hasattr(self, 'image_open')

    @property
    def dmlc_img_is_loaded(self):
        """Status of DMLC image."""
        return hasattr(self, 'image_dmlc')

    @property
    def passed(self):
        """Returns whether all the segment ratios passed the given tolerance."""
        return self.segments.passed

    def plot_analyzed_image(self, show=True):
        """Plot the analyzed images. Shows the open and dmlc images with the segments drawn; also plots the median
        profiles of the two images for visual comparison.

        Parameters
        ----------
        show : bool
            Whether to actually show the image.
        """
        fig, axes = plt.subplots(ncols=3, sharex=True)
        subimages = (OPEN, DMLC, PROFILE)
        titles = ('Open Image', 'DMLC Image', 'Median Profiles')
        for subimage, axis, title in zip(subimages, axes, titles):
            self.plot_analyzed_subimage(subimage=subimage, ax=axis, show=False)
            axis.set_title(title)

        if show:
            plt.show()

    def save_analyzed_image(self, filename, **kwargs):
        """Save the analyzed images as a png file.

        Parameters
        ----------
        filename : str, file-object
            Where to save the file to.
        kwargs
            Passed to matplotlib.
        """
        self.plot_analyzed_image(show=False)
        plt.savefig(filename, **kwargs)

    @value_accept(subimage=(DMLC, OPEN, PROFILE))
    def plot_analyzed_subimage(self, subimage='dmlc', show=True, ax=None):
        """Plot an individual piece of the VMAT analysis.

        Parameters
        ----------
        subimage : str
            Specifies which image to plot.
        show : bool
            Whether to actually plot the image.
        ax : matplotlib Axes, None
            If None (default), creates a new figure to plot to, otherwise plots to the given axes.
        """
        if ax is None:
            fig, ax = plt.subplots()

        # plot DMLC or OPEN image
        if subimage in (DMLC, OPEN):
            if subimage == DMLC:
                img = self.image_dmlc
            elif subimage == OPEN:
                img = self.image_open
            ax.imshow(img, cmap=plt.cm.Greys)
            self.segments.draw(ax)
            plt.sca(ax)
            plt.axis('off')

        # plot profile
        elif subimage == PROFILE:
            dmlc_prof, open_prof = self.median_profiles
            ax.plot(dmlc_prof.values)
            ax.plot(open_prof.values)
            ax.autoscale(axis='x', tight=True)
            ax.grid()

        if show:
            plt.show()

    def save_analyzed_subimage(self, filename, subimage='dmlc', interactive=False, **kwargs):
        """Save a subplot to file.

        Parameters
        ----------
        filename : str, file-object
            Where to save the file to.
        subimage : str
            Which subplot to save.
        interactive : bool
            Only applicable for the profile subplot. If False, saves as a .png image, else saves as an interactive
            HTML file.
        kwargs
            Passed to matplotlib.
        """
        self.plot_analyzed_subimage(subimage, show=False)
        if interactive and (subimage == PROFILE):
            mpld3 = import_mpld3()
            mpld3.save_html(plt.gcf(), filename)
        else:
            plt.savefig(filename, **kwargs)
        if isinstance(filename, str):
            print("VMAT subimage figure saved to {}".format(osp.abspath(filename)))

    @property
    def median_profiles(self):
        """Return two median profiles from the open and dmlc image. For visual comparison."""
        # dmlc median profile
        dmlc_prof = SingleProfile(np.median(self.image_dmlc, axis=0))
        dmlc_prof.stretch()
        # open median profile
        open_prof = SingleProfile(np.median(self.image_open, axis=0))
        open_prof.stretch()

        # normalize the profiles to near the same value
        norm_val = np.percentile(dmlc_prof.values, 90)
        dmlc_prof.normalize(norm_val)
        norm_val = np.percentile(open_prof.values, 90)
        open_prof.normalize(norm_val)

        return dmlc_prof, open_prof

    def return_results(self):
        """A string of the summary of the analysis results.

        Returns
        -------
        str
            The results string showing the overall result and deviation statistics by segment.
        """
        mlcs_seg_strs = ('1.7cm/s', '2.0cm/s', '1.0cm/s', '0.5cm/s')  # See TB RA procedure for values; != to Ling et al.
        DRGS_seg_strs = ('105MU/min', '210MU/min', '314MU/min', '417MU/min', '524MU/min', '592MU/min', '600MU/min')

        if self.passed:
            passfail_str = 'PASS'
        else:
            passfail_str = 'FAIL'

        if _test_is_drgs(self.settings.test_type):
            string = ('Dose Rate & Gantry Speed \nTest Results (Tol. +/-%2.1f%%): %s\n' %
                      (self.settings.tolerance * 100, passfail_str))
        elif _test_is_mlcs(self.settings.test_type):
            string = ('Dose Rate & MLC Speed \nTest Results (Tol. +/-%2.1f%%): %s\n' %
                      (self.settings.tolerance * 100, passfail_str))

        string += ('Max Deviation: %4.3f%%\n'
                   'Absolute Mean Deviation: %4.3f%%' %
                   (self.max_r_deviation, self.avg_abs_r_deviation))

        return string


class SegmentManager:
    """A class to handle the VMAT segments and perform actions across all segments."""

    def __init__(self, open_img, dmlc_img, settings):
        self.image_open = open_img
        self.image_dmlc = dmlc_img
        self.settings = settings
        self.segments = []

        points = self._construct_segment_centers(settings.test_type)
        self._construct_segments(points)

    def __len__(self):
        return len(self.segments)

    def __getitem__(self, item):
        return self.segments[item]

    def draw(self, plot):
        """Draw the segments onto a plot.

        Parameters
        ----------
        plot : matplotlib.axes.Axes
            The plot to draw the objects on.
        """
        for segment in self.segments:
            color = segment.get_bg_color()
            segment.add_to_axes(plot.axes, edgecolor=color)

    def _construct_segment_centers(self, test):
        """Construct the center points of the segments, depending on the test.

        Parameters
        ----------
        test : {'drgs', 'mlcs'}
        """
        if _test_is_drgs(test):
            num = 7
            offsets = DRGS_SETTINGS['X-plane offsets (mm)']
        else:
            num = 4
            offsets = MLCS_SETTINGS['X-plane offsets (mm)']

        points = []
        for seg, offset in zip(range(num), offsets):
            y = self.image_dmlc.center.y + self.settings.y_offset
            x_offset = offset * self.image_dmlc.dpmm + self.settings.x_offset
            x = self.image_dmlc.center.x + x_offset
            points.append(Point(x, y))
        return points

    def _construct_segments(self, center_points):
        """Construct Segment instances at the center point locations."""
        for point in center_points:
            segment = Segment(point, self.image_open, self.image_dmlc, self.settings.tolerance)
            self.segments.append(segment)
        self._update_r_corrs()

    def _update_r_corrs(self):
        """After the Segment constructions, the R_corr must be set for each segment."""
        for segment in self.segments:
            segment.r_dev = ((segment.r_corr / self._avg_r_corr) * 100) - 100

    @property
    def r_devs(self):
        """Return the deviations of all segments as an array."""
        return np.array([segment.r_dev for segment in self.segments])

    @property
    def _avg_r_corr(self):
        """Return the average R_corr of the segments."""
        return np.array([segment.r_corr for segment in self.segments]).mean()

    @property
    def avg_abs_r_deviation(self):
        """Return the average of the absolute R_deviation values."""
        return np.abs(self.r_devs).mean()

    @property
    def avg_r_deviation(self):
        """Return the average of the R_deviation values, including the sign."""
        return self.r_devs.mean()

    @property
    def max_r_deviation(self):
        """Return the value of the maximum R_deviation segment."""
        max_idx = np.abs(self.r_devs).argmax()
        return self.r_devs[max_idx]

    @property
    def passed(self):
        """Return whether all the images passed."""
        return all(segment.passed for segment in self.segments)


class Segment(Rectangle):
    """A class for holding and analyzing segment data of VMAT tests.

    For VMAT tests, there are either 4 or 7 'segments', which represents a section of the image that received
    radiation under the same conditions.

    Attributes
    ----------
    r_dev : float
            The reading deviation (R_dev) from the average readings of all the segments. See RTD for equation info.
    r_corr : float
        The corrected reading (R_corr) of the pixel values. See RTD for explanation and equation info.
    passed : boolean
        Specifies where the segment reading deviation was under tolerance.
    """
    # width of the segment (i.e. parallel to MLC motion) in pixels under reference conditions
    _nominal_width_mm = 5
    _nominal_height_mm = 100

    def __init__(self, center_point, open_image, dmlc_image, tolerance):
        self.r_dev = None  # is assigned after all segments constructed
        self._tolerance = tolerance
        self._open_image = open_image
        self._dmlc_image = dmlc_image
        width = self._nominal_width_mm * dmlc_image.dpmm
        height = self._nominal_height_mm * dmlc_image.dpmm
        super().__init__(width, height, center=center_point, as_int=True)

    @property
    def r_corr(self):
        """Return the ratio of the mean pixel values of DMLC/OPEN images."""
        dmlc_value = self._dmlc_image.array[self.bl_corner.y:self.bl_corner.y + self.height,
                     self.bl_corner.x: self.bl_corner.x + self.width].mean()
        open_value = self._open_image.array[self.bl_corner.y:self.bl_corner.y + self.height,
                     self.bl_corner.x: self.bl_corner.x + self.width].mean()
        ratio = (dmlc_value / open_value) * 100
        return ratio

    @property
    def passed(self):
        """Return whether the segment passed or failed."""
        return self.r_dev < self._tolerance * 100

    def get_bg_color(self):
        """Get the background color of the segment when plotted, based on the pass/fail status."""
        color = 'blue' if self.passed else 'red'
        return color


class Settings:
    """A helper class to hold analysis settings.

    Attributes
    ----------
    test_type : str
        The test type being done; either `drgs` or `mlcs`
    tolerance : float
        The tolerance value in percentage; e.g. 0.03 == 3%
    x_offset : int
        Offset in pixels to apply to all segments.
    y_offset : int
        Offset in pixels to apply to all segments.
    """
    test_type = typed_property('test_type', str)
    tolerance = typed_property('tolerance', (int, float))
    x_offset = typed_property('x_offset', int)
    y_offset = typed_property('y_offset', int)

    def __init__(self, test_type, tolerance):
        self.test_type = test_type
        self.tolerance = tolerance
        self.x_offset = 0
        self.y_offset = 0


def _is_open_type(image_type):
    return _x_in_y(image_type, OPEN)

def _is_dmlc_type(image_type):
    return _x_in_y(image_type, DMLC)

def _test_is_drgs(test):
    return _x_in_y(test, DRGS)

def _test_is_mlcs(test):
    return _x_in_y(test, DRMLC)

def _x_in_y(x, y):
    if x.lower() in y:
        return True
    else:
        return False


# -------------------
# VMAT demo
# -------------------
if __name__ == '__main__':
    VMAT().run_demo_drmlc()
