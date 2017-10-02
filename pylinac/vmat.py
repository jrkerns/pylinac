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
import os.path as osp
import io
from itertools import zip_longest

import matplotlib.pyplot as plt
import numpy as np

from .core import image
from .core.decorators import value_accept, type_accept
from .core.geometry import Point, Rectangle
from .core.io import get_url, TemporaryZipDirectory, retrieve_demo_file
from .core import pdf
from .core.profile import SingleProfile
from .core.utilities import typed_property, import_mpld3
from .settings import get_dicom_cmap

# test types
DRGS = 'drgs'
DRMLC = 'drmlc'
# delivery types
DMLC = 'dmlc'
OPEN = 'open'
# other constants
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
            >>> VMAT.run_demo_drgs()

        Run the DRMLC demo:
            >>> VMAT.run_demo_drmlc()

        Load images directly:
            >>> VMAT(images=('myvmat1.dcm', 'myvmat2.dcm'), delivery_types=['open', 'dmlc'])

        Or load from a ZIP archive:
            >>> VMAT.from_zip('myvmat.zip')

        Or load from a URL:
            >>> VMAT.from_url('http://myserver.com/vmat.zip')

        A typical use case:
            >>> open_img = "C:/QA Folder/VMAT/open_field.dcm"
            >>> dmlc_img = "C:/QA Folder/VMAT/dmlc_field.dcm"
            >>> myvmat = VMAT(images=(open_img, dmlc_img))
            >>> myvmat.analyze(test='drmlc', tolerance=1.5)
            >>> print(myvmat.return_results())
            >>> myvmat.plot_analyzed_image()
    """
    def __init__(self, images, delivery_types=None):
        """
        Parameters
        ----------
        images : iterable (list, tuple, etc)
            A sequence of paths to the image files.
        delivery_types : iterable, None
            If None (default), the image paths names must follow the `Naming Convention <http://pylinac.readthedocs.org/en/latest/vmat_docs.html#naming-convention>`_.
            If the image paths do not follow the naming convention, a 2-element string sequence for ``delivery_types`` must be passed in. E.g. ``['open', 'dmlc']``.
        """
        self.settings = Settings('', 1.5)

        # error checks
        if len(images) != 2:
            raise ValueError("Exactly 2 images (open, DMLC) must be passed")
        if delivery_types and len(delivery_types) != 2:
            raise ValueError("Delivery types must be 2 elements long")
        if delivery_types is None:
            delivery_types = []

        # walk over images and load the open and DMLC fields
        for img, deliv in zip_longest(images, delivery_types, fillvalue=''):
            if OPEN in img.lower() or OPEN == deliv.lower():
                self.image_open = image.load(img)
            elif DMLC in img.lower() or DMLC == deliv.lower():
                self.image_dmlc = image.load(img)
            else:
                raise ValueError("Image file names must follow the naming convention (e.g. 'DRGS_open.dcm'), or the delivery types must be passed explicitly")

        # try to determine test type
        if all(DRGS in osp.basename(img).lower() for img in images):
            self.settings.test_type = DRGS
        elif all(DRMLC in osp.basename(img).lower() for img in images):
            self.settings.test_type = DRMLC

    @classmethod
    @value_accept(type=(DRGS, DRMLC))
    def from_demo_images(cls, type=DRGS):
        """Construct a VMAT instance using the demo images.

        Parameters
        ----------
        type : {'drgs', 'drmlc'}
            Test type of images to load.
        """
        if type.lower() == DRMLC:
            url_type = 'drmlc.zip'
        else:
            url_type = 'drgs.zip'
        demo_file = retrieve_demo_file(url=url_type)
        return cls.from_zip(demo_file)

    @classmethod
    def from_zip(cls, path):
        """Load VMAT images from a ZIP file that contains both images. Must follow the naming convention.

        Parameters
        ----------
        path : str
            Path to the ZIP archive which holds the VMAT image files.
        """
        with TemporaryZipDirectory(path) as tmpzip:
            image_files = image.retrieve_image_files(tmpzip)
            return cls(images=image_files)

    @classmethod
    def from_url(cls, url):
        """Load a ZIP archive from a URL.  Must follow the naming convention.

        Parameters
        ----------
        url : str
            Must point to a valid URL that is a ZIP archive of two VMAT images.
        """
        zfile = get_url(url)
        return cls.from_zip(zfile)

    @staticmethod
    def run_demo_drgs(tolerance=1.5):
        """Run the VMAT demo for the Dose Rate & Gantry Speed test."""
        vmat = VMAT.from_demo_images(DRGS)
        vmat.analyze(test=DRGS, tolerance=tolerance, x_offset=20)  # old images (rev1, not new rev2's), which are offset
        print(vmat.return_results())
        vmat.plot_analyzed_image()

    @staticmethod
    def run_demo_drmlc(tolerance=1.5):
        """Run the VMAT demo for the MLC leaf speed test."""
        vmat = VMAT.from_demo_images(DRMLC)
        vmat.analyze(test=DRMLC, tolerance=tolerance)
        print(vmat.return_results())
        vmat.plot_analyzed_image()

    @type_accept(test=str, x_offset=int)
    @value_accept(test=(DRMLC, DRGS), tolerance=(0, 8))
    def analyze(self, test=None, tolerance=1.5, x_offset=0):
        """Analyze the open and DMLC field VMAT images, according to 1 of 2 possible tests.

        Parameters
        ----------
        test : {'drgs', 'mlcs', None}
            The test to perform, either Dose Rate Gantry Speed ('drgs') or MLC Speed ('mlcs').
            If None (default), the images passed in must have followed the clear `Naming Convention <http://pylinac.readthedocs.org/en/latest/vmat_docs.html#naming-convention>`_.
            If the images did not follow the naming convention then a string specifying the test type must be passed in.
        tolerance : float, int, optional
            The tolerance of the sample deviations in percent. Default is 1.5.
            Must be between 0 and 8.
        x_offset : int

            .. versionadded:: 0.8

            If the VMAT segments aren't aligned to the CAX (older VMAT images), the segments need a shift.
            Positive moves the segments right, negative to the left.
        """
        self.settings.tolerance = tolerance/100
        if test is None and self.settings.test_type is '':
            raise ValueError("Test parameter must be given, either 'drgs' or 'drmlc'")
        elif test is not None:
            self.settings.test_type = test.lower()
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
    def plot_analyzed_subimage(self, subimage=DMLC, show=True, ax=None):
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
            ax.imshow(img, cmap=get_dicom_cmap())
            self.segments.draw(ax)
            plt.sca(ax)
            plt.axis('off')
            plt.tight_layout()

        # plot profile
        elif subimage == PROFILE:
            dmlc_prof, open_prof = self.median_profiles
            ax.plot(dmlc_prof.values)
            ax.plot(open_prof.values)
            ax.autoscale(axis='x', tight=True)
            ax.grid()

        if show:
            plt.show()

    @value_accept(subimage=(DMLC, OPEN, PROFILE))
    def save_analyzed_subimage(self, filename, subimage=DMLC, interactive=False, **kwargs):
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
            print("VMAT subimage figure saved to {0}".format(osp.abspath(filename)))

    @property
    def median_profiles(self):
        """Return two median profiles from the open and dmlc image. For visual comparison."""
        # dmlc median profile
        dmlc_prof = SingleProfile(np.mean(self.image_dmlc, axis=0))
        dmlc_prof.stretch()
        # open median profile
        open_prof = SingleProfile(np.mean(self.image_open, axis=0))
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
        if self.passed:
            passfail_str = 'PASS'
        else:
            passfail_str = 'FAIL'

        if self.settings.test_type == DRGS:
            string = ('Dose Rate & Gantry Speed \nTest Results (Tol. +/-{0:2.1f}%): {1!s}\n'.format(self.settings.tolerance * 100, passfail_str))
        elif self.settings.test_type == DRMLC:
            string = ('Dose Rate & MLC Speed \nTest Results (Tol. +/-{0:2.1f}%): {1!s}\n'.format(self.settings.tolerance * 100, passfail_str))

        string += 'Max Deviation: {:2.3f}%\nAbsolute Mean Deviation: {:2.3f}%'.format(self.max_r_deviation, self.avg_abs_r_deviation)

        return string

    def publish_pdf(self, filename, author=None, unit=None, notes=None, open_file=False):
        """Publish (print) a PDF containing the analysis and quantitative results.

        Parameters
        ----------
        filename : (str, file-like object}
            The file to write the results to.
        """
        from reportlab.lib.units import cm
        canvas = pdf.create_pylinac_page_template(filename, file_name=osp.basename(self.image_open.path) + ", " + osp.basename(self.image_dmlc.path),
                                                  analysis_title='{} VMAT Analysis'.format(self.settings.test_type.upper()), author=author, unit=unit)
        for y, x, width, img in zip((10, 10, 1), (1, 11, 4), (9, 9, 16), (OPEN, DMLC, PROFILE)):
            data = io.BytesIO()
            self.save_analyzed_subimage(data, subimage=img)
            img = pdf.create_stream_image(data)
            canvas.drawImage(img, x * cm, y * cm, width=width * cm, height=18 * cm, preserveAspectRatio=True)
        text = ['{} VMAT results:'.format(self.settings.test_type.upper()),
                'Source-to-Image Distance (mm): {:2.0f}'.format(self.image_open.sid),
                'Tolerance (mm): {:2.1f}'.format(self.settings.tolerance),
                'X-offset applied (pixels): {:2.0f}'.format(self.settings.x_offset),
                'Absolute mean deviation (%): {:2.2f}'.format(self.avg_abs_r_deviation),
                'Maximum deviation (%): {:2.2f}'.format(self.max_r_deviation),
                ]
        pdf.draw_text(canvas, x=10 * cm, y=25.5 * cm, text=text)
        if notes is not None:
            pdf.draw_text(canvas, x=1 * cm, y=5.5 * cm, fontsize=14, text="Notes:")
            pdf.draw_text(canvas, x=1 * cm, y=5 * cm, text=notes)
        pdf.finish(canvas, open_file=open_file, filename=filename)


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
            segment.plot2axes(plot.axes, edgecolor=color)

    def _construct_segment_centers(self, test):
        """Construct the center points of the segments, depending on the test.

        Parameters
        ----------
        test : {'drgs', 'mlcs'}
        """
        if test == DRGS:
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
        return 'blue' if self.passed else 'red'


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
