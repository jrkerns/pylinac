# -*- coding: utf-8 -*-
"""The VMAT module consists of the class VMAT, which is capable of loading an EPID DICOM Open field image and MLC field image and analyzing the
images according to the `Jorgensen et al. <http://dx.doi.org/10.1118/1.3552922>`_ tests, specifically the Dose-Rate GantryAxis-Speed (DRGS) and Dose-Rate MLC (DRMLC) tests.
"""
import os.path as osp

import numpy as np
import matplotlib.pyplot as plt

from pylinac.core.decorators import value_accept, type_accept
from pylinac.core.image import Image
from pylinac.core.geometry import Point, Rectangle
from pylinac.core.io import get_filepath_UI


test_types = {'DRGS': 'drgs', 'MLCS': 'mlcs'}
im_types = {'OPEN': 'open', 'DMLC': 'dmlc'}

DRGS_SETTINGS = {'X-plane offsets (cm)': (-6, -4, -2, 0, 2, 4, 6)}
MLCS_SETTINGS = {'X-plane offsets (cm)': (-4.5, -1.5, 1.5, 4.5)}

class VMAT:
    """The VMAT class analyzes two DICOM images acquired via a linac's EPID and analyzes
        regions of interest (segments) based on the paper by `Jorgensen et al <http://dx.doi.org/10.1118/1.3552922>`_,
        specifically, the Dose Rate & Gantry Speed (DRGS) and Dose Rate & MLC speed (DRMLC) tests.

        Examples
        --------
        Run the DRGS demo:
            >>> VMAT().run_demo_drgs()

        Run the DRMLC demo:
            >>> VMAT().run_demo_drmlc()

        A typical use case:
            >>> open_img = "C:/QA Folder/VMAT/open_field.dcm"
            >>> dmlc_img = "C:/QA Folder/VMAT/dmlc_field.dcm"
            >>> myvmat = VMAT()
            >>> myvmat.load_image(open_img, im_type='open')
            >>> myvmat.load_image(dmlc_img, im_type='dmlc')
            >>> myvmat.analyze(test='mlcs', tolerance=1.5)
            >>> print(myvmat.return_results())
            >>> myvmat.plot_analyzed_image()

        Attributes
        ----------
        image_open : :class:`ImageObj`
            The open-field image object.
        image_dmlc : core.image.ImageObj
            The dmlc-field image object.
        segments : list
            A list containing :class:`Segment` instances.
        settings : :class:`Settings`
            Settings for analysis.
    """
    def __init__(self):
        self.settings = Settings('', 1.5)
        self.segments = []

    @value_accept(im_type=im_types)
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

    @value_accept(im_type=im_types)
    def load_image(self, file_path, im_type='open'):
        """Load the image directly by the file path.

        Parameters
        ----------
        file_path : str
            The path to the DICOM image.
        im_type : {'open', 'mlcs'}
            Specifies which file/image type is being loaded in.
        """
        img = Image(file_path)
        if _is_open_type(im_type):
            self.image_open = img
        elif _is_dmlc_type(im_type):
            self.image_dmlc = img

    @value_accept(type=test_types)
    def load_demo_image(self, type='drgs'):
        """Load the demo DICOM images from demo files folder.

        Parameters
        ----------
        type : {'drgs', 'mlcs'}
            Test type of images to load.
        """
        demo_folder = osp.join(osp.dirname(osp.abspath(__file__)), "demo_files", 'vmat')
        if _test_is_mlcs(type):
            im_open_path = osp.join(demo_folder, "DRMLC_open.dcm")
            im_dmlc_path = osp.join(demo_folder, 'DRMLC_dmlc.dcm')
        else:
            im_open_path = osp.join(demo_folder, "DRGS_open.dcm")
            im_dmlc_path = osp.join(demo_folder, 'DRGS_dmlc.dcm')

        self.load_image(im_open_path, im_type=im_types['OPEN'])
        self.load_image(im_dmlc_path, im_type=im_types['DMLC'])

    @type_accept(tolerance=(int, float))
    def run_demo_drgs(self, tolerance=1.5):
        """Run the VMAT demo for the Dose Rate & Gantry Speed test."""
        self.load_demo_image('drgs')
        self.analyze(test='drgs', tolerance=tolerance)  # set tolerance to 2 to show some failures
        print(self.return_results())
        self.plot_analyzed_image()

    @type_accept(tolerance=(int, float))
    def run_demo_mlcs(self, tolerance=1.5):
        """Run the VMAT demo for the Dose Rate & MLC speed test."""
        self.load_demo_image('mlcs')
        self.analyze(test='mlcs', tolerance=tolerance)
        print(self.return_results())
        self.plot_analyzed_image()

    @type_accept(test=str)
    @value_accept(test=test_types, tolerance=(0.1, 8))
    def analyze(self, test, tolerance=1.5):
        """Analyze the open and DMLC field VMAT images, according to 1 of 2 possible tests.

        Parameters
        ----------
        test : {'drgs', 'mlcs'}
            The test to perform, either Dose Rate Gantry Speed ('drgs') or MLC Speed ('mlcs').
        tolerance : float, int, optional
            The tolerance of the sample deviations in percent. Default is 1.5.
            Must be between 0.1 and 8.
        """
        if not self.open_img_is_loaded or not self.dmlc_img_is_loaded:
            raise AttributeError("Open or MLC Image not loaded yet. Use .load_image() or .load_image_UI()")

        self.settings.test_type = test
        self.settings.tolerance = tolerance/100

        """Pre-Analysis"""
        self._check_img_inversion()

        """Analysis"""
        center_points = self._construct_segment_centers(test)
        self._construct_segments(center_points)

    def _construct_segment_centers(self, test):
        """Construct the center points of the segments, depending on the test.

        Parameters
        ----------
        test : {'drgs', 'mlcs'}
        """
        if _test_is_drgs(test):
            num = 7
            offsets = DRGS_SETTINGS['X-plane offsets (cm)']
        else:
            num = 4
            offsets = MLCS_SETTINGS['X-plane offsets (cm)']

        points = []
        for seg, offset in zip(range(num), offsets):
            y = self.image_dmlc.center.y + self.settings.y_offset
            x_offset = offset * 10 * self.image_dmlc.dpmm + self.settings.x_offset
            x = self.image_dmlc.center.x + x_offset
            points.append(Point(x, y))
        return points

    def _construct_segments(self, center_points):
        """Construct Segment instances at the center point locations."""
        for point in center_points:
            segment = Segment(point, self.image_open, self.image_dmlc, self.settings.tolerance)
            self.segments.append(segment)

        # go back and assign the deviation of each segment based on the average
        avg_r_corr = self._avg_r_corr
        for segment in self.segments:
            segment.r_dev = ((segment.r_corr / avg_r_corr) * 100) - 100

    def _check_img_inversion(self):
        """Check that the images are correctly inverted.

        Pixel value should increase with dose. This is ensured by
        sampling a corner and comparing to the mean image value.
        """
        top_corner = self.image_open.array[:20,:20].mean()
        img_mean = self.image_open.array.mean()
        if top_corner > img_mean:
            self.image_open.invert()
            self.image_dmlc.invert()

    @property
    def _avg_r_corr(self):
        """Return the average R_corr of the segments."""
        return np.array([segment.r_corr for segment in self.segments]).mean()

    @property
    def avg_abs_r_deviation(self):
        """Return the average of the absolute R_deviation values."""
        return np.abs(self._get_r_dev_array).mean()

    @property
    def avg_r_deviation(self):
        """Return the average of the R_deviation values, including the sign."""
        return self._get_r_dev_array.mean()

    @property
    def max_r_deviation(self):
        """Return the value of the maximum R_deviation segment."""
        max_idx = np.abs(self._get_r_dev_array).argmax()
        return self._get_r_dev_array[max_idx]

    @property
    def min_r_deviation(self):
        """Return the value of the minimum R_deviation segment."""
        min_idx = np.abs(self._get_r_dev_array).argmin()
        return self._get_r_dev_array[min_idx]

    @property
    def _get_r_dev_array(self):
        """Helper method to get an array of the R_deviation values of all segments."""
        return np.array([segment.r_dev for segment in self.segments])

    @property
    def open_img_is_loaded(self):
        """Status of open image."""
        if hasattr(self, 'image_open'):
            return True
        else:
            return False

    @property
    def dmlc_img_is_loaded(self):
        """Status of DMLC image."""
        if hasattr(self, 'image_dmlc'):
            return True
        else:
            return False

    @property
    def passed(self):
        """Returns whether all the segment ratios passed the given tolerance."""
        for segment in self.segments:
            if not segment.passed:
                return False
        return True

    def _draw_objects(self, plot):
        """Draw lines onto the matplotlib widget/figure showing the ROIS.

        Parameters
        ----------
        plot : matplotlib.axes.Axes
            The plot to draw the objects on.
        """
        for segment_num, segment in enumerate(self.segments):
            if segment.passed:
                color = 'blue'
            else:
                color = 'red'

            segment.add_to_axes(plot.axes, edgecolor=color)

    def plot_analyzed_image(self, plot1=None, plot2=None):
        """Create 1 figure with 2 plots showing the open and MLC images
            with the samples and results drawn on.

        Parameters
        ----------
        plot1 : matplotlib.axes.plot
            If passed, plots to this plot. If None, will create a new figure.
        plot2 : matplotlib.axes.plot
            Same as above; if plot1 is supplied but plot2 left as None, will put images into
            one figure.
        """
        if plot1 is None and plot2 is None:
            fig, (ax1, ax2) = plt.subplots(1,2)
        if plot1 is not None:
            ax1 = plot1.axes
        if plot2 is not None:
            ax2 = plot2.axes

        try:
            # http://stackoverflow.com/questions/3823752/display-image-as-grayscale-using-matplotlib
            ax2.imshow(self.image_open.array, cmap=plt.cm.Greys_r)
            ax2.set_title("Open Field Image")
        except:
            pass  # axis2 wasn't defined
        ax1.imshow(self.image_dmlc.array, cmap=plt.cm.Greys_r)
        ax1.set_title("MLC Field Image")
        self._draw_objects(ax1)

        # Finally, show it all
        if plot1 is None and plot2 is None:
            plt.show()
        if plot1 is not None:
            plot1.draw()
            plot1.axes.hold(False)
        if plot2 is not None:
            plot2.draw()
            plot2.axes.hold(False)

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

        string += ('\nOverall Results:\n'
                   'Max Positive Deviation: %4.3f%%\n'
                   'Max Negative Deviation: %4.3f%%\n'
                   'Absolute Mean Deviation: %4.3f%%' %
                   (self.max_r_deviation, self.min_r_deviation, self.avg_abs_r_deviation))

        return string


class Segment(Rectangle):
    """A class for holding and analyzing segment data of VMAT tests.

    For VMAT tests, there are either 4 or 7 'segments', which represents a section of the image that received
    radiation under the same conditions.

    Attributes
    ----------
    r_dev : float
            The reading deviation (R_dev) from the average readings of all the segments.
    r_corr : float
        The corrected reading of the pixel values defined as :math:`R_{corr}(x) = \frac{R_{DRGS}(x)}{R_{open}(x)} * 100`.
    passed : boolean
        Specifies where the segment reading deviation was under tolerance.
    """
    # width of the segment (i.e. parallel to MLC motion) in pixels under reference conditions
    _nominal_width_mm = 5
    _nominal_height_mm = 100

    def __init__(self, center_point, open_image, dmlc_image, tolerance):
        """
        """
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
        if self.r_dev > self._tolerance*100:
            return False
        else:
            return True


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
    def __init__(self, test_type, tolerance):
        self.test_type = test_type
        self.tolerance = tolerance
        self.x_offset = 0
        self.y_offset = 0


def _is_open_type(image_type):
    if image_type.lower() == 'open':
        return True
    else:
        return False

def _is_dmlc_type(image_type):
    if image_type.lower() == 'dmlc':
        return True
    else:
        return False

def _test_is_drgs(test):
    if test.lower() in 'drgs':
        return True
    else:
        return False

def _test_is_mlcs(test):
    if test.lower() in 'mlcs':
        return True
    else:
        return False

# -------------------
# VMAT demo
# -------------------
if __name__ == '__main__':
    vmat = VMAT()
    vmat.settings.x_offset = 0
    vmat.run_demo_mlcs()
    # VMAT().run_demo_drmlc()  # uncomment to run MLCS demo
