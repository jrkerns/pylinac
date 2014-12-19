# -*- coding: utf-8 -*-
"""The VMAT module consists of the class VMAT, which is capable of loading an EPID DICOM Open field image and MLC field image and analyzing the
images according to the `Jorgensen et al. <http://dx.doi.org/10.1118/1.3552922>`_ tests, specifically the Dose-Rate GantryAxis-Speed (DRGS) and Dose-Rate MLC (DRMLC) tests.

"""
# Full documentation of this module and how to use it can be found at this repositories' Read the Docs site: pylinac.rtfd.org

from __future__ import print_function, division, absolute_import
from builtins import range

from future import standard_library

standard_library.install_aliases()
import os.path as osp

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Rectangle

from pylinac.core.decorators import value_accept, type_accept, lazyproperty
from pylinac.core.image_classes import ImageObj, AnalysisModule
from pylinac.core.geometry import Point, Box


# Common constants
DRMLC_seg_strs = ('1.6cm/s', '2.4cm/s', '0.8cm/s', '0.4cm/s')
DRGS_seg_strs = ('105MU/min', '210MU/min', '314MU/min', '417MU/min', '524MU/min', '592MU/min', '600MU/min')
test_types = {'DRGS': 'drgs', 'DRMLC': 'drmlc'}
im_types = {'OPEN': 'open', 'DMLC': 'dmlc'}

class VMAT(AnalysisModule):
    """The VMAT class analyzes two DICOM images acquired via a linac's EPID and analyzes
        regions of interest (segments) based on the paper by `Jorgensen et al <http://dx.doi.org/10.1118/1.3552922>`_, specifically,
        the Dose Rate & GantryAxis Speed (DRGS) and Dose Rate & MLC speed (DRMLC) tests.
    """

    def __init__(self):
        super().__init__()
        self.image_open = ImageObj()  # the Open field image
        self.image_dmlc = ImageObj()  # the MLC field image
        self._test_type = ''  # the test to perform
        self.tolerance = 3 # default of 3% tolerance as Jorgensen recommends
        # the following are dicts holding the individual segment max, min and mean deviation
        self.segments = []  # this is a list which will hold Segment objects (either 4 or 7)

    @value_accept(im_type=im_types)
    def load_image_UI(self, im_type='open'):
        """Open a UI file browser to load dicom image for given VMAT image.

        :param im_type: 'open' or 'mlc'; specifies which file/image is being loaded
        :type im_type: str
        """
        if  im_type == im_types['OPEN']:  # open
            caption = "Select Open Field EPID Image..."
            img = self.image_open
        else:  # dmlc
            caption = "Select MLC Field EPID Image..."
            img = self.image_dmlc

        fs = img._get_imagepath_UI(self, caption=caption)
        if fs:  # if user didn't hit cancel
            self.load_image(fs, im_type=im_type)

    @value_accept(im_type=im_types)
    def load_image(self, file_path, im_type='open'):
        """Load the image directly by the file path (i.e. non-interactively).

        :param file_path: The absolute path to the DICOM image
        :type file_path: str
        :param im_type: Specifies whether the image is the Open ('open') or DMLC ('dmlc') field.
        :type im_type: str
        """
        if im_type == im_types['OPEN']:
            image = self.image_open
        elif im_type == im_types['DMLC']:
            image = self.image_dmlc
        image.load_image(file_path)

    @value_accept(test_type=test_types)
    def load_demo_image(self, test_type='drgs'):
        """Load the demo DICOM images from demo files folder.

        :param test_type: Test type of images to load; 'drgs' or 'drmlc'
        :type test_type: str
        """
        demo_folder = osp.join(osp.split(osp.abspath(__file__))[0], "demo_files")
        if test_type == test_types['DRMLC']:  # DRMLC
            im_open_path = osp.join(demo_folder, "DRMLCopen-example.dcm")
            im_dmlc_path = osp.join(demo_folder, 'DRMLCmlc-example.dcm')
        else:
            im_open_path = osp.join(demo_folder, "DRGSopen-example.dcm")
            im_dmlc_path = osp.join(demo_folder, 'DRGSmlc-example.dcm')

        self.load_image(im_open_path, im_type=im_types['OPEN'])
        self.load_image(im_dmlc_path, im_type=im_types['DMLC'])

    def run_demo_drgs(self):
        """Run the demo of the module for the Dose Rate & GantryAxis Speed test.

        :param number: Which demo image set to run; there are currently 2 demos.
        :type number: int
        """
        self.load_demo_image('drgs')
        self.analyze(test='drgs', tolerance=3)  # tolerance at 2 to show some failures
        print(self.get_string_results())
        self.plot_analyzed_image()

    def run_demo_drmlc(self):
        """Run the demo of the module for the Dose Rate & MLC speed test.

        :param number: Which demo image set to run; there are currently 2 demos.
        :type number: int
        """
        self.load_demo_image('drmlc')
        self.analyze(test='drmlc', tolerance=3)
        print(self.get_string_results())
        self.plot_analyzed_image()

    def _calc_im_scaling_factors(self, SID=None):
        """Determine the image scaling factors; factors are relative to reference values from images of size 384x512 taken at 150cm SID."""
        # Image size scaling
        y_scale = self.image_dmlc.pixel_array.shape[0] / 384.0
        x_scale = self.image_dmlc.pixel_array.shape[1] / 512.0
        scale = Point(x_scale, y_scale)

        # SID scaling
        if SID:
            SID_scale = SID / 150
        else:
            try:
                SID_scale = self.image_open.properties['SID mm'] / 1500.0
            except:
                SID_scale = 1

        return SID_scale, scale

    def construct_segments(self, test, scale, SID_scale, HDMLC):
        """Construct the 4 or 7 Segments of the test."""

        # DRMLC segment x-direction offsets
        if test == test_types['DRMLC']:
            segment_offsets = (-86, -27, 32, 89)  # These are the x-direction offsets of the segments.
        # ditto for DRGS
        elif test == test_types['DRGS']:
            segment_offsets = (-97, -59, -21, 18, 55, 95, 134)

        # init and append the Segments to the segments list
        for seg_offset in segment_offsets:
            seg_center = Point(self.image_dmlc.center.x + seg_offset, self.image_dmlc.center.y)
            segment = Segment(test, seg_center, scale, SID_scale, HDMLC)
            self.segments.append(segment)

    def _extract_sample_ratios(self):
        """Extract the mean of the pixel values of the samples from both Open and DMLC images and ratio them."""
        for segment in self.segments:
            segment.extract_sample_ratios(self.image_open, self.image_dmlc)

    def _normalize_samples(self):
        """Normalize the samples for each image respectively to the mean of the central segment."""

        if self._test_type == test_types['DRGS']:
            normalization_segment = 3
        elif self._test_type == test_types['DRMLC']:
            normalization_segment = 2

        normalization_value = self.segments[normalization_segment].mean_ratio

        for segment in self.segments:
            segment.normalize_sample_ratios(normalization_value)

    def _calc_deviations(self):
        """Calculate the deviations of each sample using Jorgensen's deviation equation."""
        # preallocation
        deviation = np.zeros((self.num_samples, self.num_segments))
        # deviation calculation
        for sample_row in range(len(self.segments[0].samples)):
            # calculate mean of all 4 or 7 samples of a given sample row
            samples = np.array([segment.samples[sample_row].ratio for segment in self.segments])
            # the row mean is the denominator of deviation equation
            row_mean = samples.mean()
            for segment, sample in enumerate(samples):
                # each sample value in the row is the numerator in the deviation equation
                deviation[sample_row, segment] = (sample / row_mean) - 1

        deviation *= 100  # convert to percentage

        # set deviation values to segment property
        for segment_num, segment in enumerate(self.segments):
            segment.deviations = deviation[:, segment_num]

    @type_accept(test=str, tolerance=(float, int), SID=(None, int, float))
    @value_accept(test=test_types, tolerance=(0.3, 8), SID=(0, 180))
    def analyze(self, test, tolerance=3, SID=None, HDMLC=False):
        """Analyze 2 VMAT images, the open field image and DMLC field image, according to 1 of 2 possible tests.

        :param test: The test to perform, either Dose Rate GantryAxis Speed ('drgs') or Dose Rate MLC Speed ('drmlc').
        :type test: str
        :param tolerance: The tolerance of the sample deviations in percent. Default is 3, as Jorgensen recommends.
        :type tolerance: float, int
        :param SID: The Source to Image (detector) distance in cm. Usually doesn't need to be passed for EPID DICOM images. This argument
            will override any automatically derived value however. If left as None and no SID was determined, it will assume 150cm.
        :type SID: None, int
        """
        # error checking
        if self.image_open.pixel_array.size == 0  or self.image_dmlc.pixel_array.size == 0:
            raise AttributeError("Open or MLC Image not loaded yet. Use .load_image() or .load_image_UI()")

        self._check_img_inversion()

        self.tolerance = tolerance / 100.0
        self._test_type = test

        # get the image scaling factors and center pixels; this corrects for the SID
        SID_scale, scale = self._calc_im_scaling_factors(SID)

        # set up pixel bounds of test
        self.construct_segments(test, scale, SID_scale, HDMLC)

        # Extract the samples from each image
        self._extract_sample_ratios()

        # normalize the samples for each image respectively
        self._normalize_samples()

        # calculate deviations as per Jorgensen equation
        self._calc_deviations()

    def _check_img_inversion(self):
        """Check that the images are correctly inverted (pixel value increases with dose) by
        sampling a corner and comparing to the image mean."""
        top_corner = self.image_open.pixel_array[:20,:20].mean()
        img_mean = self.image_open.pixel_array.mean()
        if top_corner > img_mean:
            self.image_open.invert_array()
            self.image_dmlc.invert_array()

    @property
    def sample_pass_matrix(self):
        """Return a num_segments-by-num_samples boolean matrix indicating pass/fail of the sample ratio."""
        sample_passfail_matrix = np.zeros((self.num_samples, self.num_segments), dtype=bool)

        # for each sample, if ratio is within tolerance, set to True
        for seg_num, segment in enumerate(self.segments):
            for sam_num, sample in enumerate(segment.samples):
                if sample.ratio < 1 + self.tolerance and sample.ratio > 1 - self.tolerance:
                    sample_passfail_matrix[sam_num, seg_num] = True

        return sample_passfail_matrix

    @property
    def num_samples(self):
        return len(self.segments[0].samples)

    @property
    def num_segments(self):
        return len(self.segments)

    @property
    def overall_max_deviation(self):
        return np.max([segment.max_dev for segment in self.segments])

    @property
    def overall_min_deviation(self):
        return np.min([segment.min_dev for segment in self.segments])

    @property
    def overall_abs_mean_deviation(self):
        return np.mean([segment.abs_mean_dev for segment in self.segments])

    @property
    def passed(self):
        """Returns whether all the sample ratios passed the given tolerance."""
        if self.sample_pass_matrix.all():
            return True
        else:
            return False

    def _draw_objects(self, plot, draw=False):
        """
        This method draws lines on the matplotlib widget/figure showing an outline of the ROIs used in
        the calculation and text of the results.

        :type plot: matplotlib.axes.Axes
        """

        # clear images of any existing objects that are on it
        # if hasattr(self, 'roi_handles'):
        #     for handle in self.roi_handles:
        #         handle.remove()
        #     if draw:
        #         plot.draw()

        #plot ROI lines on image
        self.roi_handles = [[None for _ in range(self.num_segments)] for _ in range(self.num_samples)]

        pass_matrix = self.sample_pass_matrix

        for segment_num, segment in enumerate(self.segments):
            for sample_num, sample in enumerate(segment.samples):
                if pass_matrix[sample_num, segment_num]:
                    color = 'blue'
                else:
                    color = 'red'

                self.roi_handles[sample_num][segment_num] = Rectangle((sample.bl_corner.x,
                                                                       sample.bl_corner.y),
                                                                      height=sample.height, width=sample.width,
                                                                      fill=False, edgecolor=color, linewidth=1)
                plot.axes.add_patch(self.roi_handles[sample_num][segment_num])

        if draw:
            plot.draw()

    def plot_analyzed_image(self, plot1=None, plot2=None):
        """Create 1 figure of 2 plots showing the open and MLC images with the samples and results drawn on.

        :param plot1: If None, will create a new figure.
        :param plot2: Only for PyQA use.
        """
        if plot1 is None and plot2 is None:
            fig, (ax1, ax2) = plt.subplots(1,2)
        if plot1 is not None:
            ax1 = plot1.axes
        if plot2 is not None:
            ax2 = plot2.axes

        try:
            # http://stackoverflow.com/questions/3823752/display-image-as-grayscale-using-matplotlib
            ax2.imshow(self.image_open.pixel_array, cmap=cm.Greys_r)
            ax2.set_title("Open Field Image")
        except:
            pass  # axis2 wasn't defined
        ax1.imshow(self.image_dmlc.pixel_array, cmap=cm.Greys_r)
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

    def get_string_results(self):
        """Create a string of the summary of the analysis results.

        :returns: str -- The results string showing the overall result and deviation stats by segment.
        """

        if self.passed:
            passfail_str = 'PASS'
        else:
            # TODO: count how many failures
            passfail_str = 'FAIL'

        if self._test_type == test_types['DRGS']:
            string = ('Dose Rate & Gantry Speed \nTest Results (Tol. +/-%2.1f%%): %s\n' %
                      (self.tolerance * 100, passfail_str))
        elif self._test_type == test_types['DRMLC']:
            string = ('Dose Rate & MLC Speed \nTest Results (Tol. +/-%2.1f%%): %s\n' %
                      (self.tolerance * 100, passfail_str))

        string += ('\nOverall Results:\n'
                   'Max Positive Deviation: %4.3f%%\n'
                   'Max Negative Deviation: %4.3f%%\n'
                   'Absolute Mean Deviation: %4.3f%%' %
                   (self.overall_max_deviation, self.overall_min_deviation, self.overall_abs_mean_deviation))

        if self._test_type == test_types['DRGS']:
            seg_str = DRGS_seg_strs
        else:
            seg_str = DRMLC_seg_strs

        for seg_num, segment in enumerate(self.segments):
            string += ('\n\n%s Segment:\n'
                       'Max Positive Deviation: %4.3f%%\n'
                       'Max Negative Deviation: %4.3f%%\n'
                       'Absolute Mean Deviation: %4.3f%%' %
                       (seg_str[seg_num], segment.max_dev, segment.min_dev, segment.abs_mean_dev))

        return string


class Sample(Box):
    """A class representing a single 'sample' of a VMAT segment. A sample is an ROI of the radiation one MLC pair
        produces in one segment.
    """
    def __init__(self, width, height, center):
        super().__init__(width, height, center, as_int=True)
        self.ratio = 0

    def extract_ratio(self, open_image, dmlc_image):
        dmlc_value = dmlc_image.pixel_array[self.bl_corner.y - self.height:self.bl_corner.y,
                                self.bl_corner.x: self.bl_corner.x + self.width].mean()
        open_value = open_image.pixel_array[self.bl_corner.y - self.height:self.bl_corner.y,
                                self.bl_corner.x: self.bl_corner.x + self.width].mean()
        self.ratio = dmlc_value/open_value

    def normalize_ratio(self, norm_value):
        self.ratio /= norm_value


class Segment(object):
    """A class for holding and analyzing segment data of VMAT tests.

    For VMAT tests, there are either 4 or 7 'segments', which represents a section of the image that received
    radiation under the same conditions. Segments hold a number of Samples.
    """

    def __init__(self, test, center_point, scaling, SID_scale, HDMLC):
        self.center = center_point
        if test == test_types['DRGS']:  #DRGS
            width = np.round(24 * scaling.x * SID_scale)  # length of the sample parallel to MLC movement.
        else:
            width = np.round(42 * scaling.x * SID_scale)

        num_leaves = 38  # cut off last leaves on either end as per Jorgensen.
        leaf_spacing = 9.566 * scaling.y * SID_scale # the width from one leaf to the next leaf (perpendicular to the MLC movement).
        # 0.784mm/pixel@100cm / 1.5 = 0.52267mm/pixel@150cm. Solving for 5mm gives 9.566. This must be corrected for image scaling.
        sample_spacing = 6 * scaling.y * SID_scale # the width in pixels of a given MLC leaf sample perpendicular to MLC movement at 150cm SID.
        # This must be corrected for image scaling
        offset = leaf_spacing*num_leaves/2
        if HDMLC:
            num_leaves *= 2
            leaf_spacing /= 2
            sample_spacing /= 2


        self.samples = []  # Will be a list of Samples within segment
        for leaf in range(num_leaves):
            # the sample center has the same x-coordinate, but an offset y-coordinate.
            # The y-coord is calculated by offsetting by the total range/2 (offset) plus half a leaf (center point is in middle of leaf),
            # and finally by changing the offset by a leaf's width for each sample (leaf*leaf_spacing).
            sample_center = Point(center_point.x,
                                  (center_point.y + leaf*leaf_spacing + leaf_spacing/2 - offset))
            sample = Sample(width, height=sample_spacing, center=sample_center)
            self.samples.append(sample)

    def extract_sample_ratios(self, open_image, dmlc_image):
        """Extract the ratios of open to dmlc fields for all the sample ROIs."""
        for sample in self.samples:
            sample.extract_ratio(open_image, dmlc_image)

    def normalize_sample_ratios(self, norm_value):
        """Normalize the sample ratios by the normalization value passed."""
        for sample in self.samples:
            sample.normalize_ratio(norm_value)

    @property
    def deviations(self):
        return self._deviations

    @deviations.setter
    def deviations(self, value):
        if len(value) != len(self.samples):
            raise ValueError("Deviation matrix is not correctly sized")
        self._deviations = value

    @lazyproperty
    def abs_mean_dev(self):
        absolute_deviations = np.abs(self.deviations)
        return absolute_deviations.mean()

    @property
    def mean_ratio(self):
        return self._ratios.mean()

    @lazyproperty
    def min_dev(self):
        return self.deviations.min()

    @lazyproperty
    def max_dev(self):
        return self.deviations.max()

    @property
    def _ratios(self):
        return np.array([sample.ratio for sample in self.samples])



#---------------------------------------------------------------------------------------------------------------------
# VMAT demo.
#---------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # VMAT().run_demo_drgs()
    VMAT().run_demo_drmlc()  # uncomment to run MLCS demo
