# -*- coding: utf-8 -*-
"""The VMAT module consists of the class VMAT, which is capable of loading an EPID DICOM Open field image and MLC field image and analyzing the
images according to the `Jorgensen et al. <http://dx.doi.org/10.1118/1.3552922>`_ tests, specifically the Dose-Rate GantryAxis-Speed (DRGS) and Dose-Rate MLC (DRMLC) tests.

"""
# Full documentation of this module and how to use it can be found at this repositories' Read the Docs site: pylinac.rtfd.org

import os.path as osp

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from pylinac.core.decorators import value_accept, type_accept, lazyproperty
from pylinac.core.image import ImageObj
from pylinac.core.geometry import Point, Rectangle
from pylinac.core.io import get_filepath_UI

DRMLC_seg_strs = ('1.6cm/s', '2.4cm/s', '0.8cm/s', '0.4cm/s')
DRGS_seg_strs = ('105MU/min', '210MU/min', '314MU/min', '417MU/min', '524MU/min', '592MU/min', '600MU/min')
test_types = {'DRGS': 'drgs', 'DRMLC': 'drmlc'}
im_types = {'OPEN': 'open', 'DMLC': 'dmlc'}

class VMAT:
    """The VMAT class analyzes two DICOM images acquired via a linac's EPID and analyzes
        regions of interest (segments) based on the paper by `Jorgensen et al <http://dx.doi.org/10.1118/1.3552922>`_,
        specifically, the Dose Rate & Gantry Speed (DRGS) and Dose Rate & MLC speed (DRMLC) tests.
    """

    def __init__(self):
        super().__init__()
        self.image_open = ImageObj()  # the Open field image
        self.image_dmlc = ImageObj()  # the MLC field image
        self._test_type = ''  # the test to perform
        self.tolerance = 3 # default of 3% tolerance as Jorgensen recommends
        self.segments = []  # a list which will hold Segment objects (either 4 or 7)

    @value_accept(im_type=im_types)
    def load_image_UI(self, im_type='open'):
        """Open a UI file browser to load dicom image.

        Parameters
        ----------
        im_type : {'open', 'dmlc'}
            Specifies which file/image type is being loaded in.
        """
        if  im_type == im_types['OPEN']:  # open
            caption = "Select Open Field EPID Image..."
        else:  # dmlc
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
        im_type : {'open', 'dmlc'}
            Specifies which file/image type is being loaded in.
        """
        if im_type == im_types['OPEN']:
            image = self.image_open
        elif im_type == im_types['DMLC']:
            image = self.image_dmlc
        image.load_image(file_path)

    @value_accept(test_type=test_types)
    def load_demo_image(self, test_type='drgs'):
        """Load the demo DICOM images from demo files folder.

        Parameters
        ----------
        test_type : {'drgs', 'drmlc'}
            Test type of images to load.
        """
        demo_folder = osp.join(osp.dirname(osp.abspath(__file__)), "demo_files", 'vmat')
        if test_type == test_types['DRMLC']:  # DRMLC
            im_open_path = osp.join(demo_folder, "DRMLC_open.dcm")
            im_dmlc_path = osp.join(demo_folder, 'DRMLC_dmlc.dcm')
        else:
            im_open_path = osp.join(demo_folder, "DRGS_open.dcm")
            im_dmlc_path = osp.join(demo_folder, 'DRGS_dmlc.dcm')

        self.load_image(im_open_path, im_type=im_types['OPEN'])
        self.load_image(im_dmlc_path, im_type=im_types['DMLC'])

    def run_demo_drgs(self):
        """Run the VMAT demo for the Dose Rate & Gantry Speed test."""
        self.load_demo_image('drgs')
        self.analyze(test='drgs', tolerance=3)  # tolerance at 2 to show some failures
        print(self.get_string_results())
        self.plot_analyzed_image()

    def run_demo_drmlc(self):
        """Run the VMAT demo for the Dose Rate & MLC speed test."""
        self.load_demo_image('drmlc')
        self.analyze(test='drmlc', tolerance=3)
        print(self.get_string_results())
        self.plot_analyzed_image()

    def _calc_im_scaling_factors(self, SID=None):
        """Determine image scaling factors.

         Factors are relative to reference values from images of size 384x512 taken at 150cm SID.

         See Also
         --------
         analyze : Further parameter info.
         """
        # Image size scaling
        y_scale = self.image_dmlc.pixel_array.shape[0] / 384.0
        x_scale = self.image_dmlc.pixel_array.shape[1] / 512.0
        scale = Point(x_scale, y_scale)

        # SID scaling
        if self.image_open.SID:
            SID_scale = self.image_open.SID / 150.0
        elif SID is not None:
            SID_scale = SID / 150.0
        else:
            SID_scale = 1

        return SID_scale, scale

    def construct_segments(self, test, scale, SID_scale, HDMLC):
        """Construct the 4 or 7 Segments of the test, depending on the test.

        See Also
        --------
        Segment : Further parameter info.
        """

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
        for sample_row in range(self.num_samples):
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

    @type_accept(test=str)
    @value_accept(test=test_types, tolerance=(0.3, 8), SID=(0, 180))
    def analyze(self, test, tolerance=3, SID=None, HDMLC=False):
        """Analyze the open and DMLC field VMAT images, according to 1 of 2 possible tests.

        Parameters
        ----------
        test : {'drgs', 'drmlc'}
            The test to perform, either Dose Rate Gantry Speed ('drgs') or Dose Rate MLC Speed ('drmlc').
        tolerance : number
            The tolerance of the sample deviations in percent. Default is 3, as Jorgensen recommends.
            Must be between 0.3 and 8.
        SID : number
            The Source to Image (detector) distance in cm. Usually doesn't need to be passed for EPID DICOM images. This argument
            will override any automatically derived value however. If left as None and no SID was determined, it will assume 150cm.
        """
        # error checking
        if not self.open_img_is_loaded or not self.dmlc_img_is_loaded:
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
        """Check that the images are correctly inverted.

        Pixel value should increase with dose. This is ensured by
        sampling a corner and comparing to the mean image value.
        """
        top_corner = self.image_open.pixel_array[:20,:20].mean()
        img_mean = self.image_open.pixel_array.mean()
        if top_corner > img_mean:
            self.image_open.invert_array()
            self.image_dmlc.invert_array()

    @property
    def sample_pass_matrix(self):
        """Return a num_segments-by-num_samples boolean matrix indicating pass/fail status of all the samples."""
        sample_passfail_matrix = np.zeros((self.num_samples, self.num_segments), dtype=bool)

        # for each sample, if ratio is within tolerance, set to True
        # TODO: probably replacable with np.where()
        for seg_num, segment in enumerate(self.segments):
            for sam_num, sample in enumerate(segment.samples):
                if sample.ratio < 1 + self.tolerance and sample.ratio > 1 - self.tolerance:
                    sample_passfail_matrix[sam_num, seg_num] = True

        return sample_passfail_matrix

    @property
    def open_img_is_loaded(self):
        """Status of open image."""
        if self.image_open.pixel_array.size != 0:
            return True
        else:
            return False

    @property
    def dmlc_img_is_loaded(self):
        """Status of DMLC image."""
        if self.image_dmlc.pixel_array.size != 0:
            return True
        else:
            return False

    @property
    def num_samples(self):
        """The number of 'samples' the algorithm created."""
        return len(self.segments[0].samples)

    @property
    def num_segments(self):
        """The number of segments the algorithm created."""
        return len(self.segments)

    @property
    def overall_max_deviation(self):
        """The maximum deviation of any sample."""
        return np.max([segment.max_dev for segment in self.segments])

    @property
    def overall_min_deviation(self):
        """The minimum deviation of any sample."""
        return np.min([segment.min_dev for segment in self.segments])

    @property
    def overall_abs_mean_deviation(self):
        """The mean absolute deviation of all samples."""
        return np.mean([segment.abs_mean_dev for segment in self.segments])

    @property
    def passed(self):
        """Returns whether *all* the sample ratios passed the given tolerance."""
        if self.sample_pass_matrix.all():
            return True
        else:
            return False

    def _draw_objects(self, plot):
        """Draw lines onto the matplotlib widget/figure showing the ROIS.

        Parameters
        ----------
        plot : matplotlib.axes.Axes
            The plot to draw the objects tol
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

                self.roi_handles[sample_num][segment_num] = sample.add_to_axes(plot.axes, edgecolor=color)

    def plot_analyzed_image(self, plot1=None, plot2=None):
        """Create 1 figure with 2 plots showing the open and MLC images
            with the samples and results drawn on.

        Parameters
        ----------
        plot1 : matplotlib.axes.plot
            If passed, plots to this plot. If None, will create a new figure.
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
        """A string of the summary of the analysis results.

        Returns
        -------
        str
            The results string showing the overall result and deviation statistics by segment.
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


class Sample(Rectangle):
    """Represents a single 'sample' of a VMAT segment. A sample is an ROI of the radiation one MLC pair
        produces in one VMAT segment.
    """
    def __init__(self, width, height, center):
        """
        Parameters
        ----------
        width : int
            Width of the sample in pixels.
        height : int
            Height of the sample in pixels.
        center : geometry.Point
            Center Point of the sample.
        """
        super().__init__(width, height, center, as_int=True)
        self.ratio = 0

    def extract_ratio(self, open_image, dmlc_image):
        """Extract the ratio (DMLC/open) of the pixel values in the sample.

        Parameters
        ----------
        open_image : numpy.ndarray
            The open image array
        dmlc_image : numpy.ndarray
            The DMLC image array
        """
        dmlc_value = dmlc_image.pixel_array[self.bl_corner.y - self.height:self.bl_corner.y,
                                self.bl_corner.x: self.bl_corner.x + self.width].mean()
        open_value = open_image.pixel_array[self.bl_corner.y - self.height:self.bl_corner.y,
                                self.bl_corner.x: self.bl_corner.x + self.width].mean()
        self.ratio = dmlc_value/open_value

    def normalize_ratio(self, norm_value):
        self.ratio /= norm_value


class Segment:
    """A class for holding and analyzing segment data of VMAT tests.

    For VMAT tests, there are either 4 or 7 'segments', which represents a section of the image that received
    radiation under the same conditions. Segments hold a number of Samples.
    """
    DRGS_width = 24  # width of the segment (i.e. parallel to MLC motion) in pixels under reference conditions
    DRMLC_width = 42  # ditto for DRMLC

    def __init__(self, test, center_point, scaling, SID_scale, HDMLC):
        """
        Parameters
        ----------
        test : {'drgs', 'drmlc'}
            The test being performed, either Dose-Rate & Gantry Speed or Dose-Rate and MLC Speed.
        center_point : geometry.Point
            The center point of the Segment.
        scaling :
            The scaling of the image compared to reference values.
        SID_scale : float
            The ratio of the actual SID distance to reference SID distance (150cm).
        HDMLC : bool
            If False (default), assumes a standard Millenium 120 leaf MLC (5mm leaves).
            If True, assumes an HD Millenium MLC which has a leaf width value of 2.5mm.

        """
        self.center = center_point
        if test == test_types['DRGS']:  #DRGS
            width = np.round(self.DRGS_width * scaling.x * SID_scale)  # length of the sample parallel to MLC movement.
        else:
            width = np.round(self.DRMLC_width * scaling.x * SID_scale)

        num_leaves = 38  # cut off last leaves on either end as per Jorgensen.
        leaf_spacing = 9.566 * scaling.y * SID_scale # the width from one leaf to the next leaf (perpendicular to the MLC movement).
        # 0.784mm/pixel@100cm / 1.5 = 0.52267mm/pixel@150cm. Solving for 5mm gives 9.566. This must be corrected for image scaling.
        sample_spacing = 6 * scaling.y * SID_scale # the width in pixels of a given MLC leaf sample perpendicular to MLC movement at 150cm SID.
        # This must be corrected for image scaling
        if HDMLC:
            num_leaves *= 2
            leaf_spacing /= 2
            sample_spacing /= 2

        offset = leaf_spacing*num_leaves/2

        self.samples = []  # Will be a list of Samples within segment
        for leaf in range(num_leaves):
            # the sample center has the same x-coordinate, but an offset y-coordinate.
            # The y-coord is calculated by offsetting by the total range/2 (offset) plus half a leaf (center point is in middle of leaf),
            # and finally by changing the offset by a leaf's width for each sample (leaf*leaf_spacing).
            sample_center = Point(center_point.x,
                                  (center_point.y + leaf*leaf_spacing + leaf_spacing/2 - offset))
            sample = Sample(width=width, height=sample_spacing, center=sample_center)
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
        """Returns a matrix of len num_samples with the sample deviation values."""
        return self._deviations

    @deviations.setter
    def deviations(self, value):
        """Set the deviations matrix. Must be length of the num_samples"""
        if len(value) != len(self.samples):
            raise ValueError("Deviation matrix is not correctly sized")
        self._deviations = value

    @lazyproperty
    def abs_mean_dev(self):
        """Return the mean absolute deviation of the samples."""
        absolute_deviations = np.abs(self.deviations)
        return absolute_deviations.mean()

    @property
    def mean_ratio(self):
        """Return the mean ratio (DMLC/open) of the samples."""
        return self._ratios.mean()

    @lazyproperty
    def min_dev(self):
        """Return the min ratio of the samples."""
        return self.deviations.min()

    @lazyproperty
    def max_dev(self):
        """Return the max deviation of the samples."""
        return self.deviations.max()

    @property
    def _ratios(self):
        """Return an array of sample ratios the length of num_samples."""
        return np.array([sample.ratio for sample in self.samples])



#---------------------------------------------------------------------------------------------------------------------
# VMAT demo.
#---------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # VMAT().run_demo_drgs()
    VMAT().run_demo_drmlc()  # uncomment to run MLCS demo
