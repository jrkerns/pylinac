# -*- coding: utf-8 -*-
"""The VMAT module consists of the class VMAT, which is capable of loading an EPID DICOM Open field image and MLC field image and analyzing the
images according to the `Jorgensen et al. <http://dx.doi.org/10.1118/1.3552922>`_ tests, specifically the Dose-Rate Gantry-Speed (DRGS) and Dose-Rate MLC (DRMLC) tests.

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

from pylinac.common.decorators import value_accept
from pylinac.common.image_classes import SingleImageObject


# Common constants
DRMLC_seg_strs = ('1.6cm/s', '2.4cm/s', '0.8cm/s', '0.4cm/s')
DRGS_seg_strs = ('105MU/min', '210MU/min', '314MU/min', '417MU/min', '524MU/min', '592MU/min', '600MU/min')
test_types = ('drgs', 'drmlc')
im_types = ('open', 'mlc')

class VMAT(SingleImageObject):
    """The VMAT class analyzes two DICOM images acquired via a linac's EPID and analyzes
        regions of interest (segments) based on the paper by `Jorgensen et al <http://dx.doi.org/10.1118/1.3552922>`_, specifically,
        the Dose Rate & Gantry Speed (DRGS) and Dose Rate & MLC speed (DRMLC) tests.
        """

    def __init__(self):
        SingleImageObject.__init__(self)
        delattr(self,'image')  # delete the SIO attr 'image'. This is because we'll have two images (see below)
            # and it'd be good not to be confused by the presence of image in the class.
        self.image_open = np.array([], dtype=float)  # the Open field image
        self.image_dmlc = np.array([], dtype=float)  # the MLC field image
        self._test_type = ''  # the test to perform
        self._tolerance = 3 # default of 3% tolerance as Jorgensen recommends
        self._sample_ratios = np.array([], dtype=float)  # a float array holding the ratios (DMLC/Open) of MLC pair samples
        self._samples_passed = np.array([], dtype=bool)  # array of booleans specifying whether each sample passed or failed
        self._passed_analysis = False  # a boolean specifying a final flag about about whether or not the test passed completely.
        # the following are dicts holding the individual segment max, min and mean deviation
        self._seg_dev_max = {}
        self._seg_dev_min = {}
        self._seg_dev_mean = {}

    @value_accept(im_type=im_types)
    def load_image_UI(self, im_type='open'):
        """Open a UI file browser to load dicom image for given VMAT image.

        :param im_type: 'open' or 'mlc'; specifies which file/image is being loaded
        :type im_type: str
        """

        if  im_type == im_types[0]:  # open
            caption = "Select Open Field EPID Image..."
        else:  # dmlc
            caption = "Select MLC Field EPID Image..."

        fs = SingleImageObject._get_imagepath_UI(self, caption=caption)
        if fs:  # if user didn't hit cancel
            self.load_image(fs, im_type=im_type)

    @value_accept(im_types=im_types)
    def load_image(self, filepath, im_type='open'):
        """Load the image directly by the file path (i.e. non-interactively).

        :param filepath: The absolute path to the DICOM image
        :type filepath: str
        :param im_type: Specifies whether the image is the Open ('open') or DMLC ('dmlc') field.
        :type im_type: str
        """

        img, props = SingleImageObject.load_image(self, filepath, return_it=True)
        if im_type == im_types[0]:  # open
            self.image_open = img
        else:  # dmlc
            self.image_dmlc = img
        self.im_props = props

    @value_accept(test_types=test_types, number=(1, 2))
    def load_demo_image(self, test_type='drgs', number=1):
        """Load the demo DICOM images from demo files folder.

        :param test_type: Test type of images to load; 'drgs' or 'drmlc'
        :type test_type: str
        :param number: Which demo images to load; there are currently 2 for each test type.
        :type number: int
        """
        demo_folder = osp.join(osp.split(osp.abspath(__file__))[0], "demo_files")
        if test_type == test_types[1]:  # DRMLC
            if number == 1:
                im_open_path = osp.join(demo_folder, "DRMLCopen-example.dcm")
                im_dmlc_path = osp.join(demo_folder, 'DRMLCmlc-example.dcm')
            else:
                im_open_path = osp.join(demo_folder, "DRMLCopen-150-example.dcm")
                im_dmlc_path = osp.join(demo_folder, 'DRMLCmlc-150-example.dcm')
        else:
            if number == 1:
                im_open_path = osp.join(demo_folder, "DRGSopen-example.dcm")
                im_dmlc_path = osp.join(demo_folder, 'DRGSmlc-example.dcm')
            else:
                im_open_path = osp.join(demo_folder, "DRGSopen-150-example.dcm")
                im_dmlc_path = osp.join(demo_folder, 'DRGSmlc-150-example.dcm')

        self.load_image(im_open_path, im_type='open')
        self.load_image(im_dmlc_path, im_type='mlc')

    @value_accept(number=(1, 2))
    def run_demo_drgs(self, number=1):
        """Run the demo of the module for the Dose Rate & Gantry Speed test.

        :param number: Which demo image set to run; there are currently 2 demos.
        :type number: int
        """
        self.load_demo_image('drgs', number)
        self.analyze(test='', tolerance=3)  # tolerance at 2 to show some failures
        print(self.get_string_results())
        self.plot_analyzed_image()

    @value_accept(number=(1, 2))
    def run_demo_drmlc(self, number=1):
        """Run the demo of the module for the Dose Rate & MLC speed test.

        :param number: Which demo image set to run; there are currently 2 demos.
        :type number: int
        """
        self.load_demo_image('drmlc', number)
        self.analyze(test='drmlc', tolerance=3)
        print(self.get_string_results())
        self.plot_analyzed_image()

    def _get_im_scaling_factors(self, SID):
        """Determine the image scaling factors; factors are relative to reference values from images of size 384x512 taken at 150cm SID."""
        # Image size scaling
        y_scale = np.size(self.image_dmlc, 0) / 384.0
        x_scale = np.size(self.image_dmlc, 1) / 512.0
        # center pixel points
        x_im_center = np.size(self.image_dmlc, 1) / 2.0 + 0.5
        y_im_center = np.size(self.image_dmlc, 0) / 2.0 + 0.5
        # SID scaling
        if SID:
            SID_scale = SID / 150
        else:
            try:
                SID_scale = self.im_props['SID mm'] / 1500.0
            except:
                SID_scale = 1

        return SID_scale, x_im_center, x_scale, y_im_center, y_scale

    def _calc_sample_bounds(self, SID_scale, test, x_im_center, x_scale, y_im_center, y_scale):
        """Calculate the x and y pixel bounds of the samples to be extracted."""
        self._num_leaves = 38  # TODO: add if statement for HDMLCs
        sample_width_pix = 6  # the width in pixels of a given MLC leaf sample perpendicular to MLC movement at 150cm SID.
        leaf_width = 9.57  # the width from one sample to the next perpendicular to the MLC movement (MLC leaf to the next MLC leaf).
            # Calculation using DPmm of reference image gives 9.57
        y_pixel_bounds = [(np.round((12 + leaf * leaf_width - 192) * y_scale * SID_scale + y_im_center),
                           np.round((12 + leaf * leaf_width - 192 + sample_width_pix) * y_scale * SID_scale + y_im_center))
                          for leaf in np.arange(self._num_leaves)]

        # set up DRMLC sample points
        if test == test_types[1]:
            width = np.round(42 * x_scale * SID_scale)  # length of the sample parallel to MLC movement.
            num_segments = 4
            x_pixel_offsets = (-107, -49, 9, 68)  # offsets from the central x-pixel which are left-side starting points for the ROIs
            x_pixel_starts = [np.round(x_im_center + (offset * x_scale * SID_scale)) for offset in
                              x_pixel_offsets]  # multiply by scaling factor
            x_pixel_bounds = [(item, item + width) for item in x_pixel_starts]  # create a list of the left/right x-points
        # set up DRGS sample points
        else:
            width = np.round(24 * x_scale * SID_scale)  # length of the sample parallel to MLC movement.
            num_segments = 7
            x_pixel_offsets = (-109, -71, -33, 6, 43, 83, 122)
            x_pixel_starts = [np.round(x_im_center + (offset * x_scale * SID_scale)) for offset in
                              x_pixel_offsets]  # multiply by scaling factors
            x_pixel_bounds = [(x_start, x_start + width) for x_start in x_pixel_starts]

        self._num_segments = num_segments
        self._width = width
        self._sample_height = sample_width_pix*y_scale*SID_scale
        self._y_pixel_bounds = y_pixel_bounds
        self._x_pixel_bounds = x_pixel_bounds

    def _extract_samples(self):
        """Extract the mean of the pixel values of the samples from both Open and DMLC images based on pixel bounds."""
        # preallocation
        dmlc_samples = np.zeros((self._num_leaves, self._num_segments))
        open_samples = np.zeros((self._num_leaves, self._num_segments))
        # this pulls the samples from the Open and MLC image into new arrays
        # which are analyzed for the mean and standard deviation.
        for segment in np.arange(self._num_segments):
            for sample in np.arange(self._num_leaves):
                # extract ROI values
                dmlc_samples[sample, segment] = self.image_dmlc[self._y_pixel_bounds[sample][0]:self._y_pixel_bounds[sample][1],
                                                                self._x_pixel_bounds[segment][0]:self._x_pixel_bounds[segment][1]].mean()
                open_samples[sample, segment] = self.image_open[self._y_pixel_bounds[sample][0]:self._y_pixel_bounds[sample][1],
                                                                self._x_pixel_bounds[segment][0]:self._x_pixel_bounds[segment][1]].mean()

        return dmlc_samples, open_samples

    def _ratio_samples(self, dmlc_samples, open_samples):
        """Ratio the open and dmlc sample values to each other, i.e. the DMLC to open field ratio."""
        sample_ratios = dmlc_samples / open_samples
        self._sample_ratios = sample_ratios

    def _normalize_samples(self, dmlc_samples, open_samples):
        """Normalize the samples for each image respectively to the mean of the central segment."""
        # DRGS
        if self._test_type == test_types[0]:
            open_norm = open_samples[:, 3].mean()
            dmlc_norm = dmlc_samples[:, 3].mean()
        # MLCS
        else:
            open_norm = open_samples[:, 1:2].mean()
            dmlc_norm = dmlc_samples[:, 1:2].mean()
        open_samples /= open_norm
        dmlc_samples /= dmlc_norm
        return dmlc_samples, open_samples

    def _calc_deviations(self):
        """Calculate the deviations of each sample using Jorgensen's deviation equation."""
        # preallocation
        dev_matrix = np.zeros(self._sample_ratios.shape)
        # deviation calculation
        for sample in np.arange(self._num_leaves):
            denom = self._sample_ratios[sample, :].mean()
            for segment in np.arange(self._num_segments):
                numer = self._sample_ratios[sample, segment]
                dev_matrix[sample, segment] = numer / denom - 1
        # convert to percentages
        dev_matrix *= 100
        # saving the important values:
        # overall deviations
        self._dev_max = dev_matrix.max()
        self._dev_min = dev_matrix.min()
        self._dev_mean = dev_matrix.__abs__().mean()
        # segment deviations
        for segment in np.arange(self._num_segments):
            self._seg_dev_max[segment] = dev_matrix[:, segment].max()
            self._seg_dev_min[segment] = dev_matrix[:, segment].min()
            self._seg_dev_mean[segment] = dev_matrix[:, segment].__abs__().mean()

    @value_accept(test=test_types, tolerance=(0.3, 8), SID=(0, 180))
    def analyze(self, test, tolerance=3, SID=None):
        """Analyze 2 VMAT images, the open field image and DMLC field image, according to 1 of 2 possible tests.

        :param test: The test to perform, either Dose Rate Gantry Speed ('drgs') or Dose Rate MLC Speed ('drmlc').
        :type test: str
        :param tolerance: The tolerance of the sample deviations in percent. Default is 3, as Jorgensen recommends.
        :type tolerance: float, int
        :param SID: The Source to Image (detector) distance in cm. Usually doesn't need to be passed for EPID DICOM images. This argument
            will override any automatically derived value however. If left as None and no SID was determined, it will assume 150cm.
        :type SID: None, int
        """
        # error checking
        if len(self.image_open) == 0  or len(self.image_dmlc) == 0:
            raise AttributeError("Open or MLC Image not loaded yet. Use .load_image() or .load_image_UI()")

        self._tolerance = tolerance / 100.0
        self._test_type = test

        # get the image scaling factors and center pixels; this corrects for the SID
        SID_scale, x_im_center, x_scale, y_im_center, y_scale = self._get_im_scaling_factors(SID)

        # set up pixel bounds of test
        self._calc_sample_bounds(SID_scale, test, x_im_center, x_scale, y_im_center, y_scale)

        # Extract the samples from each image
        dmlc_samples, open_samples = self._extract_samples()

        # normalize the samples for each image respectively
        dmlc_samples, open_samples = self._normalize_samples(dmlc_samples, open_samples)

        # ratio the samples
        self._ratio_samples(dmlc_samples, open_samples)

        # calculate deviations as per Jorgensen equation
        self._calc_deviations()

        # run a pass/fail test on results
        self._run_passfail()

    def _run_passfail(self):
        """Checks whether the sample ratios passed the given tolerance."""

        # create boolean array the size of the number of samples; default to fail until pass
        self._samples_passed = np.zeros((self._num_leaves, self._num_segments), dtype=bool)

        # for each sample, if ratio is within tolerance, set to True
        #TODO: probably simple equality expression that's simpler.
        for segment in np.arange(self._num_segments):
            for sample in np.arange(self._num_leaves):
                if self._sample_ratios[sample,segment] < 1 + self._tolerance and self._sample_ratios[sample,segment] > 1 - self._tolerance:
                    self._samples_passed[sample,segment] = True

        # if all segments passed, set overall flag to pass
        if self._samples_passed.all():
            self._passed_analysis = True

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
        self.roi_handles = [[None for _ in range(self._num_segments)] for _ in range(self._num_leaves)]

        for segment in np.arange(self._num_segments):
            for sample in np.arange(self._num_leaves):
                if self._samples_passed[sample, segment]:
                    color = 'blue'
                else:
                    color = 'red'

                self.roi_handles[sample][segment] = (Rectangle((self._x_pixel_bounds[segment][0], self._y_pixel_bounds[sample][0]),
                                                  self._width, self._sample_height, fill=False, edgecolor=color, linewidth=1))
                plot.axes.add_patch(self.roi_handles[sample][segment])

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
            ax2.imshow(self.image_open, cmap=cm.Greys_r)
            ax2.set_title("Open Field Image")
        except:
            pass  # axis2 wasn't defined
        ax1.imshow(self.image_dmlc, cmap=cm.Greys_r)
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

        if self._passed_analysis:
            passfail_str = 'PASS'
        else:
            passfail_str = 'FAIL'

        if self._test_type == test_types[0]:  # DRGS
            string = ('Dose Rate & Gantry Speed \nTest Results (Tol. +/-%2.1f%%): %s\n' %
                      (self._tolerance * 100, passfail_str))
        else:  # MLC Speed
            string = ('Dose Rate & MLC Speed \nTest Results (Tol. +/-%2.1f%%): %s\n' %
                      (self._tolerance * 100, passfail_str))

        string += ('\nOverall Results:\n'
                   'Max Positive Deviation: %4.3f%%\n'
                   'Max Negative Deviation: %4.3f%%\n'
                   'Absolute Mean Deviation: %4.3f%%' %
                   (self._dev_max, self._dev_min, self._dev_mean))

        if self._test_type == test_types[0]:  # DRGS
            seg_str = DRGS_seg_strs
        else:
            seg_str = DRMLC_seg_strs
        for segment in np.arange(self._num_segments):
            string += ('\n\n%s Segment:\n'
                       'Max Positive Deviation: %4.3f%%\n'
                       'Max Negative Deviation: %4.3f%%\n'
                       'Absolute Mean Deviation: %4.3f%%' %
                       (seg_str[segment], self._seg_dev_max[segment], self._seg_dev_min[segment], self._seg_dev_mean[segment]))

        return string


#---------------------------------------------------------------------------------------------------------------------
# VMAT demo.
#---------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    v = VMAT()
    # v.load_image_UI()
    # v.load_image_UI('dmlc')
    v.load_demo_image('drmlc')
    v.analyze('drmlc', 2)
    print(v.get_string_results())
    v.plot_analyzed_image()
    # VMAT().run_demo_drmlc()  # uncomment to run MLCS demo