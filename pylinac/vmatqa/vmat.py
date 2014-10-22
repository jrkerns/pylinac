# -*- coding: utf-8 -*-
"""The VMAT module consists of the class VMAT, which is capable of loading an EPID DICOM Open field image and MLC field image and analyzing the
images according to the `Jorgensen et al. <http://dx.doi.org/10.1118/1.3552922>`_ tests, specifically the Dose-Rate Gantry-Speed (DRGS) and Dose-Rate MLC (DRMLC) tests.

"""
# Full documentation of this module and how to use it can be found at this repositories' Read the Docs site: pylinac.rtfd.org

# builtins
from __future__ import print_function, division, absolute_import
import os.path as osp

# 3rd party
from numpy import mean, arange, array, zeros
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Rectangle

# internal
from pylinac.common.image_classes import SingleImageObject


# Common constants
DRMLC_seg_strs = ('1.6cm/s', '2.4cm/s', '0.8cm/s', '0.4cm/s')
DRGS_seg_strs = ('105MU/min', '210MU/min', '314MU/min', '417MU/min', '524MU/min', '592MU/min', '600MU/min')

class VMAT(SingleImageObject):
    """
        The VMAT class analyzes two DICOM images acquired via a linac's EPID and analyzes
        regions of interest (segments) based on the paper by Jorgensen et al (http://dx.doi.org/10.1118/1.3552922), specifically,
        the dose-rate gantry speed (DRGS) and MLC speed (DRMLC) tests.
        """

    def __init__(self):
        SingleImageObject.__init__(self)
        delattr(self,'image')  # delete the SIO attr 'image'. This is because we'll have two images (see below)
            # and it'd be good not to be confused by the presence of image in the class.
        self.image_open = None # the Open field image
        self.image_mlc = None # the MLC field image
        self._test_type = None  # the test to perform
        self._tolerance = 2 # default of 3% tolerance as Jorgensen recommends
        self._sample_ratios = np.array([])
        self._samples_passed = np.array([])  # list of booleans specifying whether each ROI of the test passed or failed
        self._test_types = ('drgs', 'drmlc')
        self._im_types = ('open', 'mlc')
        self._passed_analysis = False  # a boolean specifying a final flag about about whether or not the test passed completely
        self._segment_medians = np.array([])
        self._segment_sts = np.array([])

    def load_image_UI(self, im_type='open'):
        """Open a Qt UI file browser to load dicom image for given VMAT image.

        :param im_type: 'open' or 'mlc'; specifies which file/image is being loaded
        :type im_type: str
        """

        if  im_type == self._im_types[0]:  # open
            caption = "Select Open Field EPID Image..."
        elif im_type == self._im_types[1]:  # dmlc
            caption = "Select MLC Field EPID Image..."
        else:
            raise NameError("im_type input string {s} not valid".format(im_type))

        fs = SingleImageObject.get_imagepath_UI(self, caption=caption)
        if fs:  # if user didn't hit cancel
            self.load_image(fs, im_type=im_type)

    def load_image(self, filepath, im_type='open'):
        """Load the image directly by the file path (i.e. non-interactively).
        :param filepath: The absolute path to the DICOM image
        :type filepath: str
        :param im_type: Specifies whether the image is the Open ('open') or MLC ('mlc') field.
        :type im_type: str
        """

        img, props = SingleImageObject.load_image(self, filepath, return_it=True)
        if im_type == self._im_types[0]:  # open
            self.image_open = img
        elif im_type == self._im_types[1]:  # dmlc
            self.image_mlc = img
        else:
            raise NameError("im_type input string {s} not valid".format(im_type))
        self.im_props = props

    def load_demo_image(self, test_type='drgs', number=1):
        """
        Load the demo DICOM images from demo files folder
        """
        assert test_type in self._test_types, "test_type input string was not valid. Must be 'drmlc', or 'drgs'"

        if test_type == self._test_types[1]:  # DRMLC
            if number == 1:
                im_open_path = osp.join(osp.split(osp.abspath(__file__))[0], "demo files", "DRMLCopen-example.dcm")
                im_dmlc_path = osp.join(osp.split(osp.abspath(__file__))[0], 'demo files', 'DRMLCmlc-example.dcm')
            elif number == 2:
                im_open_path = osp.join(osp.split(osp.abspath(__file__))[0], "demo files", "DRMLCopen-150-example.dcm")
                im_dmlc_path = osp.join(osp.split(osp.abspath(__file__))[0], 'demo files', 'DRMLCmlc-150-example.dcm')
        else:
            if number == 1:
                im_open_path = osp.join(osp.split(osp.abspath(__file__))[0], "demo files", "DRGSopen-example.dcm")
                im_dmlc_path = osp.join(osp.split(osp.abspath(__file__))[0], 'demo files', 'DRGSmlc-example.dcm')
            elif number == 2:
                im_open_path = osp.join(osp.split(osp.abspath(__file__))[0], "demo files", "DRGSopen-150-example.dcm")
                im_dmlc_path = osp.join(osp.split(osp.abspath(__file__))[0], 'demo files', 'DRGSmlc-150-example.dcm')

        self.load_image(im_open_path, im_type='open')
        self.load_image(im_dmlc_path, im_type='mlc')

    def run_demo_drgs(self, number=1):
        """Run the demo of the module for the Dose Rate & Gantry Speed test."""
        self.load_demo_image('drgs', number)
        self.analyze(test='drgs', tolerance=3)  # tolerance at 2 to show some failures
        print(self.get_string_results())
        self.show_img_results()

    def run_demo_drmlc(self, number=1):
        """Run the demo of the module for the MLC speed test."""
        self.load_demo_image('drmlc', number)
        self.analyze(test='drmlc', tolerance=3)
        print(self.get_string_results())
        self.show_img_results()

    def _get_im_scaling_factors(self):
        """Determine the image scaling factors; initial values were based on 384x512 images @ 150cm SID."""
        # Image size scaling
        y_scale = np.size(self.image_mlc, 0) / 384.0
        x_scale = np.size(self.image_mlc, 1) / 512.0
        # center pixel points
        x_im_center = np.size(self.image_mlc, 1) / 2.0 + 0.5
        y_im_center = np.size(self.image_mlc, 0) / 2.0 + 0.5
        # SID scaling
        try:
            SID_scale = self.im_props['SID mm'] / 1500.0
        except:
            SID_scale = 1

        return SID_scale, x_im_center, x_scale, y_im_center, y_scale

    def analyze(self, test, tolerance=3):
        """Analyze 2 VMAT images, the open field image and MLC field image, according to 1 of 2 possible tests.
        """
        # error checking that images are loaded. If using in PyTG142, warn user. If using standalone, assert error
        if self.image_open is None or self.image_mlc is None:
            raise AttributeError("Open or MLC Image not loaded yet. Use .load_image()")
        if test.lower() not in self._test_types:
            raise NameError("Test input string was not valid. Must be 'drmlc', or 'drgs'")
        if (type(tolerance) != float and type(tolerance) != int)or tolerance < 0.1:
            raise ValueError("Tolerance must be a float or int greater than 0.1")

        self._tolerance = tolerance / 100.0
        self._test_type = test

        # get the image scaling factors and center pixels
        SID_scale, x_im_center, x_scale, y_im_center, y_scale = self._get_im_scaling_factors()

        # set up ROI points used in both tests
        num_leaves = 38
        sample_width_pix = 6 # the width in pixels of a given MLC leaf sample
        leaf_width = 9.57 # (368-12)/38  # the width from sample to sample (MLC leaf to the next MLC leaf). Calculation using DPmm gives
        # 9.57

        y_pixel_bounds = [(np.round((12+leaf*leaf_width-192)*y_scale*SID_scale+y_im_center),
                           np.round((12+leaf*leaf_width-192+sample_width_pix)*y_scale*SID_scale+y_im_center))
                            for leaf in arange(num_leaves)]

        # set up MLC speed sample points
        if test == self._test_types[1]:
            width, num_segments = np.round(42 * x_scale * SID_scale), 4
            x_pixel_offsets = (-107, -49, 9, 68)  # offsets from the central x-pixel which are left-side starting points for the ROIs
            x_pixel_starts = [np.round(x_im_center + (offset*x_scale*SID_scale)) for offset in x_pixel_offsets] # multiply by scaling factor
            x_pixel_bounds = [(item, item + width) for item in x_pixel_starts]  # create a list of the left/right x-points

        # set up Dose-Rate Gantry-Speed sample points
        else:
            width, num_segments = np.round(24 * x_scale * SID_scale), 7
            x_pixel_offsets = np.array((147, 185, 223, 262, 299, 339, 378)) - 256
            x_pixel_starts = [np.round(x_im_center + (item * x_scale*SID_scale)) for item in x_pixel_offsets]  # multiply by scaling factors
            x_pixel_bounds = [(item, item + width) for item in x_pixel_starts]

        #preallocation/instantiation
        dmlc_samples = zeros((num_leaves, num_segments))
        open_samples = zeros((num_leaves, num_segments))

        # this pulls the samples from the Open and MLC image into new arrays
        # which are analyzed for the mean and standard deviation.
        for segment in arange(num_segments):
            for sample in arange(num_leaves):
                # extract ROI values
                dmlc_samples[sample, segment] = mean(array(self.image_mlc[y_pixel_bounds[sample][0]:y_pixel_bounds[sample][1],
                                                           x_pixel_bounds[segment][0]:x_pixel_bounds[segment][1]]))
                open_samples[sample, segment] = mean(array(self.image_open[y_pixel_bounds[sample][0]:y_pixel_bounds[sample][1],
                                                           x_pixel_bounds[segment][0]:x_pixel_bounds[segment][1]]))

        # normalize the samples to that of segment 4 for DRGS or the mean of segments 2 and 3 for MLCS of their own respective image
        if self._test_type == self._test_types[0]:  # DRGS
            open_norm = mean(open_samples[:,3])
            dmlc_norm = mean(dmlc_samples[:,3])
        else:  # MLCS
            open_norm = mean(open_samples[:, 1:2])
            dmlc_norm = mean(dmlc_samples[:, 1:2])

        open_samples /= open_norm
        dmlc_samples /= dmlc_norm

        # now ratio the open and mlc sample values to each other, i.e. the MLC to open field difference
        sample_ratios = dmlc_samples / open_samples

        # calculate deviations as per Jorgensen equation
        dev_matrix = np.zeros(sample_ratios.shape)
        for sample in arange(num_leaves):
            denom = sample_ratios[sample,:].mean()
            for segment in arange(num_segments):
                numer = sample_ratios[sample, segment]
                dev_matrix[sample, segment] = numer/denom - 1
        # convert to percentages
        dev_matrix *= 100
        self._dev_max = dev_matrix.max()
        self._dev_min = dev_matrix.min()
        self._dev_mean = dev_matrix.__abs__().mean()

        # attach other sample values to self
        self._sample_height = sample_width_pix*y_scale*SID_scale
        self._sample_ratios = sample_ratios
        self._width = width
        self._y_pixel_bounds = y_pixel_bounds
        self._x_pixel_bounds = x_pixel_bounds
        self._num_segments = num_segments
        self._num_leaves = num_leaves

        # run a pass/fail test on results
        self._run_passfail()

    def _run_passfail(self):
        """Checks whether the sample ratios passed the given tolerance.
        """

        # create boolean array the size of the number of samples; default to fail until pass
        self._samples_passed = np.zeros((self._num_leaves, self._num_segments), dtype=bool)

        # for each sample, if ratio is within tolerance, set to True
        #TODO: probably simple equality expression that's simpler.
        for segment in arange(self._num_segments):
            for sample in arange(self._num_leaves):
                if self._sample_ratios[sample,segment] < 1 + self._tolerance and self._sample_ratios[sample,segment] > 1 - self._tolerance:
                    self._samples_passed[sample,segment] = True

        # if all segments passed, set overall flag to pass
        if self._samples_passed.all():
            self._passed_analysis = True

    def _draw_objects(self, plot, draw=False):
        """
        This method draws lines on the matplotlib widget/figure showing an outline of the ROIs used in
        the calcuation and text of the results

        @type plot: matplotlib.axes.Axes
        """

        # clear images of any existing objects that are on it
        if hasattr(self, 'roi_handles'):
            for handle in self.roi_handles:
                handle.remove()
            if draw:
                plot.draw()

        #plot ROI lines on image
        self.roi_handles = [[None for _ in range(self._num_segments)] for _ in range(self._num_leaves)]

        for segment in arange(self._num_segments):
            for sample in arange(self._num_leaves):
                if self._samples_passed[sample, segment]:
                    color = 'blue'
                else:
                    color = 'red'

                self.roi_handles[sample][segment] = (Rectangle((self._x_pixel_bounds[segment][0], self._y_pixel_bounds[sample][0]),
                                                  self._width, self._sample_height, fill=False, edgecolor=color, linewidth=1))
                plot.axes.add_patch(self.roi_handles[sample][segment])

        if draw:
            plot.draw()

    def show_img_results(self):
        """Create 2 plots showing the open and MLC images with the samples and results drawn on.
        """
        fig, (ax1, ax2) = plt.subplots(1,2)
        ax1.imshow(self.image_open, cmap=cm.Greys_r)  # http://stackoverflow.com/questions/3823752/display-image-as-grayscale-using-matplotlib
        ax2.imshow(self.image_mlc, cmap=cm.Greys_r)
        ax1.set_title("Open Field Image")
        ax2.set_title("MLC Field Image")
        self._draw_objects(ax2)
        plt.show()

    def get_string_results(self):
        """Create a string of the summary of the analysis results."""

        if self._passed_analysis:
            passfail_str = 'PASS'
        else:
            passfail_str = 'FAIL'

        if self._test_type == self._test_types[0]:  # DRGS
            string = ('Dose Rate & Gantry Speed \nTest Results (Tol. +/-%2.1f%%): %s\n'
                      'Max Positive Deviation: %4.3f%%\n'
                      'Max Negative Deviation: %4.3f%%\n'
                      'Absolute Mean Deviation: %4.3f%%' %
                      (self._tolerance * 100, passfail_str, self._dev_max, self._dev_min, self._dev_mean))
        else:  # MLC Speed
            string = ('Dose Rate & MLC Speed \nTest Results (Tol. +/-%2.1f%%): %s\n'
                      'Max Positive Deviation: %4.3f%%\n'
                      'Max Negative Deviation: %4.3f%%\n'
                      'Absolute Mean Deviation: %4.3f%%' %
                      (self._tolerance * 100, passfail_str, self._dev_max, self._dev_min, self._dev_mean))

        return string


#---------------------------------------------------------------------------------------------------------------------
# VMAT demo.
#---------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    VMAT().run_demo_drgs(1)
    # VMAT().run_demo_drmlc()  # uncomment to run MLCS demo