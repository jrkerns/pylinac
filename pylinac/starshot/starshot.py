# -*- coding: utf-8 -*-
"""
The Starshot module analyses a starshot film or multiple superimposed EPID images that measures the wobble of the
radiation spokes, whether gantry, collimator, or couch. It is based on ideas from `Depuydt et al <http://iopscience.iop.org/0031-9155/57/10/2997>`_
and `Gonzalez et al <http://dx.doi.org/10.1118/1.1755491>`_ and evolutionary optimization.
"""
from __future__ import division, print_function, absolute_import
from future import standard_library


standard_library.install_aliases()
import os
import os.path as osp
import zipfile as zp

import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from matplotlib.patches import Circle as mpl_Circle

from pylinac.core.common_functions import Prof_Penum, point2edge_min, peak_detect
from pylinac.core.decorators import value_accept
from pylinac.core.image_classes import ImageObj, AnalysisModule
from pylinac.core.geometry import Point, Line, Circle


class Starshot(AnalysisModule):
    """Creates a Starshot instance for determining the wobble in a gantry, collimator,
    couch or MLC starshot image pattern.
    """
    def __init__(self):
        super().__init__()
        self.image = ImageObj()  # The image array and property class
        self._algo_startpoint = Point()  # Point which specifies the algorithm starting point for search algorithm
        self.circle_profile = CircleProfile()
        self.lines = []  # list which will hold Line instances representing radiation lines.
        self.wobble = Wobble()
        self.tolerance = 1  # tolerance limit of the radiation wobble
        self.tolerance_unit = 'pixels'  # tolerance units are initially pixels. Will be converted to 'mm' if conversion
        # information available in image properties

    def load_demo_image(self):
        """Load the starshot demo image."""
        #TODO: see about using Python's temporary file/folder module
        demos_folder = osp.join(osp.split(osp.abspath(__file__))[0], 'demo_files')

        im_zip_path = osp.join(demos_folder, "demo_starshot_1.zip")

        # extract file from the zip file and put it in the demos folder
        zp.ZipFile(im_zip_path).extractall(demos_folder)
        # rename file path to the extracted one
        file_path = im_zip_path.replace('.zip', '.tif')
        # load image
        self.image.load_image(file_path)
        # delete extracted file to save space
        try:
            os.remove(file_path)
        except:
            print("Extracted demo image was not able to be deleted and remains in the demo directory")

    def load_image(self, filepath):
        """Load the image directly by passed the filepath."""
        self.image.load_image(filepath)

    def load_image_UI(self):
        """Load the image by using a UI dialog box."""
        self.image.load_image_UI()

    def set_start_point(self, point, warn_if_far_away=True):
        """Set the algorithm starting point manually.

        If the point given is far away (1% image width or 15 pixels) from the point found automatically a warning is given.

        :param point: The starting point for the search algorithm.
        :type point: Point
        :param warn_if_far_away: If the point is far away from the automatic determination, warn user.
        :type warn_if_far_away: boolean
        """
        if warn_if_far_away:
            # determine automatic start point if not yet set.
            if self._algo_startpoint.x == 0:
                self._auto_set_start_point()
            tolerance = max(min(self.image.pixel_array.shape)/100, 15)  # 1% image width of smalling dimension, or 15 pixels
            auto_y_upper = self._algo_startpoint.y - tolerance
            auto_y_lower = self._algo_startpoint.y + tolerance
            auto_x_left = self._algo_startpoint.x - tolerance
            auto_x_right = self._algo_startpoint.x + tolerance
            if (point.y < auto_y_upper or point.y > auto_y_lower) \
                    or (point.x < auto_x_left or point.x > auto_x_right):
                print("Warning: The point you've set is far away from the automatic calculation.\n" +
                      " The algorithm may not calculate correctly if you continue. \nUse method .clear_start_point" +
                      " to reset if need be or don't set the starting point manually.")

        self._algo_startpoint.x = point.x
        self._algo_startpoint.y = point.y

    def clear_start_point(self):
        """Clear/reset the algorithm starting point."""
        self._algo_startpoint = Point()

    def _check_image_inversion(self, allow_inversion=True):
        """Check the image for proper inversion, i.e. that pixel value increases with dose.

        Inversion is checked by the following:
        - Summing the image along both horizontal and vertical directions.
        - If the maximum point of both horizontal and vertical is in the middle 1/3, the image is assumed to be correct.
        - Otherwise, invert the image.
        """
        if not allow_inversion:
            return

        # sum the image along each axis
        x_sum = np.sum(self.image.pixel_array, 0)
        y_sum = np.sum(self.image.pixel_array, 1)

        # determine the point of max value for each sum profile
        xmaxind = np.argmax(x_sum)
        ymaxind = np.argmax(y_sum)

        # If that maximum point isn't near the center (central 1/3), invert image.
        center_in_central_third = ((xmaxind > len(x_sum) / 3 and xmaxind < len(x_sum) * 2 / 3) and
                               (ymaxind > len(y_sum) / 3 and ymaxind < len(y_sum) * 2 / 3))

        if not center_in_central_third:
            self.image.invert_array()

    def _auto_set_start_point(self):
        """Set the algorithm starting point automatically.

        The determination of an automatic start point is accomplished by finding the Full-Width-80%-Max.
        Finding the maximum pixel does not consistently work, esp. in the presence of a pin prick. The
        FW80M is a more consistent metric for finding a good start point.
        """
        # sum the image along each axis
        x_prof = np.sum(self.image.pixel_array, 0)
        y_prof = np.sum(self.image.pixel_array, 1)

        # Calculate Full-Width, 80% Maximum
        x_point = Prof_Penum(x_prof).get_FWXM_center(80)
        y_point = Prof_Penum(y_prof).get_FWXM_center(80)
        center_point = Point(x_point, y_point)

        self.set_start_point(center_point, warn_if_far_away=False)


    @value_accept(radius=(5, 95), min_peak_height=(5, 80), SID=(0, 180))
    def analyze(self, allow_inversion=True, radius=50, min_peak_height=30, SID=0):
        """Analyze the starshot image.

         Analyze finds the minimum radius and center of a circle that touches all the lines
         (i.e. the wobble circle diameter and wobble center)

         :param allow_inversion: Specifies whether to let the algorithm automatically check the image for proper inversion. Analysis will
            likely fail without proper inversion. Use .invert_image() to manually invert.
         :type allow_inversion: boolean
         :param radius: Distance in % between starting point and closest image edge.
         :type radius: int, float
         :param min_peak_height: The percentage minimum height a peak must be to be considered a valid peak. A lower value catches
            radiation peaks that vary in magnitude (e.g. different MU delivered), but also could pick up noise. Raise if pixel values of
            strips are similar but noise is getting caught. Also try changing radius if noise is a problem.
         :type min_peak_height: int
         :param SID: The source to image distance in cm. If passed in, results will be scaled to 100cm. E.g. a wobble of
            3 pixels at an SID of 150cm will be presented as 2 pixels [3/(150/100)].
         :type SID: int
        """
        # error checking
        if self.image.pixel_array.size == 0:
            raise AttributeError("Starshot image not yet loaded")

        # check inversion
        self._check_image_inversion(allow_inversion)

        # set starting point automatically if not yet set
        if self._algo_startpoint.x == 0:
            self._auto_set_start_point()

        # extract the circle profile
        self.circle_profile.get_profile(self.image.pixel_array, radius, self._algo_startpoint)

        # find the radiation lines using the peaks of the profile
        self.lines = self.circle_profile.find_rad_lines(min_peak_height)

        # find the wobble
        self._find_wobble_2step(SID)

    def _scale_wobble(self, SID):
        # convert wobble to mm if possible
        if self.image.properties['DPmm'] != 0:
            self.tolerance_unit = 'mm'
            self.wobble.radius = self.wobble.radius_pix / self.image.properties['DPmm']
        if SID:
            self.wobble.radius /= SID / 100
            self.wobble.radius_pix /= SID / 100

    def _find_wobble_2step(self, SID):
        """Find the smallest radius ("wobble") and center of a circle that touches all the star lines.

        This is accomplished by two rounds of searching. The first round finds the radius and center down to
        the nearest pixel. The second round finds the center and radius down to sub-pixel precision using parameter scale.
        This methodology is faster than one round of searching at sub-pixel precision.
        """
        sp = self._algo_startpoint

        # first round of searching; this finds the circle to the nearest pixel
        normal_tolerance, normal_scale = 0.05, 1.0
        self._find_wobble(normal_tolerance, sp, normal_scale)

        # second round of searching; this finds the circle down to sub-pixel precision
        small_tolerance, small_scale = 0.001, 10.0
        self._find_wobble(small_tolerance, self.wobble.center, small_scale)

        # scale the wobble based on the SID
        self._scale_wobble(SID)

    def _find_wobble(self, tolerance, start_point, scale):
        """An iterative method that moves pixel by pixel to the point of minimum distance to all radiation lines.

        :param tolerance: The value the "outside" pixels must be within compared to the center pixel to stop the algorithm
        :type tolerance: float
        :param start_point: The starting point for the search algorithm.
        :type start_point: Point
        :param scale: The scale of the search in pixels. E.g. 0.1 searches to 0.1 pixel precision.
        :type scale: float, int
        """
        #TODO: use an optimization function instead of evolutionary search
        sp = start_point
        # init conditions; initialize a 3x3 "ones" matrix and make corner value 0 to start minimum distance search.
        distmax = np.ones((3, 3))
        distmax[0, 0] = 0

        # find min point within the given tolerance
        while np.any(distmax < distmax[1, 1] - tolerance):  # while any edge pixel value + tolerance is less than the center one...
            # find which pixel that is lower than center pixel
            min_idx = np.unravel_index(distmax.argmin(),distmax.shape)
            # set new starting point to min dist index point
            sp.y += (min_idx[0] - 1)/scale
            sp.x += (min_idx[1] - 1)/scale
            for x in np.arange(-1,2):
                for y in np.arange(-1,2):
                    point = Point(y=sp.y+(y/scale), x=sp.x+(x/scale))
                    distmax[y+1, x+1] = np.max([line.distance_to_point(point) for line in self.lines])
                    # distmax[y+1, x+1] = np.max([point_to_2point_line_dist(point, line) for line in self.lines])

        self.wobble.radius_pix = distmax[1,1]
        self.wobble.center = sp

    @property
    def passed(self):
        """Boolean specifying whether overall test passed or failed."""
        if self.wobble.radius * 2 < self.tolerance:
            return True
        else:
            return False

    def get_string_results(self):
        """Return the results of the analysis.

        :return string: A string with a statement of the minimum circle.
        """
        if self.passed:
            passfailstr = 'PASS'
        else:
            passfailstr = 'FAIL'

        string = ('\nResult: %s \n\n'
                  'The minimum circle that touches all the star lines has a radius of %4.3g %s. \n\n'
                  'The center of the minimum circle is at %4.1f, %4.1f') % (passfailstr, self.wobble.radius, self.tolerance_unit,
                                                                            self.wobble.center.x, self.wobble.center.y)
        return string

    def plot_analyzed_image(self, plot=None):
        """Draw the star lines, profile circle, and wobble circle on a matplotlib figure.

        :param plot: The plot to draw on. If None, will create a new one.
        :type plot: matplotlib.image.AxesImage
        """
        # plot image
        if plot is None:
            imgplot = plt.imshow(self.image.pixel_array)
        else:
            plot.axes.imshow(self.image.pixel_array)
            plot.axes.hold(True)
            imgplot = plot

        # plot radiation lines
        for line in self.lines:
            line.add_to_axes(imgplot.axes)

        # plot wobble circle
        self.wobble.add_to_axes(imgplot.axes)

        # plot profile circle
        self.circle_profile.add_to_figure(imgplot)

        # tighten plot around image
        imgplot.axes.autoscale(tight=True)

        # Finally, show it all
        if plot is None:
            plt.show()
        else:
            plot.draw()
            plot.axes.hold(False)

    def run_demo(self):
        """Run the Starshot module demo."""
        self.load_demo_image()
        self.analyze()
        print(self.get_string_results())
        self.plot_analyzed_image()


class Wobble(Circle):
    """A class that holds the wobble information of the Starshot analysis."""
    def __init__(self):

        super().__init__()
        self.radius_pix = 0  # The radius of the circle in pixels. For proper drawing of the circle on the plot.

    # def add_to_figure(self, fig, color='black'):
    #     """Plot the wobble circle to the figure."""
    #     fig.axes.add_patch(Circle((self.center.x, self.center.y), edgecolor=color, radius=self.radius_pix, fill=False))


class CircleProfile(object):
    """Holds and analyzes the circular profile found during starshot analysis."""
    def __init__(self):
        self.profile = np.ndarray  # the 1D circular profile
        self.x = np.ndarray  # the x-values of the circle
        self.y = np.ndarray  # the y-values of the circle
        self.center = Point()  # The center point of the circular profile
        self.radius_pix = 0  # the radius of the circle in pixels
        self.radius_perc = 0  # the radius of the circle in percentage of center to closest image edge

    def get_profile(self, image_array, radius, start_point):
        """Extracts values of a circular profile around the isocenter point atop the image matrix,
        later to be searched for peaks and such.
        """

        # find minimum distance from starting point to image edges
        mindist = point2edge_min(image_array, start_point)

        # create index and cos, sin points which will be the circle's rectilinear coordinates
        deg = np.arange(0, 360 - 0.01, 0.01)
        x = np.cos(np.deg2rad(deg)) * radius / 100 * mindist + start_point.x
        y = np.sin(np.deg2rad(deg)) * radius / 100 * mindist + start_point.y

        # Pull the values of the image along the y,x points defined above, creating a circular profile
        raw_profile = ndimage.map_coordinates(image_array, [y, x], order=0)
        filtered_profile = ndimage.median_filter(raw_profile, size=100)
        normalized_profile = filtered_profile - np.min(filtered_profile)

        # Roll the profile if needed
        # --------------------------
        # In order to properly find the peaks, the bounds of the circular profile must not be near a radiation strip.
        # If the profile's edge (0-index) is in the middle of a radiation strip, move it over so that it's not
        zero_ind = np.where(normalized_profile == 0)
        prof = np.roll(normalized_profile, -zero_ind[0][0])
        x = np.roll(x, -zero_ind[0][0])
        y = np.roll(y, -zero_ind[0][0])

        self.profile = prof
        self.x = x
        self.y = y
        self.center = start_point
        self.radius_pix = radius / 100 * mindist
        self.radius_perc = radius / 100

    def _roughly_find_peaks(self, min_peak_height):
        """Find the peaks using max-value search. Enforces an even number of peaks found."""

        # Find the positions of the max values
        min_peak_height = min_peak_height / 100
        min_peak_distance = len(self.profile) / 100 * 3  # 3-degree minimum distance
        max_vals, max_idxs = peak_detect(self.profile, threshold=min_peak_height,
                                         min_peak_width=min_peak_distance)
        # ensure the # of peaks found was even; every radiation "strip" should result in two peaks, one on either side of the isocenter.
        if len(max_vals) % 2 != 0 or len(max_vals) == 0:
            raise ValueError(
                "The algorithm found zero or an uneven number of radiation peaks. Ensure that the starting " \
                "point is correct and/or change the search radius. Sorry.")

        return max_idxs, max_vals

    def _finely_find_peaks(self, max_idxs, max_vals):
        """Find the "peaks" be calculating the FWHM-C, using the roughly-found peaks as a starting point."""
        # create a zero-array called strip_limits that holds the indices of the minimum between peaks.
        # In this way, we search the full-width half-max within the indices between any two indices of strip_limits
        # The first index of strip_limits is always 0 and the last is always 36,000 (or whatever the length of
        # self._circleprofile is).
        strip_limits = np.zeros(len(max_vals) + 1).astype(int)
        for i in np.arange(len(max_vals) - 1):
            strip_limits[i + 1] = (max_idxs[i + 1] - max_idxs[i]) / 2 + max_idxs[i]
        strip_limits[-1] = len(self.profile)
        # Now, create and fill an array called center_indices that will be the index of _circleprofile that the FWHM is at.
        center_indices = np.zeros(len(max_vals))
        # Determine the FWHM of each peak
        for i in np.arange(len(max_vals)):
            prof = Prof_Penum(self.profile[strip_limits[i]:strip_limits[i + 1]],
                              np.arange(strip_limits[i], strip_limits[i + 1]))
            center_indices[i] = prof.get_FWXM_center()
        center_indices = np.round(center_indices).astype(int)  # round to the nearest pixel

        # create a list of Points from the peak indices
        peak_list = [Point(self.x[i], self.y[i]) for i in center_indices]

        return peak_list

    def find_rad_lines(self, min_peak_height):
        """Find and match the positions of peaks in the circle profile and map their positions to the starshot image."""

        # find the peaks using simple max-height search
        max_idxs, max_vals = self._roughly_find_peaks(min_peak_height)

        # Using the above peaks as a starting point, find the FWHM-C
        peak_list = self._finely_find_peaks(max_idxs, max_vals)

        # Match the peaks to form lines
        lines = self._match_peaks(peak_list)
        return lines

    def _match_peaks(self, peak_list):
        """Match the peaks found to the same radiation lines.

        Peaks are matched by connecting the existing peaks based on an offset of peaks. E.g. if there are
        12 peaks, there must be 6 radiation lines. Furthermore, assuming star lines go all the way across the CAX,
        the 7th peak will be the opposite peak of the 1st peak, forming a line. This method is robust to
        starting points far away from the real center.
        """
        line_list = []
        offset = int(len(peak_list)/2)
        num_rad_lines = np.arange(len(peak_list) // 2)
        for line in num_rad_lines:
            line_list.append(Line(peak_list[line], peak_list[line + offset]))
        return line_list

    def add_to_figure(self, fig=None, color='green'):
        """Plot the circle profile that was extracted.

        :param fig: Figure to plot the profile to. If None, will create a new 1D figure plot.
        """
        if fig is None:
            plt.plot(self.profile)
            plt.show()
        else:
            fig.axes.add_patch(mpl_Circle((self.center.x, self.center.y), edgecolor=color, radius=self.radius_pix, fill=False))

# ----------------------------
# Starshot demo
# ----------------------------
if __name__ == '__main__':
    Starshot().run_demo()
    pass