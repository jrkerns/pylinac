# -*- coding: utf-8 -*-
"""
The Starshot module analyses a starshot film or multiple superimposed EPID images that measures the wobble of the
radiation spokes, whether gantry, collimator, or couch. It is based on ideas from `Depuydt et al <http://iopscience.iop.org/0031-9155/57/10/2997>`_
and `Gonzalez et al <http://dx.doi.org/10.1118/1.1755491>`_ and evolutionary optimization.
"""

import os
import os.path as osp
import zipfile as zp

import numpy as np
import matplotlib.pyplot as plt

from pylinac.core.decorators import value_accept
from pylinac.core.geometry import Point, Line, Circle
from pylinac.core.image import ImageObj
from pylinac.core.analysis import AnalysisModule
from pylinac.core.profile import CircleProfile, SingleProfile


class Starshot(AnalysisModule):
    """Creates a Starshot instance for determining the wobble in a gantry, collimator,
    couch or MLC starshot image pattern.
    """
    def __init__(self):
        super().__init__()
        self.image = ImageObj()  # The image array and property class
        self.circle_profile = StarProfile()
        self.lines = []  # list which will hold Line instances representing radiation lines.
        self.wobble = Wobble()
        self.tolerance = 1  # tolerance limit of the radiation wobble
        self.tolerance_unit = 'pixels'  # tolerance units are initially pixels. Will be converted to 'mm' if conversion
        # information available in image properties

    def load_demo_image(self):
        """Load the starshot demo image."""
        #TODO: see about using Python's temporary file/folder module
        demos_folder = osp.join(osp.dirname(__file__), 'demo_files', 'starshot')

        im_zip_path = osp.join(demos_folder, "starshot_gantry.zip")

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

    # @property
    # def start_point(self):
    #     return self.circle_profile.center

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
            if self.circle_profile.center.x == 0:
                self._auto_set_start_point()
            tolerance = max(min(self.image.pixel_array.shape)/100, 15)  # 1% image width of smalling dimension, or 15 pixels
            auto_y_upper = self.circle_profile.center.y - tolerance
            auto_y_lower = self.circle_profile.center.y + tolerance
            auto_x_left = self.circle_profile.center.x - tolerance
            auto_x_right = self.circle_profile.center.x + tolerance
            if (point.y < auto_y_upper or point.y > auto_y_lower) \
                    or (point.x < auto_x_left or point.x > auto_x_right):
                print("Warning: The point you've set is far away from the automatic calculation.\n" +
                      " The algorithm may not calculate correctly if you continue. \nUse method .clear_start_point" +
                      " to reset if need be or don't set the starting point manually.")

        self.circle_profile.center = point

    def clear_start_point(self):
        """Clear/reset the algorithm starting point."""
        self.circle_profile.center = Point()

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
        x_sum = np.sum(self.image.pixel_array, 0)
        y_sum = np.sum(self.image.pixel_array, 1)

        # Calculate Full-Width, 80% Maximum
        x_point = SingleProfile(x_sum).get_FWXM_center(80)
        y_point = SingleProfile(y_sum).get_FWXM_center(80)
        center_point = Point(x_point, y_point)

        self.set_start_point(center_point, warn_if_far_away=False)

    @value_accept(radius=(0.05, 0.95), min_peak_height=(5, 80), SID=(0, 180))
    def analyze(self, allow_inversion=True, radius=0.5, min_peak_height=0.25, SID=0):
        """Analyze the starshot image.

         Analyze finds the minimum radius and center of a circle that touches all the lines
         (i.e. the wobble circle diameter and wobble center)

         :param allow_inversion: Specifies whether to let the algorithm automatically check the image for proper inversion. Analysis will
            likely fail without proper inversion. Use .invert_image() to manually invert.
         :type allow_inversion: boolean
         :param radius: Distance in % between starting point and closest image edge. Must be between 0.05 and 0.95.
         :type radius: float
         :param min_peak_height: The percentage minimum height a peak must be to be considered a valid peak. A lower value catches
            radiation peaks that vary in magnitude (e.g. different MU delivered), but also could pick up noise. Raise if pixel values of
            strips are similar but noise is getting caught. Also try changing radius if noise is a problem.
         :type min_peak_height: int
         :param SID: The source to image distance in cm. If passed in, results will be scaled to 100cm. E.g. a wobble of
            3 pixels at an SID of 150cm will be presented as 2 pixels [3/(150/100)].
         :type SID: int
        """
        # error checking
        if not self.image_is_loaded:
            raise AttributeError("Starshot image not yet loaded")

        # check inversion
        self._check_image_inversion(allow_inversion)

        # set starting point automatically if not yet set
        if not self.start_point_is_set:
            self._auto_set_start_point()

        # set profile extraction radius
        self.circle_profile.radius = self.convert_radius_perc2pix(radius)
        # extract the circle profile
        self.circle_profile.get_profile(self.image.pixel_array)
        # find the radiation lines using the peaks of the profile
        self.lines = self.circle_profile.find_rad_lines(min_peak_height)

        # find the wobble
        self._find_wobble_2step(SID)

    def convert_radius_perc2pix(self, radius):
        """Convert a radius in percent (e.g. 50) to distance in pixels, based on the distance from center point to image edge."""
        dist = self.image.dist2edge_min(self.circle_profile.center)
        return dist*radius

    @property
    def image_is_loaded(self):
        if self.image.pixel_array.size == 0:
            return False
        else:
            return True

    @property
    def start_point_is_set(self):
        if self.circle_profile.center.x == 0:
            return False
        else:
            return True

    def _scale_wobble(self, SID):
        # convert wobble to mm if possible
        if self.image.properties['DPmm'] != 0:
            self.tolerance_unit = 'mm'
            self.wobble.radius_mm = self.wobble.radius / self.image.properties['DPmm']
        else:
            self.tolerance_unit = 'pixels'
            self.wobble.radius_mm = self.wobble.radius
        if SID:
            self.wobble.radius /= SID / 100
            self.wobble.radius_mm /= SID / 100

    def _find_wobble_2step(self, SID):
        """Find the smallest radius ("wobble") and center of a circle that touches all the star lines.

        This is accomplished by two rounds of searching. The first round finds the radius and center down to
        the nearest pixel. The second round finds the center and radius down to sub-pixel precision using parameter scale.
        This methodology is faster than one round of searching at sub-pixel precision.
        """
        sp = self.circle_profile.center

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
                    distmax[y+1, x+1] = np.max([line.distance_to(point) for line in self.lines])

        self.wobble.radius = distmax[1,1]
        self.wobble.center = sp
        self.circle_profile.center = sp

    @property
    def passed(self):
        """Boolean specifying whether overall test passed or failed."""
        if self.wobble.radius_mm * 2 < self.tolerance:
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
                  'The center of the minimum circle is at %4.1f, %4.1f') % (passfailstr, self.wobble.radius_mm, self.tolerance_unit,
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
        self.circle_profile.add_to_axes(imgplot.axes, edgecolor='green')

        # tighten plot around image
        imgplot.axes.autoscale(tight=True)

        # Finally, show it all
        if plot is None:
            plt.show()
        else:
            plot.draw()
            plot.axes.hold(False)

    def run_demo(self):
        """Run the Starshot module using the demo image."""
        self.load_demo_image()
        self.analyze()
        print(self.get_string_results())
        self.plot_analyzed_image()


class Wobble(Circle):
    """A class that holds the wobble information of the Starshot analysis."""
    def __init__(self, center_point=None, radius=None):
        super().__init__(center_point=center_point, radius=radius)
        self.radius_mm = 0  # The radius of the wobble in mm; as opposed to pixels.


class StarProfile(CircleProfile):
    """Holds and analyzes the circular profile used during starshot analysis."""
    def __init__(self):
        super().__init__()

    def get_profile(self, image_array):
        """Overload to also correct for profile positioning."""
        super().get_profile(image_array)
        self._roll_prof_to_midvalley()

    def _roll_prof_to_midvalley(self):
        """Roll the circle profile so that its edges are not near a radiation line. This is a prerequisite for properly finding star
        lines."""
        # Find index of the min value(s)
        min_idx = np.where(self.y_values == self.y_values.min())
        # Roll the profile and x and y coordinates
        self.y_values = np.roll(self.y_values, -min_idx[0][0])
        self.x_locs = np.roll(self.x_locs, -min_idx[0][0])
        self.y_locs = np.roll(self.y_locs, -min_idx[0][0])

    def find_rad_lines(self, min_peak_height):
        """Find and match the positions of peaks in the circle profile and map their positions to the starshot image."""

        # find the FWHM-C
        min_peak_distance = len(self.y_values) / 100 * 3  # 3-degree minimum distance between peaks
        self.find_peaks(min_peak_height, min_peak_distance)
        # self.find_FWHM_peaks(min_peak_height, min_peak_distance)

        # map 1D peaks to the image; i.e. add x and y coords of peaks
        self.map_peaks()

        # Match the peaks to form lines
        lines = self._match_peaks()
        return lines

    def _match_peaks(self):
        """Match the peaks found to the same radiation lines.

        Peaks are matched by connecting the existing peaks based on an offset of peaks. E.g. if there are
        12 peaks, there must be 6 radiation lines. Furthermore, assuming star lines go all the way across the CAX,
        the 7th peak will be the opposite peak of the 1st peak, forming a line. This method is robust to
        starting points far away from the real center.
        """
        # line_list = []
        offset = int(len(self.peaks)/2)
        num_rad_lines = len(self.peaks) // 2
        # for line in range(num_rad_lines):
        #     line_list.append(Line(self.peaks[line], self.peaks[line+offset]))
        line_list = [Line(self.peaks[line], self.peaks[line+offset]) for line in range(num_rad_lines)]
        return line_list

# ----------------------------
# Starshot demo
# ----------------------------
if __name__ == '__main__':
    Starshot().run_demo()
    pass