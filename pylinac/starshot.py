# -*- coding: utf-8 -*-
"""
The Starshot module analyses a starshot film or multiple superimposed EPID images that measures the wobble of the
radiation spokes, whether gantry, collimator, or couch. It is based on ideas from
`Depuydt et al <http://iopscience.iop.org/0031-9155/57/10/2997>`_
and `Gonzalez et al <http://dx.doi.org/10.1118/1.1755491>`_ and evolutionary optimization.
"""

import os.path as osp

import numpy as np
import matplotlib.pyplot as plt

from pylinac.core.decorators import value_accept
from pylinac.core.geometry import Point, Line, Circle
from pylinac.core.image import Image
from pylinac.core.analysis import AnalysisModule
from pylinac.core.profile import CircleProfile, SingleProfile


class Starshot(AnalysisModule):
    """Class that can determine the wobble in a "starshot" image, be it gantry, collimator,
        couch or MLC. The image can be DICOM or a scanned film (TIF, JPG, etc).

    Attributes
    ----------
    image : core.image.Image
    circle_profile : StarProfile
    lines : list of Line instances
    wobble : Wobble

    Examples
    --------
    Run the demo:
        >>> Starshot().run_demo()

    Typical session:
        >>> img_path = r"C:/QA/Starshots/Coll"
        >>> mystar = Starshot()
        >>> mystar.load_image(img_path)
        >>> mystar.analyze()
        >>> print(mystar.return_results())
        >>> mystar.plot_analyzed_image()
    """
    def __init__(self):
        super().__init__()
        # self.image = Image  # The image array and image property structure
        self.circle_profile = StarProfile()  # a circular profile which will detect radiation line locations
        self.lines = []  # a list which will hold Line instances representing radiation lines.
        self.wobble = Wobble()  # A Circle representing the radiation wobble
        self._tolerance = 1  # tolerance limit of the radiation wobble
        self._tolerance_unit = 'pixels'  # tolerance units are initially pixels. Will be converted to 'mm' if conversion
        # information available in image properties

    def load_demo_image(self):
        """Load the starshot demo image.

        The Pylinac package comes with compressed demo images.
        When called, the function unpacks the demo image and loads it.

        Parameters
        ----------
        cleanup : boolean
            If True (default), the extracted demo file is deleted (but not the compressed version).
            If False, leaves the extracted file alone after loading. Useful when using the demo image a lot,
            or you don't mind using the extra space.
        """
        demo_folder = osp.join(osp.dirname(__file__), 'demo_files', 'starshot')
        demo_file = osp.join(demo_folder, '10X_collimator.tif')

        self.image = Image(demo_file)

    def load_image(self, filepath):
        """Load the image via the file path.

        Parameters
        ----------
        filepath : str
            Path to the file to be loaded.

        Notes
        -----
        Wrapper for pylinac.image.ImageObj.load_image()
        """
        self.image = Image(filepath)

    def load_image_UI(self):
        """Load the image by using a UI dialog box."""
        self.image = Image.open_UI()

    @property
    def start_point(self):
        """The start point of the wobble search algorithm.

        After analysis this point is the wobble center.
        """
        return self.circle_profile.center

    def set_start_point(self, point, warn_if_far_away=True):
        """Set the algorithm starting point manually.

        Parameters
        ----------
        point : geometry.Point
            The starting point desired for the search algorithm.
        warn_if_far_away : boolean
            Flag to warn user if the starting point is far away
            (1% image width or 15 pixels) from the automatic determination.
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

    def _check_image_inversion(self):
        """Check the image for proper inversion, i.e. that pixel value increases with dose.

        Notes
        -----
        Inversion is checked by the following:
        - Summing the image along both horizontal and vertical directions.
        - If the maximum point of both horizontal and vertical is in the middle 1/3, the image is assumed to be correct.
        - Otherwise, invert the image.
        """

        # sum the image along each axis
        x_sum = np.sum(self.image.array, 0)
        y_sum = np.sum(self.image.array, 1)

        # determine the point of max value for each sum profile
        xmaxind = np.argmax(x_sum)
        ymaxind = np.argmax(y_sum)

        # If that maximum point isn't near the center (central 1/3), invert image.
        center_in_central_third = ((xmaxind > len(x_sum) / 3 and xmaxind < len(x_sum) * 2 / 3) and
                               (ymaxind > len(y_sum) / 3 and ymaxind < len(y_sum) * 2 / 3))

        if not center_in_central_third:
            self.image.invert()

    def _auto_set_start_point(self):
        """Set the algorithm starting point automatically.

        Notes
        -----
        The determination of an automatic start point is accomplished by finding the Full-Width-80%-Max.
        Finding the maximum pixel does not consistently work, esp. in the presence of a pin prick. The
        FW80M is a more consistent metric for finding a good start point.
        """
        # sum the image along each axis within the central 1/3 (avoids outlier influence from say, gantry shots)
        top_third = int(self.image.array.shape[0]/3)
        bottom_third = int(top_third * 2)
        left_third = int(self.image.array.shape[1]/3)
        right_third = int(left_third * 2)
        x_sum = np.sum(self.image.array[top_third:bottom_third, left_third:right_third], 0)
        y_sum = np.sum(self.image.array[top_third:bottom_third, left_third:right_third], 1)

        # Calculate Full-Width, 80% Maximum
        x_point = SingleProfile(x_sum).get_FWXM_center(80)
        y_point = SingleProfile(y_sum).get_FWXM_center(80)
        center_point = Point(x_point+left_third, y_point+top_third)

        self.set_start_point(center_point, warn_if_far_away=False)

    @value_accept(radius=(0.05, 0.95), min_peak_height=(0.1, 0.9), SID=(0, 180))
    def analyze(self, radius=0.5, min_peak_height=0.25, SID=100):
        """Analyze the starshot image.

        Analyze finds the minimum radius and center of a circle that touches all the lines
        (i.e. the wobble circle diameter and wobble center).

        Parameters
        ----------
        radius : float, optional
            Distance in % between starting point and closest image edge; used to build the circular profile which finds
            the radiation lines. Must be between 0.05 and 0.95.
        min_peak_height : float, optional
            The percentage minimum height a peak must be to be considered a valid peak. A lower value catches
            radiation peaks that vary in magnitude (e.g. different MU delivered), but could also pick up noise.
            Increase value for noisy images.
        SID : int, float, optional
            The source-to-image distance in cm. If a value != 100 is passed in, results will be scaled to 100cm. E.g. a wobble of
            3.0 pixels at an SID of 150cm will calculate to 2.0 pixels [3 / (150/100)].

        Raises
        ------
        AttributeError
            If an image has not yet been loaded.
        """
        # error checking
        if not self.image_is_loaded:
            raise AttributeError("Starshot image not yet loaded")

        # check inversion
        self._check_image_inversion()

        # set starting point automatically if not yet set
        if not self.start_point_is_set:
            self._auto_set_start_point()

        # set profile extraction radius
        self.circle_profile.radius = self._convert_radius_perc2pix(radius)
        # extract the circle profile
        self.circle_profile.get_profile(self.image.array)
        # find the radiation lines using the peaks of the profile
        self.lines = self.circle_profile.find_rad_lines(min_peak_height)
        # find the wobble
        self._find_wobble_2step(SID)

    def _convert_radius_perc2pix(self, radius):
        """Convert a percent radius to distance in pixels, based on the distance from center point to image
            edge.

        Parameters
        ----------
        radius : float
            The radius ratio (e.g. 0.5).
        """
        dist = self.image.dist2edge_min(self.circle_profile.center)
        return dist*radius

    @property
    def image_is_loaded(self):
        """Boolean property specifying if an image has been loaded."""
        if self.image.array.size == 0:
            return False
        else:
            return True

    @property
    def start_point_is_set(self):
        """Boolean specifying if a start point has been set."""
        if self.circle_profile.center.x == 0:
            return False
        else:
            return True

    def _scale_wobble(self, SID):
        """Scale the determined wobble by the SID.

        Parameters
        ----------
        SID : int, float
            Source to image distance in cm.
        """
        # convert wobble to mm if possible
        if self.image.dpmm != 0:
            self._tolerance_unit = 'mm'
            self.wobble.radius_mm = self.wobble.radius / self.image.dpmm
        else:
            self._tolerance_unit = 'pixels'
            self.wobble.radius_mm = self.wobble.radius

        self.wobble.radius /= SID / 100
        self.wobble.radius_mm /= SID / 100

    def _find_wobble_2step(self, SID):
        """Find the smallest radius ("wobble") and center of a circle that touches all the star lines.

        Notes
        -----
        Wobble determination is accomplished by two rounds of searching. The first round finds the radius and center down to
        the nearest pixel. The second round finds the center and radius down to sub-pixel precision using parameter scale.
        This methodology is faster than one round of searching at sub-pixel precision.

        See Also
        --------
        analyze : Further parameter info.
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
        """An iterative method that moves element by element to the point of minimum distance to all radiation lines.

        Parameters
        ----------
        tolerance : float
            The value differential between the outside elements and center element to stop the algorithm.
        start_point : geometry.Point
            The starting point for the algorithm.
        scale : int, float
            The scale of the search in pixels. E.g. 0.1 searches at 0.1 pixel precision.
        """
        # TODO: use an optimization function instead of evolutionary search
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
        """Boolean specifying whether the determined wobble was within tolerance."""
        if self.wobble.radius_mm * 2 < self._tolerance:
            return True
        else:
            return False

    def return_results(self):
        """Return the results of the analysis.

        Returns
        -------
        string
            A string with a statement of the minimum circle.
        """
        if self.passed:
            passfailstr = 'PASS'
        else:
            passfailstr = 'FAIL'

        string = ('\nResult: %s \n\n'
                  'The minimum circle that touches all the star lines has a diameter of %4.3g %s. \n\n'
                  'The center of the minimum circle is at %4.1f, %4.1f') % (passfailstr, self.wobble.radius_mm*2, self._tolerance_unit,
                                                                            self.wobble.center.x, self.wobble.center.y)
        return string

    def plot_analyzed_image(self, plot=None, show=True):
        """Draw the star lines, profile circle, and wobble circle on a matplotlib figure.

        Parameters
        ----------
        plot : matplotlib.image.AxesImage, optional
            The plot to draw on. If None, will create a new one.
        """
        # plot image
        if plot is None:
            imgplot = plt.imshow(self.image.array)
        else:
            plot.axes.imshow(self.image.array)
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
        if show:
            if plot is None:
                plt.show()
            else:
                plot.draw()
                plot.axes.hold(False)

    def run_demo(self, show=True):
        """Demonstrate the Starshot module using the demo image."""
        self.load_demo_image()
        self.analyze()
        print(self.return_results())
        self.plot_analyzed_image(show=show)


class Wobble(Circle):
    """A class that holds the wobble information of the Starshot analysis.

    Attributes
    ----------
    radius_mm : The radius of the Circle in **mm**.
    """
    def __init__(self, center_point=None, radius=None):
        super().__init__(center_point=center_point, radius=radius)
        self.radius_mm = 0  # The radius of the wobble in mm; as opposed to pixels.


class StarProfile(CircleProfile):
    """Class that holds and analyzes the circular profile which finds the radiation lines."""
    def __init__(self):
        super().__init__()

    def get_profile(self, image_array):
        """Take the profile over the image array. Overloads to also correct for profile positioning.

        See Also
        --------
        core.profile.CircleProfile.get_profile : Further parameter info
        """
        super().get_profile(image_array)
        self._roll_prof_to_midvalley()

    def _roll_prof_to_midvalley(self):
        """Roll the circle profile so that its edges are not near a radiation line.
            This is a prerequisite for properly finding star lines.
        """
        # Find index of the min value(s)
        min_idx = np.where(self.y_values == self.y_values.min())
        # Roll the profile and x and y coordinates
        self.y_values = np.roll(self.y_values, -min_idx[0][0])
        self.x_locs = np.roll(self.x_locs, -min_idx[0][0])
        self.y_locs = np.roll(self.y_locs, -min_idx[0][0])

    def find_rad_lines(self, min_peak_height, min_peak_distance=0.03):
        """Find and match the positions of peaks in the circle profile (radiation lines)
            and map their positions to the starshot image.

        Radiation lines are found by finding the FWHM of the radiation spokes, then matching them
        to form lines.

        Returns
        -------
        lines : list
            A list of Lines (radiation lines) found.

        See Also
        --------
        Starshot.analyze() : min_peak_height parameter info
        core.profile.CircleProfile.find_FWXM_peaks : min_peak_distance parameter info.
        geometry.Line : returning object
        """
        # find the FWHM-C
        self.find_FWXM_peaks(min_peak_height=min_peak_height, min_peak_distance=min_peak_distance)

        # Match the peaks to form lines
        lines = self.match_peaks()
        return lines

    def match_peaks(self):
        """Match the peaks found to the same radiation lines.

        Peaks are matched by connecting the existing peaks based on an offset of peaks. E.g. if there are
        12 peaks, there must be 6 radiation lines. Furthermore, assuming star lines go all the way across the CAX,
        the 7th peak will be the opposite peak of the 1st peak, forming a line. This method is robust to
        starting points far away from the real center.
        """
        num_rad_lines = int(len(self.peaks) / 2)
        offset = num_rad_lines
        line_list = [Line(self.peaks[line], self.peaks[line+offset]) for line in range(num_rad_lines)]
        return line_list

# ----------------------------
# Starshot demo
# ----------------------------
if __name__ == '__main__':
    Starshot().run_demo()
    # star = Starshot()
    # star.load_demo_image()
    # star.analyze(radius=0.9)
