# -*- coding: utf-8 -*-
"""
The Starshot module analyses a starshot film or multiple superimposed EPID images that measures the wobble of the
radiation spokes, whether gantry, collimator, or couch. It is based on ideas from
`Depuydt et al <http://iopscience.iop.org/0031-9155/57/10/2997>`_
and `Gonzalez et al <http://dx.doi.org/10.1118/1.1755491>`_ and evolutionary optimization.
"""

import os.path as osp
import copy

import numpy as np
import matplotlib.pyplot as plt

from pylinac.core.decorators import value_accept
from pylinac.core.geometry import Point, Line, Circle
from pylinac.core.image import Image
from pylinac.core.io import get_filepath_UI, get_filenames_UI
from pylinac.core.profile import CircleProfile, SingleProfile


class Starshot:
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
        # demo_file = osp.join(demo_folder, 'DHMC_starshot.dcm')
        self.load_image(demo_file)

    def load_image(self, filepath):
        """Load the image via the file path.

        Parameters
        ----------
        filepath : str
            Path to the file to be loaded.
        """
        self.image = Image(filepath)
        # apply filter if it's a large image to reduce noise
        if self.image.shape[0] > 1100:
            self.image.median_filter(0.002)

    def load_multiple_images(self, filepath_list):
        """Load multiple images via the file path.

        .. versionadded:: 0.5.1

        Parameters
        ----------
        filepath_list : sequence
            An iterable sequence of filepath locations.
        """
        self.image = Image.combine_multiples(filepath_list)

    def load_multiple_images_UI(self):
        """Load multiple images via a dialog box.

        .. versionadded:: 0.5.1
        """
        path_list = get_filenames_UI()
        if path_list:
            self.load_multiple_images(path_list)

    def load_image_UI(self):
        """Load the image by using a UI dialog box."""
        path = get_filepath_UI()
        if path:
            self.load_image(path)

    @property
    def start_point(self):
        """The start point of the wobble search algorithm.

        After analysis this point is the wobble center.
        """
        return self.circle_profile.center

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
        central_array = self.image.array[top_third:bottom_third, left_third:right_third]

        x_sum = np.sum(central_array, 0)
        y_sum = np.sum(central_array, 1)

        # Calculate Full-Width, 80% Maximum
        fwxm_x_point = SingleProfile(x_sum).get_FWXM_center(80) + left_third
        fwxm_y_point = SingleProfile(y_sum).get_FWXM_center(80) + top_third

        # find maximum points
        x_max = np.unravel_index(np.argmax(central_array), central_array.shape)[1]  + left_third
        y_max = np.unravel_index(np.argmax(central_array), central_array.shape)[0] + top_third

        # which one is closer to the center
        fwxm_dist = Point(fwxm_x_point, fwxm_y_point).dist_to(self.image.center)
        max_dist = Point(x_max, y_max).dist_to(self.image.center)

        if fwxm_dist < max_dist:
            center_point = Point(fwxm_x_point, fwxm_y_point)
        else:
            center_point = Point(x_max, y_max)

        self.circle_profile.center = center_point

    @value_accept(radius=(0.2, 0.95), min_peak_height=(0.1, 0.9), SID=(40, 400))
    def analyze(self, radius=0.95, min_peak_height=0.25, SID=100, fwhm=True, recursive=True):
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

            .. note::
                For EPID images (e.g. superimposed collimator shots), the SID is in the DICOM file, and this
                value will always be used if it can be found, otherwise the passed value will be used.

        fwhm : bool
            If True (defualt), the center of the FWHM of the spokes will be determined.
            If False, the peak value location is used as the spoke center.
            .. note:: In practice, this ends up being a very small difference. Set to false if behavior is unexpected.
        recursive : bool
            If True (default), will recursively search for a "reasonable" wobble, meaning the wobble radius is
            <5mm, and the wobble location is somewhere near the starting point. If the wobble found was unreasonable,
            the minimum peak height is lowered incrementally. If at that point the wobble is still unreasonable, the
            radius is lowered (closer to center) and the search (including minimum peak height)
            If False, will not

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

        wobble_unreasonable = True
        orig_peak_height = copy.copy(min_peak_height)
        while wobble_unreasonable:
            # set profile extraction radius
            self.circle_profile.radius = self._convert_radius_perc2pix(radius)
            # extract the circle profile
            self.circle_profile.get_median_profile(self.image.array)
            # find the radiation lines using the peaks of the profile
            self.lines = self.circle_profile.find_rad_lines(min_peak_height, fwhm=fwhm)
            # find the wobble
            self._find_wobble_2step(SID)
            if not recursive:
                wobble_unreasonable = False
            else:
                if self.wobble.radius_mm < 5 and self.wobble.center.dist_to(self.start_point) < 30:
                    wobble_unreasonable = False
                else:
                    if min_peak_height > 0.15:
                        min_peak_height -= 0.07
                    elif radius > 0.3:
                        min_peak_height = orig_peak_height
                        radius -= 0.05
                    else:
                        raise RuntimeError("The algorithm was unable to determine a reasonable wobble. Try setting"
                                           "recursive to False")

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
        try:
            self.image.size
            return True
        except AttributeError:
            return False

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
        if self.image.dpmm is not None:
            self._tolerance_unit = 'mm'
            self.wobble.radius_mm = self.wobble.radius / self.image.dpmm
        else:
            self._tolerance_unit = 'pixels'
            self.wobble.radius_mm = self.wobble.radius

        if self.image.SID is not None:
            self.wobble.radius /= self.image.SID / 100
            self.wobble.radius_mm /= self.image.SID / 100
        else:
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
        sp = copy.copy(self.circle_profile.center)

        # first round of searching; this finds the circle to the nearest pixel
        normal_tolerance, normal_scale = 0.05, 1.0
        self._find_wobble(normal_tolerance, sp, normal_scale)

        # second round of searching; this finds the circle down to sub-pixel precision
        small_tolerance, small_scale = 0.0001, 100.0
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
        show : bool
            Whether to actually show the image.
        """
        plt.clf()
        imgplot = plt.imshow(self.image.array, cmap=plt.cm.Greys)

        # plot radiation lines
        for line in self.lines:
            line.add_to_axes(imgplot.axes, color='blue')

        # plot wobble circle
        self.wobble.add_to_axes(imgplot.axes, edgecolor='green')

        # plot profile circle
        self.circle_profile.add_to_axes(imgplot.axes, edgecolor='green')

        # tighten plot around image
        imgplot.axes.autoscale(tight=True)

        imgplot.axes.axis('off')

        # Finally, show it all
        if show:
            plt.show()

    def save_analyzed_image(self, filename, **kwargs):
        """Save the analyzed image plot to a file.

        Parameters
        ----------
        filename : str, IO stream
            The filename to save as. Format is deduced from string extention, if there is one. E.g. 'mystar.png' will
            produce a PNG image.

        kwargs
            All other kwargs are passed to plt.savefig().
        """
        self.plot_analyzed_image(show=False)

        plt.savefig(filename, **kwargs)

    def run_demo(self):
        """Demonstrate the Starshot module using the demo image."""
        self.load_demo_image()
        self.analyze()
        print(self.return_results())
        self.plot_analyzed_image()


class Wobble(Circle):
    """A class that holds the wobble information of the Starshot analysis.

    Attributes
    ----------
    radius_mm : The radius of the Circle in **mm**.
    """
    def __init__(self, center_point=None, radius=None):
        super().__init__(center_point=center_point, radius=radius)
        self.radius_mm = 0  # The radius of the wobble in mm; as opposed to pixels.

    @property
    def diameter_mm(self):
        return self.radius_mm*2


class StarProfile(CircleProfile):
    """Class that holds and analyzes the circular profile which finds the radiation lines."""
    def __init__(self, center=None, radius=None):
        super().__init__(center=center, radius=radius)

    def get_median_profile(self, image_array):
        """Take the profile over the image array. Overloads to also correct for profile positioning.

        See Also
        --------
        core.profile.CircleProfile.get_profile : Further parameter info
        """
        prof_size = 4*self.radius*np.pi
        super().get_profile(image_array, prof_size)

        mean_prof = np.zeros(prof_size)
        rrange = np.linspace(start=0.9, stop=1, num=int(self.radius*0.1))
        for rad in rrange[::-1]:
            prof = CircleProfile(self.center, self.radius*rad)
            prof.get_profile(image_array, prof_size)
            mean_prof += prof.y_values
        mean_prof /= len(rrange)

        self.y_values = mean_prof
        self._roll_prof_to_midvalley()
        self.ground()

    def _roll_prof_to_midvalley(self):
        """Roll the circle profile so that its edges are not near a radiation line.
            This is a prerequisite for properly finding star lines.
        """
        # vals, indx = self.find_peaks(return_it=True)
        # dindx = np.diff(indx)
        # roll_amount = int(indx[0] - np.median(dindx)/2)
        # Find index of the min value(s)
        roll_amount = np.where(self.y_values == self.y_values.min())[0][0]
        # Roll the profile and x and y coordinates
        self.y_values = np.roll(self.y_values, -roll_amount)
        self.x_locs = np.roll(self.x_locs, -roll_amount)
        self.y_locs = np.roll(self.y_locs, -roll_amount)
        return roll_amount

    def find_rad_lines(self, min_peak_height, min_peak_distance=0.02, fwhm=True):
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
        if fwhm:
            # find the FWHM-C
            self.find_FWXM_peaks(fwxm=70, min_peak_height=min_peak_height, min_peak_distance=min_peak_distance)
        else:
            self.find_peaks(min_peak_height, min_peak_distance)

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
    pass
    # Starshot().run_demo()
    star = Starshot()
    # star.load_image_UI()
    star.load_demo_image()
    star.analyze(radius=0.95, min_peak_height=0.25, fwhm=True)
    # star.analyze(recursive=True)
    # print(star.return_results())
    star.plot_analyzed_image()
    # star.save_analyzed_image('tester.png', bbox_inches='tight', pad_inches=0)
