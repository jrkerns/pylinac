"""The kV module analyzes phantom images taken with the kV imager; for example, the Leeds phantom."""
from functools import lru_cache
import os.path as osp

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

from pylinac.cbct import DiskROI, LowContrastDiskROI
from pylinac.core.geometry import Point
from pylinac.core.image import Image
from pylinac.core.io import get_url
from pylinac.core.mask import bounding_box, filled_area_ratio, square_ratio
from pylinac.core.profile import CollapsedCircleProfile


class LeedsTOR:
    """Class that analyzes Leeds TOR phantom planar kV images for kV QA."""

    def __init__(self, filepath):
        """
        Parameters
        ----------
        filepath : str
            Path to the local file.
        """
        self.image = Image.load(filepath)

    @classmethod
    def from_url(cls, url):
        """
        Parameters
        ----------
        url : str
            The URL to the image.
        """
        image_file = get_url(url)
        return cls(image_file)

    def _find_phantom(self):
        """Find the phantom center and radius within the image. This searches for the thin metal band around
        the edge of the Leeds phantom.

        Returns
        -------
        center : Point
            The center point of the phantom
        radius : float
            The radius of the phantom to the metal ring
        labeled array : ndarray
            A numpy array with labeled ROIs
        """
        found = False
        thres_pct = 40

        # threshold image iteratively until the ring structure is found
        while not found:
            # threshold and label the binary image
            threshold = np.percentile(self.image, thres_pct)
            bin_img = self.image.as_binary(threshold)
            labeled_arr, num_roi = ndimage.label(bin_img)
            roi_sizes, bin_edges = np.histogram(labeled_arr, bins=num_roi + 1)
            # search through the ROIs for a ring-like structure
            for roi in range(num_roi):
                if roi_sizes[roi] < 200:
                    continue
                roi_img = np.where(labeled_arr == roi, 1, 0)
                ymin, ymax, xmin, xmax = bounding_box(roi_img)
                fill_area = filled_area_ratio(roi_img)
                is_squarelike = 0.97 < square_ratio(roi_img) < 1.03
                is_hollow = fill_area < 0.07
                is_mediumsized = (ymax - ymin) * (xmax - xmin) > 0.3 * self.image.shape[0] * self.image.shape[1]

                if is_squarelike and is_hollow and is_mediumsized:
                    found = True
                    break

            thres_pct += 7

        # determine center of ROI, which is the center of phantom
        x_cen = (xmax - xmin)/2 + xmin
        y_cen = (ymax - ymin)/2 + ymin
        radius = (xmax - xmin)/2
        return Point(x_cen, y_cen), radius, labeled_arr

    def _find_angle(self, labeled_array, radius, center):
        """Find the angle of the phantom from the center to the lead square.

        Returns
        -------
        angle : float
            The angle of the phantom in radians.
        """
        expected_size = (radius ** 2) * 0.116
        roi_sizes, bin_edges = np.histogram(labeled_array, bins=labeled_array.max() + 1)

        # search for the square-like ROI
        found = False
        roi = 0
        while not found:
            if roi >= len(roi_sizes):
                raise ValueError("Could not find the Leeds phantom within the image.")
            if roi_sizes[roi] < 100:
                roi += 1
                continue
            roi_img = np.where(labeled_array == roi, 1, 0)
            fill_area = filled_area_ratio(roi_img)
            is_squarelike = 0.98 < square_ratio(roi_img) < 1.02
            is_solid = fill_area > 0.5
            is_right_size = expected_size * 1.04 > roi_sizes[roi] > expected_size * 0.96

            if is_squarelike and is_solid and is_right_size:
                found = True
            else:
                roi += 1

        # now that the ROI is found, determine the angle from the center
        com = ndimage.measurements.center_of_mass(np.where(labeled_array == roi, 1, 0))
        adjacent = com[1] - center.x
        opposite = com[0] - center.y
        angle = np.arctan2(opposite, adjacent)
        return angle

    def _is_clockwise(self, center, radius, angle):
        """Determine if the low-contrast bubbles go from high to low clockwise or counter-clockwise.

        Returns
        -------
        boolean
        """
        circle = CollapsedCircleProfile(center, radius * 0.8, self.image, angle, width_ratio=0.05)
        first_set = circle.find_peaks(search_region=(0, 0.5), threshold=0.1)
        second_set = circle.find_peaks(search_region=(0.5, 1.0), threshold=0.1)
        return len(first_set) > len(second_set)

    def _low_contrast(self, radius, center, angle):
        """Perform the low-contrast analysis. This samples the bubbles and a background bubble just beneath it to
        determine contrast and contrast-to-noise.

        Returns
        -------
        contrast ROIs : list
            LeedsLowContrastDistROI instances
        reference ROIs : list
            Reference ROIs that determine the background pixel values.
        """
        angle = np.degrees(angle)
        # bubble_angles = [60, 45, 30, 15, 0, -15, -30, -45, -60, -120, -135, -150, -165, -180, 165, 150, 135, 120]
        bubble_angles = list(range(30, 151, 15))
        bubble_angles += list(range(210, 331, 15))
        bubble_radius = 0.03 * radius

        # sample contrast ROIs
        bubble_dist = 0.8 * radius
        crois = []
        for angle_delta in bubble_angles:
            roi = LeedsLowContrastDistROI(self.image, angle - angle_delta , bubble_radius, bubble_dist, center, self.low_contrast_threshold)
            crois.append(roi)

        # sample reference ROIs
        bubble_dist = 0.65 * radius
        rrois = []
        for idx, angle_delta in enumerate(bubble_angles):
            roi = DiskROI(self.image, angle - angle_delta, bubble_radius, bubble_dist, center)
            crois[idx].background = roi.pixel_value
            rrois.append(roi)

        return crois, rrois

    def _high_contrast(self, radius, angle, center):
        """Perform high-contrast analysis. This samples disks within the line-pair region and calculates
        relative MTF from the min and max values.

        Returns
        -------
        contrast ROIs : list
            HighContrastDiskROI
        reference ROIs : list
            Reference ROIs that determine the background pixel values.
        """
        angle = np.degrees(angle)

        # sample ROIs of the reference areas
        ref_angles = [303, 270]
        ref_dists = [0.3 * radius, 0.25 * radius]
        ref_radius = 0.04 * radius
        rrois = []
        for nominal_angle, dist in zip(ref_angles, ref_dists):
            roi = HighContrastDiskROI(self.image, angle - nominal_angle, ref_radius, dist, center,
                                      self.hi_contrast_threshold)
            rrois.append(roi)
        mtf_norm_val = (rrois[0].pixel_value - rrois[1].pixel_value) / (rrois[0].pixel_value + rrois[1].pixel_value)

        # sample ROIs of each line pair region
        contrast_angles = [-144.8, -115.1, -62.5, -169.7, -153.4, -23.7, 169.7, 151.6, 30.2]
        contrast_dists = np.array([0.3, 0.187, 0.187, 0.252, 0.092, 0.094, 0.252, 0.094, 0.0958]) * radius
        contrast_radii = np.array([0.04, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.018, 0.018, 0.015, 0.015, 0.012]) * radius
        crois = []
        for nominal_angle, dist, cradius in zip(contrast_angles, contrast_dists, contrast_radii):
            roi = HighContrastDiskROI(self.image, angle + nominal_angle + 90, cradius, dist, center, self.hi_contrast_threshold, mtf_norm=mtf_norm_val)
            crois.append(roi)

        return crois, rrois

    @classmethod
    def from_demo_image(cls):
        """Instantiate and load the demo image."""
        demo_file = osp.join(osp.dirname(__file__), 'demo_files', 'kv', 'leeds.dcm')
        return cls(demo_file)

    def analyze(self, low_contrast_threshold=0.005, hi_contrast_threshold=0.4):
        """Analyze the image.

        Parameters
        ----------
        low_contrast_threshold : float
            The threshold for the low-contrast bubbles to be "seen".
        hi_contrast_threshold : float
            The threshold percentage that the relative MTF must be above to be "seen". Must be between 0 and 1.
        """
        self.low_contrast_threshold = low_contrast_threshold
        self.hi_contrast_threshold = hi_contrast_threshold
        center, radius, binary_img = self._find_phantom()
        angle = self._find_angle(binary_img, radius, center)
        if not self._is_clockwise(center, radius, angle):
            self.image.array = np.fliplr(self.image.array)
            center, radius, binary_img = self._find_phantom()
            angle = self._find_angle(binary_img, radius, center)
        self.lcrois, self.lcrrois = self._low_contrast(radius, center, angle)
        self.hcrois, self.hcrrois = self._high_contrast(radius, angle, center)

    def plot_analyzed_image(self, image=True, lowcontrast=True, highcontrast=True, show=True):
        """Plot the analyzed image, which includes the original image with ROIs marked, low-contrast plots
        and high-contrast plots.

        Parameters
        ----------
        show : boolean
            Whether to actually show the image when called.
        """
        num_plots = sum((image, lowcontrast, highcontrast))
        fig, axes = plt.subplots(1, num_plots)
        fig.subplots_adjust(wspace=0.4)
        axes = iter(axes)

        if image:
            img_ax = next(axes)
            self.image.plot(ax=img_ax, show=False)
            img_ax.axis('off')
            img_ax.set_title('Leeds TOR Phantom Analysis')

            # plot the low contrast ROIs
            for roi in self.lcrois:
                roi.plot2axes(img_ax, edgecolor=roi.plot_color)
            for roi in self.lcrrois:
                roi.plot2axes(img_ax, edgecolor='g')
            # plot the high-contrast ROIs
            for roi in self.hcrois:
                roi.plot2axes(img_ax, edgecolor=roi.plot_color)
            for roi in self.hcrrois:
                roi.plot2axes(img_ax, edgecolor='g')

        # plot the low contrast values
        if lowcontrast:
            lowcon_ax = next(axes)
            line1, = lowcon_ax.plot([roi.contrast for roi in self.lcrois], marker='o', color='m', label='Contrast')
            lowcon_ax.axhline(self.low_contrast_threshold, color='k')
            lowcon_ax.grid('on')
            lowcon_ax.set_title('Low-frequency Contrast')
            lowcon_ax.set_xlabel('ROI #')
            lowcon_ax.set_ylabel('Contrast')
            lowcon_ax2 = lowcon_ax.twinx()
            line2, = lowcon_ax2.plot([roi.contrast_to_noise for roi in self.lcrois], marker='^', label='CNR')
            plt.legend(handles=[line1, line2])

        # plot the high contrast MTF
        if highcontrast:
            hicon_ax = next(axes)
            hc_rois = [roi.mtf for roi in self.hcrois]
            hc_rois.insert(0, 1)
            hicon_ax.plot(hc_rois, marker='*')
            hicon_ax.axhline(self.hi_contrast_threshold, color='k')
            hicon_ax.grid('on')
            hicon_ax.set_title('High-frequency rMTF')
            hicon_ax.set_xlabel('Line pair region #')
            hicon_ax.set_ylabel('relative MTF')

        if show:
            plt.show()

    def save_analyzed_image(self, filename, **kwargs):
        """Save the analyzed image to a file."""
        self.plot_analyzed_image(show=False)
        plt.savefig(filename, **kwargs)


class LeedsLowContrastDistROI(LowContrastDiskROI):
    """A low-contrast ROI class that uses the actual contrast value for pass/fail status rather than the contrast constant."""

    @property
    def passed(self):
        return self.contrast > self.contrast_threshold


class HighContrastDiskROI(DiskROI):
    """A class for analyzing the high-contrast disks."""

    def __init__(self, array, angle, roi_radius, dist_from_center, phantom_center, contrast_threshold, mtf_norm=None):
        """
        Parameters
        ----------
        contrast_threshold : float, int
            The threshold for considering a bubble to be "seen".
        """
        super().__init__(array, angle, roi_radius, dist_from_center, phantom_center)
        self.contrast_threshold = contrast_threshold
        self.mtf_norm = mtf_norm

    @property
    def mtf(self):
        """The contrast of the bubble compared to background: (ROI - backg) / (ROI + backg)."""
        mtf = (self.max - self.min) / (self.max + self.min)
        if self.mtf_norm is not None:
            mtf /= self.mtf_norm
        return mtf

    @property
    def passed(self):
        """Boolean specifying if ROI pixel value was within tolerance of the nominal value."""
        return self.mtf > self.contrast_threshold

    @property
    def plot_color(self):
        """Return one of two colors depending on if ROI passed."""
        return 'blue' if self.passed else 'red'

    @property
    @lru_cache()
    def max(self):
        """The max pixel value of the ROI."""
        masked_img = self._get_roi_mask()
        return np.nanmax(masked_img)

    @property
    @lru_cache()
    def min(self):
        """The min pixel value of the ROI."""
        masked_img = self._get_roi_mask()
        return np.nanmin(masked_img)
