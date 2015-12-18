"""The planar imaging module analyzes phantom images taken with the kV or MV imager in 2D; for example, the Leeds and PipsPro phantoms."""
import copy
from functools import lru_cache
import os.path as osp

import matplotlib.pyplot as plt
import numpy as np
from skimage import feature, measure

from pylinac.cbct import DiskROI, LowContrastDiskROI as LCDiskROI
from pylinac.core.geometry import Point
from pylinac.core.image import Image
from pylinac.core.io import get_url
from pylinac.core.profile import CollapsedCircleProfile


class PipsPro:
    """Class for analyzing high and low contrast of a PipsPro MV phantom."""

    def __init__(self, filepath):
        """
        Parameters
        ----------
        filepath : str
            Path to the image file.
        """
        self.image = Image.load(filepath)
        self.image.invert()
        self.image.ground()

    @classmethod
    def from_demo_image(cls):
        """Instantiate and load the demo image."""
        demo_file = osp.join(osp.dirname(__file__), 'demo_files', 'planar_imaging', 'pipspro.dcm')
        return cls(demo_file)

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

    @staticmethod
    def run_demo():
        """Run the PipsPro QC-3 phantom analysis demonstration."""
        pp = PipsPro.from_demo_image()
        pp.analyze()
        pp.plot_analyzed_image()

    @property
    @lru_cache()
    def _phan_region(self):
        """The skimage region of the phantom outline."""
        regions = _get_canny_regions(self.image)
        blobs = []
        for phantom_idx, region in enumerate(regions):
            if region.area < 50:
                continue
            semi_round = 0.7 > region.eccentricity > 0.3
            hollow = region.extent < 0.025
            angled = region.orientation > 0.2 or region.orientation < -0.2
            if semi_round and hollow and angled:
                blobs.append(phantom_idx)

        if len(blobs) < 1:
            raise ValueError("Unable to find the PipsPro phantom in the image.")

        # find the biggest ROI and call that the phantom outline
        big_roi_idx = np.argsort([regions[phan].major_axis_length for phan in blobs])[-1]
        phantom_idx = blobs[big_roi_idx]

        return regions[phantom_idx]

    @property
    def _phan_radius(self):
        """The radius of the phantom in pixels; the value itself doesn't matter, it's just
        used for relative distances to ROIs.

        Returns
        -------
        radius : float
        """
        return self._phan_region.major_axis_length / 14

    @property
    def _phan_angle(self):
        """The angle of the phantom.

        Returns
        -------
        angle : float
            The angle in degrees.
        """
        return -np.rad2deg(self._phan_region.orientation)

    @property
    def _phan_center(self):
        """The center point of the phantom.

        Returns
        -------
        center : Point
        """
        bbox = self._phan_region.bbox
        y = abs(bbox[0] - bbox[2]) / 2 + min(bbox[0], bbox[2])
        x = abs(bbox[1] - bbox[3]) / 2 + min(bbox[1], bbox[3])
        return Point(x, y)

    def _low_contrast(self):
        """Sample the detail contrast regions."""
        angles = np.array([90, -90, 55, -55, 128, -128]) + self._phan_angle
        dists = np.array([2, 2, 2.4, 2.4, 2.4, 2.4]) * self._phan_radius
        rrois = []

        # background ROI
        bg_roi = LowContrastDiskROI(self.image, angles[0], 0.5 * self._phan_radius, dists[0], self._phan_center, 0.05)

        for dist, angle in zip(dists[1:], angles[1:]):
            roi = LowContrastDiskROI(self.image, angle, 0.5 * self._phan_radius, dist, self._phan_center,
                                      0.05, background=bg_roi.pixel_value)
            rrois.append(roi)
        return bg_roi, rrois

    def _high_contrast(self):
        """Sample the high-contrast line pair regions."""
        dists = np.array([2.8, -2.8, 1.45, -1.45, 0]) * self._phan_radius
        rrois = []
        for dist in dists:
            roi = HighContrastDiskROI(self.image, self._phan_angle, 0.5*self._phan_radius, dist, self._phan_center,
                                      0.05)
            rrois.append(roi)
        return rrois

    def analyze(self, low_contrast_threshold=0.005, hi_contrast_threshold=0.5, invert=False):
        """Analyze the PipsPro phantom.

        Parameters
        ----------
        low_contrast_threshold : float
            The threshold for the low-contrast bubbles to be "seen".
        hi_contrast_threshold : float
            The threshold percentage that the relative MTF must be above to be "seen". Must be between 0 and 1.
        invert : bool
            Whether to force an inversion of the image. Pylinac tries to infer the correct inversion but uneven
            backgrounds can cause this analysis to fail. If the contrasts/MTF ROIs appear correctly located but the
            plots are wonky, try setting this to True.
        """
        self.image.check_inversion(box_size=30, offset=int(0.05 * max(self.image.shape)))
        if invert:
            self.image.invert()
        self.low_contrast_threshold = low_contrast_threshold
        self.hi_contrast_threshold = hi_contrast_threshold

        self._lcbgroi, self._lcrois = self._low_contrast()
        self._hcrois = self._high_contrast()

    def plot_analyzed_image(self, image=True, lowcontrast=True, highcontrast=True, show=True):
        """Plot the analyzed image."""
        num_plots = sum((image, lowcontrast, highcontrast))
        if num_plots < 1:
            return
        # set up axes and make axes iterable
        fig, axes = plt.subplots(1, num_plots)
        fig.subplots_adjust(wspace=0.4)
        if num_plots < 2:
            axes = (axes,)
        axes = iter(axes)

        # plot the marked image
        if image:
            img_ax = next(axes)
            self.image.plot(ax=img_ax, show=False)
            img_ax.axis('off')
            img_ax.set_title('PipsPro Phantom Analysis')

            # plot the low contrast ROIs
            self._lcbgroi.plot2axes(img_ax, edgecolor='b')
            for roi in self._lcrois:
                roi.plot2axes(img_ax, edgecolor=roi.plot_color)
            # plot the high-contrast ROIs
            for roi in self._hcrois:
                roi.plot2axes(img_ax, edgecolor=roi.plot_color)

        # plot the low contrast values
        if lowcontrast:
            lowcon_ax = next(axes)
            _plot_lowcontrast(lowcon_ax, self._lcrois, self.low_contrast_threshold)

        # plot the high contrast MTF
        if highcontrast:
            hicon_ax = next(axes)
            mtfs = [roi.mtf for roi in self._hcrois]
            mtfs /= mtfs[0]
            _plot_highcontrast(hicon_ax, mtfs, self.hi_contrast_threshold)

        if show:
            plt.show()

    def save_analyzed_image(self, filename, **kwargs):
        """Save the analyzed image to a file.

        Parameters
        ----------
        filename : str
            The location and filename to save to.
        kwargs
            Keyword arguments are passed to plt.savefig().
        """
        self.plot_analyzed_image(show=False)
        plt.savefig(filename, **kwargs)


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
        self.image.invert()
        self.image.ground()

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

    @property
    @lru_cache()
    def _blobs(self):
        """The indices of the regions that were significant; i.e. a phantom circle outline or lead/copper square."""
        blobs = []
        for idx, region in enumerate(self._regions):
            if region.area < 100:
                continue
            round = region.eccentricity < 0.25
            if round:
                blobs.append(idx)
        if not blobs:
            raise ValueError("Could not find the phantom in the image.")
        return blobs

    @property
    @lru_cache()
    def _regions(self):
        """All the regions of the canny image that were labeled."""
        regions = _get_canny_regions(self.image)
        return regions

    def _determine_phantom_center(self):
        """Determine the phantom center.

        This is done by searching for circular ROIs of the canny image. Those that are circular and roughly the
        same size as the biggest circle ROI are all sampled for the center of the bounding box. The values are
        averaged over all the detected circles to give a more robust value.

        Returns
        -------
        center : Point
        """
        circles = [roi for roi in self._blobs if
                   np.isclose(self._regions[roi].major_axis_length, self._determine_radius() * 3.35, rtol=0.3)]

        # get average center of all circles
        bboxs = [self._regions[roi].bbox for roi in circles]
        y = np.mean([abs(bbox[0] - bbox[2]) / 2 + min(bbox[0], bbox[2]) for bbox in bboxs])
        x = np.mean([abs(bbox[1] - bbox[3]) / 2 + min(bbox[1], bbox[3]) for bbox in bboxs])
        return Point(x, y)

    def _determine_phantom_angle(self, center):
        """Determine the angle of the phantom.

        This is done by searching for square-like boxes of the canny image. There are usually two: one lead and
        one copper. The box with the highest intensity (lead) is identified. The angle from the center of the lead
        square bounding box and the phantom center determines the phantom angle.

        Parameters
        ----------
        center : Point
            The center point of the phantom

        Returns
        -------
        angle : float
            The angle in radians.
        """
        expected_length = self._determine_radius() * 0.52
        square_rois = [roi for roi in self._blobs if np.isclose(self._regions[roi].major_axis_length, expected_length, rtol=0.2)]
        if not square_rois:
            raise ValueError("Could not find the angle of the image.")
        regions = self._regions
        lead_idx = np.argsort([regions[roi].mean_intensity for roi in square_rois])[-1]
        lead_roi = regions[square_rois[lead_idx]]
        y = abs(lead_roi.bbox[0] - lead_roi.bbox[2]) / 2 + min(lead_roi.bbox[0], lead_roi.bbox[2])
        x = abs(lead_roi.bbox[1] - lead_roi.bbox[3]) / 2 + min(lead_roi.bbox[1], lead_roi.bbox[3])
        lead_center = Point(x, y)

        adjacent = lead_center.x - center.x
        opposite = lead_center.y - center.y
        angle = np.arctan2(opposite, adjacent)
        return angle

    @lru_cache()
    def _determine_radius(self):
        """Determine the radius of the phantom.

        The radius is determined by finding the largest of the detected blobs of the canny image and taking
        its major axis length.

        Returns
        -------
        radius : float
            The radius of the phantom in pixels. The actual value is not important; it is used for scaling the
            distances to the low and high contrast ROIs.
        """
        big_circle_idx = np.argsort([self._regions[roi].major_axis_length for roi in self._blobs])[-1]
        circle_roi = self._regions[self._blobs[big_circle_idx]]
        radius = circle_roi.major_axis_length / 3.35
        return radius

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
            :class:`~pylinac.planar_imaging.LowContrastDistROI` instances of the contrast ROIs.
        reference ROIs : list
            :class:`~pylinac.planar_imaging.LowContrastDistROI` instances of the reference ROIs;
            pixel values of the reference ROIs determines the background for the contrast ROIs.
        """
        angle = np.degrees(angle)
        bubble_angles1 = np.linspace(30, 149, num=9)
        bubble_angles2 = np.linspace(209, 331, num=9)
        bubble_angles = np.concatenate((bubble_angles1, bubble_angles2))
        bubble_radius = 0.025 * radius

        # sample the contrast ROIs
        bubble_dist = 0.785 * radius
        crois = []
        for angle_delta in bubble_angles:
            roi = LowContrastDiskROI(self.image, angle - angle_delta, bubble_radius, bubble_dist, center, self.low_contrast_threshold)
            crois.append(roi)

        # sample the reference ROIs
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
            :class:`~pylinac.planar_imaging.HighContrastDiskROI` instances of the line pairs.
        reference ROIs : list
            :class:`~pylinac.planar_imaging.HighContrastDiskROI` instances of the solid ROIs that
            determine the normalization value for MTF.
        """
        angle = np.degrees(angle)

        # sample ROIs of the reference areas
        ref_angles = [303, 271]
        ref_dists = [0.3 * radius, 0.25 * radius]
        ref_radius = 0.04 * radius
        rrois = []
        for nominal_angle, dist in zip(ref_angles, ref_dists):
            roi = HighContrastDiskROI(self.image, angle - nominal_angle, ref_radius, dist, center,
                                      self.hi_contrast_threshold)
            rrois.append(roi)
        mtf_norm_val = (rrois[0].pixel_value - rrois[1].pixel_value) / (rrois[0].pixel_value + rrois[1].pixel_value)

        # sample ROIs of each line pair region
        # ordering goes from the "biggest" line pair region downward
        contrast_angles = [-144.8, -115.1, -62.5, -169.7, -153.4, -25, 169.7, 151.6, 27]
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
        demo_file = osp.join(osp.dirname(__file__), 'demo_files', 'planar_imaging', 'leeds.dcm')
        return cls(demo_file)

    @staticmethod
    def run_demo():
        """Run the Leeds TOR phantom analysis demonstration."""
        leeds = LeedsTOR.from_demo_image()
        leeds.analyze()
        leeds.plot_analyzed_image()

    def analyze(self, low_contrast_threshold=0.005, hi_contrast_threshold=0.4, invert=False):
        """Analyze the image.

        Parameters
        ----------
        low_contrast_threshold : float
            The threshold for the low-contrast bubbles to be "seen".
        hi_contrast_threshold : float
            The threshold percentage that the relative MTF must be above to be "seen". Must be between 0 and 1.
        invert : bool
            Whether to force an inversion of the image. Pylinac tries to infer the correct inversion but uneven
            backgrounds can cause this analysis to fail. If the contrasts/MTF ROIs appear correctly located but the
            plots are wonky, try setting this to True.
        """
        self.image.check_inversion(box_size=30, offset=int(0.05*max(self.image.shape)))
        if invert:
            self.image.invert()
        self.low_contrast_threshold = low_contrast_threshold
        self.hi_contrast_threshold = hi_contrast_threshold

        radius = self._determine_radius()
        center = self._determine_phantom_center()
        angle = self._determine_phantom_angle(center)
        if not self._is_clockwise(center, radius, angle):
            center, angle = self._flip_image_data(center, angle)
        self.lcrois, self.lcrrois = self._low_contrast(radius, center, angle)
        self.hcrois, self.hcrrois = self._high_contrast(radius, angle, center)

    def _flip_image_data(self, center, angle):
        """Flip the image left->right and invert the center, and angle as appropriate."""
        self.image.array = np.fliplr(self.image.array)
        new_x = self.image.shape[1] - center.x
        new_center = Point(new_x, center.y)
        new_angle = np.pi - angle
        return new_center, new_angle

    def plot_analyzed_image(self, image=True, low_contrast=True, high_contrast=True, show=True):
        """Plot the analyzed image, which includes the original image with ROIs marked, low-contrast plots
        and high-contrast plots.

        Parameters
        ----------
        image : bool
            Show the image.
        low_contrast : bool
            Show the low contrast values plot.
        high_contrast : bool
            Show the high contrast values plot.
        show : bool
            Whether to actually show the image when called.
        """
        num_plots = sum((image, low_contrast, high_contrast))
        fig, axes = plt.subplots(1, num_plots)
        fig.subplots_adjust(wspace=0.4)
        if num_plots < 2:
            axes = (axes,)
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
        if low_contrast:
            lowcon_ax = next(axes)
            _plot_lowcontrast(lowcon_ax, self.lcrois, self.low_contrast_threshold)

        # plot the high contrast MTF
        if high_contrast:
            hicon_ax = next(axes)
            hc_rois = [roi.mtf for roi in self.hcrois]
            hc_rois.insert(0, 1)
            _plot_highcontrast(hicon_ax, hc_rois, self.hi_contrast_threshold)

        if show:
            plt.show()

    def save_analyzed_image(self, filename, **kwargs):
        """Save the analyzed image to a file.

        Parameters
        ----------
        filename : str
            The location and filename to save to.
        kwargs
            Keyword arguments are passed to plt.savefig().
        """
        self.plot_analyzed_image(show=False)
        plt.savefig(filename, **kwargs)


class LowContrastDiskROI(LCDiskROI):
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


def _get_canny_regions(image, sigma=2, percentiles=(0.001, 0.01)):
    img_copy = copy.copy(image)
    img_copy.filter(kind='gaussian', size=sigma)
    img_copy.ground()
    lo_th, hi_th = np.percentile(img_copy, percentiles)
    c = feature.canny(img_copy, low_threshold=lo_th, high_threshold=hi_th)
    labeled = measure.label(c)
    regions = measure.regionprops(labeled, intensity_image=img_copy)
    return regions


def _plot_lowcontrast(axes, rois, threshold):
    line1, = axes.plot([roi.contrast for roi in rois], marker='o', color='m', label='Contrast')
    axes.axhline(threshold, color='k')
    axes.grid('on')
    axes.set_title('Low-frequency Contrast')
    axes.set_xlabel('ROI #')
    axes.set_ylabel('Contrast')
    axes2 = axes.twinx()
    line2, = axes2.plot([roi.contrast_to_noise for roi in rois], marker='^', label='CNR')
    axes.legend(handles=[line1, line2])


def _plot_highcontrast(axes, rois, threshold):
    axes.plot(rois, marker='*')
    axes.axhline(threshold, color='k')
    axes.grid('on')
    axes.set_title('High-frequency rMTF')
    axes.set_xlabel('Line pair region #')
    axes.set_ylabel('relative MTF')
