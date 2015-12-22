"""The planar imaging module analyzes phantom images taken with the kV or MV imager in 2D; for example, the Leeds and PipsPro phantoms."""
import copy
from functools import lru_cache
import os.path as osp

import matplotlib.pyplot as plt
import numpy as np
from skimage import feature, measure

from .core.roi import LowContrastDiskROI, HighContrastDiskROI, DiskROI, bbox_center
from .core.geometry import Point
from .core.image import Image
from .core.io import get_url
from .core.profile import CollapsedCircleProfile


class ImagePhantomBase:
    """Base class for planar phantom classes."""
    _demo_filename = ''

    def __init__(self, filepath):
        """
        Parameters
        ----------
        filepath : str
            Path to the image file.
        """
        self.image = Image.load(filepath)
        self.image.invert()

    @classmethod
    def from_demo_image(cls):
        """Instantiate and load the demo image."""
        demo_file = osp.join(osp.dirname(__file__), 'demo_files', 'planar_imaging', cls._demo_filename)
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

    def _get_canny_regions(self, sigma=2, percentiles=(0.001, 0.01)):
        """Compute the canny edges of the image and return the connected regions found."""
        # copy, filter, and ground the image
        img_copy = copy.copy(self.image)
        img_copy.filter(kind='gaussian', size=sigma)
        img_copy.ground()

        # computer the canny edges with very low thresholds (detects nearly everything)
        lo_th, hi_th = np.percentile(img_copy, percentiles)
        c = feature.canny(img_copy, low_threshold=lo_th, high_threshold=hi_th)

        # label the canny edge regions
        labeled = measure.label(c)
        regions = measure.regionprops(labeled, intensity_image=img_copy)
        return regions

    def _plot_lowcontrast(self, axes, rois, threshold):
        """Plot the low contrast ROIs to an axes."""
        line1, = axes.plot([roi.contrast for roi in rois], marker='o', color='m', label='Contrast')
        axes.axhline(threshold, color='k')
        axes.grid('on')
        axes.set_title('Low-frequency Contrast')
        axes.set_xlabel('ROI #')
        axes.set_ylabel('Contrast')
        axes2 = axes.twinx()
        line2, = axes2.plot([roi.contrast_to_noise for roi in rois], marker='^', label='CNR')
        axes.legend(handles=[line1, line2])

    def _plot_highcontrast(self, axes, rois, threshold):
        """Plot the high contrast ROIs to an axes."""
        axes.plot(rois, marker='*')
        axes.axhline(threshold, color='k')
        axes.grid('on')
        axes.set_title('High-frequency rMTF')
        axes.set_xlabel('Line pair region #')
        axes.set_ylabel('relative MTF')


class PipsProQC3(ImagePhantomBase):
    """Class for analyzing high and low contrast of the PipsPro QC-3 MV phantom.

    Attributes
    ----------
    lc_rois : list
        :class:`~pylinac.core.roi.LowContrastDiskROI` instances of the low
        contrast ROIs, other than the reference ROI (below).
    lc_ref_rois : list
        :class:`~pylinac.core.roi.LowContrastDiskROI` instance of the low
        contrast reference ROI (15mm PVC).
    hc_rois : list
        :class:`~pylinac.core.roi.HighContrastDiskROI` instances of the
        high contrast line pair regions.
    """
    _demo_filename = 'pipspro.dcm'

    @staticmethod
    def run_demo():
        """Run the PipsPro QC-3 phantom analysis demonstration."""
        pp = PipsProQC3.from_demo_image()
        pp.analyze()
        pp.plot_analyzed_image()

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

        self.lc_ref_rois, self.lc_rois = self._low_contrast()
        self.hc_rois = self._high_contrast()

    def plot_analyzed_image(self, image=True, low_contrast=True, high_contrast=True, show=True):
        """Plot the analyzed image.

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
            self.lc_ref_rois.plot2axes(img_ax, edgecolor='b')
            for roi in self.lc_rois:
                roi.plot2axes(img_ax, edgecolor=roi.plot_color_constant)
            # plot the high-contrast ROIs
            for roi in self.hc_rois:
                roi.plot2axes(img_ax, edgecolor='b')

        # plot the low contrast values
        if low_contrast:
            lowcon_ax = next(axes)
            self._plot_lowcontrast(lowcon_ax, self.lc_rois, self.low_contrast_threshold)

        # plot the high contrast MTF
        if high_contrast:
            hicon_ax = next(axes)
            mtfs = [roi.mtf for roi in self.hc_rois]
            mtfs /= mtfs[0]
            self._plot_highcontrast(hicon_ax, mtfs, self.hi_contrast_threshold)

        if show:
            plt.show()

    @property
    @lru_cache()
    def _phan_region(self):
        """The skimage region of the phantom outline."""
        regions = self._get_canny_regions()
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
        return bbox_center(self._phan_region)

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


class LeedsTOR(ImagePhantomBase):
    """Class that analyzes Leeds TOR phantom planar kV images for kV QA.

    Attributes
    ----------
    lc_rois : list
        :class:`~pylinac.core.roi.LowContrastDiskROI` instances of the low
        contrast ROIs.
    lc_ref_rois : list
        :class:`~pylinac.core.roi.LowContrastDiskROI` instances of the low
        contrast reference ROIs, which are placed just inside each contrast ROI.
    hc_rois : list
        :class:`~pylinac.core.roi.HighContrastDiskROI` instances of the
        high contrast line pair regions.
    hc_ref_rois : list
        :class:`~pylinac.core.roi.HighContrastDiskROI` instances of the
        2 solid areas beside the high contrast line pair regions, which determine
        the normalized MTF value.
    """
    _demo_filename = 'leeds.dcm'

    @property
    @lru_cache()
    def _blobs(self):
        """The indices of the regions that were significant; i.e. a phantom circle outline or lead/copper square."""
        blobs = []
        for idx, region in enumerate(self._regions):
            if region.area < 100:
                continue
            round = region.eccentricity < 0.3
            if round:
                blobs.append(idx)
        if not blobs:
            raise ValueError("Could not find the phantom in the image.")
        return blobs

    @property
    @lru_cache()
    def _regions(self):
        """All the regions of the canny image that were labeled."""
        regions = self._get_canny_regions()
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
        circle_rois = [self._regions[roi] for roi in circles]
        y = np.mean([bbox_center(roi).y for roi in circle_rois])
        x = np.mean([bbox_center(roi).x for roi in circle_rois])
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
        lead_center = bbox_center(lead_roi)

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
            :class:`~pylinac.core.roi.LowContrastDistROI` instances of the contrast ROIs.
        reference ROIs : list
            :class:`~pylinac.core.roi.LowContrastDistROI` instances of the reference ROIs;
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
            :class:`~pylinac.core.roi.HighContrastDiskROI` instances of the line pairs.
        reference ROIs : list
            :class:`~pylinac.core.roi.HighContrastDiskROI` instances of the solid ROIs that
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
        self.lc_rois, self.lc_ref_rois = self._low_contrast(radius, center, angle)
        self.hc_rois, self.hc_ref_rois = self._high_contrast(radius, angle, center)

    def _flip_image_data(self, center, angle):
        """Flip the image left->right and invert the center, and angle as appropriate.

        Sometimes the Leeds phantom is set upside down on the imaging panel. Pylinac's
        analysis goes counter-clockwise, so this method flips the image and coordinates to
        make the image ccw. Quicker than flipping the image and reanalyzing.
        """
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
        if num_plots < 1:
            return
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
            for roi in self.lc_rois:
                roi.plot2axes(img_ax, edgecolor=roi.plot_color)
            for roi in self.lc_ref_rois:
                roi.plot2axes(img_ax, edgecolor='g')
            # plot the high-contrast ROIs
            for roi in self.hc_rois:
                roi.plot2axes(img_ax, edgecolor=roi.plot_color)
            for roi in self.hc_ref_rois:
                roi.plot2axes(img_ax, edgecolor='g')

        # plot the low contrast values
        if low_contrast:
            lowcon_ax = next(axes)
            self._plot_lowcontrast(lowcon_ax, self.lc_rois, self.low_contrast_threshold)

        # plot the high contrast MTF
        if high_contrast:
            hicon_ax = next(axes)
            hc_rois = [roi.mtf for roi in self.hc_rois]
            hc_rois.insert(0, 1)
            self._plot_highcontrast(hicon_ax, hc_rois, self.hi_contrast_threshold)

        if show:
            plt.show()
