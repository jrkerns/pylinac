"""The planar imaging module analyzes phantom images taken with the kV or MV imager in 2D.
The following phantoms are supported:

* Leeds TOR 18
* Standard Imaging QC-3
* Las Vegas

Features:

* **Automatic phantom localization** - Set up your phantom any way you like; automatic positioning,
  angle, and inversion correction mean you can set up how you like, nor will setup variations give you headache.
* **High and low contrast determination** - Analyze both low and high contrast ROIs. Set thresholds
  as you see fit.
"""
import copy
from functools import lru_cache
import io

import matplotlib.pyplot as plt
import numpy as np
from pylinac.core.geometry import Circle
from scipy.interpolate.interpolate import interp1d
from skimage import feature, measure

from .core.utilities import open_path
from .core import image
from .core.geometry import Point
from .core.io import get_url, retrieve_demo_file
from .core.profile import CollapsedCircleProfile
from .core.roi import LowContrastDiskROI, HighContrastDiskROI, DiskROI, bbox_center
from .core import pdf


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
        self.image = image.load(filepath)
        self.image.invert()

    @classmethod
    def from_demo_image(cls):
        """Instantiate and load the demo image."""
        demo_file = retrieve_demo_file(url=cls._demo_filename)
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
        self.plot_analyzed_image(show=False, **kwargs)
        plt.savefig(filename, **kwargs)

    def _get_canny_regions(self, sigma=2, percentiles=(0.001, 0.01)):
        """Compute the canny edges of the image and return the connected regions found."""
        # copy, filter, and ground the image
        img_copy = copy.copy(self.image)
        img_copy.filter(kind='gaussian', size=sigma)
        img_copy.ground()

        # compute the canny edges with very low thresholds (detects nearly everything)
        lo_th, hi_th = np.percentile(img_copy, percentiles)
        c = feature.canny(img_copy, low_threshold=lo_th, high_threshold=hi_th)

        # label the canny edge regions
        labeled = measure.label(c)
        regions = measure.regionprops(labeled, intensity_image=img_copy, coordinates='rc')
        return regions

    @staticmethod
    def _plot_lowcontrast(axes, rois, threshold):
        """Plot the low contrast ROIs to an axes."""
        # rois.sort(key=lambda x: x.contrast_constant, reverse=True)
        line1, = axes.plot([roi.contrast_constant for roi in rois], marker='o', color='m', label='Contrast Constant')
        axes.axhline(threshold, color='k')
        axes.grid(True)
        axes.set_title('Low-frequency Contrast')
        axes.set_xlabel('ROI #')
        axes.set_ylabel('Contrast Constant')
        axes2 = axes.twinx()
        line2, = axes2.plot([roi.contrast_to_noise for roi in rois], marker='^', label='CNR')
        axes2.set_ylabel('CNR')
        axes.legend(handles=[line1, line2])

    @staticmethod
    def _plot_highcontrast(axes, rois, threshold):
        """Plot the high contrast ROIs to an axes."""
        axes.plot(rois, marker='*')
        axes.axhline(threshold, color='k')
        axes.grid(True)
        axes.set_title('High-frequency rMTF')
        axes.set_xlabel('Line pair region #')
        axes.set_ylabel('relative MTF')

    def phantom_center(self):
        pass

    def phantom_angle(self):
        pass

    def phantom_ski_region(self):
        pass

    def phantom_radius(self):
        pass

    def _mtf(self, x=50):
        norm = max(roi.mtf for roi in self.hc_rois)
        ys = [roi.mtf / norm for roi in self.hc_rois]
        xs = np.arange(len(ys))
        f = interp1d(ys, xs)
        try:
            mtf = f(x / 100)
        except ValueError:
            mtf = min(ys)
        return float(mtf)


class LasVegas(ImagePhantomBase):
    """Class for analyzing low contrast of the Las Vegas MV phantom.

    Attributes
    ----------
    lc_rois : list
        :class:`~pylinac.core.roi.LowContrastDiskROI` instances of the low
        contrast ROIs, other than the reference ROI (below).
    bg_rois : list
        :class:`~pylinac.core.roi.LowContrastDiskROI` instance of the low
        contrast reference ROI (15mm PVC).
    """
    _demo_filename = 'lasvegas.dcm'

    def __init__(self, filepath):
        super().__init__(filepath)
        self._phantom_ski_region = None
        self.image.check_inversion_by_histogram()

    @staticmethod
    def run_demo():
        """Run the Las Vegas phantom analysis demonstration."""
        lv = LasVegas.from_demo_image()
        lv.analyze()
        lv.plot_analyzed_image()

    def analyze(self, low_contrast_threshold=0.1, invert=False):
        """Analyze the Las Vegas phantom.

        Parameters
        ----------
        low_contrast_threshold : float
            The threshold for the low-contrast bubbles to be "seen".
        invert : bool
            Whether to force an inversion of the image. Pylinac tries to infer the correct inversion but uneven
            backgrounds can cause this analysis to fail. If the contrasts/MTF ROIs appear correctly located but the
            plots are wonky, try setting this to True.
        """
        self.threshold = low_contrast_threshold
        if invert:
            self.image.invert()
        self._check_direction()
        self._determine_low_contrast()

    def _check_direction(self):
        """Check that the phantom is facing the right direction and if not perform a left-right flip of the array."""
        circle = CollapsedCircleProfile(self.phantom_center, self.phantom_radius * 0.175, self.image, ccw=False,
                                        width_ratio=0.16, num_profiles=5)
        circle.filter(size=0.015, kind='median')
        valleys = circle.find_valleys(max_number=2, kind='value')
        if valleys[1] > valleys[0]:
            self.image.array = np.fliplr(self.image.array)
            self._phantom_ski_region = None

    def _determine_low_contrast(self):
        """Sample the detail contrast regions."""
        # create 4 ROIs on each side of the phantom to determine the average background
        angles = np.array([-10, 80, 170, 260]) + self.phantom_angle - 4
        dists = np.array([0.24, 0.24, 0.24, 0.24]) * self.phantom_radius
        bg_rois = []
        for dist, angle in zip(dists, angles):
            roi = LowContrastDiskROI(self.image, angle, self.phantom_radius*0.03, dist, self.phantom_center,
                                     0.05)
            bg_rois.append(roi)
        avg_bg = np.mean([roi.pixel_value for roi in bg_rois])

        # create X ROIs to sample the low contrast holes
        angles = np.array([77, 116, 134.5, 0, 13, 77, 142, 153, -21, -29, -107, 182, 174, -37, -55, -105, 206.5, 189.5, -48.1, -67.8]) + self.phantom_angle
        dists = np.array([0.107, 0.141, 0.205, 0.179, 0.095, 0.042, 0.097, 0.178, 0.174, 0.088, 0.024, 0.091, 0.179, 0.189, 0.113, 0.0745, 0.115, 0.191, 0.2085, 0.146]) * self.phantom_radius
        roi_radii = np.array([0.028, 0.028, 0.028, 0.016, 0.016, 0.016, 0.016, 0.016, 0.012, 0.012, 0.012, 0.012, 0.012, 0.007, 0.007, 0.007, 0.007, 0.007, 0.003, 0.003])
        rois = []
        for dist, angle, radius in zip(dists, angles, roi_radii):
            roi = LowContrastDiskROI(self.image, angle, self.phantom_radius*radius, dist, self.phantom_center,
                                     self.threshold, avg_bg)
            rois.append(roi)

        # normalize the threshold
        self.threshold *= max(roi.contrast_constant for roi in rois)
        for roi in rois:
            roi.contrast_threshold = self.threshold
        self.bg_rois = bg_rois
        self.lc_rois = rois

    def plot_analyzed_image(self, image=True, low_contrast=True, show=True):
        """Plot the analyzed image, which includes the original image with ROIs marked and low-contrast plots.

        Parameters
        ----------
        image : bool
            Show the image.
        low_contrast : bool
            Show the low contrast values plot.
        show : bool
            Whether to actually show the image when called.
        """
        num_plots = sum((image, low_contrast))
        if num_plots < 1:
            return
        fig, axes = plt.subplots(1, num_plots)
        fig.subplots_adjust(wspace=0.4)
        if num_plots < 2:
            axes = (axes,)
        axes = iter(axes)

        if image:
            img_ax = next(axes)
            # self.image.plot(ax=img_ax, show=False, vmin=self.bg_rois[0].pixel_value*0.92, vmax=self.bg_rois[0].pixel_value*1.08)
            self.image.plot(ax=img_ax, show=False)
            img_ax.axis('off')
            img_ax.set_title('Las Vegas Phantom Analysis')

            # plot the low contrast ROIs
            for roi in self.lc_rois:
                roi.plot2axes(img_ax, edgecolor=roi.plot_color_constant)
            for roi in self.bg_rois:
                roi.plot2axes(img_ax, edgecolor='g')
            # c = Circle(self.phantom_center, radius=5)
            # c.plot2axes(img_ax, edgecolor='r')

        # plot the low contrast values
        if low_contrast:
            lowcon_ax = next(axes)
            self._plot_lowcontrast(lowcon_ax, self.lc_rois, self.threshold)

        if show:
            plt.show()

    @property
    def phantom_center(self):
        return bbox_center(self.phantom_ski_region)

    @property
    def phantom_radius(self):
        return self.phantom_ski_region.major_axis_length

    @property
    @lru_cache(maxsize=1)
    def phantom_angle(self):
        """
        Sample all sides of the phantom, searching for the low contrast circle

        Returns
        -------

        """
        circle = CollapsedCircleProfile(self.phantom_center, self.phantom_radius*0.17, self.image, ccw=False,
                                        width_ratio=0.2, num_profiles=5)
        circle.filter(size=0.01, kind='gaussian')
        angle = circle.find_valleys(max_number=2).max()
        return angle/len(circle.values) * 360

    @property
    def phantom_ski_region(self):
        """The skimage region of the phantom outline."""
        if self._phantom_ski_region is not None:
            return self._phantom_ski_region
        else:
            regions = self._get_canny_regions()
            blobs = []
            for phantom_idx, region in enumerate(regions):
                if region.area < 50:
                    continue
                hollow = region.extent < 0.02
                sized = (130 * self.image.dpmm) < region.major_axis_length < (270 * self.image.dpmm)
                if hollow and sized:
                    blobs.append(phantom_idx)

            if not blobs:
                raise ValueError("Unable to find the Las Vegas phantom in the image.")

            # find the biggest ROI and call that the phantom outline
            big_roi_idx = np.argsort([regions[phan].major_axis_length for phan in blobs])[-1]
            phantom_idx = blobs[big_roi_idx]

            self._phantom_ski_region = regions[phantom_idx]
            return regions[phantom_idx]

    def publish_pdf(self, filename, notes=None, open_file=False, metadata=None):
        """Publish (print) a PDF containing the analysis, images, and quantitative results.

        Parameters
        ----------
        filename : (str, file-like object}
            The file to write the results to.
        notes : str, list of strings
            Text; if str, prints single line.
            If list of strings, each list item is printed on its own line.
        open_file : bool
            Whether to open the file using the default program after creation.
        metadata : dict
            Extra data to be passed and shown in the PDF. The key and value will be shown with a colon.
            E.g. passing {'Author': 'James', 'Unit': 'TrueBeam'} would result in text in the PDF like:
            --------------
            Author: James
            Unit: TrueBeam
            --------------
        """
        canvas = pdf.PylinacCanvas(filename, page_title='Las Vegas Analysis', metadata=metadata)
        for (img, lo), (w, l) in zip(((True, False), (False, True)), ((5, 12), (5, 2))):
            data = io.BytesIO()
            self.save_analyzed_image(data, image=img, low_contrast=lo)
            canvas.add_image(data, location=(w, l), dimensions=(13, 13))
        text = ['Las Vegas results:',
                f'Median Contrast: {np.median([roi.contrast for roi in self.lc_rois]):2.2f}',
                f'Median CNR: {np.median([roi.contrast_to_noise for roi in self.lc_rois]):2.1f}',
                f'ROIs "seen": {sum(roi.passed_contrast_constant for roi in self.lc_rois):2.0f}',
                ]
        canvas.add_text(text=text, location=(10, 24.5))
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 5.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 5))

        canvas.finish()
        if open_file:
            open_path(filename)


class StandardImagingQC3(ImagePhantomBase):
    """Class for analyzing high and low contrast of the Standard Imaging QC-3 MV phantom.

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
    _demo_filename = 'qc3.dcm'

    @staticmethod
    def run_demo():
        """Run the Standard Imaging QC-3 phantom analysis demonstration."""
        qc3 = StandardImagingQC3.from_demo_image()
        qc3.analyze()
        qc3.plot_analyzed_image()

    def analyze(self, low_contrast_threshold=0.005, hi_contrast_threshold=0.5, invert=False):
        """Analyze the QC-3 phantom.

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
        self.image.check_inversion_by_histogram()
        if invert:
            self.image.invert()
        self.low_contrast_threshold = low_contrast_threshold
        self.hi_contrast_threshold = hi_contrast_threshold

        self.lc_ref_rois, self.lc_rois = self._low_contrast()
        self.hc_rois = self._high_contrast()

    def publish_pdf(self, filename, notes=None, open_file=False, metadata=None):
        """Publish (print) a PDF containing the analysis, images, and quantitative results.

        Parameters
        ----------
        filename : (str, file-like object}
            The file to write the results to.
        notes : str, list of strings
            Text; if str, prints single line.
            If list of strings, each list item is printed on its own line.
        open_file : bool
            Whether to open the file using the default program after creation.
        metadata : dict
            Extra data to be passed and shown in the PDF. The key and value will be shown with a colon.
            E.g. passing {'Author': 'James', 'Unit': 'TrueBeam'} would result in text in the PDF like:
            --------------
            Author: James
            Unit: TrueBeam
            --------------
        """
        canvas = pdf.PylinacCanvas(filename, page_title='QC-3 Analysis', metadata=metadata)
        for (img, lo, hi), (w, l) in zip(((True, False, False), (False, True, False), (False, False, True)),
                                         ((5, 12), (1, 4), (11, 4))):
            data = io.BytesIO()
            self.save_analyzed_image(data, image=img, low_contrast=lo)
            canvas.add_image(data, location=(w, l), dimensions=(13, 13))
        text = ['QC-3 results:',
                f'MTF 80% (lp/mm): {self._mtf(80):2.2f}',
                f'MTF 50% (lp/mm): {self._mtf(50):2.2f}',
                f'Median Contrast: {np.median([roi.contrast for roi in self.lc_rois]):2.2f}',
                f'Median CNR: {np.median([roi.contrast_to_noise for roi in self.lc_rois]):2.1f}',
                ]
        canvas.add_text(text=text, location=(10, 25.5))
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 5.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 5))

        canvas.finish()
        if open_file:
            open_path(filename)

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
            img_ax.set_title('QC-3 Phantom Analysis')

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
            mtfs /= max(mtfs)
            self._plot_highcontrast(hicon_ax, mtfs, self.hi_contrast_threshold)

        if show:
            plt.show()

    @property
    @lru_cache()
    def phantom_ski_region(self):
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

        if not blobs:
            raise ValueError("Unable to find the QC-3 phantom in the image.")

        # find the biggest ROI and call that the phantom outline
        big_roi_idx = np.argsort([regions[phan].major_axis_length for phan in blobs])[-1]
        phantom_idx = blobs[big_roi_idx]

        return regions[phantom_idx]

    @property
    def phantom_radius(self):
        """The radius of the phantom in pixels; the value itself doesn't matter, it's just
        used for relative distances to ROIs.

        Returns
        -------
        radius : float
        """
        return self.phantom_ski_region.major_axis_length / 14

    @property
    def phantom_angle(self):
        """The angle of the phantom.

        Returns
        -------
        angle : float
            The angle in degrees.
        """
        return -np.rad2deg(self.phantom_ski_region.orientation)

    @property
    def phantom_center(self):
        """The center point of the phantom.

        Returns
        -------
        center : Point
        """
        return bbox_center(self.phantom_ski_region)

    def _low_contrast(self):
        """Sample the detail contrast regions."""
        angles = np.array([90, -90, 55, -55, 128, -128]) + self.phantom_angle
        dists = np.array([2, 2, 2.4, 2.4, 2.4, 2.4]) * self.phantom_radius
        rrois = []

        # background ROI
        bg_roi = LowContrastDiskROI(self.image, angles[0], 0.5 * self.phantom_radius, dists[0], self.phantom_center, 0.05)

        for dist, angle in zip(dists[1:], angles[1:]):
            roi = LowContrastDiskROI(self.image, angle, 0.5 * self.phantom_radius, dist, self.phantom_center,
                                     0.05, background=bg_roi.pixel_value)
            rrois.append(roi)
        return bg_roi, rrois

    def _high_contrast(self):
        """Sample the high-contrast line pair regions."""
        dists = np.array([2.8, -2.8, 1.45, -1.45, 0]) * self.phantom_radius
        rrois = []
        for dist in dists:
            roi = HighContrastDiskROI(self.image, self.phantom_angle, 0.5 * self.phantom_radius, dist, self.phantom_center,
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
    _phantom_angle = None
    _phantom_center = None

    @property
    @lru_cache(1)
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
    @lru_cache(1)
    def _regions(self):
        """All the regions of the canny image that were labeled."""
        return self._get_canny_regions()

    @property
    def phantom_center(self):
        """Determine the phantom center.

        This is done by searching for circular ROIs of the canny image. Those that are circular and roughly the
        same size as the biggest circle ROI are all sampled for the center of the bounding box. The values are
        averaged over all the detected circles to give a more robust value.

        Returns
        -------
        center : Point
        """
        if self._phantom_center is not None:
            return self._phantom_center
        circles = [roi for roi in self._blobs if
                   np.isclose(self._regions[roi].major_axis_length, self.phantom_radius * 3.35, rtol=0.3)]

        # get average center of all circles
        circle_rois = [self._regions[roi] for roi in circles]
        y = np.mean([bbox_center(roi).y for roi in circle_rois])
        x = np.mean([bbox_center(roi).x for roi in circle_rois])
        return Point(x, y)

    @property
    @lru_cache()
    def phantom_angle(self):
        """Determine the angle of the phantom.

        This is done by searching for square-like boxes of the canny image. There are usually two: one lead and
        one copper. The box with the highest intensity (lead) is identified. The angle from the center of the lead
        square bounding box and the phantom center determines the phantom angle.

        Returns
        -------
        angle : float
            The angle in radians.
        """
        circle = CollapsedCircleProfile(self.phantom_center, self.phantom_radius * 0.79, self.image,
                                        width_ratio=0.04, ccw=True)
        circle.ground()
        circle.filter(size=0.01)
        peak_idx = circle.find_fwxm_peaks(threshold=0.6, max_number=1)[0]
        shift_percent = peak_idx / len(circle.values)
        shift_radians = shift_percent * 2 * np.pi
        shift_radians_corrected = 2*np.pi - shift_radians
        return shift_radians_corrected

    @property
    @lru_cache(1)
    def phantom_radius(self):
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

    def _is_counter_clockwise(self):
        """Determine if the low-contrast bubbles go from high to low clockwise or counter-clockwise.

        Returns
        -------
        boolean
        """
        circle = CollapsedCircleProfile(self.phantom_center, self.phantom_radius * 0.79, self.image, width_ratio=0.04, ccw=True)
        circle.ground()
        circle.filter(size=0.01)
        circle.values = np.roll(circle.values, -circle.values.argmax())
        first_set = circle.find_peaks(search_region=(0.05, 0.45), threshold=0, min_distance=0.025, kind='value', max_number=9)
        second_set = circle.find_peaks(search_region=(0.55, 0.95), threshold=0, min_distance=0.025, kind='value', max_number=9)
        return max(first_set) > max(second_set)

    def _low_contrast(self, angle_offset):
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
        phantom_angle = np.degrees(self.phantom_angle) - angle_offset
        roi_set1_angles = np.linspace(30, 150, num=9)
        roi_set2_angles = np.linspace(210, 330, num=9)
        roi_nominal_angles = np.concatenate((roi_set1_angles, roi_set2_angles))
        roi_determined_angles = phantom_angle + roi_nominal_angles

        bubble_radius = 0.025 * self.phantom_radius
        bubble_dist = 0.785 * self.phantom_radius

        # sample the Contrast ROIs
        crois = []
        for angle in roi_determined_angles:
            roi = LowContrastDiskROI(self.image, angle, bubble_radius, bubble_dist, self.phantom_center, self.low_contrast_threshold)
            crois.append(roi)

        # sample the reference ROIs (closer to the center of phantom)
        bubble_dist = 0.65 * self.phantom_radius
        rrois = []
        for idx, angle in enumerate(roi_determined_angles):
            roi = DiskROI(self.image, angle, bubble_radius, bubble_dist, self.phantom_center)
            crois[idx].background = roi.pixel_value
            rrois.append(roi)

        return crois, rrois

    def _high_contrast(self, angle_offset):
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
        phantom_angle = np.degrees(self.phantom_angle) - angle_offset

        # sample ROIs of the reference areas
        # sampling these will give a reference min and max to normalize the line pairs to
        ref_angles = [57, 89]
        ref_dists = [0.3 * self.phantom_radius, 0.25 * self.phantom_radius]
        ref_radius = 0.04 * self.phantom_radius
        rrois = []
        for nominal_angle, dist in zip(ref_angles, ref_dists):
            roi = HighContrastDiskROI(self.image, phantom_angle - nominal_angle, ref_radius, dist, self.phantom_center,
                                      self.hi_contrast_threshold)
            rrois.append(roi)
        mtf_norm_val = (rrois[0].max - rrois[1].min) / (rrois[0].max + rrois[1].min)

        # sample ROIs of each line pair region
        # ordering goes from the thickest line pair region downward
        contrast_angles = np.array([-144.8, -115.1, -62.5, -169.7, -153.4, -25, 173, 156, 27]) + 90
        contrast_dists = np.array([0.3, 0.187, 0.187, 0.252, 0.092, 0.094, 0.252, 0.094, 0.0958]) * self.phantom_radius
        contrast_radii = np.array([0.04, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.018, 0.018, 0.015, 0.015, 0.012]) * self.phantom_radius
        crois = []
        for nominal_angle, dist, cradius in zip(contrast_angles, contrast_dists, contrast_radii):
            roi = HighContrastDiskROI(self.image, phantom_angle - nominal_angle, cradius, dist, self.phantom_center, self.hi_contrast_threshold, mtf_norm=mtf_norm_val)
            crois.append(roi)

        return crois, rrois

    @staticmethod
    def run_demo():
        """Run the Leeds TOR phantom analysis demonstration."""
        leeds = LeedsTOR.from_demo_image()
        leeds.analyze()
        leeds.plot_analyzed_image()

    def analyze(self, low_contrast_threshold=0.005, hi_contrast_threshold=0.4, invert=False,
                angle_offset=0):
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
        angle_offset : int, float
            Some LeedsTOR phantoms have the low contrast regions slightly offset from phantom to phantom. 
            This parameter lets the user correct for any consistent angle effects of the phantom. The offset 
            is in degrees and moves counter-clockwise. Use this if the low contrast ROI samples are offset from the real
            ROIs.
        """
        self.image.check_inversion_by_histogram()
        if invert:
            self.image.invert()
        self.low_contrast_threshold = low_contrast_threshold
        self.hi_contrast_threshold = hi_contrast_threshold

        if self._is_counter_clockwise():
            self._flip_image_data()
        self.lc_rois, self.lc_ref_rois = self._low_contrast(angle_offset)
        self.hc_rois, self.hc_ref_rois = self._high_contrast(angle_offset)

    def _flip_image_data(self):
        """Flip the image left->right and invert the center, and angle as appropriate.

        Sometimes the Leeds phantom is set upside down on the imaging panel. Pylinac's
        analysis goes counter-clockwise, so this method flips the image and coordinates to
        make the image ccw. Quicker than flipping the image and reanalyzing.
        """
        self.image.array = np.fliplr(self.image.array)
        new_x = self.image.shape[1] - self.phantom_center.x
        self._phantom_center = Point(new_x, self.phantom_center.y)

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

            # plot the phantom location circle; helps users identify if the phantom location was identified correctly
            Circle(center_point=self.phantom_center, radius=self.phantom_radius).plot2axes(img_ax, edgecolor='blue')

            # plot the low contrast ROIs
            for roi in self.lc_rois:
                roi.plot2axes(img_ax, edgecolor=roi.plot_color_constant)
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

    def publish_pdf(self, filename, notes=None, open_file=False, metadata=None):
        """Publish (print) a PDF containing the analysis, images, and quantitative results.

        Parameters
        ----------
        filename : (str, file-like object}
            The file to write the results to.
        notes : str, list of strings
            Text; if str, prints single line.
            If list of strings, each list item is printed on its own line.
        open_file : bool
            Whether to open the file using the default program after creation.
        metadata : dict
            Extra data to be passed and shown in the PDF. The key and value will be shown with a colon.
            E.g. passing {'Author': 'James', 'Unit': 'TrueBeam'} would result in text in the PDF like:
            --------------
            Author: James
            Unit: TrueBeam
            --------------
        """
        canvas = pdf.PylinacCanvas(filename, page_title='Leeds TOR-18 Analysis', metadata=metadata)
        for (img, lo, hi), (w, l) in zip(((True, False, False), (False, True, False), (False, False, True)),
                                         ((5, 12), (1, 4), (11, 4))):
            data = io.BytesIO()
            self.save_analyzed_image(data, image=img, low_contrast=lo)
            canvas.add_image(data, location=(w, l), dimensions=(10, 10))
        text = ['Leeds TOR18 results:',
                f'MTF 80% (lp/mm): {self._mtf(80):2.2f}',
                f'Median Contrast: {np.median([roi.contrast for roi in self.lc_rois]):2.2f}',
                f'Median CNR: {np.median([roi.contrast_to_noise for roi in self.lc_rois]):2.1f}',
                ]
        canvas.add_text(text=text, location=(10, 25.5))
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 5.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 5))

        canvas.finish()
        if open_file:
            open_path(filename)
