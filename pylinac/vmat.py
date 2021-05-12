# -*- coding: utf-8 -*-
"""The VMAT module consists of the class VMAT, which is capable of loading an EPID DICOM Open field image and MLC field image and analyzing the
images according to the Varian RapidArc QA tests and procedures, specifically the Dose-Rate & Gantry-Speed (DRGS)
and Dose-Rate & MLC speed (DRMLC) tests.

Features:

* **Do both tests** - Pylinac can handle either DRGS or DRMLC tests.
* **Automatic offset correction** - Older VMAT tests had the ROIs offset, newer ones are centered. No worries, pylinac finds the ROIs automatically.
* **Automatic open/DMLC identification** - Pass in both images--don't worry about naming. Pylinac will automatically identify the right images.
"""
from io import BytesIO
from typing import Union, List, Tuple, Sequence, Optional

import argue
import matplotlib.pyplot as plt
import numpy as np

from .core import image
from .core.geometry import Point, Rectangle
from .core.image import ImageLike
from .core.io import get_url, TemporaryZipDirectory, retrieve_demo_file
from .core.pdf import PylinacCanvas
from .core.profile import SingleProfile
from .core.utilities import open_path
from .settings import get_dicom_cmap

DMLC = 'dmlc'
OPEN = 'open'
PROFILE = 'profile'


class VMATBase:
    _url_suffix: str
    _result_header: str
    _result_short_header: str
    SEGMENT_X_POSITIONS_MM: Tuple
    dmlc_image: image.DicomImage
    open_image: image.DicomImage
    segments: List
    _tolerance: float

    def __init__(self, image_paths: Sequence[str]):
        """
        Parameters
        ----------
        image_paths : iterable (list, tuple, etc)
            A sequence of paths to the image files.
        """
        if len(image_paths) != 2:
            raise ValueError("Exactly 2 images (open, DMLC) must be passed")
        image1, image2 = self._load_images(image_paths)
        image1, image2 = self._check_img_inversion(image1, image2)
        self._identify_images(image1, image2)
        self.segments = []
        self._tolerance = 0

    @classmethod
    def from_url(cls, url: str):
        """Load a ZIP archive from a URL.  Must follow the naming convention.

        Parameters
        ----------
        url : str
            Must point to a valid URL that is a ZIP archive of two VMAT images.
        """
        zfile = get_url(url)
        return cls.from_zip(zfile)

    @classmethod
    def from_zip(cls, path: str):
        """Load VMAT images from a ZIP file that contains both images. Must follow the naming convention.

        Parameters
        ----------
        path : str
            Path to the ZIP archive which holds the VMAT image files.
        """
        with TemporaryZipDirectory(path) as tmpzip:
            image_files = image.retrieve_image_files(tmpzip)
            return cls(image_paths=image_files)

    @classmethod
    def from_demo_images(cls):
        """Construct a VMAT instance using the demo images."""
        demo_file = retrieve_demo_file(url=cls._url_suffix)
        return cls.from_zip(demo_file)

    @argue.bounds(tolerance=(0, 8))
    def analyze(self, tolerance: Union[float, int] = 1.5, segment_size_mm: Tuple = (5, 100)):
        """Analyze the open and DMLC field VMAT images, according to 1 of 2 possible tests.

        Parameters
        ----------
        tolerance : float, int, optional
            The tolerance of the sample deviations in percent. Default is 1.5.
            Must be between 0 and 8.
        segment_size_mm : tuple(int, int)
            The (width, height) of the ROI segments in mm.
        """
        self._tolerance = tolerance/100

        """Analysis"""
        points = self._calculate_segment_centers()
        Segment._nominal_width_mm = segment_size_mm[0]
        Segment._nominal_height_mm = segment_size_mm[1]
        self._construct_segments(points)

    @staticmethod
    def _load_images(image_paths: list) -> Tuple[ImageLike, ImageLike]:
        image1 = image.load(image_paths[0])
        image2 = image.load(image_paths[1])
        image1.ground()
        image2.ground()
        return image1, image2

    @staticmethod
    def _check_img_inversion(image1: ImageLike, image2: ImageLike) -> Tuple[ImageLike, ImageLike]:
        """Check that the images are correctly inverted."""
        for image in [image1, image2]:
            image.check_inversion()
        return image1, image2

    def _identify_images(self, image1: ImageLike, image2: ImageLike):
        """Identify which image is the DMLC and which is the open field."""
        profile1, profile2 = self._median_profiles((image1, image2))
        field_profile1 = profile1.field_values()
        field_profile2 = profile2.field_values()
        if np.std(field_profile1) > np.std(field_profile2):
            self.dmlc_image = image1
            self.open_image = image2
        else:
            self.dmlc_image = image2
            self.open_image = image1

    def results(self) -> str:
        """A string of the summary of the analysis results.

        Returns
        -------
        str
            The results string showing the overall result and deviation statistics by segment.
        """
        if self.passed:
            passfail_str = 'PASS'
        else:
            passfail_str = 'FAIL'

        string = f'{self._result_header}\nTest Results (Tol. +/-{self._tolerance*100:2.2}%): {passfail_str}\n'

        string += f'Max Deviation: {self.max_r_deviation:2.3}%\nAbsolute Mean Deviation: {self.avg_abs_r_deviation:2.3}%'
        return string

    def _calculate_segment_centers(self) -> List[Point]:
        """Construct the center points of the segments based on the field center and known x-offsets."""
        points = []
        dmlc_prof, _ = self._median_profiles((self.dmlc_image, self.open_image))
        x_field_center, _ = dmlc_prof.fwxm_center(x=30)
        for x_offset_mm in self.SEGMENT_X_POSITIONS_MM:
            y = self.open_image.center.y
            x_offset_pixels = x_offset_mm * self.open_image.dpmm
            x = x_field_center + x_offset_pixels
            points.append(Point(x, y))
        return points

    def _construct_segments(self, points: List[Point]):
        for point in points:
            segment = Segment(point, self.open_image, self.dmlc_image, self._tolerance)
            self.segments.append(segment)
        # post-analysis to update R_corr values
        self._update_r_corrs()

    def _update_r_corrs(self):
        """After the Segment constructions, the R_corr must be set for each segment."""
        avg_r_corr = np.array([segment.r_corr for segment in self.segments]).mean()
        for segment in self.segments:
            segment.r_dev = ((segment.r_corr / avg_r_corr) * 100) - 100

    @property
    def passed(self) -> bool:
        return all(segment.passed for segment in self.segments)

    @property
    def r_devs(self) -> np.ndarray:
        """Return the deviations of all segments as an array."""
        return np.array([segment.r_dev for segment in self.segments])

    @property
    def avg_abs_r_deviation(self) -> float:
        """Return the average of the absolute R_deviation values."""
        return np.abs(self.r_devs).mean()

    @property
    def avg_r_deviation(self) -> float:
        """Return the average of the R_deviation values, including the sign."""
        return self.r_devs.mean()

    @property
    def max_r_deviation(self) -> float:
        """Return the value of the maximum R_deviation segment."""
        return np.max(np.abs(self.r_devs))

    def plot_analyzed_image(self, show: bool=True):
        """Plot the analyzed images. Shows the open and dmlc images with the segments drawn; also plots the median
        profiles of the two images for visual comparison.

        Parameters
        ----------
        show : bool
            Whether to actually show the image.
        """
        fig, axes = plt.subplots(ncols=3, sharex=True)
        subimages = (OPEN, DMLC, PROFILE)
        titles = ('Open', 'DMLC', 'Median Profiles')
        for subimage, axis, title in zip(subimages, axes, titles):
            self._plot_analyzed_subimage(subimage=subimage, ax=axis, show=False)
            axis.set_title(title)
        axis.set_ylabel('Normalized Response')
        axis.legend(loc='lower center')

        if show:
            plt.tight_layout(h_pad=1.5)
            plt.show()

    @argue.options(subimage=(DMLC, OPEN, PROFILE))
    def _save_analyzed_subimage(self, filename: str, subimage: str, **kwargs):
        """Save the analyzed images as a png file.

        Parameters
        ----------
        filename : str, file-object
            Where to save the file to.
        kwargs
            Passed to matplotlib.
        """
        self._plot_analyzed_subimage(subimage=subimage, show=False)
        plt.savefig(filename, **kwargs)

    @argue.options(subimage=(DMLC, OPEN, PROFILE))
    def _plot_analyzed_subimage(self, subimage: str, show: bool=True, ax: Optional[plt.Axes]=None):
        """Plot an individual piece of the VMAT analysis.

        Parameters
        ----------
        subimage : str
            Specifies which image to plot.
        show : bool
            Whether to actually plot the image.
        ax : matplotlib Axes, None
            If None (default), creates a new figure to plot to, otherwise plots to the given axes.
        """
        plt.ioff()
        if ax is None:
            fig, ax = plt.subplots()

        # plot DMLC or OPEN image
        if subimage in (DMLC, OPEN):
            if subimage == DMLC:
                img = self.dmlc_image
            elif subimage == OPEN:
                img = self.open_image
            ax.imshow(img, cmap=get_dicom_cmap())
            self._draw_segments(ax)
            plt.sca(ax)
            plt.axis('off')
            plt.tight_layout()

        # plot profile
        elif subimage == PROFILE:
            dmlc_prof, open_prof = self._median_profiles((self.dmlc_image, self.open_image))
            ax.plot(dmlc_prof.values, label='DMLC')
            ax.plot(open_prof.values, label='Open')
            ax.autoscale(axis='x', tight=True)
            ax.legend(loc=8, fontsize='large')
            ax.grid()

        if show:
            plt.show()

    def _draw_segments(self, axis: plt.Axes):
        """Draw the segments onto a plot.

        Parameters
        ----------
        axis : matplotlib.axes.Axes
            The plot to draw the objects on.
        """
        for segment in self.segments:
            color = segment.get_bg_color()
            segment.plot2axes(axis, edgecolor=color)

    @staticmethod
    def _median_profiles(images) -> Tuple[SingleProfile, SingleProfile]:
        """Return two median profiles from the open and dmlc image. For visual comparison."""
        profile1 = SingleProfile(np.mean(images[0], axis=0))
        profile1.stretch()
        profile2 = SingleProfile(np.mean(images[1], axis=0))
        profile2.stretch()

        # normalize the profiles to approximately the same value
        norm_val = np.percentile(profile1.values, 90)
        profile1.normalize(norm_val)
        norm_val = np.percentile(profile2.values, 90)
        profile2.normalize(norm_val)

        return profile1, profile2

    def publish_pdf(self, filename: str, notes: str=None, open_file: bool=False, metadata: Optional[dict]=None):
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
        canvas = PylinacCanvas(filename=filename, page_title=f"{self._result_short_header} VMAT Analysis", metadata=metadata)
        for y, x, width, img in zip((9, 9, -2), (1, 11, 3), (9, 9, 14), (OPEN, DMLC, PROFILE)):
            data = BytesIO()
            self._save_analyzed_subimage(data, subimage=img)
            canvas.add_image(data, location=(x, y), dimensions=(width, 18))
            # canvas.add_text(text=f"{img} Image", location=(x + 2, y + 10), font_size=18)
        canvas.add_text(text='Open Image', location=(4, 22), font_size=18)
        canvas.add_text(text=f'{self.open_image.base_path}', location=(4, 21.5))
        canvas.add_text(text='DMLC Image', location=(14, 22), font_size=18)
        canvas.add_text(text=f'{self.dmlc_image.base_path}', location=(14, 21.5))
        canvas.add_text(text='Median profiles', location=(8, 12), font_size=18)
        text = [f'{self._result_header} VMAT results:',
                f'Source-to-Image Distance (mm): {self.open_image.sid:2.0f}',
                f'Tolerance (%): {self._tolerance*100:2.1f}',
                f'Absolute mean deviation (%): {self.avg_abs_r_deviation:2.2f}',
                f'Maximum deviation (%): {self.max_r_deviation:2.2f}',
                ]
        canvas.add_text(text=text, location=(10, 25.5))
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 5.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 5))

        canvas.finish()

        if open_file:
            open_path(filename)


class DRGS(VMATBase):
    """Class representing a Dose-Rate, Gantry-speed VMAT test. Will accept, analyze, and return the results."""
    _url_suffix = 'drgs.zip'
    _result_header = 'Dose Rate & Gantry Speed'
    _result_short_header = 'DR/GS'
    SEGMENT_X_POSITIONS_MM = (-60, -40, -20, 0, 20, 40, 60)

    @staticmethod
    def run_demo():
        """Run the demo for the Dose Rate & Gantry Speed test."""
        vmat = DRGS.from_demo_images()
        vmat.analyze()  # old images (rev1, not new rev2's), which are offset
        print(vmat.results())
        vmat.plot_analyzed_image()


class DRMLC(VMATBase):
    """Class representing a Dose-Rate, MLC speed VMAT test. Will accept, analyze, and return the results."""
    _url_suffix = 'drmlc.zip'
    _result_header = 'Dose Rate & MLC Speed'
    _result_short_header = 'DR/MLCS'
    SEGMENT_X_POSITIONS_MM = (-45, -15, 15, 45)

    @staticmethod
    def run_demo():
        """Run the demo for the MLC leaf speed test."""
        vmat = DRMLC.from_demo_images()
        vmat.analyze()
        print(vmat.results())
        vmat.plot_analyzed_image()


class Segment(Rectangle):
    """A class for holding and analyzing segment data of VMAT tests.

    For VMAT tests, there are either 4 or 7 'segments', which represents a section of the image that received
    radiation under the same conditions.

    Attributes
    ----------
    r_dev : float
            The reading deviation (R_dev) from the average readings of all the segments. See RTD for equation info.
    r_corr : float
        The corrected reading (R_corr) of the pixel values. See RTD for explanation and equation info.
    passed : boolean
        Specifies where the segment reading deviation was under tolerance.
    """
    # width of the segment (i.e. parallel to MLC motion) in pixels under reference conditions
    _nominal_width_mm: int
    _nominal_height_mm: int

    def __init__(self, center_point: Point, open_image: image.DicomImage, dmlc_image: image.DicomImage,
                 tolerance: Union[float, int]):
        self.r_dev: float = 0.0  # is assigned after all segments constructed
        self._tolerance = tolerance
        self._open_image = open_image
        self._dmlc_image = dmlc_image
        width = self._nominal_width_mm * dmlc_image.dpmm
        height = self._nominal_height_mm * dmlc_image.dpmm
        super().__init__(width, height, center=center_point, as_int=True)

    @property
    def r_corr(self) -> float:
        """Return the ratio of the mean pixel values of DMLC/OPEN images."""
        dmlc_value = self._dmlc_image.array[self.bl_corner.y:self.bl_corner.y + self.height,
                     self.bl_corner.x: self.bl_corner.x + self.width].mean()
        open_value = self._open_image.array[self.bl_corner.y:self.bl_corner.y + self.height,
                     self.bl_corner.x: self.bl_corner.x + self.width].mean()
        ratio = (dmlc_value / open_value) * 100
        return ratio

    @property
    def passed(self) -> bool:
        """Return whether the segment passed or failed."""
        return abs(self.r_dev) < self._tolerance * 100

    def get_bg_color(self) -> str:
        """Get the background color of the segment when plotted, based on the pass/fail status."""
        return 'blue' if self.passed else 'red'
