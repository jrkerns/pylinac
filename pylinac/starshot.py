"""
The Starshot module analyses a starshot image made of radiation spokes, whether gantry, collimator, MLC or couch.
It is based on ideas from `Depuydt et al <http://iopscience.iop.org/0031-9155/57/10/2997>`_
and `Gonzalez et al <http://dx.doi.org/10.1118/1.1755491>`_.

Features:

* **Analyze scanned film images, single EPID images, or a set of EPID images** -
  Any image that you can load in can be analyzed, including 1 or a set of EPID DICOM images and
  films that have been digitally scanned.
* **Any image size** - Have machines with different EPIDs? Scanned your film at different resolutions? No problem.
* **Dose/OD can be inverted** - Whether your device/image views dose as an increase in value or a decrease, pylinac
  will detect it and invert if necessary.
* **Automatic noise detection & correction** - Sometimes there's dirt on the scanned film; sometimes there's a dead pixel on the EPID.
  Pylinac will detect these spurious noise signals and can avoid or account for them.
* **Accurate, FWHM star line detection** - Pylinac uses not simply the maximum value to find the center of a star line,
  but analyzes the entire star profile to determine the center of the FWHM, ensuring small noise or maximum value bias is avoided.
* **Adaptive searching** - If you passed pylinac a set of parameters and a good result wasn't found, pylinac can recover and
  do an adaptive search by adjusting parameters to find a "reasonable" wobble.
"""

from __future__ import annotations

import copy
import io
import math
import webbrowser
from itertools import product
from pathlib import Path
from typing import BinaryIO

import argue
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from pydantic import Field
from scipy import optimize

from .core import image, pdf
from .core.geometry import Circle, Line, Point
from .core.io import TemporaryZipDirectory, get_url, retrieve_demo_file
from .core.plotly_utils import set_axis_range
from .core.profile import CollapsedCircleProfile, FWXMProfile
from .core.utilities import QuaacDatum, QuaacMixin, ResultBase, ResultsDataMixin
from .core.warnings import capture_warnings


@capture_warnings
class StarshotResults(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes.
    """

    tolerance_mm: float = Field(
        description="The tolerance used for the analysis in mm."
    )
    circle_diameter_mm: float = Field(
        description="The diameter of the minimum circle that touches all the star lines in mm.",
        title="Diameter of fitted circle (mm)",
    )
    circle_radius_mm: float = Field(
        description="The radius of the minimum circle that touches all the star lines in mm.",
        title="Radius of fitted circle (mm)",
    )
    circle_center_x_y: tuple[float, float] = Field(
        description="The center position of the minimum circle in pixels.",
        title="Circle center pixel (X, Y)",
    )
    angles: list[float] = Field(
        description="The angles of the radiation lines in degrees. The angles are relative to the vertical axis and range from +/- 90 degrees.",
        title="Radiation line angles (degrees)",
    )
    passed: bool = Field(description="Whether the analysis passed or failed.")


class Starshot(ResultsDataMixin[StarshotResults], QuaacMixin):
    """Class that can determine the wobble in a "starshot" image, be it gantry, collimator,
    couch or MLC. The image can be a scanned film (TIF, JPG, etc) or a sequence of EPID DICOM images.

    Attributes
    ----------
    image : :class:`~pylinac.core.image.Image`
    circle_profile : :class:`~pylinac.starshot.StarProfile`
    lines : :class:`~pylinac.starshot.LineManager`
    wobble : :class:`~pylinac.starshot.Wobble`
    tolerance : :class:`~pylinac.starshot.Tolerance`

    Examples
    --------
    Run the demo:
        >>> Starshot.run_demo()

    Typical session:
        >>> img_path = r"C:/QA/Starshots/Coll.jpeg"
        >>> mystar = Starshot(img_path, dpi=105, sid=1000)
        >>> mystar.analyze()
        >>> print(mystar.results())
        >>> mystar.plot_analyzed_image()
    """

    angles: list[float]

    def __init__(self, filepath: str | BinaryIO, **kwargs):
        """
        Parameters
        ----------
        filepath
            The path to the image file.
        kwargs
            Passed to :func:`~pylinac.core.image.load`.
        """
        self.image = image.load(filepath, **kwargs)
        self.wobble = Wobble()
        self.tolerance = 1
        if self.image.dpmm is None:
            raise ValueError(
                "DPI was not a tag in the image nor was it passed in. Please pass a DPI value"
            )
        if self.image.sid is None:
            raise ValueError(
                "Source-to-Image distance was not an image tag and was not passed in. Please pass an SID value."
            )

    @classmethod
    def from_url(cls, url: str, **kwargs):
        """Instantiate from a URL.

        Parameters
        ----------
        url : str
            URL of the raw file.
        kwargs
            Passed to :func:`~pylinac.core.image.load`.
        """
        filename = get_url(url)
        return cls(filename, **kwargs)

    @classmethod
    def from_demo_image(cls):
        """Construct a Starshot instance and load the demo image."""
        demo_file = retrieve_demo_file(name="starshot.tif")
        return cls(demo_file, sid=1000)

    @classmethod
    def from_multiple_images(
        cls,
        filepath_list: list,
        stretch_each: bool = True,
        method: str = "sum",
        **kwargs,
    ):
        """Construct a Starshot instance and load in and combine multiple images.

        Parameters
        ----------
        filepath_list : iterable
            An iterable of file paths to starshot images that are to be superimposed.
        stretch_each : bool
            Whether to stretch each image individually before combining. See ``load_multiples``.
        method : {'sum', 'mean'}
            The method to combine the images. See ``load_multiples``.
        kwargs
            Passed to :func:`~pylinac.core.image.load_multiples`.
        """
        with io.BytesIO() as stream:
            img = image.load_multiples(
                filepath_list, stretch_each=stretch_each, method=method, **kwargs
            )
            img.save(stream)
            stream.seek(0)
            return cls(stream, **kwargs)

    @classmethod
    def from_zip(cls, zip_file: str, **kwargs):
        """Construct a Starshot instance from a ZIP archive.

        Parameters
        ----------
        zip_file : str
            Points to the ZIP archive. Can contain a single or multiple images. If multiple images
            the images are combined and thus should be from the same test sequence.
        kwargs
            Passed to :func:`~pylinac.core.image.load_multiples`.
        """
        with TemporaryZipDirectory(zip_file) as tmpdir:
            image_files = image.retrieve_image_files(tmpdir)
            if not image_files:
                raise IndexError(f"No valid starshot images were found in {zip_file}")
            if len(image_files) > 1:
                return cls.from_multiple_images(image_files, **kwargs)
            else:
                return cls(image_files[0], **kwargs)

    def _get_reasonable_start_point(self) -> (Point, float):
        """Set the algorithm starting point automatically.

        Notes
        -----
        The determination of an automatic start point is accomplished by finding the Full-Width-80%-Max.
        Finding the maximum pixel does not consistently work, esp. in the presence of a pin prick. The
        FW80M is a more consistent metric for finding a good start point.
        """
        # sum the image along each axis within the central 1/3 (avoids outlier influence from say, gantry shots)
        top_third = int(self.image.array.shape[0] / 3)
        bottom_third = int(top_third * 2)
        left_third = int(self.image.array.shape[1] / 3)
        right_third = int(left_third * 2)
        central_array = self.image.array[top_third:bottom_third, left_third:right_third]

        # use max because there is an edge case where multiple spokes
        # that are nearly vertical or horizontal can cause the sum to be higher than the center
        # allowing for a dip that can affect the initial center point
        x_sum = np.max(central_array, 0)
        y_sum = np.max(central_array, 1)

        # Calculate Full-Width, 80% Maximum center
        fwxm_x_point = (
            round(FWXMProfile(values=x_sum, fwxm_height=80).center_idx) + left_third
        )
        fwxm_y_point = (
            round(FWXMProfile(values=y_sum, fwxm_height=80).center_idx) + top_third
        )
        center_point = Point(fwxm_x_point, fwxm_y_point)
        return center_point, np.percentile(central_array, 90)

    @argue.bounds(radius=(0.2, 0.95), min_peak_height=(0.05, 0.95))
    def analyze(
        self,
        radius: float = 0.85,
        min_peak_height: float = 0.25,
        max_wobble_diameter: float = 2.0,
        tolerance: float = 1.0,
        start_point: Point | tuple | None = None,
        fwhm: bool = True,
        recursive: bool = True,
        invert: bool = False,
    ):
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
            radiation peaks that vary in magnitude (e.g. different MU delivered or gantry shot), but could also pick up noise.
            If necessary, lower value for gantry shots and increase for noisy images.
        max_wobble_diameter: float, optional
            Max wobble diameter in mm. Used only when ``recursive=True``, and in that case a larger value results in an unreasonable wobble.
        tolerance : int, float, optional
            The tolerance in mm to test against for a pass/fail result.
        start_point : 2-element iterable, optional
            The point where the algorithm should center the circle profile, given as (x-value, y-value).
            If None (default), will search for a reasonable maximum point nearest the center of the image.
        fwhm : bool
            If True (default), the center of the FWHM of the spokes will be determined.
            If False, the peak value location is used as the spoke center.

            .. note:: In practice, this ends up being a very small difference. Set to false if peak locations are offset or unexpected.
        recursive : bool
            If True (default), will recursively search for a "reasonable" wobble, meaning the wobble radius is less than max_wobble_diameter
            If the wobble found was unreasonable, the minimum peak height is iteratively adjusted from low to high at the passed radius.
            If for all peak heights at the given radius the wobble is still unreasonable, the
            radius is then iterated over from most distant inward, iterating over minimum peak heights at each radius.
            If False, will simply return the first determined value or raise error if a reasonable wobble could not be determined.

            .. warning:: It is strongly recommended to leave this setting at True.

        invert : bool
            Whether to force invert the image values. This should be set to True if the automatically-determined
            pylinac inversion is incorrect.

        Raises
        ------
        RuntimeError
            If a reasonable wobble value was not found.
        """
        self.tolerance = tolerance
        self.image.check_inversion_by_histogram(percentiles=[4, 50, 96])
        self.image.ground()
        if invert:
            self.image.invert()

        auto_point, local_max = self._get_reasonable_start_point()
        if start_point is None:
            start_point = auto_point

        self._get_reasonable_wobble(
            start_point,
            fwhm,
            min_peak_height,
            radius,
            recursive,
            local_max,
            max_wobble_diameter,
        )
        self.angles = calculate_angles(self.lines)

    def _get_reasonable_wobble(
        self,
        start_point,
        fwhm,
        min_peak_height,
        radius,
        recursive,
        local_max,
        max_wobble_diameter,
    ):
        """Determine a wobble that is "reasonable". If recursive is false, the first iteration will be passed,
        otherwise the parameters will be tweaked to search for a reasonable wobble."""
        wobble_reasonable = False
        focus_point = copy.copy(start_point)
        peak_height_recursive_candidates = np.append(
            min_peak_height, np.linspace(0.05, 0.95, 10)
        )
        radius_recursive_candidates = np.append(radius, np.linspace(0.95, 0.1, 10))
        radius_and_peak_gen = product(
            radius_recursive_candidates, peak_height_recursive_candidates
        )

        while not wobble_reasonable:
            try:
                min_height = min_peak_height * local_max
                self.circle_profile = StarProfile(
                    self.image, focus_point, radius, min_height, fwhm
                )
                # must have at least 3 spokes and we must detect an even number (e.g. gantry starshot spokes can fade off)
                if (len(self.circle_profile.peaks) < 6) or (
                    len(self.circle_profile.peaks) % 2 != 0
                ):
                    if not recursive:
                        raise RuntimeError(
                            "The algorithm was unable to properly detect the radiation lines. Try setting "
                            "recursive to True or lower the minimum peak height"
                        )
                    else:
                        raise ValueError

                self.lines = LineManager(
                    self.circle_profile.peaks,
                    focus_point=focus_point,
                    dpmm=self.image.dpmm,
                )
                self._find_wobble_minimize()

                # see if the wobble is reasonable and near the center
                # if so, we can assume we found the wobble
                focus_near_center = (
                    self.wobble.center.distance_to(focus_point) < 10 * self.image.dpmm
                )
                if (
                    (
                        self.wobble.diameter_mm < max_wobble_diameter
                        and focus_near_center
                    )
                    or not recursive  # stop after first iteration if not recursive
                ):
                    wobble_reasonable = True
                else:
                    raise ValueError

            except ValueError:
                try:
                    radius, min_peak_height = next(radius_and_peak_gen)
                except StopIteration:
                    raise RuntimeError(
                        "The algorithm was unable to determine a reasonable wobble. Try setting "
                        "recursive to False and manually adjusting algorithm parameters"
                    )

    def _find_wobble_minimize(self) -> None:
        """Find the minimum distance wobble location and radius to all radiation lines.

        The minimum is found using a scipy minimization function.
        """
        # starting point
        sp = copy.copy(self.circle_profile.center)

        def distance(p, lines):
            """Calculate the maximum distance to any line from the given point."""
            return max(line.distance_to(Point(p[0], p[1])) for line in lines)

        res = optimize.minimize(
            distance,
            sp.as_array(),
            args=(self.lines,),
            method="Nelder-Mead",
            options={"fatol": 0.001},
        )
        # res = optimize.least_squares(distance, sp.as_array(), args=(self.lines,), ftol=0.001)

        self.wobble.radius = res.fun
        self.wobble.radius_mm = res.fun / self.image.dpmm
        self.wobble.center = Point(res.x[0], res.x[1])

    @property
    def passed(self) -> bool:
        """Boolean specifying whether the determined wobble was within tolerance."""
        return bool(self.wobble.radius_mm * 2 < self.tolerance)

    @property
    def _passfail_str(self) -> str:
        """Return a pass/fail string."""
        return "PASS" if self.passed else "FAIL"

    def results(self, as_list: bool = False) -> str | list[str]:
        """Return the results of the analysis.

        Parameters
        ----------
        as_list : bool
            Whether to return as a list of strings vs single string. Pretty much for internal usage.
        """
        results = [
            " - Starshot Results - ",
            f"Result: {self._passfail_str}",
            f"The minimum circle that touches all the star lines has a diameter of {self.wobble.radius_mm*2:2.3f} mm.",
            f"The center of the minimum circle is at {self.wobble.center.x:3.1f}, {self.wobble.center.y:3.1f}",
        ]
        if not as_list:
            results = "\n".join(results)
        return results

    def _generate_results_data(self) -> StarshotResults:
        """Present the results data and metadata as a dataclass or dict.
        The default return type is a dataclass."""
        return StarshotResults(
            tolerance_mm=self.tolerance,
            circle_diameter_mm=self.wobble.radius_mm * 2,
            circle_radius_mm=self.wobble.radius_mm,
            circle_center_x_y=(self.wobble.center.x, self.wobble.center.y),
            angles=self.angles,
            passed=self.passed,
        )

    def _quaac_datapoints(self) -> dict[str, QuaacDatum]:
        """Return the data points to be saved to the QuAAC file."""
        results_data = self.results_data()
        return {
            "Circle diameter": QuaacDatum(
                value=results_data.circle_diameter_mm,
                unit="mm",
                description="The diameter of the fitted circle representing isocenter.",
            ),
        }

    def plotly_analyzed_images(
        self,
        show: bool = True,
        show_colorbar: bool = True,
        show_legend: bool = True,
        **kwargs,
    ) -> dict[str, go.Figure]:
        """Plot the analyzed set of images to Plotly figures. Will plot a zoomed-out image and a zoomed-in image.


        Parameters
        ----------
        show : bool
            Whether to show the plot.
        show_colorbar : bool
            Whether to show the colorbar on the plot.
        show_legend : bool
            Whether to show the legend on the plot.
        kwargs
            Additional keyword arguments to pass to the plot.

        Returns
        -------
        dict
            A dictionary of the Plotly figures where the key is the name of the
            image and the value is the figure.
        """
        figs = {}
        for name, zoom in zip(("Image", "Wobble"), (False, True)):
            fig = self.image.plotly(
                title="Starshot Analysis",
                show=False,
                show_legend=show_legend,
                show_colorbar=show_colorbar,
                **kwargs,
            )
            for idx, line in enumerate(self.lines):
                line.plotly(
                    fig,
                    color="blue",
                    showlegend=show_legend,
                    name=f"Line {idx} ({self.angles[idx]:2.2f}°)",
                )
            self.wobble.plotly(
                fig,
                line_color="green",
                name=f"Wobble Circle {self.wobble.diameter_mm:2.2f}mm",
                hoverinfo="text",
                hovertext=f"Wobble diameter: {self.wobble.diameter_mm:2.2f} mm",
            )
            if zoom:
                set_axis_range(
                    fig=fig,
                    x=[
                        self.wobble.center.x - self.wobble.diameter,
                        self.wobble.center.x + self.wobble.diameter,
                    ],
                    y=[
                        self.wobble.center.y - self.wobble.diameter,
                        self.wobble.center.y + self.wobble.diameter,
                    ],
                )

            figs[name] = fig
        if show:
            for f in figs.values():
                f.show()
        return figs

    def plot_analyzed_image(self, show: bool = True, **plt_kwargs: dict):
        """Draw the star lines, profile circle, and wobble circle on a matplotlib figure.

        Parameters
        ----------
        show : bool
            Whether to actually show the image.
        plt_kwargs : dict
            Keyword args passed to the plt.subplots() method. Allows one to set things like figure size.
        """
        fig, axes = plt.subplots(ncols=2, **plt_kwargs)
        subimages = ("whole", "wobble")
        titles = ("Analyzed Image", "Wobble Circle")

        # show images
        for ax, subimage, title in zip(axes, subimages, titles):
            self.plot_analyzed_subimage(ax=ax, show=False, subimage=subimage)
            ax.set_title(title)

        if show:
            plt.show()

    def plot_analyzed_subimage(
        self,
        subimage: str = "wobble",
        ax: plt.Axes | None = None,
        show: bool = True,
        **plt_kwargs: dict,
    ):
        """Plot a subimage of the starshot analysis. Current options are the zoomed out image and the zoomed in image.

        Parameters
        ----------
        subimage : str
            If 'wobble', will show a zoomed in plot of the wobble circle.
            Any other string will show the zoomed out plot.
        ax : None, matplotlib Axes
            If None (default), will create a new figure to plot on, otherwise plot to the passed axes.
        show : bool
            Whether to actually show the image.
        plt_kwargs : dict
            Keyword args passed to the plt.figure() method. Allows one to set things like figure size. Only used if ax is not passed.
        """
        if ax is None:
            fig, ax = plt.subplots(**plt_kwargs)
        # show analyzed image
        self.image.plot(ax=ax, show=False)
        self.lines.plot(ax)
        self.wobble.plot2axes(ax, edgecolor="green")
        self.circle_profile.plot2axes(ax, edgecolor="green")
        ax.autoscale(tight=True)
        ax.axis("off")

        # zoom in if wobble plot
        if subimage == "wobble":
            xlims = [
                self.wobble.center.x + self.wobble.diameter,
                self.wobble.center.x - self.wobble.diameter,
            ]
            ylims = [
                self.wobble.center.y + self.wobble.diameter,
                self.wobble.center.y - self.wobble.diameter,
            ]
            ax.set_xlim(xlims)
            ax.set_ylim(ylims)
            ax.axis("on")

        if show:
            plt.show()

    def save_analyzed_image(self, filename: str, **kwargs):
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

    def save_analyzed_subimage(self, filename: str, subimage: str = "wobble", **kwargs):
        """Save the analyzed subimage to a file.

        Parameters
        ----------
        filename : str, file-object
            Where to save the file to.
        subimage : str
            If 'wobble', will show a zoomed in plot of the wobble circle.
            Any other string will show the zoomed out plot.
        kwargs
            Passed to matplotlib.
        """
        self.plot_analyzed_subimage(subimage=subimage, show=False)
        plt.savefig(filename, **kwargs)

    def publish_pdf(
        self,
        filename: str | BinaryIO,
        notes: str | list[str] | None = None,
        open_file: bool = False,
        metadata: dict | None = None,
        logo: Path | str | None = None,
    ):
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
        logo: Path, str
            A custom logo to use in the PDF report. If nothing is passed, the default pylinac logo is used.
        """
        canvas = pdf.PylinacCanvas(
            filename, page_title="Starshot Analysis", metadata=metadata, logo=logo
        )
        for img, height in zip(("wobble", "asdf"), (2, 11.5)):
            data = io.BytesIO()
            self.save_analyzed_subimage(data, img)
            canvas.add_image(data, location=(4, height), dimensions=(13, 13))
        text = [
            "Starshot results:",
            f"Source-to-Image Distance (mm): {self.image.sid:2.0f}",
            f"Tolerance (mm): {self.tolerance:2.1f}",
            f"Minimum circle diameter (mm): {self.wobble.radius_mm*2:2.2f}",
        ]
        canvas.add_text(text=text, location=(10, 25.5), font_size=12)
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 5.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 5))
        canvas.finish()

        if open_file:
            webbrowser.open(filename)

    @staticmethod
    def run_demo():
        """Demonstrate the Starshot module using the demo image."""
        star = Starshot.from_demo_image()
        star.analyze()
        print(star.results())
        star.plot_analyzed_image()


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
    def diameter_mm(self) -> float:
        """Diameter of the wobble in mm."""
        return self.radius_mm * 2


class LineManager:
    """Manages the radiation lines found."""

    def __init__(self, points: list[Point], focus_point: Point, dpmm: float):
        """
        Parameters
        ----------
        points :
            The peak points found by the StarProfile
        """
        self.lines = []
        self.focus_point = focus_point
        self.dpmm = dpmm
        self.construct_rad_lines(points)

    def __getitem__(self, item):
        return self.lines[item]

    def __len__(self):
        return len(self.lines)

    def construct_rad_lines(self, points: list[Point]):
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
        self.match_points(points)
        # check that the lines pass near the focus point. If not, we can
        # still be in a situation where we're missing spoke-halves like a gantry starshot.
        for line in self.lines:
            if line.distance_to(self.focus_point) > 10 * self.dpmm:
                raise ValueError(
                    "The radiation lines are not near the center of the image. "
                    "This could be due to missing spoke halves, such as in a gantry starshot."
                )

    def match_points(self, points: list[Point]):
        """Match the peaks found to the same radiation lines.

        Peaks are matched by connecting the existing peaks based on an offset of peaks. E.g. if there are
        12 peaks, there must be 6 radiation lines. Furthermore, assuming star lines go all the way across the CAX,
        the 7th peak will be the opposite peak of the 1st peak, forming a line. This method is robust to
        starting points far away from the real center.
        """
        num_rad_lines = int(len(points) / 2)
        offset = num_rad_lines
        self.lines = [
            Line(points[line], points[line + offset]) for line in range(num_rad_lines)
        ]

    def plot(self, axis: plt.Axes):
        """Plot the lines to the axis."""
        for line in self.lines:
            line.plot2axes(axis, color="blue")


class StarProfile(CollapsedCircleProfile):
    """Class that holds and analyzes the circular profile which finds the radiation lines."""

    def __init__(self, image, start_point, radius, min_peak_height, fwhm):
        radius = self._convert_radius_perc2pix(image, start_point, radius)
        super().__init__(
            center=start_point,
            radius=radius,
            image_array=image.array,
            width_ratio=0.1,
            sampling_ratio=3,
        )
        self.get_peaks(min_peak_height, fwhm=fwhm)

    @staticmethod
    def _convert_radius_perc2pix(image, start_point, radius):
        """Convert a percent radius to distance in pixels, based on the distance from center point to image
            edge.

        Parameters
        ----------
        radius : float
            The radius ratio (e.g. 0.5).
        """
        return image.dist2edge_min(start_point) * radius

    def _roll_prof_to_midvalley(self) -> int:
        """Roll the circle profile so that its edges are not near a radiation line.
        This is a prerequisite for properly finding star lines.
        """
        roll_amount = np.where(self.values == self.values.min())[0][0]
        self.roll(roll_amount)
        return roll_amount

    def get_peaks(self, min_peak_height, min_peak_distance=0.02, fwhm=True):
        """Determine the peaks of the profile."""
        self._roll_prof_to_midvalley()
        self.filter(size=0.003, kind="gaussian")
        self.ground()
        if fwhm:
            self.find_fwxm_peaks(
                threshold=min_peak_height, min_distance=min_peak_distance
            )
        else:
            self.find_peaks(min_peak_height, min_peak_distance)


def calculate_angles(lines: list[Line]) -> list[float]:
    """Calculate the angles of the starshot spokes. What makes this
    somewhat annoying is that the zero-angle is defined as pointing up (vs right for a unit angle)
    and that we display the image with the y-axis increasing downward (vs upward)."""
    angles = []
    for line in lines:
        try:
            phi_rad = math.atan(line.m)
            phi_deg = math.degrees(phi_rad) - 90
            # Normalize the angle to be within (-90, +90) degrees
            if phi_deg > 90:
                phi_deg -= 180
            elif phi_deg <= -90:
                phi_deg += 180
        except ZeroDivisionError:
            phi_deg = 90
        angles.append(phi_deg)
    return angles
