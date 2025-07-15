from __future__ import annotations

import copy
import io
import webbrowser
from collections.abc import Sequence
from pathlib import Path
from typing import Literal, Union

import matplotlib.pyplot as plt
import numpy as np
from pydantic import Field

from . import Centering, Edge, Normalization
from .core import image, pdf
from .core.exceptions import NotAnalyzed
from .core.geometry import Point, Rectangle
from .core.image import ImageLike
from .core.io import retrieve_demo_file
from .core.profile import (
    FWXMProfilePhysical,
    HillProfilePhysical,
    InflectionDerivativeProfilePhysical,
)
from .core.roi import RectangleROI
from .core.utilities import ResultBase, ResultsDataMixin, convert_to_enum
from .core.warnings import capture_warnings
from .metrics.profile import (
    CAXToLeftEdgeMetric,
    CAXToRightEdgeMetric,
    FlatnessDifferenceMetric,
    PenumbraLeftMetric,
    PenumbraRightMetric,
    ProfileMetric,
    SymmetryPointDifferenceMetric,
)


class FieldProfileResult(ResultBase):
    x_metrics: dict = Field(
        title="X Metrics",
        description="""
      The metrics computed from the X/crossplane profile.
      The items included depend on the metrics given to the analysis. See: :ref:`profile_builtin_plugins`. The
      key will be the name of the metric and the value will be what's returned from
      the ``calculate`` method of the metric. In addition, the following items
      are always given:

      * ``Field Width (mm)``: The width of the field in mm.
      * ``values``: A list of pixel values of the profile that the metrics were calculated from.""",
    )
    y_metrics: dict = Field(
        title="Y Metrics",
        description="""
      The metrics computed from the X/crossplane profile.
      The items included depend on the metrics given to the analysis. See: :ref:`profile_builtin_plugins`. The
      key will be the name of the metric and the value will be what's returned from
      the ``calculate`` method of the metric. In addition, the following items
      are always given:

      * ``Field Width (mm)``: The width of the field in mm.
      * ``values``: A list of pixel values of the profile that the metrics were calculated from.""",
    )
    center: dict = Field(title="Center ROI", description="The center ROI of the field")
    normalization: str = Field(
        title="Normalization", description="The normalization method used."
    )
    edge_type: str = Field(
        title="Edge Type", description="The edge detection method used."
    )
    centering: str = Field(title="Centering", description="The centering method used.")


DEFAULT_METRICS = (
    FlatnessDifferenceMetric(),
    SymmetryPointDifferenceMetric(),
    PenumbraRightMetric(),
    PenumbraLeftMetric(),
    CAXToLeftEdgeMetric(),
    CAXToRightEdgeMetric(),
)
PROFILES = {
    Edge.FWHM: FWXMProfilePhysical,
    Edge.INFLECTION_HILL: HillProfilePhysical,
    Edge.INFLECTION_DERIVATIVE: InflectionDerivativeProfilePhysical,
}
PROFILE_TYPE = Union[
    FWXMProfilePhysical, HillProfilePhysical, InflectionDerivativeProfilePhysical
]


@capture_warnings
class FieldProfileAnalysis(ResultsDataMixin[FieldProfileResult]):
    x_profile: PROFILE_TYPE
    y_profile: PROFILE_TYPE
    center_rect: RectangleROI
    x_rect: Rectangle
    y_rect: Rectangle
    _normalization: Normalization
    _edge_type: Edge
    _is_analyzed: bool = False

    def __init__(self, path: str | Path, **kwargs):
        """Field analysis of a radiation field via profiles.

        Parameters
        ----------
        path
            The path to the image to analyze.
        kwargs
            Keyword arguments to pass to the image loader.
        """
        super().__init__()
        self.image: ImageLike = image.load(path, **kwargs)
        self.image.check_inversion_by_histogram()

    @classmethod
    def from_demo_image(cls):
        """Load the demo image into an instance."""
        demo_file = retrieve_demo_file(name="flatsym_demo.dcm")
        return cls(demo_file)

    def analyze(
        self,
        centering: Centering = Centering.BEAM_CENTER,
        position: tuple[float, float] = (0.5, 0.5),
        x_width: float | int = 0.0,
        y_width: float | int = 0.0,
        normalization: Normalization = Normalization.NONE,
        edge_type: Edge = Edge.INFLECTION_DERIVATIVE,
        invert: bool = False,
        ground: bool = True,
        metrics: Sequence[ProfileMetric] = DEFAULT_METRICS,
        **kwargs,
    ):
        """Analyze the field by pulling out profiles at the specified position and width. The profiles are then analyzed
        for the given metrics.

        Parameters
        ----------
        centering
            The centering method to use.
        position
            The relative position of the field to pull the profile in (height, width).
            E.g. (0.3, 0.6) will pull profiles at 30% of the image height (from the top)
            and 60% of the image width (from the left).

            .. note::

                This is only used if centering is set to 'manual'.

        x_width
            The x-width ratio of the field to extract into the profile. Must be between 0 and 1. The pixels are averaged
            along the width to create the profile.
        y_width
            The y-width ratio of the field to extract into the profile. Must be between 0 and 1. The pixels are averaged
            along the height to create the profile.
        normalization
            The normalization method to use.
        edge_type
            The type of profile to use. Options are 'FWXM', 'Hill', and 'Inflection'.
        invert
            Whether to invert the image before analyzing.
        ground
            Whether to "ground" the profile by subtracting the minimum value from all values.
        metrics
            The metrics to compute on the profiles. Should be a sequence of individual metrics.
        """
        if invert:
            self.image.invert()
        self._normalization = convert_to_enum(normalization, Normalization)
        self._edge_type = convert_to_enum(edge_type, Edge)
        self._centering = convert_to_enum(centering, Centering)

        x_values, y_values = self._get_profile_values(position, x_width, y_width)

        self.x_profile = PROFILES[self._edge_type](
            values=x_values,
            dpmm=self.image.dpmm,
            normalization=normalization,
            ground=ground,
            **kwargs,
        )
        self.x_profile.compute(metrics=metrics)
        self.y_profile = PROFILES[self._edge_type](
            values=y_values,
            dpmm=self.image.dpmm,
            normalization=normalization,
            ground=ground,
            **kwargs,
        )
        # we have to deep copy otherwise the metrics get re-used and the y-metrics override the x-metrics
        self.y_profile.compute(metrics=copy.deepcopy(metrics))
        self._is_analyzed = True

    def _generate_results_data(self) -> FieldProfileResult:
        if not self._is_analyzed:
            raise NotAnalyzed("Image is not analyzed yet. Use analyze() first.")
        return FieldProfileResult(
            edge_type=self._edge_type,
            normalization=self._normalization,
            centering=self._centering,
            x_metrics=self.x_profile.metric_values
            | {
                "Field Width (mm)": self.x_profile.field_width_mm,
                "values": self.x_profile.values.tolist(),
            },
            y_metrics=self.y_profile.metric_values
            | {
                "Field Width (mm)": self.y_profile.field_width_mm,
                "values": self.y_profile.values.tolist(),
            },
            center={
                "mean": self.center_rect.mean,
                "stdev": self.center_rect.std,
                "min": self.center_rect.min,
                "max": self.center_rect.max,
            },
        )

    def results(self) -> str:
        """Return a string representation of the results."""
        d = self.results_data(by_alias=True, as_dict=True)
        s = ""
        for key, value in d.items():
            if isinstance(value, dict):
                s = s + f"{key}:\n"
                for k, v in value.items():
                    if not isinstance(v, list):
                        s = s + f"{k}: {v}\n"
            else:
                s = s + f"{key}: {value}\n"
        return s

    def plot_analyzed_images(
        self,
        show: bool = True,
        show_field_edges: bool = True,
        show_center: bool = True,
        show_grid: bool = True,
        mirror: Literal["beam", "geometry"] | None = None,
        **kwargs,
    ) -> list[plt.Figure]:
        """Plot the analyzed image and x and y profiles.

        Parameters
        ----------
        show
            Whether to show the plots.
        show_field_edges
            Whether to show the edges of the field of the profiles as vertical lines.
        show_center
            Whether to show the center of the field as a vertical line on the profiles.
        show_grid
            Whether to display a grid for the profiles.
        mirror
            Whether to mirror the image. Options are 'beam', 'geometry', or None.
        kwargs
            Additional keyword arguments to pass to the plt.subplots() method.

        Returns
        -------
        list[plt.Figure]
            The figures of the x profile, y profile, and image.
        """
        if not self._is_analyzed:
            raise NotAnalyzed("Image is not analyzed yet. Use analyze() first.")
        xfig, xax = plt.subplots(**kwargs)
        xax.set_title("X Profile")
        self.x_profile.plot(
            axis=xax,
            show=False,
            show_field_edges=show_field_edges,
            show_center=show_center,
            show_grid=show_grid,
            mirror=mirror,
        )
        xfig.tight_layout()
        yfig, yax = plt.subplots(**kwargs)
        self.y_profile.plot(
            axis=yax,
            show=False,
            show_field_edges=show_field_edges,
            show_center=show_center,
            show_grid=show_grid,
            mirror=mirror,
        )
        yax.set_title("Y Profile")
        yfig.tight_layout()
        ifig, ax = plt.subplots(**kwargs)
        self.image.plot(ax=ax, show=False)
        ax.set_title("Image")
        self.x_rect.plot2axes(
            ax, edgecolor="b", fill=True, alpha=0.3, facecolor="b", label="X Profile"
        )
        self.y_rect.plot2axes(
            ax, edgecolor="g", fill=True, alpha=0.3, facecolor="g", label="Y Profile"
        )
        self.center_rect.plot2axes(
            ax, edgecolor="r", fill=False, alpha=0.3, facecolor="b", label="Center ROI"
        )
        ax.legend()
        if show:
            plt.show()
        return [xfig, yfig, ifig]

    def _get_profile_values(
        self, position: tuple[float, float], x_width: float, y_width: float
    ) -> tuple[np.ndarray, np.ndarray]:
        """Get the array values of the x and y profiles. Accounts for
        width, edges, position, etc"""
        x, y = self._get_x_y_position(position)
        if x_width > 1 or x_width < 0 or y_width > 1 or y_width < 0:
            raise ValueError("Width must be between 0 and 1")
        # The +2, -1 is to ensure we have at least 2 pixel width
        # The +2 is because the indexing is exclusive
        top = round(y - self.image.shape[0] * x_width / 2 - 1)
        bottom = round(max(y + self.image.shape[0] * x_width / 2, top + 2))
        left = round(x - self.image.shape[1] * y_width / 2 - 1)
        right = round(max(x + self.image.shape[1] * y_width / 2, left + 2))
        x_box = self.image[top:bottom, :]
        y_box = self.image[:, left:right]
        # we multiply by 2 because the center can be offset and falling off the edge;
        # for this purpose is okay. 2x will cover the image no matter where the
        # center is.
        self.x_rect = Rectangle(
            width=x_box.shape[1] * 2, height=x_box.shape[0], center=(x, y)
        )
        self.y_rect = Rectangle(
            width=y_box.shape[1], height=y_box.shape[0] * 2, center=(x, y)
        )
        self.center_rect = RectangleROI(
            array=self.image.array,
            width=right - left,
            height=bottom - top,
            center=Point(x, y),
        )
        x_values = x_box.mean(axis=0)
        y_values = y_box.mean(axis=1)
        return x_values, y_values

    def _get_x_y_position(self, position: tuple[float, float]) -> (float, float):
        """Get the x and y position of the field"""
        centering = self._centering
        if centering != Centering.MANUAL:
            v_sum = self.image.array.sum(axis=0)
            h_sum = self.image.array.sum(axis=1)
            v_p: FWXMProfilePhysical = PROFILES[self._edge_type](
                values=v_sum, dpmm=self.image.dpmm
            )
            h_p: FWXMProfilePhysical = PROFILES[self._edge_type](
                values=h_sum, dpmm=self.image.dpmm
            )
            if centering == Centering.BEAM_CENTER:
                x = v_p.center_idx
                y = h_p.center_idx
            elif centering == Centering.GEOMETRIC_CENTER:  # i.e. cax
                x = v_p.cax_index
                y = h_p.cax_index
            return x, y
        else:
            # manually-specified position
            if len(position) != 2:
                raise ValueError("Position must be a tuple of two values")
            if any([pos < 0 or pos > 1 for pos in position]):
                raise ValueError("Position values must be between 0 and 1")
            # shape swapped so we get col (x), row (y).
            return self.image.shape[1] * position[1], self.image.shape[0] * position[0]

    def publish_pdf(
        self,
        filename: str,
        notes: str | list[str] = None,
        open_file: bool = False,
        metadata: dict = None,
        logo: Path | str | None = None,
        plot_kwargs: dict | None = None,
    ) -> None:
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
            Extra stream to be passed and shown in the PDF. The key and value will be shown with a colon.
            E.g. passing {'Author': 'James', 'Unit': 'TrueBeam'} would result in text in the PDF like:
            --------------
            Author: James
            Unit: TrueBeam
            --------------
        logo: Path, str
            A custom logo to use in the PDF report. If nothing is passed, the default pylinac logo is used.
        plot_kwargs : dict
            Additional keyword arguments to pass to the plot_analyzed_image method. E.g. ``show_grid=False``.
        """
        plt.ioff()
        if not self._is_analyzed:
            raise NotAnalyzed("Image is not analyzed yet. Use analyze() first.")
        canvas = pdf.PylinacCanvas(
            filename,
            page_title="Field Analysis",
            metadata=metadata,
            metadata_location=(2, 5),
            logo=logo,
        )
        # draw result text
        data = self.results_data(
            as_dict=True, by_alias=True, exclude={"pylinac_version"}
        )
        # drop values; that's too much info
        data["x_metrics"].pop("values")
        data["y_metrics"].pop("values")
        offset = 0
        for i, (key, value) in enumerate(data.items()):
            if isinstance(value, str):
                canvas.add_text(
                    text=f"{key}: {value}", location=(1, 25 - offset), font_size=12
                )
                offset += 0.75
            elif isinstance(value, dict):
                canvas.add_text(text=f"{key}:", location=(1, 25 - offset), font_size=12)
                offset += 0.75
                for j, (subkey, subvalue) in enumerate(value.items()):
                    canvas.add_text(
                        text=f"{subkey}: {subvalue:.3f}",
                        location=(2, 25 - offset),
                        font_size=12,
                    )
                    offset += 0.75

        # generate figures
        plot_kwargs = plot_kwargs or {}
        figs = self.plot_analyzed_images(show=False, **plot_kwargs)
        for fig in figs[::-1]:
            canvas.add_new_page()
            with io.BytesIO() as stream:
                fig.savefig(stream, format="png")
                stream.seek(0)
                canvas.add_image(stream, location=(-4, 13), dimensions=(28, 12))

        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 5.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 5))
        canvas.finish()

        if open_file:
            webbrowser.open(filename)
