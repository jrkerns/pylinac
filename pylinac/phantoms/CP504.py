"""CatPhan 503 and 504 phantom module classes.

Physical modules:
  CTP404 — HU linearity, geometry, slice thickness
  CTP486 — Image uniformity
  CTP528 — 21 lp/cm high-resolution gauge (line-pair MTF)
  CTP515 — Low contrast (supra-slice and sub-slice targets)
"""

from __future__ import annotations

from functools import cached_property

import matplotlib.pyplot as plt
import numpy as np
from plotly import graph_objects as go
from py_linq import Enumerable

from ..core.mtf import MTF, BeadMTF
from ..core.nps import (
    average_power,
    max_frequency,
    noise_power_spectrum_1d,
    noise_power_spectrum_2d,
)
from ..core.plotly_utils import add_title
from ..core.profile import CollapsedCircleProfile
from ..core.roi import AnnularROI, LowContrastDiskROI, RectangleROI
from ..core.warnings import capture_warnings
from ..ct import (
    AIR,
    ACRYLIC,
    DELRIN,
    LDPE,
    POLY,
    PMP,
    RAMP_ANGLE_RATIO,
    TEFLON,
    WATER,
    CatPhanBase,
    CatPhanModule,
    CTP528 as _CTP528Base,
    GeometricLine,
    Point,
    Slice,
    ThicknessROI,
    get_regions,
)


# ---------------------------------------------------------------------------
# CTP404 — HU linearity, geometry, slice thickness
# ---------------------------------------------------------------------------

class CTP404(CatPhanModule):
    """HU linearity, geometry, and slice thickness for the CTP404 module (CatPhan 503/504)."""

    attr_name = "ctp404"
    common_name = "HU Linearity"
    roi_dist_mm = 58.7
    roi_radius_mm = 5
    roi_settings = {
        "Air": {
            "value": AIR,
            "angle": -90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "PMP": {
            "value": PMP,
            "angle": -120,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "LDPE": {
            "value": LDPE,
            "angle": 180,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Poly": {
            "value": POLY,
            "angle": 120,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Acrylic": {
            "value": ACRYLIC,
            "angle": 60,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Delrin": {
            "value": DELRIN,
            "angle": 0,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Teflon": {
            "value": TEFLON,
            "angle": -60,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
    }
    background_roi_settings = {
        "1": {"angle": -30, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "2": {"angle": -150, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "3": {"angle": -210, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "4": {"angle": 30, "distance": roi_dist_mm, "radius": roi_radius_mm},
    }
    # thickness
    thickness_roi_height = 40
    thickness_roi_width = 10
    thickness_roi_distance_mm = 38
    thickness_roi_settings = {
        "Left": {
            "angle": 180,
            "width": thickness_roi_width,
            "height": thickness_roi_height,
            "distance": thickness_roi_distance_mm,
        },
        "Bottom": {
            "angle": 90,
            "width": thickness_roi_height,
            "height": thickness_roi_width,
            "distance": thickness_roi_distance_mm,
        },
        "Right": {
            "angle": 0,
            "width": thickness_roi_width,
            "height": thickness_roi_height,
            "distance": thickness_roi_distance_mm,
        },
        "Top": {
            "angle": -90,
            "width": thickness_roi_height,
            "height": thickness_roi_width,
            "distance": thickness_roi_distance_mm,
        },
    }
    # geometry
    geometry_roi_size_mm = 35
    geometry_roi_settings = {
        "Top-Horizontal": (0, 1),
        "Bottom-Horizontal": (2, 3),
        "Left-Vertical": (0, 2),
        "Right-Vertical": (1, 3),
    }
    pad: str | int
    thickness_image: Slice

    def __init__(
        self,
        catphan,
        offset: int,
        hu_tolerance: float,
        thickness_tolerance: float,
        scaling_tolerance: float,
        clear_borders: bool = True,
        thickness_slice_straddle: str | int = "auto",
        expected_hu_values: dict[str, float | int] | None = None,
    ):
        self.mm_per_pixel = catphan.mm_per_pixel
        self.hu_tolerance = hu_tolerance
        self.thickness_tolerance = thickness_tolerance
        self.scaling_tolerance = scaling_tolerance
        self.thickness_rois = {}
        self.lines = {}
        self.thickness_slice_straddle = thickness_slice_straddle
        self.expected_hu_values = expected_hu_values
        super().__init__(
            catphan, tolerance=hu_tolerance, offset=offset, clear_borders=clear_borders
        )

    def preprocess(self, catphan) -> None:
        if (
            isinstance(self.thickness_slice_straddle, str)
            and self.thickness_slice_straddle.lower() == "auto"
        ):
            if float(catphan.dicom_stack.metadata.SliceThickness) < 3.5:
                self.pad = 1
            else:
                self.pad = 0
        else:
            self.pad = self.thickness_slice_straddle
        self.thickness_image = Slice(
            catphan,
            combine_method="mean",
            num_slices=self.num_slices + self.pad,
            slice_num=self.slice_num,
            clear_borders=self.clear_borders,
        ).image

    def _replace_hu_values(self):
        if self.expected_hu_values is not None:
            for name, value in self.expected_hu_values.items():
                if name in self.roi_settings:
                    self.roi_settings[name]["value"] = value

    def _setup_rois(self) -> None:
        self._replace_hu_values()
        super()._setup_rois()
        self._setup_thickness_rois()
        if len(self.geometry_roi_settings) > 0:
            self._setup_geometry_rois()

    def _setup_thickness_rois(self) -> None:
        for name, setting in self.thickness_roi_settings.items():
            self.thickness_rois[name] = ThicknessROI.from_phantom_center(
                self.thickness_image,
                setting["width_pixels"],
                setting["height_pixels"],
                setting["angle_corrected"],
                setting["distance_pixels"],
                self.phan_center,
            )

    def _setup_geometry_rois(self) -> None:
        boxsize = self.geometry_roi_size_mm / self.mm_per_pixel
        xbounds = (int(self.phan_center.x - boxsize), int(self.phan_center.x + boxsize))
        ybounds = (int(self.phan_center.y - boxsize), int(self.phan_center.y + boxsize))
        geo_img = self.image[ybounds[0] : ybounds[1], xbounds[0] : xbounds[1]].copy()
        geo_img -= np.median(geo_img)
        nearest_extreme = min(abs(geo_img.max()), abs(geo_img.min()))
        geo_clipped = np.clip(geo_img, a_min=-nearest_extreme, a_max=nearest_extreme)
        geo_clipped_abs = np.abs(geo_clipped)
        larr, regionprops, num_roi = get_regions(
            geo_clipped_abs, fill_holes=True, clear_borders=False
        )
        if num_roi < 4:
            raise ValueError("Unable to locate the Geometric nodes")
        elif num_roi > 4:
            regionprops = sorted(
                regionprops, key=lambda x: x.filled_area, reverse=True
            )[:4]
        sorted_regions = sorted(
            regionprops, key=lambda x: 2 * x.centroid[0] + x.centroid[1]
        )
        centers = [
            Point(
                r.weighted_centroid[1] + xbounds[0], r.weighted_centroid[0] + ybounds[0]
            )
            for r in sorted_regions
        ]
        for name, order in self.geometry_roi_settings.items():
            self.lines[name] = GeometricLine(
                centers[order[0]],
                centers[order[1]],
                self.mm_per_pixel,
                self.scaling_tolerance,
            )

    @property
    def lcv(self) -> float:
        """The low-contrast visibility"""
        return (
            2
            * abs(self.rois["LDPE"].pixel_value - self.rois["Poly"].pixel_value)
            / (self.rois["LDPE"].std + self.rois["Poly"].std)
        )

    def plotly_linearity(
        self, plot_delta: bool = True, show_legend: bool = True
    ) -> go.Figure:
        fig = go.Figure()
        nominal_x_values = [roi.nominal_val for roi in self.rois.values()]
        if plot_delta:
            values = [roi.value_diff for roi in self.rois.values()]
            nominal_measurements = [0] * len(values)
            ylabel = "HU Delta"
        else:
            values = [roi.pixel_value for roi in self.rois.values()]
            nominal_measurements = nominal_x_values
            ylabel = "Measured Values"
        fig.add_scatter(
            x=nominal_x_values,
            y=values,
            name="Measured values",
            mode="markers",
            marker=dict(color="green", size=10, symbol="cross", line=dict(width=1)),
        )
        fig.add_scatter(
            x=nominal_x_values,
            y=nominal_measurements,
            mode="lines",
            name="Nominal Values",
            marker_color="green",
        )
        fig.add_scatter(
            x=nominal_x_values,
            y=np.array(nominal_measurements) + self.hu_tolerance,
            mode="lines",
            name="Upper Tolerance",
            line=dict(dash="dash", color="red"),
        )
        fig.add_scatter(
            x=nominal_x_values,
            y=np.array(nominal_measurements) - self.hu_tolerance,
            mode="lines",
            name="Lower Tolerance",
            line=dict(dash="dash", color="red"),
        )
        fig.update_layout(
            xaxis_title="Nominal Values", yaxis_title=ylabel, showlegend=show_legend
        )
        add_title(fig, "HU Linearity")
        return fig

    def plot_linearity(
        self, axis: plt.Axes | None = None, plot_delta: bool = True
    ) -> tuple:
        """Plot the HU linearity values to an axis."""
        nominal_x_values = [roi.nominal_val for roi in self.rois.values()]
        if axis is None:
            fig, axis = plt.subplots()
        if plot_delta:
            values = [roi.value_diff for roi in self.rois.values()]
            nominal_measurements = [0] * len(values)
            ylabel = "HU Delta"
        else:
            values = [roi.pixel_value for roi in self.rois.values()]
            nominal_measurements = nominal_x_values
            ylabel = "Measured Values"
        points = axis.plot(nominal_x_values, values, "g+", markersize=15, mew=2)
        axis.plot(nominal_x_values, nominal_measurements)
        axis.plot(
            nominal_x_values, np.array(nominal_measurements) + self.hu_tolerance, "r--"
        )
        axis.plot(
            nominal_x_values, np.array(nominal_measurements) - self.hu_tolerance, "r--"
        )
        axis.margins(0.05)
        axis.grid(True)
        axis.set_xlabel("Nominal Values")
        axis.set_ylabel(ylabel)
        axis.set_title("HU linearity")
        return points

    @property
    def passed_hu(self) -> bool:
        """Boolean specifying whether all the ROIs passed within tolerance."""
        return all(roi.passed for roi in self.rois.values())

    def plotly_rois(self, fig: go.Figure) -> None:
        super().plotly_rois(fig)
        for name, roi in self.thickness_rois.items():
            roi.plotly(fig, line_color="blue", name=f"Ramp {name}")
        for name, line in self.lines.items():
            line.plotly(fig, color=line.pass_fail_color, name=f"Geometry {name}")

    def plot_rois(self, axis: plt.Axes) -> None:
        super().plot_rois(axis)
        for roi in self.thickness_rois.values():
            roi.plot2axes(axis, edgecolor="blue")
        for line in self.lines.values():
            line.plot2axes(axis, color=line.pass_fail_color)

    @property
    def passed_thickness(self) -> bool:
        """Whether the slice thickness was within tolerance from nominal."""
        return (
            self.slice_thickness - self.thickness_tolerance
            < self.meas_slice_thickness
            < self.slice_thickness + self.thickness_tolerance
        )

    @property
    def meas_slice_thickness(self) -> float:
        """The average slice thickness for the wire measurements in mm."""
        return np.mean(
            sorted(
                roi.wire_fwhm * self.mm_per_pixel * RAMP_ANGLE_RATIO
                for roi in self.thickness_rois.values()
            )
        ) / (1 + 2 * self.pad)

    @property
    def avg_line_length(self) -> float:
        return float(np.mean([line.length_mm for line in self.lines.values()]))

    @property
    def passed_geometry(self) -> bool:
        """Returns whether all the line lengths were within tolerance."""
        return all(line.passed for line in self.lines.values())


# CTP404CP503 is the same physical module as CTP404 — no geometric difference.
class CTP404CP503(CTP404):
    """CTP404 module as found in CatPhan 503 (alias for backward compatibility)."""

    pass


# ---------------------------------------------------------------------------
# CTP486 — Image uniformity
# ---------------------------------------------------------------------------

class CTP486(CatPhanModule):
    """Image uniformity module (CTP486), shared across CatPhan 503/504/600/604/700."""

    attr_name = "ctp486"
    common_name = "HU Uniformity"
    roi_dist_mm = 53
    roi_radius_mm = 10
    nominal_value = 0
    nps_rois: dict[str, RectangleROI]
    roi_settings = {
        "Top": {
            "value": nominal_value,
            "angle": -90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Right": {
            "value": nominal_value,
            "angle": 0,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Bottom": {
            "value": nominal_value,
            "angle": 90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Left": {
            "value": nominal_value,
            "angle": 180,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Center": {
            "value": nominal_value,
            "angle": 0,
            "distance": 0,
            "radius": roi_radius_mm,
        },
    }

    def plot_profiles(self, axis: plt.Axes | None = None) -> None:
        """Plot the horizontal and vertical profiles of the Uniformity slice."""
        if axis is None:
            fig, axis = plt.subplots()
        horiz_data = self.image[int(self.phan_center.y), :]
        vert_data = self.image[:, int(self.phan_center.x)]
        axis.plot(horiz_data, "g", label="Horizontal")
        axis.plot(vert_data, "b", label="Vertical")
        axis.autoscale(tight=True)
        axis.axhline(self.nominal_value + self.tolerance, color="r", linewidth=3)
        axis.axhline(self.nominal_value - self.tolerance, color="r", linewidth=3)
        axis.grid(True)
        axis.set_ylabel("HU")
        axis.legend(loc=8, fontsize="small", title="")
        axis.set_title("Uniformity Profiles")

    def _setup_rois(self) -> None:
        super()._setup_rois()
        self.nps_rois = {}
        for name, setting in self.roi_settings.items():
            self.nps_rois[name] = RectangleROI.from_phantom_center(
                array=self.image,
                width=setting["radius_pixels"] * 2,
                height=setting["radius_pixels"] * 2,
                angle=setting["angle_corrected"],
                dist_from_center=setting["distance_pixels"],
                phantom_center=self.phan_center,
            )

    def plot(self, axis: plt.Axes):
        for nps_roi in self.nps_rois.values():
            nps_roi.plot2axes(axis, edgecolor="green", linestyle="-.")
        super().plot(axis)

    def plotly(self, **kwargs) -> go.Figure:
        fig = super().plotly(**kwargs)
        for name, nps_roi in self.nps_rois.items():
            nps_roi.plotly(
                fig, line_color="green", line_dash="dash", name=f"NPS {name}"
            )
        return fig

    @property
    def overall_passed(self) -> bool:
        return all(roi.passed for roi in self.rois.values())

    @property
    def uniformity_index(self) -> float:
        """The Uniformity Index. Elstrom et al equation 2."""
        center = self.rois["Center"]
        uis = [
            100 * ((roi.pixel_value - center.pixel_value) / (center.pixel_value + 1000))
            for roi in self.rois.values()
        ]
        abs_uis = np.abs(uis)
        return uis[np.argmax(abs_uis)]

    @property
    def integral_non_uniformity(self) -> float:
        """The Integral Non-Uniformity. Elstrom et al equation 1."""
        maxhu = max(roi.pixel_value for roi in self.rois.values())
        minhu = min(roi.pixel_value for roi in self.rois.values())
        return (maxhu - minhu) / (maxhu + minhu + 2000)

    # NPS region sizes (mm); converted to pixels at runtime via mm_per_pixel.
    # big_roi_size: central ROI from which sub-ROIs are extracted (IEC 62220 method).
    # small_roi_size: size of each overlapping sub-ROI.
    nps_big_roi_mm: float = 100.0
    nps_small_roi_mm: float = 20.0

    @cached_property
    def power_spectrum_2d(self) -> tuple[np.ndarray, float]:
        big_px = max(4, int(self.nps_big_roi_mm / self.mm_per_pixel))
        small_px = max(2, int(self.nps_small_roi_mm / self.mm_per_pixel))
        return noise_power_spectrum_2d(
            image_array=self.image.array.astype(float),
            pixel_size=self.mm_per_pixel,
            big_roi_size=big_px,
            small_roi_size=small_px,
        )

    @cached_property
    def power_spectrum_1d(self) -> np.ndarray:
        nps2d, _ = self.power_spectrum_2d
        return noise_power_spectrum_1d(nps2d)

    @property
    def avg_noise_power(self) -> float:
        return average_power(self.power_spectrum_1d)

    @property
    def max_noise_power_frequency(self) -> float:
        """The frequency of the maximum noise power. 0 means no pattern."""
        return max_frequency(self.power_spectrum_1d)


# ---------------------------------------------------------------------------
# CTP528 — 21 lp/cm high-resolution line-pair gauge
# ---------------------------------------------------------------------------

class CTP528CP504(_CTP528Base):
    """21 lp/cm high-resolution spatial resolution module for CatPhan 504.

    MTF is computed from the point-source bead embedded in the module insert,
    using a background-subtracted 2D PSF → radial average → Hanning → FFT
    pipeline (see :class:`~pylinac.core.mtf.BeadMTF`).  The line-pair gauge
    is retained for visual reference (``circle_profile``), but is no longer
    used for the MTF value.
    """

    attr_name: str = "ctp528"
    common_name: str = "Spatial Resolution"
    radius2linepairs_mm = 47
    combine_method: str = "mean"
    num_slices: int = 3
    boundaries: tuple[float, ...] = (
        0,
        0.107,
        0.173,
        0.236,
        0.286,
        0.335,
        0.387,
        0.434,
        0.479,
    )
    start_angle: float = np.pi
    ccw: bool = True

    # --- Bead PSF settings ---
    BEAD_DIAMETER_MM: float = 0.18
    PATCH_RADIUS_MM: float = 3.0
    BACKGROUND_INNER_MM: float = 3.5
    BACKGROUND_OUTER_MM: float = 6.0
    # Exclude the line-pair gauge ring and the outer shell; bead sits inside.
    BEAD_SEARCH_MAX_RADIUS_MM: float = 40.0

    roi_settings = {
        "region 1": {
            "start": boundaries[0],
            "end": boundaries[1],
            "num peaks": 2,
            "num valleys": 1,
            "peak spacing": 0.021,
            "gap size (cm)": 0.5,
            "lp/mm": 0.1,
        },
        "region 2": {
            "start": boundaries[1],
            "end": boundaries[2],
            "num peaks": 3,
            "num valleys": 2,
            "peak spacing": 0.01,
            "gap size (cm)": 0.25,
            "lp/mm": 0.2,
        },
        "region 3": {
            "start": boundaries[2],
            "end": boundaries[3],
            "num peaks": 4,
            "num valleys": 3,
            "peak spacing": 0.006,
            "gap size (cm)": 0.167,
            "lp/mm": 0.3,
        },
        "region 4": {
            "start": boundaries[3],
            "end": boundaries[4],
            "num peaks": 4,
            "num valleys": 3,
            "peak spacing": 0.00557,
            "gap size (cm)": 0.125,
            "lp/mm": 0.4,
        },
        "region 5": {
            "start": boundaries[4],
            "end": boundaries[5],
            "num peaks": 4,
            "num valleys": 3,
            "peak spacing": 0.004777,
            "gap size (cm)": 0.1,
            "lp/mm": 0.5,
        },
        "region 6": {
            "start": boundaries[5],
            "end": boundaries[6],
            "num peaks": 5,
            "num valleys": 4,
            "peak spacing": 0.00398,
            "gap size (cm)": 0.083,
            "lp/mm": 0.6,
        },
        "region 7": {
            "start": boundaries[6],
            "end": boundaries[7],
            "num peaks": 5,
            "num valleys": 4,
            "peak spacing": 0.00358,
            "gap size (cm)": 0.071,
            "lp/mm": 0.7,
        },
        "region 8": {
            "start": boundaries[7],
            "end": boundaries[8],
            "num peaks": 5,
            "num valleys": 4,
            "peak spacing": 0.0027866,
            "gap size (cm)": 0.063,
            "lp/mm": 0.8,
        },
    }

    def _setup_rois(self):
        pass

    def _convert_units_in_settings(self):
        pass

    @cached_property
    def bead_center(self) -> Point:
        """Sub-pixel location of the PSF bead inside the CTP528 insert.

        Searches the inner region of the module image (radius <
        ``BEAD_SEARCH_MAX_RADIUS_MM``) to avoid the line-pair gauge ring
        at ~47 mm from the phantom centre.
        """
        arr = self.image.array.astype(float)
        y_idx, x_idx = np.indices(arr.shape)
        r_mm = (
            np.sqrt(
                (x_idx - self.phan_center.x) ** 2
                + (y_idx - self.phan_center.y) ** 2
            )
            * self.mm_per_pixel
        )
        search = arr.copy()
        search[r_mm > self.BEAD_SEARCH_MAX_RADIUS_MM] = -np.inf

        if np.all(np.isinf(search)):
            raise ValueError(
                "Could not locate the PSF bead: search region is fully masked. "
                "Check BEAD_SEARCH_MAX_RADIUS_MM or phantom alignment."
            )

        peak_y, peak_x = np.unravel_index(np.argmax(search), search.shape)

        half_win = max(3, int(np.ceil(self.PATCH_RADIUS_MM / self.mm_per_pixel / 2)))
        y0 = max(0, peak_y - half_win)
        y1 = min(arr.shape[0], peak_y + half_win + 1)
        x0 = max(0, peak_x - half_win)
        x1 = min(arr.shape[1], peak_x + half_win + 1)
        win = arr[y0:y1, x0:x1].copy()
        win_pos = np.clip(win - win.min(), 0.0, None)
        total = win_pos.sum()
        if total > 0:
            ys, xs = np.indices(win.shape)
            cx = float((xs * win_pos).sum() / total) + x0
            cy = float((ys * win_pos).sum() / total) + y0
        else:
            cx, cy = float(peak_x), float(peak_y)

        return Point(x=cx, y=cy)

    @cached_property
    def _psf_data(self) -> tuple[np.ndarray, tuple[float, float]]:
        """Background-subtracted 2D PSF patch and sub-pixel bead position."""
        arr = self.image.array.astype(float)
        cx, cy = self.bead_center.x, self.bead_center.y

        patch_px = int(np.ceil(self.PATCH_RADIUS_MM / self.mm_per_pixel))
        iy, ix = int(round(cy)), int(round(cx))
        y0 = max(0, iy - patch_px)
        y1 = min(arr.shape[0], iy + patch_px + 1)
        x0 = max(0, ix - patch_px)
        x1 = min(arr.shape[1], ix + patch_px + 1)

        patch = arr[y0:y1, x0:x1].copy()

        y_full, x_full = np.indices(arr.shape)
        r_from_bead = (
            np.sqrt((x_full - cx) ** 2 + (y_full - cy) ** 2) * self.mm_per_pixel
        )
        bg_mask = (r_from_bead >= self.BACKGROUND_INNER_MM) & (
            r_from_bead <= self.BACKGROUND_OUTER_MM
        )
        background = float(arr[bg_mask].mean()) if bg_mask.any() else 0.0
        patch -= background
        return patch, (cy - y0, cx - x0)

    @property
    def psf_patch(self) -> np.ndarray:
        """Background-subtracted 2D PSF patch centred on the bead."""
        return self._psf_data[0]

    @cached_property
    def mtf(self) -> BeadMTF:
        """Bead PSF MTF."""
        patch, _ = self._psf_data
        return BeadMTF(patch, self.mm_per_pixel)

    @property
    def radius2linepairs(self) -> float:
        return self.radius2linepairs_mm * self.scaling_factor / self.mm_per_pixel

    def plotly_rois(self, fig: go.Figure) -> None:
        theta = np.linspace(0, 2 * np.pi, 100)
        r_px = self.PATCH_RADIUS_MM / self.mm_per_pixel
        fig.add_scatter(
            x=(self.bead_center.x + r_px * np.cos(theta)).tolist(),
            y=(self.bead_center.y + r_px * np.sin(theta)).tolist(),
            mode="lines",
            line=dict(color="blue"),
            name="PSF patch",
            showlegend=False,
        )

    def plot_rois(self, axis: plt.Axes) -> None:
        circle = plt.Circle(
            (self.bead_center.x, self.bead_center.y),
            radius=self.PATCH_RADIUS_MM / self.mm_per_pixel,
            color="blue",
            fill=False,
            linewidth=1.5,
        )
        axis.add_patch(circle)

    @cached_property
    def circle_profile(self) -> CollapsedCircleProfile:
        """Collapsed circle profile through the line-pair gauge (for visual reference)."""
        circle_profile = CollapsedCircleProfile(
            self.phan_center,
            self.radius2linepairs,
            image_array=self.image,
            start_angle=self.start_angle + np.deg2rad(self.catphan_roll),
            width_ratio=0.04 * self.roi_size_factor,
            sampling_ratio=2,
            ccw=self.ccw,
        )
        circle_profile.filter(0.001, kind="gaussian")
        circle_profile.ground()
        return circle_profile


class CTP528CP503(CTP528CP504):
    """CTP528 with CatPhan 503-specific line-pair boundary positions."""

    start_angle = 0
    ccw = False
    boundaries = (0, 0.111, 0.176, 0.240, 0.289, 0.339, 0.390, 0.436, 0.481)


# ---------------------------------------------------------------------------
# CTP515 — Low contrast
# ---------------------------------------------------------------------------

class CTP515(CatPhanModule):
    """Low contrast module (CTP515) for CatPhan 504.

    Measures supra-slice and sub-slice contrast targets. Contrast is multiplied
    by the target diameter to approximate human detectability (Rose model).
    """

    attr_name = "ctp515"
    common_name = "Low Contrast"
    num_slices = 1
    roi_dist_mm = 50
    roi_radius_mm = [6, 3.5, 3, 2.5, 2, 1.5]
    roi_angles = [-87.4, -69.1, -52.7, -38.5, -25.1, -12.9]
    roi_settings = {
        "15": {
            "angle": roi_angles[0],
            "distance": roi_dist_mm,
            "radius": roi_radius_mm[0],
        },
        "9": {
            "angle": roi_angles[1],
            "distance": roi_dist_mm,
            "radius": roi_radius_mm[1],
        },
        "8": {
            "angle": roi_angles[2],
            "distance": roi_dist_mm,
            "radius": roi_radius_mm[2],
        },
        "7": {
            "angle": roi_angles[3],
            "distance": roi_dist_mm,
            "radius": roi_radius_mm[3],
        },
        "6": {
            "angle": roi_angles[4],
            "distance": roi_dist_mm,
            "radius": roi_radius_mm[4],
        },
        "5": {
            "angle": roi_angles[5],
            "distance": roi_dist_mm,
            "radius": roi_radius_mm[5],
        },
    }
    background_roi_dist_ratio = 0.75
    background_roi_radius_mm = 4
    WINDOW_SIZE = 50

    def __init__(
        self,
        catphan,
        tolerance: float,
        cnr_threshold: float,
        offset: int,
        contrast_method: str,
        visibility_threshold: float,
        clear_borders: bool = True,
    ):
        self.cnr_threshold = cnr_threshold
        self.contrast_method = contrast_method
        self.visibility_threshold = visibility_threshold
        super().__init__(
            catphan, tolerance=tolerance, offset=offset, clear_borders=clear_borders
        )

    def _setup_rois(self) -> None:
        # Two annular background crowns at the same radii as the former per-rod
        # disk pairs (ratio 0.75 × and 1.25 × rod distance), with radial width
        # equal to the former disk diameter.
        dist_px = next(iter(self.roi_settings.values()))["distance_pixels"]
        half_width_px = self.background_roi_radius_mm / self.mm_per_pixel

        inner_r = dist_px * self.background_roi_dist_ratio
        outer_r = dist_px * (2 - self.background_roi_dist_ratio)

        self.background_rois["inner"] = AnnularROI(
            array=self.image,
            center=self.phan_center,
            inner_radius=inner_r - half_width_px,
            outer_radius=inner_r + half_width_px,
        )
        self.background_rois["outer"] = AnnularROI(
            array=self.image,
            center=self.phan_center,
            inner_radius=outer_r - half_width_px,
            outer_radius=outer_r + half_width_px,
        )

        background_val = float(np.mean([
            self.background_rois["inner"].pixel_value,
            self.background_rois["outer"].pixel_value,
        ]))

        for name, setting in self.roi_settings.items():
            self.rois[name] = LowContrastDiskROI.from_phantom_center(
                self.image,
                setting["angle_corrected"],
                setting["radius_pixels"],
                setting["distance_pixels"],
                self.phan_center,
                contrast_reference=background_val,
                cnr_threshold=self.cnr_threshold,
                contrast_method=self.contrast_method,
                visibility_threshold=self.visibility_threshold,
            )

    @property
    def rois_visible(self) -> int:
        return sum(roi.passed_visibility for roi in self.rois.values())

    @property
    def window_min(self) -> float:
        return (
            Enumerable(self.background_rois.values()).min(lambda r: r.pixel_value)
            - self.WINDOW_SIZE
        )

    @property
    def window_max(self) -> float:
        return (
            Enumerable(self.rois.values()).max(lambda r: r.pixel_value)
            + self.WINDOW_SIZE
        )


# ---------------------------------------------------------------------------
# CatPhan phantom classes
# ---------------------------------------------------------------------------

@capture_warnings
class CatPhan503(CatPhanBase):
    """CatPhan 503. Analyzes: HU linearity (CTP404), spatial resolution (CTP528), uniformity (CTP486)."""

    _demo_url = "CatPhan503.zip"
    _model = "503"
    catphan_radius_mm = 97
    modules = {
        CTP404CP503: {"offset": 0},
        CTP486: {"offset": -110},
        CTP528CP503: {"offset": -30},
    }

    @staticmethod
    def run_demo(show: bool = True):
        cbct = CatPhan503.from_demo_images()
        cbct.analyze()
        print(cbct.results())
        cbct.plot_analyzed_image(show)


@capture_warnings
class CatPhan504(CatPhanBase):
    """CatPhan 504. Analyzes: HU linearity (CTP404), spatial resolution (CTP528), uniformity (CTP486), low contrast (CTP515)."""

    _demo_url = "CatPhan504.zip"
    _model = "504"
    catphan_radius_mm = 101
    modules = {
        CTP404: {"offset": 0},
        CTP486: {"offset": -65},
        CTP528CP504: {"offset": 30},
        CTP515: {"offset": -30},
    }

    @staticmethod
    def run_demo(show: bool = True):
        cbct = CatPhan504.from_demo_images()
        cbct.analyze()
        print(cbct.results())
        cbct.plot_analyzed_image(show)


# ---------------------------------------------------------------------------
# Backward-compatibility aliases (old pylinac names)
# ---------------------------------------------------------------------------

#: Old name for :class:`CTP404`.
CTP404CP504 = CTP404

#: Old name for :class:`CTP528CP504`.
CTP528 = CTP528CP504
