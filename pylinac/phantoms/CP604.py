"""CatPhan 604 phantom module classes.

Physical modules:
  CTP732 — HU linearity (adds 50%/20% Bone inserts vs CTP404)
  CTP729 — Image uniformity (identical to CTP486)
  CTP528CP604 — Slanted-wire MTF (tilted tungsten wire, 5° to z-axis)
  CTP730 — Low contrast (27 rods, 3 contrast groups)
"""

from __future__ import annotations

from functools import cached_property

import matplotlib.pyplot as plt
import numpy as np
from plotly import graph_objects as go

from ..core.mtf import BeadMTF, SlantedWireMTF
from ..core.warnings import capture_warnings
from ..ct import (
    BONE_20,
    BONE_50,
    AIR,
    ACRYLIC,
    DELRIN,
    LDPE,
    POLY,
    PMP,
    TEFLON,
    CatPhanBase,
    CTP528 as _CTP528Base,
    Point,
    Slice,
    ThicknessROI,
)
from . import CP504


# ---------------------------------------------------------------------------
# CTP732 — HU linearity (CatPhan 604)
# ---------------------------------------------------------------------------

class CTP732(CP504.CTP404):
    """HU linearity module for CatPhan 604 (CTP732).

    Adds 50% Bone (at −150°) and 20% Bone (at +30°) inserts compared with the
    standard CTP404.  Only two background ROIs instead of four.

    References:
        `CP604 Manual, p.10 <CP604Manual.pdf#page=10>`_
    """

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
        "50% Bone": {
            "value": BONE_50,
            "angle": -150,
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
        "20% Bone": {
            "value": BONE_20,
            "angle": 30,
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
        "2": {"angle": -210, "distance": roi_dist_mm, "radius": roi_radius_mm},
    }


# ---------------------------------------------------------------------------
# CTP729 — Image uniformity (identical hardware to CTP486)
# ---------------------------------------------------------------------------

class CTP729(CP504.CTP486):
    """Image uniformity module for CatPhan 604 (CTP729).

    Functionally and geometrically identical to CTP486.

    References:
        `CP604 Manual, p.26 <CP604Manual.pdf#page=26>`_
    """

    pass


# ---------------------------------------------------------------------------
# CTP528CP604 — Slanted-wire MTF
# ---------------------------------------------------------------------------

class CTP528CP604(_CTP528Base):
    """Slanted-wire MTF for the CatPhan 604.

    Section 2 of the 604 phantom contains a 50 µm tungsten wire tilted ~5°
    to the z-axis.  The wire's cross-section drifts laterally by
    ``tan(5°) × slice_spacing`` per slice, providing sub-pixel oversampling
    of the PSF.  This class collects ``WIRE_SLICES`` slices on each side of
    the reference slice, registers patches at sub-pixel precision using the
    known tilt-angle shift, averages them into a high-SNR 2D PSF, then
    computes the MTF via :class:`~pylinac.core.mtf.SlantedWireMTF`.

    References:
        `CP604 Manual, p.10 <CP604Manual.pdf#page=10>`_
    """

    attr_name: str = "ctp528"
    common_name: str = "Spatial Resolution"
    combine_method: str = "mean"
    num_slices: int = 0
    start_angle = None
    roi_settings: dict = {}

    PATCH_RADIUS_MM: float = 3.0
    BACKGROUND_INNER_MM: float = 3.5
    BACKGROUND_OUTER_MM: float = 6.0
    HU_PLUG_DIST_MM: float = 58.7
    HU_PLUG_EXCLUSION_MM: float = 12.0
    TILT_ANGLE_DEG: float = 5.0
    WIRE_SLICES: int = 7  # slices on each side of the reference slice

    def preprocess(self, catphan) -> None:
        self._dicom_stack = catphan.dicom_stack

    def _setup_rois(self) -> None:
        pass

    def _convert_units_in_settings(self) -> None:
        pass

    @cached_property
    def wire_center(self) -> Point:
        """Sub-pixel location of the tungsten wire in the reference slice.

        Searches the interior of the phantom, excluding the HU plug annulus
        at ~58.7 mm from centre.
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
        phantom_radius_mm = np.sqrt(self.catphan_size / np.pi) * self.mm_per_pixel

        in_plug_annulus = np.abs(r_mm - self.HU_PLUG_DIST_MM) < self.HU_PLUG_EXCLUSION_MM
        outside_phantom = r_mm > phantom_radius_mm * 0.95

        search = arr.copy()
        search[in_plug_annulus | outside_phantom] = -np.inf

        if np.all(np.isinf(search)):
            raise ValueError(
                "Could not locate the wire: the entire search region was masked. "
                "Verify the phantom model and slice orientation."
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
        """Single-slice background-subtracted PSF patch (for display / diagnostics)."""
        arr = self.image.array.astype(float)
        cx, cy = self.wire_center.x, self.wire_center.y

        patch_px = int(np.ceil(self.PATCH_RADIUS_MM / self.mm_per_pixel))
        iy, ix = int(round(cy)), int(round(cx))
        y0 = max(0, iy - patch_px)
        y1 = min(arr.shape[0], iy + patch_px + 1)
        x0 = max(0, ix - patch_px)
        x1 = min(arr.shape[1], ix + patch_px + 1)

        patch = arr[y0:y1, x0:x1].copy()

        y_full, x_full = np.indices(arr.shape)
        r_from_wire = (
            np.sqrt((x_full - cx) ** 2 + (y_full - cy) ** 2) * self.mm_per_pixel
        )
        bg_mask = (r_from_wire >= self.BACKGROUND_INNER_MM) & (
            r_from_wire <= self.BACKGROUND_OUTER_MM
        )
        background = float(arr[bg_mask].mean()) if bg_mask.any() else 0.0
        patch -= background
        return patch, (cy - y0, cx - x0)

    @property
    def psf_patch(self) -> np.ndarray:
        """Single-slice background-subtracted PSF patch centred on the wire."""
        return self._psf_data[0]

    @cached_property
    def mtf(self) -> SlantedWireMTF:
        """Slanted-wire MTF built from a multi-slice PSF stack."""
        stack = [
            self._dicom_stack[self.slice_num + i].array.astype(float)
            for i in range(-self.WIRE_SLICES, self.WIRE_SLICES + 1)
            if 0 <= self.slice_num + i < len(self._dicom_stack)
        ]
        return SlantedWireMTF(
            image_stack=stack,
            wire_center=(self.wire_center.x, self.wire_center.y),
            pixel_size_mm=self.mm_per_pixel,
            slice_spacing_mm=self.slice_spacing,
            tilt_angle_deg=self.TILT_ANGLE_DEG,
            patch_radius_px=int(np.ceil(self.PATCH_RADIUS_MM / self.mm_per_pixel)),
        )

    def plot_rois(self, axis: plt.Axes) -> None:
        circle = plt.Circle(
            (self.wire_center.x, self.wire_center.y),
            radius=self.PATCH_RADIUS_MM / self.mm_per_pixel,
            color="blue",
            fill=False,
            linewidth=1.5,
        )
        axis.add_patch(circle)

    def plotly_rois(self, fig: go.Figure) -> None:
        theta = np.linspace(0, 2 * np.pi, 100)
        r_px = self.PATCH_RADIUS_MM / self.mm_per_pixel
        fig.add_scatter(
            x=(self.wire_center.x + r_px * np.cos(theta)).tolist(),
            y=(self.wire_center.y + r_px * np.sin(theta)).tolist(),
            mode="lines",
            line=dict(color="blue"),
            name="Wire PSF patch",
            showlegend=False,
        )


# ---------------------------------------------------------------------------
# CTP730 — Low contrast (27 rods, 3 contrast groups)
# ---------------------------------------------------------------------------

class CTP730(CP504.CTP515):
    """Low-contrast module for CatPhan 604 (CTP730).

    Contains three contrast groups (0.3 %, 0.5 %, 1.0 %) each with 9 rod
    diameters (2, 3, 4, 5, 6, 7, 8, 9, 15 mm).  All rods are POSITIVE
    contrast (denser than the LDPE background — they appear brighter in CT).

    Angles follow pylinac image convention: 0°=right, −90°=top, +90°=bottom.
    Positions calibrated from a 41-slice average of a real CatPhan 604 scan
    (0.625 mm/slice, 0.78 mm/pixel).  Group 3 smaller rods (≤9 mm) are near
    the noise floor; positions are estimated from fine radial scans.

    Background reference: two annular crowns centred on the phantom — inner at
    r≈37.5 mm and outer at r≈62.5 mm (0.75× and 1.25× the 50 mm rod radius), each
    8 mm wide — mirror the standard CTP515 approach and sample LDPE on both
    sides of the rod ring.

    References:
        `CP604 Manual, p.24 <CP604Manual.pdf#page=24>`_
    """

    attr_name = "ctp515"
    common_name = "Low Contrast (CTP730)"
    roi_dist_mm = 50   # all rods sit at r≈50 mm from centre
    num_slices = 20    # average ±20 slices for SNR improvement
    roi_settings = {
        # 1.0 % contrast group — positive contrast (Anchor: -90.00)
        "1pct_15mm": {"angle":  -90.00, "distance": roi_dist_mm, "radius": 6.0},
        "1pct_9mm":  {"angle": -108.00, "distance": roi_dist_mm, "radius": 3.5},
        "1pct_8mm":  {"angle": -124.29, "distance": roi_dist_mm, "radius": 3.0},
        "1pct_7mm":  {"angle": -138.86, "distance": roi_dist_mm, "radius": 2.5},
        "1pct_6mm":  {"angle": -151.71, "distance": roi_dist_mm, "radius": 2.0},
        "1pct_5mm":  {"angle": -162.86, "distance": roi_dist_mm, "radius": 1.5},
        "1pct_4mm":  {"angle": -172.29, "distance": roi_dist_mm, "radius": 1.2},
        "1pct_3mm":  {"angle": -180.00, "distance": roi_dist_mm, "radius": 1.0},
        "1pct_2mm":  {"angle":  174.00, "distance": roi_dist_mm, "radius": 0.8},

        # 0.5 % contrast group — positive contrast (Anchor: 150.00)
        "05pct_15mm": {"angle": 150.00, "distance": roi_dist_mm, "radius": 6.0},
        "05pct_9mm":  {"angle": 132.00, "distance": roi_dist_mm, "radius": 3.5},
        "05pct_8mm":  {"angle": 115.71, "distance": roi_dist_mm, "radius": 3.0},
        "05pct_7mm":  {"angle": 101.14, "distance": roi_dist_mm, "radius": 2.5},
        "05pct_6mm":  {"angle":  88.29, "distance": roi_dist_mm, "radius": 2.0},
        "05pct_5mm":  {"angle":  77.14, "distance": roi_dist_mm, "radius": 1.5},
        "05pct_4mm":  {"angle":  67.71, "distance": roi_dist_mm, "radius": 1.2},
        "05pct_3mm":  {"angle":  60.00, "distance": roi_dist_mm, "radius": 1.0},
        "05pct_2mm":  {"angle":  54.00, "distance": roi_dist_mm, "radius": 0.8},

        # 0.3 % contrast group — positive contrast (Anchor: 30.00)
        "03pct_15mm": {"angle":  30.00, "distance": roi_dist_mm, "radius": 6.0},
        "03pct_9mm":  {"angle":  12.00, "distance": roi_dist_mm, "radius": 3.5},
        "03pct_8mm":  {"angle":  -4.29, "distance": roi_dist_mm, "radius": 3.0},
        "03pct_7mm":  {"angle": -18.86, "distance": roi_dist_mm, "radius": 2.5},
        "03pct_6mm":  {"angle": -31.71, "distance": roi_dist_mm, "radius": 2.0},
        "03pct_5mm":  {"angle": -42.86, "distance": roi_dist_mm, "radius": 1.5},
        "03pct_4mm":  {"angle": -52.29, "distance": roi_dist_mm, "radius": 1.2},
        "03pct_3mm":  {"angle": -60.00, "distance": roi_dist_mm, "radius": 1.0},
        "03pct_2mm":  {"angle": -66.00, "distance": roi_dist_mm, "radius": 0.8},
}
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
        # Replace the global default (15) with the Rose-criterion value for CTP730.
        # Callers who pass an explicit threshold other than 15 keep their override.
        if cnr_threshold == 15:
            cnr_threshold = 1.0
        super().__init__(
            catphan,
            tolerance=tolerance,
            cnr_threshold=cnr_threshold,
            offset=offset,
            contrast_method=contrast_method,
            visibility_threshold=visibility_threshold,
            clear_borders=clear_borders,
        )

    # ------------------------------------------------------------------ helpers

    @property
    def rois_1pct(self) -> dict:
        """ROIs belonging to the 1.0 % contrast group."""
        return {k: v for k, v in self.rois.items() if k.startswith("1pct_")}

    @property
    def rois_05pct(self) -> dict:
        """ROIs belonging to the 0.5 % contrast group."""
        return {k: v for k, v in self.rois.items() if k.startswith("05pct_")}

    @property
    def rois_03pct(self) -> dict:
        """ROIs belonging to the 0.3 % contrast group."""
        return {k: v for k, v in self.rois.items() if k.startswith("03pct_")}

    def _bg_stats(self, _: str) -> tuple[float, float]:
        """Return (bg_median_HU, bg_std_HU) from the two annular background crowns."""
        inner = self.background_rois["inner"]
        outer = self.background_rois["outer"]
        bg_median = float(np.mean([inner.pixel_value, outer.pixel_value]))
        bg_std = float(np.mean([inner.std, outer.std]))
        return bg_median, bg_std

    def _scoring_table(self) -> list[dict]:
        """Per-rod metrics: contrast_pct, CNR, detectability."""
        rows: list[dict] = []
        for name, roi in self.rois.items():
            group, size_str = name.split("_", 1)
            diam_mm = int(size_str.replace("mm", ""))
            bg_median, bg_std = self._bg_stats(name)
            delta_hu    = float(roi.pixel_value) - bg_median
            contrast_pct = abs(delta_hu) / (abs(bg_median) + 1000.0) * 100.0
            cnr = abs(delta_hu) / bg_std if bg_std > 0 else 0.0
            rows.append({
                "group":        group,
                "diameter_mm":  diam_mm,
                "rod_hu":       float(roi.pixel_value),
                "bg_hu":        bg_median,
                "delta_hu":     delta_hu,
                "contrast_pct": contrast_pct,
                "cnr":          cnr,
                "detectability": contrast_pct * diam_mm,
            })
        return rows

    # ------------------------------------------------------------------ results

    _GROUP_LABELS: dict[str, str] = {
        "1pct":  "1.0 %",
        "05pct": "0.5 %",
        "03pct": "0.3 %",
    }

    def results(self) -> str:
        """Text summary table of low-contrast metrics per rod group."""
        rows = self._scoring_table()
        hdr = f"{'Diam':>6}  {'Rod HU':>8}  {'Bg HU':>8}  {'Delta':>7}  {'Ctr%':>7}  {'CNR':>6}  {'Detect':>7}"
        sep = "-" * len(hdr)
        lines: list[str] = [
            f"CTP730 Low-Contrast Results  (slice {self.slice_num}, "
            f"roll {self.catphan_roll:.2f} deg)",
            "=" * len(hdr),
        ]
        for gkey in ("1pct", "05pct", "03pct"):
            glabel = self._GROUP_LABELS[gkey]
            group_rows = sorted(
                [r for r in rows if r["group"] == gkey],
                key=lambda r: -r["diameter_mm"],
            )
            lines += [f"\n{glabel} contrast group:", hdr, sep]
            for r in group_rows:
                lines.append(
                    f"{r['diameter_mm']:>6}  {r['rod_hu']:>8.1f}  {r['bg_hu']:>8.1f}  "
                    f"{r['delta_hu']:>7.2f}  {r['contrast_pct']:>7.3f}  "
                    f"{r['cnr']:>6.2f}  {r['detectability']:>7.2f}"
                )
            visible = [r for r in group_rows if r["cnr"] >= self.cnr_threshold]
            if visible:
                sv = min(visible, key=lambda r: r["diameter_mm"])
                lines.append(
                    f"  => Smallest visible: {sv['diameter_mm']} mm  "
                    f"(CNR >= {self.cnr_threshold})"
                )
            else:
                lines.append(f"  => No rods detected  (CNR threshold {self.cnr_threshold})")
        return "\n".join(lines)

    def results_data(self) -> dict:
        """Dict export of low-contrast results."""
        return {
            "module": "CTP730",
            "slice_number": self.slice_num,
            "catphan_roll_deg": round(self.catphan_roll, 3),
            "rois": self._scoring_table(),
        }

    # ------------------------------------------------------------------ plots

    def plotly_analyzed_image(self) -> go.Figure:
        """Interactive Plotly figure: ROI overlay + contrast-detail curve.

        Returns a ``plotly.graph_objects.Figure`` with two panels:
        * Left (65 %): phantom slice with ROI circles, coloured by group.
        * Right (35 %): contrast-detail scatter (diameter vs contrast %).
        """
        from plotly.subplots import make_subplots

        rows   = self._scoring_table()
        arr    = self.image.array.astype(float)
        wl, ww = 50.0, 40.0

        _GROUP_COLORS = {"1pct": "#00ff00", "05pct": "#00ffff", "03pct": "#ffff00"}

        fig = make_subplots(
            rows=1, cols=2,
            column_widths=[0.65, 0.35],
            subplot_titles=[
                f"CTP730 ROI overlay — slice {self.slice_num}",
                "Contrast-detail curve",
            ],
        )

        # Phantom image
        fig.add_trace(
            go.Heatmap(
                z=arr,
                colorscale="gray",
                zmin=wl - ww / 2,
                zmax=wl + ww / 2,
                showscale=False,
                hoverinfo="skip",
            ),
            row=1, col=1,
        )

        # ROI circles as layout shapes
        shapes: list[dict] = []
        for name, roi in self.rois.items():
            group = name.split("_")[0]
            c = _GROUP_COLORS.get(group, "#ffffff")
            diam_label = name.split("_")[1]
            r = roi.radius
            shapes.append(
                dict(
                    type="circle",
                    xref="x", yref="y",
                    x0=roi.center.x - r, x1=roi.center.x + r,
                    y0=roi.center.y - r, y1=roi.center.y + r,
                    line=dict(color=c, width=1.5),
                    name=name,
                )
            )
            # diameter label annotation
            fig.add_annotation(
                x=roi.center.x, y=roi.center.y,
                text=diam_label.replace("mm", ""),
                showarrow=False,
                font=dict(size=7, color=c),
                xref="x", yref="y",
                row=1, col=1,
            )
        fig.update_layout(shapes=shapes)

        # Contrast-detail scatter traces
        for gkey in ("1pct", "05pct", "03pct"):
            glabel = self._GROUP_LABELS[gkey]
            c      = _GROUP_COLORS[gkey]
            gdata  = sorted(
                [r for r in rows if r["group"] == gkey],
                key=lambda r: r["diameter_mm"],
            )
            if not gdata:
                continue
            diams     = [r["diameter_mm"]  for r in gdata]
            contrasts = [r["contrast_pct"] for r in gdata]
            hover     = [
                (
                    f"<b>{glabel} — {r['diameter_mm']} mm</b><br>"
                    f"Contrast: {r['contrast_pct']:.3f} %<br>"
                    f"CNR: {r['cnr']:.2f}<br>"
                    f"Detectability: {r['detectability']:.2f}<br>"
                    f"Rod HU: {r['rod_hu']:.1f}  Bg HU: {r['bg_hu']:.1f}"
                )
                for r in gdata
            ]
            fig.add_trace(
                go.Scatter(
                    x=diams,
                    y=contrasts,
                    mode="lines+markers",
                    name=glabel,
                    line=dict(color=c, width=2),
                    marker=dict(size=8, color=c),
                    hovertext=hover,
                    hoverinfo="text",
                ),
                row=1, col=2,
            )

        fig.update_xaxes(title_text="Diameter (mm)", row=1, col=2)
        fig.update_yaxes(title_text="Contrast (%)", row=1, col=2)
        fig.update_layout(
            title=dict(
                text=(
                    f"CTP730 Low-Contrast Analysis — "
                    f"slice {self.slice_num}, "
                    f"roll {self.catphan_roll:.2f} °"
                ),
                font=dict(size=14),
            ),
            height=620,
            plot_bgcolor="#111111",
            paper_bgcolor="#1a1a1a",
            font=dict(color="#eeeeee"),
            legend=dict(
                yanchor="top", y=0.99,
                xanchor="right", x=0.99,
                bgcolor="rgba(30,30,30,0.7)",
            ),
        )
        return fig

    def save_analyzed_image(self, filename: str, **kwargs) -> None:
        """Save a static matplotlib figure: ROI overlay + contrast-detail chart."""
        rows = self._scoring_table()
        arr  = self.image.array.astype(float)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
        fig.patch.set_facecolor("#111111")
        for ax in (ax1, ax2):
            ax.set_facecolor("#1a1a1a")

        _GROUP_COLORS_MPL = {"1pct": "lime", "05pct": "cyan", "03pct": "yellow"}

        # Panel 1 — ROI overlay
        ax1.imshow(arr, cmap="gray", vmin=30.0, vmax=70.0, origin="upper")
        for name, roi in self.rois.items():
            group = name.split("_")[0]
            c = _GROUP_COLORS_MPL.get(group, "white")
            ax1.add_patch(
                plt.Circle((roi.center.x, roi.center.y), roi.radius,
                            color=c, fill=False, lw=1.2)
            )
            ax1.text(
                roi.center.x, roi.center.y,
                name.split("_")[1].replace("mm", ""),
                color=c, fontsize=5, ha="center", va="center",
            )
        ax1.set_title("CTP730 ROI overlay", color="white", fontsize=11)
        ax1.axis("off")

        # Panel 2 — contrast-detail
        for gkey in ("1pct", "05pct", "03pct"):
            glabel = self._GROUP_LABELS[gkey]
            c      = _GROUP_COLORS_MPL[gkey]
            gdata  = sorted(
                [r for r in rows if r["group"] == gkey],
                key=lambda r: r["diameter_mm"],
            )
            if gdata:
                ax2.plot(
                    [r["diameter_mm"]  for r in gdata],
                    [r["contrast_pct"] for r in gdata],
                    "o-", color=c, label=glabel, lw=2, markersize=6,
                )
        ax2.set_xlabel("Diameter (mm)", color="white")
        ax2.set_ylabel("Contrast (%)", color="white")
        ax2.set_title("Contrast-detail curve", color="white", fontsize=11)
        ax2.tick_params(colors="white")
        ax2.spines[:].set_color("#555555")
        ax2.grid(alpha=0.25, color="#444444")
        ax2.legend(facecolor="#2a2a2a", labelcolor="white")

        fig.suptitle(
            f"CTP730 — slice {self.slice_num}, roll {self.catphan_roll:.2f}°",
            color="white", fontsize=12,
        )
        fig.savefig(filename, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
        plt.close(fig)


# ---------------------------------------------------------------------------
# CatPhan 604
# ---------------------------------------------------------------------------

@capture_warnings
class CatPhan604(CatPhanBase):
    """CatPhan 604.

    Analyzes: HU linearity (CTP732), uniformity (CTP729), spatial resolution
    via bead MTF (CTP528CP604), and low contrast (CTP730).
    """

    _demo_url = "CatPhan604.zip"
    _model = "604"
    catphan_radius_mm = 101
    modules = {
        CTP732: {"offset": 0},
        CTP729: {"offset": 80},
        CTP528CP604: {"offset": 0},   # co-located with HU plug slice (Section 2)
        CTP730: {"offset": 40},       # Section 3 is at +40 mm from Section 2
    }

    @staticmethod
    def run_demo(show: bool = True):
        cbct = CatPhan604.from_demo_images()
        cbct.analyze()
        print(cbct.results())
        cbct.plot_analyzed_image(show)

    def refine_origin_slice(self, initial_slice_num: int) -> int:
        """Find the HU-linearity slice with the least wire-ramp angle.

        Scans ±5 slices around ``initial_slice_num``, measures the horizontal
        and vertical wire angles on each, and returns the slice with the
        smallest mean angle.  This is robust to the RM R1-4 jig being present.
        """
        angles = []
        ctp404, offset = self._get_module(CTP732, raise_empty=True)
        original_setup = ctp404._setup_rois
        ctp404._setup_rois = lambda x: x
        ctp = ctp404(
            self,
            offset=offset,
            clear_borders=self.clear_borders,
            hu_tolerance=0,
            scaling_tolerance=0,
            thickness_tolerance=0,
        )
        ctp404._setup_rois = original_setup
        for slice_num in range(initial_slice_num - 5, initial_slice_num + 5):
            sl = Slice(self, slice_num, clear_borders=self.clear_borders)
            troi = {}
            for name, setting in ctp.thickness_roi_settings.items():
                troi[name] = ThicknessROI.from_phantom_center(
                    sl.image,
                    setting["width_pixels"],
                    setting["height_pixels"],
                    setting["angle_corrected"],
                    setting["distance_pixels"],
                    sl.phan_center,
                )
            left_wire = troi["Left"].long_profile.center_idx
            right_wire = troi["Right"].long_profile.center_idx
            h_angle = abs(left_wire - right_wire)
            top_wire = troi["Top"].long_profile.center_idx
            bottom_wire = troi["Bottom"].long_profile.center_idx
            v_angle = abs(top_wire - bottom_wire)
            angle = (h_angle + v_angle) / 2
            angles.append(
                {
                    "slice": slice_num,
                    "angle": angle,
                    "left width": troi["Left"].long_profile.field_width_px,
                    "right width": troi["Right"].long_profile.field_width_px,
                    "left center": troi["Left"].long_profile.y_at_x(left_wire),
                    "right center": troi["Right"].long_profile.y_at_x(right_wire),
                    "left profile": troi["Left"].long_profile.values,
                    "right profile": troi["Right"].long_profile.values,
                }
            )

        median_width_l = np.median([a["left width"] for a in angles])
        median_width_r = np.median([a["right width"] for a in angles])
        median_width = (median_width_l + median_width_r) / 2
        median_left_pixel_val = np.median(
            np.concatenate([a["left profile"] for a in angles])
        )
        median_right_pixel_val = np.median(
            np.concatenate([a["right profile"] for a in angles])
        )
        median_pixel_val = (median_left_pixel_val + median_right_pixel_val) / 2
        max_left_pixel_val = np.max(
            np.concatenate([a["left profile"] for a in angles])
        )
        max_right_pixel_val = np.max(
            np.concatenate([a["right profile"] for a in angles])
        )
        max_pixel_val = (max_left_pixel_val + max_right_pixel_val) / 2

        for angle_set in angles.copy():
            if (
                angle_set["left width"] < median_width * 0.7
                or angle_set["right width"] < median_width * 0.7
            ):
                angles.remove(angle_set)
                continue
            fwxm_pixel = np.mean(
                (angle_set["left center"], angle_set["right center"])
            )
            delta_median = abs(median_pixel_val - fwxm_pixel)
            delta_max = abs(max_pixel_val - fwxm_pixel)
            if delta_median < delta_max:
                angles.remove(angle_set)

        m_slice_num = np.argsort([a["angle"] - self.catphan_roll for a in angles])
        return angles[m_slice_num[0]]["slice"]


# ---------------------------------------------------------------------------
# Backward-compatibility aliases (old pylinac / upstream names)
# ---------------------------------------------------------------------------

#: Old name for :class:`CTP732`.
CTP404CP604 = CTP732

#: Old name for :class:`CTP729`.
CTP486 = CTP729

#: CTP730 dispatches via the CTP515 ABC; re-export the old name.
CTP515 = CTP730
