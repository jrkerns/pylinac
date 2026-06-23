"""CatPhan 700 phantom module classes.

Physical modules:
  CTP682 — HU linearity (10 inserts including Lung #7112; replaces CTP404)
  CTP712 — Image uniformity (identical to CTP486)
  CTP714 — 30 lp/cm high-resolution rectangular ROIs (replaces CTP528)
  CTP515CP700 — Low contrast (same as CTP515CP600)
"""

from __future__ import annotations

from functools import cached_property

import numpy as np
from skimage.transform import EuclideanTransform

from ..core.mtf import MTF, PointSourceMTF
from ..core.roi import RectangleROI
from ..core.warnings import capture_warnings
from ..ct import (
    AIR,
    ACRYLIC,
    BONE_20,
    BONE_50,
    DELRIN,
    LDPE,
    LUNG_7112,
    POLY,
    PMP,
    TEFLON,
    WATER,
    CatPhanBase,
    CTP528 as _CTP528Base,
    SpatialResolutionROI,
)
from . import CP504, CP600


# ---------------------------------------------------------------------------
# CTP682 — HU linearity (CatPhan 700)
# ---------------------------------------------------------------------------

class CTP682(CP504.CTP404):
    """HU linearity module for CatPhan 700 (CTP682).

    Contains 10 inserts including a Lung #7112 plug (−868 HU) and an optional
    water vial.  Only Top and Bottom wire-ramp ROIs (no Left/Right — those
    positions hold bead ramps).  No geometry nodes (empty dict).

    Angles use the 180-x formula compared with the standard CTP404 layout.
    """

    roi_dist_mm = 58.7
    roi_radius_mm = 5
    roi_settings = {
        "Air": {
            "value": AIR,
            "angle": 180 - -90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "PMP": {
            "value": PMP,
            "angle": 180 - -120,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Lung #7112": {
            "value": LUNG_7112,
            "angle": 180 - -165,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Delrin": {
            "value": DELRIN,
            "angle": 180 - 165,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Poly": {
            "value": POLY,
            "angle": 180 - 120,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Teflon": {
            "value": TEFLON,
            "angle": 180 - 90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Bone 20%": {
            "value": BONE_20,
            "angle": 180 - 60,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "LDPE": {
            "value": LDPE,
            "angle": 180 - 15,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Bone 50%": {
            "value": BONE_50,
            "angle": 180 - -15,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Acrylic": {
            "value": ACRYLIC,
            "angle": 180 - -60,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Vial": {
            "value": WATER,
            "angle": 180 - -135,
            "distance": 28,
            "radius": roi_radius_mm,
        },
    }
    background_roi_settings = {
        "1": {
            "angle": -37.5,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "2": {
            "angle": -142.5,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "3": {
            "angle": 142.5,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "4": {
            "angle": 37.5,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
    }
    # CatPhan 700 has wire ramps only on top and bottom.
    # Left and right positions are bead ramps (different evaluation, not yet implemented).
    thickness_roi_height = 40
    thickness_roi_width = 10
    thickness_roi_distance_mm = 40
    thickness_roi_settings = {
        "Bottom": {
            "angle": 90,
            "width": thickness_roi_height,
            "height": thickness_roi_width,
            "distance": thickness_roi_distance_mm,
        },
        "Top": {
            "angle": -90,
            "width": thickness_roi_height,
            "height": thickness_roi_width,
            "distance": thickness_roi_distance_mm,
        },
    }
    geometry_roi_settings = {}  # no geometry nodes on CTP682

    def _setup_rois(self) -> None:
        """Remove the optional water vial ROI if no water is detected."""
        super()._setup_rois()
        if self.rois["Vial"].pixel_value < -500:  # closer to air than water
            self.rois.pop("Vial")


# ---------------------------------------------------------------------------
# CTP712 — Image uniformity (identical to CTP486)
# ---------------------------------------------------------------------------

class CTP712(CP504.CTP486):
    """Image uniformity module for CatPhan 700 (CTP712).

    Functionally and geometrically identical to CTP486.
    """

    pass


# ---------------------------------------------------------------------------
# CTP714 — 30 lp/cm high-resolution module (rectangular ROIs)
# ---------------------------------------------------------------------------

class CTP714(_CTP528Base):
    """Spatial resolution module for CatPhan 700 (CTP714).

    MTF is derived from the non-tilted tungsten wire embedded in the insert.
    The wire runs parallel to the z-axis, so its cross-section occupies the
    same pixel position in every slice; ``num_slices`` adjacent slices are
    averaged (via ``combine_method="mean"``) to improve SNR before the PSF is
    extracted.  The wire cross-section is treated as a point source and fed to
    :class:`~pylinac.core.mtf.PointSourceMTF` (no sinc correction needed for a
    50 µm wire at CT resolution).

    The rectangular line-pair gauge ROIs are retained for visual reference
    (``_setup_rois`` / ``plot_rois``), but are no longer used for the MTF value.
    """

    rois: dict[str, RectangleROI]
    attr_name: str = "ctp528"
    common_name: str = "Spatial Resolution"
    combine_method: str = "mean"
    num_slices: int = 3
    start_angle = None  # not applicable; kept for serialisation compatibility

    # --- Wire PSF settings ---
    PATCH_RADIUS_MM: float = 3.0
    BACKGROUND_INNER_MM: float = 3.5
    BACKGROUND_OUTER_MM: float = 6.0
    # Wire sits inside the line-pair gauge ring; exclude outer region.
    WIRE_SEARCH_MAX_RADIUS_MM: float = 40.0

    roi_settings = {
        # Regions listed clockwise from top-left on CT image
        "region 1": {
            "lp/mm": 0.1,
            "radial_distance": 50,
            "transversal_distance": -7,
            "rotation": -90,
            "width": 3,
            "height": 11,
        },
        "region 2": {
            "lp/mm": 0.2,
            "radial_distance": 50,
            "transversal_distance": 11,
            "rotation": -90,
            "width": 3,
            "height": 11,
        },
        "region 3": {
            "lp/mm": 0.3,
            "radial_distance": 50,
            "transversal_distance": -5.5,
            "rotation": -45,
            "width": 3,
            "height": 10,
        },
        "region 4": {
            "lp/mm": 0.4,
            "radial_distance": 50,
            "transversal_distance": 9.5,
            "rotation": -45,
            "width": 3,
            "height": 8.5,
        },
        "region 5": {
            "lp/mm": 0.5,
            "radial_distance": 50,
            "transversal_distance": -9,
            "rotation": 0,
            "width": 3,
            "height": 8,
        },
        "region 6": {
            "lp/mm": 0.6,
            "radial_distance": 50,
            "transversal_distance": 2,
            "rotation": 0,
            "width": 3,
            "height": 7,
        },
        "region 7": {
            "lp/mm": 0.7,
            "radial_distance": 50,
            "transversal_distance": 12,
            "rotation": 0,
            "width": 3,
            "height": 6,
        },
        "region 8": {
            "lp/mm": 0.8,
            "radial_distance": 50,
            "transversal_distance": -10.5,
            "rotation": 45,
            "width": 3,
            "height": 4,
        },
    }

    def _setup_rois(self) -> None:
        tform_catphan_global = EuclideanTransform(
            rotation=np.deg2rad(self.catphan_roll),
            translation=[self.phan_center.x, self.phan_center.y],
        )
        for name, setting in self.roi_settings.items():
            tform_roi_catphan = EuclideanTransform(
                translation=[
                    setting["radial_distance_pixels"],
                    setting["transversal_distance_pixels"],
                ]
            ) + EuclideanTransform(rotation=np.deg2rad(setting["rotation"]))
            tform_roi_global = tform_roi_catphan + tform_catphan_global
            self.rois[name] = SpatialResolutionROI(
                array=self.image.array,
                width=setting["width_pixels"],
                height=setting["height_pixels"],
                center=tform_roi_global.translation,
                rotation=np.rad2deg(tform_roi_global.rotation),
            )

    @cached_property
    def wire_center(self) -> tuple[float, float]:
        """Sub-pixel (x, y) centre of the tungsten wire in the combined image.

        Searches within ``WIRE_SEARCH_MAX_RADIUS_MM`` of the phantom centre to
        avoid the line-pair gauge patterns at the periphery.
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
        search[r_mm > self.WIRE_SEARCH_MAX_RADIUS_MM] = -np.inf

        if np.all(np.isinf(search)):
            raise ValueError(
                "Could not locate the wire: search region fully masked. "
                "Check WIRE_SEARCH_MAX_RADIUS_MM or phantom alignment."
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
        return cx, cy

    @cached_property
    def mtf(self) -> PointSourceMTF:
        """Wire cross-section PSF MTF (non-tilted wire, parallel to z)."""
        arr = self.image.array.astype(float)
        cx, cy = self.wire_center

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

        return PointSourceMTF(patch, self.mm_per_pixel, source_diameter_mm=0.0)


# ---------------------------------------------------------------------------
# CTP515CP700 — Low contrast (identical to CTP515CP600)
# ---------------------------------------------------------------------------

class CTP515CP700(CP600.CTP515CP600):
    """Low contrast module for CatPhan 700 (CTP515CP700).

    Identical to the CatPhan 600 layout.
    """

    pass


# ---------------------------------------------------------------------------
# CatPhan 700
# ---------------------------------------------------------------------------

@capture_warnings
class CatPhan700(CatPhanBase):
    """CatPhan 700.

    Analyzes: HU linearity (CTP682), low contrast (CTP515), uniformity
    (CTP712), and spatial resolution (CTP714).
    """

    _model = "700"
    catphan_radius_mm = 101
    modules = {
        CTP682: {"offset": 0},
        CTP515CP700: {"offset": -80},
        CTP712: {"offset": -160},
        CTP714: {"offset": -40},
    }


# ---------------------------------------------------------------------------
# Backward-compatibility aliases (old pylinac / upstream names)
# ---------------------------------------------------------------------------

#: Old name for :class:`CTP682`.
CTP404CP700 = CTP682

#: Old name for :class:`CTP714`.
CTP528CP700 = CTP714

#: Old name for :class:`CTP712`.
CTP486 = CTP712
