"""CatPhan 600 phantom module classes.

Physical modules:
  CTP763 — HU linearity, geometry, slice thickness (urethane; replaces CTP404)
  CTP764 — 21 lp/cm high-resolution gauge (replaces CTP528)
  CTP515CP600 — Low contrast (same targets as CTP515, rotated 180°)
  CTP591 — Bead geometry / slice thickness (new for 600)

CTP591 contains three bead-ramp pairs (2 coarse × 1mm z-spacing, 1 fine ×
0.25mm z-spacing) plus individual beads and a 50 µm tungsten wire for MTF.
Physical offset from CTP763: +32.5 mm (pylinac convention: offset=-32).
"""

from __future__ import annotations

import numpy as np
from plotly import graph_objects as go
from scipy.ndimage import maximum_filter

from ..core.warnings import capture_warnings
from ..ct import (
    AIR,
    ACRYLIC,
    DELRIN,
    LDPE,
    POLY,
    PMP,
    TEFLON,
    WATER,
    CatPhanBase,
    CatPhanModule,
    Point,
)
from . import CP504


# ---------------------------------------------------------------------------
# CTP763 — HU linearity module (urethane body, replaces CTP404)
# ---------------------------------------------------------------------------

class CTP763(CP504.CTP404):
    """HU linearity module for CatPhan 600 (CTP763, cast in urethane).

    Angles are mirrored vs CTP404/504.  The top slot is an optional water vial;
    it is removed from analysis if air is detected there instead of water.
    """

    roi_dist_mm = 58.7
    roi_radius_mm = 5
    roi_settings = {
        "Air": {
            "value": AIR,
            "angle": 90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "PMP": {
            "value": PMP,
            "angle": 60,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "LDPE": {
            "value": LDPE,
            "angle": 0,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Poly": {
            "value": POLY,
            "angle": -60,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Acrylic": {
            "value": ACRYLIC,
            "angle": -120,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Delrin": {
            "value": DELRIN,
            "angle": -180,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Teflon": {
            "value": TEFLON,
            "angle": 120,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Vial": {
            "value": WATER,
            "angle": -90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm - 1,  # vial sits inside the ROI slot; needs clearance
        },
    }

    def _setup_rois(self) -> None:
        """Remove the optional water vial ROI if no water is detected."""
        super()._setup_rois()
        if self.rois["Vial"].pixel_value < -500:  # closer to air than water
            self.rois.pop("Vial")


# ---------------------------------------------------------------------------
# CTP764 — 21 lp/cm high-resolution gauge (line-pair MTF)
# ---------------------------------------------------------------------------

class CTP764(CP504.CTP528CP504):
    """21 lp/cm high-resolution spatial resolution module for CatPhan 600 (CTP764).

    Functionally identical to CTP528 but cast in urethane.  Different boundary
    positions and direction compared with the 504 module.

    MTF is computed from the embedded point-source bead via
    :class:`~pylinac.core.mtf.BeadMTF` (inherited from
    :class:`~pylinac.phantoms.CP504.CTP528CP504`).
    """

    start_angle = np.pi - 0.1
    ccw = False
    boundaries = (0, 0.127, 0.195, 0.255, 0.304, 0.354, 0.405, 0.453, 0.496)
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


# ---------------------------------------------------------------------------
# CTP515CP600 — Low contrast (same hardware as CTP515, rotated 180°)
# ---------------------------------------------------------------------------

class CTP515CP600(CP504.CTP515):
    """Low contrast module for CatPhan 600.

    Same physical targets as CTP515 but the angular positions are rotated 180°
    relative to the 504 layout.
    """

    roi_angles = [
        -87.4 + 180,
        -69.1 + 180,
        -52.7 + 180,
        -38.5 + 180,
        -25.1 + 180,
        -12.9 + 180,
    ]
    roi_dist_mm = 50
    roi_radius_mm = [6, 3.5, 3, 2.5, 2, 1.5]
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


# ---------------------------------------------------------------------------
# CTP591 — Bead geometry module (NEW)
# ---------------------------------------------------------------------------

class CTP591(CatPhanModule):
    """Bead geometry module for CatPhan 600 (CTP591).

    Physical layout (CTP500/600 manual, section 19–25):

    * 2 coarse bead ramps — 0.28 mm beads, 1.0 mm z-spacing
    * 1 fine bead ramp   — 0.18 mm beads, 0.25 mm z-spacing
    * 2 individual beads (0.28 mm and 0.18 mm)
    * 1 × 50 µm tungsten wire (for PSF/MTF)

    Physical offset: +32.5 mm from CTP763.
    pylinac convention: offset = -32 (images increasing in z).

    **Implemented**
    ~~~~~~~~~~~~~~~
    * Bright-spot detection — locates up to ``N_BRIGHTEST`` candidate bead/wire
      positions in the axial slice using non-maximum suppression (local max
      filter), reporting sub-pixel (x, y) coordinates and HU values.

    **Not yet implemented**
    ~~~~~~~~~~~~~~~~~~~~~~~
    * Slice-thickness from bead ramps — requires building a z-profile for each
      ramp column across the DICOM stack and counting peaks at FWHM::

          slice_thickness = N_peaks_at_FWHM × bead_z_spacing

      The two ramp spacings (1.0 mm coarse, 0.25 mm fine) give independent
      estimates; the fine ramp is the more precise one.
    * Wire PSF/MTF — requires isolating the tungsten wire response and computing
      the ESF/LSF along the wire axis.  The wire sits at ``WIRE_DIST_MM`` from
      the phantom centre and can be distinguished from beads by its elongated
      shape across slices.
    """

    attr_name = "ctp591"
    common_name = "Bead Geometry (CTP591)"
    num_slices = 1

    # Bead ramp geometry from the CTP500/600 manual
    FINE_BEAD_SPACING_MM: float = 0.25
    COARSE_BEAD_SPACING_MM: float = 1.0
    WIRE_DIST_MM: float = 60.0  # tungsten wire distance from phantom centre

    # Detection parameters
    N_BRIGHTEST: int = 10         # max bright spots to report
    MIN_DISTANCE_MM: float = 2.0  # minimum separation between detected spots

    def _setup_rois(self) -> None:
        """Detect the N brightest pixels as candidate bead/wire locations."""
        arr = self.image.array.astype(float)
        min_dist_px = max(2, int(self.MIN_DISTANCE_MM / self.mm_per_pixel))
        local_max = arr == maximum_filter(arr, size=min_dist_px * 2 + 1)
        ys, xs = np.where(local_max)
        vals = arr[ys, xs]
        order = np.argsort(vals)[::-1][: self.N_BRIGHTEST]
        self.bead_centers: list[Point] = [
            Point(x=float(xs[i]), y=float(ys[i])) for i in order
        ]
        self.bead_values: list[float] = [float(vals[i]) for i in order]

    def _convert_units_in_settings(self) -> None:
        pass

    def plotly_rois(self, fig: go.Figure) -> None:
        for pt, val in zip(self.bead_centers, self.bead_values):
            r_px = self.MIN_DISTANCE_MM / self.mm_per_pixel
            theta = np.linspace(0, 2 * np.pi, 36)
            fig.add_scatter(
                x=(pt.x + r_px * np.cos(theta)).tolist(),
                y=(pt.y + r_px * np.sin(theta)).tolist(),
                mode="lines",
                line=dict(color="cyan", width=1),
                name=f"Bead {val:.0f} HU",
                showlegend=False,
            )

    def plot_rois(self, axis) -> None:
        import matplotlib.pyplot as plt

        for pt in self.bead_centers:
            r_px = self.MIN_DISTANCE_MM / self.mm_per_pixel
            circle = plt.Circle(
                (pt.x, pt.y), radius=r_px, color="cyan", fill=False, linewidth=1
            )
            axis.add_patch(circle)


# ---------------------------------------------------------------------------
# CatPhan 600
# ---------------------------------------------------------------------------

@capture_warnings
class CatPhan600(CatPhanBase):
    """CatPhan 600.

    Analyzes: HU linearity (CTP763), uniformity (CTP486), low contrast
    (CTP515), spatial resolution (CTP764), and bead geometry (CTP591).
    """

    _demo_url = "CatPhan600.zip"
    _model = "600"
    catphan_radius_mm = 101
    modules = {
        CTP763: {"offset": 0},
        CP504.CTP486: {"offset": -160},
        CTP515CP600: {"offset": -110},
        CTP764: {"offset": -70},
        CTP591: {"offset": -32},
    }

    @staticmethod
    def run_demo(show: bool = True):
        cbct = CatPhan600.from_demo_images()
        cbct.analyze()
        print(cbct.results())
        cbct.plot_analyzed_image(show)

    def analyze(self, *args, **kwargs) -> None:
        """Analyze the CatPhan 600; extends base analysis to include CTP591."""
        super().analyze(*args, **kwargs)
        if self._has_module(CTP591):
            ctp591_cls, offset = self._get_module(CTP591)
            self.ctp591 = ctp591_cls(
                self, offset=offset, tolerance=None, clear_borders=self.clear_borders
            )

    def find_phantom_roll(self, func=None) -> float:
        """With the CatPhan 600, we have to consider that the top air ROI
        has a water vial in it (see pg 12 of the manual). If so, the top air ROI won't be detected.
        Rather, the default algorithm will find the bottom air ROI and teflon to the left.
        It may also find the top air ROI if the water vial isn't there.
        We use the below lambda to select the bottom air and teflon ROIs consistently.
        These two ROIs are at 75 degrees from cardinal. We thus offset the default outcome by 75.

        HOWEVER, for direct density scans, the Teflon ROI might not register because of the reduced
        HU. We make a best guess depending on the detected roll. If it's ~75 degrees,
        we have caught the bottom Air and Teflon. If it's near zero, we have
        caught the top and bottom Air ROIs.
        """
        angle = super().find_phantom_roll(lambda x: -x.centroid[0])
        if abs(angle) < 10:
            return angle
        return angle + 75


# ---------------------------------------------------------------------------
# Backward-compatibility aliases (old pylinac / upstream names)
# ---------------------------------------------------------------------------

#: Old name for :class:`CTP763`.
CTP404CP600 = CTP763

#: Old name for :class:`CTP764`.
CTP528CP600 = CTP764
