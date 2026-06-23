# Catphan QA — Pylinac Audit & Implementation Notes

This document captures all findings from reading the Catphan product manuals and auditing the
installed Pylinac source (`pylinac.ct`). It exists so that a new conversation or a new project
tackling the Pylinac rework can pick up immediately without re-deriving this context.

Status legend: ✅ Done · ⚠️ Partial / known issue · ❌ Not started

---

## 1. Available phantom manuals

Located at `src/qa/ct_analysis/phantoms/`:

| File | Phantom |
|---|---|
| `CTP604Manual.pdf` | Catphan 604 — primary target, fully audited below |
| `CTP504Manual.pdf` | Catphan 504 |
| `CTP503Manual.pdf` | Catphan 503 |
| `CTP500-600.pdf` | Catphan 500 / 600 series |
| `CTP700Manual.pdf` | Catphan 700 |

Only the 604 has been read and audited in detail. The others are present for when
multi-phantom support is added.

---

## 2. Catphan 604 — Physical layout

The phantom is 20 cm in diameter. The reference section is **Section 2**; all distances are
measured from its center along the z-axis.

| z offset | Section | Module ID | Contents |
|---|---|---|---|
| 0 mm | Section 2 | *(no CTP number)* | HU sensitometry plugs · MTF bead · MTF wire · 23° wire ramps · 3 mm pixel-size holes |
| −40 mm | Section 1 | *(no CTP number)* | 15 lp/cm visual line-pair gauge |
| +40 mm | Section 3 | **CTP730** | Low-contrast rods |
| +80 mm | Section 4 | *(no CTP number)* | Solid uniformity field |

> **Sign convention note (CatPhan 604, this scanner):** In the local implementation the
> physical z-direction matches the manual sign. `CTP730` at Section 3 (+40 mm manual) uses
> `offset: +40` in `CatPhan604.modules`. The upstream pylinac default used `offset: -40`,
> which placed the slice on the wrong side; this was corrected during Phase 3 work (2026-06-23).
> Other models (504, 503, 600, 700) have not been audited for the same issue.

---

## 3. Section 2 — Sensitometry / HU accuracy ✅ Done

The 9 sensitometry inserts (all in Section 2) and their nominal HU ranges at 120 kVp per
The Phantom Laboratory manuals. Ranges are **non-symmetric and per-material** — not the
uniform ±40 HU that Pylinac applies.

| Material | Formula / Notes | Specific Gravity | HU range |
|---|---|---|---|
| Air | — | 0.00 | −1046 to −986 |
| PMP (polymethylpentene) | C₆H₁₂(CH₂) | 0.83 | −220 to −172 |
| LDPE (low-density polyethylene) | C₂H₄ | 0.92 | −121 to −87 |
| Polystyrene | C₈H₈ | 1.03 | −65 to −29 |
| Acrylic | C₅H₈O₂ | 1.18 | 92 to 137 |
| Bone 20% | Ca/C/N/O/P mix | 1.14 | 211 to 263 |
| Delrin | Proprietary (acetal) | 1.42 | 344 to 387 |
| Bone 50% | Ca/C/N/O/P mix | 1.40 | 667 to 783 |
| Teflon | CF₂ | 2.16 | 941 to 1060 |

All other Catphan models share the same ranges for common materials; model-specific
inserts (Vial/Water, Lung #7112) are included in `MANUAL_HU_RANGES` in `orientation.py`.

### What was wrong in Pylinac and what was fixed

**Problem 1 — ROI angles 180° off (8 of 9 inserts).**
Pylinac's `CTP404CP604` has hardcoded angles that were mirror-images of the actual physical
insert positions in this scanner's phantom. Only Air (at −90°/+90°) was correct; all others
were negated. Fixed in `CTP404CP604Fixed` in `orientation.py` with angles confirmed by
live HU cross-check.

**Problem 2 — Uniform ±40 HU tolerance for all materials.**
Pylinac passes a single scalar `hu_tolerance` to every `HUDiskROI`. The manufacturer
specifies different, non-symmetric ranges per material (see table above). Fixed by
`MANUAL_HU_RANGES` dict in `orientation.py`; `_extract_results()` in `catphan_analysis.py`
uses these for pass/fail and display, bypassing `roi.passed`.

**Problem 3 — Phantom orientation not detected.**
Pylinac assumes a fixed phantom orientation. If the phantom is placed the other way
(flipped along z), every non-CTP404 module lands on the wrong slice. Fixed by
`OrientationMixin` in `orientation.py` using a 4-filter content-fingerprinting pipeline
that identifies the CTP404 slice (max HU std) and determines flip direction by counting
valid phantom slices on each side.

### Known HU anomalies observed in live scans

See `HU_ANOMALY_OBSERVATIONS.md` for full analysis. In brief:
- **Air ~−984 HU** (2 HU above upper limit): likely scanner air calibration baseline drift.
- **Teflon and 50% Bone**: σ 3–4× higher than other materials + values near lower limit: beam hardening artifact (expected physics).

---

## 4. Section 2 — MTF (bead and wire) ✅ Done

### 4a. What the phantom has

Section 2 contains **two** point-source MTF targets:

1. **Tungsten carbide bead** — diameter 0.18 mm. Subpixel at typical CT resolution;
   no size-correction needed in most cases (some software optionally corrects for it).
   Recommended by the manual for thicker slices.

2. **Angled tungsten wire** — diameter 50 μm, tilted 5° to the z-axis.
   Because of the 5° angle, the wire crosses voxel boundaries as you step through slices.
   On any single slice the wire may occupy 1 or more voxels depending on sub-voxel position.
   **Correct method:** use several adjacent slices to build a high-SNR PSF via sub-pixel
   oversampling, then derive MTF from the averaged PSF.
   Recommended by the manual for thin slices.

### 4b. MTF algorithm (implemented)

The pipeline is shared across all phantom models via a class hierarchy in `pylinac/core/mtf.py`:

```
PointSourceMTF (base)
│  Accepts a 2D PSF patch + pixel_size_mm + optional source_diameter_mm.
│  Pipeline: radial average → symmetric mirror → Hanning window → FFT → normalise.
│  If source_diameter_mm > 0: divide by sinc(π·d·f) to correct finite source size.
│
├── BeadMTF (subclass)
│   SOURCE_DIAMETER_MM = 0.18 mm  → sinc correction applied automatically.
│   Used by: CTP528CP504 (503/504), CTP764 (600), CTP764-via-CTP528CP504 (600 free inheritance).
│
└── SlantedWireMTF (subclass)
    Accepts image_stack + wire_center + pixel_size_mm + slice_spacing_mm + tilt_angle_deg.
    Builds PSF by registering WIRE_SLICES patches at tan(θ)×dz/pixel_size lateral drift per slice,
    averaging into a high-SNR 2D PSF, then feeds into PointSourceMTF pipeline.
    source_diameter_mm = 0.0 (50 µm wire — no sinc correction needed at CT resolution).
    Used by: CTP528CP604 (604 slanted wire, ±7 slices, TILT_ANGLE_DEG=5°).
```

MTF strategy per model:

| Model | Module | Source | Class |
|---|---|---|---|
| 503 / 504 | `CTP528CP504` | 0.18 mm bead (Section 2) | `BeadMTF` |
| 600 | `CTP764` | 0.18 mm bead (inherited) | `BeadMTF` (via `CTP528CP504`) |
| 604 | `CTP528CP604` | 50 µm wire, 5° tilt | `SlantedWireMTF` |
| 700 | `CTP714` | 50 µm wire, parallel | `PointSourceMTF` (source_diameter=0) |

### 4c. Key implementation details

**Bead detection** (`CTP528CP504.bead_center`):
- Search within `BEAD_SEARCH_MAX_RADIUS_MM=40.0` mm to exclude the line-pair gauge ring at ~47 mm.
- Find peak pixel, refine with weighted centroid in a small window.

**Wire detection** (`CTP528CP604.wire_center`):
- Search excludes HU plug annulus at 58.7 mm (`HU_PLUG_DIST_MM`) ± 12 mm.
- `preprocess(catphan)` hook stores `catphan.dicom_stack` as `self._dicom_stack` (the only
  way to access individual slices, since `CatPhanModule.__init__` does not keep a catphan ref).

**Sinc correction**: `MTF(f) /= sinc(π·d·f)` where d = bead diameter. Applied for bead (0.18 mm),
not for wire (50 µm — correction is < 0.1 % at Nyquist for typical CT pixel sizes).

**Line-pair `circle_profile`**: retained in `CTP528CP504` as a `circle_profile` property for
visual display of the line-pair gauge, but is no longer used to compute the `mtf` value.

**Testing**: CatPhan503, CatPhan504, CatPhan604 all smoke-tested successfully (2026-06-23).
CatPhan600 and CatPhan700 have no local DICOM data to test against.

---

## 5. Section 3 — Low contrast (CTP730) ✅ Done

### 5a. What the phantom has (CTP730)

Rod diameters: **2, 3, 4, 5, 6, 7, 8, 9, 15 mm**
Nominal contrast levels: **0.3%, 0.5%, 1.0%**
Rod length: 40 mm

The three contrast groups are physically separated in the image. Rods of the same contrast
are cast from a single batch, so intra-group contrast uniformity is guaranteed.

Scoring metrics used by the manual:
```
contrast% = |rod_HU − bg_HU| / (bg_HU + 1000) × 100
CNR       = |rod_HU − bg_HU| / bg_std
detect    = contrast% × diameter_mm    (Rose model figure of merit)
```

### 5b. What is implemented (`pylinac/phantoms/CP604.py::CTP730`)

`CTP730` is a full implementation of the low-contrast module:
- All 27 ROIs defined: 9 rods × 3 contrast groups (1.0 %, 0.5 %, 0.3 %)
- ROIs grouped via `rois_1pct`, `rois_05pct`, `rois_03pct` convenience properties
- Inherits from `CTP515` for circular-ROI infrastructure (`LowContrastDiskROI`)
- `num_slices = 20` → 41-slice mean (±20 from centre, ~25 mm, ~3× noise reduction)
- `roi_dist_mm = 58` (measured from a 41-slice profile scan; original spec was 55 mm)
- Offset `+40` (correct for this scanner's z-direction; fixed from the wrong upstream `−40`)
- Background ROIs at 0.75× and 1.25× rod distance (inner / outer pair per rod)

**Per-rod scoring** — `_scoring_table()` returns a list of dicts with:
  - `rod_hu`, `bg_hu`, `delta_hu`, `contrast_pct`, `cnr`, `detectability`

**Output methods:**
| Method | Returns |
|---|---|
| `results()` | Formatted text table, one block per contrast group |
| `results_data()` | Dict with `module`, `slice_number`, `catphan_roll_deg`, `rois` list |
| `plotly_analyzed_image()` | Interactive Plotly figure: ROI overlay (65 %) + contrast-detail scatter (35 %) |
| `save_analyzed_image(filename)` | Static matplotlib figure — same two-panel layout |

### 5c. Angular position calibration

Group positions were calibrated from a 41-slice multi-scan profile at r = 58 mm:

| Group | Detected 15 mm rod | Implemented range |
|---|---|---|
| 1.0 % (most visible) | 106.7° | 107° → 171° (8° steps) |
| 0.5 % | −110° (= 250°) | −110° → −174° (8° steps) |
| 0.3 % | 20.6° | 21° → −43° (8° steps) |

**Known limitation**: only the 15 mm rod anchor per group was reliably detected.
The within-group diameter order (which angle → 15 mm, 9 mm, …, 2 mm) is assumed
monotonic from the 15 mm anchor and uses 8° uniform spacing. A precise calibration
requires comparing against the manufacturer's p.24 diagram or a higher-SNR scan.

### 5d. Pending

- **CTP515 alignment** ⏳ — The CTP515 low-contrast module (CatPhan 504) uses a
  different scoring schema (no per-rod detectability). Once the CTP730 output is
  audited in the private app, align CTP515 to the same `results()` / `results_data()`
  / `plotly_analyzed_image()` / `save_analyzed_image()` interface.

---

## 6. Section 4 — Uniformity ✅ Done

Solid cast phantom material; CT number within 20 HU of water (typical range 5–18 HU).

Measurements:
- Mean HU at centre ROI and 4 peripheral ROIs
- Noise = standard deviation within each ROI
- Integral non-uniformity = (CTmax − CTmin) / (CTmax + CTmin) along a cross-section profile

**Pylinac `CTP486`** implements exactly this. The module offset in Pylinac is `−80` (Section 4 is
at +80 mm in the manual; sign flip as described in §2). **This module is correct and usable as-is.**

---

## 7. Pylinac audit summary

| Pylinac class | Pylinac offset | Physical target | Verdict | Status |
|---|---|---|---|---|
| `CTP732` (`CTP404CP604`) | 0 | Section 2 — HU plugs | ✅ Correct geometry; angles & tolerances fixed in `CTP404CP604Fixed` | ✅ Done |
| `CTP729` (`CTP486`) | −80 | Section 4 — Uniformity | ✅ Correct, use as-is | ✅ Done |
| `CTP528CP604` | 0 | Section 2 — slanted-wire MTF | ✅ Replaced: `SlantedWireMTF` (±7 slices, 5° tilt, PSF stack) | ✅ Done |
| `CTP528CP504` / `CTP764` | varies | Section 2 — bead MTF (503/504/600) | ✅ Replaced: `BeadMTF` (0.18 mm bead, sinc-corrected PSF) | ✅ Done |
| `CTP730` | +40 | Section 3 — CTP730 low contrast | ✅ Full: 27 ROIs, contrast%/CNR/detectability, `results()`/`results_data()`/plotly/matplotlib | ✅ Done |

---

## 8. Implementation roadmap

### Phase 1 — Scaffolding & working results ✅ Done
- ✅ CT Analysis registered in catalog and workflow
- ✅ Zip and multi-DICOM upload
- ✅ Phantom model selection (604, 504, 503, 600, 700)
- ✅ HU accuracy: `CTP404CP604Fixed` with corrected angles + `MANUAL_HU_RANGES` for per-material pass/fail
- ✅ Orientation detection: `OrientationMixin` (4-filter fingerprinting, count-based flip detection)
- ✅ Phantom roll: custom Air-insert method in `_compute_phantom_roll()`
- ✅ Geometry metrics: pixel size, slice thickness, phantom roll, CTP404 origin slice
- ✅ Uniformity & noise: `CTP486`, 5-position ROI table + bar chart
- ✅ Analyzed image tabs (side view, CTP404, CTP486, CTP528, CTP730)
- ✅ CSV export

### Phase 2 — Custom MTF ✅ Done (2026-06-23)
- ✅ `PointSourceMTF` base class in `pylinac/core/mtf.py` (shared PSF → Hanning → FFT pipeline)
- ✅ `BeadMTF(PointSourceMTF)` — 0.18 mm bead, sinc correction, backward-compatible API
- ✅ `SlantedWireMTF(PointSourceMTF)` — multi-slice PSF stack, tilt-angle registration
- ✅ `CTP528CP504` → `BeadMTF` (bead detection within 40 mm radius; 503/504 + 600 via inheritance)
- ✅ `CTP528CP604` → `SlantedWireMTF` (±7 slices, `preprocess()` hook captures dicom_stack)
- ✅ `CTP714` → `PointSourceMTF` (non-tilted wire, single-slice PSF, no sinc)
- Smoke-tested: CatPhan503, CatPhan504, CatPhan604 all pass; 600/700 untested (no local DICOM)

### Phase 3 — CTP730 low contrast ✅ Done (2026-06-23)
- ✅ `CTP730` class in `pylinac/phantoms/CP604.py` with 27 ROIs (9 rods × 3 contrast groups)
- ✅ ROIs grouped by contrast level (`rois_1pct`, `rois_05pct`, `rois_03pct`)
- ✅ Angular positions calibrated from 41-slice profile scan at r=58 mm; `roi_dist_mm=58`, offset=+40
- ✅ `num_slices=20` (41-slice mean, ~3× noise reduction; std 10.4→3.3 HU)
- ✅ Per-rod metrics: contrast%, CNR, detectability (`_scoring_table()`)
- ✅ Background ROI pair (inner + outer) per rod via `_bg_stats()`
- ✅ `results()` — formatted text table per contrast group
- ✅ `results_data()` — dict export (module, slice, roll, rois list)
- ✅ `plotly_analyzed_image()` — interactive ROI overlay + contrast-detail scatter
- ✅ `save_analyzed_image()` — static matplotlib two-panel figure
- ✅ Smoke-tested end-to-end on 408-slice CatPhan 604 DICOM (0.625 mm/slice), 2026-06-23
- ⏳ Within-group intra-rod angle ordering not verified vs. physical manual diagram (p.24)
- ⏳ `CTP515` alignment (CatPhan 504) — same `results()`/`plotly` interface not yet applied

---

## 9. File structure

The phantom classes have been extracted from `pylinac/ct.py` into a dedicated sub-package.
All original names are re-exported from `ct.py` for backward compatibility.

```
pylinac/
├── ct.py                        # CatPhanBase, dispatch ABCs, re-exports from phantoms/
├── core/
│   └── mtf.py                   # PointSourceMTF, BeadMTF, SlantedWireMTF, MTF (legacy)
└── phantoms/
    ├── __init__.py              # imports CP504/CP600/CP604/CP700 sub-modules
    ├── CP504.py                 # CTP404, CTP486, CTP528CP504/503, CTP515
    ├── CP600.py                 # CTP763, CTP764, CTP515CP600, CTP591, CatPhan600
    ├── CP604.py                 # CTP732, CTP729, CTP528CP604, CTP730, CatPhan604
    ├── CP700.py                 # CTP682, CTP712, CTP714, CTP515CP700, CatPhan700
    └── CATPHAN_PYLINAC_NOTES.md # this file — audit + roadmap

Private app (separate repo):
src/qa/ct_analysis/
├── __init__.py
├── catphan_analysis.py          # main display + analysis entry point
├── orientation.py               # OrientationMixin, CTP404CP604Fixed, MANUAL_HU_RANGES
├── HU_ANOMALY_OBSERVATIONS.md  # observed HU anomalies and physics explanations
└── phantoms/
    ├── CTP604Manual.pdf
    ├── CTP504Manual.pdf
    ├── CTP503Manual.pdf
    ├── CTP500-600.pdf
    └── CTP700Manual.pdf
```

---

## 10. Key class aliases (backward compatibility)

All original upstream names remain importable from `pylinac.ct`:

| New physical name | Old upstream name | File |
|---|---|---|
| `CTP404` | — | CP504.py |
| `CTP528CP504` | — | CP504.py |
| `CTP486` | `CTP486` | CP504.py |
| `CTP515` | `CTP515` | CP504.py |
| `CTP732` | `CTP404CP604` | CP604.py |
| `CTP729` | `CTP486` (in 604 ns) | CP604.py |
| `CTP730` | `CTP515` (in 604 ns) | CP604.py |
| `CTP763` | `CTP404CP600` | CP600.py |
| `CTP764` | `CTP528CP600` | CP600.py |
| `CTP682` | `CTP404CP700` | CP700.py |
| `CTP714` | `CTP528CP700` | CP700.py |
| `CTP712` | `CTP486` (in 700 ns) | CP700.py |
