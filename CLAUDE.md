# avalanche-release-areas — Claude Code context

## What this repo is

Standalone package for physically-motivated avalanche release zone delineation.
Inputs: pre-computed SNOWPACK cluster features (CSVs + rasters).
Outputs: release polygons, scenario summaries, comparison figures.

No SNOWPACK, no zarr, no NWP pipeline — just this package and its data.

## Package layout

```
src/release_areas/
  arrest_indices.py    — Meloche et al. (2025) crack-arrest scaling laws (pure numpy)
  release_geometry.py  — BFS flood-fill propagation from trigger cluster
  snowpack_features.py — per-cluster feature helpers (slab, WL, profile)
  generate_scenarios.py — main CLI entry point (filter chain → BFS → output)
  plot_release.py      — comparison figure helpers
  scenario_writer.py   — writes GeoJSON + depth raster per scenario
  config.py            — all tunable constants (single source of truth)
tests/
  test_arrest_indices.py    — 42 unit tests for arrest_indices.py
  test_snowpack_features.py — 12 tests for element geometry helpers
  test_layered_slab.py      — 74 tests for layered sigma_t/E + ligament
  test_terrain_direction.py — 33 tests for aspect, downslope walk, directional θ
  test_theta_estimators.py  — 14 tests for the knn / plane_fit θ estimators
docs/
  release_area_methods.md — full derivation, calibration, Jan 18 application
data/little_prof/      — Little Professor path, Jan 18 2026 event (CC BY 4.0)
```

## Included data (Little Professor path, CDOT I-70)

| Path | Description |
|------|-------------|
| `features/all_start_zone_features_2026-01-18.csv` | Per-cluster WL/slab features — **reference**, current defaults |
| `features/meloche_features_all_2026-01-18.csv` | Meloche Π₁, A_ca, θ per cluster — **reference**, current defaults |
| `features/*_2026-01-18_v1.csv` | The pre-2026-10-09 baseline pair, kept for reproducing older numbers. Generated with plain element means, the undercounted `wl_thickness`, `project_fit` E, `SCH2004_EXP = 2.0` and `THETA_ESTIMATOR = 'knn'`. |
| `spatial/cluster_map.tif` | Cluster ID raster (1 m, EPSG:6342) |
| `dem_1m.tif` | 1 m DEM (EPSG:6342 + NAVD88) |
| `boundaries/start_zone.kml` | Topographic start zone boundary |
| `boundaries/domain.kml` | Full simulation domain |
| `boundaries/avalanche_release_area_20260118.geojson` | Observed crown, Jan 18 2026 (EPSG:6342, 6 934 m²) |
| `boundaries/20260118_avalanche_boundaries.geojson` | Source mapping in CRS84 (lon/lat), pre-reprojection |

## Install

```bash
pip install -e .          # runtime
pip install -e ".[dev]"   # + pytest
```

Or with uv (preferred):

```bash
uv pip install -e ".[dev]"
```

## Run the Jan 18 2026 scenario (from repo root)

```bash
python -m release_areas.generate_scenarios \
    --features-csv  data/little_prof/features/all_start_zone_features_2026-01-18.csv \
    --meloche-csv   data/little_prof/features/meloche_features_all_2026-01-18.csv \
    --cluster-map   data/little_prof/spatial/cluster_map.tif \
    --dem           data/little_prof/dem_1m.tif \
    --start-zone    data/little_prof/boundaries/start_zone.kml \
    --release-poly  data/little_prof/boundaries/avalanche_release_area_20260118.geojson \
    --out-dir       outputs/little_prof
```

Expected output (verified 2026-10-09, reference CSVs):
- Start zone clusters with features: 1608
- Filter chain, one count per stage: 490 (tau_g/slope/Sk38) → 406 (slab
  thickness) → 203 (elevation ≥ P50) → 102 (Pi1 ≥ median 45.81)
- Top-3 triggers: [2859, 5656, 6191] — Sk38 all 0.72, i.e. `N_TOP_TRIGGERS = 3`
  cuts *inside* a three-way tie, so the selection is decided by tie ordering,
  not by a margin
- Best IoU vs observed crown: **0.669** at the default `--max-clusters 500`
  (1 of 3 polygons is cap-bound), **0.677** at `--max-clusters 2000`
- At `--max-clusters 2000`: mean IoU 0.562, mean area/observed 0.63

`N_TOP_TRIGGERS` dropped 5 → 3 on 2026-10-09; five overlapping polygons made
the comparison figure unreadable. Per-scenario geometry is **unchanged** — the
three are bit-identical to the first three of the five-trigger run — but
aggregates are not, because the dropped scenarios 5817 and 348 were the two
largest (area ratios 1.03 and 1.20). **The mean area/observed 0.82 that matches
the 0.83 mask-limited ceiling, and which the methods §6 2×2 uses as its scale
check, is a five-trigger figure**: at three triggers it is 0.63. Reproduce the
2×2 with `--n-triggers 5`. Previous five-trigger output for reference:
[2859, 5656, 6191, 5817, 348], 3 of 5 cap-bound, mean IoU 0.574.

The reference CSVs are regenerated with the current defaults by
`examples/regenerate_jan18_reference_csvs.py` (needs avachain's interpreter for
zarr). **Do not delete the `_v1` pair**: that script reads it to pin the cluster
set and `group` labels, and it cannot be regenerated — it came from arithmetic
three commits back (0 of 3540 rows satisfy the post-fix D_wl identity).

Pointing both `--*-csv` flags at the `_v1` pair under **current** code gives
493 → 408 → 204 → 102, triggers [2859, 5656, 2858, 348, 1068], best IoU
**0.626** at the default cap and **0.669** at `--max-clusters 2000`, mean
area/observed 0.70 / 0.77. The default-cap figure was 0.639 before the
cross-slope θ sampling fix of 2026-10-09; it moved because those polygons are
cap-bound, so a changed `d_lat` alters which clusters fill the budget. The
`--max-clusters 2000` figure is unchanged. Old CSVs do **not** mean old
numbers — the geometry fixes live in the code, not the data.

Note `generate_scenarios` reads the CSVs and never calls `profile_features`, so
changing `THETA_ESTIMATOR` or any feature-generation default has **no effect on
a scenario run until the CSVs are regenerated**.

## Element weighting in profile_features

SNOWPACK's `z` is exactly `height_above_ground - HS`, so `np.diff` of sorted
`z` with the ground prepended at `-HS` recovers each element's own thickness
and those thicknesses sum to `HS` (verified to machine precision against the
zarr `height` variable; `element_thickness()`).

Two consequences, both fixed 2026-10-07 in this repo and in avachain's
`snowpack_analysis.py`:

- Bulk slab means are **thickness-weighted**, not element means. Element
  thicknesses span 3–38 mm within one profile, so a plain mean over-weights
  thin layers. Measured effect on Jan 18: ρ ×1.013, E ×1.033, τg ×1.004.
- A layer's thickness is the **sum of its elements' dz**, not the span of their
  `z`, which silently omits the basal element. `wl_thickness` was undercounting
  D_wl by a median 23% (p95 49%). Since the basal WL starts at the ground,
  `D_wl == hs - slab_thickness` exactly — useful for patching old CSVs.
  `crust_thickness` had the same defect, and also counted non-crust elements
  sandwiched between two crusts.

Λ ∝ √D_wl, so Λ rose ×1.13 and A_ca fell ×0.86. Measured end-to-end on Jan 18
with everything else held byte-identical: best IoU **0.667 → 0.636** at
`--max-clusters 2000`, mean area/observed 0.813 → 0.782. The fit got *worse* —
an inflated A_ca was partly compensating for the under-covered crown. That is
not a reason to revert it and **not** a licence to re-tune anything else.

**Vindicated 2026-10-09.** That call was right, and for the reason given. Once
θ was also corrected to the L_ss scale (which pushes A_ca the *other* way,
×2.9 against this fix's ×0.86 and the E relation's ×0.69), the fully-corrected
configuration became the best of all four: best IoU 0.677 and mean
area/observed 0.82. Each correction alone looks wrong; together they are right.
Keep this in mind before judging any single correction by Jan 18 IoU — see the
2×2 table in methods §6.

`wl_shear_strength` (τp) is deliberately still a plain element mean: weighting
it is a modelling choice (a crack runs in the weakest sublayer, so the minimum
may be the right reduction), and it cascades into θ and the trigger ranking.

## Material relations are selectable — check which is active

Every relation below is a *named option*; nothing is hardcoded any more, and
the defaults changed on 2026-10-07. All live in `arrest_indices.py`.

| What | Default | Options | Source of the default |
|------|---------|---------|-----------------------|
| `slab_modulus` | `vanherwijnen2016` | + `project_fit` | van Herwijnen et al. (2016) *J. Glaciol.* 62(236) Eq. 8: E = 0.93 ρ^2.8 Pa |
| `slab_tensile_strength` | (only one) | — | **this project** — σt = (ρ/300)^1.4 · 5 kPa, no published source |
| `k_ic` | `schweizer2004` | + `kirchner2000`; `borstad2013` **gated** | Schweizer, Michot & Kirchner (2004) *Ann. Glaciol.* 38 Eq. 8: K_Ic = 350 Pa·m · (ρ/917)² / √d_max |

- **`project_fit` is ours, not van Herwijnen's.** It is a hand-fit through the
  *range* that paper reports and sits 1.6–2.2× below their actual regression.
  `C_FIT` was calibrated at E = 4 MPa, ρ = 250 (Meloche Table 1 fixes ρ = 250
  and sweeps E as a constant — **there is no E(ρ) in Meloche**). At ρ = 250
  `project_fit` gives 2.54 MPa, `vanherwijnen2016` gives 4.82 MPa, so the new
  default is the one consistent with the calibration point.
- **`kirchner2000` is a LOWER BOUND on K_Ic** — apparent toughness from small
  cantilever beams, so it carries the small-specimen size effect.
- **The ρ exponent in `schweizer2004` is the paper's printed 1.9 since
  2026-10-09** (`SCH2004_EXP = SCH2004_EXP_PAPER`). Eq. 8 carries 1.9; the
  abstract says "about 2", which is where the old 2.0 default came from.
  Verified against the PDF 2026-10-08. 2.0 is 11% low at ρ = 300 (1.18 vs
  1.32 kPa·m^0.5) and 17% low at ρ = 150; `SCH2004_EXP_ROUND` reproduces it.
  Changed on fidelity grounds — nothing was ever measured to justify the
  rounding — and note it **cannot** move any release polygon: `k_ic()` reaches
  only the ligament columns, which the BFS never reads, so this was a *null*
  test against Jan 18 IoU and had to be decided from the source.
- **`schweizer2004` is only valid to 300 kg m⁻³**, below our median slab
  density (~330). `k_ic()` clamps ρ into the fitted range *for the K_Ic
  evaluation only* (never elsewhere) and `n_rho_clamped` counts it: on Jan 18
  that is a median 46 of 71 layers per cell, and over the full space-time run a
  median 46 of 81 with 94.3% of slabs carrying ≥1 clamped layer. K_Ic is
  effectively saturated over most of the slab. `clamp=False` restores
  NaN-outside-range.
- **`borstad2013` raises `NotImplementedError`** — bronze OA at Wiley only, no
  repository copy, so the regression has never been sourced.
- `d_max = DMAX_FACTOR · grain_size`, **both in metres**. Our zarr
  `grain_size` is already in m (median 5.9e-4); raw SNOWPACK `.pro` is in mm.

### Ligament (edge-crack) bound — hypothesis, `USE_LIGAMENT_BOUND = True`

The weakest-link controlling layer is treated as a crack of length a = its own
thickness in a slab of thickness h; K_Ic and σt come from the **intact** layers
(thickness-weighted over the slab minus that layer). El Haddad short-crack
correction is the default:

    sigma_c = K_Ic / (F(a/h) sqrt(pi (a + a0))),  a0 = (K_Ic/(F sigma_t_lig))^2/pi

`F` is the **same** `F(a/h)` in both places, so `a → 0` recovers `sigma_t_lig`
exactly and no iteration is needed. `self_consistent_a0=True` instead solves
`a0` at `F(a0/h)` iteratively (kept as an option; its `a → 0` limit is only
`σt·F(a0/h)/F(0)`). `model='lefm'` gives the uncorrected form, which overshoots
`sigma_t_lig` on thin controlling layers. A bound is withheld (NaN) when
`a0 ≥ h` — the intrinsic flaw exceeds the slab and the LEFM geometry is void.

**Not validated.** No K_Ic for snow slabs in tension has been verified against
a source we hold; the ligament bound has never been fed into the BFS or an IoU.

**Measured 2026-10-08, full space-time run** — all 6 565 clusters × all 501
six-hourly steps = 3 289 065 profiles, 2 169 737 (66.0%) with a resolved slab.
Medians: σ_c/σt_lig_intact **0.816** (`schweizer2004` d=1), **0.706** (d=2),
**0.595** (`kirchner2000`); σt_mean/σt(bulk) 0.875; σt_wl/σt_mean 0.650;
controlling a/h 0.0118. A_ca(`vanherwijnen2016`)/A_ca(`project_fit`) = **0.692**
— and that ratio is the *exact* identity √(E_pf/E_vh), because A_ca depends on E
only through Λ ∝ √E (L_t, τg, θ, σt are E-free; checked against `evaluate()` to
rtol 1e-12). The `a0 ≥ h` guard **never fired** in any variant, so it is
untested by data rather than confirmed by it. Full tables in methods §4.

**Gotcha for any space-time run:** the zarr `location` axis is
59 280 = 6 565 unique clusters repeated **9× at stride 6 565** (verified
byte-identical). Dedupe with `.isel(location=slice(0, 6565))`; counts computed
without it are inflated 9-fold. Chunking is `(25, 126, 85)`, so loading one
timestep costs the same as loading all 126 in its time chunk — iterate over
location blocks with all times loaded, not over timesteps.

## θ is lag-dependent — `THETA_ESTIMATOR` (default `plane_fit` since 2026-10-09)

θ is a property of the **lag it is measured at**, not of the snowpack. Jan 18
start-zone clusters are 3.0 m across and the k = 6 neighbourhood `knn` uses
averages **3.3 m**, i.e. inside the 0.5–10 m band Meloche et al. treat as local
noise in Appendix E. The τp variogram gives

    theta(d) = 3.87 + 72/d        (Pa/m, d in m)

cross-checked twice: the trend term implies a ramp of 6.07 Pa/m against 4.98
from an independent planar fit, and the noise term 72 Pa against 78 Pa
predicted from the 69 Pa nugget. So **at 3.3 m, θ is 15% slope-scale trend and
85% local noise**; at Meloche's L_ss = 20 m it is 7.5 Pa/m.

| option | what it is |
|---|---|
| `plane_fit` **(default)** | local least-squares plane through τp within `THETA_FIT_RADIUS_M` (= `L_SS`), θ = \|∇τp\|. Estimates the ramp Meloche's θ represents, at the scale it was calibrated over. |
| `knn` | mean \|Δτp\|/d over `K_NEIGHBORS` nearest centroids. Reproduces every pre-2026-10-09 CSV and published number. |

Measured on the shipped Jan 18 τp field (1 308 start-zone clusters): θ median
25.19 (`knn`) vs 7.28 (`plane_fit`) Pa/m, per-cluster ratio median **3.30**. As
A_ca ∝ 1/θ, `plane_fit` multiplies A_ca by ~3.3× — matching the variogram's
independent 3.4× and giving a **calibration-free** explanation of the
"A_ca is systematically too small" pattern in the v2 comparison. It implicates
the θ *estimator*, not the slab elastic relation and not any Meloche constant.

**Why `plane_fit` is the default.** The argument is the calibration scale, not
the fit: Meloche's θ is the gradient of a linear ramp over L_ss, so the
estimator's neighbourhood should be L_ss. `THETA_FIT_RADIUS_M = L_SS` is
asserted in the tests for exactly that reason.

Two things to keep in view, neither of which the switch resolves:
- **97.2%** of start-zone clusters fall below `THETA_VALID_MIN` (20 Pa/m) under
  `plane_fit`, against 36.2% under `knn`. At the trend scale almost the whole
  start zone is outside the paper's calibrated regime — read literally, the
  Jan 18 weak layer has no slope-scale ramp steep enough to arrest within L_ss.
  This is the strongest open objection to the whole θ-based arrest route and is
  **not** evidence against `plane_fit`; `knn` only looked compliant because it
  was measuring noise.
- θ is NaN for **15.8%** of start-zone clusters, which lack
  `THETA_FIT_MIN_NEIGHBOURS` within the radius. End to end this costs almost
  nothing (overall NaN rate 31.8% → 32.2%), because those clusters were already
  NaN for other reasons.

`plane_fit` also emits `theta_grad_east` / `theta_grad_north` (NaN under
`knn`), so a θ split is available from the estimator itself if a future caller
wants one.

**The Gaume width ratio no longer depends on the CSV θ.** `estimate_cross_slope_width`
measures *both* halves itself via `_sector_theta` — θ_along about the fall
line, θ_cross about the cross axis — from the same τp field in the same
separation band. Before 2026-10-09 the numerator was the CSV's isotropic k-NN
mean at a ~3.3 m lag while the denominator was directional at 5–50 m; given
θ(d) ≈ 3.87 + 72/d that mismatch alone inflated the ratio several-fold and
pushed it into `GAUME_ASPECT_CAP`. Fixing it cut the Jan 18 Gaume widths from
55/44/50/155/38 m to 39/34/34/112/20 m (trigger 348 is no longer
cap-saturated). The CSV θ survives only as a fallback numerator when no
cluster lies in the along-slope sector.

This change is **inert on the default Jan 18 configuration** — IoU is
bit-identical — because `d_lat = min(gaume_width, A_ca · 0.712)` and the mode
III speed cap binds on all 5 triggers, so the Gaume width is not what sets
lateral extent. It matters only with `USE_MODE3_SPEED_CAP = False` or where the
cap does not bind. Measured on the `_v1` pair with the cap off, `d_lat` becomes
the Gaume width (39/34/34/112/20 m) and best IoU is 0.646, against 0.669 with
the cap on. Worth remembering that the Gaume path is effectively dormant on
this event.

**Do not pick the estimator by Jan 18 IoU.** Changing θ moves Π₁, hence the
trigger ranking and the whole filter chain. Argue it from the propagation scale
the scaling law was calibrated at (L_ss). Full derivation in methods §5.

**Adding a τp Gaussian random field is *not* the indicated next step** — a
reading of Appendix E that earlier notes got wrong. Appendix E is a sensitivity
analysis whose stated result is that local noise "mainly affects the crack
speed"; the paper's arrest-relevant heterogeneity *is* the linear ramp. Our
measured σ_local is 69 Pa against Appendix E's 500 Pa, and superimposing that
field would inject 55–169 Pa/m of spurious θ (2–7× the real median), shrinking
A_ca further in the wrong direction.

## Terrain direction — two fixed bugs, and the lesson

`compute_slope_aspect` is **correct**: verified to 0.0000° against a map-frame
gradient over real easting/northing arrays and a least-squares plane fit. The
start zone is genuinely **ESE 121.5°** (observed crown axis 119.0°).

Two bugs confused map direction with raster indexing, both fixed 2026-10-09:

- `find_stauchwall` stepped `dr = sign(-cos(asp)), dc = sign(sin(asp))`.
  `np.sign` is ±1 for any nonzero component, so only the four bearings
  45/135/225/315 were reachable — a pure diagonal on **100.00%** of start-zone
  pixels, median **+13.9°** clockwise error, and the walk was a straight ray
  that never followed terrain. Now carries sub-pixel position.
- `estimate_cross_slope_width` picked "abeam" neighbours by raster row/col,
  which assumes a north–south fall line. Selected neighbours sat a mean **30°**
  off the true cross axis (one within 1.2° of the fall line), so `theta_cross`
  absorbed the along-slope gradient. Now samples in the map frame within
  ±`THETA_CROSS_SECTOR_DEG` of the cross axis.

**The lesson generalises:** this package converts between pixel `(row, col)`
and map `(east, north)` constantly, and **row increases southward**. Any new
code that derives a direction must be checked against the map frame, never
assumed from array axes. `tests/test_terrain_direction.py` covers both helpers
with planar-DEM fixtures at known bearings; neither had coverage before.

Side effect worth knowing: with `theta_cross` sampled correctly the Gaume width
saturates at `GAUME_ASPECT_CAP` more often, so the **mode III speed cap now
binds on all 5 Jan 18 triggers rather than 2** and alone sets lateral extent.

Residual polygon lean is the **start-zone KML itself** — its own principal axis
is 136.4° SE — i.e. the same boundary that caps IoU at 0.830.

## Run tests

```bash
pytest tests/ -v
```

All 42 tests in `test_arrest_indices.py`, the 12 in
`test_snowpack_features.py`, the 74 in `test_layered_slab.py`, the 33 in
`test_terrain_direction.py` and the 14 in `test_theta_estimators.py` must pass
(175 total). Do not modify tolerances to make
failing tests pass — fix the underlying formula or inputs.

## Key physical constraints

- **Never tune parameters to the Jan 18 2026 event.** It is the only validation event.
  Calibration targets come from Meloche et al. (2025) Table 1 / Fig 8.
  This still holds after the 2026-10-09 default changes, which are *not* an
  exception to it. `THETA_ESTIMATOR` and `SCH2004_EXP` were each chosen between
  two **named, sourced options** on source grounds — the L_ss calibration scale
  and the paper's printed exponent — not by fitting a free parameter. Jan 18 was
  run afterwards as a check, and the record (methods §6) states explicitly that
  the IoU margin is ~1% and too thin to have decided anything. If you are ever
  tempted to pick a value *because* Jan 18 likes it, note what the 2×2 showed:
  single-factor IoU selection would have rejected all three corrections that
  are, together, the best configuration. On n = 1 the metric is not just weak,
  it is actively misleading.
- `R_FIT` and `C` in `arrest_indices.py` are Meloche calibration constants — do not change.
- `arrest_indices.evaluate()` is the primary public API — keep its signature stable.
- `config.py` is the single source of truth for every tunable constant. Do not
  reintroduce duplicate literals in other modules.
- Release depth comes from `slab_thickness` (failure plane → surface), never
  from `hs` (total snowpack height) — `hs` would imply a full-depth release.
- `config.MAX_BFS_CLUSTERS` is a *safety cap*, not physics. When it binds, the
  run prints a `[CAP-BOUND]` warning and the polygon is not a physical result.
- `config.USE_MELOCHE_ARREST` is **False**: the BFS arrests on Λ/thickness
  discontinuity heuristics and the distance caps, not on the Meloche
  per-direction A_ca criterion. Enabling it changes every polygon.
- The three absolute floors — `TAU_G_ABS_FLOOR` (350 Pa),
  `MIN_PROPAGATION_SLAB` (0.5 m), `MIN_PROPAGATION_LAMBDA` (0.1 m) — are
  **inert on the Jan 18 data**: measured 0 arrests each. 350 Pa is at P0.0 of
  start-zone τg (min 553 Pa). Do not spend effort calibrating them; the
  Λ/thickness discontinuity factors (~32% of rejections) are the real target.
- **IoU is capped at 0.830 on Jan 18**: 17% of the observed crown (1 180 of
  6 934 m²) lies outside the start-zone KML, and the BFS hard-rejects clusters
  outside the mask. 29.4% of all arrests are that mask. Report IoU against
  this ceiling, not against 1.0.

## Filter chain in generate_scenarios.py

Order matters; do not reorder:
1. Spatial: intersect `features_df.index` with `cluster_map[start_zone_mask]`
2. tau_g ≥ 40 Pa, slope ≥ (stauchwall_deg + 2°), Sk38 < 1.0
3. Slab thickness 0.5–max_slab_thickness m (default 2.0; use 1.5 for skier triggers)
4. Elevation ≥ P50 of remaining candidates
5. Pi1_elastic ≥ median of remaining pool

## CRS handling

All rasters and the cluster map are **EPSG:6342** (NAD83(2011) / UTM zone 13N;
the DEM is a compound CRS with NAVD88 height). This is *not* EPSG:32613 (WGS84
UTM 13N) — earlier revisions of this file said 32613, which was wrong.

`config.RASTER_EPSG` is the single source of truth; `load_kml_mask`,
`geojson_to_mask` and `load_observed_polygon(s)` all default to it.

GeoJSON inputs may be CRS84/WGS84 or UTM — CRS is auto-detected via regex
`re.search(r'EPSG:+(\d+)', ...)` to handle both `EPSG:6342` and `EPSG::6342`.
A `urn:ogc:def:crs:OGC:1.3:CRS84` declaration does not match the regex and
correctly falls through to the EPSG:4326 default.

## Mode III Lambda TODO

`elastic_length_cross()` in `arrest_indices.py` currently uses
`Λ_III = √(G h D_wl / G_wl)` giving `Λ_III/Λ_II ≈ 0.59` (ν=0.3).
Update when Johan Gaume's mode III (antiplane) formula is received.

A `--mode3-scale` flag used to be threaded through `make_release_polygon_2d` →
`propagate_crack` but was never read; removed rather than left as a no-op.

**The code is prepared for the swap.** `propagate_crack._qualifies()` splits
arrest evaluation into along-slope vs cross-slope branches, and
`release_geometry.directional_lambda(props, is_cross)` is the single place the
two diverge. Cross-slope has its own continuity thresholds
(`LAMBDA_CROSS_DROP_FACTOR` / `..._RISE_FACTOR`, currently equal to the
along-slope pair). `compute_meloche_features` writes a `Lambda_cross` column.

To land the real formula: (1) update `arrest_indices.elastic_length_cross()`,
(2) set `config.USE_MODE3_LAMBDA = True`, (3) regenerate the meloche CSV — the
shipped Jan 18 one predates the `Lambda_cross` column. Nothing else should need
to change; if it does, the directional split has drifted.

**Caveat:** a constant Λ_III/Λ_II ratio changes nothing, because the continuity
gate is scale-invariant (relative change) and `MIN_PROPAGATION_LAMBDA` = 0.1 m
is never approached. Only a Λ_III with different *spatial structure* will move
the flanks. See methods §3.

## Mode III crack-speed cap (implemented)

Mode III cracks cap at `c_s` where mode II runs supershear at ~1.6 `c_s`
(Broberg 1989). A slower crack builds slab tension faster per unit advance, so
the cross-slope distance to first slab fracture is shorter:
`L_dyn,III/L_dyn,II = 0.712` at ν=0.3. E and ρ cancel — the ratio depends only
on ν and the two speed ratios (`arrest_indices.mode3_length_ratio`).

Unlike Λ_III this is an **absolute distance cap**, so it does change results.
`propagate_crack` composes it with the Gaume width as a minimum:
`d_lat = min(gaume_width, A_ca · 0.712) · size_factor`.

Flags: `USE_MODE3_SPEED_CAP` (default True), `MODE2_SPEED_RATIO` = 1.6,
`MODE3_SPEED_RATIO` = 1.0. Turning it off is bit-for-bit identical to pre-cap.

**Known asymmetry — read before changing this.** Only the *restrictive* half of
mode III is coded. Because G < E′, an equal-strength slab has a *larger* energy
cap cross-slope (`slab_energy_cap_cross`, G_slab,III = τ_flank²h/2G), which
pushes the other way — but it needs a flank strength τ_flank that has no
parameterisation, so it is absent. On Jan 18 the cap binds on 2 of 5 triggers
and makes the aggregate fit slightly worse (mean IoU 0.567 → 0.561, area ratio
0.92 → 0.78). That is consistent with implementing one side of a two-sided
effect; it is left on because Broberg's limit holds independently of this event.
Do not tune it against Jan 18. See methods §3.

Lateral *distance* extent is separately set by `estimate_cross_slope_width()`
(Gaume 2015 / θ ratio, capped at `GAUME_ASPECT_CAP`).

## Commit style

- Short imperative subject line
- End with: `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`
