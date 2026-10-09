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
  test_layered_slab.py      — 72 tests for layered sigma_t/E + ligament
docs/
  release_area_methods.md — full derivation, calibration, Jan 18 application
data/little_prof/      — Little Professor path, Jan 18 2026 event (CC BY 4.0)
```

## Included data (Little Professor path, CDOT I-70)

| Path | Description |
|------|-------------|
| `features/all_start_zone_features_2026-01-18.csv` | Per-cluster WL/slab features |
| `features/meloche_features_all_2026-01-18.csv` | Meloche Π₁, A_ca, θ per cluster |
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

Expected output (verified):
- Start zone clusters with features: 1608
- Filter chain, one count per stage: 493 (tau_g/slope/Sk38) → 408 (slab
  thickness) → 204 (elevation ≥ P50) → 102 (Pi1 ≥ median)
- Top-5 triggers: [2859, 5656, 2858, 348, 1068]
- Best IoU vs observed crown: 0.639 at the default `--max-clusters 500`
  (3 of 5 polygons are cap-bound), 0.669 at `--max-clusters 2000`
  (0.667 before the `find_stauchwall` direction fix of 2026-10-09; the
  default-cap figure is unchanged because those polygons are cap-bound)

**The shipped CSVs predate the element-weighting fix** (see below), so these
numbers are the pre-fix baseline. They still reproduce exactly, because
`generate_scenarios` reads the CSVs and never calls `profile_features`.

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
- **The ρ exponent we use in `schweizer2004` is 2.0; the paper's printed Eq. 8
  carries 1.9** (its abstract says "about 2"). Verified against the PDF
  2026-10-08. 2.0 is 11% low at ρ = 300 and 17% low at ρ = 150;
  `SCH2004_EXP_PAPER` restores 1.9. The default stays 2.0 because that is what
  was specified and tested — this is our approximation, not Schweizer et al.'s.
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

## Run tests

```bash
pytest tests/ -v
```

All 42 tests in `test_arrest_indices.py`, the 12 in
`test_snowpack_features.py` and the 72 in `test_layered_slab.py` must pass
(126 total). Do not modify tolerances to make
failing tests pass — fix the underlying formula or inputs.

## Key physical constraints

- **Never tune parameters to the Jan 18 2026 event.** It is the only validation event.
  Calibration targets come from Meloche et al. (2025) Table 1 / Fig 8.
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
