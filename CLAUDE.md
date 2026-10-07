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
  test_arrest_indices.py — 34 unit tests for arrest_indices.py (all must pass)
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
  (3 of 5 polygons are cap-bound), 0.667 at `--max-clusters 2000`

## Run tests

```bash
pytest tests/ -v
```

All 34 tests in `test_arrest_indices.py` must pass. Do not modify tolerances to make
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
- End with: `Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>`
