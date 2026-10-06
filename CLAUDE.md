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
  config.py            — SlopeConfig dataclass
tests/
  test_arrest_indices.py — 34 unit tests for arrest_indices.py (all must pass)
docs/
  crack_arrest_indices_reference.md — full derivation, calibration, Jan 18 application
data/little_prof/      — Little Professor path, Jan 18 2026 event (CC BY 4.0)
```

## Included data (Little Professor path, CDOT I-70)

| Path | Description |
|------|-------------|
| `features/all_start_zone_features_2026-01-18.csv` | Per-cluster WL/slab features |
| `features/meloche_features_all_2026-01-18.csv` | Meloche Π₁, A_ca, θ per cluster |
| `spatial/cluster_map.tif` | Cluster ID raster (1 m, EPSG:32613) |
| `dem_1m.tif` | 1 m DEM (EPSG:32613) |
| `boundaries/start_zone.kml` | Topographic start zone boundary |
| `boundaries/domain.kml` | Full simulation domain |
| `boundaries/avalanche_release_area_20260118.geojson` | Observed crown, Jan 18 2026 (UTM) |

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

Expected output (verified against avachain pipeline):
- Start zone clusters: 1608
- Candidates after tau_g/slope/Sk38: 493 → 408 → 204 → 102
- Top-5 triggers: [2859, 5656, 2858, 348, 1068]

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

## Filter chain in generate_scenarios.py

Order matters; do not reorder:
1. Spatial: intersect `features_df.index` with `cluster_map[start_zone_mask]`
2. tau_g ≥ 40 Pa, slope ≥ (stauchwall_deg + 2°), Sk38 < 1.0
3. Slab thickness 0.5–max_slab_thickness m (default 2.0; use 1.5 for skier triggers)
4. Elevation ≥ P50 of remaining candidates
5. Pi1_elastic ≥ median of remaining pool

## CRS handling

All rasters and the cluster map are EPSG:32613 (UTM zone 13N).
GeoJSON inputs may be WGS84 or UTM — CRS is auto-detected via regex
`re.search(r'EPSG:+(\d+)', ...)` to handle both `EPSG:32613` and `EPSG::32613`.

## Mode III Lambda TODO

`elastic_length_cross()` in `arrest_indices.py` currently uses
`Λ_III = √(G h D_wl / G_wl)` giving `Λ_III/Λ_II ≈ 0.59` (ν=0.3).
Update when Johan Gaume's mode III (antiplane) formula is received.

## Commit style

- Short imperative subject line
- End with: `Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>`
