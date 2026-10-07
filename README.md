# avalanche-release-areas

Physically-motivated avalanche release zone delineation using SNOWPACK weak-layer features and the Meloche et al. (2025) crack-arrest scaling laws.

## What this does

Given pre-computed SNOWPACK cluster features (slab density, WL shear strength, slab thickness, elastic length Λ), this package:

1. Computes Meloche (2025) crack-arrest indices Π₁, Π₂, A_ca per cluster
2. Runs spatial crack propagation from each candidate trigger cluster (outward flood-fill with per-direction arrest criteria)
3. Outputs one release polygon per trigger by default (the most likely scenario); optional sweep over size factors and depth percentiles
4. Compares against an observed release polygon (IoU metric)

No zarr, no SNOWPACK, no NWP — just CSVs + rasters.

## Included data

`data/little_prof/` contains the Little Professor avalanche path (CDOT I-70 corridor, CO):

| File | Description |
|------|-------------|
| `features/all_start_zone_features_2026-01-18.csv` | Per-cluster WL/slab features at Jan 18, 2026 |
| `features/meloche_features_all_2026-01-18.csv` | Meloche Π₁, A_ca, θ per cluster |
| `spatial/cluster_map.tif` | Cluster ID raster (1 m, EPSG:6342 — NAD83(2011) / UTM 13N) |
| `dem_1m.tif` | 1 m DEM (EPSG:6342 + NAVD88 height) |
| `boundaries/avalanche_release_area_20260118.geojson` | Observed release polygon, Jan 18 2026 (EPSG:6342, 6 934 m²) |
| `boundaries/20260118_avalanche_boundaries.geojson` | Source mapping in CRS84, before reprojection |
| `boundaries/start_zone.kml` | Start zone KML boundary |
| `boundaries/domain.kml` | Full simulation domain |

## Installation

```bash
pip install -e .          # runtime only
pip install -e ".[dev]"   # + pytest for running tests
```

## Quick start

Run from the repo root. All paths are relative to it.

**Step 1 — generate release scenarios** (one polygon per trigger by default):

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

Outputs written to `--out-dir`:

| File | Description |
|------|-------------|
| `scenarios/scenario_001/release.geojson` | Release polygon for trigger 1 (GeoJSON, EPSG:6342) |
| `scenarios/scenario_001/params.json` | Trigger cluster, A_ca, release area, mean depth, volume, density, Voellmy μ/ξ |
| `scenarios/scenario_001/depth.tif` | Release depth (m) from slab thickness, NaN outside the polygon |
| `scenario_summary.csv` | One row per scenario — area, IoU, trigger cluster |
| `release_comparison.png` | Map overlay: modelled vs observed release |

**Step 2 — re-plot from saved scenarios** (reads GeoJSONs, no recomputation):

```bash
python -m release_areas.plot_release \
    --scenario-dir  outputs/little_prof/scenarios \
    --dem           data/little_prof/dem_1m.tif \
    --start-zone    data/little_prof/boundaries/start_zone.kml \
    --release-poly  data/little_prof/boundaries/avalanche_release_area_20260118.geojson \
    --out           outputs/little_prof/release_comparison.png
```

After `pip install -e .` the same commands are also available as `generate-scenarios` and `plot-release` console scripts.

## Building the input CSVs from SNOWPACK / xsnow

The quick-start above uses pre-computed CSVs from `data/little_prof/features/`. For a new date or path, you need to generate them from a SNOWPACK simulation run loaded into [xsnow](https://xsnow.readthedocs.io).

See [`examples/compute_indices_from_xsnow.py`](examples/compute_indices_from_xsnow.py) for the complete workflow:

1. Open a distributed SNOWPACK run as an xsnow Dataset (from zarr store or `.pro` files)
2. Select the analysis timestep and call `profile_features()` for each cluster profile — extracts slab ρ/h/E/σt, WL shear strength τp, grain type, Sk38, and elastic quantities
3. Call `compute_meloche_features()` across all clusters — adds spatial θ gradient, Π₁, A_ca, R0, and the full arrest-index suite
4. Write the two CSVs that `generate_scenarios` reads directly

**Environment note:** xsnow, xarray, and dask are not in this package's dependencies — they live in the avachain simulation environment. Install this package there with `pip install -e .` and run the example from that environment.

```bash
# From avachain environment, repo root:
python examples/compute_indices_from_xsnow.py
# → data/little_prof/features/all_start_zone_features_2026-01-18.csv
# → data/little_prof/features/meloche_features_all_2026-01-18.csv
```

The script prints a sanity check (τg ≥ 40 Pa candidate count, R0 median) that should match the filter-chain output from `generate_scenarios`.

## Run tests

```bash
pytest tests/ -v
```

## Physical model

See [`docs/release_area_methods.md`](docs/release_area_methods.md) for full derivations, calibration, and the Jan 18 2026 application.

Crack-arrest scaling law (Meloche et al. 2025, JGR Earth Surface):

```
Π₁   = τ_g / (θ · Λ · √(1+δ))
A_ca = C · L_t · Π₁ · √(σ_t / τ_g)
```

Where:
- τ_g = ρ g h sin ψ — gravitational driving shear stress on WL
- θ = |∇τ_p| — WL shear-strength spatial gradient (Pa/m), from k=6 nearest cluster centroids
- Λ = √(E′h/K_wl) — elastic length of the slab-WL system (m)
- δ = 1.0 — softening coefficient (Meloche Table 1)
- L_t = σ_t / k_f — tensile length (distance to first slab fracture, m)
- σ_t — slab tensile strength (Pa)
- C = 0.045 — two-run fit to Meloche et al. (2025) Fig. 8; **do not omit**, it sets the scale of A_ca

## Reference

Meloche, F., Bobillier, G., Guillet, L., Gauthier, F., Langlois, A., & Gaume, J. (2025).
Modeling crack arrest in snow slab avalanches: Toward estimating avalanche release sizes.
*JGR Earth Surface*, 130(12), e2025JF008470.
doi:10.1029/2025JF008470

## License

CC BY 4.0 — see [LICENSE](LICENSE).
