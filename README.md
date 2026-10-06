# avalanche-release-areas

Physically-motivated avalanche release zone delineation using SNOWPACK weak-layer features and the Meloche et al. (2025) crack-arrest scaling laws.

## What this does

Given pre-computed SNOWPACK cluster features (slab density, WL shear strength, slab thickness, elastic length Λ), this package:

1. Computes Meloche (2025) crack-arrest indices Π₁, Π₂, A_ca per cluster
2. Runs a BFS crack-propagation flood-fill from each candidate trigger cluster
3. Outputs release polygons for each scenario (trigger × size factor × depth percentile)
4. Compares against an observed release polygon (IoU metric)

No zarr, no SNOWPACK, no NWP — just CSVs + rasters.

## Included data

`data/little_prof/` contains the Little Professor avalanche path (CDOT I-70 corridor, CO):

| File | Description |
|------|-------------|
| `features/all_start_zone_features_2026-01-18.csv` | Per-cluster WL/slab features at Jan 18, 2026 |
| `features/meloche_features_all_2026-01-18.csv` | Meloche Π₁, A_ca, θ per cluster |
| `spatial/cluster_map.tif` | Cluster ID raster (1 m resolution, UTM zone 13N) |
| `dem_1m.tif` | 1 m DEM (UTM zone 13N) |
| `boundaries/avalanche_release_area_20260118.geojson` | Observed release polygon, Jan 18 2026 |
| `boundaries/start_zone.kml` | Start zone KML boundary |
| `boundaries/domain.kml` | Full simulation domain |

## Installation

```bash
pip install -e ".[dev]"
```

## Quick start

```bash
python scripts/run_scenarios.py \
    --features-csv  data/little_prof/features/all_start_zone_features_2026-01-18.csv \
    --meloche-csv   data/little_prof/features/meloche_features_all_2026-01-18.csv \
    --cluster-map   data/little_prof/spatial/cluster_map.tif \
    --dem           data/little_prof/dem_1m.tif \
    --start-zone    data/little_prof/boundaries/start_zone.kml \
    --release-poly  data/little_prof/boundaries/avalanche_release_area_20260118.geojson \
    --out-dir       outputs/little_prof
```

```bash
python scripts/plot_arrest_indices.py \
    --meloche-csv   data/little_prof/features/meloche_features_all_2026-01-18.csv \
    --cluster-map   data/little_prof/spatial/cluster_map.tif \
    --dem           data/little_prof/dem_1m.tif \
    --release-poly  data/little_prof/boundaries/avalanche_release_area_20260118.geojson \
    --column        A_ca_brittle \
    --out           outputs/arrest_index_map.png
```

## Run tests

```bash
pytest tests/ -v
```

## Physical model

Crack-arrest scaling law (Meloche et al. 2025, JGR Earth Surface):

```
Π₁ = τ_g / (θ · Λ · √(1+δ))
A_ca = L_t · Π₁ · √(σ_t / τ_g)
```

Where:
- τ_g = ρ g h sin ψ — gravitational driving shear stress on WL
- θ = |∇τ_p| — WL shear-strength spatial gradient (Pa/m)
- Λ = √(E′h/K_wl) — elastic length of the slab-WL system (m)
- δ = 1.0 — softening coefficient (Meloche Table 1)
- L_t — characteristic length from fracture toughness
- σ_t — slab tensile strength (Pa)

## Reference

Meloche, J., Gaume, J., Gauthier, D., Hendrikx, J., & Simenhois, R. (2025).
Spatial crack arrest in weak snow layers: field validation of a scaling law.
*Journal of Geophysical Research: Earth Surface*.
doi:10.1029/2025JF008470

## License

CC BY 4.0 — see [LICENSE](LICENSE).
