"""
generate_scenarios — Generate release-zone scenarios for one slope and date.

Inputs (all pre-computed, no zarr/SNOWPACK dependency):
  --features-csv   all_start_zone_features_<DATE>.csv   (profile_features output)
  --meloche-csv    meloche_features_all_<DATE>.csv       (compute_meloche_features output)
  --cluster-map    cluster_map.tif
  --dem            dem_1m.tif
  --start-zone     start_zone.kml  (or .geojson)
  --release-poly   avalanche_release_area_<DATE>.geojson  (observed, for IoU)
  --out-dir        output directory

Outputs:
  <out-dir>/scenarios/scenario_NNN/   one directory per scenario (release.geojson,
                                       depth.tif, depth.asc, params.json, density.json)
  <out-dir>/scenario_summary.csv
  <out-dir>/release_comparison.png

Usage
-----
  # One (most likely) scenario per trigger — default:
  generate-scenarios \\
      --features-csv  data/little_prof/features/all_start_zone_features_2026-01-18.csv \\
      --meloche-csv   data/little_prof/features/meloche_features_all_2026-01-18.csv \\
      --cluster-map   data/little_prof/spatial/cluster_map.tif \\
      --dem           data/little_prof/dem_1m.tif \\
      --start-zone    data/little_prof/boundaries/start_zone.kml \\
      --release-poly  data/little_prof/boundaries/avalanche_release_area_20260118.geojson \\
      --out-dir       outputs/little_prof

  # Full sensitivity sweep (5 triggers × 5 sizes × 3 depths = 75 scenarios):
  generate-scenarios ... \\
      --size-factors 0.75 0.875 1.0 1.125 1.25 \\
      --depth-pcts   0.25 0.5 0.75
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

from release_areas import config
from release_areas.release_geometry import (
    make_release_polygon_2d,
    plot_release_comparison,
    load_observed_polygon,
)
from release_areas.scenario_writer import write_scenario
from release_areas.snowpack_features import geojson_to_mask


def load_kml_mask(kml_path: Path, dem_shape, transform,
                  dst_epsg: int = 32613) -> np.ndarray:
    """Parse a KML polygon using stdlib XML (no fastkml dependency)."""
    import xml.etree.ElementTree as ET
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    from pyproj import Transformer
    import rasterio.features

    tree = ET.parse(str(kml_path))
    root = tree.getroot()

    KML_NS = 'http://www.opengis.net/kml/2.2'
    rings = []
    for coords_el in root.iter(f'{{{KML_NS}}}coordinates'):
        pts = []
        for triple in coords_el.text.strip().split():
            parts = triple.split(',')
            pts.append((float(parts[0]), float(parts[1])))  # lon, lat
        if len(pts) >= 3:
            rings.append(pts)

    if not rings:
        return np.zeros(dem_shape, dtype=bool)

    tr = Transformer.from_crs('EPSG:4326', f'EPSG:{dst_epsg}', always_xy=True)
    polys = []
    for ring in rings:
        lons, lats = zip(*ring)
        xs, ys = tr.transform(lons, lats)
        polys.append(Polygon(zip(xs, ys)))

    merged = unary_union(polys)
    mask   = rasterio.features.geometry_mask(
        [merged.__geo_interface__], out_shape=dem_shape,
        transform=transform, invert=True)
    return mask


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--features-csv',  required=True)
    ap.add_argument('--meloche-csv',   required=True)
    ap.add_argument('--cluster-map',   required=True)
    ap.add_argument('--dem',           required=True)
    ap.add_argument('--start-zone',    required=True)
    ap.add_argument('--release-poly',  default=None)
    ap.add_argument('--out-dir',       required=True)
    ap.add_argument('--n-triggers',         type=int,   default=config.N_TOP_TRIGGERS)
    ap.add_argument('--size-factors',       type=float, nargs='+', default=[1.0],
                    help='Release size multipliers (default: [1.0] — most likely only)')
    ap.add_argument('--depth-pcts',         type=float, nargs='+', default=[0.5],
                    help='Depth percentiles (default: [0.5] — most likely only)')
    ap.add_argument('--stauchwall-deg',     type=float, default=28.0,
                    help='Slope threshold for stauchwall arrest (default: 28°)')
    ap.add_argument('--max-slab-thickness', type=float, default=2.0,
                    help='Maximum slab thickness for trigger candidates (m). '
                         'Use 1.5 for skier-triggered scenarios, 2.0 for natural/large-slab events.')
    ap.add_argument('--mode3-scale',        type=float, default=1.5,
                    help='Mode III lateral arrest multiplier (default: 1.5)')
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    (out_dir / 'scenarios').mkdir(parents=True, exist_ok=True)

    # --- Load inputs ---
    features_df = pd.read_csv(args.features_csv, index_col=0)
    meloche_df  = pd.read_csv(args.meloche_csv,  index_col=0)

    with rasterio.open(args.cluster_map) as src:
        cluster_map = src.read(1).astype(float)
        transform   = src.transform
        crs         = src.crs
        crs_wkt     = crs.to_wkt()

    with rasterio.open(args.dem) as src:
        dem = src.read(1).astype(float)
        dem[dem <= -9999] = np.nan

    sz_path = Path(args.start_zone)
    if sz_path.suffix.lower() == '.kml':
        start_zone_mask = load_kml_mask(sz_path, dem.shape, transform)
    else:
        start_zone_mask = geojson_to_mask(sz_path, dem.shape, transform)

    observed_polygon = None
    if args.release_poly:
        try:
            observed_polygon = load_observed_polygon(Path(args.release_poly))
            print(f"Observed polygon: {observed_polygon.area:.0f} m2")
        except Exception as e:
            print(f"Warning: could not load observed polygon: {e}")

    # --- Select top-N trigger clusters (lowest Sk38 = most unstable) ---
    sk38_col = 'min_sk38' if 'min_sk38' in features_df.columns else 'sk38_min'

    MIN_TAU_G          = 40.0
    MAX_SK38           = 1.0
    MIN_SLAB_THICKNESS = 0.5
    MAX_SLAB_THICKNESS = args.max_slab_thickness
    MIN_SLOPE_TRIGGER  = args.stauchwall_deg + 2.0

    sz_cids = set(int(c) for c in np.unique(cluster_map[start_zone_mask]) if c > 0)
    candidate_ids = [cid for cid in features_df.index if cid in sz_cids]
    print(f"  Start zone clusters with features: {len(candidate_ids)}")

    def _scalar(df, cid, col):
        if cid not in df.index or col not in df.columns:
            return np.nan
        v = df.loc[cid, col]
        if isinstance(v, pd.Series):
            v = v.iloc[0]
        return float(v)

    # Filter 1: tau_g, slope, Sk38
    if not meloche_df.empty and 'tau_g' in meloche_df.columns:
        candidate_ids = [
            cid for cid in candidate_ids
            if cid in meloche_df.index
            and _scalar(meloche_df, cid, 'tau_g') >= MIN_TAU_G
            and _scalar(meloche_df, cid, 'slope_angle') >= MIN_SLOPE_TRIGGER
            and _scalar(features_df, cid, sk38_col) < MAX_SK38
        ]
        print(f"  Candidates after tau_g >= {MIN_TAU_G} Pa, "
              f"slope >= {MIN_SLOPE_TRIGGER:.1f}°, "
              f"Sk38 < {MAX_SK38} filter: {len(candidate_ids)}")

    # Filter 2: slab thickness
    candidate_ids = [
        cid for cid in candidate_ids
        if np.isnan(_scalar(features_df, cid, 'slab_thickness'))
        or (MIN_SLAB_THICKNESS
            <= _scalar(features_df, cid, 'slab_thickness')
            <= MAX_SLAB_THICKNESS)
    ]
    print(f"  Candidates after slab_thickness {MIN_SLAB_THICKNESS}–{MAX_SLAB_THICKNESS} m filter: "
          f"{len(candidate_ids)}")

    # Filter 3: elevation — keep upper P50 of candidate zone
    if candidate_ids:
        elev_pairs = []
        for cid in candidate_ids:
            pxs = np.argwhere(cluster_map == cid)
            if len(pxs):
                elev_pairs.append((cid, float(dem[pxs[:, 0], pxs[:, 1]].mean())))
        if elev_pairs:
            elev_threshold = np.percentile([e for _, e in elev_pairs], 50)
            candidate_ids = [cid for cid, e in elev_pairs if e >= elev_threshold]
            print(f"  Candidates after elevation >= P50 "
                  f"({elev_threshold:.0f} m) filter: {len(candidate_ids)}")

    if not candidate_ids:
        print("ERROR: no candidate trigger clusters survived filters.")
        sys.exit(1)

    cand_df = features_df.loc[candidate_ids].copy()

    # Filter 4: Pi1 propagation gate (median)
    if not meloche_df.empty and 'Pi1_elastic' in meloche_df.columns:
        pi1_vals = {}
        for cid in cand_df.index:
            if cid in meloche_df.index:
                v = meloche_df.loc[cid, 'Pi1_elastic']
                pi1_vals[cid] = float(v.iloc[0]) if isinstance(v, pd.Series) else float(v)
        cand_df['Pi1_elastic'] = pd.Series(pi1_vals)
        pi1_valid = cand_df['Pi1_elastic'].dropna()
        if len(pi1_valid) > 0:
            pi1_median = float(pi1_valid.median())
            pre_gate = cand_df.dropna(subset=[sk38_col]).shape[0]
            cand_df = cand_df[cand_df['Pi1_elastic'] >= pi1_median]
            print(f"  Pi1 propagation gate (>= median {pi1_median:.2f}): "
                  f"{pre_gate} → {cand_df.shape[0]} candidates")

    triggers_df = (cand_df.dropna(subset=[sk38_col])
                          .nsmallest(args.n_triggers, sk38_col))
    trigger_ids = list(triggers_df.index)
    sk38_vals   = list(triggers_df[sk38_col].round(3))
    print(f"Top-{args.n_triggers} triggers (lowest Sk38): {trigger_ids}  Sk38={sk38_vals}")

    # --- Build depth raster from HS values in features CSV ---
    depth_raster = np.full(dem.shape, np.nan, dtype=np.float32)
    if 'hs' in features_df.columns:
        for cid in features_df.index:
            hs_val = _scalar(features_df, cid, 'hs')
            if np.isfinite(hs_val):
                depth_raster[cluster_map == cid] = float(hs_val)

    # --- Scenario sweep ---
    total_scenarios = len(trigger_ids) * len(args.size_factors) * len(args.depth_pcts)
    scenario_weight = 1.0 / total_scenarios if total_scenarios > 0 else 1.0

    scenarios   = []
    mel_polys   = []
    trig_cents  = []
    trig_labels = []
    trig_ious   = []
    scenario_n  = 0
    poly        = None  # ensure defined for mel_polys fallback

    for trig_id in trigger_ids:
        A_ca_row = meloche_df.loc[trig_id]
        A_ca_val = float(A_ca_row['A_ca_brittle']
                         if not isinstance(A_ca_row['A_ca_brittle'], pd.Series)
                         else A_ca_row['A_ca_brittle'].iloc[0])

        # Trigger centroid
        pxs = np.argwhere(cluster_map == trig_id)
        if len(pxs):
            r, c = pxs.mean(axis=0).astype(int)
            tx = transform.c + c * transform.a + r * transform.b
            ty = transform.f + c * transform.d + r * transform.e
            trig_cents.append((tx, ty))
        else:
            trig_cents.append(None)

        # Density for this trigger (used by write_scenario for params.json / density.json)
        dens_col = 'slab_density'
        if trig_id in features_df.index and dens_col in features_df.columns:
            d_series = pd.to_numeric(features_df[dens_col], errors='coerce')
            density_mean = float(d_series.loc[trig_id]) if not pd.isna(d_series.get(trig_id)) else float(d_series.median())
            density_std  = float(d_series.dropna().std()) if d_series.dropna().size > 1 else 0.0
        else:
            density_mean, density_std = 250.0, 0.0
        if not np.isfinite(density_mean):
            density_mean = 250.0

        scenario_start = len(scenarios)

        for sf in args.size_factors:
            for dpct in args.depth_pcts:
                scenario_n += 1
                sid = f"scenario_{scenario_n:03d}"

                poly = make_release_polygon_2d(
                    trigger_cluster_id = trig_id,
                    A_ca               = A_ca_val,
                    meloche_df         = meloche_df,
                    cluster_map        = cluster_map,
                    dem                = dem,
                    transform          = transform,
                    start_zone_mask    = start_zone_mask,
                    snap_features      = features_df,
                    size_factor        = sf,
                    stauchwall_deg     = args.stauchwall_deg,
                    mode3_scale        = args.mode3_scale,
                )

                iou = 0.0
                if poly and observed_polygon:
                    inter = poly.intersection(observed_polygon).area
                    union = poly.union(observed_polygon).area
                    iou   = inter / union if union > 0 else 0.0

                if poly:
                    row = write_scenario(
                        scenario_dir       = out_dir / 'scenarios',
                        scenario_id        = sid,
                        release_polygon    = poly,
                        depth_raster       = depth_raster,
                        dem_shape          = dem.shape,
                        transform          = transform,
                        crs_wkt            = crs_wkt,
                        density_mean       = density_mean,
                        density_std        = density_std,
                        trigger_cluster_id = int(trig_id),
                        A_ca               = A_ca_val,
                        size_factor        = sf,
                        depth_percentile   = int(round(dpct * 100)),
                        weight             = scenario_weight,
                    )
                    row['iou'] = iou
                else:
                    row = {
                        'scenario_id':      sid,
                        'trigger_cluster':  int(trig_id),
                        'A_ca_m':           round(A_ca_val, 2),
                        'size_factor':      sf,
                        'depth_percentile': int(round(dpct * 100)),
                        'release_area_m2':  0.0,
                        'iou':              0.0,
                    }
                scenarios.append(row)

        # Best IoU for this trigger across all size/depth combinations
        trig_best_iou = max(
            (r.get('iou', 0.0) for r in scenarios[scenario_start:]),
            default=0.0)
        trig_ious.append(trig_best_iou)

        mel_polys.append((poly, sf))
        sk38_trig = _scalar(features_df, trig_id, sk38_col)
        trig_labels.append(
            f"T{trigger_ids.index(trig_id)+1} cid={trig_id} "
            f"Sk38={sk38_trig:.2f} A_ca={A_ca_val:.0f}m"
        )

    # --- Summary CSV ---
    summary = pd.DataFrame(scenarios)
    summary_path = out_dir / 'scenario_summary.csv'
    summary.to_csv(str(summary_path))
    print(f"\nScenario summary: {summary_path}")
    if not summary.empty:
        areas = summary['release_area_m2'].dropna()
        print(f"  {len(summary)} scenarios  |  "
              f"area range: {areas.min():.0f}–{areas.max():.0f} m²")
        if summary['iou'].any():
            best = summary.loc[summary['iou'].idxmax()]
            print(f"  Best IoU: {best['iou']:.3f}  ({best['scenario_id']}, "
                  f"sf={best['size_factor']}, area={best['release_area_m2']:.0f} m²)")

    # --- Comparison plot ---
    fig_path = out_dir / 'release_comparison.png'
    plot_release_comparison(
        meloche_polygons  = mel_polys,
        observed_polygon  = observed_polygon,
        dem               = dem,
        transform         = transform,
        start_zone_mask   = start_zone_mask,
        trigger_labels    = trig_labels,
        trigger_centroids = trig_cents,
        trigger_ious      = trig_ious,
        out_path          = str(fig_path),
        title             = "Release zone scenarios vs observed",
    )
    print(f"Comparison figure: {fig_path}")


if __name__ == '__main__':
    main()
