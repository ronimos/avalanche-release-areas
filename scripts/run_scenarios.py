#!/usr/bin/env python3
"""
run_scenarios.py — Run release-zone scenario sweep for one slope and date.

Inputs (all pre-computed, no zarr/SNOWPACK dependency):
  --features-csv   all_start_zone_features_<DATE>.csv   (profile_features output)
  --meloche-csv    meloche_features_all_<DATE>.csv       (compute_meloche_features output)
  --cluster-map    cluster_map.tif  (or .npy + companion .tif for transform)
  --dem            dem_1m.tif
  --start-zone     start_zone.kml  (or .geojson)
  --release-poly   avalanche_release_area_<DATE>.geojson  (observed, for IoU)
  --out-dir        output directory

Outputs:
  <out-dir>/scenarios/scenario_NNN.geojson   one per scenario
  <out-dir>/scenario_summary.csv
  <out-dir>/release_comparison.png

Usage
-----
  python scripts/run_scenarios.py \\
      --features-csv  data/little_prof/features/all_start_zone_features_2026-01-18.csv \\
      --meloche-csv   data/little_prof/features/meloche_features_all_2026-01-18.csv \\
      --cluster-map   data/little_prof/spatial/cluster_map.tif \\
      --dem           data/little_prof/dem_1m.tif \\
      --start-zone    data/little_prof/boundaries/start_zone.kml \\
      --release-poly  data/little_prof/boundaries/avalanche_release_area_20260118.geojson \\
      --out-dir       outputs/little_prof \\
      --n-triggers 5 \\
      --size-factors 0.75 0.875 1.0 1.125 1.25 \\
      --depth-pcts   0.25 0.5 0.75
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from release_areas.release_geometry import (
    make_release_polygon_2d,
    plot_release_comparison,
)
from release_areas.scenario_writer import write_scenario


def _load_kml_mask(kml_path: Path, dem_shape, transform) -> np.ndarray:
    from fastkml import kml
    from shapely.geometry import shape
    from shapely.ops import unary_union
    import rasterio.features

    k = kml.KML()
    with open(str(kml_path), 'rb') as f:
        k.from_string(f.read())

    polys = []
    for doc in k.features():
        for feat in doc.features():
            try:
                polys.append(shape(feat.geometry))
            except Exception:
                pass

    if not polys:
        return np.zeros(dem_shape, dtype=bool)

    merged = unary_union(polys)
    mask   = rasterio.features.geometry_mask(
        [merged.__geo_interface__], out_shape=dem_shape,
        transform=transform, invert=True)
    return mask


def _load_geojson_mask(path: Path, dem_shape, transform) -> np.ndarray:
    from release_areas.features import geojson_to_mask
    return geojson_to_mask(path, dem_shape, transform)


def _load_observed_polygon(path: Path):
    import json
    import re
    from shapely.geometry import shape
    from shapely.ops import unary_union
    from pyproj import Transformer

    with open(str(path)) as f:
        gj = json.load(f)

    polys = [shape(feat['geometry']) for feat in gj['features']]
    merged = unary_union(polys)

    src_epsg = 4326
    crs_node = gj.get('crs', {}).get('properties', {}).get('name', '')
    m = re.search(r'EPSG:+(\d+)', crs_node, re.IGNORECASE)
    if m:
        src_epsg = int(m.group(1))

    if src_epsg != 32613:
        tr = Transformer.from_crs(f'EPSG:{src_epsg}', 'EPSG:32613', always_xy=True)
        from shapely.geometry import Polygon
        def reproj(poly):
            if poly.geom_type == 'Polygon':
                xs, ys = zip(*poly.exterior.coords)
                xs2, ys2 = tr.transform(xs, ys)
                return Polygon(zip(xs2, ys2))
            from shapely.geometry import MultiPolygon
            return MultiPolygon([reproj(p) for p in poly.geoms])
        merged = reproj(merged)

    return merged


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
    ap.add_argument('--n-triggers',    type=int,   default=5)
    ap.add_argument('--size-factors',  type=float, nargs='+',
                    default=[0.75, 0.875, 1.0, 1.125, 1.25])
    ap.add_argument('--depth-pcts',    type=float, nargs='+',
                    default=[0.25, 0.50, 0.75])
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

    with rasterio.open(args.dem) as src:
        dem = src.read(1).astype(float)
        dem[dem <= -9999] = np.nan

    sz_path = Path(args.start_zone)
    if sz_path.suffix.lower() == '.kml':
        start_zone_mask = _load_kml_mask(sz_path, dem.shape, transform)
    else:
        start_zone_mask = _load_geojson_mask(sz_path, dem.shape, transform)

    observed_polygon = None
    if args.release_poly:
        try:
            observed_polygon = _load_observed_polygon(Path(args.release_poly))
            print(f"Observed polygon: {observed_polygon.area:.0f} m2")
        except Exception as e:
            print(f"Warning: could not load observed polygon: {e}")

    # --- Select top-N trigger clusters ---
    if 'A_ca_brittle' not in meloche_df.columns:
        print("ERROR: meloche CSV missing 'A_ca_brittle'. Check compute_meloche_features output.")
        sys.exit(1)

    sorted_mel = (meloche_df['A_ca_brittle']
                  .dropna()
                  .sort_values(ascending=False))
    trigger_ids = list(sorted_mel.index[:args.n_triggers])
    print(f"Top-{args.n_triggers} triggers: {trigger_ids}")

    # --- Scenario sweep ---
    scenarios   = []
    mel_polys   = []
    trig_cents  = []
    trig_labels = []
    scenario_n  = 0

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
                )

                area_m2 = float(poly.area) if poly else 0.0
                iou     = 0.0
                if poly and observed_polygon:
                    inter = poly.intersection(observed_polygon).area
                    union = poly.union(observed_polygon).area
                    iou   = inter / union if union > 0 else 0.0

                scenarios.append({
                    'scenario_id':      sid,
                    'trigger_cluster':  int(trig_id),
                    'size_factor':      sf,
                    'depth_pct':        dpct,
                    'area_m2':          area_m2,
                    'A_ca':             A_ca_val,
                    'iou':              iou,
                })

                if poly:
                    out_path = out_dir / 'scenarios' / f'{sid}.geojson'
                    write_scenario(
                        scenario_id         = sid,
                        trigger_cluster_id  = trig_id,
                        release_polygon     = poly,
                        depth_grid          = np.full(dem.shape, np.nan),
                        transform           = transform,
                        meloche_row         = meloche_df.loc[trig_id],
                        out_path            = out_path,
                    )

        mel_polys.append((poly, sf))
        trig_labels.append(f"T{trigger_ids.index(trig_id)+1} (cid={trig_id})")

    # --- Summary CSV ---
    summary = pd.DataFrame(scenarios)
    summary_path = out_dir / 'scenario_summary.csv'
    summary.to_csv(str(summary_path))
    print(f"\nScenario summary: {summary_path}")
    if not summary.empty:
        print(f"  {len(summary)} scenarios  |  "
              f"area range: {summary['area_m2'].min():.0f}–{summary['area_m2'].max():.0f} m²")
        best = summary.loc[summary['iou'].idxmax()]
        print(f"  Best IoU: {best['iou']:.3f}  ({best['scenario_id']}, "
              f"sf={best['size_factor']}, area={best['area_m2']:.0f} m²)")

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
        out_path          = str(fig_path),
        title             = "Release zone scenarios vs observed",
    )
    print(f"Comparison figure: {fig_path}")


if __name__ == '__main__':
    main()
