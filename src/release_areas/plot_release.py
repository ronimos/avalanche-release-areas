"""
plot_release — Plot modeled release polygons vs. observed and compute IoU.

Reads pre-computed scenario GeoJSONs from a generate-scenarios output directory,
overlays them on hillshade/DEM, and annotates each polygon with its IoU score.

Usage
-----
  plot-release \\
      --scenario-dir  outputs/little_prof/scenarios \\
      --dem           data/little_prof/dem_1m.tif \\
      --start-zone    data/little_prof/boundaries/start_zone.kml \\
      --release-poly  data/little_prof/boundaries/avalanche_release_area_20260118.geojson \\
      --out           outputs/little_prof/release_comparison.png

  # Filter to specific scenario IDs:
  plot-release ... --scenario-ids scenario_001 scenario_002
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import rasterio
from shapely.geometry import shape

from release_areas.release_geometry import plot_release_comparison, load_observed_polygons
from release_areas.snowpack_features import geojson_to_mask
from release_areas.generate_scenarios import load_kml_mask


def _load_scenario_polygons(scenario_dir: Path, ids: list[str] | None = None):
    """Yield (polygon, size_factor, label) for each release.geojson found."""
    polys = []
    for release_json in sorted(scenario_dir.glob("*/release.geojson")):
        sid = release_json.parent.name
        if ids and sid not in ids:
            continue
        with open(release_json) as f:
            gj = json.load(f)
        features = gj.get("features", [gj] if gj.get("type") == "Feature" else [])
        if not features:
            continue
        geom = shape(features[0]["geometry"])
        props = features[0].get("properties", {})
        sf = props.get("size_factor", 1.0)
        polys.append((geom, sf, sid))
    return polys


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--scenario-dir',  required=True,
                    help='Directory containing scenario_NNN/ subdirectories')
    ap.add_argument('--dem',           required=True)
    ap.add_argument('--start-zone',    required=True)
    ap.add_argument('--release-poly',  default=None,
                    help='Observed release polygon GeoJSON (for IoU)')
    ap.add_argument('--out',           required=True, help='Output PNG path')
    ap.add_argument('--scenario-ids',  nargs='+', default=None,
                    help='Limit to specific scenario IDs (e.g. scenario_001)')
    ap.add_argument('--title',         default='Release zone comparison')
    args = ap.parse_args()

    scenario_dir = Path(args.scenario_dir)

    with rasterio.open(args.dem) as src:
        dem       = src.read(1).astype(float)
        transform = src.transform
        dem[dem <= -9999] = np.nan

    sz_path = Path(args.start_zone)
    if sz_path.suffix.lower() == '.kml':
        start_zone_mask = load_kml_mask(sz_path, dem.shape, transform)
    else:
        start_zone_mask = geojson_to_mask(sz_path, dem.shape, transform)

    observed_polygons = []
    observed_polygon = None
    if args.release_poly:
        try:
            observed_polygons = load_observed_polygons(Path(args.release_poly))
            for poly, lbl in observed_polygons:
                print(f"Observed polygon ({lbl}): {poly.area:.0f} m2")
            observed_polygon = observed_polygons[0][0] if observed_polygons else None
        except Exception as e:
            print(f"Warning: could not load observed polygons: {e}")

    entries = _load_scenario_polygons(scenario_dir, args.scenario_ids)
    if not entries:
        print(f"No scenario GeoJSONs found in {scenario_dir}")
        return

    mel_polys   = [(p, sf) for p, sf, _ in entries]
    trig_labels = [sid for _, _, sid in entries]
    trig_ious   = []
    trig_cents  = []

    for poly, _, _ in entries:
        cx, cy = poly.centroid.x, poly.centroid.y
        trig_cents.append((cx, cy))

        iou = 0.0
        if observed_polygon:
            inter = poly.intersection(observed_polygon).area
            union = poly.union(observed_polygon).area
            iou   = inter / union if union > 0 else 0.0
        trig_ious.append(iou)

    for lbl, iou in zip(trig_labels, trig_ious):
        print(f"  {lbl}: IoU={iou:.3f}")

    plot_release_comparison(
        meloche_polygons  = mel_polys,
        observed_polygons = observed_polygons,
        dem               = dem,
        transform         = transform,
        start_zone_mask   = start_zone_mask,
        trigger_labels    = trig_labels,
        trigger_centroids = trig_cents,
        trigger_ious      = trig_ious,
        out_path          = args.out,
        title             = args.title,
    )
    print(f"Saved: {args.out}")


if __name__ == '__main__':
    main()
