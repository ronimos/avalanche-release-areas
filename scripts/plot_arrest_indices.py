#!/usr/bin/env python3
"""
plot_arrest_indices.py — Spatial map of Meloche crack-arrest indices.

Reads pre-computed meloche_features CSV + cluster_map.tif and produces
a figure with per-cluster A_ca_brittle (or any other column) coloured
on the DEM hillshade, with the observed release polygon overlaid.

Usage
-----
  python scripts/plot_arrest_indices.py \\
      --meloche-csv data/little_prof/features/meloche_features_all_2026-01-18.csv \\
      --cluster-map data/little_prof/spatial/cluster_map.tif \\
      --dem         data/little_prof/dem_1m.tif \\
      --release-poly data/little_prof/boundaries/avalanche_release_area_20260118.geojson \\
      --column      A_ca_brittle \\
      --out         outputs/arrest_index_map.png
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import rasterio
from matplotlib.colors import LightSource

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))


def _load_observed_polygon(path: Path):
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
        from shapely.geometry import Polygon, MultiPolygon
        def reproj(poly):
            if poly.geom_type == 'Polygon':
                xs, ys = zip(*poly.exterior.coords)
                xs2, ys2 = tr.transform(xs, ys)
                return Polygon(zip(xs2, ys2))
            return MultiPolygon([reproj(p) for p in poly.geoms])
        merged = reproj(merged)

    return merged


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--meloche-csv',  required=True)
    ap.add_argument('--cluster-map',  required=True)
    ap.add_argument('--dem',          required=True)
    ap.add_argument('--release-poly', default=None)
    ap.add_argument('--column',       default='A_ca_brittle',
                    help='Meloche column to map (default: A_ca_brittle)')
    ap.add_argument('--out',          default='arrest_index_map.png')
    args = ap.parse_args()

    meloche_df = pd.read_csv(args.meloche_csv, index_col=0)

    with rasterio.open(args.cluster_map) as src:
        cluster_map = src.read(1).astype(float)
        transform   = src.transform

    with rasterio.open(args.dem) as src:
        dem = src.read(1).astype(float)
        dem[dem <= -9999] = np.nan

    if args.column not in meloche_df.columns:
        print(f"ERROR: column '{args.column}' not in meloche CSV. "
              f"Available: {list(meloche_df.columns)}")
        sys.exit(1)

    # Build per-pixel colour raster from per-cluster values
    col_vals   = meloche_df[args.column].dropna()
    vmin, vmax = float(col_vals.quantile(0.05)), float(col_vals.quantile(0.95))

    color_raster = np.full(cluster_map.shape, np.nan)
    for cid, val in col_vals.items():
        color_raster[cluster_map == cid] = float(val)

    # Figure
    fig, ax = plt.subplots(figsize=(10, 10))

    fill_dem  = np.where(np.isnan(dem), np.nanmean(dem), dem)
    hillshade = LightSource(azdeg=315, altdeg=45).hillshade(fill_dem, dx=1.0, dy=1.0)
    nrows, ncols = dem.shape
    extent = [transform.c,
              transform.c + ncols * transform.a,
              transform.f + nrows * transform.e,
              transform.f]

    ax.imshow(hillshade, cmap='gray', extent=extent, alpha=0.6,
              aspect='auto', origin='upper')

    cmap = plt.cm.RdYlGn
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    im   = ax.imshow(color_raster, cmap=cmap, norm=norm, extent=extent,
                     alpha=0.65, aspect='auto', origin='upper')
    plt.colorbar(im, ax=ax, shrink=0.6, label=args.column)

    if args.release_poly:
        try:
            obs = _load_observed_polygon(Path(args.release_poly))
            xs, ys = obs.exterior.xy
            ax.fill(xs, ys, alpha=0.20, color='red', zorder=3)
            ax.plot(xs, ys, color='red', linewidth=2.5, zorder=4,
                    label=f"Observed ({obs.area:.0f} m²)")
        except Exception as e:
            print(f"Warning: could not load release polygon: {e}")

    ax.legend(loc='upper right', fontsize=9)
    ax.set_title(f"{args.column}  —  per-cluster map", fontsize=11)
    ax.set_xlabel('Easting (m UTM)')
    ax.set_ylabel('Northing (m UTM)')
    ax.ticklabel_format(style='plain', axis='both')
    plt.tight_layout()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), dpi=150, bbox_inches='tight')
    print(f"Saved: {out_path}")
    plt.close(fig)


if __name__ == '__main__':
    main()
