"""Anisotropic BFS rings vs the lateral distance caps — the six-run comparison.

Reproduces the table in methods §3, "Anisotropic ring growth": each cell is one
combination of {anisotropic rings} x {lateral constraint}, scored on area,
IoU and the crown aspect (L/W and major-axis bearing of the minimum rotated
rectangle), against the observed Jan 18 2026 crown.

    rings  lateral constraint
    ----------------------------------------------------------------
    off    d_lat = min(gaume, A_ca*0.712)   pre-2026-10-10 reference
    on     d_lat = min(gaume, A_ca*0.712)
    off    d_lat = gaume width              <- current default
    on     d_lat = gaume width
    off    no lateral cap at all
    on     no lateral cap at all

The last two cells remove `d_lat` entirely, which no CLI flag does, by patching
`estimate_cross_slope_width` to return a huge width. That is the only way to ask
whether ring growth can *substitute* for the lateral cap rather than merely
compose with it — the question methods §3 answers "substitute, and whichever
constraint is tighter binds".

Trigger selection precedes propagation, so all six cells use the same triggers
and the comparison is automatically trigger-matched.

RE-RUN THIS WHEN Λ_III ARRIVES FROM JOHAN GAUME. Every number in §3 was measured
with Λ_cross == Λ_along: `config.USE_MODE3_LAMBDA` is False and the interim
Λ_III = √(G h D_wl/G_wl) is a constant 0.59 ratio, which the scale-invariant
continuity gate cannot see. Flank continuity — the one input to lateral extent
these runs cannot vary — was therefore governed by the along-slope Λ throughout.
Once the real formula lands: update `arrest_indices.elastic_length_cross()`, set
`USE_MODE3_LAMBDA = True`, **regenerate the meloche CSV** (the shipped Jan 18
one predates the `Lambda_cross` column, and `generate_scenarios` only ever reads
CSVs), then re-run this with `--lambda-cross` to get the 2x2 of
{rings} x {Λ_III}. Judge it on crown aspect, not IoU: the observed crown is
L/W 3.22 at 120.5°, and per the no-tuning rule a ~0.01 IoU margin on n = 1
decides nothing.

Runs in this repo's venv (no zarr needed), from the repo root:

    .venv/bin/python examples/rings_experiment.py
    .venv/bin/python examples/rings_experiment.py --analyse-only
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from shapely.geometry import shape

from release_areas import config, release_geometry
from release_areas import generate_scenarios as gs
from release_areas.release_geometry import load_observed_polygon

DATA = Path('data/little_prof')
OBSERVED = DATA / 'boundaries/avalanche_release_area_20260118.geojson'

# (label, rings, speed cap, drop d_lat entirely)
CELLS = [
    ('rings off  cap on   gaume', False, True,  False),
    ('rings ON   cap on   gaume', True,  True,  False),
    ('rings off  cap off  gaume', False, False, False),   # current default
    ('rings ON   cap off  gaume', True,  False, False),
    ('rings off  no lateral cap', False, False, True),
    ('rings ON   no lateral cap', True,  False, True),
]


def axes(poly):
    """Length, width and major-axis bearing (deg from N) of the min rotated rect.

    Bearing is taken modulo 180: an axis has no head or tail, so 119 and 299
    are the same orientation. Reported alongside L/W because the crown aspect
    is a sharper test of flank physics than IoU — IoU rewards area agreement,
    which a round polygon of the right size can fake.
    """
    rect = poly.minimum_rotated_rectangle
    x, y = rect.exterior.coords.xy
    seg = sorted(((math.hypot(x[i + 1] - x[i], y[i + 1] - y[i]),
                   math.degrees(math.atan2(x[i + 1] - x[i],
                                           y[i + 1] - y[i])) % 180.0)
                  for i in range(4)), reverse=True)
    length, bearing = seg[0]
    width = seg[2][0] if seg[2][0] > 0 else seg[1][0]
    return length, width, (length / width if width else float('nan')), bearing


def run_cell(out_dir: Path, rings: bool, speed_cap: bool, no_lateral_cap: bool,
             lambda_cross: bool, max_clusters: int, n_triggers: int) -> None:
    """One scenario run with the switches set. Mutates config, which is what
    generate_scenarios does with its own CLI flags, so this is the same path."""
    original = release_geometry.estimate_cross_slope_width
    if no_lateral_cap:
        release_geometry.estimate_cross_slope_width = lambda *a, **k: 1e6
    config.USE_MODE3_LAMBDA = lambda_cross
    try:
        sys.argv = [
            'generate_scenarios',
            '--features-csv', str(DATA / 'features/all_start_zone_features_2026-01-18.csv'),
            '--meloche-csv',  str(DATA / 'features/meloche_features_all_2026-01-18.csv'),
            '--cluster-map',  str(DATA / 'spatial/cluster_map.tif'),
            '--dem',          str(DATA / 'dem_1m.tif'),
            '--start-zone',   str(DATA / 'boundaries/start_zone.kml'),
            '--release-poly', str(OBSERVED),
            '--out-dir',      str(out_dir),
            '--max-clusters', str(max_clusters),
            '--n-triggers',   str(n_triggers),
            '--anisotropic-rings' if rings else '--no-anisotropic-rings',
            '--mode3-speed-cap' if speed_cap else '--no-mode3-speed-cap',
        ]
        gs.main()
    finally:
        release_geometry.estimate_cross_slope_width = original


def analyse(out_dir: Path, observed):
    """Area / IoU / aspect per scenario polygon written under out_dir."""
    rows = []
    for path in sorted(out_dir.glob('scenarios/*/release.geojson')):
        feats = json.loads(path.read_text())['features']
        poly = shape(feats[0]['geometry'])
        length, width, lw, bearing = axes(poly)
        rows.append(dict(
            area=poly.area,
            iou=poly.intersection(observed).area / poly.union(observed).area,
            lw=lw, axis=bearing))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out-root', default='outputs/rings_experiment', type=Path)
    ap.add_argument('--max-clusters', type=int, default=config.MAX_BFS_CLUSTERS,
                    help='Keep this at the default: a cap-bound polygon is not '
                         'a physical result, and region size saturates by ~1000.')
    ap.add_argument('--n-triggers', type=int, default=config.N_TOP_TRIGGERS)
    ap.add_argument('--lambda-cross', action='store_true',
                    help='Set USE_MODE3_LAMBDA for every cell — use once '
                         "Gaume's antiplane formula is in elastic_length_cross() "
                         'AND the meloche CSV has been regenerated.')
    ap.add_argument('--analyse-only', action='store_true',
                    help='Re-score existing run directories without re-running.')
    args = ap.parse_args()

    observed = load_observed_polygon(str(OBSERVED))
    length, width, lw, bearing = axes(observed)
    print(f"\nOBSERVED crown  area {observed.area:6.0f} m2   "
          f"L {length:5.1f} m  W {width:5.1f} m  L/W {lw:4.2f}  "
          f"axis {bearing:5.1f} deg\n")

    results = []
    for i, (label, rings, speed_cap, no_lat) in enumerate(CELLS, start=1):
        out_dir = args.out_root / f"{i}_{label.replace(' ', '_')}"
        if not args.analyse_only:
            print(f"\n{'=' * 72}\n=== {label}\n{'=' * 72}")
            run_cell(out_dir, rings, speed_cap, no_lat, args.lambda_cross,
                     args.max_clusters, args.n_triggers)
        if out_dir.exists():
            results.append((label, analyse(out_dir, observed)))

    print(f"\n{'=' * 96}")
    print(f"{'cell':28s} {'area/obs':>9s} {'IoU best':>9s} {'IoU mean':>9s}  "
          f"{'L/W per polygon':<22s} major axis")
    print('-' * 96)
    for label, rows in results:
        if not rows:
            print(f"{label:28s}  (no polygons)")
            continue
        area = np.array([r['area'] for r in rows])
        iou = np.array([r['iou'] for r in rows])
        lws = ' '.join(f"{r['lw']:.2f}" for r in rows)
        axs = ' '.join(f"{r['axis']:.0f}" for r in rows)
        print(f"{label:28s} {area.mean() / observed.area:9.2f} "
              f"{iou.max():9.3f} {iou.mean():9.3f}  {lws:<22s} {axs}")
    print(f"\nObserved L/W {lw:.2f}, axis {bearing:.0f} deg. IoU is capped at "
          f"0.830 on this event by the start-zone KML — report against that "
          f"ceiling, not 1.0.")


if __name__ == '__main__':
    main()
