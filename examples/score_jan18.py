#!/usr/bin/env python3
"""Score the release model on Jan 18 2026 at a calibration preset.

Port of avachain's ``scripts/inside_release_sweep.py`` and the five-trigger
``event_scoring`` (ops-2026-27), run on this repo's data/little_prof. The
release model is this repo's; the SNOWPACK features come from avachain
(``data/little_prof/features/cal_<version>``).

For each slab-crack mode:

  five-trigger   the model's own top-n triggers (select_triggers), each
                 weighted 1/Sk38: E[IoU], E[area], AUC, Brier skill
  inside         every cluster with >= half its pixels inside the observed
                 release (and the start zone), as a forced shot: released
                 count, IoU median [IQR], consensus and pairwise IoU
  ski            every cluster under an observed ski non-trigger point:
                 must not release (0/46 at 2026-11-17.1)

Forced shots have no initiation gate and no rectangle fallback. When the
trigger's own tensile gate stops a shot, the polygon is also built with that
gate removed ("ungated"); the BFS's per-cluster gate stays on.

Jan 18 is a calibration target: these are calibration numbers, never
validation.

    env -u VIRTUAL_ENV uv run --frozen python examples/score_jan18.py \\
        --modes shear opening --jobs 12 --out outputs/score_jan18
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import multiprocessing as mproc
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

from release_areas import config
from release_areas import release_geometry as _rg
from release_areas import scoring
from release_areas.generate_scenarios import (load_kml_mask, select_triggers,
                                              sk38_column)

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / 'data' / 'little_prof'
DATE = '2026-01-18'
A_CA_MAX_M, A_CA_FALLBACK_M = 300.0, 50.0     # avachain model_params

G = {}   # inputs, set in the parent before the pool forks


def apply_calibration(cal: dict) -> None:
    """Set the release-model knobs of a calibration preset on config
    (avachain calibration.apply). facet_a also acts at feature generation,
    which the preset's features_dir already carries."""
    if 'facet_a' in cal:
        config.FACET_STRENGTH_FACTOR = float(cal['facet_a'])
    if 'distance_caps' in cal:
        config.USE_DISTANCE_CAPS = bool(cal['distance_caps'])
    if 'tau_ratio_floor' in cal:
        config.TAU_RATIO_FLOOR = float(cal['tau_ratio_floor'])
    if 'slab_crack_mode' in cal:
        config.SLAB_CRACK_MODE = cal['slab_crack_mode']
    if 'downslope_slide' in cal:
        config.USE_DOWNSLOPE_SLIDE = bool(cal['downslope_slide'])
    if 'unified_slide' in cal:
        config.UNIFIED_SLIDE = bool(cal['unified_slide'])


_MIXED_DEFAULTS = {}


def set_mode(variant: str) -> None:
    """'shear' | 'opening' | 'mixed:<criterion>:<on_fail>[:<axes>]' (e.g.
    mixed:energy:stop:fall) -> config. A part left out takes the config
    default at import (MIXED_AXES 'fall', as in the 11-18 runs; it was
    'front_edge' for part of 2026-11-19)."""
    if not _MIXED_DEFAULTS:
        _MIXED_DEFAULTS.update(c=config.MIXED_CRITERION, f=config.MIXED_ON_FAIL,
                               a=config.MIXED_AXES)
    parts = variant.split(':')
    config.SLAB_CRACK_MODE = parts[0]
    if parts[0] == 'mixed':
        config.MIXED_CRITERION = parts[1] if len(parts) > 1 else _MIXED_DEFAULTS['c']
        config.MIXED_ON_FAIL = parts[2] if len(parts) > 2 else _MIXED_DEFAULTS['f']
        config.MIXED_AXES = parts[3] if len(parts) > 3 else _MIXED_DEFAULTS['a']


def load_inputs(features_dir: Path, observed: Path, ski: Path,
                data_dir: Path = DATA, date: str = DATE) -> None:
    """One event. data_dir holds the slope's spatial/cluster_map.tif,
    dem_1m.tif and boundaries/start_zone.kml (the data/little_prof layout)."""
    with rasterio.open(data_dir / 'spatial' / 'cluster_map.tif') as src:
        cluster_map = src.read(1).astype(float)
        transform = src.transform
    with rasterio.open(data_dir / 'dem_1m.tif') as src:
        dem = src.read(1).astype(float)
        dem[dem <= -9999] = np.nan
    sz = load_kml_mask(data_dir / 'boundaries' / 'start_zone.kml', dem.shape, transform)
    feats = pd.read_csv(features_dir / f'all_start_zone_features_{date}.csv', index_col=0)
    mel = pd.read_csv(features_dir / f'meloche_features_all_{date}.csv', index_col=0)
    obs_poly = _rg.load_observed_polygon(observed)
    obs_mask = scoring.rasterize([obs_poly], dem.shape, transform)[0]
    G.update(cluster_map=cluster_map, transform=transform, dem=dem, start_zone_mask=sz,
             feats=feats[~feats.index.duplicated()], mel=mel[~mel.index.duplicated()],
             observed=obs_poly, obs_mask=obs_mask,
             ski=json.loads(ski.read_text()) if ski else None,
             pixel_index=_rg.cluster_pixel_index(cluster_map))


def _a_ca(cid) -> float:
    mel = G['mel']
    if 'A_ca_brittle' in mel.columns and cid in mel.index:
        v = float(mel.loc[cid, 'A_ca_brittle'])
        if np.isfinite(v):
            return min(v, A_CA_MAX_M)
    return A_CA_FALLBACK_M


def _polygon(cid):
    state = {}
    poly = _rg.make_release_polygon_2d(
        trigger_cluster_id=cid, A_ca=_a_ca(cid), meloche_df=G['mel'],
        cluster_map=G['cluster_map'], dem=G['dem'], transform=G['transform'],
        start_zone_mask=G['start_zone_mask'], snap_features=G['feats'],
        size_factor=1.0, stauchwall_deg=config.STAUCHWALL_DEG,
        max_clusters=config.MAX_BFS_CLUSTERS, pixel_index=G['pixel_index'],
        state_out=state, rectangle_fallback=False)
    return poly, state


def run_one(job):
    """One (mode, cid, class) forced shot; polygons travel back as WKB."""
    mode, cid, cls = job
    set_mode(mode)
    out = {'mode': mode, 'cid': cid, 'class': cls}
    s_t = _rg.cluster_tensile_support(cid, G['mel'], G['feats'])
    out['trigger_s_t'] = float(s_t) if np.isfinite(s_t) else None
    with contextlib.redirect_stdout(io.StringIO()):
        poly, state = _polygon(cid)
        out['outcome'] = ('released' if poly is not None else
                          'trigger_gate' if state.get('trigger_tensile_failure') else
                          'no_propagation')
        out['wkb'] = poly.wkb if poly is not None else None
        out['wkb_ungated'] = out['wkb']
        if out['outcome'] == 'trigger_gate':
            orig = _rg.cluster_tensile_support
            _rg.cluster_tensile_support = lambda *a, **k: np.inf
            try:
                p2, _ = _polygon(cid)
            finally:
                _rg.cluster_tensile_support = orig
            out['wkb_ungated'] = p2.wkb if p2 is not None else None
    return out


def summarise(rows) -> list:
    from shapely import wkb as _wkb
    obs, shape, tf = G['observed'], G['dem'].shape, G['transform']
    for r in rows:
        for k, pk in (('wkb', 'poly'), ('wkb_ungated', 'poly_ungated')):
            r[pk] = _wkb.loads(r[k]) if r[k] else None
        r['iou'] = scoring.iou(r['poly'], obs)
        r['iou_ungated'] = scoring.iou(r['poly_ungated'], obs)
        r['area_m2'] = r['poly'].area if r['poly'] is not None else 0.0
        r['area_ungated_m2'] = r['poly_ungated'].area if r['poly_ungated'] is not None else 0.0
    out = []
    for mode in dict.fromkeys(r['mode'] for r in rows):
        for cls in ('interior', 'edge', 'inside', 'ski', 'extra'):
            sel = [r for r in rows if r['mode'] == mode and
                   (r['class'] == cls or (cls == 'inside' and r['class'] in ('interior', 'edge')))]
            if not sel:
                continue
            rec = {'mode': mode, 'class': cls, 'n': len(sel)}
            for k in ('released', 'trigger_gate', 'no_propagation'):
                rec[f'n_{k}'] = sum(r['outcome'] == k for r in sel)
            for tag, pk, ik in (('', 'poly', 'iou'), ('_ungated', 'poly_ungated', 'iou_ungated')):
                rel = [r for r in sel if r[pk] is not None]
                rec[f'n_poly{tag}'] = len(rel)
                if rel:
                    ious = np.array([r[ik] for r in rel])
                    rec[f'iou_median{tag}'] = float(np.median(ious))
                    rec[f'iou_p25{tag}'] = float(np.percentile(ious, 25))
                    rec[f'iou_p75{tag}'] = float(np.percentile(ious, 75))
                    rec[f'area_median{tag}_m2'] = float(np.median([r[pk].area for r in rel]))
                if cls in ('inside', 'ski'):
                    masks = scoring.rasterize([r[pk] for r in rel], shape, tf)
                    c, p = scoring.consensus(masks, G['obs_mask'])
                    if cls == 'inside' and tag == '' and masks:
                        G.setdefault('freqs', {})[mode] = np.mean(masks, axis=0)
                    rec[f'consensus_iou{tag}'] = c
                    rec[f'pairwise_iou{tag}'] = p
            out.append(rec)
    return out


def five_trigger(mode: str, n: int) -> dict:
    """The model's own top-n triggers, 1/Sk38-weighted. A trigger whose
    forced shot gives no polygon is dropped and its slot not refilled
    (avachain dead-trigger rule)."""
    set_mode(mode)
    feats = G['feats']
    with contextlib.redirect_stdout(io.StringIO()):
        trig = select_triggers(feats, G['mel'], G['cluster_map'], G['dem'],
                               G['start_zone_mask'], n)
        polys = {c: _polygon(c)[0] for c in trig}
    live = [c for c in trig if polys[c] is not None]
    sk = sk38_column(feats)
    w = scoring.trigger_weights([feats.loc[c, sk] for c in live])
    s = scoring.weighted_event_score([polys[c] for c in live], w, G['observed'],
                                     G['start_zone_mask'], G['dem'].shape, G['transform'])
    s.update({'mode': mode, 'triggers': [int(c) for c in trig],
              'dead_triggers': [int(c) for c in trig if polys[c] is None],
              'per_trigger': [{'cid': int(c), 'weight': wi, 'area_m2': polys[c].area,
                               'wkt': polys[c].wkt,
                               'iou': scoring.iou(polys[c], G['observed'])}
                              for c, wi in zip(live, w)]})
    return s


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--data-dir', type=Path, default=DATA,
                    help='slope data dir (cluster map, DEM, start zone); default little_prof')
    ap.add_argument('--date', default=DATE, help='event date in the feature CSV names')
    ap.add_argument('--calibration', type=Path, default=DATA / 'calibration_2026-11-17.1.json')
    ap.add_argument('--features-dir', type=Path,
                    help='default: the preset\'s features_dir under data/little_prof')
    ap.add_argument('--observed', type=Path,
                    default=DATA / 'validation' / 'survey_release_area_20260118.geojson')
    ap.add_argument('--ski', type=lambda v: None if v == 'none' else Path(v),
                    default=DATA / 'validation' / 'ski_nontriggers_20260118.json',
                    help="ski non-trigger points JSON, or 'none'")
    ap.add_argument('--modes', nargs='+', default=None,
                    help="'shear', 'opening' or 'mixed:<energy|elliptic>:<stop|intact>' "
                         "(default: the preset's slab_crack_mode)")
    ap.add_argument('--extra-triggers', nargs='*', type=int, default=[],
                    help='further clusters scored as forced shots (e.g. 3869 for the SE lobe)')
    ap.add_argument('--no-inside', action='store_true', help='five-trigger score only')
    ap.add_argument('--subsample', type=int, default=1,
                    help='use every n-th inside cluster (avachain knob sweep: 4)')
    ap.add_argument('--edge-m', type=float, default=10.0)
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args(argv)
    t0 = time.time()
    cal = json.loads(a.calibration.read_text())
    apply_calibration(cal)
    modes = a.modes or [config.SLAB_CRACK_MODE]
    fdir = a.features_dir or a.data_dir / cal['features_dir']
    load_inputs(fdir, a.observed, a.ski, a.data_dir, a.date)
    a.out.mkdir(parents=True, exist_ok=True)
    print(f"calibration {cal.get('version')} features {fdir}; observed "
          f"{G['observed'].area:.0f} m2; modes {modes}", flush=True)

    five = [five_trigger(m, int(cal.get('n_triggers', config.N_TOP_TRIGGERS))) for m in modes]
    for s in five:
        print(f"  five-trigger {s['mode']}: {s['triggers']} dead {s['dead_triggers']}  "
              f"E[IoU] {s['expected_iou']:.3f}  E[area] {s['expected_area_m2']:.0f} m2  "
              f"AUC {s['auc']:.3f}  BSS {s['brier_skill']:.3f}", flush=True)

    rows, summary = [], []
    if not a.no_inside:
        inside = scoring.inside_clusters(G['cluster_map'], G['obs_mask'], G['start_zone_mask'],
                                         G['observed'], G['transform'])
        for r in inside:
            r['class'] = 'interior' if r['edge_dist_m'] >= a.edge_m else 'edge'
        ski = (scoring.clusters_at_points(G['ski']['points_lat_lon'], G['cluster_map'],
                                          G['transform']) if G['ski'] else [])
        print(f"inside-release triggers: {len(inside)}; ski clusters: {len(ski)}", flush=True)
        inside = inside[::a.subsample]
        jobs = [(m, r['cid'], r['class']) for m in modes for r in inside]
        jobs += [(m, c, 'ski') for m in modes for c in ski]
        jobs += [(m, c, 'extra') for m in modes for c in a.extra_triggers]
        with mproc.get_context('fork').Pool(a.jobs) as pool:
            for i, r in enumerate(pool.imap_unordered(run_one, jobs, chunksize=2), 1):
                rows.append(r)
                if i % 100 == 0 or i == len(jobs):
                    print(f"  {i}/{len(jobs)}  {time.time() - t0:.0f} s", flush=True)
        rows.sort(key=lambda r: (r['mode'], r['class'], r['cid']))
        summary = summarise(rows)
        edge = {r['cid']: r['edge_dist_m'] for r in inside}
        pd.DataFrame([{k: r.get(k) for k in ('mode', 'class', 'cid', 'outcome', 'trigger_s_t',
                                             'area_m2', 'iou', 'area_ungated_m2', 'iou_ungated')}
                      | {'edge_dist_m': edge.get(r['cid'])} for r in rows]
                     ).to_csv(a.out / 'triggers.csv', index=False)

    if G.get('freqs'):
        np.savez_compressed(a.out / 'inside_frequency.npz',
                            **{m.replace(':', '_'): f for m, f in G['freqs'].items()})
    meta = {'calibration': cal, 'features_dir': str(fdir), 'observed': str(a.observed),
            'data_dir': str(a.data_dir), 'date': a.date,
            'observed_m2': G['observed'].area, 'modes': modes, 'subsample': a.subsample,
            'edge_m': a.edge_m, 'seconds': round(time.time() - t0)}
    (a.out / 'summary.json').write_text(json.dumps(
        {'meta': meta, 'five_trigger': five, 'forced_shots': summary}, indent=1))

    f = lambda x, n=3: '-' if x is None else f'{x:.{n}f}'
    lines = ['| mode | five-trigger E[IoU] | E[area] m² | AUC | Brier skill | triggers |',
             '|---|---|---|---|---|---|']
    for s in five:
        lines.append(f"| {s['mode']} | {f(s['expected_iou'])} | {f(s['expected_area_m2'], 0)} | "
                     f"{f(s['auc'])} | {f(s['brier_skill'])} | {s['triggers']} |")
    if summary:
        lines += ['', '| mode | class | n | released | trigger gate / no prop. | '
                  'IoU median [IQR] | area median m² | consensus IoU | pairwise IoU |',
                  '|---|---|---|---|---|---|---|---|---|']
        for s in summary:
            lines.append(
                f"| {s['mode']} | {s['class']} | {s['n']} | {s['n_released']} | "
                f"{s['n_trigger_gate']} / {s['n_no_propagation']} | "
                f"{f(s.get('iou_median'))} [{f(s.get('iou_p25'))}–{f(s.get('iou_p75'))}] | "
                f"{f(s.get('area_median_m2'), 0)} | {f(s.get('consensus_iou'))} | "
                f"{f(s.get('pairwise_iou'))} |")
    (a.out / 'summary.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    print(f'done in {time.time() - t0:.0f} s -> {a.out}')


if __name__ == '__main__':
    sys.exit(main())
