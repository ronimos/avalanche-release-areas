#!/usr/bin/env python3
"""Flank-line angle alpha on Jan 18: how far the release's crack front at a
flank cell is from the cell's own fall line.

'mixed' splits the relative motion across a flank into opening and shear on
axes built from the cell's fall line (t = fall line, e = cross axis). If the
weak-layer front reaches the flank on a curve (trigger near the flank), the
line the slab breaks along is tilted from the fall line by alpha, and the
same motion would split differently. This script measures alpha for every
flank step the slab-crack lines test, three ways:

  front   front normal = gradient of the BFS arrival ring, least-squares plane
          through the reached clusters within --radius m of the cell;
          alpha = asin(|normal . fall|) (0 = front tangent on the fall line)
  radial  circular front about the trigger: normal = (cell - trigger) / r
  step    the line's own step (previous cell -> cell) against the fall line

and recomputes the split on the front axes (d_open', d_shear', and the
'energy' / 'elliptic' outcome). Diagnostics only; the model is unchanged.

Five-trigger set for every mode; --inside adds every inside-release forced
shot for one mode (default mixed:energy:intact, whose lines run to their
geometric end, so every flank cell is tested).

    env -u VIRTUAL_ENV uv run --frozen python examples/alpha_jan18.py \\
        --run outputs/score_jan18/mixed_2026-11-18 --inside
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import multiprocessing as mproc
import sys
from collections import deque
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score_jan18 as sj                                     # noqa: E402
from release_areas import config                             # noqa: E402
from release_areas import release_geometry as _rg            # noqa: E402
from release_areas import scoring                            # noqa: E402
from release_areas.generate_scenarios import select_triggers  # noqa: E402

G_ACC = 9.81
CAPTURE: list = []
_ORIG = _rg.slab_crack_lines


def _wrapped(reached, cracks, trigger_cluster_id, neighbours, centroids, *a, **k):
    trace = []
    k['trace_out'] = trace
    if k.get('intact_out') is None:
        k['intact_out'] = set()
    cut = k.get('cut_out')
    if cut is None:
        k['cut_out'] = cut = set()
    reached0 = set(reached)
    severed, stats = _ORIG(reached, cracks, trigger_cluster_id, neighbours, centroids, *a, **k)
    CAPTURE.append(dict(reached=reached0, cracks=dict(cracks), trigger=trigger_cluster_id,
                        neighbours=neighbours, centroids=centroids, trace=trace,
                        intact=set(k['intact_out']), cut=set(cut), severed=set(severed),
                        stats=dict(stats)))
    return severed, stats


def rings(trigger, reached, neighbours):
    """BFS hop count from the trigger through admitted clusters: the legacy
    BFS's arrival ring (anisotropic rings are off)."""
    ring = {trigger: 0}
    q = deque([trigger])
    while q:
        c = q.popleft()
        for n in neighbours.get(c, ()):
            if n in reached and n not in ring:
                ring[n] = ring[c] + 1
                q.append(n)
    return ring


def front_normal(xy, ring, centroids, radius):
    pts = [(centroids[c][0] - xy[0], centroids[c][1] - xy[1], r)
           for c, r in ring.items() if c in centroids
           and np.hypot(centroids[c][0] - xy[0], centroids[c][1] - xy[1]) <= radius]
    if len(pts) < 4:
        return None
    A = np.array([[p[0], p[1], 1.0] for p in pts])
    b = np.array([p[2] for p in pts], float)
    if np.linalg.matrix_rank(A) < 3:
        return None
    g = np.linalg.lstsq(A, b, rcond=None)[0][:2]
    nrm = np.hypot(*g)
    return None if nrm < 1e-9 else g / nrm


def angle_to_fall(v, f):
    """alpha in deg of the line perpendicular to unit normal v, against fall f."""
    return float(np.degrees(np.arcsin(min(abs(v[0] * f[0] + v[1] * f[1]), 1.0))))


def rows_from(cap, mode, label, radii):
    out = []
    ring = rings(cap['trigger'], cap['reached'], cap['neighbours'])
    cen = cap['centroids']
    txy = np.array(cen[cap['trigger']])
    adj = {}
    for c, nb in cap['neighbours'].items():
        for u in nb:
            adj.setdefault(c, set()).add(u)
            adj.setdefault(u, set()).add(c)
    reached = cap['reached']
    for s in cap['trace']:
        n = s['cid']
        xy = np.array(cen[n])
        f = np.array([s['fall_e'], s['fall_n']])
        e = np.array([-f[1], f[0]])
        dF = np.array([s['dF_e'], s['dF_n']])
        r = xy - txy
        rec = dict(mode=mode, set=label, trigger=cap['trigger'], cid=n, sgn=s['sgn'],
                   start=s['start'], alen=s['alen'], dist_trigger=float(np.hypot(*r)),
                   ring=ring.get(n), attached=(s['n_plus'] == 0 or s['n_minus'] == 0),
                   d_open=s.get('d_open'), d_shear=s.get('d_shear'), breaks=s.get('breaks'),
                   axes=s.get('axes'),
                   intact=n in cap['intact'], kic=s['kic'])
        rec['edge'] = any(u not in reached for u in adj.get(n, ()))
        rec['alpha_radial'] = angle_to_fall(r / max(np.hypot(*r), 1e-9), f)
        st = xy - np.array(cen[s['prev']])
        rec['alpha_step'] = float(np.degrees(np.arccos(min(abs(st @ f) / max(np.hypot(*st), 1e-9), 1.0))))
        for R in radii:
            nu = front_normal(xy, ring, cen, R)
            rec[f'alpha_front_{R:g}'] = np.nan if nu is None else angle_to_fall(nu, f)
            if nu is None:
                continue
            if R == radii[0]:
                # the line's own step against the front tangent: 0 = the line
                # runs along the front (an edge), 90 = across it (interior)
                rec['line_vs_front'] = float(np.degrees(np.arcsin(
                    min(abs(st @ nu) / max(np.hypot(*st), 1e-9), 1.0))))
            nu = nu if nu @ e >= 0 else -nu               # + side, as e
            tau = np.array([nu[1], -nu[0]])
            do2, ds2 = max(dF @ nu, 0.0), abs(dF @ tau)
            rec[f'd_open_front_{R:g}'], rec[f'd_shear_front_{R:g}'] = do2, ds2
            kic = s['kic']
            if np.isfinite(kic) and kic > 0 and s['drive'] > 0:
                sc = (config.TOUGHNESS_EDGE_F * s['rho'] * G_ACC * s['L'] * s['drive'] / 2.0
                      * np.sqrt(np.pi * max(s['alen'], 0.0)))
                d_o = max(dF @ e, 0.0)
                d_s = abs(dF @ f)
                rec['brk_energy'] = (sc * d_o) ** 2 + (sc * d_s) ** 2 >= kic ** 2
                rec['brk_elliptic'] = ((sc * d_o / kic) ** 2
                                       + (sc * d_s / (config.K_IIC_OVER_K_IC * kic)) ** 2) >= 1
                rec[f'brk_energy_front_{R:g}'] = (sc * do2) ** 2 + (sc * ds2) ** 2 >= kic ** 2
                rec[f'brk_elliptic_front_{R:g}'] = (
                    (sc * do2 / kic) ** 2 + (sc * ds2 / (config.K_IIC_OVER_K_IC * kic)) ** 2) >= 1
        out.append(rec)
    return out


def shot(job):
    mode, cid, label, radii = job
    sj.set_mode(mode)
    CAPTURE.clear()
    with contextlib.redirect_stdout(io.StringIO()):
        poly, _ = sj._polygon(cid)
    rows = [r for cap in CAPTURE for r in rows_from(cap, mode, label, radii)]
    cells = [dict(mode=mode, set=label, trigger=cid,
                  cut=sorted(int(c) for c in cap['cut']),
                  intact=sorted(int(c) for c in cap['intact']),
                  severed=sorted(int(c) for c in cap['severed']),
                  cracks={int(c): o for c, o in cap['cracks'].items()},
                  stats={k: int(v) for k, v in cap['stats'].items()},
                  area=None if poly is None else poly.area) for cap in CAPTURE]
    return rows, cells


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--run', type=Path, required=True, help='score_jan18 output dir')
    ap.add_argument('--radii', type=float, nargs='+', default=[10.0, 5.0, 15.0])
    ap.add_argument('--inside', action='store_true')
    ap.add_argument('--inside-mode', default='mixed:energy:intact')
    ap.add_argument('--jobs', type=int, default=16)
    a = ap.parse_args(argv)
    summ = json.loads((a.run / 'summary.json').read_text())
    cal = summ['meta']['calibration']
    sj.apply_calibration(cal)
    meta = summ['meta']
    sj.load_inputs(Path(meta['features_dir']), Path(meta['observed']), None,
                   Path(meta.get('data_dir', sj.DATA)), meta.get('date', sj.DATE))
    _rg.slab_crack_lines = _wrapped
    modes = summ['meta']['modes']
    jobs = []
    for m in modes:
        sj.set_mode(m)
        with contextlib.redirect_stdout(io.StringIO()):
            trig = select_triggers(sj.G['feats'], sj.G['mel'], sj.G['cluster_map'], sj.G['dem'],
                                   sj.G['start_zone_mask'], int(cal.get('n_triggers', 5)))
        jobs += [(m, int(c), 'five', a.radii) for c in trig]
    if a.inside:
        inside = scoring.inside_clusters(sj.G['cluster_map'], sj.G['obs_mask'],
                                         sj.G['start_zone_mask'], sj.G['observed'],
                                         sj.G['transform'])
        jobs += [(a.inside_mode, int(r['cid']), 'inside', a.radii) for r in inside]
    rows, cells = [], []
    with mproc.get_context('fork').Pool(a.jobs) as pool:
        for r, c in pool.imap_unordered(shot, jobs, chunksize=2):
            rows += r
            cells += c
    df = pd.DataFrame(rows)
    df.to_csv(a.run / 'alpha_trace.csv', index=False)
    (a.run / 'alpha_cells.json').write_text(json.dumps(cells))
    print(f'{len(df)} flank steps from {len(jobs)} shots -> {a.run}/alpha_trace.csv')


if __name__ == '__main__':
    sys.exit(main())
