#!/usr/bin/env python3
"""Release-scale slab tension on Jan 18: would an admitted cluster's slab
break once the release has grown, and where does that happen first?

The BFS tensile gate tests each cluster once, at the initiation scale
(cantilever over its own r_c, S_t = sigma_t / (3 rho g r_c^2 / h)). Nothing
re-tests an admitted cluster as the detached slab below it grows. This script
replays the BFS rings in order and, after each ring k, computes for every
admitted cluster n the axial tension from the detached slab below it:

    sigma(n, k) = g * drive * sum(rho h dl over the slab below) / h_n
                ~ g * drive * mean(rho h)_below * L_below / h_n
    drive       = sin(theta_below) - tan(phi) cos(theta_below)   (phi = PHI_DEG)

L_below is the down-fall-line extent of the slab admitted by ring k in the
SLAB_CRACK_BAND_M strip below n (n's own fall line), as in the crown load test
of slab_crack_lines. A strip of unit width held only by the crown:
the flanks carry no load, so this is an upper bound on crown tension.
sigma >= sigma_t (SLAB_TENSILE_SIGMA_COL) is a strength exceedance.
Also reported: L* = the detached length at which sigma reaches sigma_t.

Diagnostics only; the model is unchanged.

    env -u VIRTUAL_ENV uv run --frozen python examples/release_scale_tension_jan18.py \\
        --out outputs/score_jan18/tension_2026-11-19
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import multiprocessing as mproc
import sys
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


def _wrapped(reached, cracks, trigger, neighbours, centroids, fall_xy, props, slope_of,
             cluster_map, *a, **k):
    CAPTURE.append(dict(reached=set(reached), cracks=dict(cracks), trigger=trigger,
                        neighbours=neighbours, centroids=centroids, fall_xy=fall_xy,
                        props={c: dict(props.get(c, {})) for c in reached},
                        slope={c: slope_of(c) for c in reached},
                        fall_of=dict(k.get('fall_of') or {}),
                        ring=dict(k.get('front_ring') or {})))
    return _ORIG(reached, cracks, trigger, neighbours, centroids, fall_xy, props, slope_of,
                 cluster_map, *a, **k)


def prepare(cap):
    """Arrays for one captured shot and state(k): (sigma/sigma_t, L_below,
    drive, L*) of every admitted cluster once rings <= k are admitted."""
    ids = sorted(c for c in cap['reached'] if c in cap['centroids'] and c in cap['ring'])
    n = len(ids)
    xy = np.array([cap['centroids'][c] for c in ids])
    ring = np.array([cap['ring'][c] for c in ids])
    f = np.array([cap['fall_of'].get(c, cap['fall_xy']) for c in ids])
    p = cap['props']
    get = lambda k: np.array([p[c].get(k, np.nan) for c in ids], float)
    rho, h, st = get('slab_density'), get('slab_thickness'), get(config.SLAB_TENSILE_SIGMA_COL)
    rc = get('rc_wl')
    slope = np.radians([cap['slope'][c] for c in ids])
    dx = xy[None, :, 0] - xy[:, None, 0]
    dy = xy[None, :, 1] - xy[:, None, 1]
    al = dx * f[:, None, 0] + dy * f[:, None, 1]          # row i: along i's fall line
    cr = -dx * f[:, None, 1] + dy * f[:, None, 0]
    strip = (al > 0) & (np.abs(cr) <= config.SLAB_CRACK_BAND_M)
    rh = rho * h
    ok_j = np.isfinite(rh) & np.isfinite(slope)
    tan_phi = np.tan(np.radians(config.PHI_DEG))
    t = int(np.where(np.array(ids) == cap['trigger'])[0][0])
    t_xy, t_f = xy[t], np.array(cap['fall_xy'])

    def state(k):
        m = strip & (ring[None, :] <= k) & ok_j[None, :]
        cnt = m.sum(1)
        L = np.where(cnt > 0, np.where(m, al, 0).max(1), 0.0)
        with np.errstate(invalid='ignore', divide='ignore'):
            th = np.where(cnt > 0, (np.where(m, slope[None, :], 0).sum(1)) / cnt, np.nan)
            mrh = np.where(cnt > 0, (np.where(m, rh[None, :], 0).sum(1)) / cnt, np.nan)
            drive = np.sin(th) - tan_phi * np.cos(th)
            sig = np.where((cnt > 0) & (drive > 0), G_ACC * drive * mrh * L / h, 0.0)
            lstar = np.where(drive > 0, st * h / (G_ACC * drive * mrh), np.nan)
        return sig / st, L, drive, lstar

    return dict(ids=ids, n=n, xy=xy, ring=ring, f=f, rho=rho, h=h, st=st, rc=rc,
                slope=slope, t=t, t_xy=t_xy, t_f=t_f, state=state)


def analyse(cap):
    P = prepare(cap)
    ids, n, xy, ring, rho, h, st, rc, slope = (P[k] for k in (
        'ids', 'n', 'xy', 'ring', 'rho', 'h', 'st', 'rc', 'slope'))
    t_xy, t_f, state = P['t_xy'], P['t_f'], P['state']
    rows, first = [], None
    for k in range(1, int(ring.max()) + 1):
        ratio, L, drive, _ = state(k)
        live = (ring <= k) & np.isfinite(ratio)
        hit = live & (ratio >= 1.0)
        if hit.any() and first is None:
            i = int(np.argmax(np.where(hit, ratio, -np.inf)))
            first = dict(ring=k, cid=ids[i], ratio=float(ratio[i]), L=float(L[i]), i=i)
        rows.append(dict(ring=k, admitted=int((ring <= k).sum()), exceed=int(hit.sum())))
    ratio, L, drive, lstar = state(int(ring.max()))
    gate_s = st / (3 * rho * G_ACC * np.maximum(rc, config.RC_MIN_M) ** 2 / h)

    # gate cracks: ring = 1 + the lowest ring of an admitted neighbour
    adj = {}
    for c, nb in cap['neighbours'].items():
        for u in nb:
            adj.setdefault(c, set()).add(u); adj.setdefault(u, set()).add(c)
    crk = []
    for c, o in cap['cracks'].items():
        rs = [cap['ring'][u] for u in adj.get(c, ()) if u in cap['ring']]
        if rs and c in cap['centroids']:
            crk.append((c, o, min(rs) + 1, cap['centroids'][c]))
    out = dict(trigger=cap['trigger'], n_admitted=n, n_rings=int(ring.max()),
               n_gate_cracks=len(crk),
               first_gate_ring=min((r for _, _, r, _ in crk), default=None),
               frac_exceed_final=float(np.nanmean(ratio >= 1.0)),
               ratio_median_final=float(np.nanmedian(ratio)),
               lstar_median=float(np.nanmedian(lstar)),
               lstar_p10=float(np.nanpercentile(lstar, 10)),
               lstar_p90=float(np.nanpercentile(lstar, 90)),
               gate_S_median=float(np.nanmedian(gate_s)),
               frac_no_drive=float(np.nanmean(~(drive > 0))))
    if first:
        v = xy[first['i']] - t_xy
        d = float(np.hypot(*v))
        along = float(v @ t_f)
        out.update(first_ring=first['ring'], first_cid=first['cid'], first_L=first['L'],
                   first_dist=d,
                   first_where=('upslope' if along < -0.707 * d else
                                'downslope' if along > 0.707 * d else 'cross'),
                   first_to_gate_crack=(min(float(np.hypot(*(np.array(cx) - xy[first['i']])))
                                            for _, _, _, cx in crk) if crk else None))
    cells = dict(trigger=cap['trigger'], ids=[int(c) for c in ids],
                 ratio=[None if not np.isfinite(r) else float(r) for r in ratio],
                 ring=ring.tolist(), first=out.get('first_cid'),
                 slope_deg=np.degrees(slope).tolist(),
                 cracks=[(int(c), o, int(r)) for c, o, r, _ in crk],
                 per_ring=rows)
    return out, cells


def shot(job):
    cid, label = job
    CAPTURE.clear()
    with contextlib.redirect_stdout(io.StringIO()):
        sj._polygon(cid)
    res = []
    for cap in CAPTURE:
        o, c = analyse(cap)
        o['set'] = label
        c['set'] = label
        res.append((o, c))
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--calibration', type=Path, default=sj.DATA / 'calibration_2026-11-17.1.json')
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--jobs', type=int, default=16)
    a = ap.parse_args(argv)
    cal = json.loads(a.calibration.read_text())
    sj.apply_calibration(cal)
    # 'opening' only so the BFS hands the local fall lines to the crack lines;
    # the admitted region and rings do not depend on the slab-crack mode.
    sj.set_mode('opening')
    sj.load_inputs(sj.DATA / cal['features_dir'],
                   sj.DATA / 'validation' / 'survey_release_area_20260118.geojson', None)
    _rg.slab_crack_lines = _wrapped
    with contextlib.redirect_stdout(io.StringIO()):
        trig = select_triggers(sj.G['feats'], sj.G['mel'], sj.G['cluster_map'], sj.G['dem'],
                               sj.G['start_zone_mask'], int(cal.get('n_triggers', 5)))
    inside = scoring.inside_clusters(sj.G['cluster_map'], sj.G['obs_mask'],
                                     sj.G['start_zone_mask'], sj.G['observed'], sj.G['transform'])
    jobs = [(int(c), 'five') for c in trig] + [(int(r['cid']), 'inside') for r in inside]
    outs, cells = [], []
    with mproc.get_context('fork').Pool(a.jobs) as pool:
        for res in pool.imap_unordered(shot, jobs, chunksize=2):
            for o, c in res:
                outs.append(o); cells.append(c)
    a.out.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(outs)
    df.to_csv(a.out / 'tension.csv', index=False)
    (a.out / 'tension_cells.json').write_text(json.dumps([c for c in cells if c['set'] == 'five']))
    (a.out / 'tension_cells_inside.json').write_text(json.dumps(
        [{k: c[k] for k in ('trigger', 'ids', 'ratio', 'slope_deg')}
         for c in cells if c['set'] == 'inside']))
    (a.out / 'meta.json').write_text(json.dumps({'calibration': cal, 'five': [int(c) for c in trig]}))
    print(f'{len(df)} shots -> {a.out}')


if __name__ == '__main__':
    sys.exit(main())
