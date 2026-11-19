#!/usr/bin/env python3
"""Release animation for one Jan 18 trigger: the weak-layer crack front (BFS
rings), the slab tension front and slab motion behind it, the tensile-gate
breaks, then the slab-crack lines running and the severed slab dropping out.

Three phases, all on the same map (left: fracture, right: slab tension):

  1. BFS, one frame per ring k. Left: admitted clusters coloured by arrival
     ring, the newest ring bright; tensile-gate cracks (crown red, flank
     green) appear at their ring (1 + the earliest admitted neighbour, as in
     release_scale_tension_jan18.py). Right: sigma/sigma_t from the detached
     slab below each admitted cluster at ring k (unit-width strip held by the
     crown only: an upper bound, see release_scale_tension_jan18.py); x marks
     sigma >= sigma_t (the slab would break there; the model does not re-test,
     so this is diagnostic). Arrows: the local fall line, length ~ drive =
     sin(theta) - tan(PHI_DEG) cos(theta) (crack-tip friction, 27 deg); black
     where the slope exceeds atan(SEVERED_SLIDE_MU) = 31.8 deg (the slab slides
     once it is cut free), grey where it is held.
  2. Slab-crack lines, one frame per line step (all lines step together; in
     the model they run after the BFS and their order does not matter).
     Crown lines orange, flank lines purple; 'intact' / toughness cells that
     do not break are outlined.
  3. The final release polygon (blue outline, filled); everything admitted
     but not released (severed, 'intact' / toughness cells, and whatever they
     cut off from the trigger) is light grey.

Diagnostics only; the model is unchanged. Time is the BFS ring / line step,
not seconds.

    env -u VIRTUAL_ENV uv run --frozen python examples/animate_release_jan18.py \\
        --mode opening --out outputs/score_jan18/animation_2026-11-19
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                   # noqa: E402
from matplotlib import animation, colors                          # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score_jan18 as sj                                          # noqa: E402
import release_scale_tension_jan18 as rst                         # noqa: E402
from release_areas import config                                  # noqa: E402
from release_areas import release_geometry as _rg                 # noqa: E402
from release_areas.generate_scenarios import select_triggers     # noqa: E402

LINES: list = []
OUT: dict = {}


def _wrapped(*a, **k):
    LINES.clear()
    k['lines_out'] = LINES
    res = rst._wrapped(*a, **k)
    OUT['severed'] = set(res[0])
    OUT['intact'] = set(k.get('intact_out') or ())
    return res


def run_shot(cid):
    rst.CAPTURE.clear(); OUT.clear()
    with contextlib.redirect_stdout(io.StringIO()):
        poly, state = sj._polygon(cid)
    if not rst.CAPTURE:
        return None
    return dict(cap=rst.CAPTURE[-1], lines=[dict(l, cells=list(l['cells'])) for l in LINES],
                severed=OUT.get('severed', set()), intact=OUT.get('intact', set()),
                poly=poly)


def hillshade(dem):
    gy, gx = np.gradient(np.nan_to_num(dem, nan=np.nanmean(dem)))
    sl = np.pi / 2 - np.arctan(np.hypot(gx, gy)); asp = np.arctan2(-gx, gy)
    return (np.sin(np.radians(45)) * np.sin(sl)
            + np.cos(np.radians(45)) * np.cos(sl) * np.cos(np.radians(315) - asp))


def render(shot, cid, mode, out_path, fps):
    cap = shot['cap']
    P = rst.prepare(cap)
    ids, ring, xy, f, slope = P['ids'], P['ring'], P['xy'], P['f'], P['slope']
    K = int(ring.max())
    cmap_r = sj.G['cluster_map'].astype(np.int64); tf = sj.G['transform']
    cmap_r[cmap_r < 0] = 0
    nlut = int(cmap_r.max()) + 1
    # view: the admitted region and the survey polygon, 15 m margin
    obs = sj.G['observed']
    ox0, oy0, ox1, oy1 = obs.bounds
    rr, cc = np.where(np.isin(cmap_r, ids))
    pad = 15
    r0 = max(min(rr.min(), int((oy1 - tf.f) / tf.e)) - pad, 0)
    r1 = min(max(rr.max(), int((oy0 - tf.f) / tf.e)) + pad, cmap_r.shape[0])
    c0 = max(min(cc.min(), int((ox0 - tf.c) / tf.a)) - pad, 0)
    c1 = min(max(cc.max(), int((ox1 - tf.c) / tf.a)) + pad, cmap_r.shape[1])
    crop = cmap_r[r0:r1, c0:c1]
    ext = (tf.c + c0 * tf.a, tf.c + c1 * tf.a, tf.f + r1 * tf.e, tf.f + r0 * tf.e)
    hs = hillshade(sj.G['dem'])[r0:r1, c0:c1]
    gx, gy = np.meshgrid(tf.c + (np.arange(c0, c1) + .5) * tf.a,
                         tf.f + (np.arange(r0, r1) + .5) * tf.e)
    obs_geoms = getattr(obs, 'geoms', [obs])

    # final release, rasterised on the crop: the polygon is what the model
    # keeps after severed, intact and cut-off clusters are dropped
    rel_mask = None
    if shot['poly'] is not None:
        from rasterio import features as rfeat
        from affine import Affine
        rel_mask = rfeat.geometry_mask([shot['poly']], out_shape=crop.shape,
                                       transform=tf * Affine.translation(c0, r0),
                                       invert=True)

    def img(vals: dict):
        lut = np.full(nlut, np.nan)
        for c, v in vals.items():
            if 0 < c < nlut:
                lut[c] = v
        return lut[crop]

    # gate cracks and their ring
    adj = {}
    for c, nb in cap['neighbours'].items():
        for u in nb:
            adj.setdefault(c, set()).add(u); adj.setdefault(u, set()).add(c)
    crack_ring = {}
    for c, o in cap['cracks'].items():
        rs = [cap['ring'][u] for u in adj.get(c, ()) if u in cap['ring']]
        if rs:
            crack_ring[c] = (o, min(rs) + 1)
    lines = [l for l in shot['lines'] if l['cells']]
    n_line = max((len(l['cells']) for l in lines), default=0)
    tan_mu = config.SEVERED_SLIDE_MU
    slides = slope > np.arctan(tan_mu)
    states = [P['state'](k) for k in range(K + 1)]

    frames = ([('bfs', k) for k in range(K + 1)] + [('bfs', K)] * 3
              + [('line', i) for i in range(1, n_line + 1)] + [('line', n_line)] * 2
              + [('final', 0)] * 8)

    fig = plt.figure(figsize=(13, 7.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[4.2, 1], hspace=0.18, wspace=0.04)
    axL, axR = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axT = fig.add_subplot(gs[1, :])
    ring_cmap = plt.get_cmap('Blues')
    ten_cmap = plt.get_cmap('RdBu_r').copy()
    # curves for the strip below
    adm = [int((ring <= k).sum()) for k in range(K + 1)]
    exc = [int(((ring <= k) & (states[k][0] >= 1)).sum()) for k in range(K + 1)]
    gcr = [sum(r <= k for _, r in crack_ring.values()) for k in range(K + 1)]
    xy_t = P['t_xy']

    def draw(fr):
        phase, k = fr
        kk = K if phase != 'bfs' else k
        for ax in (axL, axR):
            ax.clear()
            ax.imshow(hs, extent=ext, cmap='gray', vmin=-.2, vmax=1.1)
            ax.set_xticks([]); ax.set_yticks([])
        live = ring <= kk
        final = phase == 'final'
        sev = shot['severed'] if final else set()
        # left: rings
        vals = {c: ring[i] for i, c in enumerate(ids) if live[i] and c not in sev}
        axL.imshow(img(vals), extent=ext, cmap=ring_cmap, vmin=-.25 * K, vmax=max(K, 1),
                   interpolation='nearest', alpha=.85)
        if phase == 'bfs':
            axL.imshow(img({c: 1 for i, c in enumerate(ids) if ring[i] == kk}), extent=ext,
                       cmap=colors.ListedColormap(['#ffff33']), interpolation='nearest')
        if final:
            # admitted but not released: severed, intact, or cut off from the trigger
            dropped = np.isin(crop, ids) & ~(rel_mask if rel_mask is not None else False)
            g = np.where(dropped, 1.0, np.nan)
            axL.imshow(g, extent=ext, interpolation='nearest',
                       cmap=colors.ListedColormap(['#c8c8c8']), alpha=.95)
        cr = {c: (1 if o == 'crown' else 0) for c, (o, r) in crack_ring.items() if r <= kk}
        axL.imshow(img(cr), extent=ext, interpolation='nearest', vmin=0, vmax=1,
                   cmap=colors.ListedColormap(['#1a9641', '#d7191c']))
        if phase != 'bfs':
            step = n_line if phase == 'final' else k
            lv = {}
            for l in lines:
                for c in l['cells'][:step]:
                    lv[c] = 1 if l['orient'] == 'crown' else 0
            axL.imshow(img(lv), extent=ext, interpolation='nearest', vmin=0, vmax=1,
                       cmap=colors.ListedColormap(['#e7298a', '#fd8d3c']))
            if shot['intact']:
                m = np.isin(crop, list(shot['intact'])).astype(float)
                axL.contour(gx, gy, m, levels=[.5], colors='red', linewidths=.9)
        # right: tension
        ratio, _, drive, _ = states[kk]
        tv = {c: np.log10(max(ratio[i], 1e-3)) for i, c in enumerate(ids)
              if live[i] and np.isfinite(ratio[i]) and c not in sev}
        timg = img(tv)
        if final and rel_mask is not None:
            timg = np.where(rel_mask, timg, np.nan)
        axR.imshow(timg, extent=ext, cmap=ten_cmap, vmin=-1, vmax=1,
                   interpolation='nearest', alpha=.9)
        keep = np.ones(len(ids), bool)
        if final:
            keep = np.array([c not in shot['severed'] for c in ids])
            if rel_mask is not None:
                inside = set(np.unique(crop[rel_mask]).tolist())
                keep = np.array([c in inside for c in ids])
        live = live & keep
        brk = live & (ratio >= 1)
        axR.plot(xy[brk, 0], xy[brk, 1], 'x', color='#444444', ms=2, mew=.5)
        mv = live & (drive > 0) & np.isfinite(drive)
        sub = np.where(mv)[0][::4]
        sc = 40 * np.clip(drive[sub], 0, None)            # m per unit drive
        for sel, col in ((slides[sub], 'k'), (~slides[sub], '#777777')):
            s2 = sub[sel]
            if len(s2):
                axR.quiver(xy[s2, 0], xy[s2, 1], f[s2, 0] * sc[sel], f[s2, 1] * sc[sel],
                           color=col, angles='xy', scale_units='xy', scale=1, width=.003,
                           headwidth=3.5, headlength=4)
        for ax in (axL, axR):
            for g in obs_geoms:
                ax.plot(*g.exterior.xy, color='k', lw=1.0, ls='--')
            ax.plot(*xy_t, '*', color='gold', mec='k', ms=13)
            if phase == 'final' and shot['poly'] is not None:
                for g in getattr(shot['poly'], 'geoms', [shot['poly']]):
                    ax.plot(*g.exterior.xy, color='#0571b0', lw=2)
            ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
        title = {'bfs': f'BFS ring {kk}/{K}: weak-layer crack front',
                 'line': f'slab-crack lines, step {k}/{n_line}',
                 'final': (f"release {shot['poly'].area if shot['poly'] is not None else 0:.0f} m²"
                           f" (blue outline); light grey = admitted but dropped "
                           f"({len(shot['severed'])} severed, {len(shot['intact'])} intact, "
                           f"+ cut off)")}[phase]
        axL.set_title(title + '\nblue = admitted (darker = later ring), yellow = new ring\n'
                      'gate crack: red crown, green flank; line: orange crown, pink flank; '
                      'red outline = intact', fontsize=8)
        axR.set_title(f'slab tension σ/σt from the detached slab below (log; blue < 1 < red), '
                      f'x = σ ≥ σt\narrows: fall line × drive (φ {config.PHI_DEG:g}°), '
                      f'40 m per unit drive\nblack = slope > '
                      f'{np.degrees(np.arctan(tan_mu)):.1f}° (slides once free), grey = held',
                      fontsize=8)
        axT.clear()
        rk = np.arange(K + 1)
        axT.plot(rk, adm, color='#2c7fb8', label='admitted clusters')
        axT.plot(rk, exc, color='#d7191c', label='σ ≥ σt (release scale)')
        axT.plot(rk, np.array(gcr) * 10, color='#1a9641', label='gate cracks × 10')
        axT.axvline(kk, color='k', lw=.8)
        axT.set_xlabel('BFS ring'); axT.legend(fontsize=7.5, loc='upper left'); axT.grid(alpha=.3)
        fig.suptitle(f'Jan 18 2026, trigger {cid}, slab-crack mode {mode} '
                     f'(cal 2026-11-17.1). Time axis = BFS ring, then line step; dashed = survey',
                     fontsize=10)

    ani = animation.FuncAnimation(fig, draw, frames=frames, interval=1000 / fps)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ani.save(out_path, writer=animation.FFMpegWriter(fps=fps, bitrate=2400), dpi=110)
    draw(frames[-1]); fig.savefig(out_path.with_suffix('.png'), dpi=120)
    plt.close(fig)
    return dict(trigger=cid, rings=K, line_steps=n_line, gate_cracks=len(crack_ring),
                severed=len(shot['severed']), intact=len(shot['intact']),
                exceed_final=exc[-1], admitted=adm[-1])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--calibration', type=Path, default=sj.DATA / 'calibration_2026-11-17.1.json')
    ap.add_argument('--mode', default='opening', help='score_jan18 variant string')
    ap.add_argument('--triggers', type=int, nargs='*',
                    help='cluster ids (default: the five selected triggers)')
    ap.add_argument('--fps', type=float, default=4)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args(argv)
    cal = json.loads(a.calibration.read_text())
    sj.apply_calibration(cal)
    sj.set_mode(a.mode)
    sj.load_inputs(sj.DATA / cal['features_dir'],
                   sj.DATA / 'validation' / 'survey_release_area_20260118.geojson', None)
    _rg.slab_crack_lines = _wrapped
    trig = a.triggers
    if not trig:
        with contextlib.redirect_stdout(io.StringIO()):
            trig = [int(c) for c in select_triggers(
                sj.G['feats'], sj.G['mel'], sj.G['cluster_map'], sj.G['dem'],
                sj.G['start_zone_mask'], int(cal.get('n_triggers', 5)))]
    tag = a.mode.replace(':', '_')
    for cid in trig:
        shot = run_shot(cid)
        if shot is None:
            print(f'{cid}: no slab-crack call (no propagation)'); continue
        info = render(shot, cid, a.mode, a.out / f'release_{cid}_{tag}.mp4', a.fps)
        print(json.dumps(info))


if __name__ == '__main__':
    sys.exit(main())
