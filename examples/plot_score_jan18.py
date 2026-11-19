#!/usr/bin/env python3
"""Map of a score_jan18.py run: one column per slab-crack variant.

Row 1  five-trigger release probability (1/Sk38-weighted polygons)
Row 2  inside-release forced shots: fraction of the released shots covering
       each pixel
Row 3  row 2 minus the first variant (default: shear): blue = covered more
       often, red = less often
Row 4  (when examples/alpha_jan18.py has run) the slab-crack lines of the
       top-weighted trigger: line cells coloured by alpha, the angle between
       the BFS front's tangent and the cell's fall line; red outline = removed
       as intact, black x = a 'stop' line ended there, grey = severed; the
       trigger's final polygon in blue

Observed (survey) release in black, start zone dashed grey, on a hillshade.

    env -u VIRTUAL_ENV uv run --frozen python examples/plot_score_jan18.py \\
        outputs/score_jan18/mixed_2026-11-18
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                              # noqa: E402
import numpy as np                                           # noqa: E402
import rasterio                                              # noqa: E402
from matplotlib.colors import LinearSegmentedColormap        # noqa: E402
from shapely import wkt                                      # noqa: E402

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / 'data' / 'little_prof'

SEQ = LinearSegmentedColormap.from_list('seq', ['#e8eef7', '#9cb9df', '#3b6fb6', '#0b2e6b'])
DIV = LinearSegmentedColormap.from_list('div', ['#b2182b', '#ef8a62', '#e6e6e6', '#67a9cf', '#2166ac'])


def hillshade(dem, az=315.0, alt=45.0):
    gy, gx = np.gradient(np.nan_to_num(dem, nan=np.nanmean(dem)))
    slope = np.pi / 2 - np.arctan(np.hypot(gx, gy))
    aspect = np.arctan2(-gx, gy)
    a, z = np.radians(az), np.radians(alt)
    return np.sin(z) * np.sin(slope) + np.cos(z) * np.cos(slope) * np.cos(a - aspect)


def load_alpha(run, data=DATA):
    tr, ce = run / 'alpha_trace.csv', run / 'alpha_cells.json'
    if not (tr.exists() and ce.exists()):
        return None
    import pandas as pd
    from release_areas.release_geometry import cluster_pixel_index
    with rasterio.open(data / 'spatial' / 'cluster_map.tif') as src:
        cmap = src.read(1).astype(float)
    return dict(trace=pd.read_csv(tr, low_memory=False),
                cells=[c for c in json.loads(ce.read_text()) if c['set'] == 'five'],
                px=cluster_pixel_index(cmap))


def alpha_panel(ax, mode, s, A, hs, ext, crop, sz, obs, tf, shape):
    """Row 4: the slab-crack lines of the mode's top-weighted trigger."""
    top = max(s['per_trigger'], key=lambda p: p['weight'])
    cid = top['cid']
    cell = next((c for c in A['cells'] if c['mode'] == mode and c['trigger'] == cid), None)
    tr = A['trace']
    tr = tr[(tr['set'] == 'five') & (tr['mode'] == mode) & (tr['trigger'] == cid)]
    px = A['px']

    def paint(ids, vals=None):
        a = np.full(shape, np.nan)
        for k, c in enumerate(ids):
            p = px.get(int(c))
            if p is not None:
                a[p[0], p[1]] = 1.0 if vals is None else vals[k]
        return crop(a)

    ax.imshow(hs, extent=ext, cmap='gray', vmin=-0.2, vmax=1.1, interpolation='bilinear')
    if cell:
        sev = paint(cell['severed'])
        ax.imshow(sev, extent=ext, cmap=LinearSegmentedColormap.from_list('g', ['#9a9a9a'] * 2),
                  alpha=0.75, interpolation='nearest')
    # alpha on every tested flank cell (median over its steps); crown cells grey-blue
    fl = tr.groupby('cid')['alpha_front_10'].median()
    im = ax.imshow(paint(fl.index, fl.values), extent=ext, cmap='viridis', vmin=0, vmax=90,
                   interpolation='nearest')
    if cell:
        crowns = [c for c in cell['cut'] if c not in set(fl.index)]
        ax.imshow(paint(crowns), extent=ext, cmap=LinearSegmentedColormap.from_list(
            'c', ['#f28e2b'] * 2), interpolation='nearest')
        if cell['intact']:
            ax.contour(np.nan_to_num(paint(cell['intact'])), levels=[0.5],
                       extent=(ext[0], ext[1], ext[3], ext[2]), origin='upper',
                       colors='#d62728', linewidths=0.9)
    if mode.endswith(':stop'):
        ends = tr[tr['breaks'].astype(str) == 'False']
        xs = [tf.c + (A['px'][int(c)][1].mean() + 0.5) * tf.a for c in ends['cid'] if int(c) in px]
        ys = [tf.f + (A['px'][int(c)][0].mean() + 0.5) * tf.e for c in ends['cid'] if int(c) in px]
        ax.plot(xs, ys, 'x', color='black', ms=4.5, mew=1.1)
    poly = wkt.loads(top['wkt'])
    for g in getattr(poly, 'geoms', [poly]):
        ax.plot(*g.exterior.xy, color='#1f4e9c', lw=0.9)
    ax.plot(*obs.exterior.xy, color='black', lw=1.2)
    ax.contour(crop(sz).astype(float), levels=[0.5], extent=(ext[0], ext[1], ext[3], ext[2]),
               colors='#555555', linewidths=0.8, linestyles='--', origin='upper')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
    e = tr[tr['edge'].astype(bool)]
    n_fail = int((tr['breaks'].astype(str) == 'False').sum())
    st = cell['stats'] if cell else {}
    ax.set_title(f"trigger {cid}: {len(tr)} flank steps"
                 + (f", {n_fail} fail" if 'mixed' in mode else '') + "\n"
                 f"line cells {st.get('line_clusters', 0)} · intact {len(cell['intact']) if cell else 0}"
                 f" · severed {len(cell['severed']) if cell else 0}\n"
                 f"edge α median {np.nanmedian(e['alpha_front_10']):.0f}° "
                 f"(n={len(e)}) · area {top['area_m2']:.0f} m²", fontsize=8.5)
    return im


def main(run_dir):
    run = Path(run_dir)
    summ = json.loads((run / 'summary.json').read_text())
    freqs = dict(np.load(run / 'inside_frequency.npz'))
    data = Path(summ['meta'].get('data_dir', DATA))
    with rasterio.open(data / 'dem_1m.tif') as src:
        dem = src.read(1).astype(float)
        dem[dem <= -9999] = np.nan
        tf = src.transform
    from release_areas.generate_scenarios import load_kml_mask
    from release_areas.release_geometry import load_observed_polygon
    from release_areas import scoring
    sz = load_kml_mask(data / 'boundaries' / 'start_zone.kml', dem.shape, tf)
    obs = load_observed_polygon(Path(summ['meta']['observed']))
    rr, cc = np.where(sz | scoring.rasterize([obs], dem.shape, tf)[0])
    pad = 15
    r0, r1 = max(rr.min() - pad, 0), min(rr.max() + pad, dem.shape[0])
    c0, c1 = max(cc.min() - pad, 0), min(cc.max() + pad, dem.shape[1])
    ext = (tf.c + c0 * tf.a, tf.c + c1 * tf.a, tf.f + r1 * tf.e, tf.f + r0 * tf.e)
    hs = hillshade(dem)[r0:r1, c0:c1]
    crop = lambda a: a[r0:r1, c0:c1]

    five = {s['mode']: s for s in summ['five_trigger']}
    modes = list(five)
    ref = modes[0]
    alpha = load_alpha(run, data)
    nrow = 4 if alpha else 3
    fig, axes = plt.subplots(nrow, len(modes), figsize=(3.3 * len(modes), 3.75 * nrow),
                             squeeze=False, constrained_layout=True)
    im = {}
    for j, m in enumerate(modes):
        s = five[m]
        polys = [wkt.loads(p['wkt']) for p in s['per_trigger']]
        prob = np.zeros(dem.shape)
        for msk, p in zip(scoring.rasterize(polys, dem.shape, tf), s['per_trigger']):
            prob[msk] += p['weight']
        f = freqs.get(m.replace(':', '_'))
        fr = freqs.get(ref.replace(':', '_'))
        rows = [(np.where(prob > 0, prob, np.nan), SEQ, 0, 1, 'p'),
                (np.where(f > 0, f, np.nan) if f is not None else None, SEQ, 0, 1, 'f'),
                ((f - fr) if (f is not None and fr is not None and j > 0) else None,
                 DIV, -1, 1, 'd')]
        for i, (arr, cmap, vmin, vmax, key) in enumerate(rows):
            ax = axes[i][j]
            ax.imshow(hs, extent=ext, cmap='gray', vmin=-0.2, vmax=1.1, interpolation='bilinear')
            if arr is not None:
                a = crop(arr)
                if key == 'd':
                    a = np.where(np.abs(a) > 0.02, a, np.nan)
                im[key] = ax.imshow(a, extent=ext, cmap=cmap, vmin=vmin, vmax=vmax,
                                    alpha=0.85, interpolation='nearest')
            elif key == 'd':
                ax.text(0.5, 0.5, 'reference', transform=ax.transAxes, ha='center',
                        va='center', fontsize=10, color='#333333')
            ax.contour(crop(sz).astype(float), levels=[0.5], extent=(ext[0], ext[1], ext[3], ext[2]),
                       colors='#555555', linewidths=0.8, linestyles='--', origin='upper')
            ax.plot(*obs.exterior.xy, color='black', lw=1.4)
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
        sd = next((x for x in summ['forced_shots'] if x['mode'] == m and x['class'] == 'inside'), {})
        axes[0][j].set_title(f"{m}\nE[IoU] {s['expected_iou']:.3f} · E[area] "
                             f"{s['expected_area_m2']:.0f} m²", fontsize=9.5)
        if alpha:
            im['a'] = alpha_panel(axes[3][j], m, five[m], alpha, hs, ext, crop, sz, obs, tf,
                                  dem.shape)
        axes[1][j].set_title(f"inside median IoU {sd.get('iou_median', float('nan')):.3f}\n"
                             f"median area {sd.get('area_median_m2', float('nan')):.0f} m²",
                             fontsize=9)
    axes[0][0].set_ylabel('five-trigger P(release)', fontsize=10)
    axes[1][0].set_ylabel('inside shots: coverage', fontsize=10)
    axes[2][0].set_ylabel(f'coverage minus {ref}', fontsize=10)
    if 'p' in im:
        fig.colorbar(im['p'], ax=axes[0].tolist(), shrink=0.8, label='P(pixel in release | release)')
    if 'f' in im:
        fig.colorbar(im['f'], ax=axes[1].tolist(), shrink=0.8, label='fraction of released shots')
    if 'd' in im:
        fig.colorbar(im['d'], ax=axes[2].tolist(), shrink=0.8, label=f'difference vs {ref}')
    if 'a' in im:
        axes[3][0].set_ylabel('flank lines, top trigger', fontsize=10)
        fig.colorbar(im['a'], ax=axes[3].tolist(), shrink=0.8,
                     label='α: front tangent vs fall line (deg)')
    fig.suptitle(f"{summ['meta'].get('date', '2026-01-18')}, calibration {summ['meta']['calibration'].get('version')} — "
                 f"black: survey release ({summ['meta']['observed_m2']:.0f} m²), "
                 "dashed: start zone. Calibration data, not validation.", fontsize=11)
    out = run / 'map.png'
    fig.savefig(out, dpi=140)
    print(out)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else REPO / 'outputs' / 'score_jan18')
