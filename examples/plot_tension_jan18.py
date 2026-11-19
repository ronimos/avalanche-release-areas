#!/usr/bin/env python3
# Map of examples/release_scale_tension_jan18.py output (tension_map.png).
# env -u VIRTUAL_ENV uv run --frozen python examples/plot_tension_jan18.py
import json, sys
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, rasterio
from pathlib import Path
REPO = Path('/home/ron/avalanche-release-areas'); DATA = REPO/'data/little_prof'
run = REPO/'outputs/score_jan18/tension_2026-11-19'
cells = json.loads((run/'tension_cells.json').read_text())
from release_areas.release_geometry import cluster_pixel_index, load_observed_polygon
from release_areas.generate_scenarios import load_kml_mask
with rasterio.open(DATA/'spatial/cluster_map.tif') as s: cmap = s.read(1).astype(float); tf = s.transform
with rasterio.open(DATA/'dem_1m.tif') as s: dem = s.read(1).astype(float); dem[dem<=-9999]=np.nan
px = cluster_pixel_index(cmap)
obs = load_observed_polygon(DATA/'validation/survey_release_area_20260118.geojson')
sz = load_kml_mask(DATA/'boundaries/start_zone.kml', dem.shape, tf)
gy, gx = np.gradient(np.nan_to_num(dem, nan=np.nanmean(dem)))
sl = np.pi/2-np.arctan(np.hypot(gx,gy)); asp = np.arctan2(-gx,gy)
hs = np.sin(np.radians(45))*np.sin(sl)+np.cos(np.radians(45))*np.cos(sl)*np.cos(np.radians(315)-asp)
rr, cc = np.where(sz); pad=15
r0,r1,c0,c1 = rr.min()-pad, rr.max()+pad, cc.min()-pad, cc.max()+pad
ext = (tf.c+c0*tf.a, tf.c+c1*tf.a, tf.f+r1*tf.e, tf.f+r0*tf.e)
cen = lambda c: (tf.c+(px[c][1].mean()+.5)*tf.a, tf.f+(px[c][0].mean()+.5)*tf.e)
cells = sorted(cells, key=lambda x: x['trigger'])
fig, axes = plt.subplots(1, len(cells)+1, figsize=(3.4*(len(cells)+1), 4.3), constrained_layout=True)
for ax, x in zip(axes, cells):
    a = np.full(dem.shape, np.nan)
    for c, r in zip(x['ids'], x['ratio']):
        if r is not None and c in px: a[px[c]] = r
    ax.imshow(hs[r0:r1,c0:c1], extent=ext, cmap='gray', vmin=-.2, vmax=1.1)
    im = ax.imshow(np.log10(a[r0:r1,c0:c1]), extent=ext, cmap='RdBu_r', vmin=-1, vmax=1, interpolation='nearest')
    ax.plot(*obs.exterior.xy, 'k', lw=1.1)
    gc = [cen(c) for c,_,_ in x['cracks'] if c in px]
    ax.plot([p[0] for p in gc], [p[1] for p in gc], '.', color='#1a9641', ms=4, label='gate cracks')
    t = cen(x['trigger']); ax.plot(*t, '*', color='gold', mec='k', ms=12, label='trigger')
    if x['first'] in px: ax.plot(*cen(x['first']), 'X', color='magenta', mec='k', ms=10, label='first σ ≥ σt')
    fr = next(r['ring'] for r in x['per_ring'] if r['exceed'] > 0)
    gr = min(r for _,_,r in x['cracks'])
    ax.set_title(f"trigger {x['trigger']}\nfirst exceedance ring {fr}, first gate crack ring {gr}", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlim(ext[0],ext[1]); ax.set_ylim(ext[2],ext[3])
axes[0].legend(fontsize=7, loc='lower left')
fig.colorbar(im, ax=axes[:-1].tolist(), shrink=.8, label='log10 σ/σt at the end of the BFS')
ax = axes[-1]
for x in cells:
    pr = x['per_ring']
    ax.plot([r['ring'] for r in pr], [r['exceed']/r['admitted'] for r in pr], label=str(x['trigger']))
ax.set_xlabel('BFS ring'); ax.set_ylabel('admitted clusters with σ ≥ σt (fraction)')
ax.set_ylim(0, 1); ax.legend(fontsize=7); ax.grid(alpha=.3)
fig.suptitle('Jan 18, cal 2026-11-17.1: release-scale crown tension σ = g·drive·ρh·L_below/h (strip, flanks unloaded) vs σt', fontsize=10)
fig.savefig(run/'tension_map.png', dpi=130); print(run/'tension_map.png')
