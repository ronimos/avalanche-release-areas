#!/usr/bin/env python3
# Flank load if the slab breaks behind the front (L_block = min(L, L*)) and the
# width-based force balance (W/2); writes outputs/score_jan18/tension_2026-11-19/flank_blocks.csv.
# env -u VIRTUAL_ENV uv run --frozen python examples/flank_load_blocks_jan18.py
"""Flank load with the slab already broken behind the front: L_block = min(L, L*)."""
import contextlib, io, json, sys, multiprocessing as mp
from pathlib import Path
import numpy as np, pandas as pd
REPO = Path('/home/ron/avalanche-release-areas')
sys.path.insert(0, str(REPO / 'examples'))
import score_jan18 as sj
from release_areas import config, scoring
from release_areas import release_geometry as rg
from release_areas.generate_scenarios import select_triggers
ORIG = rg.slab_crack_lines
CAP = []
def wrap(*a, **k):
    tr = []; k['trace_out'] = tr
    out = ORIG(*a, **k); CAP.append(tr); return out
def shot(job):
    cid, lab = job
    CAP.clear()
    with contextlib.redirect_stdout(io.StringIO()):
        sj._polygon(cid)
    return [dict(t, set=lab, trigger=cid) for tr in CAP for t in tr if 'k1' in t]
if __name__ == '__main__':
    cal = json.loads((sj.DATA / 'calibration_2026-11-17.1.json').read_text())
    sj.apply_calibration(cal); sj.set_mode('mixed:energy:intact:front_edge')
    sj.load_inputs(sj.DATA / cal['features_dir'], sj.DATA / 'validation/survey_release_area_20260118.geojson', None)
    rg.slab_crack_lines = wrap
    with contextlib.redirect_stdout(io.StringIO()):
        trig = select_triggers(sj.G['feats'], sj.G['mel'], sj.G['cluster_map'], sj.G['dem'], sj.G['start_zone_mask'], 5)
    ins = scoring.inside_clusters(sj.G['cluster_map'], sj.G['obs_mask'], sj.G['start_zone_mask'], sj.G['observed'], sj.G['transform'])
    jobs = [(int(c), 'five') for c in trig] + [(int(r['cid']), 'inside') for r in ins]
    with mp.get_context('fork').Pool(24) as p:
        rows = [r for res in p.imap_unordered(shot, jobs, chunksize=2) for r in res]
    d = pd.DataFrame(rows)
    f = sj.G['feats']
    st = f[config.SLAB_TENSILE_SIGMA_COL].reindex(d.cid).values
    d['Lstar'] = st / (9.81 * d.drive * d.rho)
    d['Lblock'] = np.minimum(d.L, d.Lstar)
    s = d.Lblock / d.L
    d['brk_now'] = d.k1**2 + d.k2**2 >= d.kic**2
    d['brk_block'] = (d.k1 * s)**2 + (d.k2 * s)**2 >= d.kic**2
    q = lambda x: ' / '.join(f'{v:.2f}' for v in np.nanpercentile(x, [10, 50, 90]))
    for lab, g in d.groupby('set'):
        print(f'{lab}: {len(g)} mixed flank steps (intact, front_edge)')
        print('  L (m)            p10/50/90', q(g.L))
        print('  L* (m)           p10/50/90', q(g.Lstar))
        print('  L/L* > 1 share', round((g.L > g.Lstar).mean(), 2), '; load ratio Lblock/L p10/50/90', q(g.Lblock / g.L))
        print('  K/K_Ic now       p10/50/90', q(np.hypot(g.k1, g.k2) / g.kic))
        print('  K/K_Ic block     p10/50/90', q(np.hypot(g.k1, g.k2) * g.Lblock / g.L / g.kic))
        print(f'  steps that break: now {g.brk_now.mean():.2f}, with blocks {g.brk_block.mean():.2f}')
        for ax, h in g.groupby('axes'):
            print(f'    {ax:5s} n={len(h)} break now {h.brk_now.mean():.2f} block {h.brk_block.mean():.2f}')
    d.to_csv(REPO / 'outputs/score_jan18/tension_2026-11-19/flank_blocks.csv', index=False)
    # force balance: flank shear from the detached slab beside the line,
    # shared by two flanks: sigma = rho g (W/2) drive |d|, independent of block length
    for lab, g in d.groupby('set'):
        sW = (g.W / 2) / g.L
        print(f'{lab}: W (m) p10/50/90', q(g.W), '| W/2 over L p10/50/90', q(sW))
        kw = np.hypot(g.k1, g.k2) * sW / g.kic
        print('  K/K_Ic width-based p10/50/90', q(kw), ' break', round((kw >= 1).mean(), 2))
        for ax, h in g.groupby('axes'):
            k = np.hypot(h.k1, h.k2) * (h.W / 2) / h.L / h.kic
            print(f'    {ax:5s} break width-based {(k >= 1).mean():.2f}')
