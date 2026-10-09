"""Statistics for the full space-time ligament-bound pass.

Reads the per-block .npz files written by spacetime_ligament_pass.py and
prints the tables recorded in methods §4, "Measured parameter statistics".
Set LIG_OUT_DIR / LIG_ZARR to match the pass run.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

OUT_DIR = Path(os.environ.get('LIG_OUT_DIR', '/tmp/lig_bench/out'))
ZARR = os.environ.get(
    'LIG_ZARR',
    '/home/ron/avachain/snowpack/little_prof/output/slope_snowpack.zarr')

# Derived from config, never hardcoded: config.py is the single source of truth
# and a stale copy here would report on columns that no longer exist.
from release_areas import config as _cfg
from release_areas import arrest_indices as _ai

_NOTE = {'kirchner2000': '  (LOWER BOUND on K_Ic)'}
OPTIONS = [
    (f'{relation}, d_max = {dmax:.1f} gsz'
     + ('  [default]' if (relation == _ai.K_IC_DEFAULT
                          and dmax == _cfg.DMAX_FACTOR) else '')
     + _NOTE.get(relation, ''),
     f'__{tag}')
    for tag, relation, dmax in _cfg.LIGAMENT_VARIANTS
]


def load():
    files = sorted(OUT_DIR.glob('block_*.npz'))
    print(f'{len(files)} block files')
    cols = {}
    for f in files:
        with np.load(f) as z:
            for k in z.files:
                cols.setdefault(k, []).append(z[k])
    return pd.DataFrame({k: np.concatenate(v) for k, v in cols.items()})


def q(x, p):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.percentile(x, p)) if len(x) else np.nan


def med(x):
    return q(x, 50)


def report(df, label):
    n = len(df)
    slab = df['n_slab_layers'].to_numpy() > 0
    n_slab = int(slab.sum())
    print()
    print('=' * 74)
    print(f'{label}')
    print('=' * 74)
    print(f'profiles walked                  {n:>12,}')
    print(f'  with a resolved slab           {n_slab:>12,}  '
          f'({100*n_slab/max(n,1):.1f}%)')
    d = df[slab]
    print(f'  slab thickness h (m)           '
          f'median {med(d["slab_thickness"]):.3f}  '
          f'p5 {q(d["slab_thickness"],5):.3f}  p95 {q(d["slab_thickness"],95):.3f}')
    print(f'  slab density rho (kg/m3)       '
          f'median {med(d["slab_density"]):.1f}  '
          f'p5 {q(d["slab_density"],5):.1f}  p95 {q(d["slab_density"],95):.1f}')
    print(f'  slab layers per profile        '
          f'median {med(d["n_slab_layers"]):.0f}  '
          f'p5 {q(d["n_slab_layers"],5):.0f}  p95 {q(d["n_slab_layers"],95):.0f}')

    # --- layered vs bulk sigma_t and E (K_Ic-independent) ---
    print()
    print('-- layered aggregation vs bulk (independent of the K_Ic option) --')
    r = d['sigma_t_mean'] / d['sigma_t']
    print(f'  sigma_t_mean / sigma_t(bulk)   '
          f'median {med(r):.4f}  p5 {q(r,5):.4f}  p95 {q(r,95):.4f}')
    r = d['sigma_t_wl'] / d['sigma_t_mean']
    print(f'  sigma_t_wl / sigma_t_mean      '
          f'median {med(r):.4f}  p5 {q(r,5):.4f}  p95 {q(r,95):.4f}')
    print(f'  f_mean_weighted (facetedness)  '
          f'median {med(d["f_mean_weighted"]):.4f}  '
          f'p5 {q(d["f_mean_weighted"],5):.4f}  p95 {q(d["f_mean_weighted"],95):.4f}')
    a_h = d['wl_ctrl_thickness'] / d['slab_thickness']
    print(f'  controlling layer a (m)        '
          f'median {med(d["wl_ctrl_thickness"]):.4f}  '
          f'p95 {q(d["wl_ctrl_thickness"],95):.4f}')
    print(f'  controlling a/h                '
          f'median {med(a_h):.4f}  p95 {q(a_h,95):.4f}')

    # --- E relation / A_ca ratio ---
    print()
    print('-- E relation: A_ca(vanherwijnen2016) / A_ca(project_fit) --')
    e_vh, e_pf = d['E_eff'], d['E_eff__project_fit']
    print(f'  E_eff vanherwijnen2016 (Pa)    '
          f'median {med(e_vh):.4g}  p5 {q(e_vh,5):.4g}  p95 {q(e_vh,95):.4g}')
    print(f'  E_eff project_fit (Pa)         '
          f'median {med(e_pf):.4g}  p5 {q(e_pf,5):.4g}  p95 {q(e_pf,95):.4g}')
    re = e_pf / e_vh
    print(f'  E_pf / E_vh                    '
          f'median {med(re):.4f}  p5 {q(re,5):.4f}  p95 {q(re,95):.4f}')
    ra = np.sqrt(re)
    print(f'  A_ca ratio = sqrt(E_pf/E_vh)   '
          f'median {med(ra):.4f}  p5 {q(ra,5):.4f}  p95 {q(ra,95):.4f}')
    print('    (A_ca depends on E only through Lambda ~ sqrt(E); L_t, tau_g,')
    print('     theta and sigma_t are E-free, so this ratio is exact.)')

    # --- per K_Ic option ---
    for name, sfx in OPTIONS:
        K = d['K_Ic_ligament' + sfx]
        sc = d['sigma_t_ligament' + sfx]
        si = d['sigma_t_lig_intact' + sfx]
        a0 = d['a0_ligament' + sfx]
        lch = d['l_ch_ligament' + sfx]
        ncl = d['n_rho_clamped' + sfx]

        have_int = np.isfinite(si)
        have = np.isfinite(sc)
        withheld = have_int & ~have
        print()
        print(f'-- K_Ic option: {name} --')
        print(f'  bound returned                 {int(have.sum()):>12,}  '
              f'({100*have.sum()/max(n_slab,1):.1f}% of slabs)')
        print(f'  withheld (a0 >= h, LEFM void)  {int(withheld.sum()):>12,}  '
              f'({100*withheld.sum()/max(have_int.sum(),1):.1f}% of intact-ligament cases)')
        print(f'  K_Ic (Pa m^0.5)                '
              f'median {med(K):.1f}  p5 {q(K,5):.1f}  p95 {q(K,95):.1f}')
        print(f'  a0 intrinsic flaw (m)          '
              f'median {med(a0):.4f}  p5 {q(a0,5):.4f}  p95 {q(a0,95):.4f}')
        print(f'  l_ch = (K_Ic/sigma_t)^2 (m)    '
              f'median {med(lch):.4f}  p95 {q(lch,95):.4f}')
        print(f'  sigma_t_lig_intact (Pa)        '
              f'median {med(si):.1f}  p5 {q(si,5):.1f}  p95 {q(si,95):.1f}')
        print(f'  sigma_c ligament bound (Pa)    '
              f'median {med(sc):.1f}  p5 {q(sc,5):.1f}  p95 {q(sc,95):.1f}')
        r = sc / si
        print(f'  sigma_c / sigma_t_lig_intact   '
              f'median {med(r):.4f}  p5 {q(r,5):.4f}  p95 {q(r,95):.4f}')
        r = sc / d['sigma_t']
        print(f'  sigma_c / sigma_t(bulk)        '
              f'median {med(r):.4f}  p5 {q(r,5):.4f}  p95 {q(r,95):.4f}')
        print(f'  clamped layers per profile     '
              f'median {med(ncl):.0f}  p5 {q(ncl,5):.0f}  p95 {q(ncl,95):.0f}')
        fr = ncl / d['n_slab_layers']
        print(f'  clamped fraction of slab       '
              f'median {med(fr):.4f}  p5 {q(fr,5):.4f}  p95 {q(fr,95):.4f}')
        print(f'  profiles with >=1 clamped      '
              f'{int((ncl.to_numpy() > 0).sum()):>12,}  '
              f'({100*(ncl.to_numpy()>0).sum()/max(n_slab,1):.1f}% of slabs)')


def main():
    df = load()
    times = pd.DatetimeIndex(xr.open_zarr(ZARR).coords['time'].values)
    print(f'rows {len(df):,} = {df.cluster_id.nunique()} clusters x '
          f'{df.ti.nunique()} timesteps')
    print(f'time span {times[0]} -> {times[-1]}')

    report(df, 'FULL SPACE-TIME: all clusters, all 501 timesteps')

    jan18 = int(np.argmin(np.abs(times - pd.Timestamp('2026-01-18 12:00'))))
    report(df[df.ti == jan18],
           f'JAN 18 2026 12:00 SNAPSHOT ONLY (ti={jan18}, {times[jan18]})')


if __name__ == '__main__':
    main()
