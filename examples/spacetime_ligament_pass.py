"""Full space-time pass of profile_features over the Little Professor zarr.

Runs with config.USE_LIGAMENT_BOUND = True, so aggregate_slab computes the
ligament bound for the default relation plus every entry in
config.LIGAMENT_VARIANTS in one pass.

The zarr location axis is 59280 = 6565 x 9 exact duplicates (stride 6565,
verified byte-identical for density/grain_type/HS at several clusters and
times), so only the first 6565 are walked.

Writes one parquet per location block to OUT_DIR.
"""
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ZARR = os.environ.get(
    'LIG_ZARR',
    '/home/ron/avachain/snowpack/little_prof/output/slope_snowpack.zarr')
NLOC = 6565
MIN_DEPTH_CM = 30.0
BLOCK = 50
# Override with LIG_OUT_DIR. The default is under /tmp, so the blocks are
# ephemeral — point it somewhere persistent before a run you want to keep.
OUT_DIR = Path(os.environ.get('LIG_OUT_DIR', '/tmp/lig_bench/out'))

VARS = ['z', 'grain_type', 'HS', 'density', 'hand_hardness', 'grain_size',
        'shear_strength', 'sk38', 'ssi', 'sn38', 'stab_deformation_rate',
        'sphericity', 'dendricity']

VARIANT_KEYS = ('sigma_t_ligament', 'K_Ic_ligament', 'sigma_t_lig_intact',
                'a0_ligament', 'l_ch_ligament', 'n_rho_clamped')

# Derived from config, never hardcoded: config.py is the single source of truth,
# and a stale copy here would silently collect columns aggregate_slab no longer
# emits. '' is the default relation, which is emitted unsuffixed.
from release_areas import config as _cfg               # noqa: E402
SUFFIXES = ('',) + tuple(f'__{tag}' for tag, _, _ in _cfg.LIGAMENT_VARIANTS)

KEEP = ['slab_density', 'slab_thickness', 'slab_grain_size', 'wl_thickness',
        'hs', 'sigma_t', 'E_slab', 'Lambda', 'n_slab_layers',
        'E_eff', 'E_eff__project_fit', 'sigma_t_mean', 'sigma_t_wl',
        'f_mean_weighted', 'wl_ctrl_thickness', 'wl_ctrl_depth',
        'wl_ctrl_index'] + [k + s for s in SUFFIXES for k in VARIANT_KEYS]


def run_block(args):
    lo, hi = args
    import xarray as xr
    from release_areas.snowpack_features import profile_features

    out = OUT_DIR / f'block_{lo:05d}.npz'
    if out.exists():
        return str(out), -1

    ds = xr.open_zarr(ZARR)[VARS].isel(location=slice(lo, hi)).compute()
    loc_names = [str(x) for x in ds.coords['location'].values]
    cl_ids = np.array([int(x.split('_')[-1]) for x in loc_names])
    n_t = ds.sizes['time']

    rows = []
    for i in range(hi - lo):
        ds_i = ds.isel(location=i)
        for ti in range(n_t):
            feat = profile_features(ds_i.isel(time=ti), MIN_DEPTH_CM)
            rec = {'cluster_id': cl_ids[i], 'ti': ti}
            for k in KEEP:
                v = feat.get(k, np.nan)
                rec[k] = np.nan if v is None else v
            rows.append(rec)

    cols = {'cluster_id': np.array([r['cluster_id'] for r in rows], dtype='int32'),
            'ti': np.array([r['ti'] for r in rows], dtype='int16')}
    for k in KEEP:
        cols[k] = np.array([_f(r[k]) for r in rows], dtype='float32')
    np.savez_compressed(out, **cols)
    return str(out), len(rows)


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def main():
    import multiprocessing as mp
    import release_areas.config as cfg

    assert cfg.USE_LIGAMENT_BOUND, 'USE_LIGAMENT_BOUND must be True'
    print('config: E_RELATION           =', cfg.E_RELATION)
    print('        E_RELATIONS_EXTRA    =', cfg.E_RELATIONS_EXTRA)
    print('        K_IC_RELATION        =', cfg.K_IC_RELATION)
    print('        DMAX_FACTOR          =', cfg.DMAX_FACTOR)
    print('        LIGAMENT_MODEL       =', cfg.LIGAMENT_MODEL)
    print('        SELF_CONSISTENT_A0   =', cfg.SELF_CONSISTENT_A0)
    print('        LIGAMENT_VARIANTS    =', cfg.LIGAMENT_VARIANTS)
    print('        FACET_SP_REF         =', cfg.FACET_SP_REF)
    print(flush=True)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    blocks = [(lo, min(lo + BLOCK, NLOC)) for lo in range(0, NLOC, BLOCK)]
    n_work = int(sys.argv[1]) if len(sys.argv) > 1 else 24

    t0 = time.time()
    done = 0
    with mp.Pool(n_work, maxtasksperchild=4) as pool:
        for path, n in pool.imap_unordered(run_block, blocks):
            done += 1
            el = time.time() - t0
            print(f'[{done}/{len(blocks)}] {Path(path).name} rows={n} '
                  f'{el:.0f}s elapsed, eta {el/done*(len(blocks)-done):.0f}s',
                  flush=True)
    print(f'\nall blocks done in {time.time()-t0:.0f}s')


if __name__ == '__main__':
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    main()
