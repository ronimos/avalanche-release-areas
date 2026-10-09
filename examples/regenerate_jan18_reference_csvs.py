"""Regenerate the Jan 18 2026 reference feature CSVs with the current defaults.

Writes the canonical pair:
    all_start_zone_features_2026-01-18.csv
    meloche_features_all_2026-01-18.csv

Current defaults at the time of writing (2026-10-09): thickness-weighted
element means, summed layer thicknesses, `vanherwijnen2016` E, the layered-slab
aggregation with `schweizer2004` at the printed exponent 1.9, the ligament
bound, the mode III cross-slope columns, and `THETA_ESTIMATOR = 'plane_fit'`.

The cluster set and the `group` labels are taken from the pre-2026-10-09
baseline (now *_v1.csv) so the only difference is how parameters are generated.

Needs xarray/zarr, which are not in this repo's venv. Run with avachain's
interpreter, which also resolves `release_areas` to this working tree:

    /home/ron/avachain/.venv/bin/python -I \
        examples/regenerate_jan18_reference_csvs.py
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import xarray as xr

from release_areas import config
from release_areas import arrest_indices as ai
from release_areas.snowpack_features import (profile_features,
                                             compute_meloche_features)

REPO = Path(__file__).resolve().parents[1]
FEAT = REPO / 'data/little_prof/features'
ZARR = '/home/ron/avachain/snowpack/little_prof/output/slope_snowpack.zarr'
NLOC = 6565                      # location axis is this block repeated 9x
SNAP = pd.Timestamp('2026-01-18 12:00')
MIN_DEPTH_CM = 30.0

print(f'THETA_ESTIMATOR = {config.THETA_ESTIMATOR!r}  '
      f'(radius {config.THETA_FIT_RADIUS_M} m)')
print(f'SCH2004_EXP     = {ai.SCH2004_EXP}')
print(f'E relation      = {ai.E_DEFAULT if hasattr(ai, "E_DEFAULT") else "see config"}')

base = pd.read_csv(FEAT / 'all_start_zone_features_2026-01-18_v1.csv')
groups = dict(zip(base.cluster_id, base.group))
print(f'baseline: {len(base)} clusters, groups {base.group.value_counts().to_dict()}')

ds = xr.open_zarr(ZARR).isel(location=slice(0, NLOC))
loc_names = np.array([str(x) for x in ds.coords['location'].values])
loc_ids = np.array([int(x.split('_')[-1]) for x in loc_names])

keep = np.isin(loc_ids, base.cluster_id.to_numpy())
print(f'matched {keep.sum()} of {len(base)} baseline clusters on the location axis')

ds_sub = ds.isel(location=keep)
ids_sub = loc_ids[keep]
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    ds_t = ds_sub.sel(time=SNAP, method='nearest').compute()
print(f'snapshot {pd.Timestamp(ds_t.coords["time"].values)}')

rows = []
for i in range(ds_t.sizes['location']):
    feat = profile_features(ds_t.isel(location=i), MIN_DEPTH_CM)
    if feat:
        feat['cluster_id'] = int(ids_sub[i])
        feat['group'] = groups.get(int(ids_sub[i]), 'start_zone')
        rows.append(feat)

feats = pd.DataFrame(rows).set_index('cluster_id')
feats = feats[~feats.index.duplicated(keep='first')]
out_all = FEAT / 'all_start_zone_features_2026-01-18.csv'
feats.to_csv(out_all)
print(f'wrote {out_all.name}: {len(feats)} rows, {len(feats.columns)} columns')

# Post-fix identity check: the basal WL starts at the ground.
ok = feats.dropna(subset=['hs', 'slab_thickness', 'wl_thickness'])
hit = np.isclose(ok['wl_thickness'], ok['hs'] - ok['slab_thickness'], rtol=1e-9)
print(f'  D_wl == hs - slab_thickness: {hit.sum()}/{len(ok)}')

with rasterio.open(REPO / 'data/little_prof/spatial/cluster_map.tif') as src:
    cluster_map = src.read(1)
    transform = src.transform
with rasterio.open(REPO / 'data/little_prof/dem_1m.tif') as src:
    dem = src.read(1)

snap_data = {g: d.drop(columns=['group'])
             for g, d in feats.groupby('group')}
mel = compute_meloche_features(snap_data, cluster_map, dem, transform,
                               snap_ts=SNAP)
out_mel = FEAT / 'meloche_features_all_2026-01-18.csv'
mel.to_csv(out_mel)
print(f'wrote {out_mel.name}: {len(mel)} rows, {len(mel.columns)} columns')

th = pd.to_numeric(mel['theta'], errors='coerce')
aca = pd.to_numeric(mel['A_ca_brittle'], errors='coerce')
print(f'  theta        median {th.median():8.3f} Pa/m   NaN {th.isna().mean()*100:.1f}%')
print(f'  A_ca_brittle median {aca.median():8.1f} m      NaN {aca.isna().mean()*100:.1f}%')
