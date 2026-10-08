"""
snowpack_features.py — Per-cluster SNOWPACK feature extraction and Meloche crack-arrest parameters.

Provides:
    geojson_to_mask()         Rasterize a GeoJSON polygon to a boolean mask
    assign_cluster_groups()   Assign clusters to release / adjacent / reference groups
    split_wl_slab()           Identify WL / slab boundary via grain type
    element_thickness()       Per-element thickness dz from the z coordinate
    profile_features()        Per-cluster WL/slab/stability features (from xarray Dataset)
    compute_meloche_features() Meloche et al. (2025) crack-arrest parameters

profile_features() expects an xarray Dataset slice (one location × one timestep) as
produced by the SNOWPACK zarr pipeline.  If you have pre-computed feature CSVs, use
those directly with compute_meloche_features().
"""

from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd
import rasterio
import rasterio.features
import rasterio.transform as rt

from release_areas import config
from release_areas import arrest_indices as _ai
from release_areas.arrest_indices import evaluate as _ai_evaluate

# Slab/weak-layer material constants (config is the source of truth)
G_WL = config.G_WL   # weak-layer shear modulus (Pa)
NU   = config.NU     # slab Poisson's ratio

# SNOWPACK zarr variable holding critical cut length r_c (m)
RC_VAR = 'critical_cut_length'

GROUP_LABELS = {
    'release':   'Release zone',
    'adjacent':  'Adjacent slope (start zone)',
    'reference': 'Reference (terrain-matched)',
}


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def geojson_to_mask(path, dem_shape, transform, dst_epsg: int = config.RASTER_EPSG) -> np.ndarray:
    """Rasterize a GeoJSON polygon to a boolean mask matching dem_shape."""
    from shapely.geometry import shape
    from shapely.ops import unary_union
    from pyproj import Transformer as _T

    with open(str(path)) as f:
        gj = json.load(f)
    polys = [shape(feat['geometry']) for feat in gj['features']]
    merged = unary_union(polys)

    # Detect source CRS; handles both 'EPSG:32613' and 'urn:ogc:def:crs:EPSG::32613'
    src_epsg = 4326
    crs_node = gj.get('crs', {}).get('properties', {}).get('name', '')
    _m = re.search(r'EPSG:+(\d+)', crs_node, re.IGNORECASE)
    if _m:
        src_epsg = int(_m.group(1))

    if src_epsg != dst_epsg:
        tr = _T.from_crs(f'EPSG:{src_epsg}', f'EPSG:{dst_epsg}', always_xy=True)
        from shapely.geometry import Polygon, MultiPolygon

        def _reproj_ring(r):
            xs, ys = zip(*r.coords)
            xs2, ys2 = tr.transform(xs, ys)
            return list(zip(xs2, ys2))

        def reproj(poly):
            if poly.geom_type == 'Polygon':
                return Polygon(_reproj_ring(poly.exterior),
                               [_reproj_ring(i) for i in poly.interiors])
            return MultiPolygon([reproj(p) for p in poly.geoms])

        merged = reproj(merged)

    aff = transform if hasattr(transform, 'c') else rt.Affine(*transform[:6])
    mask = rasterio.features.geometry_mask(
        [merged.__geo_interface__], out_shape=dem_shape,
        transform=aff, invert=True)
    print(f"  Release zone mask:  {mask.sum()} cells")
    return mask


# ---------------------------------------------------------------------------
# Cluster group assignment
# ---------------------------------------------------------------------------

def assign_cluster_groups(cluster_map, release_mask, start_zone_mask,
                           dem, domain_mask) -> dict:
    """
    Assign clusters to release / adjacent / reference groups.

    Reference group is terrain-matched to release zone (similar elevation
    and slope angle).
    """
    fill_dem = np.where(np.isnan(dem), np.nanmean(dem), dem)
    dy, dx   = np.gradient(fill_dem, 1.0)
    slope    = np.degrees(np.arctan(np.sqrt(dx**2 + dy**2)))

    groups = {'release': set(), 'adjacent': set(), 'reference': set()}

    no_release = not release_mask.any()
    if no_release:
        rel_elev_mean = rel_elev_std = rel_slope_mean = None
    else:
        rel_elev_mean  = float(np.nanmean(dem[release_mask & domain_mask]))
        rel_elev_std   = float(np.nanstd(dem[release_mask & domain_mask]))
        rel_slope_mean = float(np.nanmean(slope[release_mask & domain_mask]))

    for cid in np.unique(cluster_map[domain_mask]):
        if cid <= 0:
            continue
        cells = cluster_map == cid
        n     = cells.sum()
        if n == 0:
            continue

        frac_rel   = 0.0 if no_release else (cells & release_mask).sum() / n
        frac_start = (cells & start_zone_mask).sum() / n

        if frac_rel >= 0.3:
            groups['release'].add(cid)
        elif frac_start >= 0.3:
            groups['adjacent'].add(cid)
        elif not no_release:
            c_elev  = float(np.nanmean(dem[cells]))
            c_slope = float(np.nanmean(slope[cells]))
            if (abs(c_elev - rel_elev_mean) < rel_elev_std * 1.5 and
                    abs(c_slope - rel_slope_mean) < 10.0):
                groups['reference'].add(cid)

    for g, ids in groups.items():
        print(f"  {GROUP_LABELS[g]}: {len(ids)} clusters")
    return groups


# ---------------------------------------------------------------------------
# WL / slab boundary detection
# ---------------------------------------------------------------------------

def split_wl_slab(grain_type, z):
    """
    Find the basal weak layer by scanning upward from the bottom.

    Returns (slab_mask, wl_mask, interface_z) or (None, None, None).
    """
    ok = ~np.isnan(grain_type) & ~np.isnan(z) & (z != 0) & (grain_type > 0)
    if ok.sum() < 3:
        return None, None, None

    gt_ok = np.round(grain_type[ok]).astype(int)
    z_ok  = z[ok]

    order        = np.argsort(z_ok)
    gt_sorted    = gt_ok[order]
    z_sorted     = z_ok[order]
    is_wl_sorted = (gt_sorted // 100 == 4) | (gt_sorted // 100 == 5)

    if not is_wl_sorted[0]:
        return None, None, None

    wl_top_sorted = 0
    for i in range(len(is_wl_sorted)):
        if is_wl_sorted[i]:
            wl_top_sorted = i
        else:
            break

    interface_z = float(z_sorted[wl_top_sorted])

    ok_indices = np.where(ok)[0]
    orig_order = ok_indices[order]

    slab_mask = np.zeros(len(grain_type), dtype=bool)
    wl_mask   = np.zeros(len(grain_type), dtype=bool)
    for sorted_i, orig_i in enumerate(orig_order):
        if sorted_i <= wl_top_sorted:
            wl_mask[orig_i] = True
        else:
            slab_mask[orig_i] = True

    return slab_mask, wl_mask, interface_z


def element_thickness(z, hs, valid):
    """
    Per-element thickness dz (m), in the original element index order.

    SNOWPACK's `z` is depth below the snow surface, negative downward, and is
    exactly `height_above_ground - HS` (verified to machine precision against
    the zarr `height` variable). So sorting by z and differencing gives each
    element's own thickness, with the basal element measured from the ground at
    z = -HS. `np.sum(dz) == HS` exactly when no element is excluded.

    Bulk slab properties must be weighted by dz: SNOWPACK element thicknesses
    span a factor of ~10 within a single profile (3-38 mm measured on the
    Little Professor domain), so an unweighted element mean over-weights thin
    layers. Likewise a layer's thickness is the SUM of its elements' dz, not
    the span of their z, which omits the lowest element entirely.

    valid : boolean mask of elements to include; dz is NaN elsewhere.
    """
    z = np.asarray(z, dtype=float)
    dz = np.full(z.shape, np.nan)
    idx = np.where(valid)[0]
    if len(idx) == 0 or not np.isfinite(hs):
        return dz
    order = idx[np.argsort(z[idx])]
    dz[order] = np.diff(z[order], prepend=-float(hs))
    return dz


def thickness_weighted_mean(values, mask, dz):
    """Thickness-weighted mean of `values` over `mask`; NaN if no usable element."""
    values = np.asarray(values, dtype=float)
    m = mask & np.isfinite(values) & np.isfinite(dz) & (dz > 0)
    if not m.any():
        return np.nan
    return float(np.sum(values[m] * dz[m]) / np.sum(dz[m]))


# ---------------------------------------------------------------------------
# Per-cluster feature extraction (requires xarray Dataset slice from zarr)
# ---------------------------------------------------------------------------

def profile_features(ds_t_loc, min_depth_cm: float) -> dict:
    """
    Extract slab and WL features for one cluster at one timestep.

    ds_t_loc : xarray Dataset for one location × one timestep (squeezed).
               Must contain variables: z, grain_type, HS, density,
               hand_hardness, grain_size, shear_strength, sk38, ssi, sn38.
    min_depth_cm : minimum weak-layer burial depth (cm). Profiles whose basal
               weak layer sits shallower than this carry too little slab to
               release and are returned with slab/WL features omitted.

    Returns a dict of scalar features, or {} on failure.
    """
    z      = ds_t_loc['z'].values.ravel()
    gt     = ds_t_loc['grain_type'].values.ravel()
    hs_arr = ds_t_loc['HS'].values.ravel()
    hs     = float(np.nanmean(hs_arr))

    ok = ~np.isnan(z) & ~np.isnan(gt) & (z != 0) & (gt > 0)
    if ok.sum() < 2:
        return {}

    result = {'hs': hs}
    nu = NU

    slab_m, wl_m, interface_z = split_wl_slab(gt, z)
    if slab_m is None:
        return result

    # z is negative downward from the snow surface, so the WL top sits
    # -interface_z below the surface.
    if -interface_z < min_depth_cm / 100.0:
        return result

    # Element thicknesses, used to weight every bulk mean and to sum layer
    # thicknesses. See element_thickness() for why a plain element mean and a
    # z-span are both wrong.
    dz = element_thickness(z, hs, ok)

    near_interface = ok & (np.abs(z - interface_z) <= 0.05)
    for var in ('sk38', 'ssi', 'sn38', 'stab_deformation_rate'):
        try:
            v    = ds_t_loc[var].values.ravel()
            vals = v[near_interface & ~np.isnan(v)]
            result[f'min_{var}'] = float(np.nanmin(vals)) if len(vals) else np.nan
        except Exception:
            result[f'min_{var}'] = np.nan

    if slab_m.any():
        # Thickness-weighted, so that rho * h is the true slab load per unit
        # area and E/sigma_t see a mass-weighted bulk density.
        for var in ('density', 'hand_hardness', 'grain_size'):
            try:
                v = ds_t_loc[var].values.ravel()
                result[f'slab_{var}'] = thickness_weighted_mean(v, slab_m & ok, dz)
            except Exception:
                result[f'slab_{var}'] = np.nan

        slab_z = z[slab_m & ok]
        result['slab_thickness'] = float(-interface_z) if len(slab_z) else np.nan

        gt_slab = np.round(gt[slab_m & ok]).astype(int)
        gt_slab = gt_slab[gt_slab > 0]
        if len(gt_slab):
            classes, counts = np.unique(gt_slab // 100, return_counts=True)
            result['slab_dominant_grain_class'] = int(classes[np.argmax(counts)])
        else:
            result['slab_dominant_grain_class'] = np.nan

        gt_slab_all = np.round(gt[slab_m & ok]).astype(int)
        z_slab_all  = z[slab_m & ok]
        dz_slab_all = dz[slab_m & ok]
        # MF (7xx) and IF (8xx) only. The previous `>= 7` also caught FCxr
        # (9xx, rounding facets), which is 29% of slab thickness here and
        # present in every profile — so every crust flag in the domain was a
        # false positive (zero MF/IF elements exist in it). Faceted layers are
        # handled by the layered sigma_t aggregation instead.
        is_crust    = np.isin(gt_slab_all // 100, (7, 8))
        result['has_crust']       = bool(is_crust.any())
        result['n_crust_layers']  = int(is_crust.sum())
        # Sum of crust element thicknesses. The z span this replaced both
        # omitted the lowest crust element and counted any non-crust elements
        # sandwiched between two crusts.
        result['crust_thickness'] = float(
            np.nansum(dz_slab_all[is_crust])
        ) if is_crust.any() else 0.0
        near_iface = np.abs(z_slab_all - interface_z) <= 0.10
        result['crust_at_interface'] = bool((is_crust & near_iface).any())
        result['crust_top_depth']    = float(
            -z_slab_all[is_crust].max()
        ) if is_crust.any() else np.nan
    else:
        for k in ('slab_density', 'slab_hand_hardness', 'slab_grain_size',
                  'slab_thickness', 'slab_dominant_grain_class',
                  'has_crust', 'n_crust_layers', 'crust_thickness',
                  'crust_at_interface', 'crust_top_depth'):
            result[k] = np.nan

    if wl_m.any():
        # wl_shear_strength (tau_p) is deliberately left as a plain element
        # mean. Thickness-weighting it would be a modelling change, not a bug
        # fix: a crack runs in the weakest sublayer, so the governing value is
        # arguably the minimum rather than any mean, and tau_p cascades into
        # theta, tau_p_star and the trigger ranking. Decide that separately.
        for var in ('grain_size', 'density', 'hand_hardness'):
            try:
                v = ds_t_loc[var].values.ravel()
                result[f'wl_{var}'] = thickness_weighted_mean(v, wl_m & ok, dz)
            except Exception:
                result[f'wl_{var}'] = np.nan

        try:
            v    = ds_t_loc['shear_strength'].values.ravel()
            vals = v[wl_m & ~np.isnan(v)]
            result['wl_shear_strength'] = float(np.nanmean(vals)) if len(vals) else np.nan
        except Exception:
            result['wl_shear_strength'] = np.nan

        wl_z = z[wl_m & ok]
        result['wl_burial_depth'] = float(-wl_z.max())              if len(wl_z) else np.nan
        # D_wl is the SUM of the WL elements' thicknesses. The z span this
        # replaced omitted the basal element, undercounting D_wl by a median
        # 23% on the Jan 18 2026 Little Professor profiles, which propagates as
        # Lambda ~ sqrt(D_wl) and K_wl = G_wl / D_wl.
        result['wl_thickness']    = float(np.nansum(dz[wl_m & ok])) if len(wl_z) else np.nan
        result['interface_z']     = float(interface_z)

        try:
            rc_v    = ds_t_loc[RC_VAR].values.ravel()
            rc_vals = rc_v[wl_m & ~np.isnan(rc_v)]
            result['rc_wl'] = float(np.nanmean(rc_vals)) if len(rc_vals) else np.nan
        except Exception:
            result['rc_wl'] = np.nan
    else:
        for k in ('wl_shear_strength', 'wl_grain_size', 'wl_density',
                  'wl_hand_hardness', 'wl_burial_depth', 'wl_thickness',
                  'interface_z', 'rc_wl'):
            result[k] = np.nan

    rho   = result.get('slab_density',      np.nan)
    h_m   = result.get('slab_thickness',    np.nan)
    D_wl  = result.get('wl_thickness',      np.nan)
    tau_p = result.get('wl_shear_strength', np.nan)

    if all(not np.isnan(v) for v in [rho, h_m, D_wl, tau_p]) and D_wl > 0:
        # Bulk relations, unchanged — _ai.slab_modulus / slab_tensile_strength
        # are the same expressions, lifted out so the per-layer path can reuse
        # them. See their docstrings for the (absent) provenance.
        E_slab  = float(_ai.slab_modulus(rho))
        sigma_t = float(_ai.slab_tensile_strength(rho))
        E_prime = E_slab / (1.0 - nu**2)
        K_wl    = G_WL / D_wl
        Lambda  = float(np.sqrt(E_prime * h_m / K_wl))
        result['E_slab']  = E_slab
        result['sigma_t'] = sigma_t
        result['Lambda']  = Lambda
        result['K_wl']    = K_wl
    else:
        for k in ('E_slab', 'sigma_t', 'Lambda', 'K_wl'):
            result[k] = np.nan

    # --- Layered-slab aggregation (added alongside; nothing above changes) ---
    result.update(_slab_layer_aggregate(ds_t_loc, slab_m & ok, dz, z))
    return result


def _slab_layer_aggregate(ds_t_loc, slab_sel, dz, z) -> dict:
    """Per-layer sigma_t/E aggregation over the slab. All-NaN if unavailable.

    Needs sphericity and dendricity per element; datasets without them degrade
    to NaN rather than failing, so existing outputs are unaffected.
    """
    keys = ('E_eff', 'sigma_t_mean', 'sigma_t_wl', 'f_mean_weighted',
            'wl_ctrl_index', 'wl_ctrl_depth', 'wl_ctrl_thickness',
            'sigma_t_ligament', 'K_Ic_ligament', 'sigma_t_lig_intact',
            'a0_ligament', 'l_ch_ligament', 'n_slab_layers', 'n_rho_clamped')
    try:
        rho = ds_t_loc['density'].values.ravel()[slab_sel]
        sp  = ds_t_loc['sphericity'].values.ravel()[slab_sel]
        dd  = ds_t_loc['dendricity'].values.ravel()[slab_sel]
    except Exception:
        return {k: np.nan for k in keys}
    try:
        gsz = ds_t_loc['grain_size'].values.ravel()[slab_sel]
    except Exception:
        gsz = None   # schweizer2004 needs d_max; it degrades to NaN without it
    return _ai.aggregate_slab(
        rho, dz[slab_sel], sp, dd,
        a=config.FACET_STRENGTH_FACTOR, sp_ref=config.FACET_SP_REF,
        depth_i=-z[slab_sel],
        K_Ic=config.SLAB_K_IC if config.USE_LIGAMENT_BOUND else None,
        k_ic_relation=(config.K_IC_RELATION if config.USE_LIGAMENT_BOUND
                       and config.SLAB_K_IC is None else None),
        ligament_model=config.LIGAMENT_MODEL,
        grain_size_i=gsz, dmax_factor=config.DMAX_FACTOR,
        e_relation=config.E_RELATION,
        e_relations_extra=config.E_RELATIONS_EXTRA,
        k_ic_variants=(config.LIGAMENT_VARIANTS
                       if config.USE_LIGAMENT_BOUND else ()),
        self_consistent_a0=config.SELF_CONSISTENT_A0)


# ---------------------------------------------------------------------------
# Meloche et al. (2025) crack-arrest parameters
# ---------------------------------------------------------------------------

def compute_meloche_features(snap_data: dict, cluster_map: np.ndarray,
                              dem: np.ndarray, transform,
                              snap_ts=None) -> pd.DataFrame:
    """
    Compute Meloche et al. (2025) spatial crack-arrest features.

    snap_data : dict mapping group label → DataFrame indexed by cluster_id,
                with columns from profile_features() (or equivalent CSV).
                Use {'all': features_df} when group labels don't matter.
    cluster_map : 2D int array of cluster IDs matching dem.
    dem         : 2D float array (m), same shape as cluster_map.
    transform   : rasterio Affine transform for dem. Its pixel size converts
                  centroid separations to metres, so theta comes out in Pa/m.
    snap_ts     : unused, kept for API compatibility.

    Returns DataFrame indexed by cluster_id.
    """
    from release_areas.release_geometry import cluster_pixel_index

    PHI_DEG     = config.PHI_DEG
    DELTA       = config.DELTA
    K_NEIGHBORS = config.K_NEIGHBORS
    G_GRAV      = 9.81
    px_m        = abs(transform.a) if transform is not None else 1.0

    fill_dem = np.where(np.isnan(dem), np.nanmean(dem), dem)
    dy, dx   = np.gradient(fill_dem, px_m)
    slope    = np.degrees(np.arctan(np.sqrt(dx**2 + dy**2)))

    pixel_index = cluster_pixel_index(np.where(np.isnan(dem), 0, cluster_map))

    centroids = {}
    slopes_cl = {}
    for cid, (rows, cols) in pixel_index.items():
        centroids[cid] = (float(rows.mean()), float(cols.mean()))
        slopes_cl[cid] = float(slope[rows, cols].mean())

    all_rows = []
    for grp, df in snap_data.items():
        if df.empty:
            continue
        df2 = df.copy()
        df2['group'] = grp
        all_rows.append(df2)

    if not all_rows:
        return pd.DataFrame()

    features = pd.concat(all_rows)
    features = features[~features.index.duplicated(keep='first')]
    if 'wl_shear_strength' not in features.columns:
        return pd.DataFrame()

    cids_arr = np.array([cid for cid in features.index if cid in centroids])
    if len(cids_arr) < 2:
        return pd.DataFrame()

    cents = np.array([centroids[c] for c in cids_arr])

    def _scalar(df, cid, col):
        if cid not in df.index:
            return np.nan
        val = df.loc[cid, col]
        if isinstance(val, pd.Series):
            return float(val.iloc[0])
        return float(val)

    tau_p_arr = np.array([_scalar(features, c, 'wl_shear_strength') for c in cids_arr])

    from sklearn.neighbors import NearestNeighbors
    nbrs = NearestNeighbors(n_neighbors=min(K_NEIGHBORS + 1, len(cids_arr)),
                             algorithm='ball_tree').fit(cents)
    distances, indices = nbrs.kneighbors(cents)

    rows_out = []
    for i, cid in enumerate(cids_arr):
        if cid not in features.index:
            continue
        loc_data = features.loc[cid]
        row = (loc_data.iloc[0].copy() if isinstance(loc_data, pd.DataFrame)
               else loc_data.copy())
        tau_p   = row.get('wl_shear_strength', np.nan)
        Lambda  = row.get('Lambda',            np.nan)
        E_slab  = row.get('E_slab',            np.nan)
        sigma_t = row.get('sigma_t',           np.nan)
        rho     = row.get('slab_density',      np.nan)
        h_m     = row.get('slab_thickness',    np.nan)
        psi_deg = slopes_cl.get(cid, np.nan)

        if any(np.isnan(v) for v in [tau_p, Lambda, rho, h_m, psi_deg]):
            rows_out.append({'cluster_id': cid})
            continue

        psi_rad = np.radians(psi_deg)
        sin_psi = np.sin(psi_rad)
        if sin_psi < 0.01:
            rows_out.append({'cluster_id': cid})
            continue

        tau_g_guard = rho * G_GRAV * h_m * sin_psi
        if tau_g_guard < config.TAU_G_FEATURE_FLOOR:
            rows_out.append({'cluster_id': cid, 'tau_g': tau_g_guard,
                             'slope_angle': psi_deg,
                             'note': 'tau_g below valid range for scaling law'})
            continue

        # Distances come back in pixels; px_m converts them so theta is Pa/m.
        neighbor_idx  = indices[i][1:]
        neighbor_dist = distances[i][1:] * px_m
        theta_vals = []
        for j, d in zip(neighbor_idx, neighbor_dist):
            if d < 1e-6:
                continue
            tau_j = tau_p_arr[j]
            if not np.isnan(tau_j):
                theta_vals.append(abs(tau_p - tau_j) / d)

        theta = float(np.mean(theta_vals)) if theta_vals else np.nan
        if np.isnan(theta) or theta < config.THETA_MIN:
            rows_out.append({'cluster_id': cid, 'tau_g': tau_g_guard,
                             'theta': theta, 'slope_angle': psi_deg})
            continue

        D_wl_val   = row.get('wl_thickness', np.nan)
        if np.isnan(D_wl_val) or D_wl_val <= 0:
            D_wl_val = config.D_WL_FALLBACK
        tau_p0_val = row.get('wl_shear_strength', np.nan)

        ai = _ai_evaluate(
            rho=rho, h=h_m, psi_deg=psi_deg, E=E_slab,
            sigma_t=sigma_t, D_wl=D_wl_val, G_wl=G_WL,
            theta=theta, nu=NU, phi_deg=PHI_DEG, delta=DELTA,
            tau_p0=tau_p0_val if not np.isnan(tau_p0_val) else None,
            speed_ratio_along=config.MODE2_SPEED_RATIO,
            speed_ratio_cross=config.MODE3_SPEED_RATIO,
        )

        tau_g_ai    = ai['tau_g']
        Lambda_ai   = ai['Lambda']
        denom       = theta * Lambda_ai * np.sqrt(1.0 + DELTA)
        Pi1         = tau_g_ai / denom
        Pi2         = Pi1 * np.sqrt(sigma_t / tau_g_ai)

        # Eq. 20 again under the two layered sigma_t bounds, alongside (never
        # replacing) A_ca_brittle. Only sigma_t and E change: k_f and tau_g
        # come from the bulk slab density, so they are shared.
        layered = _layered_arrest(row, ai, tau_g_ai, theta, psi_deg, h_m)

        # A_ca_elastic (JGR Eq. 19) is deliberately not emitted: the paper
        # gives it as a proportionality with no published constant, so an
        # absolute length would be uncalibrated. Pi1_elastic carries the
        # same information in dimensionless form.
        rows_out.append({
            **layered,
            'cluster_id':   cid,
            'slope_angle':  psi_deg,
            'tau_g':        tau_g_ai,
            'theta':        theta,
            'tau_p':        tau_p,
            'Pi1_elastic':  Pi1,
            'Pi2_brittle':  Pi2,
            'Lambda':       Lambda_ai,
            # Mode III elastic length from the current (unpublished)
            # derivation. Written so the data side is ready; propagation only
            # reads it when config.USE_MODE3_LAMBDA is set.
            'Lambda_cross': ai.get('Lambda_cross', np.nan),
            'L_t':          ai.get('L_t',         np.nan),
            # Distance to first slab fracture, along-slope (supershear) and
            # cross-slope (mode III capped at c_s). The ratio is the
            # crack-speed cap applied to lateral propagation.
            'L_dyn':        ai.get('L_dyn',       np.nan),
            'L_dyn_cross':  ai.get('L_dyn_cross', np.nan),
            'mode3_length_ratio': ai.get('mode3_length_ratio', np.nan),
            'A_ca_brittle': ai.get('A_ca',         np.nan),
            'rc_wl':        row.get('rc_wl',       np.nan),
            'G_slab':       ai.get('G_slab',       np.nan),
            'tau_p_star':   ai.get('tau_p_star',   np.nan),
            'R0':           ai.get('R0',            np.nan),
            'A_ca_energy':  ai.get('A_ca_energy',  np.nan),
        })

    if not rows_out:
        return pd.DataFrame()
    return pd.DataFrame(rows_out).set_index('cluster_id')


_LAYERED_KEYS = ('E_eff', 'sigma_t_mean', 'sigma_t_wl', 'sigma_t_ligament',
                 'f_mean_weighted', 'n_slab_layers', 'wl_ctrl_index',
                 'wl_ctrl_depth', 'wl_ctrl_thickness', 'K_Ic_ligament',
                 'sigma_t_lig_intact', 'a0_ligament', 'l_ch_ligament')


def _layered_arrest(row, ai, tau_g, theta, psi_deg, h_m) -> dict:
    """Eq. 20 under the layered sigma_t bounds. NaN when the columns are absent.

    Carries the aggregation columns through from profile_features and adds
    Lambda_eff (from E_eff) plus A_ca under each bound and their ratio.
    """
    # Every layered column present on the row, including the __<variant> and
    # __<E relation> suffixes, so new variants need no change here.
    keys = set(_LAYERED_KEYS)
    try:
        keys |= {k for k in row.index
                 if any(k.startswith(p + '__') for p in _LAYERED_KEYS)
                 or k.startswith('E_eff__')}
    except AttributeError:
        pass
    out = {k: _nan_get(row, k) for k in keys}
    out.update({k: np.nan for k in
                ('Lambda_eff', 'A_ca_mean', 'A_ca_wl', 'A_ca_ratio',
                 'A_ca_ligament')})

    E_eff = out['E_eff']
    D_wl  = _nan_get(row, 'wl_thickness')
    if not (np.isfinite(E_eff) and np.isfinite(D_wl) and D_wl > 0
            and np.isfinite(h_m) and np.isfinite(theta) and theta > 0):
        return out

    rho_bulk = _nan_get(row, 'slab_density')
    out['Lambda_eff'] = float(_ai.elastic_length(E_eff, h_m, D_wl, G_WL, NU))
    k_f = float(_ai.tension_gradient(rho_bulk, psi_deg, config.PHI_DEG))

    def _a_ca(sigma_t, Lam):
        if not (np.isfinite(sigma_t) and sigma_t > 0 and np.isfinite(Lam)):
            return np.nan
        L_t = float(_ai.tensile_length(sigma_t, k_f))
        if not np.isfinite(L_t):
            return np.nan
        return float(_ai.arrest_length(
            tau_g, theta, Lam, sigma_t, config.DELTA, L_t)[0])

    Lam = out['Lambda_eff']
    out['A_ca_mean'] = _a_ca(out['sigma_t_mean'], Lam)
    out['A_ca_wl'] = _a_ca(out['sigma_t_wl'], Lam)
    out['A_ca_ligament'] = _a_ca(out['sigma_t_ligament'], Lam)

    # Per-variant ligament bounds.
    for k in [k for k in out if k.startswith('sigma_t_ligament__')]:
        out['A_ca_ligament__' + k.split('__', 1)[1]] = _a_ca(out[k], Lam)

    # Alternative E relations: only Lambda changes, A_ca ~ E^-0.5.
    for k in [k for k in out if k.startswith('E_eff__')]:
        name = k.split('__', 1)[1]
        if not np.isfinite(out[k]):
            continue
        L2 = float(_ai.elastic_length(out[k], h_m, D_wl, G_WL, NU))
        out['Lambda_eff__' + name] = L2
        out['A_ca_mean__' + name] = _a_ca(out['sigma_t_mean'], L2)

    if np.isfinite(out['A_ca_mean']) and np.isfinite(out['A_ca_wl']) \
            and out['A_ca_wl'] > 0:
        out['A_ca_ratio'] = out['A_ca_mean'] / out['A_ca_wl']
    return out


def _nan_get(row, key):
    v = row.get(key, np.nan)
    if isinstance(v, pd.Series):
        v = v.iloc[0] if len(v) else np.nan
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan
