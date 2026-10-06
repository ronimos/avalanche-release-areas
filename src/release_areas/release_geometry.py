"""
release_geometry.py — Physically-motivated avalanche release polygon generation.

Provides:
    find_stauchwall()            Walk downslope from trigger until slope < threshold
    estimate_cross_slope_width() Cross-slope arrest width from Gaume (2015) / θ
    project_along_aspect()       Project a UTM point distance along slope aspect
    make_release_polygon_2d()    Full 2D release polygon from four constraints:
                                   1. Meloche A_ca (upslope)
                                   2. ~28° slope threshold (downslope / stauchwall)
                                   3. θ cross-slope gradient (flanks)
                                   4. Start zone KML (hard lateral boundary)
    propagate_crack()            Spatial crack propagation: BFS cluster flood-fill
                                   with per-direction arrest criteria
    rasterize_release_polygon()  Burn polygon to a depth raster on the DEM grid
    load_observed_polygon()      Load a GeoJSON polygon and reproject to UTM
    plot_release_comparison()    Meloche polygons vs observed release area figure

No CLI. Used by release_areas.generate_scenarios and release_areas.plot_release.

References
----------
Upslope:   Meloche et al. (2025) JGR Earth Surface, doi:10.1029/2025JF008470
Downslope: Perzl (2007) JRC; Swiss ALIP/PRA; Maggioni & Gruber;
           Bühler et al.; Veitinger et al. (2016)
Flanks:    Gaume et al. (2015) The Cryosphere, doi:10.5194/tc-9-795-2015
           + Meloche (2025) cross-slope θ gradient
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# Default physical parameters
STAUCHWALL_DEG   = 28.0   # slope threshold for downslope arrest (degrees)
FRICTION_DEG     = 27.0   # snow friction angle (degrees), Meloche Table 1
GAUME_ASPECT_CAP = 2.5    # max width/A_ca ratio (prevents runaway cross-slope)
MIN_POLYGON_AREA = 200.0  # m² — discard degenerate polygons smaller than this

USE_MELOCHE_ARREST = False
MELOCHE_DELTA      = 1.0   # softening coefficient δ (Meloche default; their Table 1)


# -----------------------------------------------------------------------
# Polygon cleanup
# -----------------------------------------------------------------------

def fill_polygon_holes(geom, max_hole_frac: float | None = None):
    """Remove interior rings (holes) from a (Multi)Polygon."""
    from shapely.geometry import Polygon, MultiPolygon

    def _fill(poly):
        if not poly.interiors:
            return poly
        if max_hole_frac is None:
            return Polygon(poly.exterior)
        ext_area = Polygon(poly.exterior).area
        keep = [r for r in poly.interiors
                if Polygon(r).area > max_hole_frac * ext_area]
        return Polygon(poly.exterior, keep)

    if geom is None or geom.is_empty:
        return geom
    if geom.geom_type == 'Polygon':
        return _fill(geom)
    if geom.geom_type == 'MultiPolygon':
        return MultiPolygon([_fill(p) for p in geom.geoms])
    return geom


# -----------------------------------------------------------------------
# Terrain helpers
# -----------------------------------------------------------------------

def compute_slope_aspect(dem: np.ndarray,
                         pixel_size: float = 1.0
                         ) -> tuple[np.ndarray, np.ndarray]:
    """Return slope (degrees) and aspect (degrees from N, clockwise) grids."""
    fill = np.where(np.isnan(dem), np.nanmean(dem), dem)
    dy, dx = np.gradient(fill, pixel_size)
    slope  = np.degrees(np.arctan(np.sqrt(dx**2 + dy**2)))
    aspect = np.degrees(np.arctan2(-dx, dy)) % 360.0
    return slope, aspect


def pixel_to_utm(row: int, col: int, transform) -> tuple[float, float]:
    """Convert raster pixel (row, col) to UTM (x, y) using Affine transform."""
    x = transform.c + col * transform.a + row * transform.b
    y = transform.f + col * transform.d + row * transform.e
    return x, y


def utm_to_pixel(x: float, y: float, transform) -> tuple[int, int]:
    """Convert UTM (x, y) to nearest raster pixel (row, col)."""
    col = int(round((x - transform.c) / transform.a))
    row = int(round((y - transform.f) / transform.e))
    return row, col


# -----------------------------------------------------------------------
# Downslope: stauchwall location
# -----------------------------------------------------------------------

def find_stauchwall(trigger_row: int,
                    trigger_col: int,
                    slope_grid: np.ndarray,
                    aspect_grid: np.ndarray,
                    transform,
                    threshold_deg: float = STAUCHWALL_DEG,
                    max_steps: int = 500
                    ) -> tuple[int, int]:
    """Walk downslope from trigger pixel until slope drops below threshold."""
    nrows, ncols = slope_grid.shape
    row, col     = trigger_row, trigger_col

    for _ in range(max_steps):
        if slope_grid[row, col] < threshold_deg:
            return row, col
        asp_rad = np.radians(aspect_grid[row, col])
        dr = int(np.sign(-np.cos(asp_rad)))
        dc = int(np.sign( np.sin(asp_rad)))
        new_row = row + dr
        new_col = col + dc
        if not (0 <= new_row < nrows and 0 <= new_col < ncols):
            break
        if np.isnan(slope_grid[new_row, new_col]):
            break
        row, col = new_row, new_col

    return row, col


# -----------------------------------------------------------------------
# Cross-slope: Gaume (2015) + θ
# -----------------------------------------------------------------------

def estimate_cross_slope_width(
        trigger_cluster_id: int,
        A_ca: float,
        meloche_df: pd.DataFrame,
        cluster_map: np.ndarray,
        transform,
        n_lateral_neighbors: int = 6,
        gaume_aspect_cap: float = GAUME_ASPECT_CAP) -> float:
    """Estimate cross-slope release width (m) using Gaume et al. (2015) / θ."""
    if meloche_df.empty or 'theta' not in meloche_df.columns:
        return A_ca

    try:
        theta_down = float(meloche_df.loc[trigger_cluster_id, 'theta'])
    except (KeyError, TypeError):
        return A_ca

    if np.isnan(theta_down) or theta_down < 1e-6:
        return A_ca

    trig_pixels = np.argwhere(cluster_map == trigger_cluster_id)
    if len(trig_pixels) == 0:
        return A_ca
    t_row, t_col = trig_pixels.mean(axis=0)

    cids_all = np.unique(cluster_map[cluster_map > 0])
    lateral  = []
    for cid in cids_all:
        if cid == trigger_cluster_id:
            continue
        pxs = np.argwhere(cluster_map == cid)
        if len(pxs) == 0:
            continue
        c_row, c_col = pxs.mean(axis=0)
        if abs(c_row - t_row) < 50 and abs(c_col - t_col) > 5:
            lateral.append((cid, abs(c_col - t_col)))

    lateral.sort(key=lambda x: x[1])
    lateral = lateral[:n_lateral_neighbors]

    if not lateral:
        return A_ca

    if 'tau_p' in meloche_df.columns:
        raw = meloche_df.loc[trigger_cluster_id, 'tau_p']
    elif 'wl_shear_strength' in meloche_df.columns:
        raw = meloche_df.loc[trigger_cluster_id, 'wl_shear_strength']
    else:
        raw = np.nan
    tau_trigger = float(pd.to_numeric(raw, errors='coerce'))
    if np.isnan(tau_trigger):
        return A_ca

    theta_cross_vals = []
    for cid, dist_px in lateral:
        if cid not in meloche_df.index:
            continue
        col_name = ('tau_p' if 'tau_p' in meloche_df.columns
                    else 'wl_shear_strength')
        tau_lat = float(pd.to_numeric(meloche_df.loc[cid, col_name], errors='coerce'))
        if not np.isnan(tau_lat) and dist_px > 0:
            theta_cross_vals.append(abs(tau_trigger - tau_lat) / dist_px)

    if not theta_cross_vals:
        return A_ca

    theta_cross = float(np.mean(theta_cross_vals))
    if theta_cross > 1e-6:
        width_factor = min(theta_down / theta_cross, gaume_aspect_cap)
    else:
        width_factor = gaume_aspect_cap

    return float(A_ca * width_factor)


# -----------------------------------------------------------------------
# Upslope projection
# -----------------------------------------------------------------------

def project_along_aspect(x: float, y: float,
                          distance: float,
                          aspect_deg: float,
                          upslope: bool = True) -> tuple[float, float]:
    """Project a UTM point by `distance` metres along (or against) aspect."""
    asp_rad = np.radians(aspect_deg)
    dx = np.sin(asp_rad) * distance
    dy = np.cos(asp_rad) * distance
    sign = -1 if upslope else 1
    return x + sign * dx, y + sign * dy


# -----------------------------------------------------------------------
# Full 2D release polygon
# -----------------------------------------------------------------------

def make_release_polygon_2d(
        trigger_cluster_id: int,
        A_ca: float,
        meloche_df: pd.DataFrame,
        cluster_map: np.ndarray,
        dem: np.ndarray,
        transform,
        start_zone_mask: Optional[np.ndarray] = None,
        snap_features: 'Optional[pd.DataFrame]' = None,
        size_factor: float = 1.0,
        stauchwall_deg: float = STAUCHWALL_DEG,
        mode3_scale: float = 1.5,
        use_propagation: bool = True):
    """
    Build a 2D release polygon.

    Primary method (use_propagation=True):
        BFS cluster flood-fill with per-direction arrest criteria.
    Fallback (use_propagation=False, or propagation returns None):
        Oriented rectangle from A_ca + stauchwall + Gaume cross-slope width.

    Returns shapely Polygon (UTM) or None.
    """
    slope_grid, _ = compute_slope_aspect(dem)

    from shapely.geometry import Polygon, MultiPolygon
    from shapely.ops import unary_union
    import rasterio.features

    sz_polygon = None
    if start_zone_mask is not None:
        shapes = list(rasterio.features.shapes(
            start_zone_mask.astype(np.uint8),
            mask=start_zone_mask, transform=transform))
        if shapes:
            sz_polygon = unary_union(
                [Polygon(s['coordinates'][0]) for s, v in shapes if v == 1])

    stauchwall_mask = slope_grid >= stauchwall_deg
    sw_shapes = list(rasterio.features.shapes(
        stauchwall_mask.astype(np.uint8),
        mask=stauchwall_mask, transform=transform))
    sw_polygon = None
    if sw_shapes:
        sw_polygon = unary_union(
            [Polygon(s['coordinates'][0]) for s, v in sw_shapes if v == 1])
        sw_polygon = sw_polygon.simplify(1.0, preserve_topology=True)

    def _clip_polygon(poly):
        if poly is None or poly.is_empty:
            return None
        if sz_polygon is not None:
            poly = poly.intersection(sz_polygon)
        if sw_polygon is not None:
            poly = poly.intersection(sw_polygon)
        if poly.is_empty:
            return None
        if poly.geom_type == 'GeometryCollection':
            polys = [g for g in poly.geoms
                     if g.geom_type in ('Polygon', 'MultiPolygon')]
            if not polys:
                return None
            poly = unary_union(polys)
        if poly.geom_type == 'MultiPolygon':
            poly = max(poly.geoms, key=lambda p: p.area)
        if poly.is_empty or poly.area < MIN_POLYGON_AREA:
            return None
        return poly

    if use_propagation and not meloche_df.empty:
        polygon, failed = propagate_crack(
            trigger_cluster_id = trigger_cluster_id,
            meloche_df         = meloche_df,
            cluster_map        = cluster_map,
            dem                = dem,
            slope_grid         = slope_grid,
            transform          = transform,
            start_zone_mask    = start_zone_mask,
            stauchwall_deg     = stauchwall_deg,
            mode3_scale        = mode3_scale,
            size_factor        = size_factor,
            A_ca               = A_ca,
            snap_features      = snap_features,
        )
        if polygon is not None:
            polygon = _clip_polygon(polygon)
            if polygon is not None:
                return polygon
        print(f"  propagate_crack returned None for cluster {trigger_cluster_id}"
              f" — falling back to rectangle")

    _, aspect_grid = compute_slope_aspect(dem)
    trig_px = np.argwhere(cluster_map == trigger_cluster_id)
    if len(trig_px) == 0:
        return None
    t_row, t_col = trig_px.mean(axis=0).astype(int)
    t_x, t_y    = pixel_to_utm(t_row, t_col, transform)
    t_aspect    = float(aspect_grid[t_row, t_col])

    up_x, up_y  = project_along_aspect(t_x, t_y, A_ca * size_factor,
                                        t_aspect, upslope=True)
    sw_row, sw_col = find_stauchwall(t_row, t_col, slope_grid, aspect_grid,
                                     transform, threshold_deg=stauchwall_deg)
    sw_x, sw_y  = pixel_to_utm(sw_row, sw_col, transform)
    half_width  = estimate_cross_slope_width(
        trigger_cluster_id, A_ca, meloche_df, cluster_map, transform
    ) * size_factor / 2.0

    asp_rad = np.radians(t_aspect)
    fall_x, fall_y   = np.sin(asp_rad), np.cos(asp_rad)
    cross_x, cross_y = -fall_y, fall_x
    corners = [
        (up_x - cross_x * half_width, up_y - cross_y * half_width),
        (up_x + cross_x * half_width, up_y + cross_y * half_width),
        (sw_x + cross_x * half_width, sw_y + cross_y * half_width),
        (sw_x - cross_x * half_width, sw_y - cross_y * half_width),
    ]
    polygon = Polygon(corners)
    if not polygon.is_valid:
        polygon = polygon.buffer(0)

    return _clip_polygon(polygon)


# -----------------------------------------------------------------------
# Spatial crack propagation
# -----------------------------------------------------------------------

def propagate_crack(
        trigger_cluster_id: int,
        meloche_df: pd.DataFrame,
        cluster_map: np.ndarray,
        dem: np.ndarray,
        slope_grid: np.ndarray,
        transform,
        start_zone_mask: Optional[np.ndarray] = None,
        stauchwall_deg: float = STAUCHWALL_DEG,
        mode3_scale: float = 2.0,
        size_factor: float = 1.0,
        A_ca: Optional[float] = None,
        snap_features: 'Optional[pd.DataFrame]' = None,
        k_neighbours: int = 8,
        max_clusters: int = 500,
        wave_callback=None):
    """
    Identify release zone as the connected region of clusters reachable
    from the trigger where crack arrest criteria are not met.

    Starting from the trigger cluster, expands outward to neighbouring
    clusters in order of proximity. A cluster is included if it passes
    distance caps (upslope A_ca, downslope stauchwall, lateral Gaume width),
    slope, slab thickness, and elastic-length continuity checks.
    """
    from collections import deque
    from shapely.geometry import MultiPolygon
    from shapely.ops import unary_union
    import rasterio.features

    if meloche_df.empty or trigger_cluster_id not in meloche_df.index:
        return None, set()

    loc = meloche_df.loc[trigger_cluster_id]
    row = loc.iloc[0] if isinstance(loc, pd.DataFrame) else loc

    col = 'Pi1_elastic' if 'Pi1_elastic' in meloche_df.columns else None
    if col is None:
        return None, set()

    pi1_raw = row[col]
    pi1_trigger = float(pi1_raw.iloc[0] if isinstance(pi1_raw, pd.Series)
                        else pi1_raw)
    if np.isnan(pi1_trigger) or pi1_trigger <= 0:
        return None, set()

    TAU_G_ABS_FLOOR = 350.0   # Pa

    print(f"    Pi1_trigger={pi1_trigger:.3f}  "
          f"tau_g_floor={TAU_G_ABS_FLOOR:.0f}Pa  "
          f"size_factor={size_factor:.2f}")

    def _scalar(cid, c):
        if cid not in meloche_df.index:
            return np.nan
        v = meloche_df.loc[cid, c]
        if isinstance(v, pd.DataFrame): v = v.iloc[0][c]
        elif isinstance(v, pd.Series):  v = v.iloc[0]
        return float(v)

    def _mean_slope(cid):
        pxs = np.argwhere(cluster_map == cid)
        if len(pxs) == 0:
            return 0.0
        return float(slope_grid[pxs[:, 0], pxs[:, 1]].mean())

    t_pxs = np.argwhere(cluster_map == trigger_cluster_id)
    if len(t_pxs) == 0:
        return None, set()
    t_row  = int(t_pxs.mean(axis=0)[0])
    t_col  = int(t_pxs.mean(axis=0)[1])
    t_x, t_y = pixel_to_utm(t_row, t_col, transform)
    _, aspect_grid = compute_slope_aspect(dem)
    t_aspect = float(aspect_grid[t_row, t_col])

    d_up = (max(A_ca, 20.0) if A_ca else 50.0) * size_factor

    sw_row, sw_col = find_stauchwall(
        t_row, t_col, slope_grid, aspect_grid, transform,
        threshold_deg=stauchwall_deg)
    sw_x, sw_y = pixel_to_utm(sw_row, sw_col, transform)
    d_down = max(float(np.sqrt((sw_x - t_x)**2 + (sw_y - t_y)**2)),
                 20.0) * size_factor

    print(f"    distance caps: up={d_up:.0f}m  down={d_down:.0f}m")

    from sklearn.neighbors import NearestNeighbors
    cids = np.array([c for c in np.unique(cluster_map) if c > 0])
    pxs  = np.array([np.argwhere(cluster_map == c).mean(axis=0) for c in cids])
    k_actual = min(k_neighbours + 1, len(cids))
    nbrs     = NearestNeighbors(n_neighbors=k_actual).fit(pxs)
    _, idxs  = nbrs.kneighbors(pxs)
    neighbours = {int(cids[i]): [int(cids[j]) for j in idxs[i][1:]]
                  for i in range(len(cids))}

    centroids = {int(c): pixel_to_utm(int(pxs[i][0]), int(pxs[i][1]), transform)
                 for i, c in enumerate(cids)}

    tau_p_trigger = None
    if snap_features is not None and trigger_cluster_id in snap_features.index:
        raw = snap_features.loc[trigger_cluster_id, 'wl_shear_strength'] \
              if 'wl_shear_strength' in snap_features.columns else np.nan
        if isinstance(raw, pd.Series): raw = raw.iloc[0]
        v = float(raw)
        tau_p_trigger = v if not np.isnan(v) and v > 0 else None
    if tau_p_trigger is None:
        print(f"    Warning: no tau_p for trigger {trigger_cluster_id}, "
              f"WL strength criterion disabled")

    d_lat = estimate_cross_slope_width(
        trigger_cluster_id, A_ca if A_ca else 50.0,
        meloche_df, cluster_map, transform) * size_factor
    d_lat = max(d_lat, 15.0)

    lambda_trigger = None
    if not meloche_df.empty and 'Lambda' in meloche_df.columns:
        lam_vals = meloche_df['Lambda'].dropna()
        if len(lam_vals):
            lam_raw = float(lam_vals.median())
            if lam_raw > 0:
                lambda_trigger = lam_raw
                print(f"    Lambda_median={lambda_trigger:.2f}m  "
                      f"(trigger own={_scalar(trigger_cluster_id, 'Lambda'):.2f}m)")

    tau_p_str = f"{tau_p_trigger:.0f} Pa" if tau_p_trigger else "N/A"
    print(f"    tau_p_trigger={tau_p_str}  d_lat={d_lat:.0f}m")

    LAMBDA_DROP_FACTOR     = 0.20
    LAMBDA_RISE_FACTOR     = 0.30
    THICKNESS_DROP_FACTOR  = 0.20
    THICKNESS_RISE_FACTOR  = 0.30
    MIN_PROPAGATION_SLAB   = 0.50
    MIN_PROPAGATION_LAMBDA = 0.1

    def _get_props(cid):
        props = {'cid': cid}
        if not meloche_df.empty and cid in meloche_df.index:
            for c in ['Lambda', 'tau_g', 'rc_wl', 'L_t']:
                if c in meloche_df.columns:
                    v = meloche_df.loc[cid, c]
                    if isinstance(v, pd.DataFrame): v = v.iloc[0][c]
                    elif isinstance(v, pd.Series):  v = v.iloc[0]
                    props[c] = float(v)
        if snap_features is not None and cid in snap_features.index:
            for c in ['slab_thickness', 'slab_density',
                      'wl_shear_strength', 'sigma_t']:
                if c in snap_features.columns:
                    v = snap_features.loc[cid, c]
                    if isinstance(v, pd.Series): v = v.iloc[0]
                    try:
                        props[c] = float(v)
                    except (ValueError, TypeError):
                        pass
        return props

    trigger_props = _get_props(trigger_cluster_id)

    def _qualifies(cid, current_props):
        if start_zone_mask is not None:
            px = np.argwhere(cluster_map == cid)
            if len(px) == 0:
                return False, 'no_data'
            r, c = int(px.mean(axis=0)[0]), int(px.mean(axis=0)[1])
            if not start_zone_mask[r, c]:
                return False, 'outside_start_zone'

        xy = centroids.get(cid)
        if xy is None:
            return False, 'no_data'
        asp_rad        = np.radians(t_aspect)
        fall_x, fall_y = np.sin(asp_rad), np.cos(asp_rad)
        dx, dy         = xy[0] - t_x, xy[1] - t_y
        dist           = np.sqrt(dx**2 + dy**2)
        dot            = (dx * fall_x + dy * fall_y) / max(dist, 1e-6)

        is_downslope = dot >  0.707
        is_upslope   = dot < -0.707

        if is_upslope   and dist > d_up:   return False, 'upslope_distance_cap'
        if is_downslope and dist > d_down: return False, 'downslope_distance_cap'
        if not is_upslope and not is_downslope and dist > d_lat:
            return False, 'lateral_distance_cap'

        if not is_upslope and _mean_slope(cid) < stauchwall_deg:
            return False, 'stauchwall_slope'

        nbr_props = _get_props(cid)
        tau_g_nbr = nbr_props.get('tau_g', np.nan)
        if not np.isnan(tau_g_nbr) and tau_g_nbr < TAU_G_ABS_FLOOR:
            return False, 'tau_g_below_floor'

        h_nbr = nbr_props.get('slab_thickness', np.nan)
        if not np.isnan(h_nbr) and h_nbr < MIN_PROPAGATION_SLAB:
            return False, 'thin_slab'

        lam_nbr_abs = nbr_props.get('Lambda', np.nan)
        if not np.isnan(lam_nbr_abs) and lam_nbr_abs < MIN_PROPAGATION_LAMBDA:
            return False, 'low_lambda'

        if USE_MELOCHE_ARREST:
            cur_cid = current_props.get('cid')
            xy_c    = centroids.get(cur_cid)
            xy_n    = centroids.get(cid)
            tau_p_c = current_props.get('wl_shear_strength', np.nan)
            tau_p_n = nbr_props.get('wl_shear_strength', np.nan)
            if (xy_c is not None and xy_n is not None
                    and not np.isnan(tau_p_c) and not np.isnan(tau_p_n)):
                dist_cn = float(np.hypot(xy_n[0] - xy_c[0], xy_n[1] - xy_c[1]))
                if dist_cn > 1e-6:
                    theta_dir = max(0.0, (tau_p_n - tau_p_c) / dist_cn)
                    if theta_dir > 0.0:
                        Lam   = nbr_props.get('Lambda',  np.nan)
                        tg    = tau_g_nbr
                        sig_t = nbr_props.get('sigma_t', np.nan)
                        L_t   = nbr_props.get('L_t',     np.nan)
                        if (not any(np.isnan(v) for v in (Lam, tg, sig_t, L_t))
                                and tg > 0.0):
                            Pi     = tg / (theta_dir * Lam
                                           * np.sqrt(1.0 + MELOCHE_DELTA))
                            A_ca_d = L_t * Pi * np.sqrt(sig_t / tg)
                            if A_ca_d < dist_cn:
                                return False, 'meloche_arrest'
            return True, 'propagated'

        lam_current = current_props.get('Lambda', np.nan)
        lam_nbr     = nbr_props.get('Lambda', np.nan)
        if (not np.isnan(lam_current) and not np.isnan(lam_nbr)
                and lam_current > 0):
            signed_change = (lam_nbr - lam_current) / lam_current
            if signed_change < 0:
                if -signed_change > LAMBDA_DROP_FACTOR * size_factor:
                    return False, 'lambda_discontinuity_drop'
            else:
                if signed_change > LAMBDA_RISE_FACTOR * size_factor:
                    return False, 'lambda_discontinuity_rise'

        h_current = current_props.get('slab_thickness', np.nan)
        if (not np.isnan(h_current) and not np.isnan(h_nbr)
                and h_current > 0):
            signed_change = (h_nbr - h_current) / h_current
            if signed_change < 0:
                if -signed_change > THICKNESS_DROP_FACTOR * size_factor:
                    return False, 'thickness_discontinuity_drop'
            else:
                if signed_change > THICKNESS_RISE_FACTOR * size_factor:
                    return False, 'thickness_discontinuity_rise'

        return True, 'propagated'

    qual_ok, qual_reason = _qualifies(trigger_cluster_id, trigger_props)
    if not qual_ok:
        pxs_t = np.argwhere(cluster_map == trigger_cluster_id)
        if len(pxs_t) == 0:
            return None, set()

    failed         = {trigger_cluster_id}
    cluster_props  = {trigger_cluster_id: trigger_props}
    visited        = {trigger_cluster_id}
    current_ring   = [trigger_cluster_id]
    ring_idx       = 0

    if wave_callback is not None:
        wave_callback(0, [trigger_cluster_id], [])

    while current_ring and len(failed) < max_clusters:
        ring_idx += 1
        next_ring = []
        ring_arrested = []

        for current in current_ring:
            current_props = cluster_props.get(current, trigger_props)
            for nbr in neighbours.get(current, []):
                if nbr in visited:
                    continue
                visited.add(nbr)
                qualifies, reason = _qualifies(nbr, current_props)
                if qualifies:
                    failed.add(nbr)
                    cluster_props[nbr] = _get_props(nbr)
                    next_ring.append(nbr)
                else:
                    ring_arrested.append({'cid': nbr, 'reason': reason})

        if not next_ring and not ring_arrested:
            break

        if wave_callback is not None:
            wave_callback(ring_idx, next_ring, ring_arrested)

        current_ring = next_ring

    dir_counts = {'upslope': 0, 'downslope': 0, 'lateral': 0, 'other': 0}
    for cid in failed:
        if cid == trigger_cluster_id:
            continue
        xy = centroids.get(cid)
        if xy is None:
            continue
        asp_rad = np.radians(t_aspect)
        fall_x, fall_y = np.sin(asp_rad), np.cos(asp_rad)
        dx, dy = xy[0] - t_x, xy[1] - t_y
        dist   = np.sqrt(dx**2 + dy**2)
        dot    = (dx * fall_x + dy * fall_y) / max(dist, 1e-6)
        if   dot < -0.707: dir_counts['upslope']  += 1
        elif dot >  0.707: dir_counts['downslope'] += 1
        else:              dir_counts['lateral']   += 1
    print(f"  Connected region: {len(failed)} clusters "
          f"(up={dir_counts['upslope']} "
          f"down={dir_counts['downslope']} "
          f"lat={dir_counts['lateral']})")

    if len(failed) < 2:
        return None, failed

    mask   = np.isin(cluster_map, list(failed)).astype(np.uint8)
    shapes = list(rasterio.features.shapes(
        mask, mask=mask.astype(bool), transform=transform))
    if not shapes:
        return None, failed

    polys   = [__import__('shapely.geometry', fromlist=['shape']).shape(s)
               for s, v in shapes if v == 1]
    polygon = unary_union(polys)
    if polygon.geom_type == 'MultiPolygon':
        polygon = max(polygon.geoms, key=lambda p: p.area)
    if polygon.is_empty or polygon.area < MIN_POLYGON_AREA:
        return None, failed

    polygon = fill_polygon_holes(polygon)

    return polygon, failed


# -----------------------------------------------------------------------
# Depth raster
# -----------------------------------------------------------------------

def rasterize_release_polygon(polygon,
                               depth_grid: np.ndarray,
                               dem_shape: tuple,
                               transform) -> np.ndarray:
    """Burn release polygon onto the depth grid. Returns float32 raster."""
    import rasterio.features

    if polygon is None or polygon.is_empty:
        return np.full(dem_shape, np.nan, dtype=np.float32)

    mask = rasterio.features.geometry_mask(
        [polygon.__geo_interface__],
        out_shape=dem_shape,
        transform=transform,
        invert=True)

    return np.where(mask & ~np.isnan(depth_grid),
                    depth_grid.astype(np.float32), np.nan)


# -----------------------------------------------------------------------
# GeoJSON loader (shared by scripts)
# -----------------------------------------------------------------------

def load_observed_polygon(path, dst_epsg: int = 32613):
    """Load a GeoJSON polygon file and reproject to UTM (default EPSG:32613)."""
    import json
    import re
    from shapely.geometry import shape, Polygon, MultiPolygon
    from shapely.ops import unary_union
    from pyproj import Transformer

    with open(str(path)) as f:
        gj = json.load(f)

    polys = [shape(feat['geometry']) for feat in gj['features']]
    merged = unary_union(polys)

    src_epsg = 4326
    crs_node = gj.get('crs', {}).get('properties', {}).get('name', '')
    m = re.search(r'EPSG:+(\d+)', crs_node, re.IGNORECASE)
    if m:
        src_epsg = int(m.group(1))

    if src_epsg != dst_epsg:
        tr = Transformer.from_crs(f'EPSG:{src_epsg}', f'EPSG:{dst_epsg}', always_xy=True)

        def _reproj(poly):
            if poly.geom_type == 'Polygon':
                xs, ys = zip(*poly.exterior.coords)
                xs2, ys2 = tr.transform(xs, ys)
                return Polygon(zip(xs2, ys2))
            return MultiPolygon([_reproj(p) for p in poly.geoms])

        merged = _reproj(merged)

    return merged


def _label_from_release_path(path) -> str:
    """Derive a human-readable label from an avalanche_release_area_YYYYMMDD.geojson path."""
    import re
    m = re.search(r'(\d{8})', Path(path).stem)
    if m:
        d = m.group(1)
        return f"Observed {d[:4]}-{d[4:6]}-{d[6:8]}"
    return Path(path).stem


def load_observed_polygons(primary_path, dst_epsg: int = 32613) -> list:
    """Load all avalanche_release_area_*.geojson files from the same directory as primary_path.

    Returns a list of (polygon, label) tuples sorted by filename so that new
    release areas added to the boundaries directory appear automatically.
    """
    if primary_path is None:
        return []
    parent = Path(primary_path).parent
    entries = []
    for p in sorted(parent.glob('avalanche_release_area_*.geojson')):
        try:
            poly = load_observed_polygon(p, dst_epsg)
            entries.append((poly, _label_from_release_path(p)))
        except Exception as e:
            print(f"Warning: could not load {p.name}: {e}")
    return entries


# -----------------------------------------------------------------------
# Comparison plot
# -----------------------------------------------------------------------

_OBS_PALETTE = ['#E31A1C', '#FF7F00', '#6A3D9A', '#8B4513']


def _poly_exterior_xy(geom):
    """Return exterior (xs, ys), falling back to largest part for MultiPolygon."""
    if geom.geom_type == 'MultiPolygon':
        geom = max(geom.geoms, key=lambda p: p.area)
    return geom.exterior.xy


def plot_release_comparison(
        meloche_polygons: list,
        observed_polygons: list,
        dem: np.ndarray,
        transform,
        start_zone_mask=None,
        trigger_labels=None,
        trigger_centroids=None,
        trigger_ious=None,
        most_likely_polygon=None,
        most_likely_label=None,
        out_path=None,
        title: str = "Release polygon comparison") -> None:
    """Plot Meloche-derived polygons vs observed release area(s).

    observed_polygons — list of (shapely_geometry, label_str) tuples.
    All matching avalanche_release_area_*.geojson files in the boundaries
    directory are loaded and rendered so new ones appear automatically.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe
    from matplotlib.colors import LightSource
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    fig, ax = plt.subplots(figsize=(10, 10))

    fill_dem  = np.where(np.isnan(dem), np.nanmean(dem), dem)
    hillshade = LightSource(azdeg=315, altdeg=45).hillshade(fill_dem, dx=1.0, dy=1.0)
    nrows, ncols = dem.shape
    extent = [transform.c,
              transform.c + ncols * transform.a,
              transform.f + nrows * transform.e,
              transform.f]
    ax.imshow(hillshade, cmap='gray', extent=extent,
              alpha=0.6, aspect='auto', origin='upper')

    if start_zone_mask is not None:
        ax.contour(start_zone_mask.astype(float), levels=[0.5],
                   colors=['limegreen'], linewidths=2.0, alpha=0.9,
                   extent=extent, origin='upper')

    obs_area = 0.0
    for j, (obs_poly, obs_label) in enumerate(observed_polygons or []):
        if obs_poly is None or obs_poly.is_empty:
            continue
        oc = _OBS_PALETTE[j % len(_OBS_PALETTE)]
        xs, ys = _poly_exterior_xy(obs_poly)
        ax.fill(xs, ys, alpha=0.30, color=oc, zorder=3)
        ax.plot(xs, ys, color=oc, linewidth=2.5, zorder=4)
        if j == 0:
            obs_area = obs_poly.area

    colors    = plt.cm.tab10.colors
    mel_areas = []

    for i, (poly, size_f) in enumerate(meloche_polygons):
        if poly is None or poly.is_empty:
            continue
        color = colors[i % len(colors)]
        xs, ys = _poly_exterior_xy(poly)
        ax.fill(xs, ys, alpha=0.20, color=color, zorder=2)
        ax.plot(xs, ys, color=color, linewidth=1.8, alpha=0.85, zorder=5)
        mel_areas.append(poly.area)

        if trigger_centroids is not None and i < len(trigger_centroids):
            cx, cy = trigger_centroids[i]
        else:
            cx, cy = poly.centroid.x, poly.centroid.y
        ax.plot(cx, cy, marker='*', markersize=14, color=color,
                markeredgecolor='white', markeredgewidth=0.8, zorder=10)

        short = f"T{i+1}"
        txt   = ax.text(cx + 8, cy + 8, short, fontsize=9,
                        color=color, fontweight='bold', zorder=11)
        txt.set_path_effects([pe.Stroke(linewidth=2.5, foreground='white'),
                               pe.Normal()])

    ml_area = 0.0
    if most_likely_polygon is not None and not most_likely_polygon.is_empty:
        ml_area = most_likely_polygon.area
        mxs, mys = _poly_exterior_xy(most_likely_polygon)
        ax.plot(mxs, mys, color='white', linewidth=4.2, zorder=7)
        ax.plot(mxs, mys, color='black', linewidth=2.4, linestyle=(0, (6, 3)),
                zorder=8)

    if mel_areas and obs_area:
        ratio = float(np.median(mel_areas)) / obs_area
        stats = ("Observed:       %6.0f m2\n"
                 "Meloche P50:    %6.0f m2\n"
                 "Ratio Mel/Obs:    %.2f") % (obs_area, float(np.median(mel_areas)), ratio)
    elif mel_areas:
        stats = "Meloche P50: %.0f m2" % float(np.median(mel_areas))
    else:
        stats = "No polygons generated"
    if ml_area:
        stats += "\nMost likely:    %6.0f m2" % ml_area
        if obs_area:
            stats += "\n  ML/Obs ratio:   %.2f" % (ml_area / obs_area)
    if trigger_ious and any(v > 0 for v in trigger_ious):
        best_i = int(np.argmax(trigger_ious))
        stats += "\nBest IoU:         %.3f  (T%d)" % (trigger_ious[best_i], best_i + 1)
    ax.text(0.02, 0.02, stats, transform=ax.transAxes, fontsize=8.5,
            verticalalignment='bottom', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='white',
                      edgecolor='gray', alpha=0.85))

    handles = []
    for j, (obs_poly, obs_label) in enumerate(observed_polygons or []):
        if obs_poly is None or obs_poly.is_empty:
            continue
        oc = _OBS_PALETTE[j % len(_OBS_PALETTE)]
        handles.append(Line2D([0], [0], color=oc, linewidth=2.5,
                              label=f"{obs_label}  ({obs_poly.area:.0f} m²)"))
    handles.append(Line2D([0],[0], color='limegreen', linewidth=2.0, label='Start zone'))
    if ml_area:
        handles.append(Line2D([0],[0], color='black', linewidth=2.4,
                              linestyle=(0, (6, 3)),
                              label=(most_likely_label
                                     or "Most likely  (%.0f m2)" % ml_area)))
    for i, (poly, size_f) in enumerate(meloche_polygons):
        if poly is None or poly.is_empty:
            continue
        color = colors[i % len(colors)]
        lbl   = trigger_labels[i] if trigger_labels else ("T%d" % (i+1))
        iou_str = ("  IoU=%.3f" % trigger_ious[i]
                   if trigger_ious and i < len(trigger_ious) else "")
        handles.append(Patch(facecolor=color, alpha=0.5, edgecolor=color,
                             label="%s  (%.0f m2)%s" % (lbl, poly.area, iou_str)))

    ax.legend(handles=handles, loc='upper right', fontsize=7.5, framealpha=0.90)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel('Easting (m UTM)')
    ax.set_ylabel('Northing (m UTM)')
    ax.ticklabel_format(style='plain', axis='both')
    plt.tight_layout()

    if out_path:
        fig.savefig(str(out_path), dpi=150, bbox_inches='tight')
        print(f"Release comparison plot: {out_path}")
        plt.close(fig)
