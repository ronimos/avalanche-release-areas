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

Mode III (cross-slope) arrest is NOT applied here. arrest_indices provides
elastic_length_cross and slab_energy_cap_cross, but no mode III multiplier is
used pending Johan Gaume's antiplane formula — see CLAUDE.md "Mode III TODO".
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from release_areas import config
from release_areas.arrest_indices import arrest_length

# Re-exported for backwards compatibility; config is the source of truth.
STAUCHWALL_DEG     = config.STAUCHWALL_DEG
FRICTION_DEG       = config.PHI_DEG
GAUME_ASPECT_CAP   = config.GAUME_ASPECT_CAP
MIN_POLYGON_AREA   = config.MIN_POLYGON_AREA
USE_MELOCHE_ARREST = config.USE_MELOCHE_ARREST
MELOCHE_DELTA      = config.DELTA


def mode3_speed_cap_ratio(trigger_cluster_id: int,
                          meloche_df: pd.DataFrame) -> Optional[float]:
    """L_dyn,III / L_dyn,II for the trigger, or None if the cap is disabled.

    Prefers the `mode3_length_ratio` column written by
    compute_meloche_features; falls back to computing it from nu so CSVs
    predating the column still get the cap. Returns None when
    config.USE_MODE3_SPEED_CAP is off, which restores pre-cap behaviour.
    """
    if not config.USE_MODE3_SPEED_CAP:
        return None

    if (not meloche_df.empty and 'mode3_length_ratio' in meloche_df.columns
            and trigger_cluster_id in meloche_df.index):
        ratio = float(first_scalar(
            meloche_df.loc[trigger_cluster_id, 'mode3_length_ratio']))
        if np.isfinite(ratio) and 0.0 < ratio <= 1.0:
            return ratio

    # E and rho cancel out of the ratio, so any positive pair works here.
    from release_areas.arrest_indices import mode3_length_ratio
    return float(mode3_length_ratio(
        4.0e6, 250.0, config.NU,
        config.MODE2_SPEED_RATIO, config.MODE3_SPEED_RATIO))


def directional_lambda(props: dict, is_cross: bool) -> float:
    """Elastic length Λ to use for propagation in a given direction.

    Along-slope (mode II) always uses `Lambda`. Cross-slope (mode III) uses
    `Lambda_cross` only when `config.USE_MODE3_LAMBDA` is set; otherwise it
    falls back to the mode II value, which is the current default because
    Gaume's antiplane formula is not published yet.

    This is the ONLY place the two directions diverge. When Λ_III lands,
    update `arrest_indices.elastic_length_cross()` and set USE_MODE3_LAMBDA —
    no other change should be needed.
    """
    if is_cross and config.USE_MODE3_LAMBDA:
        lam = props.get('Lambda_cross', np.nan)
        if not np.isnan(lam):
            return lam
    return props.get('Lambda', np.nan)


def first_scalar(value):
    """Unwrap a pandas lookup that returned a Series/DataFrame to one scalar.

    Cluster ids are not guaranteed unique across concatenated feature frames,
    so .loc[cid, col] can yield a Series rather than a scalar.
    """
    if isinstance(value, pd.DataFrame):
        return value.iloc[0, 0]
    if isinstance(value, pd.Series):
        return value.iloc[0]
    return value


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


def cluster_pixel_index(cluster_map: np.ndarray) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """Map every positive cluster id to its (rows, cols) pixel arrays.

    Built in a single sort instead of one full-raster scan per cluster, which
    is what made the per-cluster `np.argwhere(cluster_map == cid)` calls the
    dominant cost of polygon generation.
    """
    flat    = cluster_map.ravel()
    valid   = np.flatnonzero(flat > 0)
    if len(valid) == 0:
        return {}
    ids     = flat[valid].astype(np.int64)
    order   = np.argsort(ids, kind='stable')
    ids     = ids[order]
    linear  = valid[order]
    splits  = np.flatnonzero(np.diff(ids)) + 1
    ncols   = cluster_map.shape[1]
    return {int(ids[seg[0]]): (linear[seg] // ncols, linear[seg] % ncols)
            for seg in np.split(np.arange(len(ids)), splits) if len(seg)}


def cluster_centroids_px(pixel_index: dict) -> dict[int, tuple[float, float]]:
    """Per-cluster centroid in pixel (row, col) coordinates."""
    return {cid: (float(rows.mean()), float(cols.mean()))
            for cid, (rows, cols) in pixel_index.items()}


# -----------------------------------------------------------------------
# Downslope: stauchwall location
# -----------------------------------------------------------------------

def find_stauchwall(trigger_row: int,
                    trigger_col: int,
                    slope_grid: np.ndarray,
                    aspect_grid: np.ndarray,
                    threshold_deg: float = STAUCHWALL_DEG,
                    max_steps: int = 500
                    ) -> tuple[int, int]:
    """Walk downslope from trigger pixel until slope drops below threshold.

    Each step advances one pixel width along the *local* aspect, carrying
    sub-pixel position across iterations. Taking `int(np.sign(...))` of the
    direction components instead would quantise every step to a pure
    diagonal — sign() is ±1 for any nonzero component, so only the four
    bearings 45/135/225/315 are reachable and the walk becomes a straight
    ray. On the Little Professor ESE start zone (aspect 110-120 deg) that
    quantisation sent every trigger off at 135 deg, a systematic +14 deg
    median clockwise error that rotated the release polygons toward south.
    """
    nrows, ncols = slope_grid.shape
    row, col     = trigger_row, trigger_col
    row_f, col_f = float(trigger_row), float(trigger_col)

    for _ in range(max_steps):
        if slope_grid[row, col] < threshold_deg:
            return row, col
        asp_rad = np.radians(aspect_grid[row, col])
        # Downslope unit vector, map frame -> pixel frame:
        # east = +sin(aspect) -> +col, north = +cos(aspect) -> -row.
        row_f += -np.cos(asp_rad)
        col_f +=  np.sin(asp_rad)
        new_row = int(round(row_f))
        new_col = int(round(col_f))
        if not (0 <= new_row < nrows and 0 <= new_col < ncols):
            break
        if np.isnan(slope_grid[new_row, new_col]):
            break
        row, col = new_row, new_col

    return row, col


# -----------------------------------------------------------------------
# Cross-slope: Gaume (2015) + θ
# -----------------------------------------------------------------------

def _tau_p_lookup(meloche_df: pd.DataFrame,
                  snap_features: Optional[pd.DataFrame]) -> tuple:
    """Pick the frame/column holding weak-layer shear strength τp.

    compute_meloche_features() does not carry τp through to its output, so the
    value normally lives in the profile_features frame. Checking meloche_df
    first keeps hand-assembled frames that do carry it working.
    """
    for frame in (meloche_df, snap_features):
        if frame is None or frame.empty:
            continue
        for col in ('tau_p', 'wl_shear_strength'):
            if col in frame.columns:
                return frame, col
    return None, None


def _sector_theta(offsets, axis_x: float, axis_y: float,
                  tau_trigger: float, tau_frame, tau_col: str,
                  n_neighbors: int) -> float:
    """Mean |Δτp|/d over the nearest clusters lying along a given axis.

    `offsets` is a prebuilt list of (cid, dx, dy, dist_m) already filtered to
    the [MIN_SEP, MAX_SEP] separation band. A cluster qualifies when its
    separation points within THETA_CROSS_SECTOR_DEG of ±(axis_x, axis_y) —
    |d·axis|, so both halves of the axis count.

    Used for both the along-slope and cross-slope θ so the ratio that sets the
    Gaume width is formed from two estimates made the same way, at the same
    lag, from the same τp field.
    """
    cos_min = np.cos(np.radians(config.THETA_CROSS_SECTOR_DEG))
    sel = sorted(((cid, d) for cid, dx, dy, d in offsets
                  if abs(dx * axis_x + dy * axis_y) / d >= cos_min),
                 key=lambda t: t[1])[:n_neighbors]

    vals = []
    for cid, dist_m in sel:
        if cid not in tau_frame.index:
            continue
        tau_j = float(pd.to_numeric(
            first_scalar(tau_frame.loc[cid, tau_col]), errors='coerce'))
        if not np.isnan(tau_j) and dist_m > 0:
            vals.append(abs(tau_trigger - tau_j) / dist_m)
    return float(np.mean(vals)) if vals else np.nan


def estimate_cross_slope_width(
        trigger_cluster_id: int,
        A_ca: float,
        meloche_df: pd.DataFrame,
        cluster_map: np.ndarray,
        transform,
        snap_features: Optional[pd.DataFrame] = None,
        pixel_index: Optional[dict] = None,
        n_lateral_neighbors: int = 6,
        gaume_aspect_cap: float = GAUME_ASPECT_CAP,
        aspect_deg: Optional[float] = None) -> float:
    """Estimate cross-slope release width (m) using Gaume et al. (2015) / θ.

    Falls back to A_ca whenever θ, τp or `aspect_deg` is unavailable.

    `aspect_deg` is the downslope bearing at the trigger and is what defines
    "cross-slope": the lateral sample is restricted to a ±THETA_CROSS_SECTOR_DEG
    sector about the cross axis (aspect ± 90°), on either flank. Without it
    there is no cross-slope direction to sample along, so the θ widening is
    skipped rather than taken along the raster axes.

    This sets the lateral *distance cap*; the cross-slope Λ continuity test
    lives in propagate_crack via directional_lambda(). Both are mode III
    concerns and should be revisited together when Λ_III lands.

    Both halves of the width ratio are measured here, by `_sector_theta`, from
    the same τp field in the same separation band: θ_along about the fall line,
    θ_cross about the cross axis. Taking the numerator from the meloche CSV
    instead — as this did before 2026-10-09 — pairs an *isotropic* k-NN mean at
    a ~3.3 m lag with a directional mean at 5-50 m, and since
    θ(d) ≈ 3.87 + 72/d Pa/m on this slope (methods §5) the lag mismatch alone
    inflates the ratio several-fold and drives it into `gaume_aspect_cap`.
    The CSV θ is still used as a fallback numerator when no cluster lies in the
    along-slope sector.
    """
    if meloche_df.empty or 'theta' not in meloche_df.columns:
        return A_ca

    if aspect_deg is None or np.isnan(aspect_deg):
        return A_ca

    try:
        theta_down = float(meloche_df.loc[trigger_cluster_id, 'theta'])
    except (KeyError, TypeError):
        return A_ca

    if np.isnan(theta_down) or theta_down < config.THETA_MIN:
        return A_ca

    tau_frame, tau_col = _tau_p_lookup(meloche_df, snap_features)
    if tau_frame is None or trigger_cluster_id not in tau_frame.index:
        return A_ca
    tau_trigger = float(pd.to_numeric(
        first_scalar(tau_frame.loc[trigger_cluster_id, tau_col]), errors='coerce'))
    if np.isnan(tau_trigger):
        return A_ca

    if pixel_index is None:
        pixel_index = cluster_pixel_index(cluster_map)
    centroids = cluster_centroids_px(pixel_index)
    if trigger_cluster_id not in centroids:
        return A_ca
    t_row, t_col = centroids[trigger_cluster_id]
    t_x, t_y = pixel_to_utm(int(round(t_row)), int(round(t_col)), transform)

    # Separation of every other cluster from the trigger, in the map frame,
    # restricted once to the usable band and reused for both axes.
    asp_rad  = np.radians(float(aspect_deg))
    fall_x, fall_y   = np.sin(asp_rad), np.cos(asp_rad)
    cross_x, cross_y = -fall_y, fall_x

    offsets = []
    for cid, (c_row, c_col) in centroids.items():
        if cid == trigger_cluster_id:
            continue
        c_x, c_y = pixel_to_utm(int(round(c_row)), int(round(c_col)), transform)
        dx, dy   = c_x - t_x, c_y - t_y
        dist_m   = float(np.hypot(dx, dy))
        if config.THETA_CROSS_MIN_SEP_M <= dist_m <= config.THETA_CROSS_MAX_SEP_M:
            offsets.append((cid, dx, dy, dist_m))
    if not offsets:
        return A_ca

    theta_cross = _sector_theta(offsets, cross_x, cross_y, tau_trigger,
                                tau_frame, tau_col, n_lateral_neighbors)
    if np.isnan(theta_cross):
        return A_ca

    # Numerator measured the same way, about the fall line. Falling back to the
    # CSV θ reintroduces the estimator/lag mismatch, so it is a last resort.
    theta_along = _sector_theta(offsets, fall_x, fall_y, tau_trigger,
                                tau_frame, tau_col, n_lateral_neighbors)
    if np.isnan(theta_along) or theta_along < config.THETA_MIN:
        theta_along = theta_down

    if theta_cross > config.THETA_MIN:
        width_factor = min(theta_along / theta_cross, gaume_aspect_cap)
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
        max_clusters: int = config.MAX_BFS_CLUSTERS,
        pixel_index: Optional[dict] = None,
        use_propagation: bool = True):
    """
    Build a 2D release polygon.

    Primary method (use_propagation=True):
        BFS cluster flood-fill with per-direction arrest criteria.
    Fallback (use_propagation=False, or propagation returns None):
        Oriented rectangle from A_ca + stauchwall + Gaume cross-slope width.

    Returns shapely Polygon (UTM) or None.
    """
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    import rasterio.features

    px_m = abs(transform.a) if transform is not None else 1.0
    slope_grid, aspect_grid = compute_slope_aspect(dem, px_m)
    if pixel_index is None:
        pixel_index = cluster_pixel_index(cluster_map)

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
        polygon, _reached = propagate_crack(
            trigger_cluster_id = trigger_cluster_id,
            meloche_df         = meloche_df,
            cluster_map        = cluster_map,
            dem                = dem,
            slope_grid         = slope_grid,
            transform          = transform,
            start_zone_mask    = start_zone_mask,
            stauchwall_deg     = stauchwall_deg,
            size_factor        = size_factor,
            A_ca               = A_ca,
            snap_features      = snap_features,
            max_clusters       = max_clusters,
            pixel_index        = pixel_index,
            aspect_grid        = aspect_grid,
        )
        if polygon is not None:
            polygon = _clip_polygon(polygon)
            if polygon is not None:
                return polygon
        print(f"  propagate_crack returned None for cluster {trigger_cluster_id}"
              f" — falling back to rectangle")

    if trigger_cluster_id not in pixel_index:
        return None
    rows, cols  = pixel_index[trigger_cluster_id]
    t_row, t_col = int(rows.mean()), int(cols.mean())
    t_x, t_y    = pixel_to_utm(t_row, t_col, transform)
    t_aspect    = float(aspect_grid[t_row, t_col])

    up_x, up_y  = project_along_aspect(t_x, t_y, A_ca * size_factor,
                                        t_aspect, upslope=True)
    sw_row, sw_col = find_stauchwall(t_row, t_col, slope_grid, aspect_grid,
                                     threshold_deg=stauchwall_deg)
    sw_x, sw_y  = pixel_to_utm(sw_row, sw_col, transform)
    half_width  = estimate_cross_slope_width(
        trigger_cluster_id, A_ca, meloche_df, cluster_map, transform,
        snap_features=snap_features, pixel_index=pixel_index,
        aspect_deg=t_aspect,
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
        size_factor: float = 1.0,
        A_ca: Optional[float] = None,
        snap_features: 'Optional[pd.DataFrame]' = None,
        k_neighbours: int = config.BFS_K_NEIGHBOURS,
        max_clusters: int = config.MAX_BFS_CLUSTERS,
        pixel_index: Optional[dict] = None,
        aspect_grid: Optional[np.ndarray] = None,
        wave_callback=None):
    """
    Identify release zone as the connected region of clusters reachable
    from the trigger where crack arrest criteria are not met.

    Starting from the trigger cluster, expands outward to neighbouring
    clusters in order of proximity. A cluster is included if it passes
    distance caps (upslope A_ca, downslope stauchwall, lateral Gaume width),
    the stauchwall slope angle, an absolute tau_g floor, slab thickness, and
    Lambda/thickness continuity checks.

    Note: the Meloche per-direction arrest criterion is gated behind
    config.USE_MELOCHE_ARREST and is off by default — the continuity
    heuristics plus config.TAU_G_ABS_FLOOR do the arresting. See methods §5.

    Returns (polygon, reached_cluster_ids).
    """
    from shapely.ops import unary_union
    from shapely.geometry import shape as shapely_shape
    import rasterio.features

    if meloche_df.empty or trigger_cluster_id not in meloche_df.index:
        return None, set()
    if 'Pi1_elastic' not in meloche_df.columns:
        return None, set()

    pi1_trigger = float(first_scalar(meloche_df.loc[trigger_cluster_id, 'Pi1_elastic']))
    if np.isnan(pi1_trigger) or pi1_trigger <= 0:
        return None, set()

    if pixel_index is None:
        pixel_index = cluster_pixel_index(cluster_map)
    if trigger_cluster_id not in pixel_index:
        return None, set()

    px_m = abs(transform.a) if transform is not None else 1.0
    if aspect_grid is None:
        _, aspect_grid = compute_slope_aspect(dem, px_m)

    print(f"    Pi1_trigger={pi1_trigger:.3f}  "
          f"tau_g_floor={config.TAU_G_ABS_FLOOR:.0f}Pa  "
          f"size_factor={size_factor:.2f}")

    slope_cache: dict[int, float] = {}

    def _mean_slope(cid):
        if cid not in slope_cache:
            px = pixel_index.get(cid)
            slope_cache[cid] = (float(slope_grid[px[0], px[1]].mean())
                                if px is not None else 0.0)
        return slope_cache[cid]

    t_rows, t_cols = pixel_index[trigger_cluster_id]
    t_row, t_col = int(t_rows.mean()), int(t_cols.mean())
    t_x, t_y = pixel_to_utm(t_row, t_col, transform)
    t_aspect = float(aspect_grid[t_row, t_col])

    d_up = (max(A_ca, 20.0) if A_ca else 50.0) * size_factor

    sw_row, sw_col = find_stauchwall(
        t_row, t_col, slope_grid, aspect_grid, threshold_deg=stauchwall_deg)
    sw_x, sw_y = pixel_to_utm(sw_row, sw_col, transform)
    d_down = max(float(np.hypot(sw_x - t_x, sw_y - t_y)), 20.0) * size_factor

    print(f"    distance caps: up={d_up:.0f}m  down={d_down:.0f}m")

    from sklearn.neighbors import NearestNeighbors
    cents_px   = cluster_centroids_px(pixel_index)
    cids       = np.array(sorted(cents_px))
    pxs        = np.array([cents_px[c] for c in cids])
    k_actual   = min(k_neighbours + 1, len(cids))
    nbrs       = NearestNeighbors(n_neighbors=k_actual).fit(pxs)
    _, idxs    = nbrs.kneighbors(pxs)
    neighbours = {int(cids[i]): [int(cids[j]) for j in idxs[i][1:]]
                  for i in range(len(cids))}

    centroids = {int(c): pixel_to_utm(int(round(pxs[i][0])),
                                      int(round(pxs[i][1])), transform)
                 for i, c in enumerate(cids)}

    A_ca_base = A_ca if A_ca else 50.0
    d_lat_gaume = estimate_cross_slope_width(
        trigger_cluster_id, A_ca_base,
        meloche_df, cluster_map, transform,
        snap_features=snap_features, pixel_index=pixel_index,
        aspect_deg=t_aspect) * size_factor

    # Mode III crack-speed cap. Mode III cannot exceed c_s where upslope runs
    # supershear at ~1.6 c_s; the slower crack builds slab tension faster per
    # unit advance, so first slab fracture — and hence arrest — comes sooner
    # cross-slope. Composed as a min() with the Gaume width because arrest
    # happens at whichever constraint binds first.
    d_lat = d_lat_gaume
    mode3_ratio = mode3_speed_cap_ratio(trigger_cluster_id, meloche_df)
    if mode3_ratio is not None:
        d_lat_speed = A_ca_base * mode3_ratio * size_factor
        d_lat = min(d_lat, d_lat_speed)
        print(f"    d_lat={max(d_lat, 15.0):.0f}m  "
              f"(gaume={d_lat_gaume:.0f}m, mode3_speed_cap={d_lat_speed:.0f}m "
              f"@ ratio {mode3_ratio:.3f})")
    else:
        print(f"    d_lat={max(d_lat, 15.0):.0f}m  (gaume only, "
              f"mode III speed cap off)")
    d_lat = max(d_lat, 15.0)

    LAMBDA_DROP_FACTOR     = config.LAMBDA_DROP_FACTOR
    LAMBDA_RISE_FACTOR     = config.LAMBDA_RISE_FACTOR
    THICKNESS_DROP_FACTOR  = config.THICKNESS_DROP_FACTOR
    THICKNESS_RISE_FACTOR  = config.THICKNESS_RISE_FACTOR
    MIN_PROPAGATION_SLAB   = config.MIN_PROPAGATION_SLAB
    MIN_PROPAGATION_LAMBDA = config.MIN_PROPAGATION_LAMBDA
    TAU_G_ABS_FLOOR        = config.TAU_G_ABS_FLOOR

    def _get_props(cid):
        props = {'cid': cid}
        if not meloche_df.empty and cid in meloche_df.index:
            for c in ['Lambda', 'Lambda_cross', 'tau_g', 'rc_wl', 'L_t']:
                if c in meloche_df.columns:
                    props[c] = float(first_scalar(meloche_df.loc[cid, c]))
        if snap_features is not None and cid in snap_features.index:
            for c in ['slab_thickness', 'slab_density',
                      'wl_shear_strength', 'sigma_t']:
                if c in snap_features.columns:
                    try:
                        props[c] = float(first_scalar(snap_features.loc[cid, c]))
                    except (ValueError, TypeError):
                        pass
        return props

    trigger_props = _get_props(trigger_cluster_id)

    def _qualifies(cid, current_props):
        if start_zone_mask is not None:
            px = pixel_index.get(cid)
            if px is None:
                return False, 'no_data'
            r, c = int(px[0].mean()), int(px[1].mean())
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
        # Anything not within ±45° of the fall line is treated as cross-slope,
        # i.e. mode III / flank propagation.
        is_cross     = not is_upslope and not is_downslope

        if is_upslope   and dist > d_up:   return False, 'upslope_distance_cap'
        if is_downslope and dist > d_down: return False, 'downslope_distance_cap'
        if is_cross and dist > d_lat:
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

        lam_nbr_abs = directional_lambda(nbr_props, is_cross)
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
                        Lam   = directional_lambda(nbr_props, is_cross)
                        tg    = tau_g_nbr
                        sig_t = nbr_props.get('sigma_t', np.nan)
                        L_t   = nbr_props.get('L_t',     np.nan)
                        # np.isfinite, not isnan: tensile_length returns inf
                        # when k_f <= 0, and arrest_indices.evaluate() gates
                        # Eq. 20 on isfinite(L_t) for the same reason.
                        if (all(np.isfinite(v) for v in (Lam, tg, sig_t, L_t))
                                and tg > 0.0):
                            # Eq. 20 with the directional theta, from the one
                            # implementation in arrest_indices.
                            A_ca_d, _ = arrest_length(
                                tg, theta_dir, Lam, sig_t,
                                delta=MELOCHE_DELTA, L_t=L_t)
                            if A_ca_d < dist_cn:
                                return False, 'meloche_arrest'
            return True, 'propagated'

        # Λ continuity, evaluated with the elastic length and thresholds for
        # this propagation direction. Cross-slope currently resolves to the
        # same Λ and the same thresholds as along-slope, so the split is a
        # no-op until Λ_III lands — see config.USE_MODE3_LAMBDA.
        lam_current = directional_lambda(current_props, is_cross)
        lam_nbr     = directional_lambda(nbr_props,     is_cross)
        lam_drop, lam_rise = (
            (config.LAMBDA_CROSS_DROP_FACTOR, config.LAMBDA_CROSS_RISE_FACTOR)
            if is_cross else
            (LAMBDA_DROP_FACTOR, LAMBDA_RISE_FACTOR))
        if (not np.isnan(lam_current) and not np.isnan(lam_nbr)
                and lam_current > 0):
            signed_change = (lam_nbr - lam_current) / lam_current
            if signed_change < 0:
                if -signed_change > lam_drop * size_factor:
                    return False, ('lambda_cross_discontinuity_drop' if is_cross
                                   else 'lambda_discontinuity_drop')
            else:
                if signed_change > lam_rise * size_factor:
                    return False, ('lambda_cross_discontinuity_rise' if is_cross
                                   else 'lambda_discontinuity_rise')

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

    # The trigger cluster is included unconditionally: it is the seed, and it
    # already passed the generate_scenarios filter chain.
    reached        = {trigger_cluster_id}
    cluster_props  = {trigger_cluster_id: trigger_props}
    visited        = {trigger_cluster_id}
    current_ring   = [trigger_cluster_id]
    ring_idx       = 0

    if wave_callback is not None:
        wave_callback(0, [trigger_cluster_id], [])

    while current_ring and len(reached) < max_clusters:
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
                    reached.add(nbr)
                    cluster_props[nbr] = _get_props(nbr)
                    next_ring.append(nbr)
                else:
                    ring_arrested.append({'cid': nbr, 'reason': reason})

        if not next_ring and not ring_arrested:
            break

        if wave_callback is not None:
            wave_callback(ring_idx, next_ring, ring_arrested)

        current_ring = next_ring

    # Distinguish "arrest criteria stopped the crack" from "the safety cap did".
    cap_bound = len(reached) >= max_clusters and bool(current_ring)
    if cap_bound:
        print(f"  WARNING: BFS hit the max_clusters cap ({max_clusters}) with "
              f"{len(current_ring)} clusters still propagating — this polygon "
              f"is bounded by the cap, not by arrest criteria. "
              f"Raise --max-clusters to let the physics terminate it.")

    dir_counts = {'upslope': 0, 'downslope': 0, 'lateral': 0, 'other': 0}
    for cid in reached:
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
    print(f"  Connected region: {len(reached)} clusters "
          f"(up={dir_counts['upslope']} "
          f"down={dir_counts['downslope']} "
          f"lat={dir_counts['lateral']})"
          f"{'  [CAP-BOUND]' if cap_bound else ''}")

    if len(reached) < 2:
        return None, reached

    mask   = np.isin(cluster_map, list(reached)).astype(np.uint8)
    shapes = list(rasterio.features.shapes(
        mask, mask=mask.astype(bool), transform=transform))
    if not shapes:
        return None, reached

    polys   = [shapely_shape(s) for s, v in shapes if v == 1]
    polygon = unary_union(polys)
    if polygon.geom_type == 'MultiPolygon':
        polygon = max(polygon.geoms, key=lambda p: p.area)
    if polygon.is_empty or polygon.area < MIN_POLYGON_AREA:
        return None, reached

    polygon = fill_polygon_holes(polygon)

    return polygon, reached


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

def load_observed_polygon(path, dst_epsg: int = config.RASTER_EPSG):
    """Load a GeoJSON polygon file and reproject to the raster CRS."""
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


def load_observed_polygons(primary_path, dst_epsg: int = config.RASTER_EPSG) -> list:
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
    px_m      = abs(transform.a) if transform is not None else 1.0
    hillshade = LightSource(azdeg=315, altdeg=45).hillshade(fill_dem, dx=px_m, dy=px_m)
    nrows, ncols = dem.shape
    extent = [transform.c,
              transform.c + ncols * transform.a,
              transform.f + nrows * transform.e,
              transform.f]
    # aspect='equal': the extent is projected metres (EPSG:6342), so 'auto'
    # stretches east against north to fill the axes box and every bearing
    # read off the figure is wrong. On the Jan 18 domain (455 x 546 m in a
    # square figure) that was a 1.21x anisotropy, a -4.6 deg apparent rotation.
    ax.imshow(hillshade, cmap='gray', extent=extent,
              alpha=0.6, aspect='equal', origin='upper')

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
