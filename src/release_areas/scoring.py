"""
scoring — Jan 18 calibration scores, ported from avachain (ops-2026-27).

Sources: avachain ``src/avachain/event_scoring.py`` (event scores) and
``scripts/inside_release_sweep.py`` (forced shots). Kept numerically the same
so a number here can be set beside one from avachain. Runs on this repo's
``data/little_prof`` (Ron, 2026-11-18: score here, SNOWPACK features from
avachain).

Two scores:

- **Five-trigger event score.** The model's own triggers, each weighted by
  inverse Sk38 (avachain ``_trigger_weights``). P(pixel in release | release)
  is the weight sum of the polygons covering the pixel. E[IoU], Brier skill
  and AUC are taken against the observed mask over start zone OR observed.
  avachain also sums over size factors and depth percentiles, but with the
  distance caps off they share one polygon per trigger, so the trigger weights
  alone give the same numbers.
- **Forced shots.** Every cluster mostly inside the observed release
  (positives) and every cluster under an observed ski non-trigger point
  (negatives) is used as a trigger, with no initiation gate and no rectangle
  fallback.

Jan 18 is a calibration target. These are calibration numbers, never
validation.
"""
from __future__ import annotations

import numpy as np


def iou(a, b) -> float:
    if a is None or a.is_empty:
        return 0.0
    u = a.union(b).area
    return a.intersection(b).area / u if u > 0 else 0.0


def trigger_weights(sk38_values) -> list:
    """Inverse-Sk38 weights, normalised (avachain ``_trigger_weights``)."""
    inv = [1.0 / max(float(s), 0.01) for s in sk38_values]
    tot = sum(inv)
    return [v / tot for v in inv]


def rasterize(polys, shape, transform) -> list:
    """One boolean mask per polygon (None or empty -> all False)."""
    from rasterio.features import rasterize as _r
    out = []
    for p in polys:
        if p is None or p.is_empty:
            out.append(np.zeros(shape, dtype=bool))
        else:
            out.append(_r([(p, 1)], out_shape=shape, transform=transform,
                          dtype=np.uint8).astype(bool))
    return out


def auc(p: np.ndarray, o: np.ndarray) -> float:
    """Area under the ROC curve (Mann-Whitney, ties averaged); NaN if o is constant."""
    from scipy.stats import rankdata
    o = o.astype(bool)
    n1, n0 = int(o.sum()), int((~o).sum())
    if n1 == 0 or n0 == 0:
        return float('nan')
    r = rankdata(p)
    return float((r[o].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def event_score(prob: np.ndarray, start_zone: np.ndarray, observed: np.ndarray,
                pixel_area: float = 1.0) -> dict:
    """Event-day scores over the domain start_zone | observed
    (avachain ``event_scoring.score``)."""
    sz, obs = start_zone.astype(bool), observed.astype(bool)
    dom = sz | obs
    p, o = prob[dom], obs[dom].astype(float)
    clim = o.mean()
    bs, bs_ref = float(np.mean((p - o) ** 2)), float(np.mean((clim - o) ** 2))
    return {
        'domain_m2': float(dom.sum() * pixel_area),
        'expected_area_m2': float(p.sum() * pixel_area),
        'any_prob_area_m2': float((p > 0).sum() * pixel_area),
        'observed_m2': float(o.sum() * pixel_area),
        'brier': bs,
        'brier_skill': 1.0 - bs / bs_ref if bs_ref > 0 else float('nan'),
        'auc': auc(p, o),
        'expected_hit_m2': float((p * o).sum() * pixel_area),
        'expected_false_m2': float((p * (1 - o)).sum() * pixel_area),
        'observed_with_any_prob': float(np.mean(p[o > 0] > 0)),
    }


def weighted_event_score(polys, weights, observed_poly, start_zone, shape,
                         transform) -> dict:
    """E[IoU] and the event scores of a weighted polygon set."""
    masks = rasterize(polys, shape, transform)
    prob = np.zeros(shape, dtype=np.float64)
    for m, w in zip(masks, weights):
        prob[m] += w
    obs = rasterize([observed_poly], shape, transform)[0]
    out = event_score(prob, start_zone, obs, abs(transform.a * transform.e))
    out['expected_iou'] = float(sum(w * iou(p, observed_poly)
                                    for p, w in zip(polys, weights)))
    return out


def inside_clusters(cluster_map, obs_mask, start_zone_mask, observed_poly,
                    transform, min_frac: float = 0.5) -> list:
    """Clusters with >= min_frac of their pixels inside the observed polygon
    and inside the start zone, with the centroid's distance to the observed
    boundary (m). avachain ``inside_release_sweep.inside_clusters``."""
    from shapely.geometry import Point
    from release_areas.release_geometry import pixel_to_utm
    rows = []
    for cid in np.unique(cluster_map[obs_mask & (cluster_map > 0)]):
        px = cluster_map == cid
        n = px.sum()
        f_obs, f_sz = (px & obs_mask).sum() / n, (px & start_zone_mask).sum() / n
        if f_obs < min_frac or f_sz < min_frac:
            continue
        r, c = np.argwhere(px).mean(axis=0)
        x, y = pixel_to_utm(int(r), int(c), transform)
        rows.append({'cid': int(cid), 'frac_inside': round(float(f_obs), 3),
                     'edge_dist_m': round(float(observed_poly.exterior.distance(Point(x, y))), 1)})
    return rows


def clusters_at_points(points_lat_lon, cluster_map, transform, epsg=None) -> list:
    """Cluster id under each (lat, lon) WGS84 point, in order, no duplicates;
    points on no cluster are skipped (avachain ``_clusters_at_points``)."""
    from pyproj import Transformer
    from rasterio.transform import rowcol
    from release_areas import config
    epsg = config.RASTER_EPSG if epsg is None else epsg
    to_xy = Transformer.from_crs(4326, epsg, always_xy=True)
    nrow, ncol = cluster_map.shape
    cids = []
    for lat, lon in points_lat_lon:
        x, y = to_xy.transform(lon, lat)
        cid = 0
        if np.isfinite(x) and np.isfinite(y):
            r, c = rowcol(transform, x, y)
            if 0 <= r < nrow and 0 <= c < ncol:
                cid = int(cluster_map[r, c])
        if cid > 0 and cid not in cids:
            cids.append(cid)
    return cids


def consensus(masks, obs_mask):
    """(IoU of the >= half-coverage pixels vs observed, mean pairwise IoU)."""
    if not masks:
        return None, None
    stack = np.stack(masks).reshape(len(masks), -1)
    cons = stack.mean(axis=0).reshape(obs_mask.shape) >= 0.5
    cons_iou = float((cons & obs_mask).sum() / max((cons | obs_mask).sum(), 1))
    if len(masks) < 2:
        return cons_iou, None
    s = stack.astype(np.float32)
    inter = s @ s.T
    area = np.diag(inter)
    union = area[:, None] + area[None, :] - inter
    iu = np.triu_indices(len(masks), k=1)
    return cons_iou, float(np.mean(inter[iu] / np.maximum(union[iu], 1)))
