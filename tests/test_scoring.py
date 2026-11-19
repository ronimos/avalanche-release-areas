"""Unit tests for release_areas.scoring (port of avachain's Jan 18 scorers)."""
import numpy as np
import pytest
from rasterio.transform import from_origin
from shapely.geometry import box

from release_areas import scoring

TF = from_origin(0.0, 10.0, 1.0, 1.0)     # 10 x 10 grid of 1 m pixels
SHAPE = (10, 10)


def test_trigger_weights_inverse_sk38_normalised():
    w = scoring.trigger_weights([0.5, 1.0])
    assert w == pytest.approx([2 / 3, 1 / 3])
    assert sum(scoring.trigger_weights([0.79, 0.81, 0.81, 0.81, 0.82])) == pytest.approx(1.0)


def test_iou():
    assert scoring.iou(box(0, 0, 2, 2), box(1, 0, 3, 2)) == pytest.approx(1 / 3)
    assert scoring.iou(None, box(0, 0, 1, 1)) == 0.0


def test_auc_perfect_inverted_and_constant():
    o = np.array([0, 0, 1, 1], bool)
    assert scoring.auc(np.array([0., 0., 1., 1.]), o) == 1.0
    assert scoring.auc(np.array([1., 1., 0., 0.]), o) == 0.0
    assert scoring.auc(np.array([.5, .5, .5, .5]), o) == 0.5
    assert np.isnan(scoring.auc(np.zeros(4), np.zeros(4, bool)))


def test_weighted_event_score_single_exact_polygon():
    """One polygon equal to the observed: E[IoU] 1, AUC 1, Brier 0."""
    obs = box(2, 2, 6, 6)
    sz = np.ones(SHAPE, bool)
    s = scoring.weighted_event_score([obs], [1.0], obs, sz, SHAPE, TF)
    assert s['expected_iou'] == pytest.approx(1.0)
    assert s['expected_area_m2'] == pytest.approx(16.0)
    assert s['observed_m2'] == pytest.approx(16.0)
    assert s['auc'] == 1.0
    assert s['brier'] == 0.0
    assert s['brier_skill'] == pytest.approx(1.0)


def test_weighted_event_score_weights_are_probabilities():
    obs = box(0, 0, 4, 4)
    sz = np.ones(SHAPE, bool)
    s = scoring.weighted_event_score([obs, box(6, 6, 8, 8)], [0.75, 0.25], obs, sz, SHAPE, TF)
    assert s['expected_iou'] == pytest.approx(0.75)
    assert s['expected_area_m2'] == pytest.approx(0.75 * 16 + 0.25 * 4)
    assert s['expected_hit_m2'] == pytest.approx(0.75 * 16)
    assert s['expected_false_m2'] == pytest.approx(0.25 * 4)


def test_domain_includes_observed_outside_start_zone():
    """Release outside the start zone still counts as missed (avachain)."""
    obs = box(0, 0, 4, 4)
    sz = np.zeros(SHAPE, bool)
    sz[:, 5:] = True
    s = scoring.weighted_event_score([box(6, 6, 8, 8)], [1.0], obs, sz, SHAPE, TF)
    assert s['observed_m2'] == pytest.approx(16.0)
    assert s['observed_with_any_prob'] == 0.0


def test_consensus_identical_and_disjoint():
    obs = np.zeros(SHAPE, bool)
    obs[:4, :4] = True
    c, p = scoring.consensus([obs, obs.copy()], obs)
    assert (c, p) == (1.0, 1.0)
    other = np.zeros(SHAPE, bool)
    other[6:, 6:] = True
    c, p = scoring.consensus([obs, other], obs)
    assert p == 0.0
    assert c < 1.0      # >= half coverage takes both squares


def test_inside_clusters_fraction_and_edge_distance():
    cmap = np.zeros(SHAPE)
    cmap[0:4, 0:4] = 1        # rows 0-3 = northing 10..6, fully inside obs
    cmap[0:4, 4:8] = 2        # half inside
    cmap[6:10, 6:10] = 3      # outside
    obs_poly = box(0, 6, 6, 10)
    obs = scoring.rasterize([obs_poly], SHAPE, TF)[0]
    sz = np.ones(SHAPE, bool)
    rows = {r['cid']: r for r in scoring.inside_clusters(cmap, obs, sz, obs_poly, TF)}
    assert set(rows) == {1, 2}
    assert rows[1]['frac_inside'] == 1.0
    assert rows[2]['frac_inside'] == 0.5


def test_clusters_at_points_dedupes_and_skips_off_grid():
    from pyproj import Transformer
    from rasterio.transform import from_origin as fo
    from release_areas import config
    x0, y0 = 420000.0, 4400000.0
    tf = fo(x0, y0, 1.0, 1.0)
    cmap = np.zeros(SHAPE)
    cmap[0, 0], cmap[5, 5] = 7, 9
    to_ll = Transformer.from_crs(config.RASTER_EPSG, 4326, always_xy=True)
    def ll(x, y):
        lon, lat = to_ll.transform(x, y)
        return (lat, lon)
    pts = [ll(x0 + 0.5, y0 - 0.5), ll(x0 + 0.6, y0 - 0.4),      # cluster 7 twice
           ll(x0 + 5.5, y0 - 5.5),                              # cluster 9
           ll(x0 + 2.5, y0 - 2.5),                              # no cluster
           ll(x0 - 50, y0 + 50)]                                # off grid
    assert scoring.clusters_at_points(pts, cmap, tf) == [7, 9]
