"""Tests for anisotropic BFS ring growth in propagate_crack.

The hypothesis under test: if one ring of the flood-fill stands for one
increment of crack-front advance in *time*, a front running supershear
along-slope (~1.6 c_s) and capped at c_s cross-slope should already be ~1.6x
further along-slope when it arrives, so the anisotropy belongs in the growth
rather than only in a terminal distance cap (config.USE_ANISOTROPIC_RINGS).

Two invariants matter:

  1. Default OFF must be bit-identical to the isotropic flood-fill — the
     reference Jan 18 numbers depend on it.
  2. The cross-slope budget must be driven by the ALONG-SLOPE clock. Driving it
     from the global ring counter is self-defeating: flank growth keeps the
     flood-fill alive, which raises the flanks' own budget.
"""
import numpy as np
import pandas as pd
import pytest
from affine import Affine

from release_areas import config
from release_areas.release_geometry import (
    _hop_length_m,
    cluster_pixel_index,
    compute_slope_aspect,
    propagate_crack,
)

T = Affine(1.0, 0.0, 400000.0, 0.0, -1.0, 4000000.0)
BLOCK = 4          # cluster side in pixels -> 4 m hop, as on Little Professor
N = 161            # grid side in pixels; trigger at the centre


@pytest.fixture
def scene():
    """Uniform 35 deg south-facing plane tiled with 4x4 m clusters.

    Due south means the fall line runs along +row and the cross-slope axis
    along +col, so "along-slope" and "cross-slope" are the array axes and the
    expected geometry is easy to state. Every cluster carries identical
    properties, so no continuity or floor gate ever fires and the only things
    that can stop the flood-fill are the distance caps and the ring budget.
    """
    r, c = np.mgrid[0:N, 0:N]
    cluster_map = ((r // BLOCK) * (N // BLOCK + 1) + (c // BLOCK) + 1).astype(float)
    dem = -np.tan(np.radians(35.0)) * (-r.astype(float))   # falls toward +row
    slope_grid, aspect_grid = compute_slope_aspect(dem, 1.0)

    pixel_index = cluster_pixel_index(cluster_map)
    cids = sorted(pixel_index)
    trigger = int(cluster_map[N // 2, N // 2])

    meloche_df = pd.DataFrame({
        'Pi1_elastic': 50.0, 'Lambda': 2.0, 'Lambda_cross': 2.0,
        'tau_g': 800.0, 'L_t': 30.0, 'A_ca_brittle': 60.0,
    }, index=cids)
    snap_features = pd.DataFrame({
        'slab_thickness': 0.8, 'slab_density': 330.0,
        'wl_shear_strength': 2400.0, 'sigma_t': 4000.0,
    }, index=cids)

    return dict(trigger=trigger, cluster_map=cluster_map, dem=dem,
                slope_grid=slope_grid, aspect_grid=aspect_grid,
                meloche_df=meloche_df, snap_features=snap_features,
                pixel_index=pixel_index)


def _run(scene, **over):
    # A_ca is deliberately larger than the domain so neither the upslope cap
    # nor the lateral Gaume width (which falls back to A_ca without a theta
    # column) binds: the isotropic run is then bounded by the grid in every
    # direction and any anisotropy in the result is the ring budget's doing.
    kw = dict(trigger_cluster_id=scene['trigger'], meloche_df=scene['meloche_df'],
              cluster_map=scene['cluster_map'], dem=scene['dem'],
              slope_grid=scene['slope_grid'], transform=T,
              A_ca=200.0, snap_features=scene['snap_features'],
              pixel_index=scene['pixel_index'], aspect_grid=scene['aspect_grid'],
              max_clusters=20000)
    kw.update(over)
    return propagate_crack(**kw)


def _extent(scene, reached):
    """Max along-slope reach and flank width at the crown, both in metres.

    The flank width is measured only within a band of the trigger's own row,
    because the ±45° sector test classifies a corner cluster at
    (Δalong, Δcross) = (80, 80) as *downslope* — so a plain max |Δcol| over the
    whole region is set by corner leakage, not by the cross-slope budget.
    """
    px = scene['pixel_index']
    t_rows, t_cols = px[scene['trigger']]
    t_r, t_c = t_rows.mean(), t_cols.mean()
    along = cross = 0.0
    for cid in reached:
        rows, cols = px[cid]
        d_along = abs(rows.mean() - t_r)
        along = max(along, d_along)
        if d_along <= 2 * BLOCK:
            cross = max(cross, abs(cols.mean() - t_c))
    return along, cross


@pytest.fixture
def rings_on(monkeypatch):
    monkeypatch.setattr(config, 'USE_ANISOTROPIC_RINGS', True)
    monkeypatch.setattr(config, 'USE_MODE3_SPEED_CAP', False)


@pytest.fixture
def rings_off(monkeypatch):
    monkeypatch.setattr(config, 'USE_ANISOTROPIC_RINGS', False)
    monkeypatch.setattr(config, 'USE_MODE3_SPEED_CAP', False)


class TestDefaults:
    def test_rings_are_off_by_default(self):
        """The committed Jan 18 reference numbers are isotropic-ring numbers."""
        assert config.USE_ANISOTROPIC_RINGS is False

    def test_ring_step_is_derived_not_hardcoded(self):
        assert config.RING_STEP_M is None

    def test_speed_ratio_is_broberg_not_a_free_parameter(self):
        """1.6 / 1.0 comes from Broberg (1989) via arrest_indices, and the ring
        budget must reuse it rather than introduce a second constant."""
        assert config.MODE2_SPEED_RATIO == 1.6
        assert config.MODE3_SPEED_RATIO == 1.0
        assert config.MODE3_SPEED_RATIO / config.MODE2_SPEED_RATIO == 0.625


class TestHopLength:
    def test_median_edge_separation(self):
        centroids = {1: (0.0, 0.0), 2: (0.0, 3.0), 3: (0.0, 10.0)}
        neighbours = {1: [2, 3], 2: [1], 3: [1]}
        # edges: 3, 10, 3, 10 -> median 6.5
        assert _hop_length_m(neighbours, centroids) == pytest.approx(6.5)

    def test_missing_centroids_are_skipped(self):
        centroids = {1: (0.0, 0.0), 2: (0.0, 4.0)}
        neighbours = {1: [2, 99], 2: [1], 99: [1]}
        assert _hop_length_m(neighbours, centroids) == pytest.approx(4.0)

    def test_empty_graph_gives_zero(self):
        assert _hop_length_m({}, {}) == 0.0

    def test_matches_the_cluster_tiling_of_the_scene(self, scene):
        """4x4 pixel clusters on a 1 m grid -> a 4 m hop."""
        from release_areas.release_geometry import cluster_centroids_px
        from sklearn.neighbors import NearestNeighbors
        cents = cluster_centroids_px(scene['pixel_index'])
        cids = np.array(sorted(cents))
        pxs = np.array([cents[c] for c in cids])
        nb = NearestNeighbors(n_neighbors=9).fit(pxs)
        _, idx = nb.kneighbors(pxs)
        nbrs = {int(cids[i]): [int(cids[j]) for j in idx[i][1:]] for i in range(len(cids))}
        cent_m = {int(c): (cents[c][1], -cents[c][0]) for c in cids}
        assert _hop_length_m(nbrs, cent_m) == pytest.approx(BLOCK, abs=1.5)


class TestAnisotropicGrowth:
    def test_rings_shorten_the_flanks_but_not_the_fall_line(self, scene, rings_on, rings_off):
        """The whole point: cross-slope reach shrinks, along-slope reach does not."""
        # rings_off fixture applied last wins, so set explicitly per run
        config.USE_ANISOTROPIC_RINGS = False
        _, iso = _run(scene)
        config.USE_ANISOTROPIC_RINGS = True
        _, ani = _run(scene)

        iso_along, iso_cross = _extent(scene, iso)
        ani_along, ani_cross = _extent(scene, ani)

        assert ani_along == pytest.approx(iso_along, abs=BLOCK)
        assert ani_cross < iso_cross - 2 * BLOCK
        assert len(ani) < len(iso)

    def test_released_aspect_ratio_moves_toward_the_speed_ratio(self, scene, rings_on):
        """With the flanks budget-limited the region should be elongated along
        the fall line by roughly MODE2/MODE3 = 1.6, not isotropic."""
        config.USE_ANISOTROPIC_RINGS = True
        _, ani = _run(scene)
        along, cross = _extent(scene, ani)
        ratio = config.MODE2_SPEED_RATIO / config.MODE3_SPEED_RATIO
        assert along / cross == pytest.approx(ratio, rel=0.25)

    def test_isotropic_region_is_not_elongated(self, scene, rings_off):
        """Control: with no cap binding, the same scene grows round without
        rings, so the elongation above is the mechanism and not the fixture."""
        config.USE_ANISOTROPIC_RINGS = False
        _, iso = _run(scene)
        along, cross = _extent(scene, iso)
        assert along / cross == pytest.approx(1.0, abs=0.25)

    def test_flanks_stay_within_the_reported_budget(self, scene, rings_on, capsys):
        """The budget the run reports must actually bound the flank reach, to
        within one cluster width of granularity."""
        config.USE_ANISOTROPIC_RINGS = True
        _, ani = _run(scene)
        budget = float(capsys.readouterr().out
                       .split('cross-slope budget (')[1].split(' m)')[0])
        _, cross = _extent(scene, ani)
        assert cross <= budget + BLOCK

    def test_cross_budget_is_not_inflated_by_flank_growth(self, scene, rings_on, capsys):
        """Regression on the along-slope clock. Driving the budget from the
        global ring counter lets the flanks buy themselves time: measured on
        Jan 18 with the lateral cap removed, that version dropped 1 cluster at a
        132 m budget, against 43 clusters at 62 m once the clock was used.

        The clock is still not a pure fall-line measure: the ±45° sector counts
        a diagonal cluster as along-slope, so corner growth keeps the clock
        ticking after the fall-line front has stopped and the budget ends up
        somewhat above ratio x (fall-line reach). Bounded here at 2x the ideal,
        which the self-inflating version blew through by a factor of ~2.5.
        """
        config.USE_ANISOTROPIC_RINGS = True
        _, ani = _run(scene)
        along, cross = _extent(scene, ani)
        ratio = config.MODE3_SPEED_RATIO / config.MODE2_SPEED_RATIO
        assert cross < along
        assert cross <= 2.0 * ratio * along

    def test_a_zero_budget_is_refused_rather_than_hanging(self, scene, rings_on,
                                                          monkeypatch):
        """A zero step or ratio would defer every flank cluster forever."""
        monkeypatch.setattr(config, 'RING_STEP_M', 0.0)
        monkeypatch.setattr(config, 'MODE3_SPEED_RATIO', 0.0)
        config.USE_ANISOTROPIC_RINGS = True
        with pytest.raises(ValueError, match='positive ring_step'):
            _run(scene)
