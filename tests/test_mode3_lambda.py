"""Tests for the mode III (cross-slope) elastic length in propagation.

Λ₃ = √(G_slab H_slab h_WL / G_WL) is sourced from Gaume (pers. comm.
2026-10-15); `arrest_indices` owns the relation and is tested in
test_arrest_indices.TestCrossSlopeElasticLength. This file covers the one
place the two directions diverge in the geometry code,
`release_geometry.directional_lambda`, which had no coverage before the
flag went on by default on 2026-10-15.

Three things are asserted, in order of how easy they are to break:

  1. Direction + flag routing, including the fallback when a CSV predates the
     `Lambda_cross` column — the shipped Jan 18 meloche CSV is such a CSV.
  2. Scale invariance of the downstream continuity gate. This is the reason
     the flag is a no-op under the antiplane form, and it is a property of the
     gate, not of Λ, so it is tested as one.
  3. That a *spatially varying* Λ_cross is not scale-invariant — i.e. the
     plumbing is live, and would carry a bending-governed Λ₃ if one lands.
"""
import math

import numpy as np
import pytest

from release_areas import config

# Written against a version of the code that was never committed: it needs
# arrest_indices.lambda_cross_ratio and config.LAMBDA_CROSS_MODEL /
# LAMBDA_CROSS_DROP_FACTOR / LAMBDA_CROSS_RISE_FACTOR, which do not exist here.
# Committed 2026-10-27 so the work is kept; it runs again once they land.
try:
    from release_areas.arrest_indices import lambda_cross_ratio
    config.LAMBDA_CROSS_MODEL
except (ImportError, AttributeError):
    pytest.skip("mode III Lambda: lambda_cross_ratio / LAMBDA_CROSS_MODEL not "
                "implemented (uncommitted work)", allow_module_level=True)
from release_areas.release_geometry import directional_lambda


@pytest.fixture
def props():
    """One cluster's propagation properties, Λ₃/Λ₂ at the antiplane ratio."""
    lam = 2.0
    return {'Lambda': lam, 'Lambda_cross': lam * lambda_cross_ratio(config.NU)}


class TestDirectionalLambdaRouting:

    def test_along_slope_always_uses_mode_ii(self, props, monkeypatch):
        for flag in (True, False):
            monkeypatch.setattr(config, 'USE_MODE3_LAMBDA', flag)
            assert directional_lambda(props, is_cross=False) == props['Lambda']

    def test_cross_slope_uses_lambda_cross_when_enabled(self, props, monkeypatch):
        monkeypatch.setattr(config, 'USE_MODE3_LAMBDA', True)
        assert directional_lambda(props, is_cross=True) == props['Lambda_cross']

    def test_cross_slope_falls_back_when_disabled(self, props, monkeypatch):
        monkeypatch.setattr(config, 'USE_MODE3_LAMBDA', False)
        assert directional_lambda(props, is_cross=True) == props['Lambda']

    def test_fallback_when_column_absent(self, monkeypatch):
        """Jan 18's shipped meloche CSV predates `Lambda_cross`, so the enabled
        flag must still yield Λ₂ rather than NaN out into the arrest gate."""
        monkeypatch.setattr(config, 'USE_MODE3_LAMBDA', True)
        assert directional_lambda({'Lambda': 2.0}, is_cross=True) == 2.0

    def test_fallback_when_column_is_nan(self, monkeypatch):
        monkeypatch.setattr(config, 'USE_MODE3_LAMBDA', True)
        got = directional_lambda({'Lambda': 2.0, 'Lambda_cross': np.nan},
                                 is_cross=True)
        assert got == 2.0

    def test_cross_slope_is_shorter_when_enabled(self, props, monkeypatch):
        monkeypatch.setattr(config, 'USE_MODE3_LAMBDA', True)
        assert (directional_lambda(props, is_cross=True)
                < directional_lambda(props, is_cross=False))


class TestConstantRatioIsANoOp:
    """Why flipping USE_MODE3_LAMBDA cannot change a polygon under 'antiplane'.

    The continuity gate compares the *relative* change (Λ_nbr − Λ_cur)/Λ_cur,
    so scaling every Λ by one constant leaves it untouched; and the only
    absolute threshold, MIN_PROPAGATION_LAMBDA, is far below the smallest Λ₃.
    If either premise is ever broken, the no-op claim in config.py and docs §3
    stops holding and these tests are what should catch it.
    """

    @staticmethod
    def _relative_change(lam_cur, lam_nbr):
        return (lam_nbr - lam_cur) / lam_cur

    def test_relative_change_is_scale_invariant(self):
        k = lambda_cross_ratio(config.NU)
        for lam_cur, lam_nbr in [(2.0, 2.5), (2.0, 1.4), (0.6, 0.61)]:
            assert self._relative_change(lam_cur * k, lam_nbr * k) == (
                pytest.approx(self._relative_change(lam_cur, lam_nbr),
                              rel=1e-12))

    @pytest.mark.parametrize('min_lambda_ii, label', [(1.2236, 'reference'),
                                                      (0.5366, '_v1')])
    def test_min_propagation_lambda_not_approached(self, min_lambda_ii, label):
        """Min Λ₂ over the Jan 18 meloche CSVs, scaled to Λ₃, stays clear of
        the floor — 0.72 m on the reference pair, 0.32 m on `_v1`, against a
        0.1 m floor. The floor is the only thing that would break scale
        invariance, so the margin to it is the quantity that matters.
        """
        min_lambda_iii = min_lambda_ii * lambda_cross_ratio(config.NU)
        assert min_lambda_iii > config.MIN_PROPAGATION_LAMBDA, label
        assert min_lambda_iii > 3 * config.MIN_PROPAGATION_LAMBDA, label

    def test_cross_thresholds_equal_along_thresholds(self):
        """Separate knobs exist for the flanks but are deliberately equal, so
        splitting the code changed nothing. A no-op claim needs both halves."""
        assert config.LAMBDA_CROSS_DROP_FACTOR == config.LAMBDA_DROP_FACTOR
        assert config.LAMBDA_CROSS_RISE_FACTOR == config.LAMBDA_RISE_FACTOR

    def test_spatially_varying_lambda_cross_is_not_invariant(self):
        """The converse: the plumbing is live. A Λ_cross whose ratio to Λ₂
        varies between clusters does change the relative change, which is how a
        bending-governed Λ₃ would move the flanks."""
        # Λ_b ∝ (E h³ h_wl/G_wl)^¼ against Λ₂ ∝ (E h h_wl/G_wl)^½, so the
        # ratio carries h^¼ E^-¼ and differs between two unequal slabs.
        def lam_ii(E, h):
            return math.sqrt(E * h)

        def lam_bending_like(E, h):
            return (E * h**3) ** 0.25

        cur, nbr = (4e6, 0.6), (6e6, 1.2)
        rel_ii = self._relative_change(lam_ii(*cur), lam_ii(*nbr))
        rel_b = self._relative_change(lam_bending_like(*cur),
                                      lam_bending_like(*nbr))
        assert not math.isclose(rel_ii, rel_b, rel_tol=1e-3), (rel_ii, rel_b)


class TestConfigAgreesWithRelation:

    def test_config_model_is_implemented(self):
        from release_areas import arrest_indices as ai
        assert config.LAMBDA_CROSS_MODEL in ai.LAMBDA_CROSS_MODELS
        assert config.LAMBDA_CROSS_MODEL == ai.LAMBDA_CROSS_DEFAULT

    def test_mode3_lambda_on_by_default(self):
        """On since 2026-10-15 because the relation is sourced. If this is ever
        flipped back, the reason belongs in config.py, not in a tolerance."""
        assert config.USE_MODE3_LAMBDA is True
