"""Tests for the element-geometry helpers in snowpack_features.

The invariant that matters: SNOWPACK's z is exactly height_above_ground - HS,
so differencing sorted z with the ground prepended at -HS recovers each
element's own thickness and the thicknesses sum to HS.
"""
import numpy as np
import pytest

from release_areas.snowpack_features import (
    element_thickness,
    split_wl_slab,
    thickness_weighted_mean,
)

# A 4-element profile: thicknesses 0.10 / 0.20 / 0.05 / 0.15 m, HS = 0.50 m.
# Element tops above ground -> 0.10, 0.30, 0.35, 0.50; z = top - HS.
HS = 0.50
DZ_TRUE = np.array([0.10, 0.20, 0.05, 0.15])
Z = np.cumsum(DZ_TRUE) - HS          # [-0.40, -0.20, -0.15, 0.00]


class TestElementThickness:
    def test_recovers_each_element_thickness(self):
        dz = element_thickness(Z, HS, np.ones(4, bool))
        assert dz == pytest.approx(DZ_TRUE)

    def test_sums_to_hs(self):
        dz = element_thickness(Z, HS, np.ones(4, bool))
        assert float(np.sum(dz)) == pytest.approx(HS)

    def test_basal_element_measured_from_ground(self):
        # The lowest element is NOT the gap to the next one up: it spans from
        # the ground (z = -HS) to its own top. This is the element the old
        # max-min span dropped.
        dz = element_thickness(Z, HS, np.ones(4, bool))
        assert dz[0] == pytest.approx(Z[0] + HS)

    def test_returns_in_original_index_order(self):
        order = np.array([2, 0, 3, 1])
        dz = element_thickness(Z[order], HS, np.ones(4, bool))
        assert dz == pytest.approx(DZ_TRUE[order])

    def test_nan_outside_valid_mask(self):
        valid = np.array([True, True, True, False])
        dz = element_thickness(Z, HS, valid)
        assert np.isnan(dz[3])
        # Excluding the top element leaves the others untouched and the sum
        # short by exactly that element - the real effect of the (z != 0) filter.
        assert np.nansum(dz) == pytest.approx(HS - DZ_TRUE[3])
        assert dz[:3] == pytest.approx(DZ_TRUE[:3])

    def test_empty_or_bad_hs_gives_all_nan(self):
        assert np.all(np.isnan(element_thickness(Z, HS, np.zeros(4, bool))))
        assert np.all(np.isnan(element_thickness(Z, np.nan, np.ones(4, bool))))


class TestThicknessWeightedMean:
    def test_matches_plain_mean_for_uniform_elements(self):
        v = np.array([100.0, 200.0, 300.0])
        dz = np.full(3, 0.05)
        assert thickness_weighted_mean(v, np.ones(3, bool), dz) == pytest.approx(v.mean())

    def test_weights_by_thickness(self):
        # A thick dense layer must dominate a thin light one.
        v = np.array([100.0, 300.0])
        dz = np.array([0.01, 0.30])
        got = thickness_weighted_mean(v, np.ones(2, bool), dz)
        assert got == pytest.approx((100 * 0.01 + 300 * 0.30) / 0.31)
        assert got > v.mean()

    def test_ignores_masked_nan_and_nonpositive_dz(self):
        v = np.array([100.0, np.nan, 300.0, 999.0])
        dz = np.array([0.10, 0.10, 0.0, 0.10])
        mask = np.array([True, True, True, False])
        assert thickness_weighted_mean(v, mask, dz) == pytest.approx(100.0)

    def test_nan_when_nothing_usable(self):
        v = np.array([1.0, 2.0])
        assert np.isnan(thickness_weighted_mean(v, np.zeros(2, bool), np.ones(2)))
        assert np.isnan(thickness_weighted_mean(v, np.ones(2, bool), np.zeros(2)))


class TestWeakLayerThickness:
    """D_wl must be the sum of the WL elements' thicknesses."""

    # Basal WL = two depth-hoar elements (5xx), slab = two rounded-grain (3xx).
    GT = np.array([510.0, 510.0, 320.0, 320.0])

    def test_sum_equals_hs_minus_slab_thickness(self):
        slab_m, wl_m, interface_z = split_wl_slab(self.GT, Z)
        ok = np.ones(4, bool)
        dz = element_thickness(Z, HS, ok)
        d_wl = float(np.nansum(dz[wl_m & ok]))
        assert d_wl == pytest.approx(DZ_TRUE[:2].sum())
        # The identity the CSV patch relies on.
        assert d_wl == pytest.approx(HS - (-interface_z))

    def test_span_undercounts_by_the_basal_element(self):
        slab_m, wl_m, _ = split_wl_slab(self.GT, Z)
        wl_z = Z[wl_m]
        span = float(wl_z.max() - wl_z.min())      # the old, wrong form
        dz = element_thickness(Z, HS, np.ones(4, bool))
        assert span == pytest.approx(float(np.nansum(dz[wl_m])) - DZ_TRUE[0])
        assert span < float(np.nansum(dz[wl_m]))
