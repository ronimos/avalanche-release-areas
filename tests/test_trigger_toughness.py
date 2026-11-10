"""Trigger fracture toughness gate and the BFS distance-caps switch (Ron,
2026-11-10; same as ~/avachain tests/test_trigger_toughness.py). Jan 18 cannot
see the toughness gate (its five triggers pass with or without it)."""
import numpy as np
import pandas as pd
import pytest

from release_areas import arrest_indices as ai
from release_areas import config
from release_areas import release_geometry as rg


def test_toughness_limited_strength_formula():
    K, h, st = 1650.0, 1.1, 20e3
    want = K / (1.122 * np.sqrt(np.pi * 0.5 * h))
    assert ai.toughness_limited_strength(st, K, h) == pytest.approx(want)
    assert ai.toughness_limited_strength(100.0, K, h) == 100.0
    assert ai.toughness_limited_strength(st, np.nan, h) == st
    assert ai.toughness_limited_strength(st, K, np.nan) == st


def _frames(K):
    cid = 7
    snap = pd.DataFrame({config.SLAB_TENSILE_SIGMA_COL: [8000.0], 'slab_density': [300.0],
                         'slab_thickness': [1.0], config.TOUGHNESS_K_COL: [K]}, index=[cid])
    mel = pd.DataFrame({'rc_wl': [0.5]}, index=[cid])
    return cid, mel, snap


def test_trigger_support_uses_toughness(monkeypatch):
    cid, mel, snap = _frames(800.0)
    monkeypatch.setattr(config, "USE_TRIGGER_TOUGHNESS", False)
    plain = rg.cluster_tensile_support(cid, mel, snap)
    monkeypatch.setattr(config, "USE_TRIGGER_TOUGHNESS", True)
    tough = rg.cluster_tensile_support(cid, mel, snap)
    assert plain > 1.0 > tough
    s_eff = 800.0 / (1.122 * np.sqrt(np.pi * 0.5))
    assert tough == pytest.approx(plain * s_eff / 8000.0)


def test_trigger_support_without_k_column_is_plain(monkeypatch):
    cid, mel, snap = _frames(np.nan)
    monkeypatch.setattr(config, "USE_TRIGGER_TOUGHNESS", True)
    a = rg.cluster_tensile_support(cid, mel, snap)
    monkeypatch.setattr(config, "USE_TRIGGER_TOUGHNESS", False)
    assert rg.cluster_tensile_support(cid, mel, snap) == pytest.approx(a)


def test_defaults_match_avachain():
    assert config.USE_TRIGGER_TOUGHNESS is True
    assert config.USE_DISTANCE_CAPS is True
    assert config.TOUGHNESS_K_COL == 'K_Ic_ligament'
    assert (ai.TOUGHNESS_EDGE_F, ai.TOUGHNESS_A_FRAC) == (1.122, 0.5)
    assert ai.K_IC_RANGE_MODE == 'extrapolate'
