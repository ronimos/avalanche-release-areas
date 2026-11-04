"""Release density comes from the clusters inside the release polygon.

Ported from avachain (analysis_pipeline._polygon_density, 2026-11-04), one of
the behaviours where avachain is deliberately better (avachain CLAUDE.md):
this script used the trigger cluster's slab density for every scenario of
that trigger, with 250 kg/m^3 as the fallback.
"""
import numpy as np
import pandas as pd
import pytest
from rasterio.transform import from_origin
from shapely.geometry import box

from release_areas.generate_scenarios import (DEFAULT_SLAB_DENSITY,
                                               DEFAULT_SLAB_DENSITY_STD,
                                               polygon_density)

# 10 x 10 m grid, 1 m cells; clusters 1 | 2 | 3 as vertical bands
CM = np.zeros((10, 10), dtype=int)
CM[:, 0:3], CM[:, 3:6], CM[:, 6:10] = 1, 2, 3
TR = from_origin(1000.0, 2000.0, 1.0, 1.0)


def _feats(**rho):
    return pd.DataFrame({"slab_density": rho}).rename_axis("cluster_id")


def test_median_over_the_covered_clusters():
    feats = _feats(**{"1": 200.0, "2": 300.0, "3": 400.0}).rename(index=int)
    poly = box(1000.0, 1990.0, 1006.0, 2000.0)          # covers clusters 1 and 2
    mean, std = polygon_density(poly, CM, feats, TR)
    assert mean == pytest.approx(250.0)
    assert std == pytest.approx(pd.Series([200.0, 300.0]).std())


def test_not_the_trigger_cluster_alone():
    # the old rule would give cluster 1's 200 for any polygon it triggered
    feats = _feats(**{"1": 200.0, "2": 300.0, "3": 400.0}).rename(index=int)
    mean, _ = polygon_density(box(1000.0, 1990.0, 1010.0, 2000.0), CM, feats, TR)
    assert mean == pytest.approx(300.0)


def test_no_density_falls_back():
    feats = _feats(**{"1": np.nan, "2": np.nan, "3": np.nan}).rename(index=int)
    assert polygon_density(box(1000.0, 1990.0, 1010.0, 2000.0), CM, feats, TR) == \
        (DEFAULT_SLAB_DENSITY, DEFAULT_SLAB_DENSITY_STD)
    assert polygon_density(box(1000.0, 1990.0, 1010.0, 2000.0), CM,
                           pd.DataFrame(index=[1, 2, 3]), TR) == \
        (DEFAULT_SLAB_DENSITY, DEFAULT_SLAB_DENSITY_STD)


def test_one_cluster_keeps_the_default_spread():
    feats = _feats(**{"1": 210.0, "2": np.nan, "3": np.nan}).rename(index=int)
    assert polygon_density(box(1000.0, 1990.0, 1010.0, 2000.0), CM, feats, TR) == \
        (210.0, DEFAULT_SLAB_DENSITY_STD)


def test_matches_avachain_defaults():
    assert (DEFAULT_SLAB_DENSITY, DEFAULT_SLAB_DENSITY_STD) == (280.0, 20.0)
