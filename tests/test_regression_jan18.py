"""Regression guard: the default Jan 18 run must not move.

Reference is outputs/baseline/scenario_summary.csv, generated 2026-10-24 at
commit 31e5192 with the README quick-start command and all defaults.

Re-pinned 2026-10-27 for STAUCHWALL_DEG 28 -> 29 (Ron). With
--stauchwall-deg 28 the run still reproduces the old reference exactly
(2859 6019 m2 / 0.702, 5656 4911 / 0.644, 6191 5618 / 0.542). At 29 the
trigger slope gate is 31 deg and the third trigger becomes 2858.
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / 'data' / 'little_prof'

# trigger_cluster -> (release_area_m2, iou)
REFERENCE = {
    2859: (3819.0, 0.4911404831679284),
    5656: (4896.0, 0.6426230691834846),
    2858: (5360.0, 0.646098387369389),
}


def _run(out_dir, *extra):
    cmd = [
        sys.executable, '-m', 'release_areas.generate_scenarios',
        '--features-csv', str(DATA / 'features/all_start_zone_features_2026-01-18.csv'),
        '--meloche-csv',  str(DATA / 'features/meloche_features_all_2026-01-18.csv'),
        '--cluster-map',  str(DATA / 'spatial/cluster_map.tif'),
        '--dem',          str(DATA / 'dem_1m.tif'),
        '--start-zone',   str(DATA / 'boundaries/start_zone.kml'),
        '--release-poly', str(DATA / 'boundaries/avalanche_release_area_20260118.geojson'),
        '--out-dir',      str(out_dir),
        *extra,
    ]
    subprocess.run(cmd, check=True, cwd=REPO, capture_output=True)
    return pd.read_csv(out_dir / 'scenario_summary.csv')


@pytest.fixture(scope='module')
def default_summary(tmp_path_factory):
    return _run(tmp_path_factory.mktemp('jan18_default'))


def test_default_triggers(default_summary):
    assert list(default_summary['trigger_cluster']) == list(REFERENCE)


@pytest.mark.parametrize('trig', list(REFERENCE))
def test_default_area_and_iou(default_summary, trig):
    row = default_summary.set_index('trigger_cluster').loc[trig]
    area, iou = REFERENCE[trig]
    assert row['release_area_m2'] == pytest.approx(area, abs=1.0)
    assert row['iou'] == pytest.approx(iou, abs=1e-3)


def test_explicit_default_delta_is_identical(default_summary, tmp_path):
    """Passing --delta at its default value must change nothing."""
    from release_areas import config
    s = _run(tmp_path, '--delta', str(config.DELTA))
    pd.testing.assert_frame_equal(s, default_summary)
