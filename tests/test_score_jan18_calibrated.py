"""Regression guard: Jan 18 at avachain calibration 2026-11-17.1, scored here.

examples/score_jan18.py on data/little_prof/features/cal_2026-11-17.1 (avachain
analyze, facet_a 0.75, 00:00) reproduces avachain's recorded five-trigger
scores exactly (slope_config.toml [calibration] "fitted", and the 11-16 knob
sweep at facet_a 0.75 / tau_ratio_floor 0.65). Verified 2026-11-18. Any new
slab-crack mode must leave these two untouched.
"""
import importlib.util
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

# mode -> (expected_iou, expected_area_m2, auc, brier_skill)
EXPECTED = {
    'shear':   (0.662, 6233.0, 0.850, 0.368),
    'opening': (0.662, 6482.0, 0.850, 0.341),
}
TRIGGERS = [1772, 37, 1387, 6402, 565]


@pytest.fixture(scope='module')
def scorer():
    import json
    from release_areas import config
    spec = importlib.util.spec_from_file_location('score_jan18', REPO / 'examples' / 'score_jan18.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    saved = {k: getattr(config, k) for k in (
        'FACET_STRENGTH_FACTOR', 'USE_DISTANCE_CAPS', 'TAU_RATIO_FLOOR',
        'SLAB_CRACK_MODE', 'USE_DOWNSLOPE_SLIDE', 'UNIFIED_SLIDE')}
    cal = json.loads((mod.DATA / 'calibration_2026-11-17.1.json').read_text())
    mod.apply_calibration(cal)
    mod.load_inputs(mod.DATA / cal['features_dir'],
                    mod.DATA / 'validation' / 'survey_release_area_20260118.geojson', None)
    yield mod
    for k, v in saved.items():
        setattr(config, k, v)


@pytest.mark.parametrize('mode', list(EXPECTED))
def test_five_trigger_matches_avachain(scorer, mode):
    s = scorer.five_trigger(mode, 5)
    assert s['triggers'] == TRIGGERS
    assert s['dead_triggers'] == []
    e_iou, e_area, auc, bss = EXPECTED[mode]
    assert s['expected_iou'] == pytest.approx(e_iou, abs=5e-4)
    assert s['expected_area_m2'] == pytest.approx(e_area, abs=0.5)
    assert s['auc'] == pytest.approx(auc, abs=5e-4)
    assert s['brier_skill'] == pytest.approx(bss, abs=5e-4)
