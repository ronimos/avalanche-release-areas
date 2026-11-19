"""examples/evaluate_slab_crack_modes.py: the variant list covers every option
in config, the scorer parses the axes part, and pooling counts each event once."""
import importlib.util
from pathlib import Path

import pandas as pd
import pytest

from release_areas import config

REPO = Path(__file__).resolve().parents[1]


def _load(name):
    spec = importlib.util.spec_from_file_location(name, REPO / 'examples' / f'{name}.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_all_variants_cover_every_option():
    ev = _load('evaluate_slab_crack_modes')
    n_mixed = len(config.MIXED_CRITERIA) * len(config.MIXED_ON_FAILS) * len(config.MIXED_AXES_OPTIONS)
    assert len(ev.ALL_VARIANTS) == len(config.SLAB_CRACK_MODES) - 1 + n_mixed
    assert {'shear', 'opening', 'mixed:energy:intact:fall',
            'mixed:elliptic:stop:front_edge'} <= set(ev.ALL_VARIANTS)


def test_set_mode_parses_axes(monkeypatch):
    sj = _load('score_jan18')
    for k in ('SLAB_CRACK_MODE', 'MIXED_CRITERION', 'MIXED_ON_FAIL', 'MIXED_AXES'):
        monkeypatch.setattr(config, k, getattr(config, k))
    sj.set_mode('mixed:elliptic:intact:fall')
    assert (config.SLAB_CRACK_MODE, config.MIXED_CRITERION, config.MIXED_ON_FAIL,
            config.MIXED_AXES) == ('mixed', 'elliptic', 'intact', 'fall')
    sj.set_mode('mixed:energy:stop')                  # axes left out: the default
    assert config.MIXED_AXES == 'fall'


def test_pooled_counts_each_event_once():
    ev = _load('evaluate_slab_crack_modes')
    df = pd.DataFrame([
        dict(event='a', variant='x', e_iou=0.8, area_ratio=1.0, auc=0.9, bss=0.5,
             inside_iou=0.7, ski_n=10, ski_released=0),
        dict(event='a', variant='y', e_iou=0.6, area_ratio=1.2, auc=0.8, bss=0.3,
             inside_iou=0.75, ski_n=10, ski_released=1),
        dict(event='b', variant='x', e_iou=0.4, area_ratio=0.8, auc=0.7, bss=0.1,
             inside_iou=0.5, ski_n=0, ski_released=0),
        dict(event='b', variant='y', e_iou=0.5, area_ratio=0.9, auc=0.6, bss=0.2,
             inside_iou=0.4, ski_n=0, ski_released=0),
    ])
    p = ev.pooled(df).set_index('variant')
    assert p.loc['x', 'e_iou'] == pytest.approx(0.6)
    assert p.loc['x', 'rank_e_iou'] == pytest.approx(1.5)
    assert p.loc['y', 'rank_inside_iou'] == pytest.approx(1.5)
    assert p.loc['y', 'ski_released'] == 1 and p.loc['x', 'n_events'] == 2
