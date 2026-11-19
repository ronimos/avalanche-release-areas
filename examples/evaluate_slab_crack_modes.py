#!/usr/bin/env python3
"""Slab-fracture-mode evaluation across every registered avalanche.

Which slab-crack mode is right ('opening', 'shear', 'mixed' with its criterion,
failure rule and split axes) is open, and Jan 18 alone cannot decide it: it is
one event and a calibration target. This script re-runs the comparison
whenever an event is added to data/events.json, so the choice can be made on
many events instead of one.

For each event it runs examples/score_jan18.py (five-trigger score plus the
inside-release forced shots) for every variant, in its own process so one
event's calibration preset cannot leak into the next. It then pools the
results per variant: the mean over events of five-trigger E[IoU], E[area] /
observed, AUC, Brier skill and inside-release median IoU, each variant's mean
rank, and the ski non-trigger releases. Validation events are pooled apart
from calibration events; calibration scores are shown but never count as
validation.

    env -u VIRTUAL_ENV uv run --frozen python examples/evaluate_slab_crack_modes.py \\
        --out outputs/slab_mode_eval --jobs 24

Add an event: put its data under data/<slope>/ (the little_prof layout) and
append an entry to data/events.json. --reuse skips events already scored with
the same variants.
"""
from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
SCORER = REPO / 'examples' / 'score_jan18.py'

# Every slab-crack variant the model has (config.SLAB_CRACK_MODES, and for
# 'mixed' MIXED_CRITERIA x MIXED_ON_FAILS x MIXED_AXES_OPTIONS).
sys.path.insert(0, str(REPO / 'src'))
from release_areas import config  # noqa: E402

ALL_VARIANTS = [m for m in config.SLAB_CRACK_MODES if m != 'mixed'] + [
    f'mixed:{c}:{f}:{a}' for c, f, a in itertools.product(
        config.MIXED_CRITERIA, config.MIXED_ON_FAILS, config.MIXED_AXES_OPTIONS)]

METRICS = [  # (column, label, higher is better)
    ('e_iou', 'E[IoU]', True),
    ('area_ratio', 'E[area]/obs', None),
    ('auc', 'AUC', True),
    ('bss', 'Brier skill', True),
    ('inside_iou', 'inside median IoU', True),
]


def score_event(ev, variants, out, a):
    slope = REPO / 'data' / ev['slope']
    edir = out / ev['id']
    summ = edir / 'summary.json'
    if a.reuse and summ.exists():
        have = json.loads(summ.read_text())['meta'].get('modes', [])
        if set(variants) <= set(have):
            print(f"{ev['id']}: reusing {summ}", flush=True)
            return json.loads(summ.read_text())
    cmd = [sys.executable, str(SCORER), '--data-dir', str(slope), '--date', ev['date'],
           '--calibration', str(slope / ev['calibration']),
           '--observed', str(slope / ev['observed']),
           '--ski', str(slope / ev['ski']) if ev.get('ski') else 'none',
           '--modes', *variants, '--jobs', str(a.jobs), '--subsample', str(a.subsample),
           '--out', str(edir)]
    if a.no_inside:
        cmd.append('--no-inside')
    print(f"{ev['id']}: scoring {len(variants)} variants", flush=True)
    subprocess.run(cmd, check=True, stdout=None if a.verbose else subprocess.DEVNULL)
    return json.loads(summ.read_text())


def event_rows(ev, summ):
    obs = summ['meta']['observed_m2']
    inside = {s['mode']: s for s in summ['forced_shots'] if s['class'] == 'inside'}
    ski = {s['mode']: s for s in summ['forced_shots'] if s['class'] == 'ski'}
    rows = []
    for s in summ['five_trigger']:
        m = s['mode']
        rows.append({'event': ev['id'], 'role': ev.get('role', 'validation'), 'variant': m,
                     'e_iou': s['expected_iou'], 'area_ratio': s['expected_area_m2'] / obs,
                     'auc': s['auc'], 'bss': s['brier_skill'],
                     'inside_iou': inside.get(m, {}).get('iou_median', np.nan),
                     'ski_n': ski.get(m, {}).get('n', 0),
                     'ski_released': ski.get(m, {}).get('n_released', 0),
                     'dead_triggers': len(s['dead_triggers'])})
    return rows


def pooled(df):
    out = []
    for v, g in df.groupby('variant', sort=False):
        rec = {'variant': v, 'n_events': g['event'].nunique()}
        for col, _, _ in METRICS:
            rec[col] = g[col].mean()
        rec['ski_released'] = int(g['ski_released'].sum())
        rec['ski_n'] = int(g['ski_n'].sum())
        out.append(rec)
    p = pd.DataFrame(out)
    # mean rank over events (1 = best) on E[IoU] and inside median IoU
    for col in ('e_iou', 'inside_iou'):
        r = df.groupby('event')[col].rank(ascending=False, method='average')
        p[f'rank_{col}'] = df.assign(r=r).groupby('variant', sort=False)['r'].mean().values
    return p


def table(p, title, note):
    f = lambda x, n=3: '-' if pd.isna(x) else f'{x:.{n}f}'
    lines = [f'## {title}', ''] + ([note, ''] if note else []) + [
        '| variant | events | E[IoU] | rank | E[area]/obs | AUC | Brier skill | '
        'inside median IoU | rank | ski released |',
        '|---|---|---|---|---|---|---|---|---|---|']
    for _, r in p.sort_values('rank_e_iou').iterrows():
        lines.append(f"| {r['variant']} | {r['n_events']} | {f(r['e_iou'])} | "
                     f"{f(r['rank_e_iou'], 1)} | {f(r['area_ratio'], 2)} | {f(r['auc'])} | "
                     f"{f(r['bss'])} | {f(r['inside_iou'])} | {f(r['rank_inside_iou'], 1)} | "
                     f"{r['ski_released']}/{r['ski_n']} |")
    return lines + ['']


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--events', type=Path, default=REPO / 'data' / 'events.json')
    ap.add_argument('--only', nargs='*', help='event ids (default: all)')
    ap.add_argument('--variants', nargs='+', default=ALL_VARIANTS,
                    help='score_jan18 variant strings (default: all %d)' % len(ALL_VARIANTS))
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--subsample', type=int, default=1)
    ap.add_argument('--no-inside', action='store_true')
    ap.add_argument('--reuse', action='store_true')
    ap.add_argument('--verbose', action='store_true')
    a = ap.parse_args(argv)
    events = json.loads(a.events.read_text())['events']
    if a.only:
        events = [e for e in events if e['id'] in a.only]
    a.out.mkdir(parents=True, exist_ok=True)

    rows = []
    for ev in events:
        rows += event_rows(ev, score_event(ev, a.variants, a.out, a))
    df = pd.DataFrame(rows)
    df.to_csv(a.out / 'per_event.csv', index=False)

    n_val = df.loc[df.role == 'validation', 'event'].nunique()
    n_cal = df.loc[df.role == 'calibration', 'event'].nunique()
    lines = ['# Slab-fracture-mode evaluation', '',
             f'{len(events)} event(s): {n_val} validation, {n_cal} calibration. '
             f'{len(a.variants)} variants. Variant strings: '
             "`mixed:<criterion>:<on_fail>:<axes>` (config.MIXED_*). "
             'Rank = mean over events of the per-event rank (1 = best). '
             'Mean over events, not over shots: every event counts once.', '']
    if n_val:
        lines += table(pooled(df[df.role == 'validation']), 'Validation events',
                       'Events no model knob was fitted on. The comparison that counts.')
    else:
        lines += ['## Validation events', '',
                  'None yet. Every score below is from a calibration event, so it can '
                  'describe the variants but cannot choose between them.', '']
    if n_cal:
        lines += table(pooled(df[df.role == 'calibration']), 'Calibration events',
                       'A model knob was fitted on these events: description, not evidence.')
    if df['event'].nunique() > 1:
        lines += ['## Per event', '']
        for ev_id, g in df.groupby('event', sort=False):
            lines += table(pooled(g), f'{ev_id} ({g.role.iloc[0]})', '')
    (a.out / 'summary.md').write_text('\n'.join(lines))
    print('\n'.join(lines))
    print(f'-> {a.out}/summary.md, per_event.csv')


if __name__ == '__main__':
    sys.exit(main())
