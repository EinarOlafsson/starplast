"""Whole-group bootstrap intervals for aligned numerator/denominator records.

These descriptive intervals assume exchangeable supplied groups. They do not
prove biological independence or estimate per-gene confidence.
"""
from __future__ import annotations

import hashlib
import json
import numpy as np


LIMITS = ('Conditional on the supplied recorded groups and frozen predictions. '
          'Group definitions do not establish biological or source independence. '
          'Undefined resampled denominators are counted and excluded from percentiles. '
          'These intervals are not calibrated per-gene probabilities.')
PAIRED_STATE_SHA256 = '0bc3e966f6ee493185f1189735245852b8b910a6548abdcb96496e376ef061f9'


def load_paired_state(report, *, path=None, expected_sha256=PAIRED_STATE_SHA256):
    """Read the pinned interval artifact for the supported frozen paired report."""
    from pathlib import Path
    from . import functional_agreement as A, functional_results as F
    if hashlib.sha256((F._canonical(report) + '\n').encode()).hexdigest() != A.REPORT_SHA256:
        raise ValueError('Intervals require the pinned paired source report')
    path = Path(path) if path is not None else Path(__file__).with_name('data') / 'functional_agreement_uncertainty.json'
    encoded = path.read_bytes()
    if not expected_sha256 or hashlib.sha256(encoded).hexdigest() != expected_sha256:
        raise ValueError('Paired interval checksum changed')
    state = json.loads(encoded)
    return state


def inputs(groups, numerators, denominators, names):
    """Validate integer ratio records and return their canonical input identity."""
    groups = list(groups); names = list(names)
    if (not groups or any(not isinstance(g, str) or not g for g in groups)
            or not names or any(not isinstance(n, str) or not n for n in names)
            or len(set(names)) != len(names)):
        raise ValueError('Explicit nonempty groups and unique rate names required')
    a = np.asarray(numerators); b = np.asarray(denominators)
    if (a.shape != b.shape or a.shape != (len(groups), len(names))
            or a.dtype.kind not in 'iu' or b.dtype.kind not in 'iu'
            or (b < 0).any()):
        raise ValueError('Matched integer numerator/denominator matrices required')
    # Convert to Python integers before bounding, avoiding signed/unsigned overflow.
    aa, bb = a.tolist(), b.tolist()
    if any(abs(x) > y for row_a, row_b in zip(aa, bb) for x, y in zip(row_a, row_b)):
        raise ValueError('Per-record numerator magnitude exceeds denominator')
    if any(abs(x) > np.iinfo(np.int64).max // (len(groups) * len(set(groups)))
           for rows in (aa, bb) for row in rows for x in row):
        raise ValueError('Ratio count totals exceed supported integer capacity')
    encoded = json.dumps({'groups': groups, 'names': names, 'numerators': aa,
                          'denominators': bb}, sort_keys=True, separators=(',', ':')).encode()
    return groups, np.asarray(aa, dtype=np.int64), np.asarray(bb, dtype=np.int64), names, hashlib.sha256(encoded).hexdigest()


def estimate(groups, numerators, denominators, names, *, seed, draws=500, selections=None):
    """Return a JSON summary and replay arrays using the same draws for every rate.

    ``selections`` replays saved group indices; it must retain the requested shape
    and whole-group sampling unit. Zero-denominator rates remain unavailable.
    """
    groups, a, b, names, identity = inputs(groups, numerators, denominators, names)
    if type(seed) is not int or seed < 0 or type(draws) is not int or draws < 2:
        raise ValueError('Nonnegative integer seed and at least two draws required')
    unique = sorted(set(groups)); index = {g: i for i, g in enumerate(unique)}
    group_a = np.zeros((len(unique), len(names)), dtype=np.int64)
    group_b = np.zeros_like(group_a)
    for g, x, y in zip(groups, a, b):
        group_a[index[g]] += x; group_b[index[g]] += y
    sufficient = len(unique) >= 5
    shape = (draws if sufficient else 0, len(unique))
    if selections is None:
        selections = np.random.Generator(np.random.PCG64(seed)).integers(
            0, len(unique), size=shape, dtype=np.int64)
    else:
        selections = np.asarray(selections)
        if (selections.shape != shape or selections.dtype.kind not in 'iu'
                or (selections < 0).any() or (selections >= len(unique)).any()):
            raise ValueError('Replay must contain valid whole-group selections')
        selections = selections.astype(np.int64)
    boot_a = group_a[selections].sum(axis=1); boot_b = group_b[selections].sum(axis=1)
    values = np.full(boot_a.shape, np.nan)
    np.divide(boot_a, boot_b, out=values, where=boot_b != 0)
    rates = {}
    for i, name in enumerate(names):
        numerator = int(a[:, i].sum()); denominator = int(b[:, i].sum())
        valid = values[:, i][boot_b[:, i] != 0]
        available = sufficient and denominator > 0 and len(valid) >= 2
        rates[name] = {'numerator': numerator, 'denominator': denominator,
            'value': numerator / denominator if denominator else None,
            'valid_draws': len(valid), 'undefined_draws': len(values) - len(valid),
            'interval': np.quantile(valid, [.025, .975], method='linear').tolist() if available else None,
            'status': 'descriptive_group_bootstrap' if available else 'unavailable'}
    state = {'schema_version': 1, 'status': 'descriptive_group_bootstrap' if sufficient else 'unavailable',
             'reason': None if sufficient else 'Fewer than five recorded groups',
             'input_identity': identity, 'recorded_groups': len(unique), 'entities': len(groups),
             'group_order': unique, 'seed': seed, 'draws_requested': draws, 'draws_evaluated': len(values),
             'generator': 'PCG64', 'quantiles': [.025, .975], 'quantile_method': 'linear',
             'resampling_unit': 'whole recorded group', 'rates': rates, 'limits': LIMITS}
    arrays = {'selections': selections, 'numerators': boot_a, 'denominators': boot_b, 'ratios': values}
    return state, arrays


def paired_terms(report):
    """Extract paired call masks without excluding abstentions or unsupported genes."""
    rows = report['rows']; numerators = []; denominators = []
    names = ['joint_coverage', 'agreement_among_joint_calls', 'source_recovery_among_agreements',
             'correct_agreements_all_eligible', 'shared_wrong_among_joint_calls',
             'right_minus_left_recovery_all_eligible']
    for row in rows:
        x, y, t = (row[key] for key in ('left_prediction', 'right_prediction', 'truth'))
        joint = x is not None and y is not None
        agree = x is not None and x == y
        correct = agree and x == t
        wrong = joint and x != t and y != t
        numerators.append([int(joint), int(agree), int(correct), int(correct), int(wrong), int(y == t) - int(x == t)])
        denominators.append([1, int(joint), int(agree), 1, int(joint), 1])
    return [row['group'] for row in rows], numerators, denominators, names


def paired_scorecard(report, state):
    """Present verified group intervals while preserving exact paired source rates."""
    from copy import deepcopy
    from dataclasses import replace
    from . import functional_agreement as A, functional_results as F
    view = A.build_scorecard(report)
    expected, _ = estimate(*paired_terms(report), seed=state['seed'], draws=state['draws_requested'])
    if state != expected:
        raise ValueError('Group intervals do not reproduce the supplied paired records')
    for name, rate in report['rates'].items():
        if any(state['rates'][name][key] != rate[key] for key in ('numerator', 'denominator', 'value')):
            raise ValueError('Group interval point estimate differs from paired source rate')
    limits = view.source.lineage.replace('Group-aware uncertainty is not estimated in this report.', LIMITS)
    limits += ' Resampling assumes exchangeable recorded groups.'
    snapshot = deepcopy(report)
    snapshot['group_uncertainty'] = deepcopy(state)
    snapshot['interpretation'] = limits
    display_state={**state,'status':state['status'].replace('_',' '),
        'reason':state['reason'] or f"{state['draws_evaluated']} whole-group resamples; {state['recorded_groups']} recorded groups. Conditional on these groups and frozen predictions."}
    details = tuple(replace(detail, text=F._canonical(state)) if detail.key == 'group_uncertainty' else
                    replace(detail, text=F._canonical(display_state)) if detail.key == 'uncertainty' else
                    replace(detail, text=F._canonical({'status': 'unavailable',
                        'reason': 'No per-gene probability calibration; intervals describe group sampling variation'}))
                    if detail.key == 'calibration' else
                    replace(detail, text=F._canonical(limits)) if detail.key == 'gaps' else detail
                    for detail in view.details)
    evidence = {}
    for metric in view.metrics:
        rate = state['rates'][metric.key]
        bounds = rate['interval']
        interval = (f"{bounds[0]:.1%}–{bounds[1]:.1%}" if bounds is not None else 'unavailable')
        evidence[metric.label] = (f"{rate['numerator']} / {rate['denominator']} "
            f"({rate['value']:.1%}); 95% descriptive group interval: {interval}" if rate['value'] is not None else 'Unavailable')
    delta=state['rates']['right_minus_left_recovery_all_eligible']
    if delta['value'] is not None:
        label=' minus '.join(method['strategy'].replace('_',' ') for method in reversed(report['methods']))
        bounds=delta['interval']
        interval=f'{bounds[0]*100:.1f}–{bounds[1]*100:.1f} points' if bounds is not None else 'unavailable'
        evidence[label+' source recovery']=(f"{delta['numerator']} / {delta['denominator']} "
            f"({delta['value']*100:.1f} percentage points); 95% descriptive group interval: {interval}")
    details = tuple(replace(detail, text=F._canonical(evidence)) if detail.key == 'evidence'
                    else detail for detail in details)
    return replace(view, source=replace(view.source, lineage=limits), details=details,
                   metrics=tuple(replace(metric, reading=limits) for metric in view.metrics),
                   snapshot_json=F._canonical(snapshot))
