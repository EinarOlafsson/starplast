"""Frozen native conformal label sets from a separately fitted base model.

Calibration truth never fits the base model. All test sets, including empty,
multi-label and unsupported-class outcomes, remain available for efficiency and
coverage tests. Empirical coverage is not a gene-level correctness probability.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

from . import artifacts as A, strategy_learning as N
from .label_records import LabelBatch
from .splits import SplitManifest


def conformal_calls(calibration_scores, test_scores, calibration_labels, *, split, base_model, alpha=.1, per_class=True):
    """Use native quantiles/rare-class fallback on disjoint frozen calibration scores.

    ``base_model`` must be a verified train-only fitted-model artifact for native
    kNN or logistic classification. Scores have the exact complete calibration
    and outer-test populations, in frozen order, with the same native classes.
    The caller records/replays actual score generation; artifact declarations
    alone cannot establish biological source independence or exchangeability.
    """
    if not isinstance(split, SplitManifest) or not isinstance(base_model, A.Artifact):
        raise TypeError('Frozen split and verified fitted-model artifact are required')
    spec = base_model.spec
    spec.validate_split(split)
    if spec.role != 'fitted_model' or spec.strategy not in {'feature_knn', 'supervised_classifier'}:
        raise ValueError('Use a native kNN/logistic fitted-model artifact')
    if spec.fit_role != 'train' or spec.fit_entities != split.entities('train') or spec.benchmark_id != split.benchmark_id:
        raise ValueError('Base model must fit only the complete ordered training cohort')
    base_model.verify_contents()
    state = base_model.payloads.get('model_state.json', {})
    if state.get('fit_entities') != list(spec.fit_entities):
        raise ValueError('Base model payload and declared fitting population differ')
    if not isinstance(calibration_labels, pd.Series) or tuple(calibration_labels.index) != split.entities('calibration'):
        raise ValueError('Provide the complete ordered calibration-label cohort only')
    if not all(isinstance(v, str) and v for v in calibration_labels):
        raise ValueError('Calibration labels must be observed categorical values')
    split.guard_fit('calibration', calibration_labels.index)
    if isinstance(alpha, bool) or not isinstance(alpha, (int, float)) or not 0 < alpha < 1 or type(per_class) is not bool:
        raise ValueError('Declare finite conformal alpha and threshold policy')
    for frame, role in ((calibration_scores, 'calibration'), (test_scores, 'test')):
        if not isinstance(frame, pd.DataFrame) or tuple(frame.index) != split.entities(role) or frame.columns.has_duplicates:
            raise ValueError('Score population/order differs from the frozen '+role+' cohort')
        if not all(isinstance(v, str) and v for v in frame.columns):
            raise ValueError('Native class columns must be unique categorical labels')
        values = frame.to_numpy(dtype=float)
        if np.isinf(values).any() or (values[np.isfinite(values)] < 0).any() or (values[np.isfinite(values)] > 1).any():
            raise ValueError('Native class shares must be finite in [0,1] or missing')
        if values.shape[1]:
            known = np.isfinite(values).all(axis=1)
            missing = np.isnan(values).all(axis=1)
            if not (known | missing).all() or not np.allclose(values[known].sum(axis=1), 1, atol=1e-8):
                raise ValueError('Class scores need complete normalized shares or explicit unavailable rows')
    if not calibration_scores.columns.equals(test_scores.columns) or list(test_scores.columns) != state.get('native_classes'):
        raise ValueError('Native model/calibration/test class orders differ')
    classes = test_scores.columns.tolist()
    lookup = {label: j for j, label in enumerate(classes)}
    P = calibration_scores.to_numpy(dtype=float)
    nonconformity = np.array([1.0 - (P[i, lookup[label]] if label in lookup and np.isfinite(P[i, lookup[label]]) else 0.0)
        for i, label in enumerate(calibration_labels)])
    overall = N._conformal_quantile(nonconformity, alpha)
    thresholds = {label: overall for label in classes}
    fallback = []
    need = int(math.ceil(1 / alpha)) - 1
    for label in classes:
        here = nonconformity[calibration_labels.to_numpy() == label]
        if per_class and len(here) >= need:
            thresholds[label] = N._conformal_quantile(here, alpha)
        elif per_class:
            fallback.append(label)
    sets = [[label for j, label in enumerate(classes) if np.isfinite(row[j]) and 1.0-row[j] <= thresholds[label]]
        for row in test_scores.to_numpy(dtype=float)]
    calls = [members[0] if len(members) == 1 else None for members in sets]
    rows = pd.DataFrame({'entity': list(split.entities('test')), 'prediction': pd.Series(calls, dtype=object),
        'support': [None] * len(sets), 'abstained': [len(members) != 1 for members in sets],
        'calibrated_confidence': [None] * len(sets), 'prediction_set': sets,
        'native_set_text': [' | '.join(members) for members in sets], 'set_size': [len(members) for members in sets]})
    encode = lambda value: {'status': 'finite' if np.isfinite(value) else 'unbounded',
        'value': float(value) if np.isfinite(value) else None}
    calibration = {'strategy': 'conformal_calls', 'split_identity': split.identity, 'benchmark_id': split.benchmark_id,
        'base_model_identity': base_model.identity, 'fit_entities': list(spec.fit_entities), 'fit_role': 'train',
        'calibration_entities': list(calibration_labels.index), 'calibration_labels': calibration_labels.tolist(),
        'calibration_nonconformity': nonconformity.tolist(), 'alpha': alpha, 'per_class': per_class,
        'native_classes': classes, 'overall_threshold': encode(overall),
        'thresholds': {label: encode(value) for label, value in thresholds.items()},
        'rare_class_overall_fallback': fallback, 'minimum_class_calibration_size': need,
        'unsupported_calibration_classes': sorted(set(calibration_labels)-set(classes))}
    gaps = ['Set membership and raw base shares are not calibrated per-gene correctness probabilities',
        'Empirical coverage requires independent biological truth and exchangeability review',
        'Rare classes using overall fallback have no demonstrated class-specific guarantee']
    if calibration['unsupported_calibration_classes']:
        gaps.append('Calibration classes absent from the fitted model retain nonconformity one')
    if not classes:
        gaps.append('No native class scores; every test entity retained with an empty set')
    return LabelBatch(rows, test_scores.copy(), calibration, tuple(gaps))
