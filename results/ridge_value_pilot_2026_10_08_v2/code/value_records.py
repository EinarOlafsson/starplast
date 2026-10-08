"""Native numeric predictions under frozen training-only preprocessing and truth."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from . import strategy_catalog as N
from .splits import ExclusionManifest, SplitManifest


@dataclass
class ValueBatch:
    """Every numeric test outcome with explicit missingness and fitting provenance."""

    rows: pd.DataFrame
    model_state: dict
    gaps: tuple[str, ...]


def trait_ridge(features, training_values, *, split, exclusions):
    """Fit native ridge on training values with frozen rank scaling/zero imputation.

    Training ranks match native Context.matrix. Query ranks use only the training
    distribution, including average ranks for observed ties and the right ECDF
    for unseen values. Target values retain their original scale. This adapter
    never receives held-out truth or selects parameters against test outcomes.
    """
    if not isinstance(split, SplitManifest) or not isinstance(exclusions, ExclusionManifest):
        raise TypeError('Frozen split and training exclusions are required')
    if (exclusions.split_identity != split.identity or exclusions.benchmark_id != split.benchmark_id
            or exclusions.fit_entities != split.entities('train')):
        raise ValueError('Exclusions must match the complete ordered training cohort')
    ids = tuple(a.entity for a in split.assignments)
    if (not isinstance(features, pd.DataFrame) or features.index.has_duplicates or features.columns.has_duplicates
            or tuple(features.index) != ids):
        raise ValueError('Features must match the complete unique frozen population/order')
    if not isinstance(training_values, pd.Series) or tuple(training_values.index) != split.entities('train'):
        raise ValueError('Supply complete ordered training values only')
    y = pd.to_numeric(training_values, errors='raise').astype(float)
    if not np.isfinite(y).all():
        raise ValueError('Training values must be explicitly observed finite numbers')
    exclusions.guard_inputs(columns=features.columns, benchmark_id=split.benchmark_id)
    for stage in ('feature_selection', 'imputation', 'scaling', 'model_fit'):
        split.guard_fit(stage, y.index)
    numeric = features.apply(pd.to_numeric, errors='raise').astype(float)
    if np.isinf(numeric.to_numpy()).any():
        raise ValueError('Infinite features require an explicit preprocessing decision')
    train = numeric.loc[list(y.index)]
    kept = train.columns[train.notna().any()].tolist()
    withheld = train.columns[~train.notna().any()].tolist()
    distributions = {c: np.sort(train[c].dropna().to_numpy()) for c in kept}
    matrix = np.zeros((len(ids), len(kept)))
    for j, column in enumerate(kept):
        values = numeric[column].to_numpy()
        known = np.isfinite(values)
        ordered = distributions[column]
        left = np.searchsorted(ordered, values[known], side='left')
        right = np.searchsorted(ordered, values[known], side='right')
        matrix[known, j] = np.where(right > left, (left+right+1)/2, right)/len(ordered)-.5
    where = {e: i for i, e in enumerate(ids)}
    training_positions = [where[e] for e in y.index]
    test_positions = [where[e] for e in split.entities('test')]
    model = None
    predictions = np.full(len(test_positions), np.nan)
    if len(y) >= 20 and kept:
        model = N._regressor('ridge').fit(matrix[training_positions], y.to_numpy())
        predictions = model.predict(matrix[test_positions])
        if not np.isfinite(predictions).all():
            raise ValueError('Native ridge produced a nonfinite supported prediction')
    rows = pd.DataFrame({'entity': list(split.entities('test')), 'prediction': predictions,
        'abstained': ~np.isfinite(predictions), 'calibrated_confidence': [None]*len(predictions)})
    state = {'strategy': 'trait_regression', 'model': 'ridge', 'alpha': 1., 'fit_role': 'train',
        'fit_entities': y.index.tolist(), 'training_values': y.tolist(), 'split_identity': split.identity,
        'benchmark_id': split.benchmark_id, 'input_columns': features.columns.tolist(), 'kept_columns': kept,
        'withheld_all_missing_training_columns': withheld,
        'feature_transform': 'training_average_rank_percentile_centered', 'missing_feature_value': 0.,
        'unseen_feature_transform': 'right_training_ecdf_centered',
        'training_distributions': {c: a.tolist() for c, a in distributions.items()},
        'training_vectors': matrix[training_positions].tolist(),
        'coefficients': model.coef_.tolist() if model is not None else None,
        'intercept': float(model.intercept_) if model is not None else None,
        'training_target_range': [float(y.min()), float(y.max())],
        'source_exclusion_columns': list(exclusions.columns), 'source_exclusion_layers': list(exclusions.layers)}
    gaps = ['Independent biology, original measurement units and applicability require source-specific admission',
            'Training-fitted native ranks replace whole-context ranks for inductive evaluation',
            'Single fixed ridge variant; other numeric mechanisms and full outer coverage remain open']
    if model is None:
        gaps.append('Native support requires at least 20 training values and a usable feature; all test rows abstain')
    return ValueBatch(rows, state, tuple(gaps))
