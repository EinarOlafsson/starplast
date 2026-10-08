"""Frozen categorical adapters retaining native predictions and class scores.

The first adapter uses the existing kNN vote on training-fitted numeric inputs.
It never receives hidden truth, filters rare evaluation classes, tunes against
outer outcomes or labels a vote share as calibrated confidence.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from . import strategies as S
from .splits import ExclusionManifest, SplitManifest


@dataclass
class LabelBatch:
    """Native held-out calls, class-score matrix and data-only fitting provenance."""

    rows: pd.DataFrame
    class_scores: pd.DataFrame
    model_state: dict
    gaps: tuple[str, ...]


def feature_knn(features, training_labels, *, split, exclusions, k=15, min_share=0.3):
    """Predict the complete frozen test cohort using only explicitly supplied train truth.

    ``features`` is an ordered numeric dataframe for all split entities, already
    limited to permitted source columns. Unknown feature values are imputed by
    training medians; all-missing training columns are explicitly withheld.
    The native vote operates on the resulting matrix with its original k/tie
    behavior, and its raw class scores remain available after abstention.
    """
    if not isinstance(split, SplitManifest) or not isinstance(exclusions, ExclusionManifest):
        raise TypeError('Frozen split and training-derived exclusions are required')
    if exclusions.split_identity != split.identity or exclusions.benchmark_id != split.benchmark_id:
        raise ValueError('Source exclusions and split identities differ')
    if exclusions.fit_entities != split.entities('train'):
        raise ValueError('Exclusions must be frozen from the complete ordered training cohort')
    if not isinstance(features, pd.DataFrame) or features.index.has_duplicates or features.columns.has_duplicates:
        raise ValueError('Ordered unique feature rows and columns are required')
    ids = tuple(row.entity for row in split.assignments)
    if tuple(features.index) != ids:
        raise ValueError('Feature population/order differs from the frozen split')
    if not isinstance(training_labels, pd.Series) or tuple(training_labels.index) != split.entities('train'):
        raise ValueError('Provide the complete ordered training-label cohort only')
    if not all(isinstance(v, str) and v for v in training_labels):
        raise ValueError('Training labels must be observed categorical values')
    if type(k) is not int or k < 3 or not isinstance(min_share, (int, float)) or isinstance(min_share, bool) or not 0 <= min_share <= 1:
        raise ValueError('Invalid frozen kNN settings')
    exclusions.guard_inputs(columns=features.columns, benchmark_id=split.benchmark_id)
    for stage in ('feature_selection', 'imputation', 'scaling', 'model_fit'):
        split.guard_fit(stage, training_labels.index)
    numeric = features.apply(pd.to_numeric, errors='raise').astype(float)
    if np.isinf(numeric.to_numpy()).any():
        raise ValueError('Infinite features require an explicit preprocessing decision')
    train = numeric.loc[list(training_labels.index)]
    medians = train.median()
    kept = medians.index[medians.notna()].tolist()
    withheld = medians.index[medians.isna()].tolist()
    filled = numeric[kept].fillna(medians[kept])
    train_filled = filled.loc[list(training_labels.index)]
    means = train_filled.mean()
    scales = train_filled.std(ddof=0).replace(0, 1)
    X = ((filled - means) / scales).to_numpy(dtype=float)
    where = {entity: i for i, entity in enumerate(ids)}
    visible = pd.Series([None] * len(ids), dtype=object)
    for entity, label in training_labels.items():
        visible.iloc[where[entity]] = label
    query = np.asarray([where[entity] for entity in split.entities('test')], dtype=int)
    native, shares = S.knn_vote(X, visible, k, query=query)
    called = native.where(shares >= min_share)
    native_scores = S.class_scores_of(native)
    score_frame = native_scores.iloc[query].copy() if native_scores is not None else pd.DataFrame(index=range(len(query)))
    score_frame.index = list(split.entities('test'))
    predictions = called.iloc[query].astype(object).where(called.iloc[query].notna(), None).tolist()
    support = shares.iloc[query].tolist()
    rows = pd.DataFrame({'entity': list(split.entities('test')), 'prediction': predictions,
                         'support': [float(v) if np.isfinite(v) else None for v in support],
                         'abstained': [not isinstance(v, str) for v in predictions],
                         'calibrated_confidence': [None] * len(query)})
    model = {'strategy': 'feature_knn', 'split_identity': split.identity, 'benchmark_id': split.benchmark_id,
        'fit_entities': list(training_labels.index), 'fit_role': 'train', 'k': k, 'min_share': min_share,
        'input_columns': features.columns.tolist(), 'kept_columns': kept, 'withheld_all_missing_training_columns': withheld,
        'training_medians': medians[kept].to_dict(), 'training_means': means.to_dict(), 'training_scales': scales.to_dict(),
        'training_class_counts': training_labels.value_counts().to_dict(),
        'training_vectors': X[[where[e] for e in training_labels.index]].tolist(),
        'training_labels': training_labels.tolist(), 'native_classes': score_frame.columns.tolist(),
        'source_exclusion_columns': list(exclusions.columns), 'source_exclusion_layers': list(exclusions.layers)}
    gaps = ['Vote shares and class scores are method support, not calibrated probabilities',
            'Biological accuracy requires a separately admitted independent truth source']
    if withheld:
        gaps.append('All-missing training feature columns withheld explicitly')
    if not kept or len(training_labels) < 2:
        gaps.append('No supported native vote; every test entity retained as an abstention')
    return LabelBatch(rows, score_frame, model, tuple(gaps))
