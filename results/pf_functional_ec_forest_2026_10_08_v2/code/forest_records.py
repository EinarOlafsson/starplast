"""Frozen native categorical forests with exact serial tree-state replay.

This opt-in record adapter supplies the existing native forest with only ordered
training labels and training-fitted numeric ranks. It changes the constructor's
job count from four to one inside a restored local wrapper so tree probability
addition follows a reproducible order. Statistical defaults remain native.
Serializable trees preserve exact fitted state and missing predictions remain
abstentions. Scores are method support; biological accuracy, calibrated gene
confidence and deployment applicability require independent evidence.
"""
from __future__ import annotations

import json
from unittest.mock import patch

import numpy as np
import pandas as pd

from . import strategies as S, strategy_learning as SL
from .label_records import LabelBatch
from .splits import ExclusionManifest, SplitManifest


def _plain(value):
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_plain(item) for item in value]
    return value.item() if isinstance(value, np.generic) else value


def native_serial_forest(matrix, visible, query, *, trees=300, min_leaf=2, seed=0):
    """Call unchanged native forest statistics with a scoped, restored one-job constructor."""
    from sklearn.ensemble import RandomForestClassifier

    def constructor(*args, **kwargs):
        kwargs['n_jobs'] = 1
        return RandomForestClassifier(*args, **kwargs)

    with patch('sklearn.ensemble.RandomForestClassifier', constructor):
        return SL._forest(matrix, visible, query, trees, min_leaf, seed)


def _matrix(features, distributions, columns):
    matrix = np.zeros((len(features), len(columns)), dtype=float)
    for position, column in enumerate(columns):
        ordered = np.asarray(distributions[column], dtype=float)
        values = features[column].to_numpy(dtype=float)
        known = np.isfinite(values)
        left = np.searchsorted(ordered, values[known], side='left')
        right = np.searchsorted(ordered, values[known], side='right')
        matrix[known, position] = np.where(right > left, (left + right + 1) / 2, right) / len(ordered) - 0.5
    return matrix


def replay_scores(features, state):
    """Restore native Trees from finite JSON state and replay exact serial class scores."""
    from sklearn.tree import _tree
    if tuple(features.index) != tuple(state['entity_order']) or list(features.columns) != state['input_columns']:
        raise ValueError('Tree replay requires exact original feature population/order')
    matrix = _matrix(features, state['training_distributions'], state['kept_columns'])
    where = {entity: position for position, entity in enumerate(features.index)}
    query = [where[entity] for entity in state['query_entities']]
    classes = state['native_classes']
    if not state['trees']:
        return pd.DataFrame(index=state['query_entities']) if not classes else pd.DataFrame(
            index=state['query_entities'], columns=classes, dtype=float)
    values = np.asarray(matrix[query], dtype=np.float32, order='C')
    scores = np.zeros((len(query), len(classes)), dtype=np.float64)
    for saved in state['trees']:
        tree = _tree.Tree(saved['n_features'], np.asarray(saved['n_classes'], dtype=np.intp), saved['n_outputs'])
        dtype = np.dtype(saved['node_dtype_schema'])
        nodes = np.asarray([tuple(row) for row in saved['nodes']], dtype=dtype)
        tree.__setstate__({'max_depth': saved['max_depth'], 'node_count': saved['node_count'], 'nodes': nodes,
                           'values': np.asarray(saved['values'], dtype=np.float64)})
        scores += tree.predict(values)[:, :len(classes)]
    scores /= len(state['trees'])
    return pd.DataFrame(scores, index=state['query_entities'], columns=classes)


def random_forest(features, training_labels, *, split, exclusions, trees=300, min_leaf=2, seed=None):
    """Fit fixed native forest statistics on training-only ranks and retain every test entity."""
    import sklearn
    if not isinstance(split, SplitManifest) or not isinstance(exclusions, ExclusionManifest):
        raise TypeError('Typed frozen split and exclusions required')
    if (exclusions.split_identity != split.identity or exclusions.benchmark_id != split.benchmark_id or
            exclusions.fit_entities != split.entities('train')):
        raise ValueError('Training exclusion scope differs from the exact split')
    if (not isinstance(features, pd.DataFrame) or features.index.has_duplicates or features.columns.has_duplicates or
            tuple(features.index) != tuple(row.entity for row in split.assignments)):
        raise ValueError('Exact ordered unique feature population required')
    if (not isinstance(training_labels, pd.Series) or tuple(training_labels.index) != split.entities('train') or
            not all(isinstance(label, str) and label for label in training_labels)):
        raise ValueError('Exact complete ordered training labels only')
    seed = split.seed if seed is None else seed
    if type(trees) is not int or trees < 1 or type(min_leaf) is not int or min_leaf < 1 or type(seed) is not int:
        raise ValueError('Explicit integer forest settings required')
    exclusions.guard_inputs(columns=features.columns, benchmark_id=split.benchmark_id)
    for stage in ('feature_selection', 'imputation', 'scaling', 'model_fit'):
        split.guard_fit(stage, training_labels.index)
    numeric = features.apply(pd.to_numeric, errors='raise').astype(float)
    if np.isinf(numeric.to_numpy()).any():
        raise ValueError('Infinite inputs require an explicit preprocessing decision')
    train = numeric.loc[list(training_labels.index)]
    kept = train.columns[train.notna().any()].tolist()
    withheld = train.columns[~train.notna().any()].tolist()
    distributions = {column: np.sort(train[column].dropna().to_numpy()).tolist() for column in kept}
    matrix = _matrix(numeric, distributions, kept)
    where = {entity: position for position, entity in enumerate(numeric.index)}
    training = [where[entity] for entity in training_labels.index]
    query = np.asarray([where[entity] for entity in split.entities('test')], dtype=int)
    expected = (train[kept].rank(pct=True) - 0.5).fillna(0).to_numpy()
    np.testing.assert_array_equal(matrix[training], expected)
    visible = pd.Series([None] * len(numeric), dtype=object)
    visible.iloc[training] = training_labels.to_numpy()
    calls, support, model = native_serial_forest(matrix, visible, query, trees=trees, min_leaf=min_leaf, seed=seed)
    native_scores = S.class_scores_of(calls)
    scores = native_scores.iloc[query].copy() if native_scores is not None else pd.DataFrame(index=range(len(query)))
    scores.index = list(split.entities('test'))
    saved_trees = []
    for estimator in model.estimators_ if model is not None else ():
        tree = estimator.tree_
        data = tree.__getstate__()
        dtype = data['nodes'].dtype
        saved_trees.append(_plain({'n_features': tree.n_features, 'n_classes': tree.n_classes, 'n_outputs': tree.n_outputs,
            'max_depth': data['max_depth'], 'node_count': data['node_count'], 'node_dtype': data['nodes'].dtype.descr,
            'node_dtype_schema': {'names': list(dtype.names), 'formats': [dtype.fields[name][0].str for name in dtype.names],
                                  'offsets': [dtype.fields[name][1] for name in dtype.names], 'itemsize': dtype.itemsize},
            'nodes': data['nodes'].tolist(), 'values': data['values'], 'random_state': estimator.random_state}))
    state = _plain({'strategy': 'random_forest', 'benchmark_id': split.benchmark_id, 'split_identity': split.identity,
        'fit_entities': list(training_labels.index), 'fit_role': 'train', 'seed': seed, 'n_estimators': trees,
        'min_samples_leaf': min_leaf, 'max_features': 'sqrt', 'class_weight': 'balanced_subsample',
        'native_n_jobs': 4, 'execution_n_jobs': 1, 'execution_difference': 'Serial tree probability reduction; statistical parameters unchanged',
        'entity_order': list(numeric.index), 'query_entities': list(split.entities('test')), 'input_columns': list(features.columns),
        'kept_columns': kept, 'withheld_all_missing_training_columns': withheld, 'training_distributions': distributions,
        'feature_transform': 'training_average_rank_percentile_centered', 'missing_feature_value': 0.0,
        'unseen_feature_transform': 'right_training_ecdf_centered', 'training_vectors': matrix[training],
        'training_labels': training_labels.tolist(), 'training_class_counts': training_labels.value_counts().to_dict(),
        'native_classes': list(scores.columns), 'trees': saved_trees, 'sklearn_version': sklearn.__version__,
        'source_exclusion_columns': list(exclusions.columns), 'source_exclusion_layers': list(exclusions.layers)})
    restored = json.loads(json.dumps(state, allow_nan=False))
    assert restored == state
    pd.testing.assert_frame_equal(scores, replay_scores(features, restored), check_exact=True)
    predictions = calls.iloc[query].astype(object).where(calls.iloc[query].notna(), None).tolist()
    values = support.iloc[query].to_numpy(dtype=float)
    if model is not None:
        assert predictions == np.asarray(scores.columns)[scores.to_numpy().argmax(axis=1)].tolist()
        np.testing.assert_array_equal(values, scores.to_numpy().max(axis=1))
    rows = pd.DataFrame({'entity': list(split.entities('test')), 'prediction': predictions,
        'support': [float(value) if np.isfinite(value) else None for value in values],
        'abstained': [not isinstance(value, str) for value in predictions], 'calibrated_confidence': [None] * len(query)})
    gaps = ['Native forest class scores are method support, not calibrated gene confidence',
            'Native job count four replaced with one solely for exact serial probability reduction',
            'Training-only ranks replace whole-context ranking for inductive evaluation',
            'Independent biological accuracy and deployment require separately admitted evidence']
    if model is None:
        gaps.append('Unavailable native forest: fewer than two training classes or no supported numeric features')
    return LabelBatch(rows, scores, state, tuple(gaps))
