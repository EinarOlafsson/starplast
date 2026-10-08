"""Exact native forest statistics, training access guards and serialized tree replay."""
from dataclasses import replace
import json

import numpy as np
import pandas as pd
import pytest

from starplast import forest_records as F, strategies as S
from starplast.splits import Assignment, ExclusionManifest, SplitManifest


def fixture():
    ids = tuple(f'PF3D7_fixture_{i}' for i in range(12))
    split = SplitManifest('Pf', 'synthetic_forest', 'entity', tuple(Assignment(entity, entity, role)
        for entity, role in zip(ids, ['train'] * 6 + ['tune'] * 2 + ['calibration'] * 2 + ['test'] * 2)), 20261008)
    exclusions = ExclusionManifest(split.benchmark_id, ('target',), ('target', 'ec_number'), (), ids[:6], split.identity)
    features = pd.DataFrame({'x': [1., 1., 3., None, 7., 5., 99., 100., 1., 3., 2., None],
        'y': [None, 2., 2., 4., 8., 1., 100., 200., 1., 3., 3., 9.],
        'empty_train': [None] * 6 + [1., 2., 3., 4., 5., 6.]}, index=ids)
    labels = pd.Series(['["1"]', '["1"]', '["1","2"]', '["1"]', '["1","2"]', '["1","2"]'], index=ids[:6])
    return features, labels, split, exclusions


def test_catalogue_defaults_and_serial_native_constructor_statistics():
    defaults = {parameter.name: parameter.default for parameter in S.get('random_forest').params}
    assert defaults['trees'] == 300 and defaults['min_leaf'] == 2
    features, labels, split, exclusions = fixture()
    batch = F.random_forest(features, labels, split=split, exclusions=exclusions, trees=11)
    state = batch.model_state
    assert state['native_n_jobs'] == 4 and state['execution_n_jobs'] == 1
    assert state['max_features'] == 'sqrt' and state['class_weight'] == 'balanced_subsample'
    assert state['withheld_all_missing_training_columns'] == ['empty_train']
    assert state['native_classes'] == ['["1","2"]', '["1"]']
    assert len(state['trees']) == 11 and batch.rows.calibrated_confidence.isna().all()
    assert state['fit_entities'] == list(split.entities('train'))


def test_exact_native_calls_support_tree_json_and_refit():
    features, labels, split, exclusions = fixture()
    batch = F.random_forest(features, labels, split=split, exclusions=exclusions, trees=11)
    state = json.loads(json.dumps(batch.model_state, allow_nan=False))
    restored = F.replay_scores(features, state)
    pd.testing.assert_frame_equal(batch.class_scores, restored, check_exact=True)
    native_matrix = F._matrix(features, state['training_distributions'], state['kept_columns'])
    visible = pd.Series([None] * len(features), dtype=object)
    visible.iloc[:6] = labels.to_numpy()
    calls, support, model = F.native_serial_forest(native_matrix, visible, [10, 11], trees=11, min_leaf=2, seed=split.seed)
    assert model.n_jobs == 1 and model.n_estimators == 11 and model.min_samples_leaf == 2
    assert batch.rows.prediction.tolist() == calls.iloc[[10, 11]].tolist()
    np.testing.assert_array_equal(batch.rows.support.to_numpy(), support.iloc[[10, 11]].to_numpy())
    np.testing.assert_array_equal(batch.class_scores.to_numpy(), S.class_scores_of(calls).iloc[[10, 11]].to_numpy())
    repeated = F.random_forest(features, labels, split=split, exclusions=exclusions, trees=11)
    assert state == repeated.model_state
    pd.testing.assert_frame_equal(batch.rows, repeated.rows, check_exact=True)


def test_heldout_features_do_not_select_training_transform_or_tree_state():
    features, labels, split, exclusions = fixture()
    first = F.random_forest(features, labels, split=split, exclusions=exclusions, trees=11)
    features.loc[list(split.entities('test'))] = -10000.0
    second = F.random_forest(features, labels, split=split, exclusions=exclusions, trees=11)
    for key in ('training_vectors', 'training_distributions', 'training_labels', 'kept_columns', 'trees'):
        assert first.model_state[key] == second.model_state[key]


def test_source_and_hidden_role_refusals_are_explicit():
    features, labels, split, exclusions = fixture()
    with pytest.raises(ValueError, match='contamination'):
        F.random_forest(features.assign(ec_number=1.), labels, split=split, exclusions=exclusions)
    with pytest.raises(ValueError, match='training labels only'):
        F.random_forest(features, pd.concat([labels, pd.Series({'PF3D7_fixture_10': '["7"]'})]), split=split, exclusions=exclusions)
    with pytest.raises(ValueError, match='feature population'):
        F.random_forest(features.iloc[::-1], labels, split=split, exclusions=exclusions)
    with pytest.raises(ValueError, match='exclusion scope'):
        F.random_forest(features, labels, split=split, exclusions=replace(exclusions, split_identity='a' * 64))


@pytest.mark.parametrize('kind', ['no_features', 'one_class'])
def test_unavailable_native_forest_retains_every_test_entity(kind):
    features, labels, split, exclusions = fixture()
    if kind == 'no_features':
        features = features[['empty_train']]
    else:
        labels[:] = '["1"]'
    batch = F.random_forest(features, labels, split=split, exclusions=exclusions)
    assert batch.rows.entity.tolist() == list(split.entities('test'))
    assert batch.rows.abstained.all() and batch.rows.support.isna().all()
    assert not batch.model_state['trees'] and not len(batch.class_scores.columns)
    pd.testing.assert_frame_equal(batch.class_scores, F.replay_scores(features, batch.model_state), check_exact=True)


def test_serial_constructor_is_restored_after_success_and_failure(monkeypatch):
    import sklearn.ensemble as ensemble
    original = ensemble.RandomForestClassifier
    features, labels, split, exclusions = fixture()
    F.random_forest(features, labels, split=split, exclusions=exclusions, trees=3)
    assert ensemble.RandomForestClassifier is original

    def failing(*args, **kwargs):
        raise RuntimeError('sentinel native failure')

    monkeypatch.setattr(F.SL, '_forest', failing)
    with pytest.raises(RuntimeError, match='sentinel'):
        F.native_serial_forest(np.ones((2, 1)), pd.Series(['A', 'B']), [0])
    assert ensemble.RandomForestClassifier is original


@pytest.mark.parametrize('parameter', [{'trees': True}, {'trees': 0}, {'min_leaf': 0}, {'seed': True}])
def test_invalid_settings_are_refused(parameter):
    features, labels, split, exclusions = fixture()
    with pytest.raises(ValueError, match='integer forest settings'):
        F.random_forest(features, labels, split=split, exclusions=exclusions, **parameter)
