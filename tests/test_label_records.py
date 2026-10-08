"""Native-vote parity and training-only feature fitting for categorical records."""
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from starplast import label_records as L, organisms as O, strategies as S
from starplast.splits import Assignment, ExclusionManifest, SplitManifest

IDS = tuple(f'TGME49_{i:06d}' for i in range(100001, 100009))
SPLIT = SplitManifest(O.TOXOPLASMA, 'fixture', 'entity', tuple(Assignment(e, e, role)
    for role, cohort in zip(('train', 'tune', 'calibration', 'test'), (IDS[:2], IDS[2:4], IDS[4:6], IDS[6:])) for e in cohort), 3)
EXCLUSION = ExclusionManifest('fixture', ('target',), ('target',), (), IDS[:2], SPLIT.identity)


def fixture():
    return pd.DataFrame({'x': [1., 3., None, 6., 2., 3., 5., None], 'empty_in_training': [None, None, 5., 2., 1., 4., 3., 5.]}, index=IDS), pd.Series(['A', 'B'], index=IDS[:2])


def test_native_prediction_and_class_score_parity():
    features, truth = fixture()
    batch = L.feature_knn(features, truth, split=SPLIT, exclusions=EXCLUSION, k=3, min_share=.3)
    X = np.array([[0.], [.5], [0.], [.5], [0.], [.5], [.5], [0.]])
    visible = pd.Series(['A', 'B', None, None, None, None, None, None])
    native, support = S.knn_vote(X, visible, 3, query=[6, 7])
    assert batch.rows.prediction.tolist() == native.iloc[[6, 7]].tolist()
    assert batch.rows.support.tolist() == support.iloc[[6, 7]].tolist()
    np.testing.assert_allclose(batch.class_scores, S.class_scores_of(native).iloc[[6, 7]])
    assert batch.model_state['withheld_all_missing_training_columns'] == ['empty_in_training']
    assert batch.rows.entity.tolist() == list(SPLIT.entities('test'))
    assert batch.rows.calibrated_confidence.isna().all()


def test_held_out_features_cannot_change_imputation_scaling_or_feature_selection():
    features, truth = fixture()
    first = L.feature_knn(features, truth, split=SPLIT, exclusions=EXCLUSION)
    features.loc[list(IDS[2:]), 'x'] = 10000.
    second = L.feature_knn(features, truth, split=SPLIT, exclusions=EXCLUSION)
    for field in ('training_distributions', 'kept_columns', 'training_vectors'):
        assert first.model_state[field] == second.model_state[field]


def test_training_rank_transform_matches_native_with_ties_and_missingness():
    ids = tuple(f'TGME49_{i:06d}' for i in range(100001, 100011))
    roles = ['train'] * 5 + ['tune', 'calibration'] + ['test'] * 3
    assignments = tuple(Assignment(e, e, role) for e, role in zip(ids, roles))
    split = replace(SPLIT, assignments=assignments)
    exclusion = replace(EXCLUSION, fit_entities=ids[:5], split_identity=split.identity)
    features = pd.DataFrame({'x': [1., 1., 3., None, 7., 99., 99., 0., 2., 8.],
                             'y': [None, 2., 2., 4., 8., 99., 99., None, 3., 100.]}, index=ids)
    labels = pd.Series(['A', 'A', 'B', 'A', 'B'], index=ids[:5])
    batch = L.feature_knn(features, labels, split=split, exclusions=exclusion)
    native = S.Context(features.iloc[:5].reset_index(drop=True), graph={}, organism=O.TOXOPLASMA).matrix(['x', 'y'])
    np.testing.assert_array_equal(batch.model_state['training_vectors'], native)
    # Outside-support/unseen values follow the declared training ECDF, with
    # missing values fixed at the native centered median, not re-ranked.
    query = np.array([[-.5, 0.], [0., 0.], [.5, .5]])
    visible = pd.Series(['A', 'A', 'B', 'A', 'B', None, None, None])
    calls, _ = S.knn_vote(np.vstack([native, query]), visible, 15, query=[5, 6, 7])
    np.testing.assert_array_equal(batch.class_scores, S.class_scores_of(calls).iloc[[5, 6, 7]])


def test_target_and_held_out_label_leaks_are_refused():
    features, truth = fixture()
    with pytest.raises(ValueError, match='training-label cohort only'):
        L.feature_knn(features, pd.Series(['A'] * 8, index=IDS), split=SPLIT, exclusions=EXCLUSION)
    with pytest.raises(ValueError, match='contamination'):
        L.feature_knn(features.assign(target=1.), truth, split=SPLIT, exclusions=EXCLUSION)
    with pytest.raises(ValueError, match='identities differ'):
        L.feature_knn(features, truth, split=SPLIT, exclusions=replace(EXCLUSION, split_identity='a' * 64))


def test_abstention_retains_native_scores_and_no_feature_population_is_not_dropped():
    features, truth = fixture()
    batch = L.feature_knn(features, truth, split=SPLIT, exclusions=EXCLUSION, min_share=1.)
    assert batch.rows.abstained.any()
    assert batch.class_scores.notna().all().all()
    no_features = L.feature_knn(features[['empty_in_training']], truth, split=SPLIT, exclusions=EXCLUSION)
    assert len(no_features.rows) == 2 and no_features.rows.abstained.all()
    assert any('every test entity' in gap for gap in no_features.gaps)


def test_feature_order_and_nonfinite_values_require_explicit_decisions():
    features, truth = fixture()
    with pytest.raises(ValueError, match='population/order'):
        L.feature_knn(features.iloc[::-1], truth, split=SPLIT, exclusions=EXCLUSION)
    features.iloc[0, 0] = np.inf
    with pytest.raises(ValueError, match='Infinite'):
        L.feature_knn(features, truth, split=SPLIT, exclusions=EXCLUSION)
