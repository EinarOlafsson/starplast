"""Native conformal parity, rare classes and calibration/model boundaries."""
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from starplast import artifacts as A, capabilities as C, conformal_records as F, organisms as O, strategies as S, strategy_learning as N, strategy_catalog as L, scorecard as SC
from starplast.query import Query
from starplast.splits import Assignment, SplitManifest

IDS = tuple(f'TGME49_{i:06d}' for i in range(100001, 100011))
ROLES = ['train'] * 4 + ['tune'] + ['calibration'] * 3 + ['test'] * 2
SPLIT = SplitManifest(O.TOXOPLASMA, 'conformal-fixture', 'entity', tuple(Assignment(e, e, r) for e, r in zip(IDS, ROLES)), 17)
X = np.array([0., .1, .8, 1., .5, .2, .7, .9, .3, .95])[:, None]
LABELS = pd.Series(['A', 'A', 'B', 'B', None, 'A', 'B', 'C', None, None], dtype=object)


def model(tmp_path, base='kNN'):
    strategy = 'feature_knn' if base == 'kNN' else 'supervised_classifier'
    digest = 'a' * 64
    deps = tuple(A.Dependency(kind, name, value) for kind, name, value in (
        ('code', 'fixture', digest), ('table', 'entities', digest), ('exclusions', 'fixture', digest), ('split', 'nested', SPLIT.identity)))
    scope = A.canonical_object({'unit': 'gene', 'eligible_population': 4, 'truth_grade': 'synthetic_control',
        'negative_semantics': 'fixture only', 'context': {}, 'limitations': ['not biological']})
    spec = A.ArtifactSpec(Query(O.TOXOPLASMA, 'label', target='location').to_json(), strategy, SC.T_LABEL,
        C.get(strategy).outputs[0], 'location', SPLIT.entities('train'), 'fitted_model', 'train:'+SPLIT.identity,
        deps, A.canonical_object({'k': 3}), scope, 17, 'fixture', fit_entities=SPLIT.entities('train'), fit_role='train',
        benchmark_id=SPLIT.benchmark_id, confidence_kind='method_support')
    A.write_artifact(tmp_path/'base', spec, {'model_state.json': {'fit_entities': list(SPLIT.entities('train')), 'native_classes': ['A', 'B']}}, split=SPLIT)
    return A.read_artifact(tmp_path/'base', expected=spec, split=SPLIT)


def scores(base='kNN'):
    fit = LABELS.copy()
    fit.iloc[5:8] = None
    if base == 'kNN':
        predictions, _ = S.knn_vote(X, fit, 3, query=[5, 6, 7, 8, 9])
    else:
        predictions, _, _ = L._logistic(X, fit, [5, 6, 7, 8, 9], C=.5)
    raw = S.class_scores_of(predictions)
    cal, test = raw.iloc[5:8].copy(), raw.iloc[8:].copy()
    cal.index, test.index = list(SPLIT.entities('calibration')), list(SPLIT.entities('test'))
    labels = pd.Series(LABELS.iloc[5:8].to_numpy(), index=cal.index)
    return cal, test, labels


@pytest.mark.parametrize('alpha', [.1, .5])
@pytest.mark.parametrize('per_class', [True, False])
@pytest.mark.parametrize('base', ['kNN', 'logistic'])
def test_native_predictions_scores_sets_and_thresholds_on_identical_partitions(tmp_path, monkeypatch, alpha, per_class, base):
    cal, test, labels = scores(base)
    batch = F.conformal_calls(cal, test, labels, split=SPLIT, base_model=model(tmp_path, base), alpha=alpha, per_class=per_class)
    # Preserve the released algorithm; only its partition chooser is replaced
    # for parity against the already-frozen train/calibration/test roles.
    monkeypatch.setattr(N, '_group_split', lambda *args: (np.arange(4), np.arange(5, 8)))
    ctx = S.Context(pd.DataFrame({'gene_id': IDS}), graph={}, organism=O.TOXOPLASMA)
    p = {'model': base, 'k': 3, 'C': .5, 'alpha': alpha, 'thresholds': 'per class' if per_class else 'overall'}
    calls, sizes, sets, thresholds = N._conformal(ctx, X, LABELS, [8, 9], p)
    assert batch.rows.prediction.tolist() == calls.iloc[8:].astype(object).where(calls.iloc[8:].notna(), None).tolist()
    assert batch.rows.set_size.tolist() == sizes.iloc[8:].tolist()
    assert batch.rows.native_set_text.tolist() == sets.iloc[8:].tolist()
    np.testing.assert_array_equal(batch.class_scores, S.class_scores_of(calls).iloc[8:])
    for label, threshold in thresholds.items():
        assert batch.model_state['thresholds'][label] == {'status': 'finite' if np.isfinite(threshold) else 'unbounded',
            'value': threshold if np.isfinite(threshold) else None}
    assert batch.model_state['unsupported_calibration_classes'] == ['C']
    assert batch.model_state['calibration_nonconformity'][-1] == 1.


def test_calibration_truth_cannot_be_training_or_outer_truth(tmp_path):
    cal, test, labels = scores()
    base = model(tmp_path)
    with pytest.raises(ValueError, match='calibration-label cohort only'):
        F.conformal_calls(cal, test, pd.Series(['A'] * 2, index=test.index), split=SPLIT, base_model=base)
    with pytest.raises(ValueError, match='fitted-model artifact'):
        F.conformal_calls(cal, test, labels, split=SPLIT, base_model=replace(base, spec=replace(base.spec, role='reusable_base')))
    with pytest.raises(ValueError, match='class orders differ'):
        F.conformal_calls(cal, test[['B', 'A']], labels, split=SPLIT, base_model=base)
    with pytest.raises(ValueError, match='population/order'):
        F.conformal_calls(cal.iloc[::-1], test, labels, split=SPLIT, base_model=base)


def test_unavailable_scores_retain_empty_sets_and_rare_class_fallback(tmp_path):
    cal, test, labels = scores()
    test.iloc[:] = np.nan
    batch = F.conformal_calls(cal, test, labels, split=SPLIT, base_model=model(tmp_path))
    assert batch.rows.prediction_set.tolist() == [[], []]
    assert batch.rows.abstained.all() and batch.rows.calibrated_confidence.isna().all()
    assert batch.model_state['rare_class_overall_fallback'] == ['A', 'B']
    assert batch.model_state['overall_threshold']['status'] == 'unbounded'


def test_nonfinite_or_partial_or_unnormalized_scores_are_refused(tmp_path):
    cal, test, labels = scores()
    base = model(tmp_path)
    for invalid, message in (([np.inf, .2], 'finite'), ([np.nan, .2], 'complete normalized'), ([.2, .2], 'complete normalized')):
        broken = test.copy()
        broken.iloc[0] = invalid
        with pytest.raises(ValueError, match=message):
            F.conformal_calls(cal, broken, labels, split=SPLIT, base_model=base)


def test_mutated_model_vocabulary_cannot_reuse_verified_identity(tmp_path):
    cal, test, labels = scores()
    base = model(tmp_path)
    base.payloads['model_state.json']['native_classes'] = ['B', 'A']
    with pytest.raises(ValueError, match='changed after loading'):
        F.conformal_calls(cal[['B', 'A']], test[['B', 'A']], labels, split=SPLIT, base_model=base)
