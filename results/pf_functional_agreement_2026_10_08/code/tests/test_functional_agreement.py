"""Paired reporting preserves abstentions and refuses incompatible frozen scopes."""
from copy import deepcopy
import json

import pytest

from starplast import functional_agreement as A, functional_results as F
from starplast.scorecard_view import export_scorecard, render_scorecard_html


def _rows():
    calls = [(None, None), ('a', None), ('b', None), (None, 'a'), (None, 'b'),
             ('a', 'a'), ('b', 'b'), ('a', 'b'), ('b', 'a'), ('b', 'c')]
    def row(i, prediction):
        return dict(entity=str(i), truth='a', group=str(i // 2),
                    training_supported=i != 9, prediction=prediction,
                    abstained=prediction is None, calibrated_confidence=None)
    return ([row(i, x) for i, (x, y) in enumerate(calls)],
            [row(i, y) for i, (x, y) in enumerate(calls)])


def test_every_outcome_is_distinct_and_abstentions_are_not_wrong_calls():
    paired = A.paired_rows(*_rows())
    assert [row['outcome'] for row in paired] == list(A._OUTCOMES)
    summary = A._summary(paired)
    assert all(summary['counts'][key] == 1 for key in A._OUTCOMES)
    assert summary['counts']['eligible'] == 10
    assert summary['counts']['joint_answered'] == 5
    assert summary['counts']['shared_wrong'] == 2
    assert summary['counts']['unsupported'] == 1
    assert summary['rates']['source_recovery_among_agreements'] == dict(numerator=1, denominator=2, value=0.5)
    assert summary['rates']['shared_wrong_among_joint_calls'] == dict(numerator=2, denominator=5, value=2 / 5)


def test_zero_joint_answers_and_zero_agreements_are_unavailable():
    left, right = _rows()
    report = A._summary(A.paired_rows(left[:1], right[:1]))
    assert report['rates']['joint_coverage']['value'] == 0
    assert report['rates']['source_recovery_among_agreements']['value'] is None
    assert report['rates']['shared_wrong_among_joint_calls']['value'] is None
    report = A._summary(A.paired_rows(left[7:8], right[7:8]))
    assert report['rates']['agreement_among_joint_calls']['value'] == 0
    assert report['rates']['source_recovery_among_agreements']['value'] is None


@pytest.mark.parametrize('key,value', [('entity', 'another'), ('truth', 'unknown'),
    ('group', 'another'), ('training_supported', False), ('abstained', False),
    ('calibrated_confidence', 0.9), ('prediction', 2)])
def test_row_scope_and_call_semantics_are_not_silently_repaired(key, value):
    left, right = _rows()
    right[0][key] = value
    with pytest.raises(ValueError): A.paired_rows(left, right)


def test_missing_duplicate_and_reordered_population_are_refused():
    left, right = _rows()
    for candidate in ([], right[:-1], list(reversed(right))):
        with pytest.raises(ValueError): A.paired_rows(left, candidate)
    left[1] = deepcopy(left[0]); right[1] = deepcopy(right[0])
    with pytest.raises(ValueError): A.paired_rows(left, right)


@pytest.fixture(scope='module')
def pair():
    benchmarks, reason = F.shipped('Pf')
    assert not reason
    return tuple(next(b for b in benchmarks if b.metadata['strategy'] == strategy)
                 for strategy in ('feature_knn', 'random_forest'))


def _repair(benchmark):
    manifest = benchmark.metadata['source_manifest']
    manifest['key'] = F._hash(manifest['spec'])
    manifest['identity'] = F._hash({key: value for key, value in manifest.items() if key != 'identity'})
    benchmark.metadata['source_artifact_identity'] = manifest['identity']
    for name in ('summary', 'source_summary'):
        benchmark.metadata[name]['artifact_identity'] = manifest['identity']
    return benchmark


def test_shipped_pair_retains_all_rows_original_methods_and_exact_export(pair):
    original = deepcopy(pair)
    report = A.compare(*pair)
    assert report['counts']['eligible'] == 152
    assert report['counts']['unsupported'] == 4
    assert report['counts']['joint_answered'] == 143
    assert sum(report['counts'][key] for key in A._OUTCOMES) == 152
    assert report['counts']['both_abstain'] == 0
    assert report['counts']['left_only_correct'] + report['counts']['left_only_wrong'] == 0
    assert report['methods'][0]['counts']['correct'] == 40
    assert report['methods'][1]['counts']['correct'] == 59
    assert report['overlap']['same_training_entities'] is True
    assert report['overlap']['independent_evidence'] is False
    assert {'truth', 'table', 'split', 'exclusions'} <= {d['kind'] for d in report['overlap']['shared_dependencies']}
    view = A.build_scorecard(report)
    assert json.loads(export_scorecard(view))['snapshot'] == report
    assert 'Source recovery among agreements' in render_scorecard_html(view)
    report['rows'][0]['entity'] = 'changed after display'
    assert json.loads(view.snapshot_json)['rows'][0]['entity'] != report['rows'][0]['entity']
    assert pair == original


@pytest.mark.parametrize('field', ['context', 'unit', 'negative_semantics',
    'prepared_packet_manifest_sha256', 'calibration', 'deployment'])
def test_context_refusals_survive_a_consistent_artifact_envelope(pair, field):
    left, right = deepcopy(pair)
    spec = right.metadata['source_manifest']['spec']
    scope = json.loads(spec['evaluation_scope_json'])
    scope[field] = 'different'
    spec['evaluation_scope_json'] = F._canonical(scope)
    _repair(right)
    with pytest.raises(ValueError, match='context changed'): A.compare(left, right)


def test_equal_artifact_and_changed_source_context_are_refused(pair):
    with pytest.raises(ValueError, match='distinct artifacts'): A.compare(pair[0], pair[0])
    left, right = deepcopy(pair)
    right.metadata['source_context'] = 'different assay'
    with pytest.raises(ValueError, match='source scope'): A.compare(left, right)


def test_training_cohort_and_duplicate_dependencies_are_refused(pair):
    left, right = deepcopy(pair)
    spec = right.metadata['source_manifest']['spec']
    spec['fit_entities'][0] = 'different training entity'
    _repair(right)
    with pytest.raises(ValueError, match='fit_entities'): A.compare(left, right)
    left, right = deepcopy(pair)
    spec = right.metadata['source_manifest']['spec']
    spec['dependencies'].append(deepcopy(spec['dependencies'][0]))
    _repair(right)
    with pytest.raises(ValueError, match='Duplicate paired dependency'): A.compare(left, right)


@pytest.mark.parametrize('key,value', [('biological_accuracy', .9),
    ('calibrated_confidence', .9), ('biological_admission', True), ('ensemble_selection', True)])
def test_card_cannot_promote_source_recovery_or_select_a_model(pair, key, value):
    report = A.compare(*pair); report[key] = value
    with pytest.raises(ValueError, match='cannot advertise'): A.build_scorecard(report)


def test_card_rejects_changed_counts(pair):
    report = A.compare(*pair); report['counts']['agreements'] += 1
    with pytest.raises(ValueError, match='counts differ'): A.build_scorecard(report)


def test_card_rejects_recounted_but_fabricated_outcomes(pair):
    report = A.compare(*pair)
    report['rows'][0]['outcome'] = 'both_abstain'
    report.update(A._summary(report['rows']))
    with pytest.raises(ValueError, match='outcomes differ'): A.build_scorecard(report)
