"""Synthetic class cards preserve source scope, false calls and native nulls."""
from copy import deepcopy
import json

import pytest

from starplast import functional_class_view as C, functional_results as F, scorecard as S
from starplast import organisms as O
from starplast.scorecard_view import export_scorecard, render_scorecard_detail, render_scorecard_html


def fixture():
    scope = {'task': S.T_LABEL, 'organism': 'synthetic', 'target': 'recorded profiles',
        'strategy': 'fixture', 'settings': 'fixture only', 'seed': 1, 'protocol': 'fixture split',
        'partition': 'test', 'unit': 'gene', 'truth_grade': 'synthetic_control', 'gaps': ['Fixture only']}
    rows = [
        {'entity': 'g1', 'truth': '["2"]', 'prediction': '["2"]', 'abstained': False},
        {'entity': 'g2', 'truth': '["2"]', 'prediction': None, 'abstained': True},
        {'entity': 'g3', 'truth': '["1"]', 'prediction': '["2"]', 'abstained': False},
        {'entity': 'g4', 'truth': '["1"]', 'prediction': '["1"]', 'abstained': False},
    ]
    counts = {'eligible': 4, 'answered': 3, 'abstained': 1, 'correct': 2, 'wrong': 1}
    metrics = {'accuracy': .5, 'coverage': .75, 'precision_of_calls': 2 / 3}
    card = {'scope': scope, 'counts': counts, 'metrics': metrics}
    profile = {'scope': scope, 'class': '["2"]',
        'counts': {'eligible': 2, 'answered': 1, 'abstained': 1, 'correct': 1, 'wrong': 0},
        'metrics': {'accuracy': .5, 'precision_of_calls': 1., 'coverage': .5},
        'class_metrics': {'precision': .5, 'recall': .5, 'f1': .5, 'true_positive': 1,
            'false_positive': 1, 'false_negative': 1, 'evaluation_prevalence': .5},
        'extra': {'uncertainty': {'status': 'unavailable', 'reason': 'Fixture groups unresolved'}}}
    major = {'class': '2', 'eligible': 4, 'known_positive_genes': 2,
        'true_positive': 1, 'false_positive': 1, 'false_negative': 1,
        'precision': .5, 'recall': .5, 'f1': .5, 'coverage': .75, 'reference_prevalence': .5,
        'biological_precision': None, 'biological_recall': None, 'calibrated_confidence': None}
    missing = {**major, 'class': '3', 'known_positive_genes': 0, 'true_positive': 0,
        'false_positive': 0, 'false_negative': 0, 'precision': None, 'recall': None,
        'f1': None, 'reference_prevalence': 0.}
    metadata = {'organism': 'synthetic', 'biological_admission': False, 'calibrated_confidence': None,
        'source_targets': ['synthetic labels'], 'truth_grade': 'synthetic_control',
        'source_artifact_identity': 'fixture identity', 'source_context': 'fixture context',
        'summary': {'train': 3, 'tune': 1, 'calibration': 1, 'test': 4},
        'source_manifest': {'spec': {'fit_entities': ['train1', 'train2', 'train3'],
            'dependencies': [{'kind': 'table', 'name': 'installed_nodes', 'sha256': 'a' * 64}]}}}
    return F.FunctionalBenchmark(metadata, {'rows.json': rows, 'card.json': card,
        'profile_class_cards.json': [profile], 'major_class_cards.json': [major, missing],
        'baseline_cards.json': {'majority': card}})


def metrics(view):
    return {metric.key: metric for metric in view.metrics}


def snapshot(view):
    return json.loads(export_scorecard(view))['snapshot']


def test_profile_member_accuracy_and_whole_cohort_class_precision_are_distinct():
    benchmark = fixture()
    view = C.build(benchmark, 'profile:["2"]')
    values = metrics(view)
    assert dict(view.counts) == benchmark.profile_class_cards[0]['counts']
    assert values['accuracy'].value == .5
    assert values['precision_of_calls'].value == 1.
    assert values['class_precision'].value == .5
    assert values['precision_of_calls'].denominator == 'Calls among reference-profile members: 1'
    assert values['class_precision'].denominator == 'Whole-cohort calls of this class: 2'
    assert values['accuracy'].denominator == 'Reference-profile members: 2'
    assert values['class_recall'].denominator == 'Reference members in full cohort: 2'
    assert values['reference_prevalence'].denominator == 'All frozen test genes: 4'
    rows = snapshot(view)['details']['rows']
    assert [row['entity'] for row in rows] == ['g1', 'g2', 'g3']
    assert rows[1]['prediction'] is None
    assert rows[2]['truth'] != '["2"]', 'False class calls must remain available'
    assert [row['entity'] for row in snapshot(view)['details']['failures']] == ['g2', 'g3']
    assert snapshot(view)['card'] == benchmark.profile_class_cards[0]


def test_profile_custom_precision_recall_f1_prevalence_keep_native_values():
    benchmark = fixture()
    native = benchmark.profile_class_cards[0]
    view = C.build(benchmark, 'profile:["2"]')
    values = metrics(view)
    for shown, original in [('class_precision', 'precision'), ('class_recall', 'recall'),
                            ('class_f1', 'f1'), ('reference_prevalence', 'evaluation_prevalence')]:
        assert values[shown].value == native['class_metrics'][original]
    html = render_scorecard_html(view)
    assert 'Class precision, whole cohort' in html
    assert 'Reference-member accuracy among calls' in html
    assert 'Reference prevalence, whole cohort' in html
    assert 'biological precision' in html and 'Calibrated confidence unavailable' in html
    assert 'false positives' in render_scorecard_detail(view, 'class_precision')
    assert 'abstentions' in render_scorecard_detail(view, 'class_recall')
    assert 'whole frozen test cohort' in render_scorecard_detail(view, 'class_f1')


def test_major_card_uses_whole_cohort_counts_and_overall_profile_call_coverage():
    benchmark = fixture()
    native = benchmark.major_class_cards[0]
    view = C.build(benchmark, 'major:2')
    assert dict(view.counts) == {key: native[key] for key in C._COUNT_KEYS}
    values = metrics(view)
    assert set(values) == {'coverage', 'class_precision', 'class_recall', 'class_f1', 'reference_prevalence'}
    assert values['coverage'].value == .75
    assert values['coverage'].denominator == 'All frozen test genes: 4'
    assert 'not the share called this major class' in values['coverage'].definition
    assert values['class_precision'].value == native['precision']
    assert snapshot(view)['card']['native_class_card'] == native
    assert snapshot(view)['details']['native_class_card'] == native
    assert len(snapshot(view)['details']['rows']) == 3
    assert 'retained-row count is not the accuracy' in render_scorecard_detail(view, 'gaps')
    assert 'Accuracy, all eligible' not in render_scorecard_html(view)


def test_missing_major_class_truth_and_calls_keep_null_precision_recall_f1():
    view = C.build(fixture(), 'major:3')
    values = metrics(view)
    assert all(values[key].value is None for key in ('class_precision', 'class_recall', 'class_f1'))
    assert values['reference_prevalence'].value == 0.
    assert snapshot(view)['details']['rows'] == []
    assert snapshot(view)['details']['failures'] == []
    assert 'Unavailable' in render_scorecard_html(view)
    assert 'reports null' in render_scorecard_detail(view, 'class_precision')


def test_profile_zero_no_call_convention_is_preserved_without_guessing_null():
    benchmark = fixture()
    native = benchmark.profile_class_cards[0]
    native['class_metrics'].update(precision=0., recall=0., f1=0., true_positive=0,
                                  false_positive=0, false_negative=2)
    for row in benchmark.rows:
        if row['prediction'] == '["2"]':
            row.update(prediction=None, abstained=True)
    native['counts'].update(answered=0, abstained=2, correct=0)
    native['metrics'].update(accuracy=0., precision_of_calls=None, coverage=0.)
    view = C.build(benchmark, 'profile:["2"]')
    assert metrics(view)['class_precision'].value == 0.
    assert metrics(view)['precision_of_calls'].value is None
    assert 'reports zero' in render_scorecard_detail(view, 'class_precision')


def test_no_matched_class_baseline_or_biological_support_is_invented():
    view = C.build(fixture(), 'major:2')
    payload = snapshot(view)
    assert payload['details']['baseline']['status'] == 'unavailable'
    assert 'Whole-cohort controls' in payload['details']['controls']['scope']
    assert payload['details']['native_class_card']['biological_precision'] is None
    assert payload['details']['native_class_card']['calibrated_confidence'] is None
    assert 'biological precision and biological recall unavailable' in view.status
    assert view.source.grade == 'synthetic_control'


def test_detached_export_and_escaped_html_preserve_source_context_and_input():
    benchmark = fixture()
    benchmark.metadata['source_context'] = '<script>fixture()</script>'
    original = deepcopy(benchmark.payloads)
    view = C.build(benchmark, 'profile:["2"]')
    assert benchmark.payloads == original
    before = export_scorecard(view)
    benchmark.rows[0]['prediction'] = None
    benchmark.profile_class_cards[0]['class_metrics']['precision'] = .25
    assert export_scorecard(view) == before
    html = render_scorecard_html(view, expanded=True)
    assert '<script>' not in html and '&lt;script&gt;' in html
    assert 'scorecard:rows/selected' in html
    assert 'fixture identity' in render_scorecard_detail(view, 'freshness')
    assert view.source.sha256 == 'a' * 64
    assert view.source.sha256 != benchmark.metadata['source_artifact_identity']


def test_source_table_hash_and_artifact_identity_have_separate_provenance_roles():
    benchmark = fixture()
    benchmark.metadata['source_manifest']['spec']['dependencies'] = []
    view = C.build(benchmark, 'major:2')
    assert view.source.sha256 == 'Unavailable / not supplied'
    assert snapshot(view)['details']['freshness']['artifact_identity'] == 'fixture identity'
    assert snapshot(view)['details']['freshness']['source_table_sha256'] is None
    benchmark = fixture()
    dependency = benchmark.metadata['source_manifest']['spec']['dependencies'][0]
    benchmark.metadata['source_manifest']['spec']['dependencies'].append(deepcopy(dependency))
    with pytest.raises(ValueError, match='identity must be unique'):
        C.build(benchmark, 'major:2')


@pytest.mark.parametrize('address', ['label:', 'major:', 'profile:', 'major:8', 'profile:["7"]', '', None])
def test_missing_or_unsupported_class_address_is_refused(address):
    with pytest.raises(ValueError):
        C.build(fixture(), address)


def test_wrong_input_duplicate_class_and_unsupported_admission_are_refused():
    with pytest.raises(TypeError):
        C.build({}, 'major:2')
    benchmark = fixture()
    benchmark.major_class_cards.append(deepcopy(benchmark.major_class_cards[0]))
    with pytest.raises(ValueError, match='exactly one'):
        C.build(benchmark, 'major:2')
    for field, value in [('biological_admission', True), ('calibrated_confidence', .9)]:
        benchmark = fixture()
        benchmark.metadata[field] = value
        with pytest.raises(ValueError, match='admission or calibrated confidence'):
            C.build(benchmark, 'major:2')


def test_verified_existing_tg_cards_retain_native_metrics_counts_and_rows():
    benchmarks, reason = F.shipped(O.TOXOPLASMA)
    assert benchmarks and not reason
    benchmark = benchmarks[0]
    source_sha256 = next(entry['sha256'] for entry in benchmark.metadata['source_manifest']['spec']['dependencies']
                        if entry['kind'] == 'table' and entry['name'] == 'installed_nodes')
    for native in benchmark.profile_class_cards:
        value = native['class']
        view = C.build(benchmark, 'profile:' + value)
        assert view.source.sha256 == source_sha256
        shown = metrics(view)
        assert dict(view.counts) == native['counts']
        assert snapshot(view)['card'] == native
        for key, number in native['metrics'].items():
            assert shown[key].value == number
        for key in ('precision', 'recall', 'f1'):
            assert shown['class_' + key].value == native['class_metrics'][key]
        expected = [row for row in benchmark.rows if row['truth'] == value or row['prediction'] == value]
        assert snapshot(view)['details']['rows'] == expected
    for native in benchmark.major_class_cards:
        value = native['class']
        view = C.build(benchmark, 'major:' + value)
        assert view.source.sha256 == source_sha256
        shown = metrics(view)
        assert snapshot(view)['details']['native_class_card'] == native
        assert shown['coverage'].value == native['coverage']
        for key in ('precision', 'recall', 'f1'):
            assert shown['class_' + key].value == native[key]
        assert shown['reference_prevalence'].value == native['reference_prevalence']
        expected = [row for row in benchmark.rows if value in F.profile_classes(row['truth'])
            or value in (F.profile_classes(row['prediction']) or ())]
        assert snapshot(view)['details']['rows'] == expected
