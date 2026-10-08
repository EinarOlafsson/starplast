"""Synthetic presentation contracts: denominators, nulls, provenance and navigation."""
from dataclasses import FrozenInstanceError
import json

import pytest

from starplast import scorecard as SC
from starplast.scorecard_view import (
    ScorecardLink, build_scorecard_view, export_scorecard, render_scorecard_detail,
    render_scorecard_html, validate_scorecard_link,
)


def fixture(task=SC.T_LABEL):
    return {
        'scope': {'task': task, 'organism': 'fixture organism', 'strategy': 'fixture strategy',
                  'target': 'fixture target', 'settings': 'seed=3', 'protocol': 'group holdout',
                  'partition': 'test', 'truth_grade': 'synthetic_control', 'unit': 'gene',
                  'benchmark_id': 'fixture', 'negative_semantics': 'unknown is not negative',
                  'gaps': ['Fixture only; no biological validation']},
        'counts': {'eligible': 4, 'answered': 3, 'abstained': 1, 'correct': 2, 'wrong': 1},
        'metrics': {'accuracy': .5, 'precision_of_calls': 2 / 3, 'coverage': .75}
            if task == SC.T_LABEL else {SC.TASKS[task].metrics[0]: .125},
        'extra': {'uncertainty': {'status': 'unavailable', 'reason': 'Groups not supplied'}},
        'scope_identity': 'fixture scope', 'records_identity': 'fixture records',
    }


def by_key(view):
    return {metric.key: metric for metric in view.metrics}


def test_all_eligible_and_called_denominators_survive_one_model_in_both_views():
    card = fixture()
    first = build_scorecard_view(card, title='Target')
    second = build_scorecard_view(card, title='Strategy')
    assert first.metrics == second.metrics
    values = by_key(first)
    assert values['accuracy'].value == .5
    assert values['precision_of_calls'].value == 2 / 3
    assert values['accuracy'].denominator == 'All eligible: 4'
    assert values['precision_of_calls'].denominator == 'Calls made: 3'
    html = render_scorecard_html(first)
    assert 'Accuracy, all eligible' in html and 'Accuracy among calls' in html
    assert 'All eligible: 4' in html and 'Calls made: 3' in html
    assert values['accuracy'].definition == SC.METRICS['accuracy'].definition
    assert 'no biological validation' in render_scorecard_detail(first, 'gaps')


@pytest.mark.parametrize('task', list(SC.TASKS))
def test_all_six_task_structures_retain_authority_order_values_and_nulls(task):
    card = fixture(task)
    view = build_scorecard_view(card)
    assert tuple(metric.key for metric in view.metrics) == SC.TASKS[task].metrics
    for metric in view.metrics:
        assert metric.value == card['metrics'].get(metric.key)
        assert metric.definition == SC.METRICS[metric.key].definition
    assert any(metric.headline for metric in view.metrics)
    assert 'Unavailable' in render_scorecard_html(view, expanded=True)
    if task != SC.T_LABEL:
        assert 'Accuracy among calls' not in render_scorecard_html(view)


def test_no_calls_is_undefined_called_accuracy_and_zero_all_eligible_accuracy():
    card = fixture()
    card['counts'].update(answered=0, abstained=4, correct=0, wrong=0)
    card['metrics'].update(accuracy=0., precision_of_calls=None, coverage=0.)
    view = build_scorecard_view(card)
    assert by_key(view)['accuracy'].value == 0
    assert by_key(view)['precision_of_calls'].value is None
    payload = json.loads(export_scorecard(view))
    called = next(metric for metric in payload['metrics'] if metric['key'] == 'precision_of_calls')
    assert called['value'] is None
    assert payload['snapshot']['card']['metrics']['precision_of_calls'] is None
    assert 'Calls made: 0' in render_scorecard_html(view)


def test_flat_metric_contract_requires_task_and_normalizes_nan_to_null():
    with pytest.raises(ValueError, match='explicit task'):
        build_scorecard_view({'mae': 2.})
    view = build_scorecard_view({'mae': 2., 'pearson': float('nan')}, task=SC.T_VALUES)
    assert by_key(view)['mae'].value == 2.
    assert by_key(view)['pearson'].value is None
    assert 'NaN' not in export_scorecard(view)
    wrapped = build_scorecard_view({'task': SC.T_VALUES, 'mae': 2., 'status': 'Candidate'})
    assert by_key(wrapped)['mae'].value == 2.
    assert wrapped.status == 'Candidate'


def test_set_coverage_and_size_are_present_together_with_native_definitions():
    card = fixture()
    card['metrics'].update(set_coverage=.9, mean_set_size=8., singleton_share=0.)
    view = build_scorecard_view(card)
    values = by_key(view)
    assert values['set_coverage'].value == .9
    assert values['mean_set_size'].value == 8.
    html = render_scorecard_html(view)
    assert 'Prediction-set coverage' in html and 'Mean prediction-set size' in html
    assert 'An empty set is uncovered' in render_scorecard_detail(view, 'set_coverage')
    assert values['set_coverage'].definition == SC.metric_definition(SC.T_LABEL, 'set_coverage').definition


def test_unknown_metric_retains_value_with_explicit_missing_definition():
    card = fixture()
    card['metrics']['fixture_measure'] = .9
    view = build_scorecard_view(card)
    assert by_key(view)['fixture_measure'].value == .9
    assert 'Definition unavailable' in render_scorecard_detail(view, 'fixture_measure')


def test_task_aware_nominal_set_and_interval_coverage_do_not_change_native_schema():
    label = fixture()
    label['metrics'].update(promised_coverage=.9, set_efficiency=.8, empty_set_share=.25)
    value = fixture(SC.T_VALUES)
    value['metrics']['promised_coverage'] = .9
    label_view = build_scorecard_view(label)
    value_view = build_scorecard_view(value)
    assert by_key(label_view)['promised_coverage'].label == 'Nominal prediction-set coverage'
    assert by_key(value_view)['promised_coverage'].label == 'Nominal interval coverage'
    assert 'clip(set size, 1, C)' in render_scorecard_detail(label_view, 'set_efficiency')
    assert 'coverage misses' in render_scorecard_detail(label_view, 'set_efficiency')
    assert SC.METRICS['promised_coverage'].task == SC.T_VALUES
    supplemental_keys = {key for task, key in SC.SUPPLEMENTAL_METRICS if task == SC.T_LABEL}
    assert supplemental_keys.isdisjoint(SC.TASKS[SC.T_LABEL].metrics)
    assert tuple(SC._ordered(SC.T_LABEL, label['metrics'])) == SC.TASKS[SC.T_LABEL].metrics
    assert SC.metric_definition(SC.T_VALUES, 'set_coverage') is None
    assert SC.metric_definition(SC.T_LABEL, 'fixture_measure') is None
    with pytest.raises(ValueError, match='different task'):
        value['metrics']['set_coverage'] = .9
        build_scorecard_view(value)


def test_numeric_interval_coverage_never_hides_width_or_baseline_matching():
    card = fixture(SC.T_VALUES)
    card['metrics'].update(interval_coverage=.9, mean_interval_width=None,
                           unbounded_interval_share=.5, value_coverage=.75)
    baseline = {'name': 'train mean', 'eligible_rows': 4, 'matched_answered_rows': 3,
                'eligible_metrics': {'mae': 30.}, 'matched_metrics': {'mae': 2.},
                'mae_skill': .5, 'mse_skill': None}
    card['extra']['baseline_comparison'] = baseline
    view = build_scorecard_view(card)
    html = render_scorecard_html(view)
    assert 'Interval coverage, all eligible' in html
    assert 'Mean interval width' in html and 'Unavailable' in html
    assert 'matched_answered_rows' in render_scorecard_detail(view, 'baseline')
    assert json.loads(export_scorecard(view))['snapshot']['card']['extra']['baseline_comparison'] == baseline


@pytest.mark.parametrize('location', ['details', 'card', 'extra'])
def test_calibration_record_has_dedicated_escaped_route_and_exact_snapshot(location):
    card = fixture()
    calibration = {'status': 'recorded fixture only', 'population': ['c1', 'c2'],
        'source': '<script>fixture</script>', 'nominal_coverage': .9,
        'confidence_kind': 'prediction_set', 'biological_admission': False}
    details = {'calibration': calibration} if location == 'details' else {}
    if location == 'card':
        card['calibration'] = calibration
    elif location == 'extra':
        card['extra']['calibration'] = calibration
    view = build_scorecard_view(card, details=details)
    record = next(detail for detail in view.details if detail.key == 'calibration')
    assert json.loads(record.text) == calibration
    assert 'scorecard:detail/calibration' in render_scorecard_html(view)
    assert validate_scorecard_link('scorecard:detail/calibration') == 'scorecard:detail/calibration'
    popup = render_scorecard_detail(view, 'calibration')
    assert 'c1' in popup and '&lt;script&gt;' in popup and '<script>' not in popup
    exported = json.loads(export_scorecard(view))
    assert exported['snapshot'] == {'card': card, 'details': details}
    assert by_key(view)['accuracy'].value == card['metrics']['accuracy']


def test_missing_calibration_is_explicit_and_does_not_add_snapshot_fields():
    card = fixture()
    view = build_scorecard_view(card)
    popup = render_scorecard_detail(view, 'calibration')
    assert 'unavailable' in popup and 'Calibration metadata not supplied' in popup
    assert 'Calibration metadata not supplied' in render_scorecard_html(view)
    assert json.loads(export_scorecard(view))['snapshot'] == {'card': card, 'details': {}}


def test_identity_only_freshness_exposes_unknown_date_and_version_without_inference():
    card = fixture()
    freshness = {'artifact_identity': 'fixture immutable identity', 'source_table_sha256': 'a' * 64}
    details = {'freshness': freshness}
    view = build_scorecard_view(card, details=details)
    record = json.loads(next(detail.text for detail in view.details if detail.key == 'freshness'))
    assert record['artifact_identity'] == freshness['artifact_identity']
    assert record['date_availability'] == record['version_availability'] == 'Unavailable / not supplied'
    assert not any(key in record for key in ('date', 'created_at', 'source_version'))
    popup = render_scorecard_detail(view, 'freshness')
    assert 'date_availability' in popup and 'version_availability' in popup
    assert json.loads(export_scorecard(view))['snapshot'] == {'card': card, 'details': details}


def test_supplied_freshness_values_and_explicit_nulls_are_preserved():
    card = fixture()
    known = {'retrieved_at': 'fixture date', 'source_version': 'fixture version', 'artifact_identity': 'identity'}
    view = build_scorecard_view(card, details={'freshness': known})
    assert json.loads(next(detail.text for detail in view.details if detail.key == 'freshness')) == known
    unknown = {'date': None, 'source_version': 'unresolved', 'artifact_identity': 'identity'}
    view = build_scorecard_view(card, details={'freshness': unknown})
    shown = json.loads(next(detail.text for detail in view.details if detail.key == 'freshness'))
    assert shown['date'] is None and shown['source_version'] == 'unresolved'
    assert shown['date_availability'] == shown['version_availability'] == 'Unavailable / not supplied'
    assert json.loads(export_scorecard(view))['snapshot']['details']['freshness'] == unknown


def test_unknown_freshness_markers_ignore_case_and_space_preserving_raw_metadata():
    freshness = {'retrieved_at': ' Unknown ', 'source_version': ' UNRESOLVED ',
                 'artifact_identity': 'fixture identity'}
    view = build_scorecard_view(fixture(), details={'freshness': freshness})
    shown = json.loads(next(detail.text for detail in view.details if detail.key == 'freshness'))
    assert shown['retrieved_at'] == freshness['retrieved_at']
    assert shown['source_version'] == freshness['source_version']
    assert shown['date_availability'] == shown['version_availability'] == 'Unavailable / not supplied'
    assert json.loads(export_scorecard(view))['snapshot']['details']['freshness'] == freshness


def test_source_grade_context_lineage_and_unknown_negatives_remain_distinct():
    card = fixture()
    card['scope']['truth_grade'] = 'prediction'
    metadata = {'source': {'name': 'Fixture classifier output', 'grade': 'prediction',
        'version': 'v1', 'sha256': 'fixture hash', 'context': 'stored annotation only',
        'lineage': 'prediction -> transfer', 'negative_semantics': 'unmeasured is unknown'},
        'freshness': {'retrieved': 'fixture date', 'source_release': 'fixture release'}}
    view = build_scorecard_view(card, status='Candidate recovery; biological admission unresolved', details=metadata)
    assert view.source.grade == 'prediction'
    assert view.source.lineage == 'prediction -> transfer'
    assert view.source.negative_semantics == 'unmeasured is unknown'
    html = render_scorecard_html(view, expanded=True)
    assert 'biological admission unresolved' in html
    assert 'stored annotation only' in html and 'fixture release' in html
    assert 'calibrated confidence' not in html


def test_evidence_and_individual_outcome_cards_show_no_per_gene_accuracy():
    evidence = build_scorecard_view({'evidence': {'mapped': 2, 'unmeasured': None},
        'source': {'grade': 'curation'}}, kind='evidence')
    outcome = build_scorecard_view({'entity': 'g1', 'outcome': {'truth': 'A',
        'prediction': None, 'abstained': True}, 'per_gene_accuracy_probability': None}, kind='outcome')
    for view in (evidence, outcome):
        assert view.metrics == ()
        assert 'Result quality / coverage' not in render_scorecard_html(view)
    assert 'Evidence quality' in render_scorecard_html(evidence)
    assert 'Individual outcome' in render_scorecard_html(outcome)
    assert 'no per-gene accuracy probability' in render_scorecard_html(outcome)
    assert 'null' in render_scorecard_detail(outcome, 'outcome')


def test_unavailable_card_and_single_outcome_caveat_stay_visible():
    view = build_scorecard_view({'task': SC.T_LABEL, 'metrics': {}})
    assert 'Untested / unavailable' in render_scorecard_html(view)
    assert 'Unavailable / not supplied' in render_scorecard_detail(view, 'freshness')
    card = fixture()
    card['counts'] = {'eligible': 1, 'answered': 1, 'correct': 1, 'abstained': 0, 'wrong': 0}
    assert 'Single recorded outcome' in render_scorecard_html(build_scorecard_view(card))


def test_immutable_model_is_detached_and_exports_original_exact_values():
    card = fixture()
    details = {'rows': [{'entity': 'g1', 'prediction': None}], 'controls': {'fixture': True}}
    view = build_scorecard_view(card, details=details)
    before = export_scorecard(view)
    card['metrics']['accuracy'] = .1
    details['rows'][0]['prediction'] = 'changed'
    assert export_scorecard(view) == before
    with pytest.raises(FrozenInstanceError):
        view.title = 'changed'
    exported = json.loads(before)
    assert exported['snapshot']['card']['metrics']['precision_of_calls'] == 2 / 3
    assert exported['snapshot']['details']['rows'][0]['prediction'] is None
    assert exported['source']['grade'] == 'synthetic_control'


def test_rendering_and_popup_escape_all_user_metadata_and_link_labels():
    card = fixture()
    card['scope']['target'] = '<img src=x onerror=bad()>'
    view = build_scorecard_view(card, title='<script>bad()</script>',
        details={'failures': ['<b>bad</b>'], 'rows': '<script>rows()</script>'},
        links=(ScorecardLink('<i>source</i>', 'https://example.org/paper?a=1&b=2'),))
    html = render_scorecard_html(view, expanded=True)
    assert '<script>' not in html and '<img ' not in html and '<i>' not in html
    assert '&lt;script&gt;' in html and '&amp;b=2' in html
    assert '&lt;b&gt;bad&lt;/b&gt;' in render_scorecard_detail(view, 'failures')
    assert 'scorecard:detail/source' in html and 'scorecard:detail/split' in html
    assert 'scorecard:detail/controls' in html and 'scorecard:detail/rows' in html
    assert 'ALL hidden genes' in render_scorecard_detail(view, 'accuracy')
    with pytest.raises(KeyError):
        render_scorecard_detail(view, 'invented')


@pytest.mark.parametrize('url', ['javascript:alert(1)', 'file:///tmp/test', 'http://example.org',
    'https://user:secret@example.org', 'https://example.org\n/evil', 'https://example.org/%0aevil',
    'https://example.org\\evil', 'scorecard:metric/a/b', 'scorecard:metric/%2f..',
    'scorecard:delete/all', 'scorecard:metric/<script>', 'https://example.org:wrong'])
def test_invalid_navigation_is_refused_before_render_and_dispatch(url):
    with pytest.raises(ValueError):
        validate_scorecard_link(url)
    with pytest.raises(ValueError):
        ScorecardLink('destination', url)


@pytest.mark.parametrize('url', ['scorecard:expand', 'scorecard:metric/accuracy',
    'scorecard:detail/split', 'scorecard:rows/fixture', 'scorecard:outcome/gene.1',
    'https://example.org/paper#figure'])
def test_approved_routes_are_stable(url):
    assert validate_scorecard_link(url) == url


def test_contradictory_task_or_counts_and_non_numeric_rates_are_refused():
    with pytest.raises(ValueError, match='recorded task'):
        build_scorecard_view(fixture(), task=SC.T_VALUES)
    broken = fixture()
    broken['counts']['answered'] = 8
    with pytest.raises(ValueError, match='eligible population'):
        build_scorecard_view(broken)
    broken = fixture()
    broken['metrics']['accuracy'] = True
    with pytest.raises(ValueError, match='numeric'):
        build_scorecard_view(broken)
    broken = fixture()
    broken['metrics']['mae'] = 3.
    with pytest.raises(ValueError, match='different task'):
        build_scorecard_view(broken)
