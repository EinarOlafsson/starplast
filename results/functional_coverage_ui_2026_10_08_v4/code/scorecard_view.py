"""Immutable scorecard presentation shared by HTML and desktop hosts.

This module presents one supplied evaluation; it never aggregates cohorts,
admits biological truth, estimates uncertainty or converts method support into
confidence. Metric definitions come from :mod:`starplast.scorecard`.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass
from html import escape
import json
import math
from numbers import Real
import re
from urllib.parse import urlsplit

from . import scorecard as SC


@dataclass(frozen=True)
class ScorecardLink:
    """A named navigation destination; hosts validate again before dispatch."""

    label: str
    url: str

    def __post_init__(self):
        if not isinstance(self.label, str) or not self.label:
            raise ValueError('Scorecard links require a label')
        validate_scorecard_link(self.url)


@dataclass(frozen=True)
class SourceView:
    """Explicit source interpretation, independent of performance measurements."""

    grade: str
    name: str
    version: str
    sha256: str
    context: str
    lineage: str
    negative_semantics: str


@dataclass(frozen=True)
class MetricView:
    """One supplied metric with its task definition and explicit denominator."""

    key: str
    label: str
    value: float | int | None
    status: str
    definition: str
    reading: str
    range: str
    chance: str
    denominator: str
    headline: bool = False


@dataclass(frozen=True)
class DetailView:
    """An immutable full-record section reached from the compact card."""

    key: str
    label: str
    text: str


@dataclass(frozen=True)
class ScorecardView:
    """One immutable public interface for a scorecard, popup and JSON export.

    ``kind`` is performance, evidence or outcome. All mutable input is copied
    into scalar/tuple values and a canonical JSON snapshot. No numeric value is
    recomputed. ``status`` describes the supplied evaluation, not truth admission.
    """

    kind: str
    title: str
    status: str
    task: str | None
    scope: tuple[tuple[str, str], ...]
    counts: tuple[tuple[str, int | None], ...]
    source: SourceView
    metrics: tuple[MetricView, ...]
    details: tuple[DetailView, ...]
    links: tuple[ScorecardLink, ...]
    snapshot_json: str


_KINDS = {'performance', 'evidence', 'outcome'}
_HEADLINES = {
    SC.T_LABEL: ('accuracy', 'precision_of_calls', 'coverage', 'set_coverage', 'mean_set_size'),
    SC.T_RANK: ('auprc', 'auprc_lift', 'recall_at_10pct'),
    SC.T_SET: ('precision', 'recall', 'returned'),
    SC.T_CLUSTER: ('weighted_f1_clusters', 'ari', 'noise_share'),
    SC.T_VALUES: ('mae', 'rmse', 'value_coverage', 'interval_coverage', 'mean_interval_width'),
    SC.T_REPL: ('replication_rate', 'replication_lift', 'findings'),
}
_SCOPE_KEYS = ('organism', 'target', 'strategy', 'settings', 'seed', 'protocol',
               'partition', 'benchmark_id', 'unit')
_DETAIL_LABELS = {
    'source': 'Ground-truth source and interpretation', 'split': 'Split and fitted populations',
    'sizes': 'Sample sizes', 'baseline': 'Matched baseline', 'uncertainty': 'Uncertainty',
    'calibration': 'Calibration record',
    'controls': 'Controls', 'failures': 'Failure examples', 'rows': 'Row-level outcomes',
    'freshness': 'Freshness', 'gaps': 'Limitations and gaps', 'identity': 'Evaluation identity',
    'outcome': 'Individual outcome', 'evidence': 'Evidence quality',
    'parameters': 'Evaluation parameters', 'task_details': 'Task-specific test record',
}


def validate_scorecard_link(url: str) -> str:
    """Return an allowed destination or reject it before rendering/navigation.

    Internal routes are ``scorecard:expand`` and ``scorecard:<metric|detail|rows|
    outcome>/<identifier>``. HTTPS sources cannot contain credentials, whitespace,
    control characters, backslashes or encoded control characters. Arbitrary file,
    script and application schemes are refused.
    """
    if (not isinstance(url, str) or not url or re.search(r'[\s\\\x00-\x1f\x7f]', url)
            or re.search(r'%(?:0[0-9a-f]|1[0-9a-f]|7f)', url, re.I)):
        raise ValueError('Invalid scorecard destination')
    if url == 'scorecard:expand':
        return url
    if re.fullmatch(r'scorecard:(?:metric|detail|rows|outcome)/[A-Za-z0-9_.:-]+', url):
        return url
    parts = urlsplit(url)
    if (parts.scheme == 'https' and parts.hostname and parts.username is None
            and parts.password is None):
        try:
            parts.port
        except ValueError as exc:
            raise ValueError('Invalid scorecard destination port') from exc
        return url
    raise ValueError('Unsupported scorecard destination')


def _plain(value):
    """Copy JSON-like presentation metadata; preserve null and missing rates."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, Real):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) for key in value):
            raise TypeError('Scorecard metadata keys must be strings')
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(item) for item in value]
    raise TypeError('Scorecard metadata requires JSON scalars, lists and mappings')


def _text(value):
    if value is None or value == '':
        return 'Unavailable / not supplied'
    if isinstance(value, str):
        return value
    return json.dumps(value, sort_keys=True, ensure_ascii=False, indent=2, allow_nan=False)


def _number(value):
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError('Metrics require finite numeric values or null')
    return value if math.isfinite(value) else None


def _denominator(key, counts):
    if key in {'accuracy', 'coverage', 'value_coverage', 'interval_coverage',
               'interval_availability', 'set_coverage', 'singleton_share', 'empty_set_share'}:
        return 'All eligible: ' + _text(counts.get('eligible'))
    if key == 'precision_of_calls':
        return 'Calls made: ' + _text(counts.get('answered'))
    return ''


def _freshness_known(value):
    if isinstance(value, bool):
        return False
    if isinstance(value, (int, float)):
        return True
    return isinstance(value, str) and value.strip().casefold() not in {
        '', 'unknown', 'unresolved', 'unavailable', 'unavailable / not supplied', 'not supplied', 'n/a'}


def _freshness_detail(value):
    """Keep supplied metadata and mark missing date/version information explicitly."""
    metadata = dict(value) if isinstance(value, dict) else {'supplied_metadata': value}
    for status, keys in (
        ('date_availability', ('date', 'created_at', 'evaluated_at', 'generated_at',
                               'retrieved', 'retrieved_at', 'retrieval_date', 'source_release_date')),
        ('version_availability', ('version', 'source_version', 'source_release')),
    ):
        known = any(_freshness_known(metadata.get(key)) for key in keys)
        if not known:
            metadata[status] = 'Unavailable / not supplied'
    return metadata


def build_scorecard_view(card: Mapping, *, task: str | None = None,
                        kind: str = 'performance', title: str | None = None,
                        status: str | None = None, details: Mapping | None = None,
                        links: tuple[ScorecardLink, ...] = ()) -> ScorecardView:
    """Adapt a record card or a standard metric dictionary without changing values.

    ``details`` supplies source, split, sizes, baseline, uncertainty, calibration, controls,
    failures, rows, freshness and gaps when the card lacks them. ``source`` is a
    mapping of name/version/sha256/context/lineage/negative_semantics/grade. A
    source grade is provenance, never an automatic biological-admission verdict.
    Outcome cards accept the individual ``outcome`` dictionary from record views;
    they show that outcome and no per-entity performance rate.
    """
    if not isinstance(card, Mapping) or kind not in _KINDS:
        raise ValueError('One scorecard mapping and explicit card kind required')
    supplied = _plain(card)
    supplemental = _plain(details or {})
    scope = supplied.get('scope') or {}
    if not isinstance(scope, dict):
        raise ValueError('Scope must be a mapping')
    declared = scope.get('task', supplied.get('task'))
    if task is not None and declared is not None and task != declared:
        raise ValueError('Requested task differs from the recorded task')
    task = declared if task is None else task
    if task is not None and task not in SC.TASKS:
        raise ValueError('Unknown scorecard task')
    if kind == 'performance' and task is None:
        raise ValueError('Performance cards require an explicit task')
    counts = supplied.get('counts') or {}
    if not isinstance(counts, dict):
        raise ValueError('Counts must be a mapping')
    for value in counts.values():
        if value is not None and (type(value) is not int or value < 0):
            raise ValueError('Sample counts require nonnegative integers or null')
    eligible, answered, correct = (counts.get(key) for key in ('eligible', 'answered', 'correct'))
    if eligible is not None and answered is not None:
        if answered > eligible or (counts.get('abstained') is not None
                and counts['abstained'] != eligible - answered):
            raise ValueError('Answered and abstained counts contradict the eligible population')
    if answered is not None and correct is not None:
        if correct > answered or (counts.get('wrong') is not None and counts['wrong'] != answered - correct):
            raise ValueError('Correct and wrong counts contradict the called population')
    metrics = supplied.get('metrics', {key: value for key, value in supplied.items()
        if key not in {'task', 'status', 'title', 'source', 'freshness', 'uncertainty',
                       'baseline', 'controls', 'failure_examples', 'rows'}}
        if task and 'scope' not in supplied and 'counts' not in supplied else {})
    if not isinstance(metrics, dict):
        raise ValueError('Metrics must be a mapping')
    extra = supplied.get('extra') or {}
    if not isinstance(extra, dict):
        raise ValueError('Extra scorecard metadata must be a mapping')
    source = supplemental.get('source', supplied.get('source')) or {}
    if not isinstance(source, dict):
        raise ValueError('Source interpretation must be a mapping')
    source_view = SourceView(
        _text(source.get('grade', scope.get('truth_grade'))), _text(source.get('name')),
        _text(source.get('version')), _text(source.get('sha256')), _text(source.get('context')),
        _text(source.get('lineage')), _text(source.get('negative_semantics', scope.get('negative_semantics'))))
    metric_views = []
    if kind == 'performance':
        keys = list(SC.TASKS[task].metrics)
        keys.extend(key for key in metrics if key not in keys)
        for key in keys:
            if not re.fullmatch(r'[A-Za-z0-9_.:-]+', key):
                raise ValueError('Metric keys require a safe stable identifier')
            authority = SC.metric_definition(task, key)
            if authority is None and (key in SC.METRICS or any(
                    defined_key == key for _, defined_key in SC.SUPPLEMENTAL_METRICS)):
                raise ValueError('Metric belongs to a different task: ' + key)
            value = _number(metrics.get(key))
            label = {'accuracy': 'Accuracy, all eligible',
                     'precision_of_calls': 'Accuracy among calls'}.get(key,
                        authority.label if authority else key.replace('_', ' ').capitalize())
            metric_views.append(MetricView(key, label, value,
                'available' if value is not None else supplied.get('metric_status', {}).get(key, 'unavailable'),
                authority.definition if authority else 'Definition unavailable in the shared metric glossary.',
                authority.reading if authority else 'Inspect the task-specific test record and evaluation parameters.',
                authority.range if authority else 'Unspecified', authority.chance if authority else 'Unspecified',
                _denominator(key, counts), key in _HEADLINES[task]))
    available = any(metric.value is not None for metric in metric_views)
    state = status or supplied.get('status') or ('Recorded evaluation' if available else 'Untested / unavailable')
    if not isinstance(state, str):
        raise ValueError('Scorecard status must be explicit text')
    if kind == 'outcome' and status is None and 'status' not in supplied:
        state = 'Individual outcome; no per-gene accuracy probability'
    content = {
        'source': asdict(source_view),
        'split': supplemental.get('split', {key: scope.get(key) for key in ('protocol', 'partition')}),
        'sizes': supplemental.get('sizes', counts or None),
        'baseline': supplemental.get('baseline', extra.get('baseline_comparison', supplied.get('baseline'))),
        'uncertainty': supplemental.get('uncertainty', extra.get('uncertainty', supplied.get('uncertainty'))),
        'calibration': supplemental.get('calibration', supplied.get('calibration', extra.get('calibration'))),
        'controls': supplemental.get('controls', supplied.get('controls')),
        'failures': supplemental.get('failures', supplied.get('failure_examples')),
        'rows': supplemental.get('rows', supplied.get('rows')),
        'freshness': _freshness_detail(supplemental.get('freshness', supplied.get('freshness'))),
        'gaps': supplemental.get('gaps', scope.get('gaps')),
        'identity': {key: supplied.get(key) for key in ('scope_identity', 'cohort_identity',
                    'records_identity', 'parameter_identity')},
        'parameters': supplied.get('evaluation_parameters'), 'task_details': extra or None,
    }
    if content['calibration'] is None:
        content['calibration'] = {'status': 'unavailable', 'reason': 'Calibration metadata not supplied'}
    if kind == 'outcome':
        content['outcome'] = supplied.get('outcome', supplied)
    if kind == 'evidence':
        content['evidence'] = supplied.get('evidence', supplied)
    detail_views = tuple(DetailView(key, _DETAIL_LABELS[key], _text(value)) for key, value in content.items())
    link_views = tuple(links)
    if any(not isinstance(link, ScorecardLink) for link in link_views):
        raise TypeError('Links require typed ScorecardLink instances')
    snapshot = json.dumps({'card': supplied, 'details': supplemental}, sort_keys=True,
                          ensure_ascii=False, separators=(',', ':'), allow_nan=False)
    selected_title = title or {'performance': 'Performance scorecard',
        'evidence': 'Evidence quality', 'outcome': 'Individual outcome'}[kind]
    if not isinstance(selected_title, str):
        raise ValueError('Scorecard title must be text')
    return ScorecardView(kind, selected_title, state, task,
        tuple((key, _text(scope.get(key))) for key in _SCOPE_KEYS), tuple(counts.items()),
        source_view, tuple(metric_views), detail_views, link_views, snapshot)


def _value(value):
    return 'Unavailable' if value is None else format(value, '.6g')


def _anchor(label, url):
    return '<a href="' + escape(validate_scorecard_link(url), quote=True) + '">' + escape(label) + '</a>'


def _compact_detail(detail):
    """Keep the supplied full record behind its link; summarize supplied fields."""
    try:
        value = json.loads(detail.text)
    except (ValueError, TypeError):
        return detail.text
    if not isinstance(value, dict):
        return detail.text
    if detail.key == 'sizes':
        keys = ('eligible', 'answered', 'abstained', 'correct', 'wrong')
    elif detail.key in {'uncertainty', 'calibration'}:
        keys = ('status', 'reason', 'biological_groups')
    elif detail.key == 'baseline':
        keys = ('name', 'matched_answered_rows', 'mae_skill', 'mse_skill')
    else:
        return detail.text
    summary = '; '.join(key.replace('_', ' ') + ': ' + _text(value[key]) for key in keys if key in value)
    return summary or detail.text


def render_scorecard_detail(view: ScorecardView, key: str) -> str:
    """Escaped HTML for a metric definition or one full-test-record popup."""
    for metric in view.metrics:
        if metric.key == key:
            text = [metric.definition, 'Range: ' + metric.range, 'Chance: ' + metric.chance,
                    metric.reading, metric.denominator]
            return '<h3>' + escape(metric.label) + '</h3>' + ''.join(
                '<p>' + escape(part) + '</p>' for part in text if part)
    for detail in view.details:
        if detail.key == key:
            return '<h3>' + escape(detail.label) + '</h3><pre>' + escape(detail.text) + '</pre>'
    raise KeyError('Unknown scorecard detail: ' + key)


def _evidence_html(value):
    """Present supplied evidence as readable fields while keeping raw export exact."""
    if isinstance(value,dict):
        return '<table>'+''.join('<tr><td><b>'+escape(key.replace('_',' ').capitalize())+
            '</b></td><td>'+_evidence_html(item)+'</td></tr>' for key,item in value.items())+'</table>'
    if isinstance(value,list):
        return '<ul>'+''.join('<li>'+_evidence_html(item)+'</li>' for item in value)+'</ul>'
    if value is None:return 'Unavailable'
    if isinstance(value,bool):return 'Yes' if value else 'No'
    return escape(str(value)).replace('\n','<br>')


def render_scorecard_html(view: ScorecardView, *, expanded: bool = False) -> str:
    """Render compact result quality and one-action access to the full record.

    QTextBrowser hosts handle ``scorecard:expand`` by rendering with expanded=True,
    and metric/detail routes with :func:`render_scorecard_detail`.
    """
    html = ['<div class="scorecard"><h3>', escape(view.title), '</h3><p><b>',
            escape(view.kind.title()), '</b> · ', escape(view.status), '</p>']
    scope = dict(view.scope)
    html.append('<p>' + escape(' · '.join(scope[key] for key in ('organism', 'target', 'strategy'))) + '</p>')
    html.append('<p>Truth grade: ' + escape(view.source.grade) + '</p>')
    if dict(view.counts).get('eligible') == 1 and view.kind == 'performance':
        html.append('<p>Single recorded outcome; no per-gene accuracy probability.</p>')
    if view.kind == 'performance':
        selected = view.metrics if expanded else tuple(metric for metric in view.metrics if metric.headline)
        html.append('<table><tr><th>Result quality / coverage</th><th>Value</th><th>Population</th></tr>')
        for metric in selected:
            html.append('<tr><td>' + _anchor(metric.label, 'scorecard:metric/' + metric.key)
                        + '</td><td>' + escape(_value(metric.value))
                        + (' (' + escape(metric.status) + ')' if metric.value is None else '')
                        + '</td><td>' + escape(metric.denominator) + '</td></tr>')
        html.append('</table>')
    else:
        key = 'outcome' if view.kind == 'outcome' else 'evidence'
        detail=next(detail for detail in view.details if detail.key==key)
        try:value=json.loads(detail.text)
        except ValueError:value=detail.text
        html.append('<h4>'+escape(detail.label)+'</h4>'+_evidence_html(value))
    for key in ('sizes', 'uncertainty', 'calibration', 'baseline', 'freshness'):
        detail = next(detail for detail in view.details if detail.key == key)
        html.append('<p>' + _anchor(detail.label, 'scorecard:detail/' + key) + ': '
                    + escape(_compact_detail(detail)).replace('\n', '<br>') + '</p>')
    html.append('<p>' + _anchor('Full test record and definitions', 'scorecard:expand') + '</p>')
    html.append('<p>' + ' · '.join(_anchor(detail.label, 'scorecard:detail/' + detail.key)
                for detail in view.details if detail.key not in {'sizes', 'uncertainty', 'calibration', 'baseline', 'freshness'}) + '</p>')
    if view.links:
        html.append('<p>' + ' · '.join(_anchor(link.label, link.url) for link in view.links) + '</p>')
    if expanded:
        html.append('<h4>Evaluation scope</h4><pre>' + escape(_text(scope)) + '</pre>')
        html.extend(render_scorecard_detail(view, detail.key) for detail in view.details)
        html.extend(render_scorecard_detail(view, metric.key) for metric in view.metrics)
    html.append('</div>')
    return ''.join(html)


def export_scorecard(view: ScorecardView) -> str:
    """Export the exact immutable displayed model and supplied evaluation snapshot."""
    payload = asdict(view)
    payload['snapshot'] = json.loads(payload.pop('snapshot_json'))
    return json.dumps(payload, sort_keys=True, ensure_ascii=False, indent=2, allow_nan=False) + '\n'
