"""Present frozen source-profile classes with explicit evaluation populations.

This adapter reads supplied, verified records only. A profile-member aggregate
and a whole-cohort class precision answer different questions; their displayed
values and denominators remain separate. No fitting or biological admission is
performed here.
"""
from __future__ import annotations

from dataclasses import replace

from . import functional_results as F, scorecard as SC
from .scorecard_view import MetricView, ScorecardLink, ScorecardView, build_scorecard_view


_LIMIT = ('Reference-annotation recovery; biological precision and biological recall unavailable. '
          'Calibrated confidence unavailable.')
_COUNT_KEYS = ('eligible', 'known_positive_genes', 'true_positive', 'false_positive', 'false_negative')


def _class_metric(key, label, value, definition, denominator):
    return MetricView(key, label, value, 'available' if value is not None else 'unavailable',
        definition, 'This measures recovery of recorded source annotations; independent biological '
        'activity accuracy and per-gene calibrated confidence are unavailable.', '0 to 1',
        'No matched class-specific null supplied', denominator, True)


def build(benchmark: F.FunctionalBenchmark, address: str) -> ScorecardView:
    """Return an immutable shared card for one recorded profile or major class.

    ``address`` is ``profile:<stored profile>`` or ``major:<stored major class>``.
    Callers obtain ``benchmark`` from the existing verified bundle loader. Exact
    native cards are retained in the export snapshot; supplementary presentation
    labels do not change the global metric or algorithm schemas. Retained rows
    include reference members, abstentions and false class calls from other
    profiles. Their displayed row count is not substituted for a metric population.
    """
    if not isinstance(benchmark, F.FunctionalBenchmark):
        raise TypeError('A verified FunctionalBenchmark is required')
    if not isinstance(address, str) or ':' not in address:
        raise ValueError('An explicit profile or major class address is required')
    kind, value = address.split(':', 1)
    if kind not in {'profile', 'major'} or not value:
        raise ValueError('An explicit profile or major class address is required')
    metadata = benchmark.metadata
    if metadata.get('biological_admission') is not False or metadata.get('calibrated_confidence') is not None:
        raise ValueError('Class presentation cannot assert biological admission or calibrated confidence')
    candidates = benchmark.profile_class_cards if kind == 'profile' else benchmark.major_class_cards
    matching = [card for card in candidates if card['class'] == value]
    if len(matching) != 1:
        raise ValueError('Class address must identify exactly one recorded card')
    native = matching[0]
    full_scope = benchmark.card['scope']
    full_n = benchmark.card['counts']['eligible']
    if kind == 'profile':
        F.profile_classes(value)
        actual = lambda row: row['truth'] == value
        predicted = lambda row: row['prediction'] == value
        classes = native['class_metrics']
        counts = native['counts']
        tp, fp, fn = (classes[key] for key in ('true_positive', 'false_positive', 'false_negative'))
        prevalence = classes['evaluation_prevalence']
        presentation = native
        title = 'Recorded reference profile: ' + F.class_title(value)
        convention = 'The supplied profile-card convention reports zero when no class calls were made.'
    else:
        actual = lambda row: value in F.profile_classes(row['truth'])
        predicted = lambda row: value in (F.profile_classes(row['prediction']) or ())
        classes = native
        counts = {key: native[key] for key in _COUNT_KEYS}
        tp, fp, fn = (native[key] for key in ('true_positive', 'false_positive', 'false_negative'))
        prevalence = native['reference_prevalence']
        # Standard label accuracy does not describe a binary major-class card.
        # Retain its supplied overall profile-call coverage and custom class rates.
        presentation = {'scope': full_scope, 'counts': counts,
                        'metrics': {'coverage': native['coverage']}, 'native_class_card': native}
        title = 'Recorded reference class: ' + F.class_title(value)
        convention = 'The supplied major-class card reports null when no class calls were made.'
    rows = [row for row in benchmark.rows if actual(row) or predicted(row)]
    failures = [row for row in rows if actual(row) != predicted(row)]
    class_sizes = {key: classes[key] for key in ('true_positive', 'false_positive', 'false_negative')}
    class_sizes.update(whole_test_cohort=full_n, reference_members=tp + fn,
                       whole_cohort_class_calls=tp + fp, retained_outcome_rows=len(rows))
    spec = metadata.get('source_manifest', {}).get('spec', {})
    source_references = [entry for entry in spec.get('dependencies', ())
                         if entry.get('kind') == 'table' and entry.get('name') == 'installed_nodes']
    if len(source_references) > 1:
        raise ValueError('Installed source-table identity must be unique')
    source_sha256 = source_references[0]['sha256'] if source_references else None
    summary = metadata.get('summary', {})
    baseline_controls = {name: {'counts': card['counts'], 'metrics': card['metrics']}
                         for name, card in benchmark.baseline_cards.items()}
    view = build_scorecard_view(presentation, task=SC.T_LABEL, title=title, status=_LIMIT,
        details={
            'source': {'name': ', '.join(metadata['source_targets']), 'grade': metadata['truth_grade'],
                'sha256': source_sha256, 'context': metadata['source_context'],
                'lineage': 'Frozen recorded reference profiles; independent biological activity not admitted',
                'negative_semantics': 'Missing or unresolved source annotation is unknown, not biological absence'},
            'split': {'scope': full_scope, 'roles': {key: summary.get(key)
                for key in ('train', 'tune', 'calibration', 'test')}, 'fit_entities': spec.get('fit_entities')},
            'sizes': {**counts, 'class_metric_populations': class_sizes},
            'baseline': {'status': 'unavailable', 'reason': 'No matched class-specific baseline card supplied'},
            'controls': {'scope': 'Whole-cohort controls; not a matched class-specific baseline',
                         'cards': baseline_controls},
            'rows': rows, 'failures': failures,
            'freshness': {'artifact_identity': metadata['source_artifact_identity'],
                          'source_table_sha256': source_sha256},
            'gaps': list(full_scope.get('gaps', ())) + [_LIMIT,
                'Class precision includes false calls from the entire frozen test cohort. '
                'Retained rows include reference members and all calls of the selected class; '
                'the retained-row count is not the accuracy or prevalence denominator.'],
            'native_class_card': native, 'class_address': address,
        }, links=(ScorecardLink('Retained class outcomes', 'scorecard:rows/selected'),))
    metrics = []
    for metric in view.metrics:
        if kind == 'major' and metric.key != 'coverage':
            continue
        if metric.key == 'accuracy':
            metric = replace(metric, label='Reference-member accuracy, all members',
                definition='Correct complete-profile calls among genes whose recorded profile is this class, '
                    'divided by all recorded members of this profile. Abstentions count as errors.',
                denominator='Reference-profile members: ' + str(counts['eligible']))
        elif metric.key == 'precision_of_calls':
            metric = replace(metric, label='Reference-member accuracy among calls',
                definition='Correct complete-profile calls divided by calls made among recorded members '
                    'of this profile. False calls of this profile from other reference profiles are '
                    'excluded here and included in the separate whole-cohort class precision.',
                denominator='Calls among reference-profile members: ' + str(counts['answered']))
        elif metric.key == 'coverage':
            metric = replace(metric,
                label='Reference-member call coverage' if kind == 'profile' else 'Profile-call coverage, whole cohort',
                definition=('Any complete-profile calls divided by all recorded members of this profile.'
                    if kind == 'profile' else 'Any complete-profile calls divided by all frozen test genes; '
                    'this is overall caller reach, not the share called this major class.'),
                denominator=('Reference-profile members: ' + str(counts['eligible']) if kind == 'profile'
                    else 'All frozen test genes: ' + str(full_n)))
        metrics.append(metric)
    metrics.extend((
        _class_metric('class_precision', 'Class precision, whole cohort', classes['precision'],
            'TP / (TP + FP) for calls of this recorded class across the entire frozen test cohort. '
            'Calls from other reference profiles count as false positives. ' + convention,
            'Whole-cohort calls of this class: ' + str(tp + fp)),
        _class_metric('class_recall', 'Class recall, all reference members', classes['recall'],
            'TP / (TP + FN) among all recorded members of this class in the frozen test cohort. '
            'Missed calls and abstentions count as false negatives.',
            'Reference members in full cohort: ' + str(tp + fn)),
        _class_metric('class_f1', 'Class F1, whole cohort', classes['f1'],
            'The supplied class F1 balances class precision and recall using TP, FP and FN from '
            'the whole frozen test cohort: 2 * TP / (2 * TP + FP + FN) where defined.',
            'Whole-cohort class counts: 2 × TP + FP + FN'),
        _class_metric('reference_prevalence', 'Reference prevalence, whole cohort', prevalence,
            'Recorded members of this class divided by all frozen test genes. '
            'This describes source annotation prevalence, not biological presence.',
            'All frozen test genes: ' + str(full_n)),
    ))
    return replace(view, metrics=tuple(metrics))
