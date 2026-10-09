"""Descriptive comparisons of aligned, frozen source-profile recovery tests.

These counts do not estimate biological accuracy, independent evidence or
confidence for an unknown gene. No model is fitted or selected here.
"""
from __future__ import annotations

from . import organisms as _ORGANISMS

from collections import Counter
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from . import functional_results as F
from .scorecard_view import DetailView, MetricView, ScorecardView, SourceView


_OUTCOMES = (
    'both_abstain', 'left_only_correct', 'left_only_wrong',
    'right_only_correct', 'right_only_wrong', 'agree_correct', 'agree_wrong',
    'conflict_left_correct', 'conflict_right_correct', 'conflict_both_wrong',
)
_LIMIT = ('Recorded annotation recovery only. Biological accuracy and calibrated '
          'confidence are unavailable. Shared inputs prevent an independence claim. '
          'Group-aware uncertainty is not estimated in this report. '
          'This outer test cannot select an ensemble or tune a decision rule.')
REPORT_SHA256 = '940da900b3d5ba356c6e360f220787b849d7382138132ccc96f2d1684178e98b'


def find_pair(selected, benchmarks):
    """Return the two currently supported methods only at the selected address."""
    if selected is None:
        return None
    if (selected.organism, selected.namespace, selected.metadata['target']) != (
            _ORGANISMS.FALCIPARUM, 'ec_major', 'ec_direct_complete_major_profile'):
        return None
    pair = []
    for strategy in ('feature_knn', 'random_forest'):
        matches = [b for b in benchmarks if b.organism == selected.organism
                   and b.namespace == selected.namespace and b.metadata['target'] == selected.metadata['target']
                   and b.metadata['strategy'] == strategy]
        if len(matches) != 1:
            return None
        pair.append(matches[0])
    return tuple(pair)


def load_report(left, right, context, *, path=None, expected_sha256=REPORT_SHA256):
    """Load the externally pinned report and refuse stale artifacts or source tables."""
    path = Path(path) if path is not None else Path(__file__).with_name('data') / 'functional_agreement.json'
    encoded = path.read_bytes()
    if not expected_sha256 or hashlib.sha256(encoded).hexdigest() != expected_sha256:
        raise ValueError('Paired report checksum changed')
    report = json.loads(encoded)
    for benchmark in (left, right):
        F.require_context(benchmark, context)
    derived = compare(left, right)
    for key, value in derived.items():
        if key == 'overlap':
            for name, original in value.items():
                if name != 'selected_feature_overlap' and report[key].get(name) != original:
                    raise ValueError('Paired report dependency scope changed: ' + name)
        elif report.get(key) != value:
            raise ValueError('Paired report differs from current artifacts: ' + key)
    build_scorecard(report)
    return report


def paired_rows(left, right):
    """Partition every matched test gene, refusing mismatched populations or truth."""
    if not left or len(left) != len(right):
        raise ValueError('Paired rows require the same nonempty population')
    seen = set()
    result = []
    for a, b in zip(left, right):
        for key in ('entity', 'truth', 'group', 'training_supported'):
            if key not in a or key not in b or a[key] != b[key]:
                raise ValueError('Paired row alignment changed: ' + key)
        if (not isinstance(a['entity'], str) or not a['entity'] or a['entity'] in seen
                or not isinstance(a['group'], str) or not a['group']
                or any(type(row['training_supported']) is not bool for row in (a, b))):
            raise ValueError('Paired entity/group/support must remain explicit and unique')
        seen.add(a['entity'])
        truth = a['truth']
        if not isinstance(truth, str) or not truth:
            raise ValueError('Paired truth must be known')
        for row in (a, b):
            if ('prediction' not in row or type(row.get('abstained')) is not bool
                    or row['abstained'] != (row['prediction'] is None)
                    or row.get('calibrated_confidence') is not None
                    or row['prediction'] is not None and not isinstance(row['prediction'], str)):
                raise ValueError('Paired prediction/abstention/confidence semantics changed')
        x, y = a['prediction'], b['prediction']
        if x is None and y is None: outcome = 'both_abstain'
        elif y is None: outcome = 'left_only_correct' if x == truth else 'left_only_wrong'
        elif x is None: outcome = 'right_only_correct' if y == truth else 'right_only_wrong'
        elif x == y: outcome = 'agree_correct' if x == truth else 'agree_wrong'
        elif x == truth: outcome = 'conflict_left_correct'
        elif y == truth: outcome = 'conflict_right_correct'
        else: outcome = 'conflict_both_wrong'
        result.append({'entity': a['entity'], 'truth': truth, 'group': a['group'],
                       'training_supported': a['training_supported'],
                       'left_prediction': x, 'right_prediction': y, 'outcome': outcome})
    return result


def _summary(rows):
    counts = {key: 0 for key in _OUTCOMES}
    counts.update(Counter(row['outcome'] for row in rows))
    agreements = counts['agree_correct'] + counts['agree_wrong']
    conflicts = sum(counts[key] for key in _OUTCOMES if key.startswith('conflict_'))
    joint = agreements + conflicts
    counts.update(eligible=len(rows), joint_answered=joint, agreements=agreements,
                  conflicts=conflicts, shared_wrong=counts['agree_wrong'] + counts['conflict_both_wrong'],
                  unsupported=sum(not row['training_supported'] for row in rows),
                  groups=len({row['group'] for row in rows}))
    def rate(numerator, denominator):
        return {'numerator': numerator, 'denominator': denominator,
                'value': numerator / denominator if denominator else None}
    return {'counts': counts, 'rates': {
        'joint_coverage': rate(joint, len(rows)),
        'agreement_among_joint_calls': rate(agreements, joint),
        'source_recovery_among_agreements': rate(counts['agree_correct'], agreements),
        'correct_agreements_all_eligible': rate(counts['agree_correct'], len(rows)),
        'shared_wrong_among_joint_calls': rate(counts['shared_wrong'], joint)}}


def compare(left: F.FunctionalBenchmark, right: F.FunctionalBenchmark):
    """Revalidate artifacts and compare only exactly aligned, explicit scopes.

    Native EC artifacts currently supply the required context, target/cohort and
    prepared-packet identities. Other benchmarks without these fields are refused
    rather than aligned by label name alone.
    """
    left, right = (F._benchmark(deepcopy({'metadata': b.metadata, 'payloads': b.payloads}))
                   for b in (left, right))
    specs = [b.metadata['source_manifest']['spec'] for b in (left, right)]
    if left.metadata['source_artifact_identity'] == right.metadata['source_artifact_identity']:
        raise ValueError('A paired report requires distinct artifacts')
    for key in ('organism', 'target', 'truth_grade', 'source_targets', 'source_context'):
        if left.metadata[key] != right.metadata[key]:
            raise ValueError('Paired source scope changed: ' + key)
    if left.namespace != right.namespace:
        raise ValueError('Paired profile namespace changed')
    for key in ('organism', 'target', 'task', 'seed', 'protocol', 'partition',
                'benchmark_id', 'truth_grade', 'unit', 'negative_semantics'):
        if left.card['scope'][key] != right.card['scope'][key]:
            raise ValueError('Paired evaluation scope changed: ' + key)
    for key in ('query_json', 'entity_order', 'fit_entities', 'fit_role', 'role',
                'partition_id', 'feature_access', 'task', 'output',
                'calibration_scope', 'evaluation_partition', 'confidence_kind'):
        if specs[0][key] != specs[1][key]:
            raise ValueError('Paired artifact scope changed: ' + key)
    if any(spec['fit_role'] != 'train' for spec in specs):
        raise ValueError('Paired fitted state requires training-only roles')
    evaluations = [json.loads(spec['evaluation_scope_json']) for spec in specs]
    keys = ('unit', 'eligible_population', 'truth_grade', 'negative_semantics',
            'target_identity', 'cohort_identity', 'prepared_packet_manifest_sha256',
            'context', 'biological_admission', 'calibration', 'deployment')
    for key in keys:
        if (key not in evaluations[0] or key not in evaluations[1]
                or evaluations[0][key] != evaluations[1][key]):
            raise ValueError('Paired explicit evaluation context changed: ' + key)
    if any(not evaluations[0][key] for key in ('context', 'target_identity', 'cohort_identity',
                                              'prepared_packet_manifest_sha256')):
        raise ValueError('Paired context and source identities must be explicit')
    if (evaluations[0]['biological_admission'] is not False
            or evaluations[0]['calibration'] != 'unavailable'
            or evaluations[0]['deployment'] != 'unknown'):
        raise ValueError('Paired source recovery cannot claim biological admission, calibration or deployment')
    dependencies = []
    for spec in specs:
        mapping = {}
        for dep in spec['dependencies']:
            key = (dep['kind'], dep['name'])
            if key in mapping:
                raise ValueError('Duplicate paired dependency address')
            mapping[key] = dep['sha256']
        dependencies.append(mapping)
    common = set(dependencies[0]) & set(dependencies[1])
    shared = [{'kind': kind, 'name': name, 'sha256': dependencies[0][kind, name]}
              for kind, name in sorted(common)
              if dependencies[0][kind, name] == dependencies[1][kind, name]]
    changed = [{'kind': kind, 'name': name, 'left_sha256': dependencies[0][kind, name],
                'right_sha256': dependencies[1][kind, name]}
               for kind, name in sorted(common)
               if dependencies[0][kind, name] != dependencies[1][kind, name]]
    rows = paired_rows(left.rows, right.rows)
    return {'schema_version': 1, 'interpretation': _LIMIT,
            'scope': {key: deepcopy(evaluations[0][key]) for key in keys},
            'organism': left.organism, 'target': left.metadata['target'],
            'namespace': left.namespace, 'protocol': left.card['scope']['protocol'],
            'methods': [{'strategy': b.metadata['strategy'],
                         'artifact_identity': b.metadata['source_artifact_identity'],
                         'settings': spec['settings_json'], 'counts': deepcopy(b.card['counts'])}
                        for b, spec in zip((left, right), specs)],
            'overlap': {'shared_dependencies': shared, 'changed_dependencies': changed,
                        'same_training_entities': True, 'same_truth_and_groups': True,
                        'independent_evidence': False,
                        'selected_feature_overlap': 'Not measured by this compact artifact comparison',
                        'unmatched_dependency_addresses': [
                            [{'kind': kind, 'name': name, 'sha256': mapping[kind, name]}
                             for kind, name in sorted(set(mapping) - common)]
                            for mapping in dependencies]},
            **_summary(rows), 'rows': rows,
            'group_uncertainty': {'status': 'unavailable', 'reason': 'No group-aware interval estimated'},
            'biological_accuracy': None, 'calibrated_confidence': None,
            'biological_admission': False, 'ensemble_selection': False}


def build_scorecard(report) -> ScorecardView:
    """Freeze the supplied paired report into the existing card/export interface."""
    if (report.get('biological_admission') is not False
            or report.get('biological_accuracy') is not None
            or report.get('calibrated_confidence') is not None
            or report.get('ensemble_selection') is not False
            or report.get('overlap', {}).get('independent_evidence') is not False):
        raise ValueError('Paired source recovery cannot advertise confidence or independence')
    native_rows = []
    for side in ('left', 'right'):
        native_rows.append([{**{key: row[key] for key in ('entity', 'truth', 'group', 'training_supported')},
                             'prediction': row[side + '_prediction'],
                             'abstained': row[side + '_prediction'] is None,
                             'calibrated_confidence': None} for row in report['rows']])
    verified_rows = paired_rows(*native_rows)
    if verified_rows != report['rows']:
        raise ValueError('Paired scorecard outcomes differ from original calls')
    recomputed = _summary(verified_rows)
    if any(report[key] != recomputed[key] for key in ('counts', 'rates')):
        raise ValueError('Paired scorecard counts differ from row outcomes')
    labels = {'joint_coverage': 'Both methods answered',
              'agreement_among_joint_calls': 'Agreement among joint calls',
              'source_recovery_among_agreements': 'Source recovery among agreements',
              'correct_agreements_all_eligible': 'Correct agreements / all test genes',
              'shared_wrong_among_joint_calls': 'Both wrong among joint calls'}
    metrics = tuple(MetricView(key, labels[key], rate['value'],
                    'available' if rate['value'] is not None else 'unavailable',
                    labels[key] + '; complete reference profiles compared exactly.',
                    _LIMIT, '0 to 1', 'No independent null supplied',
                    f"{rate['numerator']} / {rate['denominator']}", True)
                    for key, rate in report['rates'].items())
    details = tuple(DetailView(key, label, F._canonical(report[key])) for key, label in (
        ('methods', 'Original method scopes and counts'), ('overlap', 'Shared inputs and dependencies'),
        ('group_uncertainty', 'Uncertainty'), ('rows', 'All paired gene outcomes')))
    content = {
        'evidence': ('Paired source-profile outcomes', {
            labels[key]: (f"{rate['numerator']} / {rate['denominator']} ({rate['value']:.1%})"
                          if rate['value'] is not None else 'Unavailable: no eligible denominator')
            for key, rate in report['rates'].items()}),
        'sizes': ('Sample sizes', report['counts']),
        'uncertainty': ('Uncertainty', report['group_uncertainty']),
        'calibration': ('Calibration', {'status': 'unavailable', 'reason': _LIMIT}),
        'baseline': ('Matched baseline', {'status': 'not estimated',
                     'reason': 'Original per-method control cards remain separate; no consensus control is fitted'}),
        'freshness': ('Frozen artifact identities', report['methods']),
        'gaps': ('Limitations', _LIMIT),
    }
    details += tuple(DetailView(key, label, F._canonical(value))
                     for key, (label, value) in content.items())
    return ScorecardView('evidence', 'Paired source-profile recovery', 'descriptive', None,
        tuple((key, str(report[key])) for key in ('organism', 'target', 'namespace', 'protocol'))
        + (('strategy', ' / '.join(method['strategy'] for method in report['methods'])),),
        tuple(report['counts'].items()),
        SourceView(report['scope']['truth_grade'], 'Recorded direct EC annotations', 'unresolved',
                   report['scope']['target_identity'], F._canonical(report['scope']['context']),
                   _LIMIT, report['scope']['negative_semantics']),
        metrics, details, (), F._canonical(report))


def build_gene_scorecard(report, entity) -> ScorecardView:
    """Present one observed paired outcome without a per-gene accuracy probability."""
    from dataclasses import replace
    view = build_scorecard(report)
    matches = [row for row in report['rows'] if row['entity'] == entity]
    if len(matches) != 1:
        raise ValueError('Gene is outside the frozen paired test')
    row = matches[0]
    snapshot = {'outcome': row, 'methods': report['methods'],
                'scope': report['scope'], 'paired_cohort_counts': report['counts'],
                'biological_accuracy': None, 'calibrated_confidence': None}
    details = tuple(detail for detail in view.details if detail.key not in ('rows', 'sizes', 'evidence'))
    details += (DetailView('outcome', 'Observed paired test outcome', F._canonical(row)),
                DetailView('rows', 'Selected gene outcome', F._canonical([row])),
                DetailView('sizes', 'One gene from the frozen paired cohort',
                           F._canonical({'eligible': 1, 'paired_cohort': report['counts']})))
    return replace(view, kind='outcome', title='Paired held-out outcome: ' + entity,
                   status='Observed source-profile outcome; no per-gene confidence probability',
                   counts=(('eligible', 1),), metrics=(), details=details,
                   snapshot_json=F._canonical(snapshot))
