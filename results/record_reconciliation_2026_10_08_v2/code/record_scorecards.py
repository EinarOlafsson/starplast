"""Recompute scorecards from individual outcomes, keeping cohorts and tasks separate.

Legacy track records have incomplete truth/split lineage. Their arithmetic can
be reconciled without retroactively admitting their biological benchmarks.
Repeated settings, seeds and hold-out modes do not create independent samples.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json

import numpy as np
import pandas as pd

from . import capabilities as C, organisms as O, scorecard as SC, strategies as S
from .provenance import EVIDENCE_GRADES


@dataclass(frozen=True)
class RecordScope:
    """One evaluation cohort/settings/protocol address; missing lineage stays explicit."""

    organism: str
    strategy: str
    target: str
    task: str
    settings: str
    seed: int
    protocol: str
    partition: str
    benchmark_id: str
    truth_grade: str
    unit: str
    negative_semantics: str
    gaps: tuple[str, ...] = ()

    def __post_init__(self):
        if self.organism not in O.SPACES and self.organism not in O.HOST_TABLES:
            raise ValueError('Explicit organism required')
        if self.task not in SC.TASKS or type(self.seed) is not int:
            raise ValueError('Unknown task or seed')
        if not isinstance(self.settings, str) or not isinstance(self.gaps, tuple):
            raise ValueError('Settings and gaps require immutable explicit metadata')
        if self.unit not in {'gene', 'protein', 'pair', 'finding', 'module'}:
            raise ValueError('Explicit biological evaluation unit required')
        if self.truth_grade not in EVIDENCE_GRADES or self.unit != C.get(self.strategy).benchmark_unit:
            raise ValueError('Truth grade or strategy evaluation unit is not declared')
        if not all(isinstance(v, str) and v for v in (self.strategy, self.target, self.protocol,
                self.partition, self.benchmark_id, self.truth_grade, self.negative_semantics)):
            raise ValueError('Evaluation addresses and truth interpretation must be explicit')
        if any(v == 'unresolved' for v in (self.protocol, self.benchmark_id, self.truth_grade)) and not self.gaps:
            raise ValueError('Missing lineage requires an explicit gap')

    @property
    def identity(self):
        """Stable evaluation-scope identity, distinct from its entity cohort hash."""
        return hashlib.sha256(json.dumps(asdict(self), sort_keys=True).encode()).hexdigest()


def _finite(metrics):
    return {k: float(v) if v is not None and np.isfinite(v) else None for k, v in metrics.items()}


def _required(rows, columns):
    if set(columns) - set(rows.columns):
        raise ValueError('Missing outcome columns: ' + ', '.join(sorted(set(columns) - set(rows.columns))))


def _labels(rows):
    _required(rows, ('truth', 'prediction'))
    if not all(isinstance(x, str) and x for x in rows.truth):
        raise ValueError('Every eligible label outcome needs known explicit truth')
    if any(pd.notna(x) and not isinstance(x, str) for x in rows.prediction):
        raise ValueError('Categorical predictions must be labels or missing abstentions')
    if any(isinstance(x, str) and not x for x in rows.prediction):
        raise ValueError('Empty labels must be normalized to missing abstentions')
    answered = rows.prediction.map(lambda x: isinstance(x, str) and bool(x)).to_numpy()
    correct = np.array([isinstance(p, str) and p == t for p, t in zip(rows.prediction, rows.truth)], dtype=bool)
    if 'abstained' in rows and not np.array_equal(rows.abstained.to_numpy(dtype=bool), ~answered):
        raise ValueError('Recorded abstention contradicts the actual prediction')
    if 'correct' in rows:
        observed = rows.correct.fillna(False).to_numpy(dtype=bool)
        if not np.array_equal(observed, correct):
            raise ValueError('Recorded correctness contradicts prediction and truth')
    return answered, correct


def _binary(rows, column, *, unknown=False):
    _required(rows, (column,))
    values = rows[column]
    if any(not isinstance(x, (bool, np.bool_)) and not (unknown and pd.isna(x)) for x in values):
        raise ValueError('Binary truth/decisions must be explicit booleans; unknown is not negative')
    return values


def _interval(rows, answered, correct, seed):
    if 'group' not in rows or rows.group.isna().any() or not rows.group.astype(str).str.len().gt(0).all():
        return {'status': 'unavailable', 'reason': 'Biological groups unresolved', 'biological_groups': None}
    groups = rows.group.astype(str).to_numpy()
    unique = np.unique(groups)
    if len(unique) < 5:
        return {'status': 'unavailable', 'reason': 'Fewer than five biological groups', 'biological_groups': len(unique)}
    indices = [np.flatnonzero(groups == group) for group in unique]
    rng = np.random.default_rng(seed)
    all_rates, call_rates = [], []
    for _ in range(500):
        selected = np.concatenate([indices[i] for i in rng.integers(0, len(indices), len(indices))])
        all_rates.append(float(correct[selected].mean()))
        n = int(answered[selected].sum())
        if n:
            call_rates.append(float(correct[selected].sum() / n))
    return {'status': 'descriptive_group_bootstrap', 'biological_groups': len(unique), 'draws': 500,
            'accuracy_all_hidden': np.quantile(all_rates, [0.025, 0.975]).tolist(),
            'accuracy_among_calls': np.quantile(call_rates, [0.025, 0.975]).tolist() if call_rates else None,
            'limits': 'Groups supplied by caller; not an individual-gene probability or source-independence proof'}


def aggregate(rows, scope, *, parameters=None):
    """Build a card from one complete unique cohort, using the standard metric authority.

    Known truth is required where a metric identifies correctness. Positive-only
    rankings/retrieval do not manufacture negative labels or precision. Cluster
    geometry and numerical estimates have no categorical correct/wrong count.
    """
    if not isinstance(scope, RecordScope):
        raise TypeError('A typed evaluation scope is required')
    rows = pd.DataFrame(rows).copy().reset_index(drop=True)
    _required(rows, ('entity',))
    if rows.entity.isna().any() or not rows.entity.astype(str).str.len().gt(0).all() or rows.entity.astype(str).duplicated().any():
        raise ValueError('Each eligible entity must appear exactly once in this cohort')
    parameters = dict(parameters or {})
    n = len(rows)
    if not n:
        raise ValueError('An empty cohort requires an unavailable card, not a numeric accuracy')
    answered, correct = np.zeros(n, dtype=bool), None
    metrics, extra = {}, {}
    task = scope.task
    if task == SC.T_LABEL:
        answered, correct = _labels(rows)
        scores = parameters.get('class_scores')
        if scores is not None and not pd.DataFrame(scores).index.equals(rows.index):
            raise ValueError('Class-score order must match the exact cohort')
        predictions = rows.prediction.astype(object).where(rows.prediction.notna(), None)
        metrics = SC.label_calls(predictions, rows.truth, np.arange(n), class_scores=scores)
        extra['confusion'] = [{'truth': str(truth), 'prediction': str(pred) if pd.notna(pred) else None, 'count': int(count)}
            for (truth, pred), count in rows.groupby(['truth', 'prediction'], dropna=False, observed=True).size().items()]
    elif task == SC.T_VALUES:
        _required(rows, ('truth', 'prediction'))
        truth, pred = rows.truth.to_numpy(dtype=float), rows.prediction.to_numpy(dtype=float)
        if not np.isfinite(truth).all():
            raise ValueError('Unmeasured numeric truth cannot enter an eligible cohort')
        answered = np.isfinite(pred)
        metrics = SC.values(pred, truth)
        extra['quantity_unit'] = parameters.get('quantity_unit', 'unresolved')
    elif task in {SC.T_RANK, SC.T_SET}:
        truth = _binary(rows, 'positive', unknown=True)
        verified = parameters.get('verified_negatives')
        if type(verified) is not bool:
            raise ValueError('Declare whether negatives are verified')
        if not verified and (truth.dropna() == False).any():
            raise ValueError('Unverified absence must be unknown, not false')
        if task == SC.T_RANK:
            _required(rows, ('score',))
            scores = rows.score.to_numpy(dtype=float)
            answered = np.isfinite(scores)
            if verified:
                if truth.isna().any():
                    raise ValueError('Standard ranking requires observed truth for the whole candidate universe')
                metrics = SC.ranking(scores, truth.to_numpy(dtype=bool))
            else:
                depth = parameters.get('depth')
                if type(depth) is not int or not 0 < depth <= n:
                    raise ValueError('Positive-only recovery needs a frozen valid depth')
                order = np.argsort(-np.where(answered, scores, -np.inf), kind='stable')
                positive = truth.fillna(False).to_numpy(dtype=bool)
                hits = int(positive[order[:depth]].sum())
                metrics = {k: None for k in SC.TASKS[task].metrics}
                extra.update(known_positive_recall_at_depth=hits / int(positive.sum()) if positive.any() else None,
                             depth=depth, precision_identifiable=False, negative_assumptions='unknown_is_not_negative')
        else:
            returned = _binary(rows, 'returned').to_numpy(dtype=bool)
            answered = returned
            positive = truth.fillna(False).to_numpy(dtype=bool)
            if verified:
                if truth.isna().any():
                    raise ValueError('Standard set retrieval requires observed candidate membership')
                metrics = SC.set_retrieval(np.flatnonzero(returned), np.flatnonzero(positive), np.arange(n))
            else:
                metrics = {k: None for k in SC.TASKS[task].metrics}
                metrics['returned'] = int(returned.sum())
                extra.update(known_positive_recall=int((returned & positive).sum()) / int(positive.sum()) if positive.any() else None,
                             precision_identifiable=False, negative_assumptions='unknown_is_not_negative')
    elif task == SC.T_CLUSTER:
        _required(rows, ('truth', 'cluster'))
        if not all(isinstance(x, str) and x for x in rows.truth):
            raise ValueError('Cluster recovery requires known biological labels')
        chosen = parameters.get('chosen_clusters')
        if not isinstance(chosen, dict):
            raise ValueError('Freeze label-to-cluster selection from inner training before testing')
        raw = rows.cluster.to_numpy(dtype=float)
        if not np.isfinite(raw).all() or not np.equal(raw, np.floor(raw)).all():
            raise ValueError('Cluster IDs must be finite integers; no silent truncation')
        labels = raw.astype(int)
        answered = labels != -1
        metrics = SC.cluster_recovery(labels, rows.truth.to_numpy(), chosen)
        extra.update(placed=int(answered.sum()), unplaced=int((~answered).sum()), correct_wrong_interpretation='Unavailable: geometry is not a categorical caller')
    else:
        replicated = _binary(rows, 'replicated').to_numpy(dtype=bool)
        answered = np.ones(n, dtype=bool)
        metrics = SC.replication(int(replicated.sum()), n, parameters.get('null_rates'))
        extra.update(tested_findings=n, replicated_findings=int(replicated.sum()), unreplicated_findings=int((~replicated).sum()))
    called = int(answered.sum())
    counts = {'eligible': n, 'answered': called, 'abstained': n - called,
              'correct': int(correct.sum()) if correct is not None else None,
              'wrong': called - int(correct.sum()) if correct is not None else None,
              'unique_biological_entities': int(rows.entity.nunique())}
    if correct is not None:
        extra.update(accuracy_all_hidden=float(correct.mean()), accuracy_among_calls=float(correct.sum() / called) if called else None,
            uncertainty=_interval(rows, answered, correct, scope.seed), single_outcome=n == 1)
    else:
        extra.setdefault('correct_wrong_interpretation', 'Not a categorical prediction task; use the declared task metrics')
    clean = _finite(metrics)
    return {'scope': asdict(scope), 'scope_identity': scope.identity,
            'cohort_identity': hashlib.sha256(json.dumps(rows.entity.tolist()).encode()).hexdigest(),
            'counts': counts, 'metrics': clean, 'metric_status': {k: 'available' if v is not None else 'unavailable_for_this_cohort_or_truth' for k, v in clean.items()},
            'small_sample': n < 5, 'extra': extra}


def class_cards(rows, scope):
    """Class precision includes false calls from other classes; recall includes abstentions."""
    if scope.task != SC.T_LABEL:
        raise ValueError('Class call metrics require label predictions')
    rows = pd.DataFrame(rows).reset_index(drop=True)
    _labels(rows)
    result = []
    for label in sorted(set(rows.truth)):
        actual, predicted = rows.truth == label, (rows.prediction == label).fillna(False)
        tp, fp, fn = int((actual & predicted).sum()), int((~actual & predicted).sum()), int((actual & ~predicted).sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn)
        card = aggregate(rows[actual], scope)
        card['class'] = label
        card['class_metrics'] = {'precision': precision, 'recall': recall,
            'f1': 2 * precision * recall / (precision + recall) if precision + recall else 0.0,
            'true_positive': tp, 'false_positive': fp, 'false_negative': fn,
            'evaluation_prevalence': int(actual.sum()) / len(rows)}
        result.append(card)
    return result


def views(rows, scope, *, parameters=None):
    """Reuse one record set across levels without averaging unrelated tasks or component accuracy."""
    card = aggregate(rows, scope, parameters=parameters)
    rows = pd.DataFrame(rows)
    method = []
    for technique in S.get(scope.strategy).techniques:
        method.append({'technique': technique, 'pipeline_card': card,
                       'component_accuracy': None, 'status': 'requires_direct_or_ablation_test'})
    organism = {'organism': scope.organism, 'task': scope.task, 'cards': [card],
                'unique_biological_entities': int(rows.entity.nunique()), 'pooled_accuracy': None}
    return {'target': card, 'strategy': card, 'method': method, 'organism': organism,
            'class': class_cards(rows, scope) if scope.task == SC.T_LABEL else [],
            'gene': [{'entity': str(row.entity), 'outcome': row.to_dict(), 'validation_scope': scope.identity,
                      'per_gene_accuracy_probability': None} for _, row in rows.iterrows()]}


def legacy_cohorts(ledger):
    """Yield distinct legacy cohort rows with unresolved protocol/truth lineage.

    Every setting, seed, hold-out mode and named set remains separate. Folds are
    disjoint parts of a together cohort, never additional independent samples.
    """
    columns = ['organism', 'strategy', 'target', 'setting_key', 'seed', 'mode', 'set_name']
    _required(ledger, [*columns, 'gene_id', 'truth', 'prediction', 'correct', 'abstained'])
    for keys, part in ledger.groupby(columns, observed=True, dropna=False, sort=False):
        organism, strategy, target, settings, seed, mode, set_name = keys
        name = '' if pd.isna(set_name) else str(set_name)
        scope = RecordScope(str(organism), str(strategy), str(target), SC.T_LABEL, str(settings), int(seed),
            'unresolved', str(mode) + ':' + name, 'unresolved', 'unresolved', 'gene',
            'categorical labels; not biological presence/absence',
            ('Legacy file does not pin independent truth, source exclusions or nested fit/calibration lineage',))
        rows = part.rename(columns={'gene_id': 'entity'})
        if rows.entity.duplicated().any():
            raise ValueError('Repeated legacy evaluation entities cannot be counted as independent samples')
        yield scope, rows


def organism_overview(batches):
    """Retain per-scope cards and distinct entity counts, never a cross-task accuracy."""
    cards, entities, scopes = [], {}, set()
    for scope, rows, parameters in batches:
        if scope.identity in scopes:
            raise ValueError('Duplicate evaluation scope; repeated results are not independent samples')
        scopes.add(scope.identity)
        card = aggregate(rows, scope, parameters=parameters)
        cards.append(card)
        key = (scope.organism, scope.unit)
        entities.setdefault(key, set()).update(pd.DataFrame(rows).entity.astype(str))
    return {'cards': cards, 'unique_entities_by_organism_and_unit': [
        {'organism': organism, 'unit': unit, 'entities': len(values)}
        for (organism, unit), values in sorted(entities.items())],
        'pooled_accuracy': None, 'limits': 'Unique entities are not independent studies; seeds/settings/protocols retain separate cards'}
