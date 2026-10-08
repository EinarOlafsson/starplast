"""Metric denominator reconciliation and cohort/task separation."""
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from starplast import organisms as O, record_scorecards as R, scorecard as SC


def scope(key='feature_knn', task=SC.T_LABEL, unit='gene'):
    return R.RecordScope(O.TOXOPLASMA, key, 'fixture', task, 'frozen_fixture', 3,
        'fixture_protocol', 'outer_test', 'fixture_benchmark', 'synthetic_control', unit, 'fixture truth only')


def labels():
    return pd.DataFrame({'entity': ['g1', 'g2', 'g3', 'g4'], 'truth': ['A', 'A', 'B', 'B'],
        'prediction': ['A', None, 'A', 'B']})


def test_all_hidden_and_answered_accuracy_match_standard_card():
    rows = labels()
    card = R.aggregate(rows, scope())
    assert card['counts'] == {'eligible': 4, 'answered': 3, 'abstained': 1, 'correct': 2, 'wrong': 1, 'unique_biological_entities': 4}
    assert card['extra']['accuracy_all_hidden'] == card['metrics']['accuracy'] == .5
    assert card['extra']['accuracy_among_calls'] == card['metrics']['precision_of_calls'] == pytest.approx(2 / 3)
    assert card['extra']['uncertainty']['reason'] == 'Biological groups unresolved'
    assert card['small_sample']
    assert card['metrics']['macro_auroc'] is None
    classes = R.class_cards(rows, scope())
    a = next(c for c in classes if c['class'] == 'A')
    assert a['class_metrics']['precision'] == .5, 'B called A must remain a false positive'
    assert a['class_metrics']['recall'] == .5, 'abstention must remain a miss'
    assert np.mean([c['class_metrics']['f1'] for c in classes]) == pytest.approx(card['metrics']['macro_f1'])


def test_nullable_predictions_and_all_abstention_remain_missing():
    rows = labels().astype({'prediction': 'string'})
    assert R.aggregate(rows, scope())['counts']['correct'] == 2
    rows['prediction'] = None
    card = R.aggregate(rows, scope())
    assert card['extra']['accuracy_among_calls'] is None
    assert card['extra']['accuracy_all_hidden'] == 0
    assert card['counts']['wrong'] == 0
    rows.loc[0, 'prediction'] = ''
    with pytest.raises(ValueError, match='normalized'):
        R.aggregate(rows, scope())


def test_duplicate_or_inconsistent_outcomes_are_refused():
    rows = labels()
    with pytest.raises(ValueError, match='exactly once'):
        R.aggregate(pd.concat([rows, rows]), scope())
    for column in ('correct', 'abstained'):
        broken = rows.assign(**{column: True})
        with pytest.raises(ValueError, match='contradict'):
            R.aggregate(broken, scope())
    rows.loc[0, 'truth'] = None
    with pytest.raises(ValueError, match='known explicit truth'):
        R.aggregate(rows, scope())


def test_all_other_tasks_use_the_standard_metric_functions():
    rows = pd.DataFrame({'entity': ['g1', 'g2', 'g3', 'g4'], 'truth': [1., 2., 3., 4.], 'prediction': [1.2, 2.1, 3.5, np.nan]})
    card = R.aggregate(rows, scope('trait_regression', SC.T_VALUES), parameters={'quantity_unit': 'unit'})
    assert card['metrics'] == R._finite(SC.values(rows.prediction, rows.truth))
    assert card['counts']['answered'] == 3 and card['counts']['correct'] is None
    binary = pd.DataFrame({'entity': rows.entity, 'positive': [True, False, True, False], 'score': [4., 3., 2., 1.], 'returned': [True, True, False, False]})
    assert R.aggregate(binary, scope('positive_unlabeled', SC.T_RANK), parameters={'verified_negatives': True})['metrics'] == R._finite(SC.ranking(binary.score, binary.positive))
    assert R.aggregate(binary, scope('geneset_hunt', SC.T_SET), parameters={'verified_negatives': True})['metrics'] == R._finite(SC.set_retrieval([0, 1], [0, 2], range(4)))
    cluster = labels().assign(cluster=[0, 0, 1, -1])
    assert R.aggregate(cluster, scope('holdout_search', SC.T_CLUSTER), parameters={'chosen_clusters': {'A': 0, 'B': 1}})['metrics'] == R._finite(SC.cluster_recovery(cluster.cluster, cluster.truth, {'A': 0, 'B': 1}))
    findings = pd.DataFrame({'entity': ['f1', 'f2'], 'replicated': [True, False]})
    assert R.aggregate(findings, scope('blind_battery', SC.T_REPL, 'finding'), parameters={'null_rates': [.1, .2]})['metrics'] == R._finite(SC.replication(1, 2, [.1, .2]))


def test_positive_only_truth_never_manufactures_precision_or_negatives():
    rows = pd.DataFrame({'entity': ['g1', 'g2', 'g3'], 'positive': [True, None, True], 'score': [3., 2., 1.], 'returned': [True, True, False]})
    card = R.aggregate(rows, scope('positive_unlabeled', SC.T_RANK), parameters={'verified_negatives': False, 'depth': 2})
    assert card['metrics']['auroc'] is None and card['metrics']['auprc'] is None
    assert card['extra']['known_positive_recall_at_depth'] == .5
    card = R.aggregate(rows, scope('geneset_hunt', SC.T_SET), parameters={'verified_negatives': False})
    assert card['metrics']['precision'] is None and card['extra']['known_positive_recall'] == .5
    rows.loc[1, 'positive'] = False
    with pytest.raises(ValueError, match='unknown, not false'):
        R.aggregate(rows, scope('positive_unlabeled', SC.T_RANK), parameters={'verified_negatives': False, 'depth': 2})


def test_groups_and_individual_outcomes_do_not_become_probability_claims():
    rows = pd.concat([labels().assign(entity=lambda d: d.entity + str(i), group='group' + str(i)) for i in range(5)], ignore_index=True)
    card = R.aggregate(rows, scope())
    assert card['extra']['uncertainty']['biological_groups'] == 5
    assert card['extra']['uncertainty']['status'] == 'descriptive_group_bootstrap'
    single = R.views(labels().iloc[:1], scope())
    assert single['target']['extra']['single_outcome']
    assert single['gene'][0]['per_gene_accuracy_probability'] is None
    assert all(m['component_accuracy'] is None for m in single['method'])


def test_organism_overview_keeps_seed_and_task_scope_separate():
    rows = labels()
    overview = R.organism_overview([(scope(), rows, {}), (replace(scope(), seed=4), rows, {})])
    assert len(overview['cards']) == 2
    assert overview['unique_entities_by_organism_and_unit'][0]['entities'] == 4
    assert overview['pooled_accuracy'] is None
    with pytest.raises(ValueError, match='not independent'):
        R.organism_overview([(scope(), rows, {}), (scope(), rows, {})])


def test_legacy_adapter_separates_modes_seeds_and_named_sets():
    rows = labels().rename(columns={'entity': 'gene_id'}).assign(organism=O.TOXOPLASMA,
        strategy='feature_knn', target='fixture', setting_key='k=15', seed=0, mode='together', set_name=None,
        correct=[True, None, False, True], abstained=[False, True, False, False])
    ledger = pd.concat([rows, rows.assign(seed=1), rows.assign(mode='set', set_name='whole_class')], ignore_index=True)
    cohorts = list(R.legacy_cohorts(ledger))
    assert len(cohorts) == 3
    assert all(s.truth_grade == 'unresolved' and s.gaps for s, _ in cohorts)
    with pytest.raises(ValueError, match='Repeated legacy'):
        list(R.legacy_cohorts(pd.concat([ledger, ledger])))


def test_cluster_ids_cannot_be_silently_truncated_or_filled():
    rows = labels().assign(cluster=[0., 0., 1.5, -1.])
    with pytest.raises(ValueError, match='no silent truncation'):
        R.aggregate(rows, scope('holdout_search', SC.T_CLUSTER), parameters={'chosen_clusters': {'A': 0, 'B': 1}})


def test_class_score_order_is_not_silently_reindexed():
    scores = pd.DataFrame({'A': [.8, .7, .4, .1], 'B': [.2, .3, .6, .9]}, index=[3, 2, 1, 0])
    with pytest.raises(ValueError, match='exact cohort'):
        R.aggregate(labels(), scope(), parameters={'class_scores': scores})


def test_record_and_evaluation_parameter_identities_change_with_metric_inputs():
    rows = labels()
    first = R.aggregate(rows, scope())
    changed = rows.copy()
    changed.loc[0, 'prediction'] = 'B'
    assert R.aggregate(changed, scope())['records_identity'] != first['records_identity']
    scores = pd.DataFrame({'A': [.8, .7, .4, .1], 'B': [.2, .3, .6, .9]})
    scored = R.aggregate(rows, scope(), parameters={'class_scores': scores})
    assert scored['parameter_identity'] != first['parameter_identity']
    assert scored['evaluation_parameters']['class_scores']['columns'] == ['A', 'B']
