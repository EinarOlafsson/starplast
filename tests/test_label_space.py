"""Native label scope and frozen cohort arithmetic use ordinary synthetic fixtures."""
from dataclasses import FrozenInstanceError
import json

import pandas as pd
import pytest

from starplast import organisms as O
from starplast.label_space import LabelSpace
from starplast.query import BiologicalContext

GENES = ['TGME49_100001', 'TGME49_100002', 'TGME49_100003', 'TGME49_100004']


def nodes():
    return pd.DataFrame({'gene_id': GENES, 'compartment': ['A', 'A', 'B', None],
                         'has_domain': pd.Series([False, True, pd.NA, True], dtype='boolean'),
                         'ec_number': ['1.1.1.1 (first); 2.2.2.2 (second)', '1.1.1.1', None, None],
                         'interpro_id': ['IPR000001;IPR000002', None, None, None],
                         'interpro_desc': ['original one;original two', None, None, None],
                         'pfam_id': ['PF00001.2;PF00002', None, None, None]})


def ledger(seed=1, setting='fixed'):
    return pd.DataFrame({'organism': [O.TOXOPLASMA] * 3, 'strategy': ['feature_knn'] * 3,
                         'target': ['compartment'] * 3, 'setting_key': [setting] * 3,
                         'seed': [seed] * 3, 'mode': ['together'] * 3, 'set_name': [None] * 3,
                         'gene_id': GENES[:3], 'truth': ['A', 'A', 'B'],
                         'prediction': ['A', None, 'A'], 'correct': [True, None, False],
                         'abstained': [False, True, False]})


def knn(model, value=None):
    return next(row for row in model.mechanisms('compartment', value) if row['strategy'] == 'feature_knn')


def test_exact_membership_unknown_counts_and_native_false_identity():
    model = LabelSpace(nodes(), O.TOXOPLASMA)
    label = model.labels.set_index('target').loc['compartment']
    assert label.annotated_genes == 3 and label.unannotated_genes == 1
    assert model.members('compartment', 'A').gene_id.tolist() == GENES[:2]
    false = model.classes('has_domain').set_index('value').loc['False']
    assert false.native_value is False and false.annotated_genes == 1
    assert model.members('has_domain', False).gene_id.tolist() == [GENES[0]]
    assert model.query('has_domain', False).values == ('False',)


def test_all_native_function_memberships_and_original_descriptions_are_retained():
    model = LabelSpace(nodes(), O.TOXOPLASMA)
    assert set(model.classes('ec_number').value) == {'1.1.1.1', '2.2.2.2'}
    assert len(model.members('ec_number')) == 3
    assert set(model.classes('pfam_id').value) == {'PF00001.2', 'PF00002'}
    assert set(model.classes('interpro_id').value) == {'IPR000001', 'IPR000002'}
    one = model.classes('interpro_id').set_index('value').loc['IPR000001']
    assert one.source_description == 'original one'
    assert model.members('ec_number', '2.2.2.2').description.iloc[0] == '2.2.2.2 (second)'


def test_strict_domains_do_not_create_partial_truth_from_malformed_source_cells():
    table = nodes()
    table.loc[0, 'pfam_id'] = 'PF00001.2; malformed PF00002'
    model = LabelSpace(table, O.TOXOPLASMA)
    assert model.members('pfam_id').empty
    label = model.labels.set_index('target').loc['pfam_id']
    assert label.annotated_genes == 0 and label.unannotated_genes == 4
    assert label.unparsed_annotation_genes == 1


def test_explicit_list_labels_and_function_lists_keep_every_member():
    table = nodes()
    table['protein_classes'] = pd.Series([['one', 'two'], ['one'], [], None], dtype=object)
    table['pfam_id'] = pd.Series([['PF00001.2', 'PF00002'], [], None, None], dtype=object)
    model = LabelSpace(table, O.TOXOPLASMA)
    assert set(model.classes('protein_classes').value) == {'one', 'two'}
    assert len(model.members('protein_classes')) == 3
    assert set(model.classes('pfam_id').value) == {'PF00001.2', 'PF00002'}


def test_evidence_card_native_balance_coverage_and_typed_query_are_distinct():
    context = BiologicalContext(stage='tachyzoite')
    model = LabelSpace(nodes(), O.TOXOPLASMA, context=context)
    assert model.query('compartment', 'A').context == context
    card = model.evidence_card('compartment', 'A')
    assert card.kind == 'evidence' and card.metrics == ()
    assert dict(card.counts)['whole_universe'] == 4
    assert dict(card.counts)['class_members'] == 2 and dict(card.counts)['unknown_genes'] == 1
    assert json.loads(card.snapshot_json)['card']['evidence']['query']['values'] == ['A']
    assert json.loads(card.snapshot_json)['card']['evidence']['hierarchy']['status'] == 'unavailable'
    assert json.loads(card.snapshot_json)['card']['evidence']['unmeasured_classes']['status'] == 'unavailable'
    with pytest.raises(FrozenInstanceError):
        card.status = 'changed'


def test_all_39_strategies_visible_with_null_unavailable_metrics_without_fits():
    model = LabelSpace(nodes(), O.TOXOPLASMA)
    methods = model.mechanisms('compartment')
    assert len(methods) == 39 and len({row['strategy'] for row in methods}) == 39
    assert all(row['metrics'] is None and row['card'] is None for row in methods)
    assert all(row['gaps'] for row in methods)
    assert not next(row for row in methods if row['strategy'] == 'trait_regression')['applicable']


def test_shared_class_precision_false_calls_abstentions_and_full_confusion():
    table = nodes()
    model = LabelSpace(table, O.TOXOPLASMA, ledger=ledger(), ledger_nodes=table)
    target = knn(model)
    assert target['metrics']['accuracy'] == 1 / 3
    assert target['metrics']['macro_recall'] == .25
    result = knn(model, 'A')['evaluations'][0]
    assert result['class_metrics']['precision'] == .5 and result['class_metrics']['recall'] == .5
    assert result['class_metrics']['false_positive'] == 1 and result['class_metrics']['false_negative'] == 1
    assert {'truth': 'B', 'prediction': 'A', 'count': 1} in result['confusion']
    assert {'truth': 'A', 'prediction': None, 'count': 1} in result['confusion']
    assert result['card'].kind == 'performance' and result['scope']['truth_grade'] == 'unresolved'
    saved = json.loads(result['card'].snapshot_json)
    assert saved['card']['extra']['class_metrics']['false_positive'] == 1
    assert len(saved['card']['extra']['full_cohort_confusion']) == 3
    assert model.labels.set_index('target').loc['compartment'].held_out_rows == 3


def test_seed_settings_and_modes_keep_separate_evaluations_not_pooled_rates():
    table = nodes()
    records = pd.concat([ledger(), ledger(seed=2, setting='other')], ignore_index=True)
    model = LabelSpace(table, O.TOXOPLASMA, ledger=records, ledger_nodes=table)
    method = knn(model)
    assert len(method['evaluations']) == 2 and method['metrics'] is None and method['card'] is None
    assert {row['scope']['seed'] for row in method['evaluations']} == {1, 2}


@pytest.mark.parametrize('change', ['value', 'order', 'column', 'unbound', 'context'])
def test_changed_tables_and_requested_context_cannot_borrow_archived_accuracy(change):
    original, current = nodes(), nodes()
    binding, context = original, None
    if change == 'value':
        current.loc[0, 'compartment'] = 'B'
    elif change == 'order':
        current = current.iloc[::-1]
    elif change == 'column':
        current['new'] = 0
    elif change == 'unbound':
        binding = None
    else:
        context = BiologicalContext(condition='specific assay')
    model = LabelSpace(current, O.TOXOPLASMA, ledger=ledger(), ledger_nodes=binding, context=context)
    assert knn(model)['metrics'] is None
    assert any('withheld' in gap for gap in knn(model)['gaps'])


def test_frozen_truth_mismatch_is_explicit_unavailable_even_with_equal_binding():
    table = nodes()
    records = ledger()
    records.loc[0, 'truth'] = 'B'
    model = LabelSpace(table, O.TOXOPLASMA, ledger=records, ledger_nodes=table)
    assert knn(model)['metrics'] is None
    assert any('truth/member identity' in gap for gap in knn(model)['gaps'])


def test_classes_species_and_label_text_do_not_merge_and_inputs_are_copied():
    table = nodes()
    model = LabelSpace(table, O.TOXOPLASMA)
    table.loc[0, 'compartment'] = 'changed'
    assert model.members('compartment', 'A').gene_id.tolist() == GENES[:2]
    copied = model.members('compartment')
    copied.loc[0, 'value'] = 'changed'
    assert model.members('compartment').value.iloc[0] == 'A'
    with pytest.raises(ValueError):
        model.query('compartment', 'unknown')
    with pytest.raises(ValueError):
        LabelSpace(nodes(), O.FALCIPARUM)
