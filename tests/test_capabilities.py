"""Semantic adapter refusals and declared component tests for every strategy."""
from dataclasses import replace

import pytest

from starplast import capabilities as C, organisms as O, strategies as S, techniques as T
from starplast.query import BiologicalContext, EntityRef, Query

IDS = tuple(f'TGME49_{i:06d}' for i in range(100001, 100007))
ENTITIES = tuple(EntityRef(O.TOXOPLASMA, 'gene', i) for i in IDS)
EVIDENCE = C.EvidenceState(O.TOXOPLASMA, IDS,
    columns=(('location', 'categorical'), ('growth', 'numeric'), ('baseline', 'numeric')),
    layers=(('contacts', 'physical'), ('folds', 'structural'), ('papers', 'literature')),
    permitted_feature_count=10, orthogroups=True, attention=True,
    other_organism=O.FALCIPARUM, other_columns=('other_growth',), other_orthogroups=True,
    class_values=(('location', ('rhoptry', 'nucleus')),), admitted_gene_space=True)


def request(cap):
    if cap.minimum_seeds or cap.strategy == 'neighbour_space':
        query = Query(O.TOXOPLASMA, 'gene_set', ENTITIES)
    elif cap.outputs[0].kind == 'pair_ranking':
        query = Query(O.TOXOPLASMA, 'pair', ENTITIES[:2])
    elif cap.column_parameters:
        kind = 'trait' if cap.column_parameters[0][1] == 'numeric' else 'label'
        query = Query(O.TOXOPLASMA, kind, target='growth' if kind == 'trait' else 'location')
    else:
        query = Query(O.TOXOPLASMA, 'gene', ENTITIES[:1])
    settings = {}
    for parameter, kind in cap.column_parameters:
        settings[parameter] = 'other_growth' if kind == 'other' else 'growth' if kind == 'numeric' else 'location'
    if cap.strategy == 'condition_shift':
        settings['baseline'] = 'baseline'
    if 'layer' in cap.required:
        settings['layer'] = 'folds' if cap.strategy == 'structural_homology' else 'papers' if cap.strategy == 'attention_correction' else 'contacts'
    return query, settings


@pytest.mark.parametrize('cap', C.catalog(), ids=lambda c: c.strategy)
def test_every_strategy_resolves_a_reviewed_native_plan_with_truth_gap(cap):
    query, settings = request(cap)
    plan = C.resolve(query, cap.strategy, EVIDENCE, settings=settings)
    assert plan.applicable, (plan.reason, plan.detail)
    assert plan.reason == 'applicable' and plan.query_json == query.to_json()
    assert 'Independent biological benchmark unavailable' in plan.validation_gaps
    for output in plan.outputs:
        cap.validate_output(output)
    assert all(t.status == 'declared_not_measured' for t in cap.techniques)


@pytest.mark.parametrize('cap', C.catalog(), ids=lambda c: c.strategy)
def test_missing_gene_pack_never_silently_becomes_an_answer(cap):
    query, settings = request(cap)
    assert C.resolve(query, cap.strategy, replace(EVIDENCE, admitted_gene_space=False), settings=settings).reason == 'gene_space_unavailable'


@pytest.mark.parametrize('cap', C.catalog(), ids=lambda c: c.strategy)
def test_composite_cannot_manufacture_calibrated_gene_probability(cap):
    with pytest.raises(ValueError, match='manufacture'):
        cap.validate_output(C.Output('numeric_estimates', 'gene', 'calibrated probability of a biological label'))


@pytest.mark.parametrize('kind', ['gene', 'protein', 'gene_set', 'label', 'class', 'trait', 'pair'])
@pytest.mark.parametrize('cap', C.catalog(), ids=lambda c: c.strategy)
def test_every_query_kind_has_an_answer_or_named_refusal(kind, cap):
    query = Query(O.TOXOPLASMA, kind,
        entities=ENTITIES if kind == 'gene_set' else ENTITIES[:2] if kind == 'pair' else (EntityRef(O.TOXOPLASMA, 'protein', 'P00001'),) if kind == 'protein' else ENTITIES[:1] if kind == 'gene' else (),
        target='growth' if kind == 'trait' else 'location' if kind in {'label', 'class'} else '',
        values=('rhoptry',) if kind == 'class' else ())
    _, settings = request(cap)
    plan = C.resolve(query, cap.strategy, EVIDENCE, settings=settings)
    assert plan.reason and plan.detail
    assert plan.applicable == bool(plan.outputs)
    assert plan.applicable or not plan.benchmark_ids


def test_catalogue_and_technique_extensions_require_reviewed_declarations(monkeypatch):
    added = replace(S.catalog()[0], key='unreviewed_strategy')
    monkeypatch.setattr(S, 'catalog', lambda: [added])
    with pytest.raises(ValueError, match='explicit capability'):
        C.catalog()
    with pytest.raises(ValueError, match='explicit reviewed'):
        C.technique_validation('unreviewed_technique')


def test_all_technique_roles_and_cache_dependencies_are_explicit():
    caps = C.catalog()
    assert {t.technique for c in caps for t in c.techniques} == set(T.TECHNIQUES)
    for c in caps:
        assert c.parameters == tuple(p.name for p in S.get(c.strategy).params)
        assert {'query', 'settings', 'source_exclusions', 'fit_split', 'fit_role', 'context'} <= set(c.cache_dependencies)
        assert all(t.roles and t.test for t in c.techniques)
    assert C.get('block_ablation').outputs[0].unit == 'evidence_family'
    assert C.get('recoverability_atlas').outputs[0].unit == 'class'
    assert C.get('masked_imputation').outputs[0].semantics.startswith('imputed percentile')
    for key in ('link_prediction', 'attention_correction', 'unwritten_links', 'paralog_divergence', 'network_training'):
        assert C.get(key).outputs[0].unit == 'pair'


def test_independent_organism_context_and_universe_refusals():
    q = Query(O.TOXOPLASMA, 'gene', ENTITIES[:1], target='location')
    other = C.EvidenceState(O.FALCIPARUM, ('PF3D7_0100100',), admitted_gene_space=True)
    assert C.resolve(q, 'feature_knn', other).reason == 'organism_mismatch'
    assert C.resolve(replace(q, context=BiologicalContext(stage='bradyzoite')), 'feature_knn', EVIDENCE).reason == 'context_unavailable'
    assert C.resolve(q, 'feature_knn', replace(EVIDENCE, entity_ids=IDS[1:])).reason == 'entity_not_in_universe'
    host = Query(O.HUMAN, 'gene', (EntityRef(O.HUMAN, 'gene', 'ENSG000001'),))
    state = C.EvidenceState(O.HUMAN, ('ENSG000001',), admitted_gene_space=True)
    assert C.resolve(host, 'feature_knn', state).reason == 'organism_adapter_unavailable'
    cross = Query(O.TOXOPLASMA, 'pair', (ENTITIES[0], EntityRef(O.HUMAN, 'protein', 'P00001')))
    assert C.resolve(cross, 'link_prediction', EVIDENCE).reason == 'endpoint_unit_or_organism_unsupported'


def test_target_class_and_parameter_refusals():
    q = Query(O.TOXOPLASMA, 'class', target='location', values=('rhoptry',))
    assert C.resolve(q, 'feature_knn', replace(EVIDENCE, class_values=())).reason == 'class_values_unavailable'
    assert C.resolve(replace(q, values=('absent',)), 'feature_knn', EVIDENCE).reason == 'class_value_unavailable'
    assert C.resolve(q, 'feature_knn', EVIDENCE, settings={'target': 'growth'}).reason == 'target_conflict'
    assert C.resolve(replace(q, target='growth'), 'feature_knn', EVIDENCE).reason == 'target_kind_mismatch'
    for settings in ({'alien': 1}, {'k': True}, {'k': float('nan')}, {'k': 1}, {'k': '15'}, {'target': ''}):
        plan = C.resolve(Query(O.TOXOPLASMA, 'gene', ENTITIES[:1]), 'feature_knn', EVIDENCE, settings=settings)
        assert plan.reason == 'invalid_settings'
    assert C.resolve(Query(O.TOXOPLASMA, 'gene', ENTITIES[:1]), 'feature_knn', EVIDENCE).reason == 'column_required'


def test_seed_and_required_evidence_refusals():
    small = Query(O.TOXOPLASMA, 'gene_set', ENTITIES[:2])
    assert C.resolve(small, 'positive_unlabeled', EVIDENCE).reason == 'insufficient_seeds'
    q = Query(O.TOXOPLASMA, 'gene_set', ENTITIES)
    assert C.resolve(q, 'positive_unlabeled', EVIDENCE, settings={'genes': IDS[0]}).reason == 'seed_conflict'
    label = Query(O.TOXOPLASMA, 'label', target='location')
    assert C.resolve(label, 'feature_knn', replace(EVIDENCE, permitted_feature_count=0)).reason == 'required_evidence_unavailable'
    assert C.resolve(label, 'physical_partners', replace(EVIDENCE, layers=())).reason == 'required_evidence_unavailable'
    assert C.resolve(label, 'structural_homology', EVIDENCE, settings={'layer': 'contacts'}).reason == 'layer_kind_mismatch'
    assert C.resolve(Query(O.TOXOPLASMA, 'pair', ENTITIES[:2]), 'attention_correction', EVIDENCE, settings={'layer': 'contacts', 'target': 'location'}).reason == 'layer_kind_mismatch'


def test_adaptive_ortholog_output_preserves_task():
    numeric = C.resolve(Query(O.TOXOPLASMA, 'trait', target='growth'), 'ortholog_transfer', EVIDENCE, settings={'source': 'other_growth'})
    labels = C.resolve(Query(O.TOXOPLASMA, 'label', target='location'), 'ortholog_transfer', EVIDENCE, settings={'source': 'other_growth'})
    assert numeric.applicable and labels.applicable
    assert [o.kind for o in numeric.outputs] == ['numeric_estimates']
    assert [o.kind for o in labels.outputs] == ['label_calls']
    assert C.resolve(replace(Query(O.TOXOPLASMA, 'label', target='location'), output='meta_inference'), 'feature_knn', EVIDENCE).reason == 'output_adapter_unavailable'


def test_benchmark_links_match_effective_target_and_biological_unit():
    from types import SimpleNamespace
    def benchmark(identity, target='location', unit='gene'):
        return SimpleNamespace(organism=O.TOXOPLASMA, task='label calls', target=target,
            evaluation_unit=unit, status='candidate', benchmark_id=identity)
    q = Query(O.TOXOPLASMA, 'gene', ENTITIES[:1])
    plan = C.resolve(q, 'feature_knn', EVIDENCE, settings={'target': 'location'},
        benchmarks=(benchmark('right'), benchmark('other_target', 'growth'), benchmark('wrong_unit', unit='protein')))
    assert plan.benchmark_ids == ('right',)
    assert C.resolve(Query(O.TOXOPLASMA, 'pair', ENTITIES[:2]), 'attention_correction', EVIDENCE,
        settings={'layer': 'papers'}).reason == 'column_required'
    assert C.resolve(Query(O.TOXOPLASMA, 'trait', target='growth'), 'condition_shift', EVIDENCE,
        settings={'baseline': 'growth'}).reason == 'distinct_columns_required'


def test_network_only_mode_cannot_use_features_as_network_evidence():
    query = Query(O.TOXOPLASMA, 'gene_set', ENTITIES)
    plan = C.resolve(query, 'seed_expansion', replace(EVIDENCE, layers=()), settings={'mode': 'networks'})
    assert plan.reason == 'required_evidence_unavailable'
