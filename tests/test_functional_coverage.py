"""Functional coverage keeps organisms, target encodings and evaluation meanings separate."""
from copy import deepcopy
import json

import pytest

from starplast import functional_coverage as F, functional_results as R, organisms as O


def _catalogues():
    return {O.TOXOPLASMA:[{'organism':O.TOXOPLASMA,'target':'ec_number','family':'Function',
        'genes':4,'annotated_genes':3,'unannotated_genes':1,'classes':2,'source_ids':'recorded_ec'},
        {'organism':O.TOXOPLASMA,'target':'compartment','family':'Localization',
        'genes':4,'annotated_genes':4,'unannotated_genes':0,'classes':2}],
        O.FALCIPARUM:[{'organism':O.FALCIPARUM,'target':'ec_number','family':'Function',
        'genes':7,'annotated_genes':2,'unannotated_genes':5,'classes':1,'source_ids':'other_ec'}]}


def _capabilities():
    return [{'strategy':'feature_knn','organisms':[O.TOXOPLASMA,O.FALCIPARUM],
        'query_kinds':['gene','class','label'],'benchmark_tasks':['label calls'],'benchmark_unit':'gene','required':['features']},
        {'strategy':'trait_regression','organisms':[O.TOXOPLASMA,O.FALCIPARUM],
        'query_kinds':['gene','trait'],'benchmark_tasks':['values'],'benchmark_unit':'gene','required':['features','groups']}]


def _row(document,organism,source,target,strategy='feature_knn'):
    return next(row for row in document['rows'] if row['organism']==organism and row['source_label']==source
        and row['derived_target']==target and row['strategy']==strategy)


def test_coverage_reconciles_counts_and_never_manufactures_biology_calibration_or_deployment():
    document=F.build(_catalogues(),capabilities=_capabilities())
    assert document['summary']['coverage_addresses']==len(document['rows'])
    assert sum(document['summary']['addresses_by_organism'].values())==len(document['rows'])
    assert sum(document['summary']['applicability_states'].values())==len(document['rows'])
    assert document['summary']['source_labels']==2
    assert document['summary']['pooled_accuracy'] is None and document['summary']['pooled_annotation_genes'] is None
    assert document['summary']['independent_biological_tests']==document['summary']['calibrated_targets']==document['summary']['deployment_targets']==0
    for row in document['rows']:
        assert row['independent_biological_test']=={'status':'unavailable','admitted':False,'accuracy':None}
        assert row['calibration']['calibrated_confidence'] is None
        assert row['deployment']['unknown_gene_claims'] is None
        assert row['applicable'] is None
    assert _row(document,O.TOXOPLASMA,'ec_number','ec_number')['applicability_status']=='requires_multivalued_target_adapter'
    assert _row(document,O.TOXOPLASMA,'ec_number','ec_major_classes')['applicability_status']=='declared_adapter_inputs_not_verified'
    for organism in (O.HUMAN,O.MOUSE):
        host=_row(document,organism,'ec_number','ec_major_classes')
        assert host['known_annotation'] is None and host['source_inventory_status']=='not_inventoried'
        assert host['applicability_status']=='organism_adapter_unavailable'


def test_exact_legacy_addresses_do_not_become_recovery_or_independent_truth():
    legacy=[{'organism':O.FALCIPARUM,'target':'ec_number','strategy':'feature_knn','task':'label calls','rows':10}]
    document=F.build(_catalogues(),capabilities=_capabilities(),legacy_records=legacy)
    pf=_row(document,O.FALCIPARUM,'ec_number','ec_number')
    assert pf['legacy_result']['rows']==10 and pf['legacy_result']['accuracy'] is None
    assert pf['reference_recovery']['status']=='unavailable'
    assert _row(document,O.TOXOPLASMA,'ec_number','ec_number')['legacy_result']['status']=='unavailable'
    assert _row(document,O.FALCIPARUM,'ec_number','ec_major_classes')['legacy_result']['status']=='unavailable'


def test_verified_packaged_recovery_never_migrates_to_other_source_or_organism():
    benchmarks,_=R.shipped(O.TOXOPLASMA)
    catalogues=_matching_catalogue(benchmarks[0])
    document=F.build(catalogues,capabilities=_capabilities(),benchmarks=benchmarks)
    target=_row(document,O.TOXOPLASMA,'ec_number','ec_major_classes')
    assert target['reference_recovery']['status']=='verified_reference_recovery'
    assert target['reference_recovery']['metrics']==benchmarks[0].card['metrics']
    assert target['reference_recovery']['truth_grade']=='unresolved'
    assert target['independent_biological_test']['status']=='unavailable'
    assert _row(document,O.FALCIPARUM,'ec_number','ec_major_classes')['reference_recovery']['status']=='unavailable'
    assert _row(document,O.TOXOPLASMA,'ec_number','ec_number')['reference_recovery']['status']=='unavailable'
    assert document['summary']['reference_recovery_addresses']==document['summary']['unique_reference_artifacts']==2
    domain=_row(document,O.TOXOPLASMA,'pfam_id','pfam_id_complete_profile')
    assert domain['reference_recovery']['rows']==635
    assert domain['reference_recovery']['correct']==12
    assert domain['reference_recovery']['abstained']==615
    assert domain['reference_recovery']['source_artifact_identity']!=target['reference_recovery']['source_artifact_identity']
    assert domain['independent_biological_test']['status']=='unavailable'


@pytest.mark.parametrize('mutation', ['catalogue_organism','count_gap','duplicate_source','legacy_duplicate','unknown_organism'])
def test_ambiguous_or_nonreconciling_coverage_inputs_fail_closed(mutation):
    catalogues=deepcopy(_catalogues());legacy=[];organisms=F.DEFAULT_ORGANISMS
    if mutation=='catalogue_organism':catalogues[O.TOXOPLASMA][0]['organism']=O.FALCIPARUM
    elif mutation=='count_gap':catalogues[O.TOXOPLASMA][0]['unknown_annotation_genes']=99;catalogues[O.TOXOPLASMA][0]['unannotated_genes']=99
    elif mutation=='duplicate_source':catalogues[O.TOXOPLASMA].append(deepcopy(catalogues[O.TOXOPLASMA][0]))
    elif mutation=='legacy_duplicate':legacy=[{'organism':O.TOXOPLASMA,'target':'ec_number','strategy':'feature_knn','task':'label calls','rows':10}]*2
    elif mutation=='unknown_organism':organisms=('unknown',)
    with pytest.raises(ValueError):F.build(catalogues,capabilities=_capabilities(),legacy_records=legacy,organisms=organisms)


def test_census_refuses_loose_result_json_and_nonrecovery_admission():
    with pytest.raises(TypeError):F.build(_catalogues(),capabilities=_capabilities(),benchmarks=[{'biological_admission':True}])
    benchmark=deepcopy(R.shipped(O.TOXOPLASMA)[0][0])
    benchmark.metadata['biological_admission']=True
    with pytest.raises(ValueError):F.build(_catalogues(),capabilities=_capabilities(),benchmarks=[benchmark])


def test_duplicate_recovery_scopes_are_not_pooled():
    benchmark=R.shipped(O.TOXOPLASMA)[0][0]
    with pytest.raises(ValueError,match='Multiple recovery scopes'):
        F.build(_matching_catalogue(benchmark),capabilities=_capabilities(),benchmarks=[benchmark,benchmark])


def _matching_catalogue(benchmark):
    catalogues=_catalogues()
    evaluation=json.loads(benchmark.metadata['source_manifest']['spec']['evaluation_scope_json'])
    genes=evaluation['source_annotated_genes']+evaluation['unknown_genes']
    source=catalogues[O.TOXOPLASMA][0]
    source.update(genes=genes,annotated_genes=evaluation['source_annotated_genes'],unannotated_genes=evaluation['unknown_genes'])
    return catalogues


def test_frozen_recovery_cannot_be_attached_to_an_impossible_current_source_population():
    benchmark=R.shipped(O.TOXOPLASMA)[0][0]
    with pytest.raises(ValueError,match='exceeds known source'):
        F.build(_catalogues(),capabilities=_capabilities(),benchmarks=[benchmark])
    catalogues=_matching_catalogue(benchmark)
    catalogues[O.TOXOPLASMA][0]['genes']+=1;catalogues[O.TOXOPLASMA][0]['unannotated_genes']+=1
    with pytest.raises(ValueError,match='source universe differs'):
        F.build(catalogues,capabilities=_capabilities(),benchmarks=[benchmark])
