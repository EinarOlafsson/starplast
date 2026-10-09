"""Functional evidence and strategy coverage without borrowed biological accuracy.

Addresses preserve organism, original source label, derived target, strategy and
task. Declarations, annotation membership and recorded outcomes are separate.
The current evidence contains reference-recovery tests only: independent biology,
calibrated confidence and unknown-gene deployment remain explicitly unavailable.
This census never fits methods or pools accuracy across settings or tasks.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import asdict,is_dataclass
import json

from . import organisms as O

SCHEMA_VERSION=1
DEFAULT_ORGANISMS=(O.TOXOPLASMA,O.FALCIPARUM,O.HUMAN,O.MOUSE)
_HOST_SOURCE_SLOTS=('ec_number','interpro_id','pfam_id')


def _records(value):
    if hasattr(value,'to_dict'):return value.to_dict('records')
    return list(value)


def _capability(value):
    return value.to_dict() if hasattr(value,'to_dict') else asdict(value) if is_dataclass(value) else dict(value)


def _count(value,name):
    if type(value) is not int or value<0:raise ValueError(name+' requires a nonnegative integer count')
    return value


def _annotations(catalogues,organisms):
    result={}
    for organism,frame in catalogues.items():
        if organism not in organisms:raise ValueError('Catalogue organism outside requested census')
        for row in _records(frame):
            if row.get('organism')!=organism:raise ValueError('Catalogue organism address mismatch')
            if row.get('family')!='Function':continue
            key=(organism,row['target'])
            if key in result:raise ValueError('Duplicate functional source-label catalogue address')
            genes=_count(row['genes'],'genes');known=_count(row['annotated_genes'],'annotated_genes')
            unknown=_count(row['unannotated_genes'],'unannotated_genes')
            if known+unknown!=genes:raise ValueError('Source annotation gene counts do not reconcile')
            result[key]={'genes':genes,'known_annotation_genes':known,'unknown_annotation_genes':unknown,
                'known_classes':_count(row['classes'],'classes'),'source_ids':row.get('source_ids','unresolved')}
    return result


def _legacy(records,organisms):
    result={}
    for entry in _records(records):
        if entry['organism'] not in organisms:raise ValueError('Legacy result organism outside requested census')
        key=(entry['organism'],entry['target'],entry['strategy'],entry['task'])
        if key in result:raise ValueError('Duplicate scoped legacy result summary')
        result[key]={'status':'recorded_legacy_scope_unresolved','rows':_count(entry['rows'],'legacy rows'),
            'biological_admission':False,'accuracy':None,
            'interpretation':'Historical rows exist; source/split/independence admission is not inherited'}
    return result


def _recovery(benchmarks,organisms,capabilities,annotations):
    from . import functional_results as F
    result={}
    declared={cap['strategy']:cap for cap in capabilities}
    for benchmark in benchmarks:
        if not isinstance(benchmark,F.FunctionalBenchmark):
            raise TypeError('Functional recovery requires the checksum-validated packaged result contract')
        benchmark=F._benchmark({'metadata':benchmark.metadata,'payloads':benchmark.payloads})
        meta=benchmark.metadata;card=benchmark.card;scope=card['scope']
        if meta['organism'] not in organisms:raise ValueError('Benchmark organism outside requested census')
        if meta['biological_admission'] is not False or meta['calibrated_confidence'] is not None:
            raise ValueError('Current functional census accepts verified reference recovery only')
        if (meta['organism']!=scope['organism'] or meta['target']!=scope['target']
                or meta['strategy']!=scope['strategy']):raise ValueError('Recovery scope address mismatch')
        cap=declared.get(meta['strategy'])
        if not cap or scope['task'] not in cap['benchmark_tasks'] or meta['organism'] not in cap['organisms']:
            raise ValueError('Recovery result has no declared organism/strategy/task adapter')
        for source in meta['source_targets']:
            annotation=annotations.get((meta['organism'],source))
            if annotation is not None:
                if annotation['known_annotation_genes']<meta['summary']['eligible_profiles']:
                    raise ValueError('Frozen recovery population exceeds known source annotations')
                evaluation=json.loads(meta['source_manifest']['spec']['evaluation_scope_json'])
                if {'source_annotated_genes','unknown_genes'}<=set(evaluation):
                    if annotation['genes']!=evaluation['source_annotated_genes']+evaluation['unknown_genes']:
                        raise ValueError('Frozen recovery source universe differs from the current catalogue')
            key=(meta['organism'],source,meta['target'],meta['strategy'],scope['task'])
            if key in result:raise ValueError('Multiple recovery scopes require separate coverage addresses')
            result[key]={'status':'verified_reference_recovery','rows':card['counts']['eligible'],
                'answered':card['counts']['answered'],'abstained':card['counts']['abstained'],
                'correct':card['counts']['correct'],'wrong':card['counts']['wrong'],
                'metrics':dict(card['metrics']),'source_artifact_identity':meta['source_artifact_identity'],
                'benchmark_id':scope['benchmark_id'],'truth_grade':scope['truth_grade'],
                'settings':scope['settings'],'protocol':scope['protocol'],
                'partition':scope['partition'],'biological_admission':False,
                'calibrated_confidence':None,'interpretation':'Recovery of recorded reference profiles; independent biological activity accuracy unknown'}
    return result


def build(catalogues_by_organism, *, capabilities=None, benchmarks=(), legacy_records=(), organisms=DEFAULT_ORGANISMS):
    """Return a complete, scoped matrix and counts from supplied verified evidence.

    Catalogues are current label-level rows; legacy inputs are exact
    organism/target/strategy/task row-count summaries. ``benchmarks`` must come
    from ``functional_results.load``/``shipped`` integrity validation. This first
    contract deliberately accepts no loose biological/calibration/deployment
    admission flags. Availability declarations do not assert runnable inputs.
    """
    organisms=tuple(organisms)
    if len(set(organisms))!=len(organisms) or set(organisms)-set(DEFAULT_ORGANISMS):
        raise ValueError('Explicit unique supported census organisms required')
    if capabilities is None:
        from . import capabilities as C
        capabilities=C.catalog()
    capabilities=[_capability(value) for value in capabilities]
    if len({cap['strategy'] for cap in capabilities})!=len(capabilities):
        raise ValueError('Duplicate strategy capability declaration')
    annotations=_annotations(catalogues_by_organism,organisms)
    legacy=_legacy(legacy_records,organisms)
    recovery=_recovery(benchmarks,organisms,capabilities,annotations)
    targets={(organism,source,source) for organism,source in annotations}
    for organism in organisms:
        if organism not in catalogues_by_organism:
            targets.update((organism,source,source) for source in _HOST_SOURCE_SLOTS)
    # EC source labels share a declared profile question, but results and actual
    # eligibility never migrate between organisms or orthology-derived sources.
    targets.update((organism,source,'ec_major_classes') for organism,source,target in tuple(targets)
        if source.startswith('ec_number'))
    targets.update((key[0],key[1],key[2]) for key in recovery)
    rows=[]
    for organism,source,target in sorted(targets):
        annotation=annotations.get((organism,source))
        encoding=('complete_categorical_profile' if target=='ec_major_classes' else
            'overlapping_source_memberships' if source.startswith(('interpro','pfam','ec_number')) else 'recorded_annotation_flag')
        for cap in capabilities:
            organism_declared=organism in cap['organisms']
            query_declared=bool({'label','class'}&set(cap['query_kinds']))
            for task in cap['benchmark_tasks']:
                key=(organism,source,target,cap['strategy'],task)
                status=('organism_adapter_unavailable' if not organism_declared else
                    'source_not_inventoried' if annotation is None else
                    'label_query_adapter_unavailable' if not query_declared else
                    'requires_multivalued_target_adapter' if encoding=='overlapping_source_memberships' else
                    'numeric_target_adapter_incompatible' if task=='values' else
                    'declared_adapter_inputs_not_verified')
                recorded=legacy.get((organism,target,cap['strategy'],task))
                result=recovery.get(key)
                rows.append({'organism':organism,'source_label':source,'derived_target':target,
                    'strategy':cap['strategy'],'task':task,'target_encoding':encoding,
                    'benchmark_unit':cap.get('benchmark_unit','unresolved'),
                    'truth_requirement':cap.get('truth_requirement','unresolved'),
                    'required_evidence':list(cap.get('required',())),
                    'organism_adapter_declared':organism_declared,'label_query_adapter_declared':query_declared,
                    'applicability_status':status,'applicable':None,
                    'source_inventory_status':'available' if annotation is not None else 'not_inventoried',
                    'known_annotation':dict(annotation) if annotation is not None else None,
                    'legacy_result':recorded if recorded is not None else {'status':'unavailable','rows':None,'accuracy':None},
                    'reference_recovery':result if result is not None else {'status':'unavailable','rows':None,'metrics':None},
                    'independent_biological_test':{'status':'unavailable','admitted':False,'accuracy':None},
                    'calibration':{'status':'unavailable','calibrated_confidence':None},
                    'deployment':{'status':'unavailable','unknown_gene_claims':None},
                    'negative_semantics':'Missing annotation and false/zero annotation flags do not establish biological absence'})
    addresses=[(row['organism'],row['source_label'],row['derived_target'],row['strategy'],row['task']) for row in rows]
    if len(set(addresses))!=len(addresses):raise ValueError('Duplicate functional coverage address')
    recovered=[row for row in rows if row['reference_recovery']['status']=='verified_reference_recovery']
    result_identities={row['reference_recovery']['source_artifact_identity'] for row in recovered}
    summary={'coverage_addresses':len(rows),'source_labels':len(annotations),'organisms':list(organisms),
        'addresses_by_organism':dict(Counter(row['organism'] for row in rows)),
        'applicability_states':dict(Counter(row['applicability_status'] for row in rows)),
        'reference_recovery_addresses':len(recovered),'unique_reference_artifacts':len(result_identities),
        'independent_biological_tests':0,'calibrated_targets':0,'deployment_targets':0,
        'pooled_accuracy':None,'pooled_annotation_genes':None,
        'interpretation':'Counts describe coverage addresses, not independent studies, genes or pooled method performance'}
    return {'schema_version':SCHEMA_VERSION,'rows':rows,'summary':summary,
        'source_label_inventory':[{'organism':organism,'source_label':source,**data}
            for (organism,source),data in sorted(annotations.items())]}
