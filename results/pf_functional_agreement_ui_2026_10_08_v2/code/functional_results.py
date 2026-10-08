"""Pinned offline functional recovery results, with explicit biological truth gaps.

The packaged view retains the original held-out artifact identity and hashes of
every displayed payload. It is a view of a frozen annotation-recovery test, not
a deployment model, a calibrated confidence estimate or new gene activity truth.
Corrupt, incomplete and differently scoped bundles are refused before display.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re

# Updated only when a verified immutable pilot is deliberately packaged.
FUNCTIONAL_RESULTS_SHA256 = '386c8e05f5b8cf32f4b9ba6c777d32933c73c08ce40d9e01d0b27959e9dbadda'
_PAYLOAD_NAMES = ('rows.json','card.json','profile_class_cards.json',
                  'major_class_cards.json','baseline_cards.json')
_CLASSES = {'1':'Oxidoreductases','2':'Transferases','3':'Hydrolases',
            '4':'Lyases','5':'Isomerases','6':'Ligases','7':'Translocases'}


def _canonical(value):
    return json.dumps(value,sort_keys=True,ensure_ascii=False,allow_nan=False,separators=(',',':'))


def _hash(value):
    return hashlib.sha256(_canonical(value).encode()).hexdigest()


def _invalid_constant(value):
    raise ValueError('Nonfinite functional result: '+value)


def profile_classes(value, namespace='ec_major'):
    """Read a complete supported profile without treating missing truth as absence."""
    if namespace not in {'ec_major','pfam'}:raise ValueError('Unsupported functional profile namespace')
    if value is None:return None
    members=json.loads(value)
    if (not isinstance(members,list) or not members or not all(isinstance(member,str) for member in members)
            or members!=sorted(set(members))
            or any(member not in _CLASSES if namespace=='ec_major' else re.fullmatch(r'PF\d{5}',member) is None for member in members)):
        raise ValueError('Invalid complete functional profile in namespace '+namespace)
    return tuple(members)


def class_title(value, namespace='ec_major'):
    """Display a pinned major enzyme-class name or the original profile identifier."""
    if namespace not in {'ec_major','pfam'}:raise ValueError('Unsupported functional profile namespace')
    if namespace=='pfam' and re.fullmatch(r'PF\d{5}',str(value)):
        try:
            from . import functional_domains as D
            metadata=D.shipped().lookup(value)
            if metadata.get('status')=='current_metadata' and metadata.get('name'):
                return value+' — '+metadata['name']+' (current nomenclature)'
        except (OSError,ValueError,KeyError,TypeError):pass
    return f'{value} — {_CLASSES[value]}' if namespace=='ec_major' and value in _CLASSES else str(value)


def member_cards(rows, namespace='ec_major'):
    """Derive recorded-presence cards; profile complements are never biological negatives."""
    actual=[profile_classes(row['truth'],namespace) for row in rows]
    predicted=[profile_classes(row['prediction'],namespace) or () for row in rows]
    if not rows or any(value is None for value in actual):raise ValueError('Member cards require complete observed reference profiles')
    terms=list(_CLASSES) if namespace=='ec_major' else sorted({term for members in actual+predicted for term in members})
    answered=sum(row['prediction'] is not None for row in rows)
    cards=[]
    for term in terms:
        positive=[term in members for members in actual];called=[term in members for members in predicted]
        tp=sum(t and p for t,p in zip(positive,called));fp=sum(not t and p for t,p in zip(positive,called));fn=sum(t and not p for t,p in zip(positive,called))
        cards.append({'class':term,'eligible':len(rows),'known_positive_genes':sum(positive),
            'true_positive':tp,'false_positive':fp,'false_negative':fn,
            'precision':tp/(tp+fp) if tp+fp else None,'recall':tp/sum(positive) if sum(positive) else None,
            'f1':2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None,'coverage':answered/len(rows),
            'reference_prevalence':sum(positive)/len(rows),'biological_precision':None,'biological_recall':None,
            'calibrated_confidence':None,'truth_semantics':'Recorded presence in complete source profiles',
            'complement_semantics':'Not recorded in this complete profile; not verified biological absence'})
    return cards


def require_context(benchmark, context):
    """Bind a frozen result to exact installed node bytes and current table values.

    This adapter currently supports the numeric-only installed-source pilot.
    A matching gene count cannot bind an imported or altered source to archived
    accuracy. Graph/representation adapters need their own dependency validation.
    """
    import pandas as pd
    from . import organisms as O
    if not isinstance(benchmark,FunctionalBenchmark) or context.organism!=benchmark.organism:
        raise ValueError('Frozen result organism differs from this context')
    references=[entry for entry in benchmark.metadata['source_manifest']['spec']['dependencies']
        if entry['kind']=='table' and entry['name']=='installed_nodes']
    if len(references)!=1:
        raise ValueError('Frozen result has no unique installed-node identity')
    path=Path(O.nodes_path(benchmark.organism))
    if hashlib.sha256(path.read_bytes()).hexdigest()!=references[0]['sha256']:
        raise ValueError('Installed node source differs from the frozen result')
    try:
        pd.testing.assert_frame_equal(context.nodes,pd.read_parquet(path),check_exact=True)
    except AssertionError as exc:
        raise ValueError('Current node values or entities differ from the frozen installed source') from exc
    return True


@dataclass(frozen=True)
class FunctionalBenchmark:
    """One verified held-out recovery artifact and its compact display payloads."""

    metadata: dict
    payloads: dict

    @property
    def organism(self):return self.metadata['organism']

    @property
    def namespace(self):return self.metadata.get('profile_namespace','ec_major')

    @property
    def rows(self):return self.payloads['rows.json']

    @property
    def card(self):return self.payloads['card.json']

    @property
    def major_class_cards(self):return self.payloads['major_class_cards.json']

    @property
    def profile_class_cards(self):return self.payloads['class_cards.json' if self.namespace=='pfam' else 'profile_class_cards.json']

    @property
    def baseline_cards(self):
        cards=self.payloads['baseline_cards.json']
        nested=self.namespace=='pfam' or self.metadata.get('control_format')=='native_ec_controls_v1'
        return {name:value['card'] for name,value in cards.items()} if nested else cards

    def source_matches(self, target):
        """Match an installed source label rather than pretending it is the derived target."""
        return target in self.metadata['source_targets'] or target==self.metadata['target']


def _benchmark(entry):
    if not isinstance(entry,dict) or not isinstance(entry.get('metadata'),dict) or not isinstance(entry.get('payloads'),dict):
        raise ValueError('Functional benchmark must retain explicit metadata and payloads')
    metadata=entry['metadata'];payloads=entry['payloads'];manifest=metadata['source_manifest']
    namespace=metadata.get('profile_namespace','ec_major')
    if namespace not in {'ec_major','pfam'}:raise ValueError('Unsupported functional profile namespace')
    control_format=metadata.get('control_format')
    native_ec=control_format=='native_ec_controls_v1'
    if control_format is not None and not native_ec:
        raise ValueError('Unsupported functional control format')
    if native_ec and (namespace!='ec_major' or metadata.get('organism')!='Pf'
            or metadata.get('target')!='ec_direct_complete_major_profile'
            or metadata.get('source_targets')!=['ec_number']):
        raise ValueError('Functional native EC control source/target address mismatch')
    if namespace=='pfam' and (metadata.get('source_targets')!=['pfam_id'] or metadata.get('target')!='pfam_id_complete_profile'):
        raise ValueError('Functional Pfam source/derived-target address mismatch')
    if metadata.get('biological_admission') is not False or metadata.get('calibrated_confidence') is not None:
        raise ValueError('Recovery tests cannot advertise biological admission or calibrated confidence')
    body={key:value for key,value in manifest.items() if key!='identity'}
    if _hash(body)!=manifest['identity'] or metadata['source_artifact_identity']!=manifest['identity']:
        raise ValueError('Functional source artifact identity changed')
    spec=manifest['spec'];scope=payloads['card.json']['scope']
    if _hash(spec)!=manifest['key']:
        raise ValueError('Functional source scope identity changed')
    native_names=_PAYLOAD_NAMES if namespace=='ec_major' else ('rows.json','card.json','class_cards.json','baseline_cards.json')
    for name in native_names:
        encoded=(_canonical(payloads[name])+'\n').encode()
        receipt=manifest['files'][name]
        if receipt['sha256']!=hashlib.sha256(encoded).hexdigest() or receipt['bytes']!=len(encoded):
            raise ValueError('Functional source payload changed: '+name)
    if namespace=='pfam':
        expected_derivation={'recipe':'recorded_profile_member_cards_v1','namespace':'pfam',
            'source_rows_sha256':manifest['files']['rows.json']['sha256'],
            'payload_sha256':_hash(payloads['major_class_cards.json'])}
        if metadata.get('membership_derivation')!=expected_derivation:
            raise ValueError('Functional membership derivation receipt changed')
    if (spec['role']!='held_out' or spec['evaluation_partition']!='test'
            or spec['confidence_kind']!='method_support'):
        raise ValueError('Functional recovery requires held-out method-support results')
    for key in ('strategy','target','benchmark_id'):
        if metadata[key]!=scope[key] or metadata[key]!=spec[key]:
            raise ValueError('Functional benchmark address mismatch: '+key)
    if metadata['organism']!=scope['organism'] or metadata['truth_grade']!=scope['truth_grade']:
        raise ValueError('Functional organism or truth grade mismatch')
    query=json.loads(spec['query_json'])
    if query['organism']!=metadata['organism']:
        raise ValueError('Functional query organism mismatch')
    if not metadata['source_targets'] or not metadata['source_context']:
        raise ValueError('Functional source labels/context must remain explicit')
    rows=payloads['rows.json'];entities=[row['entity'] for row in rows]
    if not rows or len(set(entities))!=len(entities) or entities!=spec['entity_order']:
        raise ValueError('Functional held-out gene population/order changed')
    if set(entities)&set(spec['fit_entities']):
        raise ValueError('Functional held-out genes overlap training')
    correct=answered=0
    for row in rows:
        if row['truth'] is None:raise ValueError('Unknown functional truth cannot enter recovery metrics')
        profile_classes(row['truth'],namespace);profile_classes(row['prediction'],namespace)
        if row.get('calibrated_confidence') is not None or row['abstained']!=(row['prediction'] is None):
            raise ValueError('Functional confidence or abstention semantics changed')
        answered+=row['prediction'] is not None
        correct+=row['prediction'] is not None and row['truth']==row['prediction']
    counts=payloads['card.json']['counts']
    expected={'eligible':len(rows),'answered':answered,'abstained':len(rows)-answered,
              'correct':correct,'wrong':answered-correct}
    if any(counts.get(key)!=value for key,value in expected.items()):
        raise ValueError('Functional row/card denominators disagree')
    summary=metadata['summary']
    if (summary['test']!=len(rows) or summary['train']!=len(spec['fit_entities'])
            or summary['artifact_identity']!=manifest['identity'] or summary['counts']!=counts
            or summary['metrics']!=payloads['card.json']['metrics']):
        raise ValueError('Functional source summary differs from the pinned artifact')
    if native_ec:
        original=metadata['source_summary'];evaluation=json.loads(spec['evaluation_scope_json'])
        if (any(summary[key]!=original[key] for key in ('train','test','artifact_identity','counts','metrics'))
                or any(type(row.get('training_supported')) is not bool for row in rows)
                or sum(row['training_supported'] is False for row in rows)!=original['unsupported_test_genes']
                or any(original[key]!=evaluation[key] for key in ('target_identity','cohort_identity'))
                or original['split_identity']!=scope['protocol']):
            raise ValueError('Functional native EC summary/capacity changed')
    if namespace=='pfam':
        original=metadata['source_summary']
        mapping={'train':'train_genes','test':'test_genes','eligible_profiles':'eligible_genes'}
        if (metadata.get('summary_mapping')!=mapping or any(summary[key]!=original[value] for key,value in mapping.items())
                or any(summary[key]!=original[key] for key in ('artifact_identity','counts','metrics'))
                or sum(row.get('training_supported') is False for row in rows)!=original['unsupported_test_genes']
                or any(type(row.get('training_supported')) is not bool for row in rows)):
            raise ValueError('Functional Pfam summary/capacity mapping changed')
        split=metadata['split_manifest'];receipt=metadata['source_split_receipt']
        split_identity=_hash(split)
        assignments=split['assignments'];roles=('train','tune','calibration','test')
        role_entities={role:[row['entity'] for row in assignments if row['role']==role] for role in roles}
        groups={}
        for row in assignments:
            if row['role'] not in roles or groups.setdefault(row['group'],row['role'])!=row['role']:
                raise ValueError('Functional split role/group boundary changed')
        if (receipt.get('identity')!=split_identity or scope['protocol']!=split_identity
                or split['organism']!=metadata['organism'] or split['benchmark_id']!=metadata['benchmark_id']
                or len({row['entity'] for row in assignments})!=len(assignments)
                or role_entities['train']!=spec['fit_entities'] or role_entities['test']!=entities
                or metadata.get('summary_role_mapping')!={role:role for role in roles}
                or any(not role_entities[role] or summary[role]!=len(role_entities[role]) for role in roles)
                or summary['eligible_profiles']!=len(assignments)
                or re.fullmatch(r'[a-f0-9]{64}',receipt.get('sha256','')) is None
                or type(receipt.get('bytes')) is not int or receipt['bytes']<=0
                or not any(dep['kind']=='split' and dep['sha256']==split_identity for dep in spec['dependencies'])):
            raise ValueError('Functional split-derived role mapping changed')
    classes=payloads['major_class_cards.json']
    if namespace=='ec_major' and [card['class'] for card in classes]!=list(_CLASSES):
        raise ValueError('Functional view must retain all seven major-class cards')
    if namespace=='pfam' and classes!=member_cards(rows,namespace):
        raise ValueError('Functional domain membership cards differ from recorded profiles')
    if any(card.get('biological_precision') is not None or card.get('biological_recall') is not None
            or card.get('calibrated_confidence') is not None or card['eligible']!=len(rows) for card in classes):
        raise ValueError('Functional class scope or biological interpretation changed')
    for card in classes if namespace=='ec_major' else ():
        term=card['class']
        actual=[term in profile_classes(row['truth'],namespace) for row in rows]
        predicted=[term in (profile_classes(row['prediction'],namespace) or ()) for row in rows]
        tp=sum(t and p for t,p in zip(actual,predicted));fp=sum(not t and p for t,p in zip(actual,predicted))
        fn=sum(t and not p for t,p in zip(actual,predicted))
        expected={'known_positive_genes':sum(actual),'true_positive':tp,'false_positive':fp,'false_negative':fn,
            'precision':tp/(tp+fp) if tp+fp else None,'recall':tp/sum(actual) if sum(actual) else None,
            'f1':2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None,'coverage':answered/len(rows),
            'reference_prevalence':sum(actual)/len(rows)}
        if any(card.get(key)!=value for key,value in expected.items()):
            raise ValueError('Functional major-class metrics differ from held-out rows')
    profiles=payloads['class_cards.json' if namespace=='pfam' else 'profile_class_cards.json']
    if [card['class'] for card in profiles]!=sorted(set(row['truth'] for row in rows)):
        raise ValueError('Functional view must retain every observed reference-profile class')
    for card in profiles:
        term=card['class'];actual=[row['truth']==term for row in rows]
        predicted=[row['prediction']==term for row in rows]
        tp=sum(t and p for t,p in zip(actual,predicted));fp=sum(not t and p for t,p in zip(actual,predicted))
        fn=sum(t and not p for t,p in zip(actual,predicted))
        precision=tp/(tp+fp) if tp+fp else 0.;recall=tp/(tp+fn)
        expected={'precision':precision,'recall':recall,'f1':2*precision*recall/(precision+recall) if precision+recall else 0.,
            'true_positive':tp,'false_positive':fp,'false_negative':fn,'evaluation_prevalence':sum(actual)/len(rows)}
        if (card['scope']!=scope or card['counts']['eligible']!=sum(actual)
                or any(card['class_metrics'].get(key)!=value for key,value in expected.items())):
            raise ValueError('Functional profile-class metrics differ from held-out rows')
    nested=namespace=='pfam' or native_ec
    controls={'training_majority','training_prevalence'} if nested else {'majority','prevalence_call'}
    if not controls<=set(payloads['baseline_cards.json']):
        raise ValueError('Functional view must retain matched training-only controls')
    for name,control in payloads['baseline_cards.json'].items():
        card=control['card'] if nested else control
        if (not nested and card['scope']!=scope) or card['counts']['eligible']!=len(rows):
            raise ValueError('Functional baseline must use the identical frozen test scope')
        if native_ec:
            extra=card.get('extra',{});control_scope=card['scope']
            comparison=metadata['source_summary']['baseline_comparisons'].get(name,{})
            if (name not in controls or any(control_scope.get(key)!=scope[key] for key in
                    ('organism','target','seed','task','truth_grade','unit','negative_semantics','protocol','benchmark_id','partition'))
                    or control_scope.get('strategy')!=name
                    or extra.get('control_name')!=name
                    or extra.get('baseline_parameters_estimated_from_training_only') is not True
                    or extra.get('classifier_fitting_performed') is not False
                    or extra.get('test_support_used_for_selection') is not False
                    or extra.get('biological_admission') is not False
                    or extra.get('biological_accuracy') is not None
                    or comparison.get('source_scope')!=control_scope
                    or comparison.get('counts')!=card['counts'] or comparison.get('metrics')!=card['metrics']
                    or comparison.get('same_entity_order_truth_groups_and_support') is not True):
                raise ValueError('Functional native EC baseline provenance/summary changed')
        if namespace=='pfam':
            extra=card.get('extra',{})
            control_scope=card['scope'];evaluation=json.loads(spec['evaluation_scope_json'])
            if (any(control_scope.get(key)!=scope[key] for key in ('organism','target','seed','task','truth_grade','unit','negative_semantics'))
                    or control_scope.get('strategy')!=name or control_scope.get('partition')!='test'
                    or extra.get('split_identity')!=scope['protocol']
                    or any(extra.get(key)!=evaluation[key] for key in ('target_identity','cohort_identity'))
                    or extra.get('baseline_parameters_estimated_from_training_only') is not True
                    or extra.get('classifier_or_feature_fitting_performed') is not False
                    or extra.get('test_support_used_for_selection') is not False
                    or extra.get('biological_accuracy') is not None
                    or metadata['source_summary']['baseline_comparisons'][name]['metrics']!=card['metrics']):
                raise ValueError('Functional baseline provenance/summary changed')
    return FunctionalBenchmark(metadata,payloads)


def load(path, *, expected_sha256):
    """Load a compact data-only view with an externally pinned complete-file identity."""
    path=Path(path)
    if path.is_symlink():raise ValueError('Functional bundle must be a regular local file')
    data=path.read_bytes()
    if not expected_sha256 or hashlib.sha256(data).hexdigest()!=expected_sha256:
        raise ValueError('Functional bundle checksum mismatch')
    document=json.loads(data,parse_constant=_invalid_constant)
    if not isinstance(document,dict) or document.get('schema_version')!=1 or not isinstance(document.get('benchmarks'),list):
        raise ValueError('Unsupported functional bundle schema')
    benchmarks=[_benchmark(entry) for entry in document['benchmarks']]
    addresses=[(b.organism,b.metadata['benchmark_id'],b.metadata['strategy']) for b in benchmarks]
    if len(set(addresses))!=len(addresses):raise ValueError('Duplicate functional benchmark address')
    return benchmarks


def shipped(organism):
    """Return offline results and an explicit unavailable reason; never fit or download."""
    path=Path(__file__).with_name('data')/'functional_results.json'
    if not path.exists():return [],'No frozen functional recovery results are packaged for this organism.'
    try:benchmarks=load(path,expected_sha256=FUNCTIONAL_RESULTS_SHA256)
    except (OSError,ValueError,KeyError,TypeError) as exc:
        return [],'Functional results unavailable: '+str(exc)
    selected=[benchmark for benchmark in benchmarks if benchmark.organism==organism]
    return selected,('' if selected else 'No frozen functional recovery results are packaged for this organism.')
