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

# Updated only when a verified immutable pilot is deliberately packaged.
FUNCTIONAL_RESULTS_SHA256 = '226da90d0df1e2de1933c672c8158579e1e46d656113a072eb03bbd3b14de849'
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


def profile_classes(value):
    """Read a complete EC major-class profile without treating missing truth as absence."""
    if value is None:return None
    members=json.loads(value)
    if (not isinstance(members,list) or not members or members!=sorted(set(members))
            or any(member not in _CLASSES for member in members)):
        raise ValueError('Invalid complete EC major-class profile')
    return tuple(members)


def class_title(value):
    """Display a pinned major enzyme-class name or the original profile identifier."""
    return f'{value} — {_CLASSES[value]}' if value in _CLASSES else str(value)


@dataclass(frozen=True)
class FunctionalBenchmark:
    """One verified held-out recovery artifact and its compact display payloads."""

    metadata: dict
    payloads: dict

    @property
    def organism(self):return self.metadata['organism']

    @property
    def rows(self):return self.payloads['rows.json']

    @property
    def card(self):return self.payloads['card.json']

    @property
    def major_class_cards(self):return self.payloads['major_class_cards.json']

    @property
    def profile_class_cards(self):return self.payloads['profile_class_cards.json']

    @property
    def baseline_cards(self):return self.payloads['baseline_cards.json']

    def source_matches(self, target):
        """Match an installed source label rather than pretending it is the derived target."""
        return target in self.metadata['source_targets'] or target==self.metadata['target']


def _benchmark(entry):
    if not isinstance(entry,dict) or not isinstance(entry.get('metadata'),dict) or not isinstance(entry.get('payloads'),dict):
        raise ValueError('Functional benchmark must retain explicit metadata and payloads')
    metadata=entry['metadata'];payloads=entry['payloads'];manifest=metadata['source_manifest']
    if metadata.get('biological_admission') is not False or metadata.get('calibrated_confidence') is not None:
        raise ValueError('Recovery tests cannot advertise biological admission or calibrated confidence')
    body={key:value for key,value in manifest.items() if key!='identity'}
    if _hash(body)!=manifest['identity'] or metadata['source_artifact_identity']!=manifest['identity']:
        raise ValueError('Functional source artifact identity changed')
    spec=manifest['spec'];scope=payloads['card.json']['scope']
    if _hash(spec)!=manifest['key']:
        raise ValueError('Functional source scope identity changed')
    for name in _PAYLOAD_NAMES:
        encoded=(_canonical(payloads[name])+'\n').encode()
        receipt=manifest['files'][name]
        if receipt['sha256']!=hashlib.sha256(encoded).hexdigest() or receipt['bytes']!=len(encoded):
            raise ValueError('Functional source payload changed: '+name)
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
        profile_classes(row['truth']);profile_classes(row['prediction'])
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
    classes=payloads['major_class_cards.json']
    if [card['class'] for card in classes]!=list(_CLASSES):
        raise ValueError('Functional view must retain all seven major-class cards')
    if any(card.get('biological_precision') is not None or card.get('biological_recall') is not None
            or card.get('calibrated_confidence') is not None or card['eligible']!=len(rows) for card in classes):
        raise ValueError('Functional class scope or biological interpretation changed')
    for card in classes:
        term=card['class']
        actual=[term in profile_classes(row['truth']) for row in rows]
        predicted=[term in (profile_classes(row['prediction']) or ()) for row in rows]
        tp=sum(t and p for t,p in zip(actual,predicted));fp=sum(not t and p for t,p in zip(actual,predicted))
        fn=sum(t and not p for t,p in zip(actual,predicted))
        expected={'known_positive_genes':sum(actual),'true_positive':tp,'false_positive':fp,'false_negative':fn,
            'precision':tp/(tp+fp) if tp+fp else None,'recall':tp/sum(actual) if sum(actual) else None,
            'f1':2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None,'coverage':answered/len(rows),
            'reference_prevalence':sum(actual)/len(rows)}
        if any(card.get(key)!=value for key,value in expected.items()):
            raise ValueError('Functional major-class metrics differ from held-out rows')
    profiles=payloads['profile_class_cards.json']
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
    if not {'majority','prevalence_call'}<=set(payloads['baseline_cards.json']):
        raise ValueError('Functional view must retain matched training-only controls')
    for card in payloads['baseline_cards.json'].values():
        if card['scope']!=scope or card['counts']['eligible']!=len(rows):
            raise ValueError('Functional baseline must use the identical frozen test scope')
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
