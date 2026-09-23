"""Versioned observations and reviewable literature assertions.

Gene-level matrices remain useful model inputs, but they are derived summaries.
This module preserves source, biological entity, condition, uncertainty and
missingness before aggregation. Parquet files use JSON scalar values so booleans,
numbers and text keep their types without a mixed-type Arrow column.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

EVIDENCE_STATUSES = {"measured", "computed", "orthology_transfer", "model_prediction", "unclassified"}
MISSING_STATES = {"observed", "not_assayed", "below_detection", "failed_qc", "ambiguous_mapping", "unknown"}
ENTITY_TYPES = {"gene", "transcript", "protein", "residue", "allele", "pair", "host_gene", "metabolite"}


@dataclass(frozen=True)
class Observation:
    """One source-specific observation; contradictory records are retained separately."""
    entity_id: str
    trait: str
    value: object
    source_id: str
    organism: str
    entity_type: str = "gene"
    entity_version: str = "unspecified"
    source_version: str = "unspecified"
    source_url: str = ""
    source_location: str = ""
    unit: str = "unspecified"
    context: dict = field(default_factory=dict)
    replicate: str = ""
    uncertainty: float | None = None
    uncertainty_kind: str = "unspecified"
    missing_state: str = "observed"
    evidence_status: str = "measured"
    derived_from: tuple = ()

    def __post_init__(self):
        if not all((self.entity_id, self.trait, self.source_id, self.organism)):
            raise ValueError("entity, trait, source and organism are required")
        if self.entity_type not in ENTITY_TYPES or self.evidence_status not in EVIDENCE_STATUSES:
            raise ValueError("invalid entity type or evidence status")
        if self.missing_state not in MISSING_STATES:
            raise ValueError("invalid missingness state")
        if self.value is None and self.missing_state == "observed":
            raise ValueError("a missing value needs an explicit missingness state")
        if self.missing_state in {"not_assayed", "failed_qc", "ambiguous_mapping", "unknown"} and self.value is not None:
            raise ValueError("unobserved values must be None; store detection limits in context")
        if self.uncertainty is not None and (not np.isfinite(self.uncertainty) or self.uncertainty < 0):
            raise ValueError("uncertainty must be finite and non-negative")
        if self.uncertainty is not None:
            object.__setattr__(self, 'uncertainty', float(self.uncertainty))
        object.__setattr__(self, 'derived_from', tuple(self.derived_from))
        json.dumps(asdict(self), allow_nan=False)

    @property
    def record_id(self):
        """Content identity distinguishes source versions and experimental conditions."""
        return hashlib.sha256(json.dumps(asdict(self),sort_keys=True).encode()).hexdigest()


def write_observations(observations, path):
    """Write a deduplicated, immutable-content snapshot and return its SHA-256."""
    rows = []
    for observation in observations:
        row = asdict(observation)
        row['record_id'] = observation.record_id
        for column in ('value','context','derived_from'):
            row[column+'_json'] = json.dumps(row.pop(column),sort_keys=True,allow_nan=False)
        rows.append(row)
    if not rows:
        raise ValueError("no observations to save")
    table = pd.DataFrame(rows).drop_duplicates('record_id').sort_values('record_id')
    path = Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary = path.with_suffix('.tmp.parquet')
    table.to_parquet(temporary,index=False)
    temporary.replace(path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    metadata = {'schema_version':1,'records':len(table),'sha256':digest,
                'created_utc':datetime.now(timezone.utc).isoformat(),
                'sources':sorted(table.source_id.unique())}
    path.with_suffix('.json').write_text(json.dumps(metadata,indent=2)+'\n')
    return digest


def read_observations(path):
    """Reconstitute typed observations and validate each stored content identity."""
    records = []
    for row in pd.read_parquet(path).to_dict('records'):
        record_id = row.pop('record_id')
        for column in ('value','context','derived_from'):
            row[column] = json.loads(row.pop(column+'_json'))
        row['derived_from'] = tuple(row['derived_from'])
        if pd.isna(row['uncertainty']): row['uncertainty'] = None
        observation = Observation(**row)
        if observation.record_id != record_id:
            raise ValueError('observation content does not match its stored identity')
        records.append(observation)
    return records


def from_gene_table(nodes, organism, source_version, columns=None, include_missing=False):
    """Yield source-linked gene summaries without inventing unavailable assay detail.

    The existing table is already aggregated. Original replicate, dose and time
    information cannot be recovered from a summary column. Those fields should
    enter as separate Observation records when raw source tables are imported.
    """
    from .datasets import provenance
    from .slots import all_slots, declared_columns
    contexts = {}
    for slot in all_slots(organism):
        if slot.unit == 'gene':
            for col in declared_columns(nodes,slot):
                contexts.setdefault(col,set()).add(slot.context)
    for column in columns or [c for c in nodes if c!='gene_id']:
        source = provenance(column)
        source_id = source.key if source else 'unregistered_column:'+column
        status = ('orthology_transfer' if column.startswith('ortholopit_') else
                  'computed' if source and (source.derived_from or 'computed' in source.kind.lower()) else 'unclassified')
        for gene,value in zip(nodes.gene_id,nodes[column]):
            if isinstance(value,np.generic): value=value.item()
            missing = value is None or (np.isscalar(value) and pd.isna(value))
            if missing and not include_missing: continue
            if not missing and not isinstance(value,(str,float,int,bool)):
                continue  # structured rows require explicit entity/trait interpretation
            yield Observation(str(gene),column,None if missing else value,source_id,organism,
                              source_version=source_version,source_url=source.url or '' if source else '',
                              source_location='column:'+column,missing_state='unknown' if missing else 'observed',
                              evidence_status=status,derived_from=source.derived_from if source else (),
                              context={'resolution':'existing_gene_summary','slot_contexts':sorted(contexts.get(column,()))})


@dataclass(frozen=True)
class LiteratureAssertion:
    """A sourced claim proposed for review, distinct from a measured observation."""
    entity_id: str
    relation: str
    object: str
    source_id: str
    source_url: str
    passage: str
    source_location: str = ""
    context: dict = field(default_factory=dict)
    negated: bool = False
    review_status: str = "pending"
    reviewer: str = ""
    reviewed_at: str = ""

    def __post_init__(self):
        if not all((self.entity_id,self.relation,self.object,self.source_id,self.passage)):
            raise ValueError('assertions require an entity, relation, object, source and passage')
        if self.review_status not in {'pending','accepted','rejected'}:
            raise ValueError('invalid review status')
        if self.review_status!='pending' and not (self.reviewer and self.reviewed_at):
            raise ValueError('a reviewed assertion must name the reviewer and review time')

    def review(self, accepted, reviewer):
        """Return a reviewed copy; retain the source passage and original claim."""
        from dataclasses import replace
        return replace(self,review_status='accepted' if accepted else 'rejected',reviewer=reviewer,
                       reviewed_at=datetime.now(timezone.utc).isoformat())

    def observation(self, organism):
        """Convert an accepted claim to a sourced literature observation; pending claims cannot enter features."""
        if self.review_status!='accepted':
            raise ValueError('review and accept a literature assertion before using it as evidence')
        return Observation(self.entity_id,'literature:'+self.relation,
                           {'object':self.object,'negated':self.negated},self.source_id,organism,
                           source_url=self.source_url,source_location=self.source_location,
                           context={**self.context,'reviewer':self.reviewer,'reviewed_at':self.reviewed_at,
                                    'evidence_form':'reviewed_literature_claim'},evidence_status='computed')


def propose_assertions(passages, gene_index):
    """Extract simple relation phrases as review candidates, never accepted facts.

    Each input supplies text, source_id, source_url and optional location/context.
    Ambiguous multi-gene sentences are omitted. This conservative rule extractor
    does not resolve coreference, infer missing contexts or establish truth.
    """
    pattern=re.compile(r'\b(locali[sz](?:es|ed)\s+(?:to|in)|required for|interacts with)\s+([^.;!?]{3,100})',re.I)
    proposed=[]
    for record in passages:
        for sentence in re.split(r'(?<=[.!?])\s+',record['text']):
            genes=sorted({gene for gene,_ in gene_index.find(sentence)})
            if len(genes)!=1: continue
            for match in pattern.finditer(sentence):
                before=sentence[max(0,match.start()-40):match.start()].lower()
                proposed.append(LiteratureAssertion(genes[0],match[1].lower(),match[2].strip(),
                    record['source_id'],record.get('source_url',''),sentence,
                    source_location=record.get('location',''),context=record.get('context',{}),
                    negated=bool(re.search(r'\b(not|never|no evidence)\b',before))))
    return proposed


def save_assertions(assertions, path):
    """Save review status and evidence passages as portable JSON lines."""
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(''.join(json.dumps(asdict(a),sort_keys=True)+'\n' for a in assertions))


def load_assertions(path):
    """Read and validate assertions without changing their review status."""
    return [LiteratureAssertion(**json.loads(line)) for line in Path(path).read_text().splitlines() if line.strip()]
