"""Conservative source and derived-input closure for frozen functional benchmarks.

This opt-in benchmark policy supplements training-only empirical closure. Domain
node annotations and domain-built edges occupy different registry families, so a
source-family match alone cannot protect the held-out target. Research attention
and unresolved functional/homology summaries are also withheld. This changes no
installed strategy, source annotation, calibration or historical frozen artifact.

Callers must declare provenance for precomputed columns/layers. An opaque fitted
representation cannot become independent merely because its name is unfamiliar.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import re

from . import datasets, discovery_labels as D, strategies as S
from .splits import make_exclusions as training_exclusions

POLICY_VERSION = 'functional-source-closure-v2'
SOURCE_COLUMNS = frozenset((*D.FUNCTION_FIELDS,'has_ec','has_domain','n_interpro','n_pfam',
    'n_domains','domain_count','orthogroup','paralog_number','has_pf_ortholog',
    'has_cp_ortholog','lineage_specific','n_publications','n_fulltext'))
SOURCE_LAYERS = frozenset(('domain','orthogroup',*S.LITERATURE_LAYERS))
RELATED_COLUMN = re.compile(r'^(?:interpro|pfam|ec_number)(?:_|$)|^(?:domain_|n_domain|n_interpro|n_pfam)'
    r'|^(?:n_publications|n_fulltext|publication_count|fulltext_count)(?:_|$)')
ATTENTION_COLUMN = re.compile(
    r'^(?:n_)?(?:publications?|full_?texts?|citations?|papers?)(?:_|$)', re.IGNORECASE)


def _address(kind, name):
    if kind not in {'column','layer'} or not isinstance(name,str) or not name or name!=name.strip():
        raise ValueError('Derived inputs require an explicit column/layer address')


@dataclass(frozen=True)
class DerivedInput:
    """One declared representation/operator input with immutable source addresses."""

    kind: str
    name: str
    parents: tuple[tuple[str,str], ...]

    def __post_init__(self):
        _address(self.kind,self.name)
        if not isinstance(self.parents,tuple) or not self.parents:
            raise ValueError('Derived input requires immutable nonempty parent addresses')
        for parent in self.parents:
            if not isinstance(parent,tuple) or len(parent)!=2:
                raise ValueError('Derived parent must be a typed address tuple')
            _address(*parent)
        if len(set(self.parents))!=len(self.parents):
            raise ValueError('Duplicate derived input parent')


def make_exclusions(nodes, targets, *, source_targets, benchmark_id, split, derived_inputs=()):
    """Bind training-only source closure and explicit functional lineage to a split.

    Source targets must name actual installed EC/InterPro/Pfam fields. ``targets``
    additionally names the encoded benchmark question. Every named field must be
    present in the training table; held-out labels cannot select exclusions.
    Known annotation/homology/literature layers are withheld regardless of their
    registry family, and declared downstream columns/layers close transitively.
    The result is the existing immutable ExclusionManifest, usable by its guards.
    """
    if isinstance(targets,str):targets=(targets,)
    if isinstance(source_targets,str):source_targets=(source_targets,)
    if not targets or any(not isinstance(target,str) or not target for target in targets):
        raise ValueError('Declare a nonempty encoded benchmark target')
    if not source_targets or not set(source_targets)<=D.FUNCTION_FIELDS:
        raise ValueError('Declare actual functional annotation source fields')
    if not isinstance(derived_inputs,tuple) or not all(isinstance(item,DerivedInput) for item in derived_inputs):
        raise ValueError('Derived input declarations must be immutable typed tuples')
    addresses=[(item.kind,item.name) for item in derived_inputs]
    if len(set(addresses))!=len(addresses):raise ValueError('Duplicate derived input address')
    named=tuple(dict.fromkeys((*source_targets,*targets)))
    manifest=training_exclusions(nodes,named,benchmark_id=benchmark_id,split=split)
    columns=set(manifest.columns)|set(SOURCE_COLUMNS)
    columns.update(column for column in nodes if S.NEVER_FEATURES.search(column)
                   or RELATED_COLUMN.search(column) or ATTENTION_COLUMN.search(column))
    # Supplemental source bans need the same registry closure as the initial target.
    columns.update(datasets.derived_dependents(columns))
    layers=set(manifest.layers)|set(SOURCE_LAYERS)
    lineage=list(derived_inputs)
    lineage.extend(DerivedInput('layer',name,tuple(('layer',parent) for parent in parents))
        for name,parents in S.DERIVED_LAYERS.items())
    parents={}
    for item in lineage:
        parents.setdefault((item.kind,item.name),set()).update(item.parents)
    terminals={('column',name) for name in nodes}|{('column',name) for name in SOURCE_COLUMNS}
    terminals|={('layer',name) for name in SOURCE_LAYERS}
    terminals|={('layer',name) for names in S.DERIVED_LAYERS.values() for name in names}
    visiting=set();visited=set()
    def validate_lineage(address):
        """Require acyclic lineage rooted in named source columns or known layers."""
        if address in visiting:raise ValueError('Cyclic derived input provenance')
        if address in visited:return
        if address not in parents:
            if address not in terminals:raise ValueError('Unresolved derived input source: '+str(address))
            return
        visiting.add(address)
        for parent in parents[address]:validate_lineage(parent)
        visiting.remove(address);visited.add(address)
    for address in parents:validate_lineage(address)
    changed=True
    while changed:
        changed=False
        for item in lineage:
            excluded=columns if item.kind=='column' else layers
            if item.name in excluded:continue
            if any(name in (columns if kind=='column' else layers) for kind,name in item.parents):
                excluded.add(item.name);changed=True
    return replace(manifest,columns=tuple(sorted(columns)),layers=tuple(sorted(layers)))
