"""Present the reviewed human source foundation as evidence quality, never accuracy.

This adapter uses only the supplied foundation summary. Gene preservation and
protein-link counts describe source handling, not measured biological accuracy.
Missing source versions, independent tests and deployment remain unavailable;
the adapter neither registers a host space nor admits data for redistribution.
"""
from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re

from . import organisms as O

COUNT_FIELDS=('source_records','canonical_genes','qualified_PAR_Y_records',
    'unmapped_genes','single_protein_genes','multiple_protein_genes',
    'source_genes_with_shared_protein','reference_mapping_records',
    'reference_gene_protein_pairs','mapping_genes_outside_source_universe',
    'mapped_proteins_in_source_universe','reference_proteins_with_multiple_genes',
    'source_measurement_cells','canonical_gene_measurement_cells',
    'protein_measurements_projected','source_records_lost')
GATE_FIELDS=('benchmark_admitted','installed_space_registered','distributable_pack_built')
SUMMARY_SHA256='dd26db5ba4bb6855e9aabde667fe6835b03ea12da73f22ff562782b22ee3df6d'


def _validate(summary):
    if not isinstance(summary,Mapping):raise TypeError('A reviewed source-foundation summary is required')
    data=deepcopy(dict(summary))
    if data.get('organism')!=O.HUMAN or data.get('unit')!='gene':
        raise ValueError('Human source-gene foundation scope is required')
    for field in COUNT_FIELDS:
        if type(data.get(field)) is not int or data[field]<0:
            raise ValueError(field+' requires a supplied nonnegative integer count')
    if data['canonical_genes']+data['qualified_PAR_Y_records']!=data['source_records']:
        raise ValueError('Canonical and qualified source records do not reconcile')
    if data['unmapped_genes']+data['single_protein_genes']+data['multiple_protein_genes']!=data['canonical_genes']:
        raise ValueError('Mapping categories do not reconcile to canonical source genes')
    if (data['source_genes_with_shared_protein']>data['canonical_genes']
            or data['canonical_gene_measurement_cells']>data['source_measurement_cells']):
        raise ValueError('Source-subset counts exceed their declared source population')
    if any(data.get(field) is not False for field in GATE_FIELDS):
        raise ValueError('This foundation-only adapter cannot assert benchmark, space or pack admission')
    if data['protein_measurements_projected']!=0 or data['source_records_lost']!=0:
        raise ValueError('Reviewed foundation scope requires no protein projection or source-record loss')
    if data.get('redistribution')!='unresolved':
        raise ValueError('Redistribution admission requires a separate reviewed contract')
    if not isinstance(data.get('source_sha256'),str) or not re.fullmatch(r'[a-f0-9]{64}',data['source_sha256']):
        raise ValueError('The original source SHA-256 must remain explicit')
    if (not isinstance(data.get('universe_scope'),str) or not data['universe_scope']
            or not isinstance(data.get('gaps'),list) or not data['gaps']
            or any(not isinstance(gap,str) or not gap for gap in data['gaps'])):
        raise ValueError('Source scope and reviewed limitations must remain explicit')
    return data


def build(summary, *, summary_sha256=None):
    """Return a shared evidence card retaining reviewed counts and unfulfilled gates.

    Only source-summary facts are used. No performance task, biological accuracy,
    confidence or mapping success percentage is invented. Future completed host
    admission needs another adapter backed by its own reviewed source contracts.
    """
    from .scorecard_view import build_scorecard_view
    data=_validate(summary)
    if summary_sha256 is not None and (not isinstance(summary_sha256,str) or not re.fullmatch(r'[a-f0-9]{64}',summary_sha256)):
        raise ValueError('Summary provenance requires an explicit SHA-256 or unknown')
    counts={field:data[field] for field in COUNT_FIELDS}
    evidence={'source_preservation':{key:data[key] for key in ('source_records','canonical_genes',
        'qualified_PAR_Y_records','source_measurement_cells','canonical_gene_measurement_cells','source_records_lost')},
        'reference_mapping':{key:data[key] for key in ('unmapped_genes','single_protein_genes',
            'multiple_protein_genes','source_genes_with_shared_protein','reference_mapping_records',
            'reference_gene_protein_pairs','mapping_genes_outside_source_universe',
            'mapped_proteins_in_source_universe','reference_proteins_with_multiple_genes')},
        'protein_measurements_projected':data['protein_measurements_projected'],
        'admission':{key:data[key] for key in (*GATE_FIELDS,'redistribution')},
        'scope':data['universe_scope'],
        'reviewed_context_and_admission_gaps':list(data['gaps']),
        'biological_accuracy':None,'independent_benchmark':None,'calibrated_confidence':None,
        'mapping_interpretation':'Reference links preserve ambiguity and shared targets; they are not biological validation or a license to project measurements',
        'source_interpretation':'Source records retained without protein projection; qualified Y records remain separate. Source gene annotation is not genome completeness'}
    source={'grade':'unresolved','name':'Reviewed human source-gene foundation',
        'version':None,'sha256':data['source_sha256'],'context':data['universe_scope'],
        'lineage':'Original gene records and mapping summary; no protein measurements projected; independent source/benchmark admission unresolved',
        'negative_semantics':'Unmapped or ambiguous reference links and missing measurements are unknown; they do not establish biological absence'}
    card={'scope':{'organism':O.HUMAN,'target':'source_gene_foundation','strategy':'evidence_audit',
        'unit':'gene','settings':'Summary-only evidence presentation; no model fitted',
        'protocol':'Reviewed source preservation and reference mapping foundation',
        'partition':'Source audit; not a held-out inference test','benchmark_id':None,
        'truth_grade':'unresolved','gaps':list(data['gaps'])},
        'counts':counts,'metrics':{},'source':source,'evidence':evidence}
    return build_scorecard_view(card,kind='evidence',title='Human source-gene foundation',
        status='Foundation reviewed; host space, benchmark and distribution unavailable',
        details={'source':source,'sizes':counts,'gaps':data['gaps'],
            'freshness':{'summary_sha256':summary_sha256,'source_version':None,'source_release_to_publication':'unresolved'},
            'split':{'status':'unavailable','reason':'This is source preservation, not an inference benchmark'},
            'baseline':None,'uncertainty':None,'controls':None,'rows':None})


def load(path, *, expected_sha256=None):
    """Read this reviewed summary only, optionally enforcing a separately pinned identity."""
    path=Path(path)
    if path.is_symlink():raise ValueError('Foundation summary must be a regular local file')
    raw=path.read_bytes();digest=hashlib.sha256(raw).hexdigest()
    if expected_sha256 is not None and expected_sha256!=digest:raise ValueError('Foundation summary checksum mismatch')
    return build(json.loads(raw),summary_sha256=digest)


def shipped():
    """Return the pinned summary-only host card or an explicit unavailable reason."""
    try:
        return load(Path(__file__).with_name('data')/'host_foundation_summary.json',
            expected_sha256=SUMMARY_SHA256),''
    except (OSError,ValueError,TypeError,KeyError) as exc:
        return None,'Host source summary unavailable: '+str(exc)
