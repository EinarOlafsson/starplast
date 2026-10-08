"""Discover every available annotation, keeping known membership separate from claims.

Functional domains and enzyme classes are multi-valued annotations, not exclusive
localizations. Missing annotation is unknown membership, not a biological negative.
This catalogue reports inference/test availability without generating confidence.
"""
from __future__ import annotations

import re

import pandas as pd

from . import datasets as D, functional_domains as FD, strategies as S

FUNCTION_FIELDS = {'interpro_id','interpro_ids','pfam_id','pfam_ids','ec_number','ec_number_orthology'}
ANNOTATION_FLAGS = re.compile(r'^(?:has_|is_)|(?:_positive|_hit|_stringent|_member|_confident|_agree|_accepted|_dependent|_newly_made)$|^lineage_specific$')
IPR = re.compile(r'\bIPR\d{6}\b')
PFAM = re.compile(r'\bPF\d{5}(?:\.\d+)?\b')
EC = re.compile(r'(?<![\d.])(?:[1-7]|-)(?:\.(?:\d+|-)){3}(?![\d.])')
FAMILY_ORDER = {'Function':0,'Phenotype':1,'Stage / expression':2,'Localization':3,'Structure':4,'Other labels':5}
CATALOGUE_COLUMNS = ('organism','target','title','family','genes','annotated_genes','unannotated_genes',
    'classes','annotation_coverage','source_annotated_genes','unparsed_annotation_genes','claims','tested_claims',
    'held_out_rows','evaluated_strategies','evaluation_status','source_ids','truth_interpretation','annotation_warning','biological_accuracy')


def _ec_entries(text):
    """Split source entries outside descriptive parentheses, retaining internal semicolons."""
    depth,start,entries = 0,0,[]
    for i,char in enumerate(text):
        if char=='(':depth+=1
        elif char==')':depth=max(0,depth-1)
        elif char==';' and depth==0:
            entries.append(text[start:i].strip());start=i+1
    return [entry for entry in (*entries,text[start:].strip()) if entry]


def label_family(target: str) -> str:
    """Describe the question represented by a column without turning it into an inference."""
    if target in FUNCTION_FIELDS or target in {'has_domain','has_ec'}:return 'Function'
    if 'phenotype' in target or target.startswith(('screen_','drug_','resistance_','enteric_')):return 'Phenotype'
    if 'stage' in target or 'cellcycle' in target:return 'Stage / expression'
    if any(part in target for part in ('compartment','lopit','location','exported')):return 'Localization'
    if target in {'dtm_class','has_structure','is_tm','has_signal_peptide'}:return 'Structure'
    return 'Other labels'


def annotation_members(ctx, target: str, *, domain_lookup=None) -> pd.DataFrame:
    """One row per organism/label/class/gene; preserve overlaps and unknown absence.

    EC identifiers are parsed independently of their optional source description.
    InterPro descriptions are paired only when source ID/name positions agree.
    Each EC description stays attached to its own assigned source entry; numbers
    mentioned in a replacement/deprecation note are not additional memberships.
    """
    columns = ['organism','target','value','description','gene_id']
    if domain_lookup is not None:
        columns += ['description_source','ontology_release','ontology_status','accession_version_status']
    if target not in ctx.nodes:return pd.DataFrame(columns=columns)
    rows = []
    values = ctx.truth(target)
    description_field = 'pfam_desc' if target.startswith('pfam') else 'interpro_desc'
    descriptions = ctx.nodes.get(description_field,pd.Series('',index=ctx.nodes.index))
    for i,text in values.dropna().items():
        if target in {'interpro_id','interpro_ids','pfam_id','pfam_ids'}:
            terms = list(dict.fromkeys((PFAM if target.startswith('pfam') else IPR).findall(str(text))))
            names = str(descriptions.iloc[i]).split(';') if pd.notna(descriptions.iloc[i]) else []
            raw_terms = str(text).split(';')
            paired = dict(zip((term.strip() for term in raw_terms),names)) if len(raw_terms)==len(names) else {}
            labels = [(term,paired.get(term,'')) for term in terms]
        elif target.startswith('ec_number'):
            labels = []
            for entry in _ec_entries(str(text)):
                entry = entry.strip()
                match = EC.match(entry)
                if match:labels.append((match.group(),entry))
        else:labels = [(str(text),'')]
        for term,description in labels:
            row = [ctx.organism,target,term,description,str(ctx.gene_ids[i])]
            if domain_lookup is not None:
                metadata = domain_lookup.lookup(term) if target in {'interpro_id','interpro_ids','pfam_id','pfam_ids'} else {}
                description_source = 'Original source annotation' if description.strip() else ''
                if not description.strip() and metadata.get('status')=='current_metadata':
                    row[3] = metadata['name']
                    description_source = metadata['source']+' current nomenclature'
                row += [description_source,metadata.get('ontology_release',''),metadata.get('status',''),
                    metadata.get('accession_version_status','')]
            rows.append(row)
    return pd.DataFrame(rows,columns=columns).drop_duplicates(['organism','target','value','gene_id']).reset_index(drop=True)


def catalogue(ctx, claims, recipes, records=None) -> tuple[pd.DataFrame,pd.DataFrame]:
    """All available categorical/function annotations plus every recorded claim recipe.

    Coverage is known annotation membership, not accuracy. Proven recipes retain
    their legacy evaluation meaning; absent or uncalibrated recipes never gain a
    score from having many annotations. Source ownership is registry-declared,
    with unrecorded lineage explicitly unresolved.
    """
    records = pd.DataFrame() if records is None else records
    try:domain_lookup = FD.shipped()
    except (OSError,ValueError,KeyError,TypeError):domain_lookup = None
    record_counts = (records.groupby('target').agg(held_out_rows=('target','size'),
        evaluated_strategies=('strategy','nunique')).to_dict('index') if len(records) else {})
    targets = (set(ctx.categorical_columns(min_labelled=1)) - {'interpro_desc','pfam_desc'}) | (FUNCTION_FIELDS & set(ctx.nodes))
    for column in ctx.nodes:
        series = ctx.nodes[column]
        if S.NEVER_FEATURES.search(column) or column in {'interpro_desc','pfam_desc'}:
            continue
        if pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series):
            if not ANNOTATION_FLAGS.search(column) or not set(series.dropna().unique()) <= {0,1}:continue
        if 0<ctx.truth(column).nunique()<=60:targets.add(column)
    for frame in (claims,recipes,records):
        if len(frame) and 'target' in frame:targets.update(frame.target.astype(str))
    members,rows = [],[]
    for target in sorted(targets,key=lambda t:(FAMILY_ORDER[label_family(t)],t)):
        membership = annotation_members(ctx,target,domain_lookup=domain_lookup)
        members.append(membership)
        inferred = claims[claims.target.astype(str).eq(target)] if len(claims) else claims
        rec = recipes[recipes.target.astype(str).eq(target)] if len(recipes) else recipes
        held = record_counts.get(target,{'held_out_rows':0,'evaluated_strategies':0})
        known = membership.gene_id.nunique()
        source_known = int(ctx.truth(target).notna().sum()) if target in ctx.nodes else 0
        owners = sorted(d.key for d in D.REGISTRY if d.organism==ctx.organism and target in d.columns)
        status = ('Recipe calibration incomplete' if len(rec) and not rec.proven.astype(bool).any()
            else 'Legacy recipe evaluated' if len(rec) else 'Inference not evaluated')
        if held['held_out_rows']:status = 'Legacy held-out records; '+status.lower()
        warning = ('Legacy EC replacement notes are descriptions, not additional class assignments; current ontology admission pending'
            if target.startswith('ec_number') and target in ctx.nodes and ctx.truth(target).str.contains('Transferred entry:',regex=False,na=False).any() else '')
        if target in {'interpro_id','interpro_ids','pfam_id','pfam_ids'}:
            if domain_lookup is None:
                warning = 'Current domain nomenclature unavailable; original identifiers and annotations retained.'
            else:
                unique = membership.drop_duplicates('value')
                current = int(unique.ontology_status.eq('current_metadata').sum())
                withdrawn = int(unique.ontology_status.eq('withdrawn').sum())
                release = ' / '.join(sorted(set(unique.ontology_release)-{''}))
                warning = (f'Current nomenclature ({release}): {current:,}/{len(unique):,} identifiers named, '
                    f'{withdrawn:,} withdrawn; remaining identifiers have unknown current status. '
                    'Original source descriptions take precedence; additional names describe identifiers, '
                    'not verified gene function or the release that assigned the gene. No replacement memberships applied.')
        rows.append({'organism':ctx.organism,'target':target,'title':target.replace('_',' '),
            'family':label_family(target),'genes':ctx.n,'annotated_genes':int(known),'unannotated_genes':int(ctx.n-known),
            'classes':int(membership.value.nunique()),'annotation_coverage':known/ctx.n if ctx.n else None,
            'source_annotated_genes':source_known,'unparsed_annotation_genes':int(source_known-known),
            'claims':int(len(inferred)),'tested_claims':int(inferred.status.astype(str).eq('tested').sum()) if len(inferred) else 0,
            'held_out_rows':int(held['held_out_rows']),'evaluated_strategies':int(held['evaluated_strategies']),
            'evaluation_status':status,'source_ids':'; '.join(owners) or 'unresolved',
            'truth_interpretation':'Source annotation membership; unannotated is unknown, not a verified negative',
            'annotation_warning':warning,
            'biological_accuracy':'Not independently admitted'})
    combined = pd.concat(members,ignore_index=True) if members else pd.DataFrame(columns=['organism','target','value','description','gene_id'])
    return pd.DataFrame(rows,columns=CATALOGUE_COLUMNS),combined


def class_summary(members: pd.DataFrame, target: str) -> pd.DataFrame:
    """Known class sizes and source descriptions; no prediction precision invented."""
    subset = members[members.target.eq(target)]
    rows = []
    for value,group in subset.groupby('value',sort=True):
        descriptions = sorted(set(group.description.dropna().astype(str))-{'','nan'})
        rows.append({'value':value,'description':' / '.join(descriptions),'annotated_genes':group.gene_id.nunique(),
            'precision':None,'recall':None,'evaluation_status':'Class inference not independently evaluated'})
    return pd.DataFrame(rows,columns=['value','description','annotated_genes','precision','recall','evaluation_status'])
