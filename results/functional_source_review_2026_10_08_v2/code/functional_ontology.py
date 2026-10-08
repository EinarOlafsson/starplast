"""Interpret pinned functional nomenclature without inventing gene activity truth.

An ENZYME entry defines a reaction class; it does not prove that a parasite gene
catalyses it. Resolve only unique, active EC transfer chains. Split replacements,
deleted or preliminary entries and missing identifiers remain explicit gaps.
Complete enzyme-class profiles preserve overlapping memberships and are suitable
for declared annotation-recovery tests, not independent biological validation.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import re

import pandas as pd

from . import discovery_labels as D

EC_ID = re.compile(r'[1-7]\.\d+\.\d+\.(?:\d+|n\d+)')


def valid_ec_source(text: str) -> bool:
    """Require balanced assigned-entry syntax; a valid numeric prefix is insufficient."""
    depth=0
    for char in text:
        if char=='(':depth+=1
        elif char==')':
            depth-=1
            if depth<0:return False
    if depth:return False
    entries=D._ec_entries(text)
    for entry in entries:
        match=D.EC.match(entry)
        if match is None:return False
        tail=entry[match.end():].strip()
        if tail and not (tail.startswith('(') and tail.endswith(')')):return False
    return bool(entries)


@dataclass(frozen=True)
class EnzymeEntry:
    """One immutable ENZYME nomenclature entry with explicit obsolete-entry semantics."""

    ec: str
    name: str
    status: str
    replacements: tuple[str,...] = ()

    def __post_init__(self):
        if not EC_ID.fullmatch(self.ec) or not self.name or self.status not in {'active','preliminary','transferred','deleted'}:
            raise ValueError('Invalid explicit ENZYME entry')
        if any(not EC_ID.fullmatch(term) for term in self.replacements):
            raise ValueError('Invalid transfer destination')
        if self.status!='transferred' and self.replacements:
            raise ValueError('Only transferred entries declare replacement identities')


def parse_enzyme(text: str) -> dict[str,EnzymeEntry]:
    """Read ID/DE records in the official flat file; reject duplicates and truncation."""
    entries,identifier,description = {},None,[]
    for line in text.splitlines():
        if line.startswith('ID   '):
            if identifier is not None:raise ValueError('Unterminated ENZYME entry')
            identifier,description = line[5:].strip(),[]
        elif identifier is not None and line.startswith('DE   '):description.append(line[5:].strip())
        elif identifier is not None and line=='//':
            name = ' '.join(description)
            status = ('deleted' if name.startswith('Deleted entry') else
                'transferred' if name.startswith('Transferred entry:') else
                'preliminary' if '.n' in identifier else 'active')
            targets = tuple(dict.fromkeys(EC_ID.findall(name.split(':',1)[1]))) if status=='transferred' else ()
            entry = EnzymeEntry(identifier,name,status,targets)
            if identifier in entries:raise ValueError('Duplicate ENZYME identity')
            entries[identifier]=entry;identifier=None
    if identifier is not None:raise ValueError('Truncated ENZYME record')
    if not entries:raise ValueError('No valid ENZYME entries')
    return entries


def resolve_ec(term: str, entries: dict[str,EnzymeEntry]) -> dict:
    """Follow a unique active replacement; never expand a split into gene assignments."""
    chain,current = [],term
    while True:
        if current in chain:return {'source_ec':term,'status':'cycle','canonical_ec':None,'chain':chain+[current],'alternatives':[]}
        chain.append(current)
        entry = entries.get(current)
        if entry is None:status='missing'
        else:status=entry.status
        if status=='active':
            return {'source_ec':term,'status':'active' if len(chain)==1 else 'unique_transfer',
                'canonical_ec':current,'chain':chain,'alternatives':[]}
        if status!='transferred':
            return {'source_ec':term,'status':status,'canonical_ec':None,'chain':chain,'alternatives':[]}
        if len(entry.replacements)!=1:
            return {'source_ec':term,'status':'ambiguous_transfer' if entry.replacements else 'unresolved_transfer',
                'canonical_ec':None,'chain':chain,'alternatives':list(entry.replacements)}
        current=entry.replacements[0]


def ec_profiles(ctx, target: str, entries: dict[str,EnzymeEntry]) -> tuple[pd.DataFrame,pd.DataFrame]:
    """Preserve every source gene and every EC, including multi-class and unresolved profiles.

    A profile is eligible only if all assigned source entries parse and resolve
    uniquely. Unannotated, deleted, missing or ambiguous genes never become a
    negative class and never acquire a partial truth profile. Original annotations
    remain unchanged; canonical membership is a separate nomenclature derivation.
    """
    if target not in {'ec_number','ec_number_orthology'} or target not in ctx.nodes:
        raise ValueError('Declare an installed EC source field')
    if len(set(ctx.gene_ids))!=ctx.n:raise ValueError('Unique source genes are required')
    members = D.annotation_members(ctx,target)
    by_gene = {str(gene):group.value.tolist() for gene,group in members.groupby('gene_id',sort=False)}
    truth = ctx.truth(target)
    rows,resolutions = [],[]
    for i,gene in enumerate(ctx.gene_ids):
        assigned = by_gene.get(str(gene),[])
        source_entries = D._ec_entries(str(truth.iloc[i])) if pd.notna(truth.iloc[i]) else []
        malformed = bool(source_entries) and not valid_ec_source(str(truth.iloc[i]))
        resolved = [resolve_ec(term,entries) for term in assigned]
        for result in resolved:resolutions.append(dict(result,gene_id=str(gene),source_target=target))
        complete = bool(assigned) and not malformed and all(r['canonical_ec'] for r in resolved)
        canonical = sorted({r['canonical_ec'] for r in resolved if r['canonical_ec']})
        major = sorted({term.split('.')[0] for term in canonical}) if complete else []
        rows.append({'organism':ctx.organism,'gene_id':str(gene),'source_target':target,
            'source_ec_codes':assigned,'canonical_ec_codes':canonical,'major_classes':major,
            'profile':json.dumps(major,separators=(',',':')) if complete else None,
            'eligible':complete,'status':'complete' if complete else 'unannotated' if not source_entries else 'unresolved',
            'malformed_source_entry':malformed})
    return pd.DataFrame(rows),pd.DataFrame(resolutions,columns=['source_ec','status','canonical_ec','chain','alternatives','gene_id','source_target'])


def profile_members(value) -> set[str]:
    """Decode a complete canonical major-class profile; reject malformed or partial labels."""
    try:terms=json.loads(value)
    except (TypeError,ValueError) as exc:raise ValueError('Explicit complete enzyme-class profile required') from exc
    if not isinstance(terms,list) or not terms or any(not isinstance(term,str) or term not in {'1','2','3','4','5','6','7'} for term in terms):
        raise ValueError('Profile must contain explicit major enzyme class identifiers')
    if terms!=sorted(set(terms)):raise ValueError('Profile must be sorted and unique')
    return set(terms)


def member_class_cards(rows) -> list[dict]:
    """Recover individual classes within complete recorded profiles, retaining abstentions.

    False-call counts refer only to recovery of these recorded profiles in the
    frozen eligible cohort. They are not adjudicated biological negatives or
    accuracy estimates for unannotated genes. Joint profile calls remain native
    categorical predictions; per-class probabilities are not inferred from them.
    """
    actual=[profile_members(value) for value in rows.truth]
    predicted=[profile_members(value) if isinstance(value,str) else set() for value in rows.prediction]
    answered=[isinstance(value,str) for value in rows.prediction]
    result=[]
    for term in ('1','2','3','4','5','6','7'):
        tp=sum(term in t and term in p for t,p in zip(actual,predicted))
        fp=sum(term not in t and term in p for t,p in zip(actual,predicted))
        fn=sum(term in t and term not in p for t,p in zip(actual,predicted))
        n=len(actual);positives=tp+fn
        precision=tp/(tp+fp) if tp+fp else None
        recall=tp/positives if positives else None
        result.append({'class':term,'eligible':n,'known_positive_genes':positives,'true_positive':tp,'false_positive':fp,'false_negative':fn,
            'precision':precision,'recall':recall,'f1':2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None,
            'coverage':sum(answered)/n if n else None,'reference_prevalence':positives/n if n else None,
            'interpretation':'Recovery of recorded complete profiles; not independent biological activity accuracy',
            'biological_precision':None,'biological_recall':None,'calibrated_confidence':None})
    return result
