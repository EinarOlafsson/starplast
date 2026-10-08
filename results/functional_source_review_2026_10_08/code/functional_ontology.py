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
        malformed = any(D.EC.match(entry) is None for entry in source_entries)
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
