#!/usr/bin/env python3
"""Turn a PRIDE deposit's own search output into one number per gene.

## What this reads, and what it refuses to

The submitter's result tables: MaxQuant `*Sites.txt` and `evidence.txt`, and mzIdentML `.mzid`. Not
the raw spectra -- see `scripts/fetch_candidates.py` for why that matters, but briefly: a `.raw` is
an instrument recording, and turning one into protein identifications needs a search engine, a
sequence database and a set of parameters. Those choices are a study of their own, and inventing
them here would produce numbers that look like measurements and are not.

So every column this produces is a count of what the AUTHORS reported, and the docstring of each
slot it fills should say so. "How many lactylation sites did this study report on this gene" is a
smaller claim than "how many lactylation sites does this gene have", and it is the one the data
supports.

## Why counts rather than intensities

Intensity is comparable within one experiment and not across two. These deposits come from
different laboratories, instruments and normalisation choices, and a column mixing them would rank
genes by whose mass spectrometer was more sensitive. A site count is coarse, robust to all of that,
and answers the question the slot actually asks -- is this gene modified this way, and heavily.

Genes absent from a deposit get NaN rather than zero, and the difference is the whole of what makes
the column honest: a gene the study never saw is not a gene with no sites. `na_policy="indicator"`
downstream then treats "not measured" as its own signal, which it is.
"""
from __future__ import annotations

import os
import re
import subprocess

import numpy as np
import pandas as pd

#: Toxoplasma gene identifiers, in either of the two strain namespaces the deposits use.
GENE = re.compile(r"\b(TGME49_\d{6}|TGGT1_\d{6})\b")

#: Files that carry per-site or per-peptide identifications. `Sites` tables are the direct answer;
#: `evidence`, `proteinGroups` and `.mzid` are the fallback when a deposit published no sites table.
#: Inside an archive the net is wider, because a submitter who zips a MaxQuant `txt` folder may
#: rename what is in it -- one of these deposits ships `txt.zip` containing a single `txt.txt`.
RESULT_FILES = re.compile(r"(Sites?\.txt|evidence\.txt|proteinGroups\.txt|\.mzid)$", re.I)
ARCHIVE_FILES = re.compile(r"\.(txt|tsv|csv|mzid)$", re.I)

#: Never counted, however many genes it names. A FASTA is the SEARCH DATABASE -- it lists every
#: gene in the proteome by construction, so counting it reports the genome rather than the
#: experiment. It inflated one verification from 128 genes to 7,701 before this line existed.
NOT_EVIDENCE = re.compile(r"\.(fasta|fa|faa)$|parameters|summary", re.I)


def _members(path: str) -> list:
    """What is inside an archive, whatever kind it is."""
    try:
        listing = subprocess.run(["bsdtar", "-tf", path], capture_output=True, text=True,
                                 timeout=300)
    except (OSError, subprocess.SubprocessError):
        return []
    return [name for name in listing.stdout.splitlines() if not name.endswith("/")]


def _extract(path: str, member: str, limit: int = 80_000_000) -> str:
    """One member of an archive, as text."""
    try:
        got = subprocess.run(["bsdtar", "-xOf", path, member], capture_output=True, timeout=600)
    except (OSError, subprocess.SubprocessError):
        return ""
    return got.stdout[:limit].decode("utf8", "replace")


def gene_counts(text: str) -> pd.Series:
    """How many times each gene is named in a result table.

    Counted from the identifiers themselves rather than from a named column, because these tables
    do not agree on what that column is called: MaxQuant writes `Proteins`, `Leading proteins` and
    `Protein`; an mzIdentML puts the accession in an attribute. What every one of them does is spell
    the gene, so that is what is counted -- one row of a sites table is one reported site.
    """
    if not text:
        return pd.Series(dtype=float)
    rows = text.splitlines()
    counts: dict[str, int] = {}
    for row in rows:
        # A contaminant row names a human or bovine protein and no Toxoplasma gene, so it drops out
        # by finding nothing rather than by being filtered -- MaxQuant's CON__ prefix is not
        # consistent across versions and matching on it would be the fragile way to do this.
        for gene in set(GENE.findall(row)):
            counts[gene] = counts.get(gene, 0) + 1
    return pd.Series(counts, dtype=float)


#: A MaxQuant `txt` folder holds one Sites table per modification the search looked for, and most
#: of them are not what the slot asked about. `Carbamidomethyl (C)Sites` is an artefact of sample
#: preparation -- every cysteine gets it -- `Deamidation (NQ)Sites` is a different modification, and
#: `allPeptides` is every peptide seen at all. Counting the folder wholesale put 90% of the proteome
#: in the S-nitrosylation slot, which is what a modification measured on nearly every gene should
#: always look like: a bug.
GENERIC_TABLES = re.compile(r"allPeptides|matchedFeatures|ms\d*Scans|msms|peptides\.txt"
                            r"|modificationSpecificPeptides|Carbamidomethyl|Deamidation"
                            r"|Oxidation", re.I)


def deposit_counts(folder: str, limit_files: int = 40, wants: str = "") -> pd.Series:
    """Sites per gene for one downloaded deposit, across the tables that measure `wants`.

    `wants` is the modification the SLOT asked for, and a deposit's own file names are what decide
    whether a table measures it. Without that, a folder of eleven Sites tables answers every slot
    with all eleven.
    """
    total = pd.Series(dtype=float)
    if not os.path.isdir(folder):
        return total
    for name in sorted(os.listdir(folder)):
        path = os.path.join(folder, name)
        if name.startswith("GSE") or not os.path.isfile(path):
            continue
        if NOT_EVIDENCE.search(name):
            continue
        texts = []
        if name.lower().endswith((".rar", ".zip", ".tar", ".tgz")):
            for member in _members(path)[:limit_files]:
                if not ARCHIVE_FILES.search(member) or NOT_EVIDENCE.search(member):
                    continue
                if GENERIC_TABLES.search(member):
                    continue
                if wants and not re.search(wants, member, re.I):
                    continue
                texts.append(_extract(path, member))
        elif name.lower().endswith(".gz"):
            # A lone .gz is one compressed file, not an archive with members: `bsdtar -tf` lists
            # nothing and the deposit reads as empty. `Results.mzid.gz` is exactly this.
            import gzip
            try:
                with gzip.open(path, "rb") as fh:
                    texts.append(fh.read(80_000_000).decode("utf8", "replace"))
            except OSError:
                pass
        elif RESULT_FILES.search(name):
            with open(path, "rb") as fh:
                texts.append(fh.read(80_000_000).decode("utf8", "replace"))
        for text in texts:
            total = total.add(gene_counts(text), fill_value=0.0)
    return total


#: The verified deposits, what each fills, and which tables inside measure that modification. The
#: `wants` pattern is the difference between a slot counting its own modification and a slot
#: counting every modification the search happened to look for.
#:
#: PXD017032 is deliberately NOT here. It was ingested as a kinase-substrate deposit and is neither:
#: it is Wang 2022's sporulated-oocyst vs tachyzoite phosphoproteome, already in the registry as
#: `phospho_quantitative` from the authors' own supplement. Counted from the raw deposit it produced
#: a column correlating at rho = 0.729 with `phospho_sites_measured` over 1,592 shared genes -- one
#: experiment entering the map twice, under a slot whose own citations name CDPK1 and CDPK7.
#:
#: Verified here means the file itself was read and found to name Toxoplasma genes AND the
#: modification -- not that a metadata field agreed. Of 102 datasets proposed from GEO, 8 survived
#: that check; the rest are in the same quarantine directory and are deliberately not listed.
DEPOSITS = (
    ("Tg", "acetylation", "PXD079431", "n_acetylation_sites", r"acetyl|GCN5|FLAG"),
    ("Tg", "interaction_proximity_labelling", "PXD059579", "n_proximity_partners", r"mzid|Results"),
    ("Tg", "S_nitrosylation", "PXD046083", "n_nitrosylation_sites", r"iodo ?TMT|nitrosyl|SNO"),
    ("Tg", "ubiquitination_SUMOylation", "PXD042937", "n_ubiquitination_sites",
     r"GlyGly|ubiquitin"),
    # `La (K)Sites` is MaxQuant's name for the lactylation search. This is the SECOND lactylation
    # deposit tried: PXD022700 ships a RAR that bsdtar cannot open, and the "0 genes" that produced
    # was a fact about the reader rather than about the study, which names 537 proteins. Both are
    # downloaded; this is the one that can be read.
    ("Tg", "lactylation", "PXD031526", "n_lactylation_sites", r"La \(K\)Sites"),
    # An AAL-lectin pulldown, so the count is peptide identifications and not sites -- named for that.
    # The same standing as the proximity-labelling column: both are "what came down in this pulldown",
    # which is a claim about enrichment rather than about a residue.
    ("Tg", "glycosylation", "PXD004426", "n_o_fucosyl_peptides", r"mzid|AAL"),
)

#: Where the fetcher puts things. Reading from quarantine is deliberate: a deposit is not promoted
#: by being downloaded, only by being checked, and `DEPOSITS` is the record of which were.
QUARANTINE = os.path.join("datasets", "quarantine", "2026_08_16_pride")


def load_all(base: str, index, log=print, resolve=None) -> pd.DataFrame:
    """Every verified deposit as columns, aligned to a gene index.

    The differentiation reporter screen rides along here rather than in `screens.crispr_screens`,
    which reads the published fitness tables from the dataset archive. This one is in quarantine and
    is computed from counts, so it belongs with the other deposits that were verified by being read.
    """
    from . import screens
    out = pd.DataFrame(index=pd.Index(index, dtype=object))
    for table in (screens.second_background_fitness(
                      os.path.join(base, QUARANTINE, "Tg",
                                   "essentiality_in_a_second_background"), log=log,
                      resolve=resolve),
                  screens.differentiation_screen(
                      os.path.join(base, QUARANTINE, "Tg",
                                   "essentiality_in_a_second_background"), log=log),
                  screens.thermal_shift(
                      os.path.join(base, "datasets", "quarantine", "2026_08_16_unverified", "Tg",
                                   "thermal_shift"), log=log, resolve=resolve),
                  screens.cyst_wall_interactome(
                      os.path.join(base, "datasets", "quarantine", "2026_08_16_unverified", "Tg",
                                   "cyst_wall"), log=log, resolve=resolve),
                  screens.oxidative_stress_screen(
                      os.path.join(base, "datasets", "quarantine", "2026_08_16_unverified", "Tg",
                                   "fitness_oxidative_stress"), log=log, resolve=resolve)):
        if table.empty:
            continue
        for column in table.columns:
            out[column] = table[column].reindex(out.index)
    for organism, folder, accession, column, wants in DEPOSITS:
        where = os.path.join(base, QUARANTINE, organism, folder)
        if not os.path.isdir(where):
            log(f"proteomics: {accession} not downloaded, {column} left out")
            continue
        got = column_for(where, column, out.index, wants=wants, resolve=resolve)
        out[column] = got[column]
        log(f"proteomics: {column} from {accession}, "
            f"{int(got[column].notna().sum()):,} genes measured")
    return out


def column_for(folder: str, column: str, index, wants: str = "",
               resolve=None) -> pd.DataFrame:
    """One deposit as one column, aligned to a gene index, NaN where the study saw nothing.

    NaN and not zero. A gene this deposit never reported is a gene nobody measured for this
    modification, and writing zero would state that it carries none -- a claim the data cannot make
    and one that would put thousands of unmeasured genes at the bottom of every ranking.

    Accessions go through `resolve` before the join. A deposit keyed on the type I strain writes
    `TGGT1_273760` where the node table says `TGME49_`, and matching those as strings joins nothing
    -- silently, because a join that matches no rows looks exactly like a study with no coverage.
    The lactylation deposit is entirely TGGT1_ and reported 0 of its 524 genes until this existed.
    """
    counts = deposit_counts(folder, wants=wants)
    if resolve is not None and not counts.empty:
        mapped = pd.Series([resolve(g) or g for g in counts.index], index=counts.index)
        counts = counts.groupby(mapped.to_numpy()).sum()
    values = pd.Series(np.nan, index=pd.Index(index, dtype=object), dtype=float)
    if counts.empty:
        return values.to_frame(column)
    seen = counts.reindex(values.index)
    values[seen.notna()] = seen[seen.notna()]
    return values.to_frame(column)
