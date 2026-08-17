#!/usr/bin/env python3
"""The host table and the bridges to it: rows are human proteins, not parasite genes.

`interaction · with host proteins` asks about pairs whose two ends live in different organisms, and
a parasite gene table cannot hold it. Instruction 39 settled the shape -- a second table keyed by the
host's own identifier, and pair slots as the bridge -- and this is the first one built.

## What identifies a host row

The UniProt entry name (`PDCD6_HUMAN`) and accession, taken from the deposit's own protein groups.
Not a gene symbol: mapping accession to symbol needs UniProt's ID-mapping service, and the entry name
is already stable, unambiguous and sufficient to join a second study.

## What the bridge is, and what it is not

One row per (parasite gene, host protein) pair, with the evidence that put it there. It is NOT an
edge in `graph.npz`: those are index pairs into the parasite node table, and neither end of a bridge
is guaranteed to be in it.

## Reading a co-immunoprecipitation honestly

Two filters do most of the work, and both were arrived at by getting it wrong first.

A protein group is a contaminant group if ANY entry in it is a contaminant. MaxQuant prefixes the
group with `CON__` only when the LEADING entry is one, so keratin arrives in the middle of a group
headed by `sp|` and survives the obvious filter. The four most enriched host proteins in this deposit
are keratins until that is fixed.

And a protein seen by one peptide in one run is an identification, not an interaction. Requiring two
unique peptides in BOTH bait replicates takes 674 host groups down to 219.
"""
from __future__ import annotations

import os
import re

import numpy as np
import pandas as pd

HOST_TABLE = "host_proteins.parquet"
BRIDGE_TABLE = "host_bridges.parquet"

#: The deposit, its bait, and the columns that make the comparison. `bait` is a parasite gene, and
#: naming it here is what makes every row of the bridge attributable.
MYR1_IP = {
    "accession": "PXD016383",
    "bait": "TGME49_254470",
    "folder": os.path.join("datasets", "quarantine", "2026_08_16_pride", "Tg", "host_interaction"),
    "file": "proteinGroups.txt",
    "bait_columns": ("LFQ intensity M1", "LFQ intensity M2"),
    "control_columns": ("LFQ intensity R1", "LFQ intensity R2"),
    "peptide_columns": ("Unique peptides M1", "Unique peptides M2"),
}

MIN_UNIQUE_PEPTIDES = 2
ENTRY = re.compile(r"\|([A-Z0-9]+_HUMAN)")
ACCESSION = re.compile(r"\b(?:sp|tr)\|([A-Z0-9]+)\|")


def read_ip(base: str, spec: dict = None) -> pd.DataFrame:
    """One immunoprecipitation as host protein rows with their enrichment over the control."""
    spec = spec or MYR1_IP
    path = os.path.join(base, spec["folder"], spec["file"])
    if not os.path.exists(path):
        return pd.DataFrame()
    d = pd.read_csv(path, sep="\t", low_memory=False)
    if "Protein IDs" not in d.columns:
        return pd.DataFrame()
    ids = d["Protein IDs"].astype(str)
    d = d[~ids.str.contains("CON__")].copy()
    needed = (*spec["bait_columns"], *spec["control_columns"], *spec["peptide_columns"])
    if not set(needed) <= set(d.columns):
        return pd.DataFrame()
    for column in needed:
        d[column] = pd.to_numeric(d[column], errors="coerce").fillna(0.0)
    keep = np.all([d[c] >= MIN_UNIQUE_PEPTIDES for c in spec["peptide_columns"]], axis=0)
    d = d[keep]
    ids = d["Protein IDs"].astype(str)
    # A group naming a parasite accession is the bait's own side of the experiment and belongs to
    # the parasite table, not this one.
    host = d[ids.str.contains("_HUMAN") & ~ids.str.contains("TGGT1_|TGME49_")].copy()
    if host.empty:
        return pd.DataFrame()
    bait = host[list(spec["bait_columns"])].mean(axis=1)
    control = host[list(spec["control_columns"])].mean(axis=1)
    out = pd.DataFrame({
        "host_id": host["Protein IDs"].astype(str).str.extract(ENTRY, expand=False),
        "host_accession": host["Protein IDs"].astype(str).str.extract(ACCESSION, expand=False),
        "host_ip_enrichment_log2": np.log2((bait + 1.0) / (control + 1.0)).to_numpy()})
    out = out.dropna(subset=["host_id"])
    return out.groupby("host_id", as_index=False).max()


def bridges(base: str, spec: dict = None) -> pd.DataFrame:
    """Parasite gene to host protein, one row per pair, with what put it there."""
    spec = spec or MYR1_IP
    host = read_ip(base, spec)
    if host.empty:
        return pd.DataFrame()
    return pd.DataFrame({
        "gene_id": spec["bait"],
        "host_id": host["host_id"],
        "host_ip_enrichment_log2": host["host_ip_enrichment_log2"],
        "evidence": f"IP-MS {spec['accession']}"})


def load(base: str, name: str = HOST_TABLE) -> pd.DataFrame:
    """A shipped host table, or an empty frame if it has not been built."""
    path = os.path.join(base, "starplast", "data", name)
    return pd.read_parquet(path) if os.path.exists(path) else pd.DataFrame()
