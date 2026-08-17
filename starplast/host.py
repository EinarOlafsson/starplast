#!/usr/bin/env python3
"""The host table and the bridges to it: rows are human proteins, not parasite genes.

`interaction · with host proteins` asks about pairs whose two ends live in different organisms, and
a parasite gene table cannot hold it. Instruction 39 settled the shape -- a second table keyed by the
host's own identifier, and pair slots as the bridge -- and this is the first one built.

## What identifies a host row

The UniProt ACCESSION (`O75340`), with whatever readable name the deposit gave beside it. The three
deposits here name host proteins three different ways -- MaxQuant entry names (`PDCD6_HUMAN`), gene
symbols (`PDCD6`) and bare accessions (`O75340`) -- and a bridge that mixed them would report the
same protein as three, which would turn "reached by three baits" into an artefact of formatting.
The accession is the one identifier all three carry.

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

#: The EAF1 and GRA35 affinity purifications, published as a differential-abundance table rather
#: than raw intensities -- so there is nothing to filter and nothing to normalise, only a contrast to
#: read. One deposit, two baits, each against the same wild-type control.
DIA_IPS = (
    {"accession": "PXD080696", "bait": "TGME49_225160", "bait_name": "EAF1",
     "logfc": "logFC_EAF1_II_vs_WT_III", "padj": "adj.P.Val_EAF1_II_vs_WT_III"},
    {"accession": "PXD080696", "bait": "TGME49_226380", "bait_name": "GRA35",
     "logfc": "logFC_GRA35_I_vs_WT_III", "padj": "adj.P.Val_GRA35_I_vs_WT_III"},
)
DIA_FILE = os.path.join("datasets", "quarantine", "2026_08_16_escrt", "PXD080696",
                        "combined_results.csv")

#: Enrichment and significance a host protein must clear to become a bridge row. A DIA contrast
#: gives both, so both are used; the MYR1 IP above gives neither and is filtered on peptides instead.
DIA_MIN_LOGFC = 1.0
DIA_MAX_PADJ = 0.05

#: The GRA64 pulldowns, published as two replicates side by side in one sheet, each with its own
#: accession column because the two are sorted differently. A protein counts when BOTH replicates
#: enrich it -- which is the same standard the MYR1 IP meets through its peptide requirement, met
#: here through replication instead.
REPLICATE_IPS = (
    {"accession": "PMC9426488", "bait": "TGME49_202620", "bait_name": "GRA64, tachyzoite",
     "archive": "PMC9426488/PMC9426488_supplementary.zip",
     "member": "mbio.01442-22-s0009.xlsx", "sheet": "IP_Tachyzoite_Summary"},
    {"accession": "PMC9426488", "bait": "TGME49_202620", "bait_name": "GRA64, bradyzoite",
     "archive": "PMC9426488/PMC9426488_supplementary.zip",
     "member": "mbio.01442-22-s0009.xlsx", "sheet": "IP_Bradyzoite_Summary"},
    {"accession": "PMC9426488", "bait": "TGME49_202620", "bait_name": "GRA64, HFF TurboID",
     "archive": "PMC9426488/PMC9426488_supplementary.zip",
     "member": "mbio.01442-22-s0010.xlsx", "sheet": "HFF TurboID Results Summary"},
)
ESCRT_ROOT = os.path.join("datasets", "quarantine", "2026_08_16_escrt")
REPLICATE_MIN_FOLD = 1.0
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
        "host_id": host["Protein IDs"].astype(str).str.extract(ACCESSION, expand=False),
        "host_name": host["Protein IDs"].astype(str).str.extract(ENTRY, expand=False),
        "host_ip_enrichment_log2": np.log2((bait + 1.0) / (control + 1.0)).to_numpy()})
    out = out.dropna(subset=["host_id"])
    return out.groupby("host_id", as_index=False).max()


def read_dia(base: str, spec: dict) -> pd.DataFrame:
    """One bait of the DIA deposit, as host protein rows that clear both thresholds.

    Gene symbols arrive quoted as `="PDCD6"` -- a spreadsheet's way of stopping Excel reading a
    symbol as a formula or a date -- and are unwrapped rather than matched.
    """
    path = os.path.join(base, DIA_FILE)
    if not os.path.exists(path):
        return pd.DataFrame()
    d = pd.read_csv(path)
    if not {spec["logfc"], spec["padj"], "organism", "gene_symbol"} <= set(d.columns):
        return pd.DataFrame()
    host = d[d["organism"].astype(str).str.contains("Homo sapiens", na=False)].copy()
    logfc = pd.to_numeric(host[spec["logfc"]], errors="coerce")
    padj = pd.to_numeric(host[spec["padj"]], errors="coerce")
    keep = (logfc >= DIA_MIN_LOGFC) & (padj < DIA_MAX_PADJ)
    host = host[keep]
    if host.empty:
        return pd.DataFrame()
    symbol = host["gene_symbol"].astype(str).str.replace(r'^="?|"?$', "", regex=True)
    out = pd.DataFrame({"host_id": host["uniprot_id"].astype(str).to_numpy(),
                        "host_name": symbol.to_numpy(),
                        "host_ip_enrichment_log2": logfc[keep].to_numpy()})
    return out[out["host_id"].str.len() > 0].groupby("host_id", as_index=False).max()


def read_replicated(base: str, spec: dict) -> pd.DataFrame:
    """One pulldown published as two replicate blocks, kept where both agree.

    The sheet holds each replicate's accessions in its own column because they are sorted
    differently, so the two blocks are read separately and intersected rather than read row-wise.
    Reading them row-wise would pair replicate 1 of one protein with replicate 2 of another.
    """
    import zipfile

    path = os.path.join(base, ESCRT_ROOT, spec["archive"])
    if not os.path.exists(path):
        return pd.DataFrame()
    with zipfile.ZipFile(path) as archive:
        if spec["member"] not in archive.namelist():
            return pd.DataFrame()
        import io as _io
        blob = _io.BytesIO(archive.read(spec["member"]))
        # Three specs share one workbook and each names its own sheet, so a workbook that carries
        # only some of them must yield nothing for the rest rather than raising.
        if spec["sheet"] not in pd.ExcelFile(blob).sheet_names:
            return pd.DataFrame()
        blob.seek(0)
        d = pd.read_excel(blob, sheet_name=spec["sheet"], header=1)
    acc = [c for c in d.columns if "Protein Accessions" in str(c)]
    fold = [c for c in d.columns if "Protein Fold Change" in str(c)]
    if len(acc) < 2 or len(fold) < 2:
        return pd.DataFrame()
    blocks = []
    for a, f in zip(acc[:2], fold[:2]):
        # A row can name several accessions when a peptide matched a protein group; the leading one
        # is the representative the search engine chose, and it is what the rest of this project
        # takes. Keeping the whole string would make `P08134; P61586` its own protein.
        block = pd.DataFrame({"host_id": d[a].astype(str).str.split(";").str[0].str.strip(),
                              "fold": pd.to_numeric(d[f], errors="coerce")}).dropna()
        blocks.append(block.groupby("host_id")["fold"].max())
    both = pd.concat(blocks, axis=1, join="inner")
    both = both[(both > REPLICATE_MIN_FOLD).all(axis=1)]
    if both.empty:
        return pd.DataFrame()
    # Parasite accessions belong to the parasite side of the experiment, not the host table.
    keep = ~both.index.str.contains("TGME49_|TGGT1_")
    both = both[keep]
    if both.empty:
        return pd.DataFrame()
    return pd.DataFrame({"host_id": both.index,
                         "host_name": both.index,
                         "host_ip_enrichment_log2": both.mean(axis=1).to_numpy()})


def bridges(base: str, spec: dict = None) -> pd.DataFrame:
    """Parasite gene to host protein, one row per pair, with what put it there."""
    spec = spec or MYR1_IP
    host = read_ip(base, spec)
    if host.empty:
        return pd.DataFrame()
    return pd.DataFrame({
        "gene_id": spec["bait"],
        "host_id": host["host_id"],
        "host_name": host["host_name"],
        "host_ip_enrichment_log2": host["host_ip_enrichment_log2"],
        "evidence": f"IP-MS {spec['accession']}"})


def all_bridges(base: str) -> pd.DataFrame:
    """Every bait in the project, as one bridge table.

    Four baits from two deposits. Keeping them in one table with the bait named per row is what lets
    a reader ask which parasite proteins reach a given host protein -- the question a single IP
    cannot answer and the reason a multi-bait bridge is worth more than the sum of its IPs.
    """
    parts = [bridges(base)]
    for spec in DIA_IPS:
        got = read_dia(base, spec)
        if got.empty:
            continue
        parts.append(pd.DataFrame({
            "gene_id": spec["bait"], "host_id": got["host_id"],
            "host_name": got["host_name"],
            "host_ip_enrichment_log2": got["host_ip_enrichment_log2"],
            "evidence": f"DIA affinity purification {spec['accession']} ({spec['bait_name']})"}))
    for spec in REPLICATE_IPS:
        got = read_replicated(base, spec)
        if got.empty:
            continue
        parts.append(pd.DataFrame({
            "gene_id": spec["bait"], "host_id": got["host_id"],
            "host_name": got["host_name"],
            "host_ip_enrichment_log2": got["host_ip_enrichment_log2"],
            "evidence": f"IP-MS {spec['accession']} ({spec['bait_name']})"}))
    parts = [p for p in parts if not p.empty]
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


#: A proximity experiment at the vacuole whose bait is the compartment rather than a named protein,
#: so it cannot be a bridge -- a bridge needs a parasite gene at one end. What it gives instead is a
#: property OF a host protein: how enriched it is at the vacuole, across three infection contexts.
PV_UPTAKE = {"accession": "PMC8700025",
             "archive": "PMC8700025/PMC8700025_supplementary.zip",
             "member": "ppat.1010138.s008.xlsx", "sheet": "Results Summary"}


def pv_enrichment(base: str) -> pd.DataFrame:
    """How enriched each host protein is at the parasitophorous vacuole, averaged over contexts.

    The sheet lists parasite and host proteins together, which is how the authors show the
    experiment worked -- the dense granule proteins top it. Only the host rows are kept here; the
    parasite rows are the positive control and belong to no host table.
    """
    import zipfile

    path = os.path.join(base, ESCRT_ROOT, PV_UPTAKE["archive"])
    if not os.path.exists(path):
        return pd.DataFrame()
    with zipfile.ZipFile(path) as archive:
        if PV_UPTAKE["member"] not in archive.namelist():
            return pd.DataFrame()
        import io as _io
        d = pd.read_excel(_io.BytesIO(archive.read(PV_UPTAKE["member"])),
                          sheet_name=PV_UPTAKE["sheet"], header=1)
    if d.shape[1] < 6:
        return pd.DataFrame()
    ident = d[d.columns[1]].astype(str).str.strip()
    host = ~ident.str.contains("TGME49_|TGGT1_")
    out = pd.DataFrame({
        "host_id": ident[host],
        "host_name": d[d.columns[0]].astype(str).str.strip()[host],
        "pv_enrichment_log2": pd.to_numeric(d[d.columns[5]], errors="coerce")[host]})
    out = out.dropna(subset=["host_id", "pv_enrichment_log2"])
    out = out[out["host_id"].str.fullmatch(
        r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")]
    return out.groupby("host_id", as_index=False).max()


def host_table(base: str) -> pd.DataFrame:
    """Every host protein the project knows, with whatever is measured about it.

    Identity from the bridges, and properties joined on. A host protein reached by a bait but never
    measured for vacuole enrichment keeps a NaN there, which is the difference between "not at the
    vacuole" and "this experiment did not see it".
    """
    bridge = all_bridges(base)
    if bridge.empty:
        return pd.DataFrame()
    out = bridge[["host_id", "host_name"]].drop_duplicates("host_id").set_index("host_id")
    pv = pv_enrichment(base)
    if not pv.empty:
        out = out.join(pv.set_index("host_id")[["pv_enrichment_log2"]], how="outer")
        out["host_name"] = out["host_name"].fillna(
            pv.set_index("host_id")["host_name"])
    return out.reset_index()


def load(base: str, name: str = HOST_TABLE) -> pd.DataFrame:
    """A shipped host table, or an empty frame if it has not been built."""
    path = os.path.join(base, "starplast", "data", name)
    return pd.read_parquet(path) if os.path.exists(path) else pd.DataFrame()
