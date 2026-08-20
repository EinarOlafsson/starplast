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

One row per (parasite gene, host protein) pair, with the evidence that put it there and the number
that evidence produced. Different evidence produces different NUMBERS -- a log2 enrichment from a
fold change, a probability from SAINT -- and they occupy different columns rather than being poured
into one. A SAINT probability of 1.0 beside a log2 of 4.11 in the same column would read as the same
quantity twice. It is NOT an
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

#: The same two affinity purifications scored the way the AP-MS field scores them -- SAINT, MiST and
#: CompPASS rather than a fold change and a p-value. SAINT gives a probability that an interaction is
#: real given the controls, which is a better instrument than the `logFC > 1` used above, so where
#: both exist these win. UNPUBLISHED: manuscript supplementary tables.
ORTHOGONAL_IPS = (
    {"bait": "TGME49_225160", "bait_name": "EAF1", "file": "orthogonal/eaf1_orthogonal_scores.csv"},
    {"bait": "TGME49_239740", "bait_name": "GRA14", "file": "orthogonal/gra14_orthogonal_scores.csv"},
)
#: SAINT's conventional cut. Below it an interaction is not called, and calling one anyway would be
#: substituting a threshold of mine for the method's own.
SAINT_MIN = 0.9


def _symbol_to_accession(base: str) -> dict:
    """Gene symbol to UniProt accession, from the deposit that carries both.

    The orthogonal tables name preys by symbol and the bridge is keyed on accession, so the mapping
    has to come from somewhere. It comes from the same experiment's own protein report rather than
    from an external service, which keeps the join inside the data.
    """
    path = os.path.join(base, DIA_FILE)
    if not os.path.exists(path):
        return {}
    d = pd.read_csv(path)
    if "gene_symbol" not in d.columns or "uniprot_id" not in d.columns:
        return {}
    symbol = d["gene_symbol"].astype(str).str.replace(r'^="?|"?$', "", regex=True)
    return dict(zip(symbol, d["uniprot_id"].astype(str)))


def read_orthogonal(base: str, spec: dict, symbols: dict = None) -> pd.DataFrame:
    """One bait scored by SAINT, keeping preys the method itself calls."""
    path = os.path.join(base, ESCRT_ROOT, spec["file"])
    if not os.path.exists(path):
        return pd.DataFrame()
    d = pd.read_csv(path)
    if not {"genesymbol", "SAINT_AvgP"} <= set(d.columns):
        return pd.DataFrame()
    symbols = _symbol_to_accession(base) if symbols is None else symbols
    name = d["genesymbol"].astype(str).str.replace(r'^="?|"?$', "", regex=True)
    saint = pd.to_numeric(d["SAINT_AvgP"], errors="coerce")
    keep = saint >= SAINT_MIN
    # SAINT is a PROBABILITY, so it goes in its own column. Putting it in the log2 column would
    # make a 1.0 sit next to a 4.11 as though they were the same quantity, which is the mistake the
    # rest of this project spends its docstrings avoiding.
    out = pd.DataFrame({"host_name": name[keep],
                        "host_id": name[keep].map(symbols),
                        "saint_avgp": saint[keep]})
    # A prey with no accession is the parasite side of the experiment -- the bait's own partners
    # among Toxoplasma proteins, which belong to no host table.
    out = out.dropna(subset=["host_id"])
    return out.groupby("host_id", as_index=False).max()
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
    symbols = _symbol_to_accession(base)
    for spec in ORTHOGONAL_IPS:
        got = read_orthogonal(base, spec, symbols)
        if got.empty:
            continue
        parts.append(pd.DataFrame({
            "gene_id": spec["bait"], "host_id": got["host_id"], "host_name": got["host_name"],
            "saint_avgp": got["saint_avgp"],
            "evidence": f"AP-MS SAINT (unpublished, {spec['bait_name']})"}))
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


# --------------------------------------------------------------------------- host tissue proteomes
#: The red blood cell the blood stage lives in, measured as two fractions of one preparation.
#: `(folder, pmid, file)`; the sheets are named for the fractions.
ERYTHROCYTE = ("erythrocyte", "41654503", "41597_2026_6792_MOESM2_ESM.xlsx")
ERYTHROCYTE_FRACTIONS = {"Membrane extract": "rbc_membrane_psms",
                         "Cytoplasmic extract": "rbc_cytoplasm_psms"}


def erythrocyte_proteome(dataset_root: str, log=print) -> pd.DataFrame:
    """Human red blood cell proteins, by fraction, keyed the way the host table is keyed.

    A host tissue reference rather than a bridge: these are the proteins present in the cell the
    parasite lives in, which is the question instruction 39's host slots ask and which no pulldown
    can answer -- a pulldown says what a bait touched, not what is there to touch.

    The two fractions are kept apart because they are different measurements: a protein in the
    membrane extract is at the surface the parasite invades through, one in the cytoplasm is in the
    haemoglobin soup around it. Self-validating in the way a fractionation should be -- spectrin
    beta heads the membrane list, haemoglobin alpha heads the cytoplasmic one.
    """
    folder, pmid, name = ERYTHROCYTE
    path = os.path.join(dataset_root, "host", folder, pmid, name)
    if not os.path.exists(path):
        return pd.DataFrame()
    book = pd.ExcelFile(path)
    frames = []
    for sheet, column in ERYTHROCYTE_FRACTIONS.items():
        if sheet not in book.sheet_names:
            continue
        d = book.parse(sheet)
        if not {"Accession", "Gene", "# PSMs"} <= set(d.columns):
            continue
        part = pd.DataFrame({
            "host_id": d["Accession"].astype(str).str.strip(),
            # A row can name several genes (`HBA1; HBA2`); the host table keys on the accession and
            # keeps the string as given rather than choosing one of them.
            "host_name": d["Gene"].astype(str).str.strip(),
            column: pd.to_numeric(d["# PSMs"], errors="coerce")})
        frames.append(part[part["host_id"].str.match(r"^[A-Z0-9]{6,10}$", na=False)])
    if not frames:
        return pd.DataFrame()
    out = frames[0]
    for part in frames[1:]:
        out = out.merge(part, on=["host_id", "host_name"], how="outer")
    out = out.drop_duplicates("host_id").reset_index(drop=True)
    found = {c: int(out[c].notna().sum()) for c in ERYTHROCYTE_FRACTIONS.values() if c in out}
    log(f"erythrocyte proteome: {len(out):,} human proteins ({found})")
    return out


#: The Cell Surface Protein Atlas, `(folder, pmid, file)`. Cell-surface capture on 41 human and 31
#: mouse cell types; only one of them is a tissue this project's parasites live in, and it is the
#: primary mouse bone-marrow-derived macrophage the tachyzoite is grown in.
CSPA = ("cspa", "25894527", "S1_File.xlsx")
#: The sheets, and the column each one calls the macrophage. They are named differently in the two
#: matrices, which is why both are written out rather than derived from the tissue name.
CSPA_SHEETS = {"annotation": "Table_A", "mouse_matrix": "Table_C", "mouse_intensity": "Table_F"}
CSPA_BMDM = {"mouse_matrix": "Macrophages (BM-der.)", "mouse_intensity": "BMd_Macrophages"}


def surface_repertoire(dataset_root: str, log=print) -> pd.DataFrame:
    """The surface proteins of a primary mouse bone-marrow macrophage, with their absences.

    A repertoire rather than a proteome, and the difference is the point: cell-surface capture
    labels what is exposed on an intact cell, so a protein here is at the surface the tachyzoite
    meets rather than merely present somewhere in the cell. The row space is every protein the atlas
    saw on the surface of ANY mouse cell type, so a `False` is a real negative -- the same capture
    ran on macrophages and did not find it -- which is what instruction 39 means by a repertoire
    slot being FILLED rather than averaged.

    The deposit's two matrices disagree for twelve proteins: nine carry a macrophage intensity in
    the abundance sheet with no mark in the detection matrix, and three are marked without one. A
    protein either sheet places on the macrophage is counted present, because both are the authors'
    own record of having measured it there, and the disagreement is logged rather than smoothed.

    Self-validating: the strongest signals are Emr1 (F4/80), Siglec1 (CD169), Itgb2, Cd47 and
    H2-K1 -- the surface a macrophage is identified by.
    """
    folder, pmid, name = CSPA
    path = os.path.join(dataset_root, "host", folder, pmid, name)
    if not os.path.exists(path):
        return pd.DataFrame()
    book = pd.ExcelFile(path)
    if not set(CSPA_SHEETS.values()) <= set(book.sheet_names):
        return pd.DataFrame()

    def _columns(frame):
        frame.columns = [str(c).strip() for c in frame.columns]
        return frame

    matrix = _columns(book.parse(CSPA_SHEETS["mouse_matrix"]))
    intensity = _columns(book.parse(CSPA_SHEETS["mouse_intensity"]))
    annotation = _columns(book.parse(CSPA_SHEETS["annotation"]))
    if CSPA_BMDM["mouse_matrix"] not in matrix.columns or \
            CSPA_BMDM["mouse_intensity"] not in intensity.columns:
        return pd.DataFrame()

    # One row per accession. The annotation sheet lists a protein once per peptide, so it is
    # collapsed before it can multiply the matrix it is joined to.
    mouse = annotation[annotation["organism"].astype(str).str.strip().str.lower() == "mouse"]
    symbols = (mouse.drop_duplicates("ID_link").set_index("ID_link")["ENTREZ gene symbol"]
               .astype(str).str.strip())
    values = (intensity.set_index("Protein")[CSPA_BMDM["mouse_intensity"]]
              .groupby(level=0).max().dropna())
    ids = matrix["ID_link"].astype(str).str.strip()
    out = pd.DataFrame({
        "host_id": ids,
        "host_name": [symbols.get(i, "") for i in ids],
        "bmdm_surface_intensity": [values.get(i, np.nan) for i in ids]})
    out["marked"] = (matrix[CSPA_BMDM["mouse_matrix"]] == 1).to_numpy()
    # Nullable boolean on purpose. This column joins a table that also holds human red cell rows,
    # where the question was never asked, and a plain numpy bool would turn those into `False` --
    # "cell-surface capture looked and did not find it" said about a cell it never touched. The
    # dtype also survives the parquet round trip as a boolean, which is what makes the atlas count
    # its TRUEs rather than its rows.
    out["bmdm_surface_detected"] = (out["marked"]
                                    | out["bmdm_surface_intensity"].notna()).astype("boolean")
    out = out[out["host_id"].str.fullmatch(
        r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")]
    out = out.drop_duplicates("host_id").reset_index(drop=True)
    disagree = int((out["marked"] != out["bmdm_surface_intensity"].notna()).sum())
    out = out.drop(columns=["marked"])
    quantified = int(out["bmdm_surface_intensity"].notna().sum())
    present = int(out["bmdm_surface_detected"].sum())
    log(f"macrophage surfaceome: {present} of {len(out):,} mouse surface proteins on BMDM "
        f"({quantified} with an intensity, {disagree} where the two sheets disagree)")
    return out


def merge_tissue(existing: pd.DataFrame, new: pd.DataFrame) -> pd.DataFrame:
    """One host table, several owners, and a shared key rather than a shared column.

    `build_graph` writes the Toxoplasma pulldown's columns into the same file a tissue reference
    writes into, and a tissue of one host species has no accession in common with a tissue of
    another. Joined on the key so both survive: a column the new frame does not carry keeps its
    old values, and a row it does not mention keeps existing.
    """
    if new.empty:
        return existing
    if not len(existing):
        return new
    # `host_id` and `host_name` are the key and its label rather than measurements, so they are
    # kept from BOTH sides and reconciled below. Dropping the old `host_name` because the new frame
    # also carries one is how 311 bridge proteins lost their symbols in the shipped table.
    keep = [c for c in existing.columns
            if c not in new.columns or c in ("host_id", "host_name")]
    merged = existing[keep].merge(new, on="host_id", how="outer")
    # `host_name` comes from both sides; the one already there wins, since it is what the bridges
    # were written against.
    if "host_name_x" in merged.columns:
        merged["host_name"] = merged["host_name_x"].fillna(merged["host_name_y"])
        merged = merged.drop(columns=["host_name_x", "host_name_y"])
    return merged


#: The red cell SURFACE, `(folder, pmid, file)`, and the sheet holding all 267 plasma-membrane
#: proteins. A different question from the fractionation above: plasma membrane profiling asks what
#: is on the outside of an intact cell, and the merozoite's receptors are exactly there.
RBC_SURFACE = ("erythrocyte", "31552303", "42003_2019_596_MOESM6_ESM.xlsx")
RBC_SURFACE_SHEET = "Data S2A"
#: Column in the sheet -> column shipped. The two donor populations are kept APART: the paper's
#: finding is that they differ, and averaging a Duffy-positive population with a Duffy-negative one
#: would erase the single best-known receptor polymorphism in malaria.
RBC_SURFACE_COLUMNS = {"Estimated copies/cell UK sample": "rbc_surface_copies_uk",
                       "Estimated copies/cell Senegal sample": "rbc_surface_copies_senegal",
                       "Found in UK sample": "rbc_surface_found_uk",
                       "Found in Senegal sample": "rbc_surface_found_senegal"}


def surface_receptors(dataset_root: str, log=print) -> pd.DataFrame:
    """The red blood cell surface, per donor population, in copies per cell.

    The invasion interface. `erythrocyte_proteome` says what is in the cell; this says what is
    reachable from outside it, which is the question a receptor slot asks and the reason the study
    was done -- the authors' stated aim is candidate Plasmodium receptors.

    Both populations ship as their own columns because the difference IS the result. A zero is the
    paper's own encoding for "not identified in this population's donors", which is why the found
    flags ship beside the counts: for ACKR1 the zero is the Duffy-negative phenotype of West Africa
    rather than a detection failure, and nothing but the flag lets a reader tell those apart.

    Self-validating against numbers measured long before mass spectrometry: band 3 comes out at
    1.3 million copies per cell and glycophorin A at 3.3 million, the two most abundant proteins of
    the red cell membrane; basigin, the receptor PfRH5 must bind, is present in both populations;
    CR1 and the Duffy antigen are far lower in the Senegalese donors, which is what the population
    genetics of malaria says they should be.
    """
    folder, pmid, name = RBC_SURFACE
    path = os.path.join(dataset_root, "host", folder, pmid, name)
    if not os.path.exists(path):
        return pd.DataFrame()
    book = pd.ExcelFile(path)
    if RBC_SURFACE_SHEET not in book.sheet_names:
        return pd.DataFrame()
    d = book.parse(RBC_SURFACE_SHEET)
    d.columns = [str(c).strip() for c in d.columns]
    if not {"Uniprot", "Name"} | set(RBC_SURFACE_COLUMNS) <= set(d.columns):
        return pd.DataFrame()
    out = pd.DataFrame({"host_id": d["Uniprot"].astype(str).str.strip(),
                        "host_name": d["Name"].astype(str).str.strip()})
    for source, column in RBC_SURFACE_COLUMNS.items():
        values = pd.to_numeric(d[source], errors="coerce")
        out[column] = values.astype("boolean") if column.startswith("rbc_surface_found") else values
    out = out[out["host_id"].str.fullmatch(
        r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")]
    out = out.drop_duplicates("host_id").reset_index(drop=True)
    both = int((out["rbc_surface_found_uk"].fillna(False)
                & out["rbc_surface_found_senegal"].fillna(False)).sum())
    log(f"red cell surface: {len(out):,} plasma-membrane proteins, {both} in both populations "
        f"({int(out['rbc_surface_found_uk'].fillna(False).sum())} UK, "
        f"{int(out['rbc_surface_found_senegal'].fillna(False).sum())} Senegal)")
    return out


#: Every tissue reference the project has loaded, newest last. A tissue is one entry, so adding one
#: is a line here plus its loader rather than an edit to a build script.
TISSUE_REFERENCES = ("erythrocyte_proteome", "surface_receptors", "surface_repertoire")


def tissue_references(dataset_root: str, log=print) -> pd.DataFrame:
    """The host tissues, merged onto one key. Empty if none of them are on disk."""
    out = pd.DataFrame()
    for name in TISSUE_REFERENCES:
        out = merge_tissue(out, globals()[name](dataset_root, log=log))
    return out
