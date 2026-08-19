#!/usr/bin/env python3
"""Relation layers for the Plasmodium table: pairs of genes, in their own index space.

The Toxoplasma graph stores edges as index pairs into `nodes.parquet`. A *falciparum* gene has no
index there, so the second species needs a second graph rather than more layers in the first -- the
same reason instruction 39 gives one table per species. Nothing is merged and nothing is compared
across the two files; an index means a row of the table it was built from and of no other.

## Three layers, built the way the Toxoplasma ones are

`orthogroup` and `domain` join genes that share a category. `coexpression` joins genes whose
abundance moves together across the life stages. The constructions are copied deliberately rather
than re-invented: if the two arms are ever compared, a difference should mean the biology differs
and not that the edges were drawn by different rules.

## The self-edge trap, inherited on purpose

InterPro reports one row per domain MATCH, so a protein carrying two copies of a repeat domain shows
up twice in a member list. Pairing a list that contains a gene twice emits (g, g) as an edge and
duplicates every other pair through it -- which is how 31 self-edges and 106 duplicated pairs once
shipped in the Toxoplasma graph. A self-edge draws as a zero-length line and a duplicate draws twice,
reading as twice the evidence. Members are collected into a SET here for that reason, and a test
asserts no edge joins a gene to itself.
"""
from __future__ import annotations

import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd

GRAPH = "pf_graph.npz"

#: Groups larger than this say only that a family is big, so they are dropped rather than drawn.
#: The same caps the Toxoplasma arm uses, for the same reason.
GROUP_CAP = 60

#: Correlation layer settings, also matched to the Toxoplasma arm.
CORRELATION = 0.95
NEIGHBOURS = 25

#: The life-stage columns co-transcription is computed across. Steady-state and polysomal are left
#: out: they are a three-point time course of one experiment, and mixing them with the seven-stage
#: series would let sampling depth decide which pairs correlate.
STAGE_COLUMNS = ("expr_ring", "expr_early_trophozoite", "expr_late_trophozoite", "expr_schizont",
                 "expr_gametocyte_ii", "expr_gametocyte_v", "expr_ookinete", "expr_asexual_blood",
                 "expr_oocyst", "expr_sporozoite")


def _pairs_within_groups(members: dict, cap: int = GROUP_CAP) -> tuple:
    """Every DISTINCT pair inside each group, weighted by how many groups the pair shares.

    Two proteins can share more than one InterPro domain, and emitting the pair once per shared
    domain draws the same edge several times -- which reads as several independent pieces of
    evidence. Plasmodium makes this unmissable: the naive construction produced 10,764 duplicate
    emissions among 34,887, almost all of them inside the var, rifin and stevor families whose
    members share whole multi-domain architectures. The count is kept as the weight instead, which
    is both correct and more informative than the 1.0 it replaces: two genes sharing four domains
    are more alike than two sharing one.

    The Toxoplasma arm builds this layer the same way and currently emits no duplicates at all --
    checked, not assumed -- so it has nothing to correct today. It would acquire the same fault the
    moment its annotation gained a pair sharing two domains.
    """
    counts = defaultdict(int)
    for _key, indices in members.items():
        ordered = sorted(indices)            # canonical order, so a < b for every pair
        if not 2 <= len(ordered) <= cap:
            continue
        for i in range(len(ordered)):
            for j in range(i + 1, len(ordered)):
                counts[(ordered[i], ordered[j])] += 1
    if not counts:
        return np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=float)
    pairs = sorted(counts)
    return (np.array([p[0] for p in pairs], dtype=int),
            np.array([p[1] for p in pairs], dtype=int),
            np.array([float(counts[p]) for p in pairs], dtype=float))


def orthogroup_edges(nodes: pd.DataFrame) -> tuple:
    """Genes sharing an OrthoMCL group."""
    if "orthogroup" not in nodes.columns:
        return np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=float)
    members = defaultdict(set)
    for index, key in enumerate(nodes["orthogroup"].fillna("").astype(str)):
        if key and key not in ("nan", "unassigned", "N/A"):
            members[key].add(index)
    return _pairs_within_groups(members)


def domain_edges(nodes: pd.DataFrame) -> tuple:
    """Genes sharing an InterPro domain. Members are a SET; see the module docstring."""
    if "interpro_ids" not in nodes.columns:
        return np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=float)
    members = defaultdict(set)
    for index, cell in enumerate(nodes["interpro_ids"].fillna("").astype(str)):
        for domain in cell.split(";"):
            domain = domain.strip()
            if domain and domain != "nan":
                members[domain].add(index)
    return _pairs_within_groups(members)


def correlation_edges(nodes: pd.DataFrame, columns=STAGE_COLUMNS,
                      threshold: float = CORRELATION, neighbours: int = NEIGHBOURS) -> tuple:
    """Genes whose abundance moves together, keeping each gene's strongest neighbours."""
    present = [c for c in columns if c in nodes.columns]
    if len(present) < 3:
        return np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=float)
    matrix = nodes[present].to_numpy(dtype=float)
    usable = ~np.isnan(matrix).any(axis=1)
    if usable.sum() < 50:
        return np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=float)
    values = matrix[usable]
    values = (values - values.mean(0)) / (values.std(0) + 1e-9)
    correlation = np.corrcoef(values)
    np.fill_diagonal(correlation, 0.0)
    rows = np.where(usable)[0]
    a, b, w = [], [], []
    for i in range(len(rows)):
        for j in np.argsort(-correlation[i])[:neighbours]:
            if correlation[i, j] < threshold or rows[i] >= rows[j]:
                continue
            a.append(rows[i])
            b.append(rows[j])
            w.append(float(correlation[i, j]))
    return np.array(a, dtype=int), np.array(b, dtype=int), np.array(w, dtype=float)


def build(nodes: pd.DataFrame, log=print, dataset_root: str | None = None) -> dict:
    """Every layer, in the `<name>__a` / `<name>__b` / `<name>__w` shape the graph file uses."""
    out = {}
    if nodes.empty:
        return out
    for label, (a, b, w) in (("orthogroup", orthogroup_edges(nodes)),
                             ("domain", domain_edges(nodes))):
        if len(a):
            out[f"{label}__a"], out[f"{label}__b"], out[f"{label}__w"] = a, b, w
            shared = f", up to {int(w.max())} shared" if w.max() > 1 else ""
            log(f"{label}: {len(a):,} edges{shared}")
    a, b, w = correlation_edges(nodes)
    if len(a):
        out["coexpression__a"], out["coexpression__b"], out["coexpression__w"] = a, b, w
        log(f"coexpression: {len(a):,} edges (top-{NEIGHBOURS} neighbours, r >= {CORRELATION})")
    if dataset_root:
        a, b, w = crosslink_edges(nodes, dataset_root, log=log)
        if len(a):
            out["xlms__a"], out["xlms__b"], out["xlms__w"] = a, b, w
        a, b, w = ip_ms_edges(nodes, dataset_root, log=log)
        if len(a):
            out["ip_ms__a"], out["ip_ms__b"], out["ip_ms__w"] = a, b, w
    return out


def save(nodes: pd.DataFrame, path: str, log=print, dataset_root: str | None = None) -> dict:
    """Build the layers AND the 3D layout, and write them together.

    The layout was missing until 2026-08-17, and its absence is why the application could not open
    this arm at all: `app.load` reads `xyz` out of the graph file, and this file had every edge layer
    and no coordinates. The arm was 41 filled slots that nothing could draw.

    Laid out by the same `build_graph.embed` the Toxoplasma table uses, so the two maps are made the
    same way -- not the same axes, which would be meaningless across species with different features,
    but the same construction, which is what makes comparing them honest.
    """
    layers = build(nodes, log=log, dataset_root=dataset_root)
    if not layers:
        return layers
    from .build_graph import embed
    layers["xyz"] = embed(nodes)
    log(f"layout: {len(layers['xyz']):,} genes placed in 3D")
    np.savez_compressed(path, **layers)
    return layers


def load(base: str):
    """The shipped Plasmodium graph, or None if it has not been built."""
    path = os.path.join(base, "starplast", "data", GRAPH)
    return np.load(path, allow_pickle=True) if os.path.exists(path) else None


# --------------------------------------------------------------------------- measured contacts
#: Crosslinking mass spectrometry of the infected erythrocyte. The layer is named `xlms` to match the
#: Toxoplasma one, which answers the same slot.
CROSSLINK = ("crosslink", "41966402", "mmc1.xlsx")
CROSSLINK_SHEET = "(D) Protein pairs"
ACCESSION = re.compile(r"PF3D7_\w+")


def _classify(cell, known: set, owners: dict) -> tuple:
    """(gene, is_parasite) for one protein in a crosslink pair.

    THREE outcomes, not two, and conflating the last two is a bug this cost a test to find. A protein
    is a known parasite gene, or a parasite protein this table does not carry -- a deprecated
    accession, a gene dropped from the current annotation -- or genuinely host. Treating the middle
    case as host inflates host degree and would put a parasite protein into a host bridge; treating it
    as parasite would index it against a row that does not exist. Either way the pair is unusable, so
    it is reported as parasite-with-no-gene and the caller skips it.
    """
    found = ACCESSION.search(str(cell))
    if found:
        return (found.group(0) if found.group(0) in known else None), True
    parts = str(cell).split("|")
    gene = owners.get(parts[1]) if len(parts) > 2 else None
    if gene:
        return (gene if gene in known else None), True
    return None, False


def crosslink_edges(nodes: pd.DataFrame, dataset_root: str, log=print) -> tuple:
    r"""Protein pairs joined by a measured crosslink, weighted by how many were found.

    Two filters matter here and neither is optional. The experiment crosslinked PARASITE INSIDE
    ERYTHROCYTE, so a third of the pairs have a human protein at one or both ends -- spectrin, band 3,
    protein 4.2. Those are real contacts and they are a host bridge, not a parasite-parasite edge, so
    they are dropped from this layer rather than quietly indexed against a table that has no row for
    them. And a pair whose two ends resolve to the same gene is a homomeric crosslink: evidence the
    protein self-associates, not an edge between two genes, and drawing it would put a zero-length
    line in the graph.

    Identifying which end is which must go through the IDENTITY INDEX and not through a pattern. The
    first version of this matched `PF3D7_\w+` in the mapping field and shipped 73 edges; the source
    writes some rows with a UniProt symbol instead -- `sp|Q6ZMA7|Pfs16` is Pfs16, which is
    PF3D7_0406200 and a parasite gene -- so six real parasite-parasite contacts were dropped as
    "host at one end", and one parasite protein was on its way into a host bridge. Resolving through
    the accession index recovers all 79. This is the mistake the Toxoplasma identity layer exists to
    prevent, made here anyway because a regex on an accession looks like resolution and is not.
    """
    folder, pmid, name = CROSSLINK
    path = os.path.join(dataset_root, "reference", "plasmodb", folder, pmid, name)
    empty = (np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=float))
    if not os.path.exists(path) or nodes.empty:
        return empty
    book = pd.ExcelFile(path)
    if CROSSLINK_SHEET not in book.sheet_names:
        return empty
    d = book.parse(CROSSLINK_SHEET)
    if not {"Protein1", "Protein2"} <= set(d.columns):
        return empty
    from .plasmodium import UNIPROT_TABLE, uniprot_index
    index = {gene: i for i, gene in enumerate(nodes["gene_id"].astype(str))}
    owners = uniprot_index(os.path.join(dataset_root, "reference", "plasmodb", UNIPROT_TABLE))

    known = set(index)

    def resolve(cell):
        """The row this protein is, or None if it is host or a parasite gene the table lacks."""
        gene, _parasite = _classify(cell, known, owners)
        return index.get(gene) if gene else None

    counts = {}
    for one, two, links in zip(d["Protein1"], d["Protein2"],
                               pd.to_numeric(d.get("Num_Crosslinks", 1), errors="coerce")):
        a, b = resolve(one), resolve(two)
        if a is None or b is None or a == b:
            continue
        pair = (min(a, b), max(a, b))
        counts[pair] = max(counts.get(pair, 0.0), float(links) if links == links else 1.0)
    if not counts:
        return empty
    pairs = sorted(counts)
    log(f"xlms: {len(pairs)} parasite-parasite pairs of {len(d)} crosslinked protein pairs "
        f"(the rest have a human protein at one end, or are homomeric)")
    return (np.array([p[0] for p in pairs], dtype=int),
            np.array([p[1] for p in pairs], dtype=int),
            np.array([counts[p] for p in pairs], dtype=float))


def ip_ms_edges(nodes: pd.DataFrame, dataset_root: str, log=print) -> tuple:
    """Bait-to-partner edges from the co-immunoprecipitations, weighted by spectra in the bait.

    Named `ip_ms` to match the Toxoplasma layer answering the same slot -- and, as everywhere on this
    arm, the shared name is why the species guard exists: these indices are positions in
    `pf_nodes.parquet` and mean different genes in the other graph.

    A bait's own row is not an edge (a protein does not bind itself into the graph), and the pairs
    are undirected here even though the experiment is not: `A pulled down B` and `B pulled down A`
    are one line between two genes, and the direction is kept in the table the loader returns rather
    than in the layer.
    """
    from .plasmodium import ip_ms
    empty = (np.array([], dtype=int),) * 2 + (np.array([], dtype=float),)
    pairs_table = ip_ms(dataset_root, log=log)
    if pairs_table.empty or nodes.empty:
        return empty
    index = {gene: i for i, gene in enumerate(nodes["gene_id"].astype(str))}
    weight = {}
    for bait, prey, spectra in zip(pairs_table["bait"], pairs_table["prey"],
                                   pairs_table["spectra_bait"]):
        a, b = index.get(bait), index.get(prey)
        if a is None or b is None or a == b:
            continue
        pair = (min(a, b), max(a, b))
        weight[pair] = max(weight.get(pair, 0.0), float(spectra))
    if not weight:
        return empty
    pairs = sorted(weight)
    log(f"ip_ms: {len(pairs)} bait-partner edges over {len(set(pairs_table['bait']))} pulldowns")
    return (np.array([p[0] for p in pairs], dtype=int),
            np.array([p[1] for p in pairs], dtype=int),
            np.array([weight[p] for p in pairs], dtype=float))


def host_bridge(nodes: pd.DataFrame, dataset_root: str, log=print) -> pd.DataFrame:
    """The crosslinks this layer throws away: parasite protein to HUMAN protein.

    The same file the `xlms` layer reads, and the other half of its filter. A crosslink between a
    parasite protein and an erythrocyte one is a measured contact across the host boundary, which is
    a bridge rather than an edge -- the pair's two ends live in different tables and a human protein
    has no index in this one, exactly as instruction 39 describes.

    Which end is which goes through the accession index, for the reason `crosslink_edges` records:
    reading it off a pattern classified Pfs16, a parasite protein, as human.
    """
    from .plasmodium import UNIPROT_TABLE, uniprot_index
    folder, pmid, name = CROSSLINK
    path = os.path.join(dataset_root, "reference", "plasmodb", folder, pmid, name)
    if not os.path.exists(path) or nodes.empty:
        return pd.DataFrame()
    book = pd.ExcelFile(path)
    if CROSSLINK_SHEET not in book.sheet_names:
        return pd.DataFrame()
    d = book.parse(CROSSLINK_SHEET)
    if not {"Protein1", "Protein2"} <= set(d.columns):
        return pd.DataFrame()
    known = set(nodes["gene_id"].astype(str))
    owners = uniprot_index(os.path.join(dataset_root, "reference", "plasmodb", UNIPROT_TABLE))

    rows = []
    for one, two, links in zip(d["Protein1"], d["Protein2"],
                               pd.to_numeric(d.get("Num_Crosslinks", 1), errors="coerce")):
        (a, a_par), (b, b_par) = _classify(one, known, owners), _classify(two, known, owners)
        # A bridge needs one known parasite gene and one protein that is definitely host. A parasite
        # protein the table does not carry is neither, and its pair is skipped.
        if not ((a and not b_par) or (b and not a_par)):
            continue
        gene, other = (a, two) if a else (b, one)
        parts = str(other).split("|")
        if len(parts) < 3:
            continue
        rows.append({"gene_id": gene, "host_id": parts[1], "host_name": parts[2],
                     "crosslinks": float(links) if links == links else 1.0,
                     "evidence": f"crosslink MS {CROSSLINK[1]}", "bridge": "host"})
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows).drop_duplicates(["gene_id", "host_id"]).reset_index(drop=True)
    log(f"host bridge (crosslink MS): {len(out)} pairs, {out.gene_id.nunique()} parasite genes, "
        f"{out.host_id.nunique()} human proteins")
    return out


def host_degree(nodes: pd.DataFrame, dataset_root: str, log=print) -> pd.Series:
    """How many host proteins each parasite gene was crosslinked to.

    Three states, and the middle one is the reason this is not a `fillna(0)`. The Toxoplasma arm sets
    this to 0 for every gene without a curated host target, which is right there because its source is
    a curated table covering the whole literature. This source is ONE experiment: a gene it never
    detected has not been shown to lack host partners, it was simply not seen. So genes that appear in
    the crosslink data get a count -- zero included, because being crosslinked only to parasite
    proteins is a real observation -- and genes absent from it stay missing.
    """
    folder, pmid, name = CROSSLINK
    path = os.path.join(dataset_root, "reference", "plasmodb", folder, pmid, name)
    if not os.path.exists(path) or nodes.empty:
        return pd.Series(dtype=float)
    book = pd.ExcelFile(path)
    if CROSSLINK_SHEET not in book.sheet_names:
        return pd.Series(dtype=float)
    d = book.parse(CROSSLINK_SHEET)
    if not {"Protein1", "Protein2"} <= set(d.columns):
        return pd.Series(dtype=float)
    from .plasmodium import UNIPROT_TABLE, uniprot_index
    known = set(nodes["gene_id"].astype(str))
    owners = uniprot_index(os.path.join(dataset_root, "reference", "plasmodb", UNIPROT_TABLE))

    seen, partners = set(), {}
    for one, two in zip(d["Protein1"], d["Protein2"]):
        (a, a_par), (b, b_par) = _classify(one, known, owners), _classify(two, known, owners)
        for gene in (a, b):
            if gene:
                seen.add(gene)
                partners.setdefault(gene, set())
        # Only a definite host protein counts as a partner. A parasite gene the table lacks does not.
        if a and not b_par:
            gene, other = a, two
        elif b and not a_par:
            gene, other = b, one
        else:
            continue
        parts = str(other).split("|")
        if len(parts) > 2:
            partners[gene].add(parts[1])
    if not seen:
        return pd.Series(dtype=float)
    counts = {gene: float(len(partners.get(gene, ()))) for gene in seen}
    out = nodes["gene_id"].astype(str).map(counts)
    log(f"host degree: {len(seen)} genes seen in the crosslink data, "
        f"{sum(1 for v in counts.values() if v)} with a host partner; the rest of the proteome "
        f"stays missing rather than zero")
    return out
