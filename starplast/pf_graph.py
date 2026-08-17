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
    return out


def save(nodes: pd.DataFrame, path: str, log=print, dataset_root: str | None = None) -> dict:
    """Build the layers and write them, or write nothing if there are none to write."""
    layers = build(nodes, log=log, dataset_root=dataset_root)
    if layers:
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


def crosslink_edges(nodes: pd.DataFrame, dataset_root: str, log=print) -> tuple:
    """Protein pairs joined by a measured crosslink, weighted by how many were found.

    Two filters matter here and neither is optional. The experiment crosslinked PARASITE INSIDE
    ERYTHROCYTE, so a third of the pairs have a human protein at one or both ends -- spectrin, band 3,
    protein 4.2. Those are real contacts and they are a host bridge, not a parasite-parasite edge, so
    they are dropped from this layer rather than quietly indexed against a table that has no row for
    them. And a pair whose two ends resolve to the same gene is a homomeric crosslink: evidence the
    protein self-associates, not an edge between two genes, and drawing it would put a zero-length
    line in the graph.
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
    index = {gene: i for i, gene in enumerate(nodes["gene_id"].astype(str))}

    def resolve(cell):
        found = ACCESSION.search(str(cell))
        return index.get(found.group(0)) if found else None

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
