#!/usr/bin/env python3
"""Build the Toxoplasma knowledge graph and its 3D embedding, once, into a cache the app loads instantly.

Design notes that matter (see ../HANDOFF.md):

* Position is a UMAP embedding of a multimodal feature matrix, NOT a force-directed layout. Force-directed
  position is aesthetic and arbitrary; UMAP over expression, the seven fitness screens, compartment,
  orthology breadth, domain content and structure confidence gives positions where proximity means
  biological similarity.
* Six edge types are stored separately and never merged. Each answers a different question and they are not
  interchangeable -- the seven CRISPR screens in particular disagree with each other substantially.
* Co-mention edges are attention-biased. We therefore store BOTH the raw co-mention count and an
  attention-corrected residual (observed minus expected given each gene's own publication count). The app
  defaults to the corrected form because the raw form is confidently misleading.

Output: data/graph.npz + data/nodes.parquet
"""
from __future__ import annotations

import itertools
import os
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from . import corpus, identity, literature

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.dirname(HERE)                      # toxoplasma_projects
DS = os.path.join(BASE, "datasets")
OUT = os.path.join(HERE, "data")
os.makedirs(OUT, exist_ok=True)

FIT = ["fit_invitro_hff", "fit_invivo_PE", "fit_invivo_lung", "fit_invivo_liver",
       "fit_invivo_spleen", "fit_naive_bmdm", "fit_ifng"]


def log(m):
    print(f"[build_graph] {m}", flush=True)


# --------------------------------------------------------------------------- nodes
def load_nodes() -> pd.DataFrame:
    n = pd.read_parquet(os.path.join(BASE, "toxonet", "data", "interim", "nodes.parquet"))
    n = n.drop_duplicates("gene_id").reset_index(drop=True)
    log(f"{len(n)} genes from the node table")

    lp = pd.read_csv(os.path.join(DS, "lopit_toxoplasma_gondii_ME49.csv"), low_memory=False)
    lp = lp[["gene_source_id", "MAP_location"]].rename(
        columns={"gene_source_id": "gene_id", "MAP_location": "compartment"})
    n = n.merge(lp.drop_duplicates("gene_id"), on="gene_id", how="left")

    prod = pd.read_csv(os.path.join(DS, "orthomcl_toxoplasma_gondii_ME49.csv"), low_memory=False)
    prod = prod[["gene_source_id", "gene_product"]].rename(columns={"gene_source_id": "gene_id"})
    n = n.merge(prod.drop_duplicates("gene_id"), on="gene_id", how="left", suffixes=("", "_p"))
    if "product" not in n.columns:
        n["product"] = n.get("gene_product")
    n["product"] = n["product"].fillna(n.get("gene_product")).fillna("unannotated")

    for c in FIT + ["mean_plddt", "paralog_number", "n_interpro", "n_phosphosites"]:
        if c in n.columns:
            n[c] = pd.to_numeric(n[c], errors="coerce")
        else:
            n[c] = np.nan

    rna = [c for c in n.columns if c.startswith("rna108740_")]
    for c in rna:
        n[c] = pd.to_numeric(n[c], errors="coerce")
    tach = [c for c in rna if "Tachyzoite" in c or "Tachyzoites" in c]
    cyst = [c for c in rna if "Tissue_cyst" in c or "Tissue_cysts" in c]
    n["expr_tachy"] = np.log2(n[tach].mean(axis=1) + 1) if tach else np.nan
    n["expr_cyst"] = np.log2(n[cyst].mean(axis=1) + 1) if cyst else np.nan
    n["expr_max"] = np.log2(n[rna].max(axis=1) + 1) if rna else np.nan

    # hyperLOPIT assignment tracks abundance, so unassigned must read as UNKNOWN, never as a compartment
    n["compartment"] = n["compartment"].fillna("unassigned")
    n["has_domain"] = n.get("has_domain", pd.Series(False, index=n.index)) \
        .astype("boolean").fillna(False).astype(int)
    n["lineage_specific"] = n.get("lineage_specific", pd.Series(False, index=n.index)) \
        .astype("boolean").fillna(False).astype(int)
    return n


# --------------------------------------------------------------------------- literature
ABSTRACTS = os.path.join(BASE, ".claude", "skills", "toxoplasma-scientist", "corpus",
                         "pubmed_toxoplasma.jsonl")
# Open-access full texts. Machine-local: on a machine without this disk the build silently falls back to
# abstracts only, which is why the committed cache is what the app actually ships.
FULLTEXTS = "/mnt/wd4tb/skill_corpora/toxoplasma-scientist"


def literature_layer(nodes: pd.DataFrame):
    """Scan both literature sources through the identity layer; return co-mention edges and node columns.

    This replaces an earlier one-pass scan that matched only current ME49 accessions and ToxoDB symbols.
    That version could not see the old TGME49_0xxxxx accessions, the GT1/VEG accessions papers use
    interchangeably, or the Tg- prefixed symbol forms, and it wrote no intermediate anyone could audit.
    See identity.py, corpus.py and literature.py for the three layers this composes.

    Abstracts and full texts stay separate all the way through: they are different populations (all of
    PubMed versus whichever papers a publisher deposited open access) and support different claims.
    """
    ix = identity.build_index(nodes.gene_id, os.path.join(OUT, "toxodb_identity.tsv"), log=log)
    identity.add_strain_accessions(
        ix, {"GT1": os.path.join(OUT, "toxodb_strain_gt1.tsv"),
             "VEG": os.path.join(OUT, "toxodb_strain_veg.tsv")}, log=log)

    docs = list(corpus.iter_abstracts(ABSTRACTS))
    if not docs:
        log(f"abstracts not found at {ABSTRACTS}; literature layer will be empty")
    n_ft = corpus.count_fulltexts(FULLTEXTS)
    if n_ft:
        log(f"{n_ft:,} open-access full texts on disk")
        docs = itertools.chain(docs, corpus.iter_fulltexts(FULLTEXTS))
    else:
        log(f"no full texts at {FULLTEXTS} (machine-local); scanning abstracts only")

    mentions, co, meta = literature.scan(docs, ix, log=log)
    if mentions.empty:
        nodes["n_publications"] = 0
        nodes["n_fulltext"] = 0
        nodes["lit_tier"] = ""
        nodes["attention_depth"] = ""
        for d in literature.DEPTH_ORDER:
            nodes[f"n_papers_{d}"] = 0
        return {}

    # The auditable intermediate: every downstream literature figure can be recomputed from this file
    # without re-scanning 40,000 documents.
    mentions.to_parquet(os.path.join(OUT, "mentions.parquet"), index=False)
    log("coverage by source and confidence tier:\n" +
        literature.coverage(mentions).to_string(index=False))

    idx = {g: i for i, g in enumerate(nodes.gene_id)}
    edges = {}
    for src, label in (("abstract", "comention"), ("fulltext", "comention_ft")):
        pairs = literature.comention_edges(co.get(src, Counter()), meta, src)
        if not pairs:
            continue
        a = np.array([idx[g1] for g1, _, _, _ in pairs])
        b = np.array([idx[g2] for _, g2, _, _ in pairs])
        w = np.array([x for _, _, x, _ in pairs])
        r = np.array([x for _, _, _, x in pairs])
        edges[label] = (a, b, w, r)
        log(f"{label}: {len(a):,} edges (raw + attention-corrected residual)")

    abstracts = literature.publication_counts(mentions, "abstract")
    fulltext = literature.publication_counts(mentions, "fulltext")
    nodes["n_publications"] = nodes.gene_id.map(abstracts).fillna(0).astype(int)
    nodes["n_fulltext"] = nodes.gene_id.map(fulltext).fillna(0).astype(int)
    # Best evidence available for each gene, so the app can say what a mention actually rests on.
    best = (mentions.assign(rank=mentions.tier.map({"accession": 0, "symbol": 1}))
            .sort_values("rank").drop_duplicates("gene_id").set_index("gene_id").tier)
    nodes["lit_tier"] = nodes.gene_id.map(best).fillna("")

    # Depth of attention: being named is not being studied (see literature.DEPTH_OF).
    depth = literature.attention_depth(mentions)
    for c in depth.columns:
        nodes[c] = nodes.gene_id.map(depth[c]).fillna(0 if c.startswith("n_") else "")
        if c.startswith("n_"):
            nodes[c] = nodes[c].astype(int)
    counts = nodes.attention_depth.value_counts()
    log("depth of attention: " + ", ".join(
        f"{counts.get(d, 0):,} {d}" for d in literature.DEPTH_ORDER) +
        f", {int((nodes.attention_depth == '').sum()):,} never named")
    return edges


# --------------------------------------------------------------------------- edges
def build_edges(nodes: pd.DataFrame):
    idx = {g: i for i, g in enumerate(nodes.gene_id)}
    edges = dict(literature_layer(nodes))

    def pair_from_groups(series, label, cap=400):
        """Edges within shared-category groups; groups larger than `cap` are skipped as uninformative."""
        a, b, w = [], [], []
        for key, grp in nodes.groupby(series):
            if not isinstance(key, str) or key in ("", "unassigned", "nan") or len(grp) < 2:
                continue
            ii = [idx[g] for g in grp.gene_id]
            if len(ii) > cap:
                continue
            for i in range(len(ii)):
                for j in range(i + 1, len(ii)):
                    a.append(ii[i]); b.append(ii[j]); w.append(1.0)
        if a:
            edges[label] = (np.array(a), np.array(b), np.array(w), np.array(w))
            log(f"{label}: {len(a):,} edges")

    if "orthogroup" in nodes.columns:
        pair_from_groups(nodes.orthogroup.fillna("").astype(str), "orthogroup", cap=60)
    pair_from_groups(nodes.compartment.astype(str), "compartment", cap=250)

    ip = os.path.join(DS, "interpro_tgon.csv")
    if os.path.exists(ip):
        d = pd.read_csv(ip, low_memory=False)
        c0 = "gene_source_id" if "gene_source_id" in d.columns else d.columns[0]
        d = d[d.interpro_id.notna()]
        dom = defaultdict(list)
        for g, i in zip(d[c0], d.interpro_id.astype(str)):
            if g in idx:
                dom[i].append(idx[g])
        a, b, w = [], [], []
        for i, ii in dom.items():
            if 2 <= len(ii) <= 60:
                for x in range(len(ii)):
                    for y in range(x + 1, len(ii)):
                        a.append(ii[x]); b.append(ii[y]); w.append(1.0)
        if a:
            edges["domain"] = (np.array(a), np.array(b), np.array(w), np.array(w))
            log(f"domain: {len(a):,} edges")

    # correlation-based edges: co-expression across stages, and co-fitness across the seven screens
    def corr_edges(cols, label, thr, kmax=25):
        M = nodes[cols].to_numpy(dtype=float)
        ok = ~np.isnan(M).any(axis=1)
        if ok.sum() < 50 or len(cols) < 3:
            return
        X = M[ok]
        X = (X - X.mean(0)) / (X.std(0) + 1e-9)
        C = np.corrcoef(X)
        np.fill_diagonal(C, 0.0)
        rows = np.where(ok)[0]
        a, b, w = [], [], []
        for i in range(len(rows)):
            order = np.argsort(-C[i])[:kmax]
            for j in order:
                if C[i, j] < thr or rows[i] >= rows[j]:
                    continue
                a.append(rows[i]); b.append(rows[j]); w.append(float(C[i, j]))
        if a:
            edges[label] = (np.array(a), np.array(b), np.array(w), np.array(w))
            log(f"{label}: {len(a):,} edges (top-{kmax} neighbours, r >= {thr})")

    rna = [c for c in nodes.columns if c.startswith("rna108740_")]
    if len(rna) >= 3:
        corr_edges(rna, "coexpression", 0.95)
    corr_edges(FIT, "cofitness", 0.90)

    # Derived last, because it is defined over the other edge types.
    hole = structural_holes(edges, nodes)
    if hole is not None:
        edges["structural_hole"] = hole
    return edges


# --------------------------------------------------------------------------- structural holes
# Independent lines of biological evidence. `orthogroup` and `domain` are ONE family, not two: paralogs
# almost always share domains, so counting them separately manufactures agreement out of one fact. That
# collapse is not cosmetic -- of the 66 sharpest candidates before it, 53 were nothing but paralogy.
#
# `compartment` is excluded entirely. Sharing one of 27 hyperLOPIT classes is real co-localisation but far
# too unspecific at 118,712 edges, and hyperLOPIT assignment tracks abundance, so it would preferentially
# link the well-expressed genes that are already well studied.
EVIDENCE_FAMILY = {"coexpression": "expression", "cofitness": "fitness",
                   "orthogroup": "homology", "domain": "homology"}

# A hole must rest on both independent *phenotype* measurements: co-expression across the stage series and
# co-fitness across the CRISPR screens. Homology can corroborate but cannot be one of the two legs.
# Allowing it as a leg admitted 291 expression+homology pairs of which 220 (76%) were same-orthogroup
# paralogs -- genes that co-express *because* they are paralogs -- and put one protein family at the top
# of every ranking. Requiring both phenotypes leaves 3 paralogs in 255 pairs.
REQUIRED_FAMILIES = frozenset({"expression", "fitness"})


def structural_holes(edges: dict, nodes: pd.DataFrame):
    """Gene pairs that two independent kinds of biology link and the literature never has.

    This is the app's governing question: not what the field says, but where its map has a gap that the
    data says should be crossed. A hole is a *derived* relation, not an observed one -- it is the absence
    of a co-mention edge across a pair the measurements agree about -- and it is labelled as such.

    Both co-mention layers count as literature, so a pair discussed anywhere, in any abstract or any
    open-access paragraph, is not a hole.
    """
    def undirected(key):
        a, b = edges[key][0], edges[key][1]
        return set(map(tuple, np.sort(np.stack([a, b], 1), axis=1))) if len(a) else set()

    fam = defaultdict(set)
    for key, family in EVIDENCE_FAMILY.items():
        if key in edges:
            for p in undirected(key):
                fam[p].add(family)
    lit = set()
    for key in ("comention", "comention_ft"):
        if key in edges:
            lit |= undirected(key)

    holes = {p: f for p, f in fam.items() if REQUIRED_FAMILIES <= f and p not in lit}
    if not holes:
        nodes["n_holes"] = 0
        return None
    a = np.array([p[0] for p in holes])
    b = np.array([p[1] for p in holes])
    w = np.array([float(len(f)) for f in holes.values()])

    deg = Counter()
    for i, j in holes:
        deg[i] += 1
        deg[j] += 1
    nodes["n_holes"] = [deg.get(i, 0) for i in range(len(nodes))]

    depth = nodes.attention_depth.to_numpy()
    studied = np.isin(depth, ["focal", "substantive"])
    sharp = int((studied[a] & studied[b]).sum())
    og = nodes.orthogroup.astype(str).to_numpy() if "orthogroup" in nodes.columns else None
    par = 0 if og is None else sum(
        1 for i, j in holes if og[i] == og[j] and og[i] not in ("", "nan", "None"))
    log(f"structural_hole: {len(a):,} pairs over {len(deg):,} genes "
        f"({sharp} with both endpoints studied -- the sharpest); "
        f"{sum(1 for f in fam.values() if REQUIRED_FAMILIES <= f) - len(holes):,} candidate pairs are "
        f"already co-mentioned and so are not holes; {par} pairs are paralogs")
    return a, b, w, w.copy()


# --------------------------------------------------------------------------- embedding
def embed(nodes: pd.DataFrame) -> np.ndarray:
    """3D UMAP of a multimodal feature matrix, so position means biological similarity."""
    feats = ["expr_tachy", "expr_cyst", "expr_max", "mean_plddt", "paralog_number",
             "n_interpro", "n_phosphosites", "has_domain", "lineage_specific"] + FIT
    feats = [f for f in feats if f in nodes.columns]
    X = nodes[feats].to_numpy(dtype=float)
    med = np.nanmedian(X, axis=0)
    X = np.where(np.isnan(X), med, X)
    X = (X - X.mean(0)) / (X.std(0) + 1e-9)
    comp = pd.get_dummies(nodes.compartment.astype(str)).to_numpy(dtype=float)
    X = np.hstack([X, comp * 0.5])          # compartment contributes, without dominating
    log(f"embedding {X.shape[0]} genes x {X.shape[1]} features")
    try:
        import umap
        Y = umap.UMAP(n_components=3, n_neighbors=25, min_dist=0.25, metric="euclidean",
                      random_state=0).fit_transform(X)
    except Exception as e:
        log(f"UMAP unavailable ({type(e).__name__}); falling back to PCA")
        from sklearn.decomposition import PCA
        Y = PCA(n_components=3, random_state=0).fit_transform(X)
    Y = np.asarray(Y, dtype=np.float32)
    Y -= Y.mean(0)
    Y /= (np.abs(Y).max() + 1e-9)
    return Y * 50.0


def main():
    nodes = load_nodes()
    edges = build_edges(nodes)
    xyz = embed(nodes)

    keep = ["gene_id", "product", "compartment", "orthogroup", "n_publications", "n_fulltext",
            "lit_tier", "attention_depth", "n_papers_focal", "n_papers_substantive",
            "n_papers_incidental", "n_holes", "has_domain",
            "lineage_specific", "paralog_number", "mean_plddt", "n_interpro", "n_phosphosites",
            "expr_tachy", "expr_cyst", "expr_max"] + FIT
    keep = [c for c in keep if c in nodes.columns]
    nodes[keep].to_parquet(os.path.join(OUT, "nodes.parquet"), index=False)

    flat = {"xyz": xyz}
    for k, (a, b, w, r) in edges.items():
        flat[f"{k}__a"] = a.astype(np.int32)
        flat[f"{k}__b"] = b.astype(np.int32)
        flat[f"{k}__w"] = w.astype(np.float32)
        flat[f"{k}__r"] = r.astype(np.float32)
    np.savez_compressed(os.path.join(OUT, "graph.npz"), **flat)
    log(f"wrote {OUT}/nodes.parquet and graph.npz "
        f"({len(nodes)} nodes, {len(edges)} edge types)")


if __name__ == "__main__":
    main()
