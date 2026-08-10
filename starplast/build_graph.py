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

import json
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

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
GENE_RX = re.compile(
    r"\b((?:GRA|ROP|RON|MIC|SRS|SAG|MYR|ASP|AP2[A-Z]*|CDPK|IMC|WNG|FIKK|SUB|CST|MORC|HDAC|BFD|TgIST|"
    r"TEEGR|HCE|MAF|PPM|RASP|NSM)\d{0,3}[A-Za-z]?)\b")


def literature_edges(nodes: pd.DataFrame):
    """Co-mention counts plus per-gene publication counts, from all 33,924 abstracts.

    Genes are matched two ways: by accession (TGME49_xxxxxx) and by the names that map to exactly one gene
    in the product annotation. Names mapping to several genes are discarded rather than guessed.
    """
    path = os.path.join(BASE, ".claude", "skills", "toxoplasma-scientist", "corpus",
                        "pubmed_toxoplasma.jsonl")
    if not os.path.exists(path):
        log("abstracts not found; co-mention edges will be empty")
        return Counter(), Counter()

    name2gene, multi = {}, set()

    def add(sym, gid):
        k = sym.upper()
        if k in name2gene and name2gene[k] != gid:
            multi.add(k)
        name2gene[k] = gid

    # ToxoDB symbols are the main source: products carry a usable symbol for only ~400 genes, and
    # matching on products alone found 234 of 8,140 genes named anywhere in the corpus.
    names_tsv = os.path.join(OUT, "toxodb_gene_names.tsv")
    if os.path.exists(names_tsv):
        # keep_default_na=False, or ToxoDB's literal "N/A" symbol is read as a missing value
        t = pd.read_csv(names_tsv, sep="\t", keep_default_na=False, dtype=str)
        t.columns = ["gene_id", "gene_name", "product", "src"][:len(t.columns)]
        for gid, sym in zip(t.gene_id, t.gene_name):
            sym = (sym or "").strip()
            if len(sym) < 3 or sym in ("N/A", "n/a"):
                continue
            # a symbol needs a digit or four capitals to be safely matchable in free text:
            # "AAP" would otherwise match prose, "GRA16" and "MORC" would not
            if not (any(c.isdigit() for c in sym) or (sym.isupper() and len(sym) >= 4)):
                continue
            add(sym, gid)
        log(f"{len(name2gene)} symbols from ToxoDB")
    else:
        log("data/toxodb_gene_names.tsv absent -- run python -m starplast.fetch_names; "
            "literature layer will badly under-count")

    for gid, prod in zip(nodes.gene_id, nodes["product"].astype(str)):
        for s in set(GENE_RX.findall(prod)):
            add(s, gid)
    for k in multi:
        name2gene.pop(k, None)
    known = set(nodes.gene_id)
    name2gene = {k: v for k, v in name2gene.items() if v in known}
    acc_rx = re.compile(r"\bTGME49_(\d{5,6}[A-Za-z]?)\b")
    sym_rx = re.compile(r"\b(" + "|".join(sorted((re.escape(k) for k in name2gene),
                                                 key=len, reverse=True)) + r")\b", re.I)
    log(f"{len(name2gene)} gene symbols resolve to exactly one gene "
        f"({len(multi)} ambiguous symbols dropped rather than guessed)")

    pubs, co = Counter(), Counter()
    n_abs = 0
    for line in open(path):
        try:
            r = json.loads(line)
        except Exception:
            continue
        n_abs += 1
        txt = (r.get("title") or "") + " " + (r.get("abstract") or "")
        hits = {f"TGME49_{m.group(1)}" for m in acc_rx.finditer(txt)}
        hits |= {name2gene[m.group(1).upper()] for m in sym_rx.finditer(txt)}
        hits &= known
        if not hits:
            continue
        for g in hits:
            pubs[g] += 1
        hs = sorted(hits)
        # abstracts naming very many genes are usually lists or screens: they inflate co-mention
        # without evidencing a relation, so they are excluded
        if len(hs) > 12:
            continue
        for i in range(len(hs)):
            for j in range(i + 1, len(hs)):
                co[(hs[i], hs[j])] += 1
    log(f"{n_abs:,} abstracts scanned; {len(pubs)} genes mentioned; {len(co):,} co-mention pairs")
    return co, pubs


# --------------------------------------------------------------------------- edges
def build_edges(nodes: pd.DataFrame):
    idx = {g: i for i, g in enumerate(nodes.gene_id)}
    edges = {}

    co, pubs = literature_edges(nodes)
    if co:
        tot = sum(pubs.values()) or 1
        a, b, raw, resid = [], [], [], []
        for (g1, g2), c in co.items():
            if c < 2:                    # a single shared abstract is not a relation
                continue
            # attention correction: expected co-mention if the two genes were mentioned independently
            exp = pubs[g1] * pubs[g2] / tot
            a.append(idx[g1]); b.append(idx[g2])
            raw.append(float(c)); resid.append(float(np.log2((c + 0.5) / (exp + 0.5))))
        edges["comention"] = (np.array(a), np.array(b), np.array(raw), np.array(resid))
        log(f"comention: {len(a):,} edges (raw + attention-corrected residual)")
        nodes["n_publications"] = nodes.gene_id.map(pubs).fillna(0).astype(int)
    else:
        nodes["n_publications"] = 0

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
    return edges


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

    keep = ["gene_id", "product", "compartment", "orthogroup", "n_publications", "has_domain",
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
