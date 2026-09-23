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

from . import (cellcycle, chromatin, codons, corpus, expression, identity, interaction_studies,
               phenotype_screen, resistance, drug_sensitivity, enteric,
               iedb, interactions, palmitome, small_tables,
               literature,
               proteomics,
               localization, screens, variation)

from . import paths

# Every one of these used to be computed from this file's location, which encoded one machine's layout:
# the dataset root was literally "the repository's parent directory". See paths.py.
BASE = os.path.dirname(paths.dataset_root())
DS = paths.dataset_root()
OUT = paths.data_dir()
os.makedirs(OUT, exist_ok=True)

FIT = ["fit_invitro_hff", "fit_invivo_PE", "fit_invivo_lung", "fit_invivo_liver",
       "fit_invivo_spleen", "fit_naive_bmdm", "fit_ifng"]


def log(m):
    """Print a build step with a timestamp, so a long rebuild shows where it is."""
    print(f"[build_graph] {m}", flush=True)


# --------------------------------------------------------------------------- nodes
def load_nodes() -> pd.DataFrame:
    """Assemble the node table: every gene and every measured column, joined and identifier-resolved.

    The authoritative assembly for the shipped data. Published supplements cite whatever accession
    was current when they were written, so everything joins through the identity layer -- without it
    a table keyed on TGGT1_ matches nothing and looks like a dataset with no coverage.
    """
    n = pd.read_parquet(os.path.join(BASE, "toxonet", "data", "interim", "nodes.parquet"))
    n = n.drop_duplicates("gene_id").reset_index(drop=True)
    log(f"{len(n)} genes from the node table")

    n = localization.lopit_labels(DS, n, log=log)

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

    # Every RNA column is kept, not just the three summaries. The second series (rna206344) is the
    # oocyst sporulation time course -- a whole life-stage axis that summarising to tachyzoite/cyst/max
    # threw away entirely.
    rna = [c for c in n.columns if c.startswith("rna")]
    for c in rna:
        n[c] = pd.to_numeric(n[c], errors="coerce")
    tach = [c for c in rna if "Tachyzoite" in c or "Tachyzoites" in c]
    cyst = [c for c in rna if "Tissue_cyst" in c or "Tissue_cysts" in c]
    spor = [c for c in rna if "Sporulat" in c or "Unsporulated" in c]
    n["expr_tachy"] = np.log2(n[tach].mean(axis=1) + 1) if tach else np.nan
    n["expr_cyst"] = np.log2(n[cyst].mean(axis=1) + 1) if cyst else np.nan
    n["expr_max"] = np.log2(n[rna].max(axis=1) + 1) if rna else np.nan
    n["expr_sporulated"] = np.log2(n[spor].mean(axis=1) + 1) if spor else np.nan

    # Published screens and mass-spec abundance, joined on gene id. Left missing where a targeted screen
    # never tested a gene -- absent from a 237-gene library is not a measurement of zero effect.
    #
    # Accessions go through the identity layer: papers cite whatever id was current when they were
    # written, and the 2019 in vivo screen uses pre-2012 ones for every single gene.
    ix = identity.build_index(n.gene_id, os.path.join(OUT, "toxodb_identity.tsv"),
                              log=lambda *a: None)
    # Strain accessions belong here too, not only in the literature path. Published supplements cite
    # TGGT1_ more often than TGME49_, so without this every dataset keyed on a type I accession joins
    # zero rows -- silently, since a join that matches nothing looks exactly like a dataset with no
    # coverage. The cell-cycle table is 964 rows of TGGT1_ and contributed nothing until this was added.
    identity.add_strain_accessions(
        ix, {"GT1": os.path.join(OUT, "toxodb_strain_gt1.tsv"),
             "VEG": os.path.join(OUT, "toxodb_strain_veg.tsv")}, log=lambda *a: None)

    global _SYMBOL_INDEX
    _SYMBOL_INDEX = ix.lookup

    def resolve(acc):
        hit = ix.lookup.get(identity.norm(acc))
        return hit[0] if hit else None

    for tbl in (screens.crispr_screens(BASE, log=log, resolve=resolve),
                screens.proteomics(BASE, log=log, resolve=resolve),
                screens.host_transcription_signatures(BASE, log=log, resolve=resolve),
                expression.load_all(BASE, resolve=resolve, log=log)):
        if tbl is not None and not tbl.empty:
            for c in tbl.columns:
                n[c] = n.gene_id.map(tbl[c])

    # Mass spectrometry, from deposits that were verified by reading them rather than by trusting a
    # metadata field. Merged here rather than inside `expression.load_all` because these are counts
    # of reported sites, not abundances, and mixing the two under one loader would invite them to be
    # normalised together.
    # `resolve` is not optional here, even though the parameter is. Two of the quarantine tables this
    # reaches are keyed on TGGT1_ -- the oxidative-stress screen and the calcium mineCETSA sheet --
    # and without the identity layer they join zero rows, which looks exactly like a dataset with no
    # coverage. That is how 7,384 and 2,348 values went missing from a build that still exited 0.
    ms = proteomics.load_all(BASE, n.gene_id.astype(str), log=log, resolve=resolve)
    for c in ms.columns:
        n[c] = ms[c].to_numpy()

    for table in (variation.strain_snps(BASE, resolve=resolve, log=log),
                  codons.codon_usage(BASE, resolve=resolve, log=log),
                  chromatin.chromatin_signals(BASE, resolve=resolve, log=log),
                  palmitome.palmitome(BASE, resolve=resolve, log=log),
                  small_tables.evidence(BASE, resolve=resolve, log=log),
                  small_tables.enzyme_classification(BASE, resolve=resolve,
                                                        log=log),
                  iedb.bcell_epitopes(BASE, resolve=resolve, log=log),
                  phenotype_screen.screen(BASE, resolve=resolve, log=log),
                  resistance.resistance(BASE, resolve=resolve, log=log),
                  drug_sensitivity.sensitivity(BASE, resolve=resolve, log=log),
                  enteric.fitness(BASE, resolve=resolve, log=log)):
        if table.empty:
            continue
        aligned = table.reindex(pd.Index(n.gene_id.astype(str)))
        for column in table.columns:
            n[column] = aligned[column].to_numpy()

    n = cellcycle.add_all(BASE, n, resolve=resolve, log=log)

    # The downloaded BioID/IP-MS corpus is ingested as auditable study membership, not promoted to
    # interaction edges. Most supplements contain complete quantification backgrounds; treating every
    # named protein as enriched would manufacture tens of thousands of bindings. Curated edges remain
    # separate, while the parsed rows and failures now ship for later per-paper curation.
    study_root = os.path.join(BASE, "datasets", "post_translation")
    members, studies = interaction_studies.parse_studies(study_root, resolve=resolve, log=log)
    if not studies.empty:
        studies = interaction_studies.guess_baits(
            studies, {key: value[0] for key, value in ix.lookup.items()}, log=log)
        studies.to_parquet(os.path.join(OUT, "interaction_studies.parquet"), index=False)
    if not members.empty:
        members.to_parquet(os.path.join(OUT, "interaction_study_members.parquet"), index=False)
    for method, column in (("BioID", "n_bioid_studies"), ("IPMS", "n_ipms_studies")):
        subset = members[members.method == method] if not members.empty else pd.DataFrame()
        counts = (subset.groupby("gene_id").pmid.nunique() if not subset.empty else pd.Series(dtype=int))
        n[column] = n.gene_id.map(counts).fillna(0).astype(int)

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
# Physical binding (StarPath XL-MS, IP-MS) and Foldseek structural similarity, already resolved to
# TGME49_ gene ids by the toxonet build.
TOXONET = os.path.join(BASE, "toxonet", "data", "interim", "edges_v3.parquet")


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
_SYMBOL_INDEX = {}


def _resolve_symbol(sym):
    """Symbol -> current gene id, using the identity index built during the node load."""
    hit = _SYMBOL_INDEX.get(identity.norm(sym))
    return hit[0] if hit else None


def build_edges(nodes: pd.DataFrame):
    """Build every edge layer, kept separate and never merged into one score.

    Merging would let a co-mention borrow the credibility of a measured crosslink, and no weighting
    recovers the difference afterwards -- see the Edges menu in the application.
    """
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
        # A set per domain, not a list. InterPro reports one row per MATCH, so a protein with two
        # copies of the same domain -- common, and the norm for repeat families -- appeared twice in the
        # member list. The pair loop then emitted (g, g) as an edge and duplicated every pair involving
        # it: 31 self-edges and 106 duplicated pairs in the shipped graph. A self-edge draws as a
        # zero-length line and a duplicate draws twice, reading as twice the evidence.
        dom = defaultdict(set)
        for g, i in zip(d[c0], d.interpro_id.astype(str)):
            if g in idx:
                dom[i].add(idx[g])
        a, b, w = [], [], []
        for i, members in dom.items():
            ii = sorted(members)                     # canonical order, so a < b for every pair
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
    # Co-translation is the same construction over ribosome footprints that co-expression is over
    # transcripts, and it is a different layer rather than a copy: 97% of its edges are not
    # co-expression edges (Jaccard 0.003). What says it is measuring co-translation is that
    # ribosomal proteins pair with each other 464 times where chance gives 3 -- they are made
    # together stoichiometrically, which is the textbook case.
    rpf = [c for c in nodes.columns if c.startswith("rpf")]
    if len(rpf) >= 3:
        corr_edges(rpf, "cotranslation", 0.95)

    # Measured physical binding and structural similarity. Already keyed by TGME49_ gene id upstream, so
    # these bypass the identity layer entirely.
    bind, _ = interactions.load_toxonet(TOXONET, nodes.gene_id, log=log)
    edges.update(bind)
    models = interactions.crosslink_models(BASE, set(nodes.gene_id), log=log)
    if not models.empty:
        models.to_parquet(os.path.join(OUT, "crosslink_models.parquet"), index=False)
    interactions.gene_attributes(edges, models, nodes)

    # Curated host targets. Deliberately the hand-curated table rather than anything mined from the 97
    # interaction supplements, which publish full quantification tables and would invent thousands.
    host = interaction_studies.host_interactions(BASE, resolve=lambda a: _resolve_symbol(a), log=log)
    if not host.empty:
        host.to_parquet(os.path.join(OUT, "host_interactions.parquet"), index=False)
        # Zero where the curated table names no partner, because that table is a reading of the whole
        # literature and its silence is an answer. That reasoning holds ONLY when the table is there.
        nodes["n_host_targets"] = nodes.gene_id.map(
            host.groupby("gene_id").host_target.nunique()).fillna(0).astype(int)
    else:
        # And when it is not, the column is left OUT rather than filled with zeros. Writing 0 for
        # every gene said "no host partner is known for any protein in this parasite" on the evidence
        # of a file that failed to load -- absence rendered as measurement, in the project whose first
        # rule forbids it. It shipped that way: the cache carried 8,140 zeros while the curated table
        # sat beside it naming 14 genes, and the slot graded A at 100% coverage on it.
        log("host targets: the curated table did not load, so the column is left out entirely")

    # Derived last, because these are defined over the other edge types.
    hole = structural_holes(edges, nodes)
    if hole is not None:
        edges["structural_hole"] = hole
    unwritten = unwritten_interactions(edges, nodes)
    if unwritten is not None:
        edges["unwritten_interaction"] = unwritten
    return edges


# --------------------------------------------------------------------------- structural holes
# Independent lines of biological evidence. `orthogroup` and `domain` are ONE family, not two: paralogs
# almost always share domains, so counting them separately manufactures agreement out of one fact. That
# collapse is not cosmetic -- of the 66 sharpest candidates before it, 53 were nothing but paralogy.
#
# `compartment` is excluded entirely. Sharing one of 27 hyperLOPIT classes is real co-localization but far
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
    of a co-mention edge across a pair the measurements agree about -- and it is labeled as such.

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


def unwritten_interactions(edges: dict, nodes: pd.DataFrame):
    """Pairs that were *measured* to interact and that no paper has ever discussed together.

    A stronger claim than `structural_hole`, and a different one. A hole says two independent phenotype
    measurements predict a relationship. This says a physical interaction was observed -- two residues
    covalently joined in a cell lysate, or a replicated pulldown -- and the literature still never put the
    two proteins in one sentence.

    This is a merge of `xlms` and `ip_ms`, and it is an explicit, labeled one rather than a silent
    blending of edge types: both remain separately toggleable.
    """
    src = [k for k in ("xlms", "ip_ms") if k in edges]
    if not src:
        return None

    def undirected(key):
        a, b = edges[key][0], edges[key][1]
        return set(map(tuple, np.sort(np.stack([a, b], 1), axis=1))) if len(a) else set()

    lit = set()
    for key in ("comention", "comention_ft"):
        if key in edges:
            lit |= undirected(key)

    measured = {}
    for key in src:
        a, b, w, _ = edges[key]
        for i, j, weight in zip(a, b, w):
            p = (min(i, j), max(i, j))
            measured[p] = max(measured.get(p, 0.0), float(weight))
    silent = {p: w for p, w in measured.items() if p not in lit}
    if not silent:
        return None

    a = np.array([p[0] for p in silent])
    b = np.array([p[1] for p in silent])
    w = np.array(list(silent.values()))
    depth = nodes.attention_depth.to_numpy()
    studied = np.isin(depth, ["focal", "substantive"])
    log(f"unwritten_interaction: {len(a):,} measured interactions never co-mentioned "
        f"({len(a) / max(len(measured), 1):.0%} of all measured pairs); "
        f"{int((studied[a] & studied[b]).sum())} join two genes that are each well studied")
    return a, b, w, w.copy()


# --------------------------------------------------------------------------- embedding
#: The features the Toxoplasma map is laid out on. Named rather than discovered, because a map whose
#: axes change when a dataset is added is a map two figures cannot be compared across.
TG_EMBED_FEATURES = ["expr_tachy", "expr_cyst", "expr_max", "mean_plddt", "paralog_number",
                     "n_interpro", "n_phosphosites", "has_domain", "lineage_specific"]


def embed_features(nodes: pd.DataFrame) -> list:
    """Which columns lay a table out, for whichever species' table this is.

    The Toxoplasma list is explicit and stays that way. A second species cannot reuse it -- only six
    of its sixteen columns exist in `pf_nodes.parquet`, because the two arms measure different things
    and their stage names differ on purpose -- so for anything else this takes the numeric columns
    measured for more than half the genes.

    Discovered rather than named for the second arm, and the reason is honest: naming a Plasmodium
    feature list would be inventing one, and a list invented without looking at what the data supports
    is worse than a rule that says what it did. The >50% floor is what keeps the layout from being
    driven by a column measured for forty genes.
    """
    named = [f for f in TG_EMBED_FEATURES + FIT if f in nodes.columns]
    if len(named) >= 8:                       # the Toxoplasma table, laid out as it always was
        return named
    numeric = [c for c in nodes.columns
               if c not in ("gene_id",) and pd.api.types.is_numeric_dtype(nodes[c])]
    return [c for c in numeric if float(nodes[c].notna().mean()) > 0.5]


def embed(nodes: pd.DataFrame, return_metadata: bool = False):
    """Build the packaged display through the interactive embedding recipe system."""
    from .embedding import embed as project, default_spec, EmbeddingCoordinates
    import hashlib
    spec = default_spec(nodes)
    from .embedding import columns_for, normalize
    if not columns_for(nodes, spec):
        log("no numeric layout features: using an explicitly uninformative line")
        values = np.zeros((len(nodes),3), dtype=np.float32)
        values[:,0] = np.arange(len(nodes))
        values = normalize(values) if len(nodes) else values
        metadata = {"requested_method":spec.method,"executed_method":"uninformative_line",
                    "backend":"numpy","fallback":"no numeric features","recipe":spec.to_dict(),
                    "features":[],"ordered_gene_ids":nodes.gene_id.astype(str).tolist(),
                    "coordinates_sha256":hashlib.sha256(values.tobytes()).hexdigest()}
        coords = EmbeddingCoordinates(values,metadata)
        return (coords,metadata) if return_metadata else coords
    coords, _features, rows, metadata = project(nodes, spec, log=log, return_metadata=True)
    if not np.all(rows):
        raise ValueError("a packaged layout must retain every gene")
    if coords.shape[1] < 3:
        values = np.pad(np.asarray(coords), ((0,0),(0,3-coords.shape[1])))
        metadata["display_padding_dimensions"] = 3-coords.shape[1]
        metadata["coordinates_sha256"] = hashlib.sha256(values.tobytes()).hexdigest()
        coords = EmbeddingCoordinates(values,metadata)
    return (coords, metadata) if return_metadata else coords


def main():
    """Rebuild the whole cache from the raw datasets. Minutes, and needs the raw tree."""
    from .structure_catalog import attach_features
    nodes = attach_features(load_nodes())
    edges = build_edges(nodes)
    xyz, layout_metadata = embed(nodes, return_metadata=True)

    # The app is meant to be standalone, so ship every column that survives the build rather than an
    # allowlist that silently drops whole assays -- an earlier version kept 3 of 18 RNA columns and 7 of
    # 8 fitness screens without saying so. Only genuinely internal scratch columns are dropped.
    # Exact duplicates carried in from the upstream table. Left in place they would inflate any
    # feature selection that picks "all localization columns", and double-weight that block in an
    # embedding. `compartment` is kept over `lopit_map` because it carries the explicit "unassigned".
    DROP = {"structure_path", "lopit_class", "lopit_posterior", "gene_product"}
    keep = [c for c in nodes.columns if c not in DROP]
    nodes[keep].to_parquet(os.path.join(OUT, "nodes.parquet"), index=False)

    import json
    flat = {"xyz": xyz, "gene_ids": nodes.gene_id.to_numpy(dtype=str),
            "layout_metadata": np.array(json.dumps(layout_metadata))}
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
