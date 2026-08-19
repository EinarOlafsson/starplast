#!/usr/bin/env python3
"""Numbers the documents quote ABOUT THE CACHE, checked against the cache.

Prose drifts silently. The methods document said the cache carried 95 per-gene columns when it carried
169, and HANDOFF said 85; both were true once. A wrong number in a methods section is a wrong number in
a paper, and nothing about reading it reveals that.

Only claims about the shipped cache are checked, and only where the sentence is unambiguously about it.
A first attempt matched every number followed by "genes" and flagged 5,574 -- the genes named in no
paper at all -- as a wrong cache size. Per-dataset coverage figures are not claims about the cache, and
a test that cannot tell them apart would force the prose to stop quoting them.

Figures from a build log or an analysis run are recorded in instructions/done/ alongside the run that
produced them, and those files say when they were measured.
"""
from __future__ import annotations

import os
import re
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from starplast import paths  # noqa: E402

DOCS = ("README.md", "HANDOFF.md", "MATERIALS_AND_METHODS.md")


@pytest.fixture(scope="module")
def facts():
    nodes = pd.read_parquet(paths.cache_file("nodes.parquet"))
    graph = np.load(paths.cache_file("graph.npz"), allow_pickle=True)
    types = sorted({k.split("__")[0] for k in graph.files if "__" in k})
    # The SHIPPED cache, which is what the documents quote: the directories the application writes
    # into while it runs are excluded. Saved embeddings and kept clusterings are the user's own work,
    # they are gitignored, and on a machine where somebody has run a search they are 69 MB — so
    # measuring them here turns "you used this application" into "the documentation is wrong", which
    # is the opposite of what this test is for.
    written_while_running = {"embeddings", "runs", "logs"}
    size = 0
    for dp, dirs, fs in os.walk(paths.data_dir()):
        dirs[:] = [d for d in dirs if d not in written_while_running]
        size += sum(os.path.getsize(os.path.join(dp, f)) for f in fs)
    # The second species has its own cache, and a document quoting its shape should be held to it
    # rather than measured against the first. Keyed by gene count, because that is what the sentence
    # says next to the column count and it is what tells the two caches apart.
    by_genes = {len(nodes): nodes.shape[1]}
    pf = paths.cache_file("pf_nodes.parquet")
    if os.path.exists(pf):
        pf_nodes = pd.read_parquet(pf)
        by_genes[len(pf_nodes)] = pf_nodes.shape[1]
    return {"genes": len(nodes), "columns": nodes.shape[1], "columns_by_genes": by_genes,
            "edge_types": len(types), "mb": size / 1e6}


def _doc(name):
    return open(os.path.join(ROOT, name), encoding="utf8").read()


def test_the_cache_column_count_is_quoted_correctly(facts):
    """The claim that was wrong: 95 in the methods document, 85 in HANDOFF, against 169 on disk."""
    for doc in ("HANDOFF.md", "MATERIALS_AND_METHODS.md"):
        text = _doc(doc)
        pairs = re.findall(
            r"(\d{2,3})\s+(?:per-gene\s+)?columns for (?:all )?([\d,]+)(?:\s+\w+)? genes", text)
        assert pairs, f"{doc} no longer states how many columns a cache carries"
        for columns, genes in pairs:
            n = int(genes.replace(",", ""))
            assert n in facts["columns_by_genes"], (
                f"{doc} quotes {columns} columns for {genes} genes, and no cache has {genes} genes")
            assert int(columns) == facts["columns_by_genes"][n], (
                f"{doc} says {columns} columns for {genes} genes; that cache has "
                f"{facts['columns_by_genes'][n]}")
        assert facts["columns"] in {int(c) for c, _g in pairs}, (
            f"{doc} no longer states the Toxoplasma cache's column count")


def test_the_cache_size_is_quoted_correctly(facts):
    """Rounded to the megabyte in prose, so exactness is not required -- but 11 against 17 is drift."""
    for doc in ("HANDOFF.md", "MATERIALS_AND_METHODS.md"):
        text = _doc(doc)
        quoted = {int(m) for m in re.findall(r"cache (?:is|\()\s*(\d{1,3})\s*MB", text)}
        assert quoted, f"{doc} no longer states the cache size"
        for q in quoted:
            assert abs(q - facts["mb"]) < facts["mb"] / 3, (
                f"{doc} says {q} MB; the cache is {facts['mb']:.0f} MB")


def test_the_gene_total_is_quoted_correctly(facts):
    """8,140 is the one gene count that describes the whole table; the others are subsets."""
    for doc in DOCS:
        assert f"{facts['genes']:,}" in _doc(doc), f"{doc} never states the gene total"


def test_the_relation_type_count_is_quoted_correctly(facts):
    for doc in DOCS:
        text = _doc(doc)
        quoted = {int(m) for m in
                  re.findall(r"(\d{1,2})\s+(?:relation types|kinds of relation|edge types)", text)}
        assert not (quoted - {facts["edge_types"]}), (
            f"{doc} says {quoted} relation types; the graph has {facts['edge_types']}")


def test_no_document_points_at_the_pre_move_cache_path():
    """The cache moved into the package, so `data/graph.npz` names a directory that no longer exists."""
    for doc in DOCS:
        stale = re.findall(r"(?<![\w/])data/(?:graph\.npz|nodes\.parquet|toxodb_)", _doc(doc))
        assert not stale, f"{doc} points at the pre-move cache path: {set(stale)}"


def test_the_superseded_figures_do_not_come_back():
    """Each of these was true once, which is exactly why it would pass a reading."""
    for doc, gone in (("MATERIALS_AND_METHODS.md", ["95 per-gene columns", "(11 MB)"]),
                      ("HANDOFF.md", ["carries 85 columns", "cache is 11 MB"])):
        text = _doc(doc)
        for phrase in gone:
            assert phrase not in text, f"{doc} has reverted to: {phrase!r}"


def test_the_explainers_name_the_features_the_shipped_map_was_actually_built_from():
    """Both texts said "expression, fitness screens, protein features and literature co-mention",
    and the shipped layout uses no literature column at all while using the measured hyperLOPIT
    compartment, one-hot at half weight, as an input.

    Getting that backwards is not a wording slip. A reader told compartment was held out reads the
    compartment coloring as a finding, when genes of one compartment sit together partly by
    construction -- and the chat text is a system prompt, so the model repeats it."""
    from starplast import build_graph
    from starplast.app import MAP_EXPLANATION
    from starplast.chat import GROUNDING

    feats = ["expr_tachy", "expr_cyst", "expr_max", "mean_plddt", "paralog_number", "n_interpro",
             "n_phosphosites", "has_domain", "lineage_specific"] + list(build_graph.FIT)
    source = open(os.path.join(ROOT, "starplast", "build_graph.py"), encoding="utf8").read()
    body = source[source.index("def embed("):source.index("def embed(") + 1500]
    assert "nodes.compartment" in body, "the shipped embedding no longer one-hots the compartment"
    assert not [f for f in feats if f.startswith(("n_publications", "n_fulltext", "n_papers"))], \
        "a literature column became an input; the explainers say there is none"

    for name, text in (("the map explainer", MAP_EXPLANATION), ("the chat grounding", GROUNDING)):
        assert "hyperLOPIT compartment" in text, f"{name} does not say compartment is an input"
        assert "no literature column is an input" in text.lower(), \
            f"{name} still implies literature is an input"
        assert "1.1%" in text, f"{name} does not say how much of the matrix compartment carries"


def test_the_compartment_share_of_the_shipped_matrix_is_what_the_explainers_say():
    """Named as an input and doing almost nothing are both true, and quoting one without the other
    misleads in a different direction each time. Measured here rather than remembered."""
    from starplast import build_graph, paths
    cache = paths.cache_file("nodes.parquet")
    if not os.path.exists(cache):
        pytest.skip("no built cache on this machine")
    nodes = pd.read_parquet(cache)
    feats = [f for f in ["expr_tachy", "expr_cyst", "expr_max", "mean_plddt", "paralog_number",
                         "n_interpro", "n_phosphosites", "has_domain", "lineage_specific"]
             + list(build_graph.FIT) if f in nodes.columns]
    X = nodes[feats].to_numpy(dtype=float)
    med = np.nanmedian(X, axis=0)
    X = np.where(np.isnan(X), np.where(np.isfinite(med), med, 0.0), X)
    X = (X - X.mean(0)) / (X.std(0) + 1e-9)
    comp = pd.get_dummies(nodes.compartment.astype(str)).to_numpy(dtype=float) * 0.5
    share = comp.var(0).sum() / (X.var(0).sum() + comp.var(0).sum())
    assert 0.005 < share < 0.02, f"compartment now carries {share:.1%}, not the 1.1% quoted"


def test_no_shipped_column_says_the_same_thing_about_every_gene():
    """A column measured for the whole proteome and holding one value is not a measurement.

    `n_host_targets` shipped as 0 for all 8,140 genes -- the build's fallback wrote a zero per gene
    when the curated table failed to load, so the cache asserted that no protein in this parasite has
    a known host partner, on the evidence of a file that did not open. The curated table was sitting
    beside it naming 14 genes, and the slot graded A at 100% coverage on the column.

    The rule is narrow on purpose: a column covering a handful of genes may legitimately hold one
    value (the curated resistance table describes one gene), so this only fires when a column claims
    to have measured EVERY gene and still says one thing.
    """
    for name in ("nodes.parquet", "pf_nodes.parquet"):
        path = paths.cache_file(name)
        if not os.path.exists(path):
            continue
        table = pd.read_parquet(path)
        flat = [c for c in table.columns
                if table[c].notna().all() and table[c].nunique(dropna=True) == 1]
        assert not flat, f"{name}: columns with one value for every gene: {flat}"
