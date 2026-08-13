#!/usr/bin/env python3
"""The whole pipeline, over the real committed cache.

Everything else in this suite tests a stage against a fixture. This tests the artefact the application
actually ships, and asserts the invariants that must hold across stages — the ones that were silently
false in the past and produced output that looked fine:

    a dataset contributing zero rows          the StringTie ids, half the in vivo data
    a block collapsing to one column          the two-row merged header, 3,000 proteins
    an id resolving to nothing                pre-2012 accessions, a whole screen
    a label restating an embedding input      compartment reporting its own twin at V = 0.96
    an edge referencing a gene that is absent indexes out of range at draw time

None of these raise. Each produces a number, and the number is wrong.

The raw dataset tree is not committed, so the rebuild from source cannot run on a clean checkout. That
is a separate opt-in test at the bottom, which SKIPS loudly when the dataset root is absent and runs
when task 08's resolver finds it — a skip nobody notices is a test that never runs while everyone
believes it does.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import datasets, identity, paths, search  # noqa: E402


@pytest.fixture(scope="module")
def nodes():
    return pd.read_parquet(paths.cache_file("nodes.parquet"))


@pytest.fixture(scope="module")
def graph():
    return np.load(paths.cache_file("graph.npz"), allow_pickle=True)


@pytest.fixture(scope="module")
def mentions():
    return pd.read_parquet(paths.cache_file("mentions.parquet"))


@pytest.fixture(scope="module")
def index(nodes):
    ix = identity.build_index(nodes.gene_id, paths.cache_file("toxodb_identity.tsv"),
                              log=lambda *_: None)
    identity.add_strain_accessions(
        ix, {"GT1": paths.cache_file("toxodb_strain_gt1.tsv"),
             "VEG": paths.cache_file("toxodb_strain_veg.tsv")}, log=lambda *_: None)
    return ix


# --------------------------------------------------------------------------- the cache itself
def test_the_cache_is_complete_and_resolves_with_no_configuration():
    ok, msg = paths.check()
    assert ok, msg


def test_the_node_table_has_one_row_per_gene(nodes):
    assert len(nodes) == nodes.gene_id.nunique()
    assert nodes.gene_id.notna().all()


def test_every_gene_id_is_a_current_me49_accession(nodes):
    """A strain or pre-2012 accession surviving into the node table would join to nothing downstream
    and look like a gene nobody has studied."""
    assert nodes.gene_id.str.fullmatch(r"TGME49_\d{5,6}").all()


# --------------------------------------------------------------------------- identity
def test_every_gene_resolves_to_itself(nodes, index):
    """The floor: if a gene's own accession does not resolve, nothing keyed on it can join."""
    unresolved = [g for g in nodes.gene_id if identity.norm(g) not in index.lookup]
    assert not unresolved, f"{len(unresolved)} genes do not resolve, e.g. {unresolved[:5]}"


def test_resolution_never_invents_a_gene_that_is_not_in_the_table(nodes, index):
    known = set(nodes.gene_id)
    assert all(gid in known for gid, _ in index.lookup.values())


def test_ambiguous_strings_are_recorded_rather_than_resolved(index):
    """Every ambiguous key must be absent from the lookup: recording it and also resolving it would be
    the worst of both."""
    assert not (set(index.ambiguous) & set(index.lookup))


def test_strain_accessions_reach_the_node_table(index):
    """Published supplements cite TGGT1_ more often than TGME49_, and this mapping was missing from the
    node-column path entirely -- so every dataset keyed on a type I accession joined zero rows."""
    hits = [k for k, (_, kind) in index.lookup.items() if kind == "accession_strain"]
    assert len(hits) > 5000


# --------------------------------------------------------------------------- literature
def test_every_mentioned_gene_exists_in_the_node_table(nodes, mentions):
    assert set(mentions.gene_id) <= set(nodes.gene_id)


def test_the_two_literature_sources_are_kept_apart(mentions):
    """All of PubMed against whichever papers a publisher deposited open access: different
    populations, different claims."""
    assert set(mentions.source) <= {"abstract", "fulltext"}
    assert mentions.source.nunique() == 2


def test_most_of_the_proteome_is_named_in_no_paper_at_all(nodes):
    """The single idea the whole project is built around. If this ever stops being true, either the
    corpus grew enormously or the matcher started matching prose."""
    never_named = int((nodes.attention_depth.fillna("") == "").sum())
    assert never_named > len(nodes) * 0.6


def test_most_coverage_is_incidental_rather_than_attention(nodes):
    """The distinction the whole literature layer exists to make: appearing in a screen's hit table is
    not the same as being studied, and counting the two together overstates attention threefold."""
    tiers = nodes.attention_depth.value_counts()
    incidental = int(tiers.get("incidental", 0))
    real_attention = int(tiers.get("focal", 0)) + int(tiers.get("substantive", 0))
    assert incidental > real_attention * 2, (
        f"incidental {incidental} vs focal+substantive {real_attention}")


def test_the_attention_tiers_are_ordered_by_strength(nodes):
    """A gene's tier is the strongest place any paper names it, so focal must be the rarest."""
    tiers = nodes.attention_depth.value_counts()
    assert tiers.get("focal", 0) < tiers.get("substantive", 0) < tiers.get("incidental", 0)


# --------------------------------------------------------------------------- datasets joined in
def test_every_registered_dataset_contributes_at_least_one_column(nodes):
    """A dataset joining zero rows looks exactly like a dataset nobody measured. The in vivo screen
    contributed 0 of 8,140 rows for exactly this reason until it was routed through identity."""
    missing = {}
    for d in datasets.REGISTRY:
        if not d.columns:
            continue
        present = [c for c in d.columns if c in nodes.columns]
        if not present:
            missing[d.key] = d.columns
    assert not missing, f"registered datasets contributing nothing: {missing}"


def test_every_declared_column_actually_carries_data(nodes):
    """A column that exists and is entirely empty is worse than an absent one: it looks measured."""
    empty = []
    for d in datasets.REGISTRY:
        for c in d.columns:
            if c in nodes.columns and nodes[c].notna().sum() == 0:
                empty.append(f"{d.key}:{c}")
    assert not empty, f"declared columns with no values at all: {empty}"


def test_no_dataset_collapsed_into_a_single_column(nodes):
    """The two-row merged header collapsed a 3,000-protein abundance block into one unnamed column,
    and the dataset appeared to contribute a single column rather than fifteen."""
    for key, expect in (("proteome_total", 5), ("oocyst_itraq", 4)):
        try:
            d = datasets.get(key)
        except Exception:
            continue
        present = [c for c in d.columns if c in nodes.columns]
        if present:
            assert len(present) >= min(expect, len(d.columns))


def test_the_measured_and_derived_stage_labels_are_separate_columns(nodes):
    """One is another laboratory's measurement, the other is a restatement of columns already here.
    Merged, the derived one would be held out against an embedding built from its own sources."""
    assert "cellcycle_phase" in nodes.columns
    assert "stage_enriched_derived" in nodes.columns
    assert nodes.cellcycle_phase.notna().sum() > 500
    assert nodes.stage_enriched_derived.notna().sum() > 500


def test_the_derived_column_declares_what_it_was_computed_from(nodes):
    sources = datasets.derived_sources("stage_enriched_derived")
    assert sources
    for c in sources:
        assert c in nodes.columns, f"{c} is declared as a source but is not in the table"


# --------------------------------------------------------------------------- the graph
def test_every_edge_references_a_gene_that_exists(nodes, graph):
    """An edge to a gene that is not in the table indexes out of range at draw time."""
    n = len(nodes)
    for key in graph.files:
        if not key.endswith(("_a", "_b")):
            continue
        idx = graph[key]
        if idx.dtype.kind not in "iu" or idx.size == 0:
            continue
        assert idx.min() >= 0 and idx.max() < n, f"{key} indexes outside the node table"


def test_no_edge_joins_a_gene_to_itself(graph):
    for key in graph.files:
        if not key.endswith("_a"):
            continue
        b = key[:-2] + "_b"
        if b not in graph.files:
            continue
        a, bb = graph[key], graph[b]
        if a.dtype.kind in "iu" and a.size:
            assert not (a == bb).any(), f"{key[:-2]} contains a self-edge"


def test_no_undirected_edge_is_stored_twice(graph):
    """Stored twice, an edge draws twice and reads as twice the evidence. The `domain` layer had 106
    duplicated pairs because InterPro reports one row per MATCH, so a protein carrying two copies of
    the same domain entered its member list twice."""
    dupes = {}
    for key in graph.files:
        if not key.endswith("_a"):
            continue
        b = key[:-2] + "_b"
        if b not in graph.files:
            continue
        a, bb = graph[key], graph[b]
        if a.dtype.kind not in "iu" or not a.size:
            continue
        lo, hi = np.minimum(a, bb), np.maximum(a, bb)
        n_unique = len(set(zip(lo.tolist(), hi.tolist())))
        if n_unique != len(a):
            dupes[key[:-2]] = len(a) - n_unique
    assert not dupes, f"duplicated undirected pairs: {dupes}"


# --------------------------------------------------------------------------- the analysis layer
def test_a_held_out_target_is_genuinely_absent_from_the_embedding(nodes):
    """Measured rather than named. `compartment` fed one embedding and the battery reported its exact
    twin and its derivations as the top discoveries at V = 0.96."""
    from starplast.embedding import EmbeddingSpec, columns_for
    for target in ("compartment", "cellcycle_phase", "stage_enriched_derived"):
        banned = search.excluded_for(nodes, target)
        assert target in banned
        spec = search._spec_without(EmbeddingSpec(), nodes, banned)
        used = {c for cols in columns_for(nodes, spec).values() for c in cols}
        assert not (used & banned), f"{target}: {sorted(used & banned)} would leak into the embedding"


def test_the_declared_derivation_is_excluded_from_its_own_embedding(nodes):
    """The joint-function case no pairwise measure can catch: each source explains ~0.6 of the label
    while the label is a deterministic function of all three."""
    banned = search.excluded_for(nodes, "stage_enriched_derived")
    for c in datasets.derived_sources("stage_enriched_derived"):
        assert c in banned


def test_an_embedding_can_be_built_from_the_shipped_table(nodes):
    from starplast.embedding import EmbeddingSpec, build_matrix
    X, names, keep = build_matrix(nodes, EmbeddingSpec(), log=lambda *_: None)
    assert np.isfinite(X).all()
    assert X.shape[0] > 1000 and X.shape[1] > 10
    assert len(names) == X.shape[1]


def test_no_column_in_the_embedding_is_constant_or_all_missing(nodes):
    """A constant column contributes nothing but takes a share of the variance budget."""
    from starplast.embedding import EmbeddingSpec, build_matrix
    X, names, _ = build_matrix(nodes, EmbeddingSpec(), log=lambda *_: None)
    var = X.var(axis=0)
    dead = [n for n, v in zip(names, var) if v <= 0]
    assert not dead, f"columns carrying no variance: {dead[:5]}"


def test_the_variance_shares_sum_to_one_and_no_block_dominates_by_column_count(nodes):
    """27 one-hot compartment columns took 56% of the matrix by being numerous, which is why blocks are
    normalized before weighting."""
    from starplast.embedding import EmbeddingSpec, variance_share
    share = variance_share(nodes, EmbeddingSpec())
    assert share.share.sum() == pytest.approx(1.0, abs=1e-6)
    assert share.share.max() < 0.75


# --------------------------------------------------------------------------- the rebuild
# Probed by landmark rather than by "does a dataset root directory exist".
#
# dataset_roots() includes the user cache, which is where ensure() stages its downloads -- so
# fetching a single dataset created that directory, and every test guarded on it stopped skipping and
# started failing on a machine holding none of the raw tree. The bookkeeping ensure() leaves behind
# is not the raw inputs. Landmark files that only the real tree carries answer the actual question.
DATASET_ROOT_PRESENT = any(datasets.local_path(k) for k in ("lopit_tgon", "orthomcl"))


@pytest.mark.skipif(not DATASET_ROOT_PRESENT,
                    reason=f"no raw dataset tree found; set ${paths.ENV_DATASETS} to run the "
                           f"rebuild check (searched: {paths._REPO}/datasets, "
                           f"{os.path.dirname(paths._REPO)}/datasets, {paths.user_cache_dir()})")
def test_the_registry_resolves_against_the_real_dataset_tree():
    """Runs only where the raw inputs exist. Measured from a clean checkout before task 08, 1 of 16
    registry paths resolved and nothing reported it."""
    total = [d.key for d in datasets.REGISTRY if d.path]
    missing = datasets.missing()
    assert len(missing) < len(total) / 2, (
        f"{len(missing)} of {len(total)} registered datasets are unreachable: {missing}")


@pytest.mark.skipif(not DATASET_ROOT_PRESENT, reason="no raw dataset tree found")
def test_the_shipped_expression_columns_still_reproduce_their_geo_source(nodes, index):
    """The verification behind task 05: rho 1.0000 and 0.9970 over 7,880 genes when last measured."""
    from starplast import verify

    def resolve(a):
        hit = index.lookup.get(identity.norm(a))
        return hit[0] if hit else None

    t = verify.verify_all(nodes, resolve=resolve, log=lambda *_: None)
    if t.empty:
        pytest.skip("no GEO primary matrix on this machine")
    assert t.agrees.all(), t.to_string()
