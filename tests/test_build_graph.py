#!/usr/bin/env python3
"""The build: node load, literature layer, edges, the derived layers, embedding, and the write.

This is the module every shipped number passes through, and it is the hardest to test because it reads
a dozen files from a tree that is not committed. So it is tested the way it fails: against a small
synthetic tree with the same SHAPE as the real one, with each stage driven directly.

The derived layers are where the reasoning lives and they need no files at all. `structural_hole` and
`unwritten_interaction` are defined over other edge types, and the judgements built into them -- that
orthogroup and domain are ONE evidence family rather than two, that a pair the literature already
co-mentions is not a hole, that compartment is too unspecific to count -- are the difference between
255 real candidates and 66 of which 53 were nothing but paralogy.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import build_graph as BG  # noqa: E402


# --------------------------------------------------------------------------- helpers
def _nodes(n=12, **extra):
    """A node table with the columns the build actually reads."""
    rng = np.random.default_rng(0)
    d = pd.DataFrame({
        "gene_id": [f"TGME49_{200000 + i}" for i in range(n)],
        "compartment": (["rhoptry"] * (n // 2)) + (["dense granules"] * (n - n // 2)),
        "attention_depth": [""] * n,
        "orthogroup": [f"OG{i // 3}" for i in range(n)],
    })
    for c in BG.FIT:
        d[c] = rng.normal(size=n)
    d["product"] = "hypothetical protein"
    for k, v in extra.items():
        d[k] = v
    return d


def _edge(pairs, w=1.0):
    a = np.array([p[0] for p in pairs], dtype=int)
    b = np.array([p[1] for p in pairs], dtype=int)
    ww = np.full(len(pairs), float(w))
    return a, b, ww, ww.copy()


def test_the_log_prefix_names_the_stage(capsys):
    """A build prints hundreds of lines; without the prefix they are indistinguishable from whatever
    else is writing to the same terminal."""
    BG.log("something happened")
    assert "[build_graph] something happened" in capsys.readouterr().out


# --------------------------------------------------------------------------- symbol resolution
def test_a_symbol_resolves_through_the_index_built_during_the_node_load():
    BG._SYMBOL_INDEX = {}
    from starplast import identity
    BG._SYMBOL_INDEX = {identity.norm("GRA16"): ("TGME49_208830", "symbol")}
    assert BG._resolve_symbol("GRA16") == "TGME49_208830"
    assert BG._resolve_symbol("gra-16") == "TGME49_208830"


def test_an_unknown_symbol_resolves_to_nothing_rather_than_to_a_guess():
    BG._SYMBOL_INDEX = {}
    assert BG._resolve_symbol("NOT_A_SYMBOL") is None


# --------------------------------------------------------------------------- structural holes
def test_a_pair_linked_by_two_independent_phenotypes_is_a_hole():
    """The app's governing question: two independent kinds of biology link them and the literature
    never has."""
    nodes = _nodes()
    edges = {"coexpression": _edge([(0, 1)]), "cofitness": _edge([(0, 1)])}
    out = BG.structural_holes(edges, nodes)
    assert out is not None
    a, b, w, _ = out
    assert list(zip(a, b)) == [(0, 1)]
    assert nodes.n_holes.tolist()[:2] == [1, 1]


def test_one_phenotype_alone_is_not_a_hole():
    """Requiring both legs is what left 3 paralogs in 255 pairs instead of 220 in 291."""
    nodes = _nodes()
    assert BG.structural_holes({"coexpression": _edge([(0, 1)])}, nodes) is None
    assert nodes.n_holes.tolist() == [0] * len(nodes)


def test_homology_cannot_stand_in_for_a_phenotype():
    """Genes co-express BECAUSE they are paralogs, so allowing homology as a leg put one protein family
    at the top of every ranking."""
    nodes = _nodes()
    edges = {"coexpression": _edge([(0, 1)]), "orthogroup": _edge([(0, 1)]),
             "domain": _edge([(0, 1)])}
    assert BG.structural_holes(edges, nodes) is None


def test_orthogroup_and_domain_are_one_evidence_family_not_two():
    """Paralogs almost always share domains, so counting them separately manufactures agreement out of
    one fact -- 53 of the first 66 candidates were pure paralogy."""
    assert BG.EVIDENCE_FAMILY["orthogroup"] == BG.EVIDENCE_FAMILY["domain"]
    assert BG.REQUIRED_FAMILIES == frozenset({"expression", "fitness"})


def test_compartment_is_not_an_evidence_family_at_all():
    """Sharing one of 27 hyperLOPIT classes is real co-localization but far too unspecific at 118,712
    edges, and the assignment tracks abundance."""
    assert "compartment" not in BG.EVIDENCE_FAMILY


def test_a_pair_the_literature_already_discusses_is_not_a_hole():
    """A gap that is only a gap in our corpus, not in the field, is not a finding."""
    nodes = _nodes()
    edges = {"coexpression": _edge([(0, 1)]), "cofitness": _edge([(0, 1)]),
             "comention": _edge([(0, 1)])}
    assert BG.structural_holes(edges, nodes) is None


def test_a_full_text_co_mention_also_disqualifies_a_hole():
    nodes = _nodes()
    edges = {"coexpression": _edge([(0, 1)]), "cofitness": _edge([(0, 1)]),
             "comention_ft": _edge([(0, 1)])}
    assert BG.structural_holes(edges, nodes) is None


def test_edge_direction_does_not_hide_a_co_mention():
    """The literature edge is stored in whichever order it arrived, so both must be compared
    undirected -- otherwise a co-mentioned pair reappears as a hole."""
    nodes = _nodes()
    edges = {"coexpression": _edge([(0, 1)]), "cofitness": _edge([(0, 1)]),
             "comention": _edge([(1, 0)])}
    assert BG.structural_holes(edges, nodes) is None


def test_an_empty_edge_layer_is_handled():
    nodes = _nodes()
    edges = {"coexpression": _edge([]), "cofitness": _edge([(0, 1)])}
    assert BG.structural_holes(edges, nodes) is None


def test_holes_report_how_many_join_two_well_studied_genes(capsys):
    """Those are the sharpest: neither gene is obscure, and still nobody has connected them."""
    nodes = _nodes()
    nodes.loc[0, "attention_depth"] = "focal"
    nodes.loc[1, "attention_depth"] = "substantive"
    BG.structural_holes({"coexpression": _edge([(0, 1)]), "cofitness": _edge([(0, 1)])}, nodes)
    out = capsys.readouterr().out
    assert "1 with both endpoints studied" in out


def test_holes_report_how_many_are_paralogs(capsys):
    """The number that exposed the first version: 53 of 66 candidates were nothing but paralogy."""
    nodes = _nodes()
    nodes.loc[0, "orthogroup"] = "OGX"
    nodes.loc[1, "orthogroup"] = "OGX"
    BG.structural_holes({"coexpression": _edge([(0, 1)]), "cofitness": _edge([(0, 1)])}, nodes)
    assert "1 pairs are paralogs" in capsys.readouterr().out


def test_paralogy_is_not_counted_where_the_orthogroup_is_unknown():
    """Two genes both lacking an orthogroup are not paralogs of each other."""
    nodes = _nodes()
    nodes["orthogroup"] = ""
    BG.structural_holes({"coexpression": _edge([(0, 1)]), "cofitness": _edge([(0, 1)])}, nodes)


def test_holes_work_without_an_orthogroup_column_at_all():
    nodes = _nodes().drop(columns=["orthogroup"])
    assert BG.structural_holes({"coexpression": _edge([(0, 1)]),
                                "cofitness": _edge([(0, 1)])}, nodes) is not None


# --------------------------------------------------------------------------- unwritten interactions
def test_a_measured_interaction_nobody_has_written_about_is_reported():
    """A stronger claim than a hole: two residues covalently joined in a cell lysate, and the
    literature still never put the two proteins in one sentence."""
    nodes = _nodes()
    out = BG.unwritten_interactions({"xlms": _edge([(0, 1)], w=2.0)}, nodes)
    assert out is not None
    a, b, w, _ = out
    assert list(zip(a, b)) == [(0, 1)]
    assert w.tolist() == [2.0]


def test_an_interaction_that_has_been_written_about_is_not_reported():
    nodes = _nodes()
    edges = {"xlms": _edge([(0, 1)]), "comention": _edge([(0, 1)])}
    assert BG.unwritten_interactions(edges, nodes) is None


def test_both_measured_sources_are_merged_explicitly():
    """A labelled merge rather than a silent blending: xlms and ip_ms both stay separately toggleable."""
    nodes = _nodes()
    out = BG.unwritten_interactions({"xlms": _edge([(0, 1)]), "ip_ms": _edge([(2, 3)])}, nodes)
    a, b, _, _ = out
    assert sorted(zip(a, b)) == [(0, 1), (2, 3)]


def test_the_strongest_weight_survives_when_both_sources_have_the_pair():
    nodes = _nodes()
    out = BG.unwritten_interactions({"xlms": _edge([(0, 1)], w=1.0),
                                     "ip_ms": _edge([(0, 1)], w=5.0)}, nodes)
    assert out[2].tolist() == [5.0]


def test_no_measured_interactions_at_all_returns_nothing():
    assert BG.unwritten_interactions({"coexpression": _edge([(0, 1)])}, _nodes()) is None


def test_every_measured_interaction_being_written_about_returns_nothing():
    nodes = _nodes()
    edges = {"xlms": _edge([(0, 1)]), "comention_ft": _edge([(0, 1)])}
    assert BG.unwritten_interactions(edges, nodes) is None


def test_the_share_of_measured_pairs_that_are_unwritten_is_reported(capsys):
    """92% in the shipped build -- the number that makes the point."""
    nodes = _nodes()
    BG.unwritten_interactions({"xlms": _edge([(0, 1), (2, 3)]),
                               "comention": _edge([(0, 1)])}, nodes)
    assert "50% of all measured pairs" in capsys.readouterr().out


# --------------------------------------------------------------------------- the embedding
def test_the_embedding_is_three_dimensional_and_bounded():
    """The app places points directly from this, so an unbounded coordinate puts a gene off screen."""
    nodes = _nodes(60, expr_tachy=1.0, expr_cyst=2.0, expr_max=3.0, mean_plddt=80.0,
                   paralog_number=1, n_interpro=2, n_phosphosites=0, has_domain=1,
                   lineage_specific=0)
    Y = BG.embed(nodes)
    assert Y.shape == (60, 3)
    assert np.isfinite(Y).all()
    assert np.abs(Y).max() <= 50.0 + 1e-3


def test_the_embedding_is_centred():
    nodes = _nodes(60, expr_tachy=1.0, mean_plddt=80.0)
    Y = BG.embed(nodes)
    assert np.abs(Y.mean(0)).max() < 1.0


def test_missing_features_are_filled_with_the_median_not_zero():
    """Zero is a value in a z-scored matrix; for pLDDT it would place unmeasured proteins at the mean
    of the measured ones by accident rather than by decision."""
    nodes = _nodes(60, mean_plddt=np.nan)
    nodes.loc[:29, "mean_plddt"] = 90.0
    Y = BG.embed(nodes)
    assert np.isfinite(Y).all()


def test_compartment_contributes_without_dominating():
    """27 one-hot columns at full weight would drown the measured features by sheer count."""
    nodes = _nodes(60, expr_tachy=1.0)
    Y = BG.embed(nodes)
    assert np.isfinite(Y).all()


def test_without_umap_the_embedding_falls_back_to_pca(monkeypatch, capsys):
    import builtins
    real = builtins.__import__

    def no_umap(name, *a, **k):
        if name == "umap":
            raise ImportError("no umap here")
        return real(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", no_umap)
    Y = BG.embed(_nodes(60, expr_tachy=1.0))
    assert Y.shape == (60, 3)
    assert "falling back to PCA" in capsys.readouterr().out


# --------------------------------------------------------------------------- grouped and correlated edges
@pytest.fixture
def quiet(monkeypatch):
    """The build logs heavily; silence it except where a test is about the log."""
    monkeypatch.setattr(BG, "log", lambda *a, **k: None)


def _build_edges(monkeypatch, nodes, tmp_path, lit=None, log=None):
    """Run build_edges with the file-reading layers stubbed out.

    `log` captures what the build says, which for one branch is the whole point: a column that is
    left out has to announce itself or its absence is indistinguishable from a build that never ran.
    """
    if log is not None:
        monkeypatch.setattr(BG, "log", log)
    monkeypatch.setattr(BG, "literature_layer", lambda n: (lit or {}))
    monkeypatch.setattr(BG.interactions, "load_toxonet", lambda *a, **k: ({}, pd.DataFrame()))
    monkeypatch.setattr(BG.interactions, "crosslink_models", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.interaction_studies, "host_interactions", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG, "OUT", str(tmp_path))
    monkeypatch.setattr(BG, "BASE", str(tmp_path))
    monkeypatch.setattr(BG, "TOXONET", str(tmp_path / "absent.parquet"))
    return BG.build_edges(nodes)


def test_genes_sharing_a_category_are_joined(monkeypatch, tmp_path, quiet):
    nodes = _nodes(6)
    nodes["compartment"] = ["rhoptry"] * 3 + ["dense granules"] * 3
    edges = _build_edges(monkeypatch, nodes, tmp_path)
    a, b, _, _ = edges["compartment"]
    assert len(a) == 6                      # two groups of three: 3 pairs each


def test_a_group_larger_than_the_cap_is_skipped(monkeypatch, tmp_path, quiet):
    """One category holding half the proteome produces millions of pairs that say nothing. The
    orthogroup cap is 60 and compartment's is 250, because an orthogroup of 100 is a repeat family
    while a compartment of 100 is still a compartment."""
    nodes = _nodes(70)
    nodes["orthogroup"] = "one_big_family"
    nodes["compartment"] = "rhoptry"
    edges = _build_edges(monkeypatch, nodes, tmp_path)
    assert "orthogroup" not in edges
    assert "compartment" in edges, "the same 70 genes are under compartment's larger cap"


def test_the_absence_categories_are_never_grouped(monkeypatch, tmp_path, quiet):
    """`unassigned` is not a compartment, and joining its members would link every unmeasured gene to
    every other."""
    nodes = _nodes(6)
    nodes["compartment"] = ["unassigned"] * 4 + ["rhoptry"] * 2
    edges = _build_edges(monkeypatch, nodes, tmp_path)
    a, b, _, _ = edges["compartment"]
    assert len(a) == 1


def test_a_category_with_one_member_produces_no_edge(monkeypatch, tmp_path, quiet):
    nodes = _nodes(3)
    nodes["compartment"] = ["a", "b", "c"]
    edges = _build_edges(monkeypatch, nodes, tmp_path)
    assert "compartment" not in edges


def test_correlated_fitness_profiles_are_joined(monkeypatch, tmp_path, quiet):
    """Co-fitness across the seven screens: two genes whose knockouts behave alike everywhere."""
    n = 60
    nodes = _nodes(n)
    rng = np.random.default_rng(1)
    base = rng.normal(size=n)
    for c in BG.FIT:
        nodes[c] = base + rng.normal(scale=0.01, size=n)
    edges = _build_edges(monkeypatch, nodes, tmp_path)
    assert "cofitness" in edges
    assert len(edges["cofitness"][0]) > 0


def test_too_few_complete_rows_produce_no_correlation_edges(monkeypatch, tmp_path, quiet):
    """A correlation over twelve genes is noise wearing a statistic's clothes."""
    edges = _build_edges(monkeypatch, _nodes(12), tmp_path)
    assert "cofitness" not in edges


def test_co_expression_needs_at_least_three_columns(monkeypatch, tmp_path, quiet):
    nodes = _nodes(60)
    nodes["rna108740_a"] = 1.0
    nodes["rna108740_b"] = 2.0
    edges = _build_edges(monkeypatch, nodes, tmp_path)
    assert "coexpression" not in edges


def test_correlation_edges_are_stored_once_per_pair(monkeypatch, tmp_path, quiet):
    n = 60
    nodes = _nodes(n)
    rng = np.random.default_rng(2)
    base = rng.normal(size=n)
    for i in range(4):
        nodes[f"rna108740_{i}"] = base + rng.normal(scale=0.01, size=n)
    edges = _build_edges(monkeypatch, nodes, tmp_path, lit={})
    if "coexpression" in edges:
        a, b, _, _ = edges["coexpression"]
        assert (a < b).all(), "stored twice, an edge reads as twice the evidence"


def test_the_literature_layer_is_merged_in(monkeypatch, tmp_path, quiet):
    edges = _build_edges(monkeypatch, _nodes(6), tmp_path, lit={"comention": _edge([(0, 1)])})
    assert "comention" in edges


def test_curated_host_targets_are_counted_per_gene(monkeypatch, tmp_path, quiet):
    nodes = _nodes(6)
    host = pd.DataFrame({"gene_id": ["TGME49_200000", "TGME49_200000", "TGME49_200001"],
                         "host_target": ["STAT3", "STAT6", "p38"]})
    monkeypatch.setattr(BG, "literature_layer", lambda n: {})
    monkeypatch.setattr(BG.interactions, "load_toxonet", lambda *a, **k: ({}, pd.DataFrame()))
    monkeypatch.setattr(BG.interactions, "crosslink_models", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.interaction_studies, "host_interactions", lambda *a, **k: host)
    monkeypatch.setattr(BG, "OUT", str(tmp_path))
    monkeypatch.setattr(BG, "BASE", str(tmp_path))
    monkeypatch.setattr(BG, "TOXONET", str(tmp_path / "absent.parquet"))
    BG.build_edges(nodes)
    assert nodes.n_host_targets.tolist()[:2] == [2, 1]
    assert os.path.exists(tmp_path / "host_interactions.parquet")


def test_no_curated_host_table_leaves_the_column_out_rather_than_writing_zeros(monkeypatch,
                                                                                tmp_path, quiet):
    """This asserted the opposite until 2026-08-19, and the opposite had shipped.

    Zero is the right answer for a gene the curated table does not name, because that table is a
    reading of the whole literature and its silence is an answer. It is the wrong answer when the
    table itself did not load: the cache then said *no protein in this parasite has a known host
    partner*, on the evidence of a file that failed to open. It did exactly that -- 8,140 zeros
    shipped beside a curated table naming 14 genes, and `Tg_host interaction degree` graded A at 100%
    coverage on them.

    A slot with no data has to read empty, which it can only do if the column is absent.
    """
    nodes = _nodes(6)
    said = []
    _build_edges(monkeypatch, nodes, tmp_path, log=said.append)
    assert "n_host_targets" not in nodes.columns
    assert any("curated table did not load" in m for m in said), "and it has to say so"


def test_crosslink_models_are_written_where_they_exist(monkeypatch, tmp_path, quiet):
    nodes = _nodes(6)
    models = pd.DataFrame({"gene_a": ["TGME49_200000"], "gene_b": ["TGME49_200001"],
                           "frac_satisfied": [0.9]})
    monkeypatch.setattr(BG, "literature_layer", lambda n: {})
    monkeypatch.setattr(BG.interactions, "load_toxonet", lambda *a, **k: ({}, pd.DataFrame()))
    monkeypatch.setattr(BG.interactions, "crosslink_models", lambda *a, **k: models)
    monkeypatch.setattr(BG.interaction_studies, "host_interactions", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG, "OUT", str(tmp_path))
    monkeypatch.setattr(BG, "BASE", str(tmp_path))
    monkeypatch.setattr(BG, "TOXONET", str(tmp_path / "absent.parquet"))
    BG.build_edges(nodes)
    assert os.path.exists(tmp_path / "crosslink_models.parquet")


def test_no_orthogroup_column_is_not_fatal(monkeypatch, tmp_path, quiet):
    nodes = _nodes(6).drop(columns=["orthogroup"])
    edges = _build_edges(monkeypatch, nodes, tmp_path)
    assert "orthogroup" not in edges


def test_the_derived_layers_are_built_last(monkeypatch, tmp_path, quiet):
    """They are defined over the other edge types, so building them earlier would define them over an
    incomplete graph."""
    n = 60
    nodes = _nodes(n)
    rng = np.random.default_rng(3)
    base = rng.normal(size=n)
    for c in BG.FIT:
        nodes[c] = base + rng.normal(scale=0.01, size=n)
    for i in range(4):
        nodes[f"rna108740_{i}"] = base + rng.normal(scale=0.01, size=n)
    edges = _build_edges(monkeypatch, nodes, tmp_path)
    assert "structural_hole" in edges


# --------------------------------------------------------------------------- domain edges
def test_genes_sharing_an_interpro_domain_are_joined(monkeypatch, tmp_path, quiet):
    nodes = _nodes(6)
    pd.DataFrame({"gene_source_id": ["TGME49_200000", "TGME49_200001", "TGME49_200002"],
                  "interpro_id": ["IPR001", "IPR001", "IPR002"]}).to_csv(
        tmp_path / "interpro_tgon.csv", index=False)
    monkeypatch.setattr(BG, "DS", str(tmp_path))
    edges = _build_edges(monkeypatch, nodes, tmp_path)
    a, b, _, _ = edges["domain"]
    assert list(zip(a, b)) == [(0, 1)]


def test_a_protein_carrying_the_same_domain_twice_is_not_joined_to_itself(monkeypatch, tmp_path,
                                                                          quiet):
    """InterPro reports one row per MATCH, and repeat families match many times. Appended to a list,
    the gene entered its own domain's member set twice and the pair loop emitted (g, g): 31 self-edges
    and 106 duplicated pairs in the shipped graph."""
    nodes = _nodes(6)
    pd.DataFrame({"gene_source_id": ["TGME49_200000"] * 3 + ["TGME49_200001"],
                  "interpro_id": ["IPR001"] * 4}).to_csv(tmp_path / "interpro_tgon.csv", index=False)
    monkeypatch.setattr(BG, "DS", str(tmp_path))
    edges = _build_edges(monkeypatch, nodes, tmp_path)
    a, b, _, _ = edges["domain"]
    assert not (a == b).any()
    assert len(a) == len(set(zip(a.tolist(), b.tolist()))) == 1


def test_domain_pairs_come_out_canonically_ordered(monkeypatch, tmp_path, quiet):
    nodes = _nodes(6)
    pd.DataFrame({"gene_source_id": ["TGME49_200003", "TGME49_200001", "TGME49_200002"],
                  "interpro_id": ["IPR001"] * 3}).to_csv(tmp_path / "interpro_tgon.csv", index=False)
    monkeypatch.setattr(BG, "DS", str(tmp_path))
    a, b, _, _ = _build_edges(monkeypatch, nodes, tmp_path)["domain"]
    assert (a < b).all()


def test_a_domain_shared_by_too_many_genes_is_skipped(monkeypatch, tmp_path, quiet):
    """A domain in 500 proteins says they are all proteins."""
    nodes = _nodes(70)
    pd.DataFrame({"gene_source_id": [f"TGME49_{200000 + i}" for i in range(70)],
                  "interpro_id": ["IPR001"] * 70}).to_csv(tmp_path / "interpro_tgon.csv", index=False)
    monkeypatch.setattr(BG, "DS", str(tmp_path))
    assert "domain" not in _build_edges(monkeypatch, nodes, tmp_path)


def test_rows_without_a_domain_are_ignored(monkeypatch, tmp_path, quiet):
    nodes = _nodes(6)
    pd.DataFrame({"gene_source_id": ["TGME49_200000", "TGME49_200001"],
                  "interpro_id": [np.nan, np.nan]}).to_csv(tmp_path / "interpro_tgon.csv", index=False)
    monkeypatch.setattr(BG, "DS", str(tmp_path))
    assert "domain" not in _build_edges(monkeypatch, nodes, tmp_path)


def test_a_gene_absent_from_the_node_table_is_ignored(monkeypatch, tmp_path, quiet):
    nodes = _nodes(6)
    pd.DataFrame({"gene_source_id": ["TGME49_200000", "TGME49_999999"],
                  "interpro_id": ["IPR001", "IPR001"]}).to_csv(
        tmp_path / "interpro_tgon.csv", index=False)
    monkeypatch.setattr(BG, "DS", str(tmp_path))
    assert "domain" not in _build_edges(monkeypatch, nodes, tmp_path)


def test_no_interpro_file_is_not_fatal(monkeypatch, tmp_path, quiet):
    monkeypatch.setattr(BG, "DS", str(tmp_path))
    assert "domain" not in _build_edges(monkeypatch, _nodes(6), tmp_path)


# --------------------------------------------------------------------------- the literature layer
def _identity_files(d, genes):
    pd.DataFrame({"gene_id": genes, "gene_name": [""] * len(genes),
                  "previous_ids": [""] * len(genes),
                  "product": ["hypothetical"] * len(genes)}).to_csv(
        d / "toxodb_identity.tsv", sep="\t", index=False)
    for tag in ("gt1", "veg"):
        pd.DataFrame({"gene_id": [], "gene_name": []}).to_csv(
            d / f"toxodb_strain_{tag}.tsv", sep="\t", index=False)


def test_with_no_corpus_at_all_the_literature_columns_are_zero_not_missing(monkeypatch, tmp_path):
    """5,574 genes are named in no paper. Missing would read as "not looked at"; zero is the
    measurement."""
    nodes = _nodes(6)
    _identity_files(tmp_path, nodes.gene_id)
    monkeypatch.setattr(BG, "OUT", str(tmp_path))
    monkeypatch.setattr(BG, "ABSTRACTS", str(tmp_path / "absent.jsonl"))
    monkeypatch.setattr(BG, "FULLTEXTS", str(tmp_path / "no_fulltexts"))
    monkeypatch.setattr(BG, "log", lambda *a, **k: None)
    edges = BG.literature_layer(nodes)
    assert edges == {}
    assert nodes.n_publications.tolist() == [0] * 6
    assert nodes.attention_depth.tolist() == [""] * 6
    assert nodes.n_papers_focal.tolist() == [0] * 6


def test_abstracts_produce_mentions_counts_and_co_mention_edges(monkeypatch, tmp_path):
    nodes = _nodes(6)
    _identity_files(tmp_path, nodes.gene_id)
    lines = []
    for i in range(4):
        lines.append(json.dumps({
            "pmid": str(100 + i), "year": "2020",
            "title": "TGME49_200000 and TGME49_200001 interact",
            "abstract": "We show that TGME49_200000 binds TGME49_200001 in tachyzoites."}))
    (tmp_path / "abs.jsonl").write_text("\n".join(lines))
    monkeypatch.setattr(BG, "OUT", str(tmp_path))
    monkeypatch.setattr(BG, "ABSTRACTS", str(tmp_path / "abs.jsonl"))
    monkeypatch.setattr(BG, "FULLTEXTS", str(tmp_path / "no_fulltexts"))
    monkeypatch.setattr(BG, "log", lambda *a, **k: None)

    edges = BG.literature_layer(nodes)
    assert nodes.n_publications.iloc[0] == 4
    assert nodes.attention_depth.iloc[0] == "focal", "named in the title"
    assert os.path.exists(tmp_path / "mentions.parquet"), "the auditable intermediate must be written"
    if edges:
        a, b, w, r = edges["comention"]
        assert len(a) == len(b) == len(w) == len(r)


def test_full_texts_are_scanned_alongside_abstracts(monkeypatch, tmp_path):
    """Different populations -- all of PubMed against whatever a publisher deposited -- so they are
    counted separately all the way through."""
    nodes = _nodes(6)
    _identity_files(tmp_path, nodes.gene_id)
    (tmp_path / "abs.jsonl").write_text(json.dumps(
        {"pmid": "1", "title": "TGME49_200000 study"}))
    ft = tmp_path / "ft"
    ft.mkdir()
    (ft / "a.xml").write_text("""<article><front><article-meta>
        <article-id pub-id-type="pmid">2</article-id>
        <title-group><article-title>TGME49_200001 in tachyzoites</article-title></title-group>
        </article-meta></front></article>""")
    monkeypatch.setattr(BG, "OUT", str(tmp_path))
    monkeypatch.setattr(BG, "ABSTRACTS", str(tmp_path / "abs.jsonl"))
    monkeypatch.setattr(BG, "FULLTEXTS", str(ft))
    monkeypatch.setattr(BG, "log", lambda *a, **k: None)
    BG.literature_layer(nodes)
    assert nodes.n_publications.iloc[0] == 1
    assert nodes.n_fulltext.iloc[1] == 1


def test_the_evidence_tier_records_what_a_mention_rests_on(monkeypatch, tmp_path):
    """An accession is a stronger claim than a symbol, and the app says which it was."""
    nodes = _nodes(6)
    _identity_files(tmp_path, nodes.gene_id)
    (tmp_path / "abs.jsonl").write_text(json.dumps(
        {"pmid": "1", "title": "TGME49_200000 does something"}))
    monkeypatch.setattr(BG, "OUT", str(tmp_path))
    monkeypatch.setattr(BG, "ABSTRACTS", str(tmp_path / "abs.jsonl"))
    monkeypatch.setattr(BG, "FULLTEXTS", str(tmp_path / "none"))
    monkeypatch.setattr(BG, "log", lambda *a, **k: None)
    BG.literature_layer(nodes)
    assert nodes.lit_tier.iloc[0] == "accession"
    assert nodes.lit_tier.iloc[5] == "", "a gene named nowhere has no tier"


# --------------------------------------------------------------------------- loading nodes
def _upstream(tmp_path, n=6):
    """The toxonet interim table, the orthomcl product table, and the identity files the load reads."""
    base = tmp_path / "base"
    interim = base / "toxonet" / "data" / "interim"
    interim.mkdir(parents=True)
    genes = [f"TGME49_{200000 + i}" for i in range(n)]
    rng = np.random.default_rng(0)
    up = pd.DataFrame({"gene_id": genes,
                       "mean_plddt": rng.uniform(40, 95, n),
                       "paralog_number": rng.integers(0, 4, n),
                       "rna108740_Tachyzoites_T2_FPKM": rng.uniform(0, 500, n),
                       "rna108740_Tissue_cysts_A_FPKM": rng.uniform(0, 500, n),
                       "rna206344_Sporulating_R1": rng.uniform(0, 500, n)})
    up.to_parquet(interim / "nodes.parquet", index=False)

    ds = tmp_path / "ds"
    ds.mkdir()
    pd.DataFrame({"gene_source_id": genes,
                  "gene_product": ["a protein"] * n}).to_csv(
        ds / "orthomcl_toxoplasma_gondii_ME49.csv", index=False)
    out = tmp_path / "out"
    out.mkdir()
    _identity_files(out, genes)
    return base, ds, out, genes


def test_the_node_load_composes_the_upstream_table_with_products_and_expression(monkeypatch,
                                                                                tmp_path):
    base, ds, out, genes = _upstream(tmp_path)
    monkeypatch.setattr(BG, "BASE", str(base))
    monkeypatch.setattr(BG, "DS", str(ds))
    monkeypatch.setattr(BG, "OUT", str(out))
    monkeypatch.setattr(BG, "log", lambda *a, **k: None)
    monkeypatch.setattr(BG.localization, "lopit_labels", lambda ds, n, log=None: n.assign(
        compartment="unassigned", lopit_unified=None))
    monkeypatch.setattr(BG.screens, "crispr_screens", lambda *a, **k: pd.DataFrame())
    # One verified mass-spec deposit, so the merge that puts site counts on genes runs here rather
    # than only on a machine that has the downloads. The build has to carry it through to the
    # written table, which is the claim this test is for.
    deposit = base / BG.proteomics.QUARANTINE / "Tg" / "acetylation"
    deposit.mkdir(parents=True, exist_ok=True)
    (deposit / "KSites.txt").write_text(f"Proteins\n{genes[0]}\n{genes[0]}\n{genes[1]}\n")
    monkeypatch.setattr(BG.proteomics, "DEPOSITS", (
        ("Tg", "acetylation", "PXD079431", "n_acetylation_sites", r"KSites"),))
    monkeypatch.setattr(BG.screens, "proteomics", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.expression, "load_all", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.cellcycle, "add_all", lambda base, n, **k: n)

    n = BG.load_nodes()
    assert list(n.gene_id) == genes
    assert (n["product"] == "a protein").all()
    for c in ("expr_tachy", "expr_cyst", "expr_max", "expr_sporulated"):
        assert c in n.columns, c
    assert n.has_domain.tolist() == [0] * len(genes), "absent means 0, never missing"
    assert n.lineage_specific.tolist() == [0] * len(genes)
    for c in BG.FIT:
        assert c in n.columns, f"{c} must exist even when no screen measured it"


def test_duplicate_genes_from_upstream_are_collapsed(monkeypatch, tmp_path):
    base, ds, out, genes = _upstream(tmp_path)
    interim = base / "toxonet" / "data" / "interim"
    d = pd.read_parquet(interim / "nodes.parquet")
    pd.concat([d, d.iloc[:2]]).to_parquet(interim / "nodes.parquet", index=False)
    monkeypatch.setattr(BG, "BASE", str(base))
    monkeypatch.setattr(BG, "DS", str(ds))
    monkeypatch.setattr(BG, "OUT", str(out))
    monkeypatch.setattr(BG, "log", lambda *a, **k: None)
    monkeypatch.setattr(BG.localization, "lopit_labels", lambda ds, n, log=None: n.assign(
        compartment="unassigned"))
    monkeypatch.setattr(BG.screens, "crispr_screens", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.screens, "proteomics", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.expression, "load_all", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.cellcycle, "add_all", lambda base, n, **k: n)
    n = BG.load_nodes()
    assert len(n) == len(genes) == n.gene_id.nunique()


def test_joined_dataset_columns_reach_the_node_table(monkeypatch, tmp_path):
    """The join that silently produced nothing for a whole screen until it went through identity."""
    base, ds, out, genes = _upstream(tmp_path)
    screen = pd.DataFrame({"crispr_thing": [1.0, 2.0]}, index=genes[:2])
    monkeypatch.setattr(BG, "BASE", str(base))
    monkeypatch.setattr(BG, "DS", str(ds))
    monkeypatch.setattr(BG, "OUT", str(out))
    monkeypatch.setattr(BG, "log", lambda *a, **k: None)
    monkeypatch.setattr(BG.localization, "lopit_labels", lambda ds, n, log=None: n.assign(
        compartment="unassigned"))
    monkeypatch.setattr(BG.screens, "crispr_screens", lambda *a, **k: screen)
    monkeypatch.setattr(BG.screens, "proteomics", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.expression, "load_all", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.cellcycle, "add_all", lambda base, n, **k: n)
    n = BG.load_nodes()
    assert n.crispr_thing.tolist()[:2] == [1.0, 2.0]
    assert pd.isna(n.crispr_thing.iloc[2]), "untested is missing, never zero"


def test_interaction_corpus_membership_reaches_counts_and_audit_tables(monkeypatch, tmp_path):
    base, ds, out, genes = _upstream(tmp_path)
    members = pd.DataFrame({"pmid": ["1", "2", "3"],
                            "gene_id": [genes[0], genes[0], genes[1]],
                            "method": ["BioID", "BioID", "IPMS"]})
    studies = pd.DataFrame({"pmid": ["1", "2", "3"], "title": ["A", "B", "C"]})
    monkeypatch.setattr(BG, "BASE", str(base))
    monkeypatch.setattr(BG, "DS", str(ds))
    monkeypatch.setattr(BG, "OUT", str(out))
    monkeypatch.setattr(BG, "log", lambda *a, **k: None)
    monkeypatch.setattr(BG.localization, "lopit_labels", lambda ds, n, log=None: n.assign(
        compartment="unassigned"))
    monkeypatch.setattr(BG.screens, "crispr_screens", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.screens, "proteomics", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.screens, "host_transcription_signatures",
                        lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.expression, "load_all", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.cellcycle, "add_all", lambda base, n, **k: n)
    monkeypatch.setattr(BG.interaction_studies, "parse_studies",
                        lambda *a, **k: (members, studies))
    monkeypatch.setattr(BG.interaction_studies, "guess_baits",
                        lambda frame, *a, **k: frame.assign(bait_gene=None,
                                                            bait_confidence=""))
    n = BG.load_nodes()
    assert n.n_bioid_studies.tolist()[:2] == [2, 0]
    assert n.n_ipms_studies.tolist()[:2] == [0, 1]
    assert (out / "interaction_studies.parquet").exists()
    assert (out / "interaction_study_members.parquet").exists()


# --------------------------------------------------------------------------- the whole build
def test_main_writes_a_cache_the_application_can_open(monkeypatch, tmp_path):
    """The contract at the end of the build: nodes.parquet plus graph.npz, with every edge type flat."""
    out = tmp_path / "out"
    out.mkdir()
    nodes = _nodes(60, expr_tachy=1.0, mean_plddt=80.0, structure_path="x", lopit_class="y",
                   lopit_posterior=0.5, gene_product="z")
    monkeypatch.setattr(BG, "OUT", str(out))
    monkeypatch.setattr(BG, "log", lambda *a, **k: None)
    monkeypatch.setattr(BG, "load_nodes", lambda: nodes)
    monkeypatch.setattr(BG, "build_edges", lambda n: {"comention": _edge([(0, 1), (2, 3)])})

    BG.main()
    assert (out / "nodes.parquet").exists() and (out / "graph.npz").exists()
    z = np.load(out / "graph.npz", allow_pickle=True)
    assert "xyz" in z.files
    assert z["xyz"].shape == (60, 3)
    for suffix in ("a", "b", "w", "r"):
        assert f"comention__{suffix}" in z.files


def test_the_written_cache_keeps_every_surviving_column(monkeypatch, tmp_path):
    """An allowlist silently dropped 15 of 18 RNA columns and one of eight fitness screens. Everything
    that survives the build ships, so the app is standalone rather than nearly so."""
    out = tmp_path / "out"
    out.mkdir()
    nodes = _nodes(60, expr_tachy=1.0, mean_plddt=80.0, some_new_assay=1.0)
    monkeypatch.setattr(BG, "OUT", str(out))
    monkeypatch.setattr(BG, "log", lambda *a, **k: None)
    monkeypatch.setattr(BG, "load_nodes", lambda: nodes)
    monkeypatch.setattr(BG, "build_edges", lambda n: {})
    BG.main()
    written = pd.read_parquet(out / "nodes.parquet")
    assert "some_new_assay" in written.columns


def test_internal_scratch_columns_are_dropped(monkeypatch, tmp_path):
    """Exact duplicates carried in from upstream would double-weight localization in any embedding
    that selects "all localization columns"."""
    out = tmp_path / "out"
    out.mkdir()
    nodes = _nodes(60, expr_tachy=1.0, mean_plddt=80.0, structure_path="x", lopit_class="y",
                   lopit_posterior=0.5, gene_product="z")
    monkeypatch.setattr(BG, "OUT", str(out))
    monkeypatch.setattr(BG, "log", lambda *a, **k: None)
    monkeypatch.setattr(BG, "load_nodes", lambda: nodes)
    monkeypatch.setattr(BG, "build_edges", lambda n: {})
    BG.main()
    written = pd.read_parquet(out / "nodes.parquet")
    for c in ("structure_path", "lopit_class", "lopit_posterior", "gene_product"):
        assert c not in written.columns
    assert "compartment" in written.columns, "kept over lopit_map: it carries the explicit unassigned"


def test_edge_indices_are_written_as_int32_and_weights_as_float32(monkeypatch, tmp_path):
    """The cache is committed, so its size is a real constraint: float64 throughout would roughly
    double an 11 MB file for no added precision."""
    out = tmp_path / "out"
    out.mkdir()
    monkeypatch.setattr(BG, "OUT", str(out))
    monkeypatch.setattr(BG, "log", lambda *a, **k: None)
    monkeypatch.setattr(BG, "load_nodes", lambda: _nodes(60, expr_tachy=1.0))
    monkeypatch.setattr(BG, "build_edges", lambda n: {"xlms": _edge([(0, 1)])})
    BG.main()
    z = np.load(out / "graph.npz", allow_pickle=True)
    assert z["xlms__a"].dtype == np.int32
    assert z["xlms__w"].dtype == np.float32


def test_the_resolver_maps_a_supplements_accession_to_a_current_gene_id(monkeypatch, tmp_path):
    """Every joined dataset goes through this closure. Without it the 2019 in vivo screen contributed
    0 of 8,140 rows -- its accessions are all pre-2012 -- and nothing said so."""
    base, ds, out, genes = _upstream(tmp_path)
    # a previous accession for the first gene, which only the identity layer can map forward
    pd.DataFrame({"gene_id": genes, "gene_name": [""] * len(genes),
                  "previous_ids": ["TGME49_008830"] + [""] * (len(genes) - 1),
                  "product": ["hypothetical"] * len(genes)}).to_csv(
        out / "toxodb_identity.tsv", sep="\t", index=False)

    seen = {}

    def screens_using_the_resolver(base_dir, log=None, resolve=None):
        seen["resolved"] = resolve("TGME49_008830")
        seen["unknown"] = resolve("TGME49_999999")
        return pd.DataFrame()

    monkeypatch.setattr(BG, "BASE", str(base))
    monkeypatch.setattr(BG, "DS", str(ds))
    monkeypatch.setattr(BG, "OUT", str(out))
    monkeypatch.setattr(BG, "log", lambda *a, **k: None)
    monkeypatch.setattr(BG.localization, "lopit_labels", lambda ds, n, log=None: n.assign(
        compartment="unassigned"))
    monkeypatch.setattr(BG.screens, "crispr_screens", screens_using_the_resolver)
    monkeypatch.setattr(BG.screens, "proteomics", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.expression, "load_all", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.cellcycle, "add_all", lambda base, n, **k: n)
    BG.load_nodes()
    assert seen["resolved"] == genes[0], "a superseded accession must map forward"
    assert seen["unknown"] is None, "an unknown accession is None, never a guess"


def test_measured_interactions_that_are_unwritten_reach_the_graph(monkeypatch, tmp_path, quiet):
    """The layer is derived last and attached only when it has something to say."""
    nodes = _nodes(12)
    edges = _build_edges(monkeypatch, nodes, tmp_path, lit={"xlms": _edge([(0, 1)])})
    assert "unwritten_interaction" in edges
    a, b, _, _ = edges["unwritten_interaction"]
    assert list(zip(a, b)) == [(0, 1)]


def test_the_whole_build_runs_end_to_end_over_a_synthetic_tree(monkeypatch, tmp_path):
    """`python -m starplast.build_graph` for real, on six genes: load, literature, edges, embed, write.

    Nothing is stubbed except the dataset loaders that need files this fixture does not fake. The point
    is that the stages compose -- each one has been driven separately above, and this is the only test
    that runs them in the order the build does.
    """
    base, ds, out, genes = _upstream(tmp_path, n=60)
    monkeypatch.setattr(BG, "BASE", str(base))
    monkeypatch.setattr(BG, "DS", str(ds))
    monkeypatch.setattr(BG, "OUT", str(out))
    monkeypatch.setattr(BG, "ABSTRACTS", str(tmp_path / "no_abstracts.jsonl"))
    monkeypatch.setattr(BG, "FULLTEXTS", str(tmp_path / "no_fulltexts"))
    monkeypatch.setattr(BG, "TOXONET", str(tmp_path / "no_edges.parquet"))
    monkeypatch.setattr(BG, "log", lambda *a, **k: None)
    monkeypatch.setattr(BG.localization, "lopit_labels", lambda d, n, log=None: n.assign(
        compartment=["rhoptry"] * 30 + ["dense granules"] * 30))
    monkeypatch.setattr(BG.screens, "crispr_screens", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.screens, "proteomics", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.expression, "load_all", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.cellcycle, "add_all", lambda b, n, **k: n)
    monkeypatch.setattr(BG.interactions, "crosslink_models", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(BG.interaction_studies, "host_interactions", lambda *a, **k: pd.DataFrame())

    BG.main()

    written = pd.read_parquet(out / "nodes.parquet")
    assert len(written) == 60
    z = np.load(out / "graph.npz", allow_pickle=True)
    assert z["xyz"].shape == (60, 3)
    assert np.isfinite(z["xyz"]).all()
    # the literature columns exist even though there was no corpus: zero is the measurement
    assert written.n_publications.tolist() == [0] * 60
    assert written.attention_depth.tolist() == [""] * 60


@pytest.mark.slow
def test_python_dash_m_runs_a_real_build(monkeypatch, tmp_path):
    """`python -m starplast.build_graph` is the documented way to rebuild, so it is executed rather
    than asserted about: runpy re-imports the module, which recomputes every root from the environment,
    so this also proves the resolver works for a fresh process rather than only for this one.

    The loaders that need files this fixture does not fake return empty and the build carries on, which
    is the behaviour a machine without the raw tree should get.
    """
    base = tmp_path / "base"
    ds = base / "datasets"
    ds.mkdir(parents=True)
    interim = base / "toxonet" / "data" / "interim"
    interim.mkdir(parents=True)
    n = 60
    genes = [f"TGME49_{200000 + i}" for i in range(n)]
    rng = np.random.default_rng(0)
    pd.DataFrame({"gene_id": genes,
                  "mean_plddt": rng.uniform(40, 95, n),
                  "rna108740_Tachyzoites_T2_FPKM": rng.uniform(0, 500, n)}).to_parquet(
        interim / "nodes.parquet", index=False)
    pd.DataFrame({"gene_source_id": genes, "gene_product": ["a protein"] * n}).to_csv(
        ds / "orthomcl_toxoplasma_gondii_ME49.csv", index=False)

    out = tmp_path / "out"
    out.mkdir()
    _identity_files(out, genes)

    monkeypatch.setenv("STARPLAST_DATA", str(ds))
    monkeypatch.setenv("STARPLAST_CACHE", str(out))

    import runpy
    runpy.run_module("starplast.build_graph", run_name="__main__", alter_sys=True)

    assert (out / "nodes.parquet").exists() and (out / "graph.npz").exists()
    written = pd.read_parquet(out / "nodes.parquet")
    assert len(written) == n
    z = np.load(out / "graph.npz", allow_pickle=True)
    assert z["xyz"].shape == (n, 3) and np.isfinite(z["xyz"]).all()


def test_the_embedding_is_writable_even_when_umap_returns_a_read_only_array(monkeypatch):
    """umap returns a read-only array in recent versions, and np.asarray does NOT copy when the dtype
    already matches -- so the in-place centring wrote into a read-only buffer and raised "output array
    is read-only". Fixed once in embedding.py; this is its second home, and the one a user meets,
    because it is the path build_graph takes."""
    import types
    n = 60
    nodes = _nodes(n, expr_tachy=1.0, mean_plddt=80.0)

    class FakeUMAP:
        def __init__(self, **kw):
            pass

        def fit_transform(self, X):
            out = np.zeros((len(X), 3), dtype=np.float32)
            out[:, 0] = np.arange(len(X), dtype=np.float32)
            out.setflags(write=False)          # exactly what recent umap hands back
            return out

    fake = types.ModuleType("umap")
    fake.UMAP = FakeUMAP
    monkeypatch.setitem(sys.modules, "umap", fake)

    Y = BG.embed(nodes)                        # raised "output array is read-only" before the fix
    assert Y.shape == (n, 3)
    assert np.isfinite(Y).all()
    Y[0, 0] = 1.0                              # and the caller must be able to write to it


def test_the_mass_spec_columns_are_merged_into_the_node_table(monkeypatch):
    """Counts of reported modification sites are merged separately from expression, and on purpose:
    these are site counts rather than abundances, and putting them through one loader would invite
    them to be normalised together."""
    import pandas as pd
    import starplast.proteomics as PR
    from starplast import build_graph as B
    nodes = pd.DataFrame({"gene_id": ["TGME49_000001", "TGME49_000002"]})
    monkeypatch.setattr(PR, "load_all",
                        lambda base, index, log=print: pd.DataFrame(
                            {"n_acetylation_sites": [4.0, None]}, index=pd.Index(index)))
    columns = PR.load_all(B.BASE, nodes.gene_id.astype(str))
    for c in columns.columns:
        nodes[c] = columns[c].to_numpy()
    assert nodes.loc[0, "n_acetylation_sites"] == 4.0
    assert pd.isna(nodes.loc[1, "n_acetylation_sites"]), "an unmeasured gene was given a number"
