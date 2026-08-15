#!/usr/bin/env python3
"""Reading the findings back as claims.

The thing being tested is a piece of writing, which makes it tempting to test nothing. What is
actually checkable is everything that matters: that the ranking puts the useful finding above the
statistically-stronger useless one, that every number quoted in a sentence is the number in the row,
that a circular finding is never silently promoted, and that no report ever states a claim as fact.
"""
from __future__ import annotations

import os
import re
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import interpret as I  # noqa: E402


def nodes(n=40, papers=0):
    return pd.DataFrame({"gene_id": [f"TGME49_{i:06d}" for i in range(n)],
                         "product": ["hypothetical protein"] * n,
                         "n_publications": [papers] * n})


def guilt_row(**kw):
    row = {"kind": "guilt", "layer": "compartment", "layer_kind": "discrete", "cluster": 5,
           "category": "IMC", "n_cluster": 58, "n_known": 58, "n_hits": 41, "purity": 0.71,
           "background": 0.06, "lift": 11.8, "p": 1e-24, "q": 1e-21, "n_predicted": 12,
           "circular": False, "genes": [f"TGME49_{i:06d}" for i in range(12)]}
    row.update(kw)
    return row


def continuous_row(**kw):
    row = {"kind": "guilt", "layer": "fit_invitro_hff", "layer_kind": "continuous", "cluster": 12,
           "category": "", "n_cluster": 30, "n_known": 23, "n_hits": 23, "mean_in": -3.4,
           "mean_out": -0.31, "effect": -2.1, "lift": 2.1, "p": 1e-9, "q": 1e-7, "n_predicted": 7,
           "circular": False, "genes": [f"TGME49_{i:06d}" for i in range(7)]}
    row.update(kw)
    return row


def split_row(**kw):
    row = {"kind": "disagreement", "layer": "compartment", "layer_kind": "discrete",
           "other": "fit_invitro_hff", "other_kind": "continuous", "cluster": 9,
           "category": "rhoptry", "n_cluster": 36, "n_known": 36, "purity": 0.88, "gap": 2.6,
           "mean_low": -3.2, "mean_high": -0.1, "n_low": 14, "n_high": 22, "p": 1e-6, "q": 1e-4,
           "n_predicted": 14, "circular": False, "genes": [f"TGME49_{i:06d}" for i in range(14)]}
    row.update(kw)
    return row


def conjunction_row(**kw):
    row = {"kind": "conjunction", "layer": "compartment", "other": "stage",
           "cluster": 41, "category": "golgi", "other_category": "bradyzoite",
           "purity": 0.78, "joint_lift": 6.2, "interaction_ratio": 3.1,
           "q": 2e-9, "n_predicted": 9, "n_known": 32, "circular": False,
           "genes": [f"TGME49_{i:06d}" for i in range(9)]}
    row.update(kw)
    return row


# --------------------------------------------------------------------------- ranking
def test_a_useful_finding_outranks_a_stronger_one_about_nothing():
    """The whole job. Statistical strength alone puts the giant obvious cluster first every time,
    and the reader already knows about it."""
    table = pd.DataFrame([
        guilt_row(cluster=1, q=1e-40, n_predicted=2),      # overwhelming, and about two genes
        guilt_row(cluster=2, q=1e-8, n_predicted=25),      # solid, and about twenty-five
    ])
    ranked = I.interest(table, nodes())
    assert ranked.iloc[0].cluster == 2, ranked[["cluster", "interest", "strength", "reach"]]


def test_a_prediction_about_well_studied_genes_ranks_below_one_about_unknowns():
    """Predicting the localisation of a gene with two hundred papers is not a prediction."""
    famous = nodes(papers=200)
    unknown = nodes(papers=0)
    assert I.novelty(unknown, unknown.gene_id[:5]) > I.novelty(famous, famous.gene_id[:5])
    hot = I.interest(pd.DataFrame([guilt_row()]), famous).iloc[0]
    cold = I.interest(pd.DataFrame([guilt_row()]), unknown).iloc[0]
    assert cold.interest > hot.interest


def test_a_circular_finding_scores_zero_however_strong_it_is():
    table = pd.DataFrame([guilt_row(q=1e-60, n_predicted=90, circular=True), guilt_row()])
    ranked = I.interest(table, nodes())
    assert ranked.iloc[0].circular == False        # noqa: E712 -- numpy bool, `is False` fails
    assert float(ranked[ranked.circular].iloc[0].interest) == 0.0


def test_novelty_is_neither_rewarded_nor_punished_with_no_literature_layer():
    bare = nodes().drop(columns=["n_publications"])
    assert I.novelty(bare, ["TGME49_000001"]) == 0.5
    assert I.novelty(nodes(), []) == 0.5
    assert I.novelty(nodes(), ["not_a_gene"]) == 0.5


def test_an_empty_table_ranks_to_nothing():
    assert I.interest(pd.DataFrame(), nodes()).empty
    assert I.interest(None, nodes()).empty


# --------------------------------------------------------------------------- the sentences
def test_the_discrete_claim_says_what_was_seen_and_what_is_predicted():
    text = I.sentence(pd.Series(guilt_row()), nodes())
    assert "41 of 58" in text and "IMC" in text
    assert "71%" in text and "6.0%" in text and "11.8x" in text
    assert "12 genes" in text and "TGME49_000000" in text
    assert "predicts" in text.lower()


def test_the_continuous_claim_gives_both_means_and_the_direction():
    text = I.sentence(pd.Series(continuous_row()), nodes())
    assert "-3.40" in text and "-0.31" in text
    assert "lower" in text and "d = -2.1" in text
    assert "7 genes" in text and "never measured" in text


def test_the_split_claim_names_both_groups_and_the_gap():
    text = I.sentence(pd.Series(split_row()), nodes())
    assert "88% rhoptry" in text and "fit_invitro_hff" in text
    assert "-3.20" in text and "-0.10" in text and "2.6 pooled SD" in text
    assert "14" in text and "22" in text
    assert "subdivision" in text


def test_a_categorical_disagreement_lists_the_groups_it_split_into():
    row = split_row(other="cellcycle_phase", other_kind="discrete",
                    groups=["G1", "S/M"], counts=[19, 14])
    text = I.sentence(pd.Series(row), nodes())
    assert "G1 (19)" in text and "S/M (14)" in text
    assert "disagrees" in text


def test_a_conjunction_names_both_layers_without_stating_causation():
    text = I.sentence(pd.Series(conjunction_row()), nodes())
    assert "golgi" in text and "bradyzoite" in text and "3.1x" in text
    assert "hypothesis" in text and "proof" in text


def test_a_product_name_is_used_where_there_is_one():
    table = nodes()
    table.loc[0, "product"] = "inner membrane complex protein IMC1"
    text = I.sentence(pd.Series(guilt_row()), table)
    assert "IMC1" in text
    assert "hypothetical protein" not in text, "a placeholder was quoted as a name"


def test_a_long_list_of_genes_is_summarised_rather_than_dumped():
    row = guilt_row(n_predicted=60, genes=[f"TGME49_{i:06d}" for i in range(60)])
    text = I.sentence(pd.Series(row), nodes(60))
    assert "and 54 more" in text


def test_gene_ids_survive_a_table_with_no_gene_id_column():
    bare = nodes().drop(columns=["gene_id"])
    assert "TGME49_000001" in I._names(bare, ["TGME49_000001"])


# --------------------------------------------------------------------------- caveats and report
def test_every_finding_carries_the_reason_to_doubt_it():
    """Never empty by design: a report of claims with no caveats reads as a report of results."""
    for row in (guilt_row(), continuous_row(), split_row()):
        out = I.caveats(pd.Series(row), nodes())
        assert out and any("hypothesis" in c for c in out)


def test_the_caveats_name_the_specific_problem_where_there_is_one():
    assert any("CIRCULAR" in c for c in I.caveats(pd.Series(guilt_row(circular=True)), nodes()))
    assert any("well published" in c
               for c in I.caveats(pd.Series(guilt_row(novelty=0.1)), nodes()))
    assert any("few points" in c for c in I.caveats(pd.Series(guilt_row(n_known=8)), nodes()))
    assert any("chance" in c for c in I.caveats(pd.Series(guilt_row(q=0.04)), nodes()))


def test_the_report_leads_with_how_much_survived_correction():
    table = pd.DataFrame([guilt_row(), continuous_row(), split_row(),
                          guilt_row(cluster=77, q=0.9)])
    text = I.report(table, nodes())
    assert "3 claims survive correction" in text
    assert "guilt by association**: 2" in text and "layer disagreement**: 1" in text
    assert text.index("### 1.") < text.index("### 2."), "the report is not in ranked order"


def test_the_report_says_so_when_there_is_nothing_to_report():
    text = I.report(pd.DataFrame(), nodes())
    assert "Nothing to report" in text
    assert "real answer" in text, "an empty result was reported as a failure"


def test_circular_findings_are_shown_separately_rather_than_dropped():
    """A reader who sees only the clean findings cannot tell whether the run was clean."""
    text = I.report(pd.DataFrame([guilt_row(), guilt_row(cluster=3, circular=True)]), nodes())
    assert "Excluded as circular" in text and "compartment" in text
    assert "1 claims survive" in text


def test_no_report_ever_states_a_claim_as_a_fact():
    """The wording is load-bearing. A report that overstates once is a report nobody can use
    unread, and these are hypotheses about genes that behave alike."""
    text = I.report(pd.DataFrame([guilt_row(), continuous_row(), split_row()]), nodes()).lower()
    for word in (" proves", " demonstrates that", " shows that these genes are"):
        assert word not in text, word
    assert "hypothesis" in text and "predict" in text


def test_the_report_is_markdown_a_notebook_can_hold():
    text = I.report(pd.DataFrame([guilt_row()]), nodes())
    assert text.startswith("## ")
    assert re.search(r"\*interest \d\.\d+ = strength", text), "the ranking is not shown inline"


def test_a_cluster_enriched_for_two_things_says_they_are_alternatives():
    """Seen on the real map the first time this ran: cluster 29 came back 30% ER and 26% golgi, and
    the same 113 unlabelled genes were offered to both. A gene is in one compartment."""
    table = pd.DataFrame([guilt_row(category="ER"), guilt_row(category="golgi")])
    ranked = I.interest(table, nodes())
    out = I.caveats(ranked.iloc[0], nodes(), ranked)
    assert any("alternatives" in c and "golgi" in c for c in out)
    assert not any("alternatives" in c for c in I.caveats(ranked.iloc[0], nodes()))
    assert "alternatives" in I.report(table, nodes())
