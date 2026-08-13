#!/usr/bin/env python3
"""The annotation store, which is mostly a set of refusals.

Writing down "this cluster is mostly dense granules, so these unlabelled members probably are too"
is a few lines of code and the easiest way for this project to start manufacturing exactly what it
exists to prevent: a candidate list looks identical whether it is 90% right or 6%, and measured on
the coarse components of this map it was 1-6%. So these tests are about what the store will not do
-- save without a validated precision, save at a precision that makes the row a coin toss, or write
anything at all into the node table.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import annotations as A  # noqa: E402


def _annotation(**kw):
    base = dict(gene_id="TGME49_208830", proposed="dense granules", target="compartment",
                cluster=3, blocks="fitness_screens", n_neighbors=15, min_dist=0.1,
                min_cluster_size=25, seed=42, cluster_size=140, cluster_frac_category=0.62,
                cluster_frac_contradicting=0.08, enrichment=6.4, precision=0.71, recall=0.44,
                folds=5, refit=False, date="2026-08-12", reasoning="dense and mostly GRAs")
    base.update(kw)
    return A.Annotation(**base)


# --------------------------------------------------------------------------- the refusals
def test_a_candidate_with_no_validated_precision_cannot_be_saved(tmp_path):
    """The one that matters. A candidate list with no error rate is a list of guesses, and stored it
    becomes a list of guesses that looks like a record."""
    store = A.AnnotationStore(str(tmp_path / "a.csv"))
    with pytest.raises(ValueError, match="no validated precision"):
        store.save(_annotation(precision=float("nan")))
    assert not os.path.exists(store.path), "a refused save still wrote a file"


def test_a_precision_that_makes_it_a_coin_toss_is_refused(tmp_path):
    """Measured on this map's coarse components, annotating from them would be right 1-6% of the
    time. A store that accepted those rows would be a machine for producing them."""
    store = A.AnnotationStore(str(tmp_path / "a.csv"))
    with pytest.raises(ValueError, match="coin toss"):
        store.save(_annotation(precision=0.04))


def test_the_refusal_says_what_to_do_rather_than_only_no(tmp_path):
    msg = _annotation(precision=float("nan")).check()
    assert "Validation tab" in msg and "compartment" in msg


def test_an_annotation_needs_a_gene_and_a_label(tmp_path):
    assert "needs a gene" in _annotation(gene_id="").check()
    assert "needs a gene" in _annotation(proposed="").check()


def test_one_unjustified_row_refuses_the_whole_save(tmp_path):
    """A partial save that dropped the bad rows silently would leave the user believing they had
    saved what they were looking at."""
    store = A.AnnotationStore(str(tmp_path / "a.csv"))
    with pytest.raises(ValueError):
        store.save([_annotation(), _annotation(gene_id="TGME49_000001", precision=float("nan"))])
    assert store.load().empty


# --------------------------------------------------------------------------- what it does save
def test_an_annotation_round_trips_with_every_number_that_justifies_it(tmp_path):
    store = A.AnnotationStore(str(tmp_path / "a.csv"))
    store.save(_annotation())
    got = store.load()
    assert len(got) == 1
    r = got.iloc[0]
    assert r.gene_id == "TGME49_208830" and r.proposed == "dense granules"
    assert r.precision == pytest.approx(0.71) and r.recall == pytest.approx(0.44)
    assert r.cluster_frac_category == pytest.approx(0.62) and r.blocks == "fitness_screens"
    assert r.reasoning == "dense and mostly GRAs" and r.date == "2026-08-12"
    assert bool(r.refit) is False, "whether it re-embedded per fold must survive the round trip"


def test_the_file_is_readable_by_a_person(tmp_path):
    """It is meant to be shared and argued with, so the gene and the label come first and the
    reasoning last."""
    store = A.AnnotationStore(str(tmp_path / "a.csv"))
    store.save(_annotation())
    head = open(store.path).readline().strip().split(",")
    assert head[:2] == ["gene_id", "proposed"] and head[-1] == "reasoning"


def test_annotating_the_same_gene_again_replaces_rather_than_duplicates(tmp_path):
    """The second opinion is the one someone just formed. Two rows for one gene is a file that
    cannot answer the question it exists to answer."""
    store = A.AnnotationStore(str(tmp_path / "a.csv"))
    store.save(_annotation(proposed="dense granules"))
    store.save(_annotation(proposed="rhoptries", precision=0.55))
    got = store.load()
    assert len(got) == 1 and got.iloc[0].proposed == "rhoptries"


def test_the_same_gene_under_a_different_target_is_a_different_annotation(tmp_path):
    """Its compartment and its cell-cycle phase are two proposals, not one changed twice."""
    store = A.AnnotationStore(str(tmp_path / "a.csv"))
    store.save(_annotation(target="compartment"))
    store.save(_annotation(target="cellcycle_phase", proposed="S"))
    assert len(store.load()) == 2


def test_a_proposal_can_be_withdrawn(tmp_path):
    """Withdrawing has to be as easy as proposing, or the file fills with things nobody believes."""
    store = A.AnnotationStore(str(tmp_path / "a.csv"))
    store.save([_annotation(), _annotation(gene_id="TGME49_000002")])
    store.remove(["TGME49_000002"])
    assert list(store.load().gene_id) == ["TGME49_208830"]
    store.remove(["TGME49_208830"], target="cellcycle_phase")
    assert len(store.load()) == 1, "removing under the wrong target must not remove anything"
    store.remove(["TGME49_208830"], target="compartment")
    assert store.load().empty


def test_removing_from_an_empty_store_is_not_an_error(tmp_path):
    assert A.AnnotationStore(str(tmp_path / "none.csv")).remove(["x"]).empty


def test_replacing_the_whole_file_is_possible_but_not_the_default(tmp_path):
    store = A.AnnotationStore(str(tmp_path / "a.csv"))
    store.save(_annotation())
    store.save(_annotation(gene_id="TGME49_000003"), replace=True)
    assert list(store.load().gene_id) == ["TGME49_000003"]


def test_an_empty_store_has_the_columns_a_full_one_has(tmp_path):
    empty = A.AnnotationStore(str(tmp_path / "none.csv")).load()
    assert list(empty.columns) == A.COLUMNS and empty.empty


def test_a_file_missing_a_column_still_loads(tmp_path):
    """These files are hand-edited. One saved by an older version, or one a collaborator trimmed in
    a spreadsheet, has to open rather than take the tab down."""
    path = tmp_path / "old.csv"
    pd.DataFrame({"gene_id": ["g1"], "proposed": ["x"]}).to_csv(path, index=False)
    got = A.AnnotationStore(str(path)).load()
    assert set(A.COLUMNS) <= set(got.columns) and len(got) == 1


def test_the_mask_says_which_genes_carry_an_annotation(tmp_path):
    """What the map's fourth colour is drawn from, defined once so "is this gene annotated" has one
    answer."""
    store = A.AnnotationStore(str(tmp_path / "a.csv"))
    store.save(_annotation(gene_id="g2"))
    m = store.mask(["g1", "g2", "g3"])
    assert list(m) == [False, True, False]
    assert not A.AnnotationStore(str(tmp_path / "none.csv")).mask(["g1"]).any()


# --------------------------------------------------------------------------- joining the numbers
def _candidates():
    return pd.DataFrame({"gene_id": ["g1", "g2"], "proposed": ["dense granules"] * 2,
                         "cluster": [3, 3], "cluster_size": [140, 140],
                         "cluster_frac_category": [0.62, 0.62],
                         "cluster_frac_contradicting": [0.08, 0.08], "enrichment": [6.4, 6.4]})


class _Scores:
    precision, recall, n_folds, refit = 0.71, 0.44, 5, True


def test_candidates_arrive_carrying_the_validation_numbers(tmp_path):
    """The join happens once, in one place: a candidate list and a precision that came from a
    different category, or a different clustering, would be worse than no number at all."""
    ann = A.from_candidates(_candidates(), "compartment", _Scores(),
                            configuration={"blocks": "fitness_screens", "n_neighbors": 15,
                                           "min_dist": 0.1, "min_cluster_size": 25, "seed": 42},
                            reasoning="looks clean", date="2026-08-12")
    assert len(ann) == 2
    assert all(a.precision == 0.71 and a.recall == 0.44 and a.refit for a in ann)
    assert all(a.blocks == "fitness_screens" and a.n_neighbors == 15 for a in ann)
    assert all(a.check() == "" for a in ann)


def test_candidates_from_an_unvalidated_run_refuse_themselves(tmp_path):
    class NoScores:
        pass

    ann = A.from_candidates(_candidates(), "compartment", NoScores())
    assert ann and all("no validated precision" in a.check() for a in ann)


def test_an_annotation_is_never_written_into_the_node_table(tmp_path):
    """The discipline the whole application rests on. The node table is measurement; an inference
    stored beside it becomes indistinguishable from one the moment anybody reads the table without
    knowing which columns are which."""
    from starplast import paths
    nodes = pd.read_parquet(paths.cache_file("nodes.parquet"))
    store = A.AnnotationStore(str(tmp_path / "a.csv"))
    store.save(_annotation(gene_id=str(nodes.gene_id.iloc[0])))
    after = pd.read_parquet(paths.cache_file("nodes.parquet"))
    assert list(after.columns) == list(nodes.columns)
    assert after.equals(nodes)


def test_the_annotation_colour_is_used_for_nothing_else():
    """Measurement, inference and absence have one each. A proposal reading as any of them is the
    failure mode -- so it is checked against the palettes rather than assumed."""
    from starplast import theme as TH
    from starplast.app import DEPTH_COLOUR, GREY, PALETTE
    taken = [tuple(np.round(c, 3)) for c in PALETTE]
    taken += [tuple(np.round(c, 3)) for c in DEPTH_COLOUR.values()]
    taken.append(tuple(np.round(GREY, 3)))
    for theme in TH.THEMES:
        taken += [tuple(np.round(c[:3], 3)) for c in TH.categorical_colours(30, theme)]
        taken.append(tuple(np.round(TH.unknown_colour(theme)[:3], 3)))
    mine = np.asarray(A.ANNOTATION_COLOUR)
    assert all(np.linalg.norm(mine - np.asarray(c)) > 0.15 for c in taken), \
        "the annotation colour is too close to one that already means something else"
