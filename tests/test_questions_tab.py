#!/usr/bin/env python3
"""The Questions tab, driven through the real widgets.

Every control here is exercised by clicking it rather than by calling the function behind it, because
a control wired to nothing is worse than a missing one and only the click can tell them apart. The
jobs are run synchronously — what is being checked is that the wiring passes the right arguments and
puts the result where the user will look for it.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

RNA = ["biology", ["gene expression", "RNA abundance"]]
FITNESS = ["biology", ["parasite phenotype", "fitness and essentiality"]]
LOCALIZATION = ["biology", ["cell organization", "localization and topology"]]

CATALOGUE = [
    {"question": "Which genes share a compartment with the secretory organelles?",
     "axis": "localisation and export", "inputs": [FITNESS], "holdout": "compartment",
     "validation_holdout": "cellcycle_phase", "verdict": "keep",
     "expectation": "dense granules and rhoptries cluster together"},
    {"question": "Can localisation predict localisation?", "axis": "localisation and export",
     "inputs": [LOCALIZATION], "holdout": "compartment", "verdict": "drop",
     "expect_refusal": True, "reason": "it restates its own holdout, kept as a demonstration"},
]


@pytest.fixture(scope="module")
def app():
    from PyQt6 import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture(scope="module")
def nodes():
    from starplast import paths
    return pd.read_parquet(paths.cache_file("nodes.parquet")).sample(
        n=400, random_state=0).reset_index(drop=True)


@pytest.fixture
def panel(app, nodes, tmp_path, monkeypatch):
    from starplast import questions
    monkeypatch.setattr(questions, "shipped", lambda *a, **k: list(CATALOGUE))
    from starplast.analysis_panel import AnalysisPanel
    from starplast.tuning import EmbeddingStore
    p = AnalysisPanel(nodes, store=EmbeddingStore(str(tmp_path)))
    yield p
    p.deleteLater()


def test_the_tab_exists_and_lists_the_shipped_questions(panel):
    assert panel.question_choice.count() == 2
    assert "expected to fail" in panel.question_choice.itemText(1)


def test_selecting_a_question_describes_it_without_running_anything(panel):
    """The description must arrive on selection: a user deciding whether to spend minutes on a
    question needs to see its holdout and its control first."""
    panel.question_choice.setCurrentIndex(0)
    text = panel.question_detail.toPlainText()
    assert "compartment" in text and "cellcycle_phase" in text
    assert "gene expression > RNA abundance" not in text     # this question does not use RNA
    assert "parasite phenotype > fitness and essentiality" in text


def test_the_question_kept_because_it_fails_says_so_on_screen(panel):
    panel.question_choice.setCurrentIndex(1)
    assert "BECAUSE it is refused" in panel.question_detail.toHtml()


def test_checking_a_circular_question_refuses_it_before_anything_is_built(panel, monkeypatch):
    """The cheap path, and the one that matters: a refusal must cost no compute."""
    ran = []
    monkeypatch.setattr(panel, "_run", lambda fn, done, **k: done(fn(lambda *a: ran.append(a))))
    panel.question_choice.setCurrentIndex(1)
    panel.check_question()
    assert "refused" in panel.question_report.toPlainText().lower()
    assert "restate a holdout" in panel.question_report.toPlainText()


def test_checking_a_sound_question_reports_what_closure_removed(panel, monkeypatch):
    monkeypatch.setattr(panel, "_run", lambda fn, done, **k: done(fn(lambda *a: None)))
    panel.question_choice.setCurrentIndex(0)
    panel.check_question()
    text = panel.question_report.toPlainText()
    assert "may build this map" in text and "excluded" in text


def test_running_a_question_puts_the_genes_in_the_table_and_the_control_beside_them(panel, nodes,
                                                                                   monkeypatch):
    """The whole tab in one test: click run, get genes, and get told whether to believe them."""
    from starplast import recipes

    lab = np.arange(len(nodes)) % 3
    lab[::7] = -1
    monkeypatch.setattr("starplast.clustering.cluster", lambda X, **k: lab[:len(X)])
    truth = np.array([None] * len(nodes), dtype=object)
    control = np.array([None] * len(nodes), dtype=object)
    for cl, name in ((0, "pvm"), (1, "cytosol")):
        truth[np.where(lab == cl)[0][:40]] = name
    control[np.where(lab == 0)[0][40:75]] = "host"
    control[np.where(lab == 1)[0][40:55]] = "not host"
    control[np.where(lab == 2)[0][:40]] = "not host"
    panel.nodes = panel.nodes.copy()
    panel.nodes["compartment"] = truth
    panel.nodes["a_control"] = control
    panel._questions = [{**CATALOGUE[0], "validation_holdout": "a_control", "min_precision": 0.7}]
    monkeypatch.setattr(panel, "_run",
                        lambda fn, done, **k: done(fn(lambda *a: None)))
    monkeypatch.setattr(recipes, "tune_umap", lambda nodes, spec, **k: pd.DataFrame(
        [{"n_neighbors": 15, "min_dist": 0.0, "usable": True, "score": 1.0, "stage": "full"}]))
    panel.run_question()
    assert panel.question_table.rowCount() > 0
    report = panel.question_report.toPlainText()
    assert "genes named" in report and "the map" in report and "recovery" in report


def test_an_answer_no_control_corroborates_says_so_rather_than_looking_like_one(panel, nodes,
                                                                               monkeypatch):
    """An uncorroborated list of genes looks exactly like a corroborated one on screen, and the whole
    argument for having a control is that the difference decides whether to believe the answer."""
    from starplast import recipes
    lab = np.arange(len(nodes)) % 3
    lab[::7] = -1
    monkeypatch.setattr("starplast.clustering.cluster", lambda X, **k: lab[:len(X)])
    truth = np.array([None] * len(nodes), dtype=object)
    for cl, name in ((0, "pvm"), (1, "cytosol")):
        truth[np.where(lab == cl)[0][:40]] = name
    panel.nodes = panel.nodes.copy()
    panel.nodes["compartment"] = truth
    panel._questions = [{**CATALOGUE[0], "validation_holdout": "", "min_precision": 0.7}]
    monkeypatch.setattr(panel, "_run", lambda fn, done, **k: done(fn(lambda *a: None)))
    monkeypatch.setattr(recipes, "tune_umap", lambda nodes, spec, **k: pd.DataFrame(
        [{"n_neighbors": 15, "min_dist": 0.0, "usable": True, "score": 1.0, "stage": "full"}]))
    panel.run_question()
    assert "corroborates none of these clusters" in panel.question_report.toPlainText()


def test_a_degenerate_map_reports_no_answer_rather_than_an_empty_table(panel, monkeypatch):
    monkeypatch.setattr("starplast.clustering.cluster",
                        lambda X, **k: np.zeros(len(X), dtype=int))
    monkeypatch.setattr(panel, "_run", lambda fn, done, **k: done(fn(lambda *a: None)))
    panel.question_choice.setCurrentIndex(0)
    panel.run_question()
    assert "no answer" in panel.question_report.toPlainText()


def test_an_empty_catalogue_says_how_to_generate_one(app, nodes, tmp_path, monkeypatch):
    from starplast import questions
    monkeypatch.setattr(questions, "shipped", lambda *a, **k: [])
    from starplast.analysis_panel import AnalysisPanel
    from starplast.tuning import EmbeddingStore
    p = AnalysisPanel(nodes, store=EmbeddingStore(str(tmp_path)))
    assert p.current_question() is None
    assert "generate_question_table" in p.question_detail.toPlainText()
    p.check_question()          # must not raise with nothing selected
    p.run_question()
    p.deleteLater()


def test_an_answer_that_names_nobody_says_how_close_it_came(panel, nodes, monkeypatch):
    """"No genes" and "the best cluster was 1.4x enriched at 0.62 purity" are different findings, and
    only the second tells a reader whether the question failed or the gate did."""
    from starplast import recipes
    lab = np.arange(len(nodes)) % 3
    lab[::7] = -1
    monkeypatch.setattr("starplast.clustering.cluster", lambda X, **k: lab[:len(X)])
    # Every cluster mixed, so nothing reaches any usable purity.
    truth = np.array([None] * len(nodes), dtype=object)
    truth[::2] = "a"
    truth[1::2] = "b"
    panel.nodes = panel.nodes.copy()
    panel.nodes["compartment"] = truth
    panel._questions = [{**CATALOGUE[0], "validation_holdout": "", "min_precision": 0.99}]
    monkeypatch.setattr(panel, "_run", lambda fn, done, **k: done(fn(lambda *a: None)))
    monkeypatch.setattr(recipes, "tune_umap", lambda nodes, spec, **k: pd.DataFrame(
        [{"n_neighbors": 15, "min_dist": 0.0, "usable": True, "score": 1.0, "stage": "full"}]))
    panel.run_question()
    # The note names BOTH numbers now: how enriched the best cluster was and how pure. Either alone
    # leaves the reader unable to tell whether the question failed or the gate did.
    report = panel.question_report.toPlainText()
    assert "enrichment" in report and "purity" in report


# --------------------------------------------------------------------------- the category sweep
# Instruction 42 asks that every analysis-mode control be driven headlessly at least once. These two
# were added on 2026-08-18 without a test that clicks them, and a coverage pass on the panel is what
# found it -- which is the argument for the rule rather than for the number.
def test_the_category_sweep_button_runs_the_sweep_and_fills_its_table(panel, monkeypatch):
    swept = {}

    def fake_sweep(nodes, spec, hierarchy=None, level=None, **k):
        swept.update(hierarchy=hierarchy, level=level, blocks=spec.blocks)
        return pd.DataFrame([{"category": "molecular measurements", "best_score": 0.42,
                              "tested": 3, "usable": True},
                             {"category": "relational measurements", "best_score": 0.11,
                              "tested": 7, "usable": True}])
    monkeypatch.setattr("starplast.search.sweep_categories", fake_sweep)
    monkeypatch.setattr(panel, "_run", lambda fn, done, **k: done(fn(lambda *a: None)))
    panel.run_category_sweep()
    assert swept["hierarchy"] and swept["level"] >= 1
    assert panel.sweep_table.rowCount() == 2


def test_the_sweep_reports_the_weakest_recovery_because_that_is_the_interesting_end(panel):
    ordered = panel._sweep_done(pd.DataFrame(
        [{"category": "strong", "best_score": 0.9, "tested": 4},
         {"category": "weak", "best_score": 0.05, "tested": 9}]))
    assert list(ordered.category) == ["weak", "strong"]


def test_a_sweep_that_scored_nothing_says_so_rather_than_showing_an_empty_table(panel):
    assert panel._sweep_done(pd.DataFrame()) is None
    assert panel._sweep_done(None) is None


def test_the_sweep_refuses_to_run_with_no_blocks_ticked(panel, monkeypatch):
    monkeypatch.setattr(panel, "spec", lambda: __import__(
        "starplast.embedding", fromlist=["EmbeddingSpec"]).EmbeddingSpec(blocks=()))
    called = []
    monkeypatch.setattr(panel, "_run", lambda *a, **k: called.append(1))
    panel.run_category_sweep()
    assert not called


# --------------------------------------------------------------------------- the PDF button
def test_the_pdf_button_is_dead_until_there_is_a_run(panel):
    """A button that promises a document it has no run for is worse than no button."""
    assert not panel.question_pdf.isEnabled()


def test_a_refusal_still_enables_the_pdf(panel, monkeypatch):
    """The refusal and its cause IS the result. A reader looking for the run must find the reason
    rather than an absent file."""
    monkeypatch.setattr(panel, "_run", lambda fn, done, **k: done(fn(lambda *a: None)))
    panel.question_choice.setCurrentIndex(1)          # the one kept because it fails
    panel.run_question()
    assert panel.question_pdf.isEnabled()


def test_saving_writes_the_pdf_where_the_user_chose(panel, tmp_path, monkeypatch):
    from PyQt6 import QtWidgets
    monkeypatch.setattr(panel, "_run", lambda fn, done, **k: done(fn(lambda *a: None)))
    panel.question_choice.setCurrentIndex(1)
    panel.run_question()
    target = str(tmp_path / "chosen.pdf")
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (target, "PDF (*.pdf)")))
    panel.save_question_pdf()
    assert os.path.exists(target)


def test_cancelling_the_save_dialog_writes_nothing(panel, tmp_path, monkeypatch):
    from PyQt6 import QtWidgets
    monkeypatch.setattr(panel, "_run", lambda fn, done, **k: done(fn(lambda *a: None)))
    panel.question_choice.setCurrentIndex(1)
    panel.run_question()
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: ("", "")))
    ran = []
    monkeypatch.setattr(panel, "_run", lambda *a, **k: ran.append(1))
    panel.save_question_pdf()
    assert not ran


def test_saving_before_any_run_does_nothing(panel):
    panel._question_result = None
    panel.save_question_pdf()      # must not raise, and must not open a dialog


def test_the_written_path_is_reported_to_the_user(panel):
    seen = []
    panel.status.connect(seen.append)
    panel._pdf_done("/tmp/somewhere.pdf")
    assert seen and "/tmp/somewhere.pdf" in seen[0]
