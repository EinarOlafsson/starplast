#!/usr/bin/env python3
"""Building the same map twice and saying what differs.

The clock is the smaller half. cuml's UMAP is a different implementation of a stochastic algorithm,
so turning the switch on does not make a map faster -- it makes a DIFFERENT map of the same data.
Measured on the shipped cache at 2,000 genes: 27x faster, and 71% of each gene's fifteen nearest
neighbours survive. Both halves of that sentence are the answer.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import benchmark as B  # noqa: E402
from starplast.embedding import EmbeddingSpec  # noqa: E402


def _nodes(n=300, seed=0):
    rng = np.random.default_rng(seed)
    d = pd.DataFrame({f"fit_{i}": rng.normal(size=n) for i in range(6)})
    d.insert(0, "gene_id", [f"TGME49_{200000 + i}" for i in range(n)])
    return d


# --------------------------------------------------------------------------- the comparison itself
def test_an_embedding_agrees_with_itself_completely():
    """The sanity check the number is worthless without."""
    Y = np.random.default_rng(0).normal(size=(200, 3))
    assert B.knn_overlap(Y, Y) == 1.0


def test_an_unrelated_embedding_agrees_with_almost_nothing():
    rng = np.random.default_rng(1)
    assert B.knn_overlap(rng.normal(size=(200, 3)), rng.normal(size=(200, 3))) < 0.15


def test_rotation_and_scale_do_not_count_as_disagreement():
    """Two embeddings can look unlike -- rotated, mirrored, differently scaled -- while preserving
    every neighbourhood, which is why coordinates cannot be compared and neighbours can."""
    Y = np.random.default_rng(2).normal(size=(200, 3))
    theta = 0.7
    R = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    assert B.knn_overlap(Y, (Y @ R) * 3.0 + 11.0) == 1.0
    assert B.knn_overlap(Y, Y * np.array([1.0, -1.0, 1.0])) == 1.0


def test_too_few_points_is_undefined_rather_than_a_number():
    assert np.isnan(B.knn_overlap(np.zeros((2, 3)), np.zeros((2, 3))))
    assert np.isnan(B.knn_overlap(np.zeros((10, 3)), np.zeros((9, 3))))


def test_k_is_bounded_by_the_points_there_are():
    """Asking for 15 neighbours of 5 points must not index past the end."""
    Y = np.random.default_rng(3).normal(size=(5, 3))
    assert B.knn_overlap(Y, Y, k=15) == 1.0


# --------------------------------------------------------------------------- running both ways
def test_with_no_cuml_it_says_so_instead_of_showing_an_empty_bar(monkeypatch):
    """torch and cupy do array work and neither implements UMAP, so with only those installed there
    is no second map to compare -- which is a sentence, not a blank."""
    from starplast import gpu
    monkeypatch.setattr(gpu, "available", lambda: {"cuml": False, "cupy": True, "torch": True,
                                                   "device": "fake"})
    out = B.compare(_nodes(), EmbeddingSpec(blocks=("fitness_screens",)), sample=120,
                    log=lambda *_: None)
    assert out["gpu_available"] is False
    assert "gpu" not in out
    assert "no cuml" in out["note"] and "starplast-gpu" in out["note"]
    assert out["cpu"]["seconds"] > 0 and out["cpu"]["coords"].shape[1] == 3


def test_both_runs_happen_and_are_compared(monkeypatch):
    """Driven with a stand-in cuml, because the wiring has to be exercised on machines without one:
    two embeddings, two clocks, one overlap."""
    import types
    fake = types.ModuleType("cuml")
    fake.manifold = types.ModuleType("cuml.manifold")
    fake.__version__ = "26.08.00"

    class FakeUMAP:
        def __init__(self, **kw):
            pass

        def fit_transform(self, X):
            # A different embedding of the same data, which is exactly the situation being reported.
            rng = np.random.default_rng(7)
            return np.asarray(X)[:, :3] + rng.normal(scale=0.3, size=(len(X), 3))

    fake.manifold.UMAP = FakeUMAP
    monkeypatch.setitem(sys.modules, "cuml", fake)
    monkeypatch.setitem(sys.modules, "cuml.manifold", fake.manifold)

    out = B.compare(_nodes(), EmbeddingSpec(blocks=("fitness_screens",)), sample=150,
                    log=lambda *_: None)
    assert out["gpu"]["backend"].startswith("cuml")
    assert out["cpu"]["backend"].startswith("umap-learn") or "pca" in out["cpu"]["backend"]
    assert out["speedup"] > 0
    assert 0.0 <= out["knn_overlap"] <= 1.0
    assert "two maps of the same data" in out["note"]


def test_the_switch_is_left_as_it_was_found(monkeypatch):
    """It forces the setting for each run; leaving it forced would change what every later map is
    built by, from a button that says it only measures."""
    from starplast import gpu
    monkeypatch.setenv(gpu.ENV_GPU, "0")
    B.compare(_nodes(), EmbeddingSpec(blocks=("fitness_screens",)), sample=120, log=lambda *_: None)
    assert os.environ[gpu.ENV_GPU] == "0"
    monkeypatch.delenv(gpu.ENV_GPU, raising=False)
    B.compare(_nodes(), EmbeddingSpec(blocks=("fitness_screens",)), sample=120, log=lambda *_: None)
    assert gpu.ENV_GPU not in os.environ, "the comparison left the GPU forced afterwards"


def test_the_sample_bounds_the_work():
    """A demonstration, not a run: waiting four minutes for the full proteome to answer 'are these
    different' is not seeing it sooner."""
    out = B.compare(_nodes(400), EmbeddingSpec(blocks=("fitness_screens",)), sample=100,
                    log=lambda *_: None)
    assert out["sample"] == 100 and len(out["cpu"]["coords"]) == 100


def test_a_table_smaller_than_the_sample_is_used_whole():
    out = B.compare(_nodes(80), EmbeddingSpec(blocks=("fitness_screens",)), sample=500,
                    log=lambda *_: None)
    assert out["sample"] == 80


# --------------------------------------------------------------------------- the dialog
def test_the_dialog_shows_both_maps_and_both_clocks(win):
    rng = np.random.default_rng(0)
    result = {"sample": 400, "gpu_available": True, "speedup": 27.5, "knn_overlap": 0.71,
              "cpu": {"coords": rng.normal(size=(400, 3)), "seconds": 14.3,
                      "backend": "umap-learn 0.5.12", "genes": 400, "features": 25},
              "gpu": {"coords": rng.normal(size=(400, 3)), "seconds": 0.52,
                      "backend": "cuml 26.08.00", "genes": 400, "features": 25},
              "note": "27.5x faster, and 71% of each gene's 15 nearest neighbours are the same"}
    d = win.build_comparison(result)
    from PyQt6 import QtWidgets
    labels = [q.text() for q in d.findChildren(QtWidgets.QLabel) if q.text()]
    assert any("umap-learn 0.5.12" in t and "14.3 s" in t for t in labels)
    assert any("cuml 26.08.00" in t and "0.5 s" in t for t in labels)
    assert any("nearest neighbours" in t for t in labels)
    pictures = [q for q in d.findChildren(QtWidgets.QLabel)
                if q.pixmap() is not None and not q.pixmap().isNull()]
    assert len(pictures) == 2, "both maps have to be visible; one is not a comparison"


def test_the_dialog_says_which_side_is_missing(win):
    """With no cuml there is one map, and the empty side has to say why rather than look broken."""
    result = {"sample": 400, "gpu_available": False, "note": "no cuml on this machine",
              "cpu": {"coords": np.random.default_rng(0).normal(size=(50, 3)), "seconds": 1.0,
                      "backend": "umap-learn 0.5.12", "genes": 50, "features": 6}}
    d = win.build_comparison(result)
    from PyQt6 import QtWidgets
    labels = [q.text() for q in d.findChildren(QtWidgets.QLabel) if q.text()]
    assert any("not available" in t for t in labels)
    assert any("no cuml" in t for t in labels)


def test_the_button_is_in_preferences(win):
    d = win.build_preferences()
    assert win.gpu_test.text().startswith("compare CPU and GPU")
    assert win.gpu_test.toolTip()
    d.close()


@pytest.fixture(scope="module")
def qapp():
    """Held for the module's lifetime. Created and dropped, it is collected out from under the
    widgets built on it, and the interpreter aborts rather than failing a test."""
    from PyQt6 import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture(scope="module")
def win(qapp):
    from starplast.app import Window
    w = Window()
    yield w
    w.close()


def test_pressing_the_button_runs_it_off_the_gui_thread_and_opens_the_result(win, monkeypatch):
    """Two UMAPs is minutes on the full proteome, so it goes through the job runner like everything
    else slow. `exec` is stubbed because a modal dialog in a test hangs forever rather than failing.
    """
    from PyQt6 import QtWidgets
    import starplast.benchmark as BM
    shown = []
    monkeypatch.setattr(QtWidgets.QDialog, "exec", lambda self: shown.append(self) or 0)
    monkeypatch.setattr(BM, "compare", lambda nodes, spec=None, sample=0, log=print: {
        "sample": 10, "gpu_available": False, "note": "stub",
        "cpu": {"coords": np.zeros((10, 3)), "seconds": 0.1, "backend": "umap-learn x",
                "genes": 10, "features": 3}})
    job = win.compare_backends(sample=10)
    win.jobs.wait(20000)
    QtWidgets.QApplication.processEvents()
    assert job.state == "done"
    assert shown, "the comparison finished and showed nothing"


def test_a_comparison_that_fails_says_so_rather_than_opening_an_empty_window(win, monkeypatch):
    from PyQt6 import QtWidgets
    import starplast.benchmark as BM

    def boom(*a, **k):
        raise RuntimeError("no map today")

    monkeypatch.setattr(BM, "compare", boom)
    monkeypatch.setattr(QtWidgets.QDialog, "exec", lambda self: pytest.fail("opened anyway"))
    win.compare_backends(sample=10)
    win.jobs.wait(20000)
    QtWidgets.QApplication.processEvents()
    assert "comparison failed" in win.statusBar().currentMessage()


def test_another_job_finishing_first_does_not_open_the_dialog(win, monkeypatch):
    """The runner's signal carries EVERY job's completion, so the handler has to check the id. The
    unrelated job is made to finish while the comparison is still running, which is the only way to
    reach that check -- afterwards the handler has already disconnected itself."""
    import threading
    from PyQt6 import QtWidgets
    import starplast.benchmark as BM
    go = threading.Event()
    opened = []
    monkeypatch.setattr(QtWidgets.QDialog, "exec", lambda self: opened.append(1) or 0)
    monkeypatch.setattr(BM, "compare", lambda *a, **k: (go.wait(20), {
        "sample": 1, "gpu_available": False, "note": "",
        "cpu": {"coords": np.zeros((5, 3)), "seconds": 0.1, "backend": "x", "genes": 5,
                "features": 3}})[1])
    win.compare_backends(sample=5)
    other = win.run_job(lambda: "something else", "unrelated")
    for _ in range(200):
        QtWidgets.QApplication.processEvents()
        if other.state in ("done", "failed"):
            break
        time.sleep(0.02)
    assert other.state == "done"
    assert opened == [], "an unrelated job opened the comparison dialog"
    go.set()
    win.jobs.wait(30000)
    QtWidgets.QApplication.processEvents()
    assert opened == [1], "the comparison's own result never arrived"


def test_a_second_comparison_can_be_run(win, monkeypatch):
    """The handler disconnects itself when its own job lands. If it did not, the second run would
    open two dialogs, then three -- and the first missed disconnect is the one nobody notices."""
    from PyQt6 import QtWidgets
    import starplast.benchmark as BM
    result = {"sample": 1, "gpu_available": False, "note": "",
              "cpu": {"coords": np.zeros((5, 3)), "seconds": 0.1, "backend": "x", "genes": 5,
                      "features": 3}}
    monkeypatch.setattr(BM, "compare", lambda *a, **k: result)
    opened = []
    monkeypatch.setattr(QtWidgets.QDialog, "exec", lambda self: opened.append(1) or 0)
    for _ in range(2):
        win.compare_backends(sample=5)
        win.jobs.wait(20000)
        QtWidgets.QApplication.processEvents()
    assert opened == [1, 1], f"two runs opened {len(opened)} dialogs"
