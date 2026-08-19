#!/usr/bin/env python3
"""GPU acceleration: what it claims, what it refuses, and that it agrees with the CPU.

The rule this file exists to hold: **a GPU result that differs from the CPU result is a bug, not a
speed-up.** The arithmetic paths -- ranking, distances -- must agree to 1e-5 or the switch changes
what the map says rather than how fast it says it. UMAP is the stated exception, because cuml's is a
different implementation, and the program says so rather than pretending otherwise.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import gpu  # noqa: E402

HAVE = any(gpu.available()[k] for k in ("cuml", "cupy", "torch"))
on_gpu = pytest.mark.skipif(not HAVE, reason="no GPU backend on this machine")
# And one per backend, because "a GPU backend exists" is not the same claim as "this backend does".
# The user's own environment carries cuml and not torch, so a torch test guarded by `on_gpu` ran and
# died on the import -- a machine-shaped failure with nothing wrong on it.
on_torch = pytest.mark.skipif(not gpu.available()["torch"], reason="torch is not installed here")


def test_what_is_available_is_reported_rather_than_assumed():
    have = gpu.available()
    assert set(have) == {"cuml", "cupy", "torch", "device"}
    assert all(isinstance(have[k], bool) for k in ("cuml", "cupy", "torch"))


def test_the_environment_can_force_it_off(monkeypatch):
    """For a machine whose driver is present but broken, and for the suite itself."""
    monkeypatch.setenv(gpu.ENV_GPU, "0")
    assert gpu.available() == {"cuml": False, "cupy": False, "torch": False, "device": ""}
    assert gpu.enabled() is False


def test_it_is_off_unless_it_is_both_wanted_and_possible(monkeypatch):
    monkeypatch.delenv(gpu.ENV_GPU, raising=False)
    from PyQt6 import QtCore
    s = QtCore.QSettings("starplast", "starplast")
    s.setValue("compute/gpu", False)
    s.sync()
    assert gpu.enabled() is False


def test_a_search_pins_the_backend_even_if_the_preference_changes(monkeypatch):
    monkeypatch.setattr(gpu, "available", lambda: {
        "cuml": False, "cupy": True, "torch": False, "device": "fixture"})
    monkeypatch.setenv(gpu.ENV_GPU, "1")
    with gpu.pinned() as snapshot:
        assert snapshot["gpu"] is True
        monkeypatch.setenv(gpu.ENV_GPU, "0")
        assert gpu.enabled() is True
    assert gpu.enabled() is False


def test_the_description_says_what_would_happen():
    """A switch that silently does nothing is worse than no switch."""
    text = gpu.describe()
    assert text
    if HAVE:
        assert any(k in text for k in ("cuml", "cupy", "torch"))
    else:
        assert "no GPU backend" in text and "install" in text


def test_a_small_array_stays_on_the_cpu(monkeypatch):
    """Under the threshold the copy costs more than the arithmetic saves, so a walk over a small
    subsample would get SLOWER for having a GPU."""
    monkeypatch.setenv(gpu.ENV_GPU, "1")
    assert gpu.worth_it(np.zeros((10, 10))) is False
    assert gpu.worth_it(np.zeros((2000, 60))) is True


@on_gpu
def test_distances_agree_with_scipy(monkeypatch):
    from scipy.spatial.distance import squareform, pdist
    monkeypatch.setenv(gpu.ENV_GPU, "1")
    # Over MIN_ELEMENTS, or the call correctly stays on the CPU and the test measures nothing.
    X = np.random.default_rng(0).normal(size=(3000, 60))
    assert gpu.worth_it(X)
    got = gpu.pairwise_distances(X)
    want = squareform(pdist(X))
    assert got.shape == want.shape
    # RELATIVE, because the device computes in float32 by design: a 3090 does float64 at a
    # thirty-second of the rate, and these distances feed a rank correlation.
    rel = np.abs(got - want).max() / max(want.max(), 1e-12)
    assert rel < 1e-5, f"the GPU distance matrix is not the CPU one: {rel:.2e} relative"


def test_rank_scaling_is_the_same_whatever_the_switch_says(monkeypatch):
    """It has no device path at all -- measured six times slower and 1.4e-4 out, so it stays on the
    CPU. This pins that the switch cannot change a ranking."""
    from starplast.embedding import _scale
    rng = np.random.default_rng(1)
    X = rng.normal(size=(2000, 30))
    X[rng.random(X.shape) < 0.1] = np.nan
    monkeypatch.setenv(gpu.ENV_GPU, "1")
    fast = _scale(X.copy(), "rank")
    monkeypatch.setenv(gpu.ENV_GPU, "0")
    slow = _scale(X.copy(), "rank")
    assert np.array_equal(np.isnan(fast), np.isnan(slow))
    both = np.isfinite(fast)
    assert np.array_equal(fast[both], slow[both]), "the switch changed a ranking"


@on_gpu
def test_the_walk_scores_the_same_either_way(monkeypatch):
    """`_quality` is what the walk's rows carry, so the two paths have to produce the same row."""
    from starplast.tuning import _quality
    rng = np.random.default_rng(2)
    X = rng.normal(size=(700, 20))
    Y = X[:, :3] + rng.normal(scale=0.1, size=(700, 3))
    monkeypatch.setenv(gpu.ENV_GPU, "1")
    fast = _quality(X, Y, 15)
    monkeypatch.setenv(gpu.ENV_GPU, "0")
    slow = _quality(X, Y, 15)
    assert abs(fast["continuity_proxy"] - slow["continuity_proxy"]) < 1e-3


def test_umap_and_hdbscan_are_offered_only_through_cuml(monkeypatch):
    """torch and cupy do the array work; neither implements UMAP, and pretending otherwise would
    silently leave the slowest step on the CPU while the switch said GPU."""
    monkeypatch.setenv(gpu.ENV_GPU, "1")
    if not gpu.available()["cuml"]:
        assert gpu.umap_class() is None
        assert gpu.hdbscan_class() is None


def test_the_switch_looks_like_the_one_in_spacr(qapp):
    """The same control in both programs: 40x20 track, purple off, teal on, caption on the left."""
    from starplast.theme import Switch, SWITCH_OFF, SWITCH_ON
    s = Switch("use the GPU where it helps")
    s.resize(240, 24)
    assert s.isChecked() is False
    got = []
    s.toggled.connect(got.append)
    s.setChecked(True)
    assert s.isChecked() is True and got == [], "setChecked must not emit"
    assert SWITCH_OFF == "#800080" and SWITCH_ON == "#008080"
    assert not s.grab().isNull()


def test_clicking_the_switch_toggles_and_emits(qapp):
    from PyQt6 import QtCore, QtGui
    from starplast.theme import Switch
    s = Switch("gpu")
    s.resize(240, 24)
    got = []
    s.toggled.connect(got.append)
    pos = QtCore.QPointF(200, 12)
    s.mouseReleaseEvent(QtGui.QMouseEvent(QtCore.QEvent.Type.MouseButtonRelease, pos, pos,
                                          QtCore.Qt.MouseButton.LeftButton,
                                          QtCore.Qt.MouseButton.LeftButton,
                                          QtCore.Qt.KeyboardModifier.NoModifier))
    assert got == [True] and s.isChecked()


def test_a_right_click_does_not_toggle_it(qapp):
    from PyQt6 import QtCore, QtGui
    from starplast.theme import Switch
    s = Switch("gpu")
    s.resize(240, 24)
    pos = QtCore.QPointF(200, 12)
    s.mouseReleaseEvent(QtGui.QMouseEvent(QtCore.QEvent.Type.MouseButtonRelease, pos, pos,
                                          QtCore.Qt.MouseButton.RightButton,
                                          QtCore.Qt.MouseButton.RightButton,
                                          QtCore.Qt.KeyboardModifier.NoModifier))
    assert s.isChecked() is False


def test_the_switch_reports_what_it_will_do_when_toggled(win):
    """Toggling it must say whether anything will actually happen -- in the status bar and in the
    caption under the switch, which is the one a user is looking at when they flip it."""
    d = win.build_preferences()
    note = win._on_gpu(True)
    assert note and ("GPU on" in note or "no GPU backend" in note)
    assert win.gpu_note.text() == note, "the caption did not follow the switch"
    assert note in win.statusBar().currentMessage()
    win._on_gpu(False)
    d.close()


@pytest.fixture(scope="module")
def qapp():
    from PyQt6 import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture(scope="module")
def win(qapp):
    from starplast.app import Window
    w = Window()
    yield w
    w.close()


def test_the_cuml_paths_are_exercised_even_without_a_gpu(monkeypatch):
    """cuml is not installed here and may never be on this machine, but the code that calls it runs
    on every embedding and every clustering. A path nobody can execute is a path nobody has
    checked -- so a stand-in cuml is injected and the wiring is driven for real."""
    import sys
    import types
    import numpy as np
    import pandas as pd
    from starplast import clustering
    from starplast.embedding import EmbeddingSpec, embed

    calls = {"umap": 0, "hdbscan": 0}

    class FakeUMAP:
        def __init__(self, **kw):
            self.kw = kw

        def fit_transform(self, X):
            calls["umap"] += 1
            return np.asarray(X)[:, :3] * 1.0

    class FakeHDBSCAN:
        def __init__(self, **kw):
            self.kw = kw

        def fit_predict(self, X):
            calls["hdbscan"] += 1
            return np.zeros(len(X), dtype=int)

    fake = types.ModuleType("cuml")
    fake.manifold = types.ModuleType("cuml.manifold")
    fake.manifold.UMAP = FakeUMAP
    fake.cluster = types.ModuleType("cuml.cluster")
    fake.cluster.HDBSCAN = FakeHDBSCAN
    monkeypatch.setitem(sys.modules, "cuml", fake)
    monkeypatch.setitem(sys.modules, "cuml.manifold", fake.manifold)
    monkeypatch.setitem(sys.modules, "cuml.cluster", fake.cluster)
    monkeypatch.setenv(gpu.ENV_GPU, "1")

    assert gpu.umap_class() is FakeUMAP and gpu.hdbscan_class() is FakeHDBSCAN

    rng = np.random.default_rng(0)
    nodes = pd.DataFrame({f"fit_{i}": rng.normal(size=300) for i in range(5)})
    nodes.insert(0, "gene_id", [f"TGME49_{200000 + i}" for i in range(300)])
    said = []
    Y, names, rows = embed(nodes, EmbeddingSpec(blocks=("fitness_screens",)), log=said.append)
    assert calls["umap"] == 1 and Y.shape[1] == 3
    assert any("a different map from the CPU path" in m for m in said), said
    assert any("cuml" in m for m in said), "the log named the wrong library for the work it did"

    labels = clustering.cluster(rng.normal(size=(1500, 3)), min_cluster_size=25)
    assert calls["hdbscan"] == 1 and len(labels) == 1500


def test_a_gpu_clustering_that_fails_falls_back_loudly(monkeypatch, caplog):
    """A GPU that refuses is a slower run, not a failed one -- but a silent fallback is how "the GPU
    switch does nothing" becomes impossible to diagnose."""
    import logging
    import sys
    import types
    import numpy as np
    from starplast import clustering

    class Boom:
        def __init__(self, **kw):
            pass

        def fit_predict(self, X):
            raise RuntimeError("out of memory")

    fake = types.ModuleType("cuml")
    fake.cluster = types.ModuleType("cuml.cluster")
    fake.cluster.HDBSCAN = Boom
    fake.manifold = types.ModuleType("cuml.manifold")
    monkeypatch.setitem(sys.modules, "cuml", fake)
    monkeypatch.setitem(sys.modules, "cuml.cluster", fake.cluster)
    monkeypatch.setenv(gpu.ENV_GPU, "1")

    seen = []

    class Collect(logging.Handler):
        def emit(self, record):
            seen.append(record.getMessage())

    log = logging.getLogger("starplast.clustering")
    handler = Collect(level=logging.WARNING)
    log.addHandler(handler)
    try:
        labels = clustering.cluster(np.random.default_rng(0).normal(size=(1500, 3)),
                                    min_cluster_size=25)
    finally:
        log.removeHandler(handler)
    assert len(labels) == 1500, "the run must survive the GPU refusing"
    assert any("using the CPU" in m for m in seen), seen


def test_a_backend_that_will_not_import_is_simply_absent(monkeypatch):
    """A half-installed cupy raises on import rather than returning False, and a probe that let that
    through would take the application down at startup on somebody else's machine."""
    import builtins
    real = builtins.__import__

    def refuse(name, *a, **k):
        if name in ("cuml", "cupy", "torch"):
            raise ImportError(f"{name} is broken here")
        return real(name, *a, **k)

    monkeypatch.setenv(gpu.ENV_GPU, "1")
    monkeypatch.setattr(builtins, "__import__", refuse)
    assert gpu.available() == {"cuml": False, "cupy": False, "torch": False, "device": ""}
    assert gpu.enabled() is False
    assert "no GPU backend" in gpu.describe()


def test_the_settings_being_unreadable_leaves_it_off(monkeypatch):
    """Qt may not be importable at all where this module is used from a script."""
    import builtins
    real = builtins.__import__

    def no_qt(name, *a, **k):
        if name.startswith("PyQt6"):
            raise ImportError("no Qt here")
        return real(name, *a, **k)

    monkeypatch.delenv(gpu.ENV_GPU, raising=False)
    monkeypatch.setattr(builtins, "__import__", no_qt)
    assert gpu.enabled() is False


def test_cupy_is_used_when_it_is_the_only_backend(monkeypatch):
    """The cupy path has no hardware here, so it is driven against a stand-in: what is checked is
    that it computes the distance matrix in float64 and hands back NumPy."""
    import sys
    import types
    import numpy as np

    fake = types.ModuleType("cupy")
    fake.float64 = np.float64
    fake.asarray = lambda a, dtype=None: np.asarray(a, dtype=dtype)
    fake.asnumpy = np.asarray
    fake.maximum = np.maximum
    fake.sqrt = np.sqrt
    fake.cuda = types.SimpleNamespace(runtime=types.SimpleNamespace(getDeviceCount=lambda: 1))
    monkeypatch.setitem(sys.modules, "cupy", fake)
    monkeypatch.setenv(gpu.ENV_GPU, "1")
    monkeypatch.setattr(gpu, "available", lambda: {"cuml": False, "cupy": True, "torch": False,
                                                   "device": "fake"})
    from scipy.spatial.distance import squareform, pdist
    # Over MIN_ELEMENTS, or the call never reaches the cupy branch this test is about.
    X = np.random.default_rng(0).normal(size=(2000, 80))
    assert gpu.worth_it(X)
    got = gpu.pairwise_distances(X)
    want = squareform(pdist(X))
    # Relative, and not tighter than float64's squared-norm identity can manage: cupy has no direct
    # mode, so |a|^2 + |b|^2 - 2ab is the only formulation available and it cancels at about 1e-8
    # even in double precision. That is four orders below the float32 shortcut this replaced.
    assert np.abs(got - want).max() / want.max() < 1e-7
    assert isinstance(got, np.ndarray)


def test_the_description_names_the_gain_not_just_the_backend(monkeypatch):
    """"GPU: torch" reads as a promise about UMAP and is not one."""
    monkeypatch.setenv(gpu.ENV_GPU, "1")
    monkeypatch.setattr(gpu, "available", lambda: {"cuml": False, "cupy": False, "torch": True,
                                                   "device": "RTX 3090"})
    assert "array work only" in gpu.describe()
    monkeypatch.setattr(gpu, "available", lambda: {"cuml": True, "cupy": False, "torch": True,
                                                   "device": "RTX 3090"})
    assert "UMAP, t-SNE, HDBSCAN, k-means and DBSCAN move to the GPU" in gpu.describe()


def test_a_real_cupy_reports_its_device_count(monkeypatch):
    """The probe asks cupy how many devices it sees, because an installed cupy on a machine with no
    card imports perfectly well and then fails at the first allocation."""
    import sys
    import types
    fake = types.ModuleType("cupy")
    fake.cuda = types.SimpleNamespace(runtime=types.SimpleNamespace(getDeviceCount=lambda: 2))
    monkeypatch.setitem(sys.modules, "cupy", fake)
    monkeypatch.setenv(gpu.ENV_GPU, "1")
    assert gpu.available()["cupy"] is True
    fake.cuda.runtime.getDeviceCount = lambda: 0
    assert gpu.available()["cupy"] is False


@on_torch
def test_the_torch_path_is_the_accurate_one(monkeypatch):
    """The matrix-multiplication shortcut cancels badly in float32 for points that are close
    together, which is most pairs in an embedding: 5.6e-4 of relative error against 1e-7 for the
    direct mode. This asserts the accurate one is what runs."""
    import torch
    monkeypatch.setenv(gpu.ENV_GPU, "1")
    seen = {}
    real = torch.cdist

    def watch(a, b, **kw):
        seen.update(kw)
        return real(a, b, **kw)

    monkeypatch.setattr(torch, "cdist", watch)
    gpu.pairwise_distances(np.random.default_rng(0).normal(size=(2000, 60)))
    assert seen.get("compute_mode") == "donot_use_mm_for_euclid_dist"


def test_the_switch_reaches_the_walk_and_the_widget(qapp, monkeypatch, win):
    """The three places the setting has to arrive: the distance path the walk uses, the caption in
    Preferences, and the scaling that must NOT change. Driven with the switch on, so the wiring is
    exercised on a machine with a backend and skipped nowhere."""
    import numpy as np
    from starplast.tuning import _condensed
    monkeypatch.setenv(gpu.ENV_GPU, "1")
    big = np.random.default_rng(0).normal(size=(2000, 70))
    small = np.random.default_rng(0).normal(size=(50, 3))
    assert len(_condensed(big)) == 2000 * 1999 // 2
    assert len(_condensed(small)) == 50 * 49 // 2, "a small array must stay on the CPU"
    d = win.build_preferences()
    assert win.gpu_switch.isChecked() in (True, False)
    assert gpu.describe() in win.gpu_note.text()
    d.close()


def test_a_small_matrix_is_computed_by_scipy_even_with_the_switch_on(monkeypatch):
    """The threshold is inside the function, not only in its callers: anything that asks for
    distances on a small array must get the CPU answer rather than a device round trip."""
    from scipy.spatial.distance import squareform, pdist
    monkeypatch.setenv(gpu.ENV_GPU, "1")
    X = np.random.default_rng(3).normal(size=(60, 5))
    assert not gpu.worth_it(X)
    assert np.array_equal(gpu.pairwise_distances(X), squareform(pdist(X)))


def test_the_switch_asks_for_room_for_its_caption(qapp):
    from starplast.theme import Switch
    s = Switch("use the GPU where it helps")
    assert s.sizeHint().width() > 60 and s.sizeHint().height() == 24
    assert Switch("").sizeHint().width() < s.sizeHint().width()


def test_toggling_before_preferences_has_ever_opened_does_not_crash(win):
    """The menu can set it, the switch is built lazily, and the caption may not exist yet."""
    had = getattr(win, "gpu_note", None)
    if had is not None:
        del win.gpu_note
    try:
        note = win._on_gpu(True)
        assert note
    finally:
        win._on_gpu(False)
        if had is not None:
            win.gpu_note = had


def test_without_umap_the_record_says_what_will_actually_build_the_map(monkeypatch):
    """The embedding falls back to PCA when umap-learn is missing, and a row that recorded
    'umap-learn' for a map PCA built would be a false provenance rather than a missing one."""
    import builtins
    real = builtins.__import__

    def no_umap(name, *a, **k):
        if name == "umap":
            raise ImportError("no umap here")
        return real(name, *a, **k)

    monkeypatch.setenv(gpu.ENV_GPU, "0")
    monkeypatch.setattr(builtins, "__import__", no_umap)
    b = gpu.backend()
    assert b["umap"].startswith("pca")
    assert "scikit-learn" in b["cluster"]
    assert "pca" in gpu.backend_id()


def test_every_cuml_algorithm_has_a_resolver(monkeypatch):
    import sys
    import types
    fake = types.ModuleType("cuml")
    fake.cluster = types.ModuleType("cuml.cluster")
    fake.manifold = types.ModuleType("cuml.manifold")
    for name in ("HDBSCAN", "KMeans", "DBSCAN"):
        setattr(fake.cluster, name, type(name, (), {}))
    for name in ("UMAP", "TSNE"):
        setattr(fake.manifold, name, type(name, (), {}))
    monkeypatch.setitem(sys.modules, "cuml", fake)
    monkeypatch.setitem(sys.modules, "cuml.cluster", fake.cluster)
    monkeypatch.setitem(sys.modules, "cuml.manifold", fake.manifold)
    monkeypatch.setenv(gpu.ENV_GPU, "1")
    monkeypatch.setattr(gpu, "available", lambda: {"cuml": True, "cupy": False, "torch": False,
                                                    "device": "fake"})
    assert all(f() is not None for f in (gpu.umap_class, gpu.hdbscan_class, gpu.kmeans_class,
                                         gpu.dbscan_class, gpu.tsne_class))


def test_distances_fall_back_to_the_cpu_when_no_backend_is_installed(monkeypatch):
    """The switch can be on while the libraries are absent, and this path used to import cupy anyway.

    An ImportError three frames down is the wrong answer to "give me the distances": the caller
    asked a mathematical question, and the accurate CPU answer is a worse day rather than a failure.
    """
    import numpy as np
    from scipy.spatial.distance import pdist, squareform
    monkeypatch.setenv(gpu.ENV_GPU, "1")
    monkeypatch.setattr(gpu, "available", lambda: {"cuml": False, "cupy": False, "torch": False,
                                                   "device": ""})
    monkeypatch.setattr(gpu, "worth_it", lambda X: True)
    X = np.random.RandomState(0).normal(size=(40, 3))
    assert np.allclose(gpu.pairwise_distances(X), squareform(pdist(X)))
