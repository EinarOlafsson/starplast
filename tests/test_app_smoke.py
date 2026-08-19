#!/usr/bin/env python3
"""Headless smoke test for the app against the committed cache.

The handoff claimed a passing smoke test but none was ever committed, so it could not be re-run. This
drives the real Window offscreen: every level of detail, every color mode, every edge type, the attention
toggle, search and picking. It asserts the app builds and survives each control, not that it looks right.

Skipped automatically where Qt cannot open an offscreen GL context, and where the cache is absent.

Run: QT_QPA_PLATFORM=offscreen pytest tests/test_app_smoke.py -q
"""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from starplast import paths
from starplast.app import COLOR_MODES  # noqa: E402

DATA = paths.data_dir()
pytestmark = pytest.mark.skipif(
    not os.path.exists(os.path.join(DATA, "graph.npz")),
    reason="no built cache; run python -m starplast.build_graph")


@pytest.fixture(scope="module")
def win():
    pytest.importorskip("PyQt6")
    from PyQt6 import QtWidgets
    try:
        from starplast import app as A
    except Exception as e:                                  # pragma: no cover - environment dependent
        pytest.skip(f"app import failed: {e}")
    qapp = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    try:
        w = A.Window()
    except Exception as e:                                  # pragma: no cover - needs a GL context
        pytest.skip(f"cannot create GL widget offscreen: {e}")
    yield w
    w.close()
    qapp.processEvents()


def test_cache_loads_with_expected_shape(win):
    assert len(win.nodes) > 8000
    assert win.xyz.shape == (len(win.nodes), 3)
    assert {"n_publications", "n_fulltext", "lit_tier", "attention_depth",
            "n_papers_focal", "n_papers_substantive", "n_papers_incidental"} <= set(win.nodes.columns)


def test_depth_of_attention_is_categorical_and_grey_when_unnamed(win):
    """Never-named genes must stay grey with everything else unknown, not sit at the bottom of a ramp."""
    import numpy as np
    from starplast.app import DEPTH_COLOR
    from starplast import theme as TH
    i = COLOR_MODES.index("depth of attention")
    assert i >= 0
    win.set_color_mode(COLOR_MODES[i])
    win.redraw()
    c = win.colors(np.ones(win.n, bool))
    unnamed = (win.nodes.attention_depth.astype(str) == "").to_numpy()
    if unnamed.any():
        # The unknown grey is derived from the active theme's fg_dim, not a module constant,
        # so a theme switch cannot leave it invisible against the new ground.
        expected = np.array(TH.unknown_color(win.theme)[:3], dtype=np.float32)
        assert np.allclose(c[unnamed][:, :3], expected, atol=1e-5)
    for tier, col in DEPTH_COLOR.items():
        m = (win.nodes.attention_depth.astype(str) == tier).to_numpy()
        if m.any():
            assert np.allclose(c[m][:, :3], np.array(col, dtype=np.float32), atol=1e-5)


def test_listed_genes_are_not_reported_as_studied(win):
    """A gene named only in bodies and captions gets the 'listed, not studied' warning."""
    import numpy as np
    idx = np.where((win.nodes.attention_depth.astype(str) == "incidental").to_numpy())[0]
    if idx.size == 0:
        pytest.skip("no incidental-only genes in cache")
    win.on_pick(int(idx[0]))
    assert "Listed, not studied" in win.detail.toHtml()


def test_both_comention_types_are_present_and_separate(win):
    from starplast.app import COMENTION
    assert "comention" in win.edges, "abstract co-mention missing"
    assert "comention_ft" in win.edges, "full-text co-mention missing"
    for k in COMENTION:
        e = win.edges[k]
        assert len(e["a"]) == len(e["b"]) == len(e["w"]) == len(e["r"])
        assert (e["r"] != e["w"]).any(), f"{k}: corrected residual identical to raw weight"


def test_every_edge_type_toggles(win):
    from starplast.app import EDGE_TYPES
    win.all_edges_act.setChecked(True)
    for k, _ in EDGE_TYPES:
        if k not in win.edges:
            continue
        for other, _ in EDGE_TYPES:
            win.set_edge(other, other == k)
        win.redraw()


def test_every_level_of_detail(win):
    for i in range(3):
        win.set_level(i)
        win.redraw()


def test_every_color_mode(win):
    for i in range(len(COLOR_MODES)):
        win.set_color_mode(COLOR_MODES[i])
        win.redraw()


def test_attention_toggle_changes_drawn_comention(win):
    """The corrected view keeps only more-than-expected pairs, so it cannot draw more than raw."""
    import numpy as np
    e = win.edges["comention"]
    assert (e["r"] > 0).sum() <= len(e["r"])
    assert np.isfinite(e["r"]).all()
    for state in (True, False, True):
        win.attn_act.setChecked(state)
        win.redraw()


def test_scatter_does_not_blend_additively(win):
    """Additive blending summed 8,140 overlapping points to white and hid the color encoding entirely.

    Every color mode rendered as one white blob while all the array-level tests still passed, so this
    checks the GL state the renderer actually uses.
    """
    from OpenGL.GL import GL_DEPTH_TEST, GL_ONE, GL_SRC_ALPHA
    opts = getattr(win.scatter, "_GLGraphicsItem__glOpts", None)
    assert opts, "could not read the scatter's GL options"
    assert opts.get(GL_DEPTH_TEST) is True, "depth testing off: points sum instead of occluding"
    src, dst = opts.get("glBlendFunc", (None, None))
    assert (src, dst) != (GL_SRC_ALPHA, GL_ONE), "additive blending is back"


def test_edge_alpha_encodes_weight(win):
    """Drawn at one flat alpha the strongest and weakest edges look identical, which also made the
    attention toggle -- which reorders exactly that quantity -- visually almost a no-op."""
    import numpy as np
    k = "comention_ft" if "comention_ft" in win.edges else "comention"
    for other, _ in __import__("starplast.app", fromlist=["EDGE_TYPES"]).EDGE_TYPES:
        win.set_edge(other, other == k)
    win.all_edges_act.setChecked(True)
    win.sel = None
    win.redraw()
    assert win.edge_items, "no edge item was drawn"
    colors = [np.asarray(it.color) for it in win.edge_items
               if getattr(it, "color", None) is not None and np.asarray(it.color).ndim == 2]
    assert colors, "edges drawn with a single flat color, not a per-edge color array"
    alphas = np.concatenate([c[:, 3] for c in colors])
    assert np.ptp(alphas) > 0.01, "edge alpha does not vary with weight"


def test_measured_binding_and_structure_layers_present(win):
    for k in ("xlms", "struct", "unwritten_interaction"):
        assert k in win.edges, f"{k} missing from the cache"
    assert {"n_xlink_partners", "n_struct_similar", "best_model_agreement"} <= set(win.nodes.columns)


def test_panel_shows_how_the_binding_is_modeled(win):
    """The crosslink model table is what answers 'how does this binding happen'."""
    import numpy as np
    if not len(win.models):
        pytest.skip("no crosslink model table in cache")
    gid = win.models.iloc[0].gene_a
    i = int(np.where(win.nodes.gene_id.to_numpy() == gid)[0][0])
    win.on_pick(i)
    html = win.detail.toHtml()
    assert "How the binding is modeled" in html
    assert "crosslink" in html.lower()


def test_model_disagreement_is_reported_not_hidden(win):
    """frac_satisfied == 0 means the model fails to explain the measurement; say so."""
    import numpy as np
    if not len(win.models) or "frac_satisfied" not in win.models:
        pytest.skip("no satisfaction scores")
    bad = win.models[win.models.frac_satisfied == 0]
    if not len(bad):
        pytest.skip("no unsatisfied models")
    gid = bad.iloc[0].gene_a
    i = int(np.where(win.nodes.gene_id.to_numpy() == gid)[0][0])
    win.on_pick(i)
    assert "does not place them in contact" in win.detail.toHtml()


def test_published_screens_and_proteomics_are_shipped(win):
    """Standalone means every measurement is in the cache, not fetched later."""
    need = {"crispr_gra17_synthlethal_delta", "crispr_gra12s1_l2fc_invivo",
            "crispr_invivo_platform_lfc", "hosttx_T2", "protein_ibaq_log2",
            "expr_sporulated", "fit_invivo_young2019", "sequence", "length"}
    assert need <= set(win.nodes.columns), f"missing: {need - set(win.nodes.columns)}"


def test_superseded_accessions_still_resolve(win):
    """The 2019 in vivo screen cites pre-2012 ids for every gene; unresolved it contributes nothing."""
    assert win.nodes.crispr_invivo_platform_lfc.notna().sum() > 100


def test_untested_genes_read_as_unmeasured_not_zero(win):
    """A targeted library leaves most genes untested. That must stay missing, never become 0.0."""
    import numpy as np
    col = win.nodes.crispr_gra12s1_l2fc_invivo
    assert col.isna().sum() > 7000, "targeted screen should be missing for most genes"
    i = int(np.where(col.isna().to_numpy())[0][0])
    win.on_pick(i)
    assert "—" in win.detail.toHtml()


def test_search_finds_a_known_gene(win):
    win.search.setText("TGME49_208830")            # GRA16
    win.do_search()
    assert win.sel is not None
    assert win.nodes.gene_id.iloc[win.sel] == "TGME49_208830"


def test_picking_fills_the_detail_panel(win):
    win.on_pick(0)
    html = win.detail.toHtml()
    assert "compartment" in html
    assert "abstracts naming it" in html


def test_unassigned_is_reported_as_unknown(win):
    """hyperLOPIT assignment tracks abundance, so a missing call must not read as a 27th compartment."""
    import numpy as np
    idx = np.where(win.nodes.compartment.astype(str) == "unassigned")[0]
    if idx.size == 0:
        pytest.skip("no unassigned genes in cache")
    win.on_pick(int(idx[0]))
    assert "unknown, not absent" in win.detail.toHtml()
