"""Species-specific displays must not silently change the calibrated inference recipe."""
import numpy as np
import pandas as pd

from starplast import embedding as E, organisms as O, strategies as S


def test_display_uses_its_own_slots_and_resolves_saved_recipes():
    nodes = pd.read_parquet(O.nodes_path(O.FALCIPARUM))
    spec = E.default_spec(nodes)
    assert spec.blocks and all(b.startswith(O.FALCIPARUM + "_") for b in spec.blocks)
    selected = E.columns_for(nodes, spec)
    columns = [c for values in selected.values() for c in values]
    assert len(columns) == len(set(columns))
    from starplast.slots import source_columns
    eligible = {c for slot in E.slot_blocks(O.FALCIPARUM).values() for c in source_columns(nodes, slot)}
    assert set(columns) == eligible
    matrix, features, rows = E.build_matrix(nodes, spec, log=lambda *a: None)
    assert len(features) > 100 and matrix.shape == (len(nodes), len(features))
    assert np.isfinite(matrix).all()
    saved = E.EmbeddingSpec.from_dict(spec.to_dict())
    again = E.build_matrix(nodes, saved, log=lambda *a: None)
    np.testing.assert_array_equal(matrix, again[0])
    assert features == again[1]
    other = pd.read_parquet(O.nodes_path(O.TOXOPLASMA))
    assert E.columns_for(other, spec) == {}


def test_inference_grouping_is_independent_of_pf_display_defaults(monkeypatch):
    ctx = S.Context.shipped(O.FALCIPARUM)
    targets = (None, *O.get(O.FALCIPARUM).targets)
    before = {t: ctx.blocks(t) for t in targets}

    def fail(*args):
        raise AssertionError("display recipe must not be consulted")

    monkeypatch.setattr(E, "default_spec", fail)
    after = S.Context.shipped(O.FALCIPARUM)
    assert {t: after.blocks(t) for t in targets} == before
    assert E.inference_spec(ctx.nodes).blocks == E.EmbeddingSpec().blocks


def test_drop_columns_reports_resolved_feature_names():
    nodes = pd.read_parquet(O.nodes_path(O.FALCIPARUM))
    spec = E.default_spec(nodes)
    spec.na_policy = "drop_columns"
    matrix, names, rows = E.build_matrix(nodes, spec, log=lambda *a: None)
    assert matrix.shape[1] == len(names) and np.isfinite(matrix).all()


def test_pf_panel_and_optimizer_offer_pf_slots(qapp):
    from starplast.analysis_panel import AnalysisPanel
    from starplast.optimize import block_pool
    nodes = pd.read_parquet(O.nodes_path(O.FALCIPARUM))
    panel = AnalysisPanel(nodes)
    try:
        assert set(panel._slot_blocks) == set(E.slot_blocks(O.FALCIPARUM))
        assert set(panel.blocks.checked()) == set(E.default_spec(nodes).blocks)
        pool = block_pool(nodes)
        assert pool and all(b.startswith(O.FALCIPARUM + "_") for b in pool)
    finally:
        panel.close()
