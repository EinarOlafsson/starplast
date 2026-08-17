#!/usr/bin/env python3
"""The Plasmodium relation layers, and the duplicate that Toxoplasma never showed.

An edge drawn twice reads as twice the evidence, and the category layers are exposed to it because
two proteins can share more than one InterPro domain. The Toxoplasma arm builds this layer the same
way and emits no duplicates at all, so the fault was invisible there for the whole life of the
project; Plasmodium produced 10,764 duplicate emissions among 34,887 the first time it ran, because
the var, rifin and stevor families share whole multi-domain architectures. The count is now the
weight, which is both correct and more informative than the 1.0 it replaced.

The other thing worth failing over is the index space. These indices point into pf_nodes.parquet and
mean a different gene in the Toxoplasma graph, which is exactly why there are two files.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import pf_graph as G  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GRAPH = os.path.join(ROOT, "starplast", "data", G.GRAPH)
NODES = os.path.join(ROOT, "starplast", "data", "pf_nodes.parquet")


def _nodes(n=4, **columns):
    frame = pd.DataFrame({"gene_id": [f"PF3D7_010{i:04d}" for i in range(n)]})
    for name, values in columns.items():
        frame[name] = values
    return frame


# --------------------------------------------------------------------------- category layers
def test_genes_sharing_an_orthogroup_are_joined():
    a, b, w = G.orthogroup_edges(_nodes(3, orthogroup=["OG6_1", "OG6_1", "OG6_2"]))
    assert list(zip(a, b)) == [(0, 1)] and list(w) == [1.0]


def test_a_group_of_one_makes_no_edge():
    a, _b, _w = G.orthogroup_edges(_nodes(2, orthogroup=["OG6_1", "OG6_2"]))
    assert len(a) == 0


def test_placeholders_are_not_treated_as_a_group():
    """`nan`, `unassigned` and an empty string are absence, and joining on absence joins everything."""
    a, _b, _w = G.orthogroup_edges(_nodes(4, orthogroup=["", "nan", "unassigned", "N/A"]))
    assert len(a) == 0


def test_a_group_larger_than_the_cap_is_dropped():
    big = _nodes(G.GROUP_CAP + 2, orthogroup=["OG6_1"] * (G.GROUP_CAP + 2))
    assert len(G.orthogroup_edges(big)[0]) == 0


def test_a_pair_sharing_several_domains_is_one_edge_weighted_by_how_many():
    """The duplicate. Two proteins with the same three domains are one relationship, not three."""
    a, b, w = G.domain_edges(_nodes(2, interpro_ids=["IPR1;IPR2;IPR3", "IPR1;IPR2;IPR3"]))
    assert list(zip(a, b)) == [(0, 1)]
    assert list(w) == [3.0]


def test_a_repeated_domain_within_one_protein_makes_no_self_edge():
    """InterPro reports one row per MATCH, so a repeat family lists the same protein twice."""
    a, b, _w = G.domain_edges(_nodes(2, interpro_ids=["IPR1;IPR1;IPR1", "IPR1"]))
    assert not any(x == y for x, y in zip(a, b))
    assert list(zip(a, b)) == [(0, 1)]


def test_a_gene_with_no_domains_joins_nothing():
    a, _b, _w = G.domain_edges(_nodes(3, interpro_ids=["", "nan", "IPR1"]))
    assert len(a) == 0


def test_missing_columns_yield_no_edges():
    assert len(G.orthogroup_edges(_nodes(2))[0]) == 0
    assert len(G.domain_edges(_nodes(2))[0]) == 0


# --------------------------------------------------------------------------- correlation layer
def _series(n, pattern):
    frame = _nodes(n)
    for i, column in enumerate(G.STAGE_COLUMNS):
        frame[column] = [row[i] for row in pattern]
    return frame


def test_correlated_genes_are_joined_and_uncorrelated_ones_are_not():
    rng = np.random.default_rng(0)
    base = rng.normal(size=len(G.STAGE_COLUMNS))
    rows = [list(base), list(base + 0.001), list(-base), list(rng.normal(size=len(G.STAGE_COLUMNS)))]
    frame = _series(4, rows)
    frame = pd.concat([frame] * 20, ignore_index=True)     # past the 50-row floor
    a, b, w = G.correlation_edges(frame)
    assert len(a) and (w >= G.CORRELATION).all()
    assert all(x != y for x, y in zip(a, b))


def test_too_few_stage_columns_yields_no_correlation():
    frame = _nodes(60)
    frame[G.STAGE_COLUMNS[0]] = 1.0
    frame[G.STAGE_COLUMNS[1]] = 2.0
    assert len(G.correlation_edges(frame)[0]) == 0


def test_too_few_complete_rows_yields_no_correlation():
    rng = np.random.default_rng(1)
    frame = _series(4, rng.normal(size=(4, len(G.STAGE_COLUMNS))).tolist())
    assert len(G.correlation_edges(frame)[0]) == 0


# --------------------------------------------------------------------------- assembly
def test_build_names_layers_the_way_the_graph_file_expects():
    frame = _nodes(3, orthogroup=["OG6_1", "OG6_1", "OG6_2"],
                   interpro_ids=["IPR1", "IPR1", ""])
    layers = G.build(frame, log=lambda *a: None)
    assert {"orthogroup__a", "orthogroup__b", "orthogroup__w"} <= set(layers)
    assert {"domain__a", "domain__b", "domain__w"} <= set(layers)


def test_an_empty_table_builds_no_layers():
    assert G.build(pd.DataFrame(), log=lambda *a: None) == {}


def test_save_writes_nothing_when_there_is_nothing_to_write(tmp_path):
    path = str(tmp_path / "empty.npz")
    assert G.save(pd.DataFrame(), path, log=lambda *a: None) == {}
    assert not os.path.exists(path)


def test_save_and_load_round_trip(tmp_path):
    frame = _nodes(3, orthogroup=["OG6_1", "OG6_1", "OG6_2"])
    data = tmp_path / "starplast" / "data"
    data.mkdir(parents=True)
    G.save(frame, str(data / G.GRAPH), log=lambda *a: None)
    loaded = G.load(str(tmp_path))
    assert loaded is not None and "orthogroup__a" in loaded.files


def test_load_returns_nothing_when_the_graph_is_not_built(tmp_path):
    assert G.load(str(tmp_path)) is None


# --------------------------------------------------------------------------- the shipped graph
@pytest.mark.skipif(not os.path.exists(GRAPH), reason="Plasmodium graph not built")
def test_no_layer_draws_an_edge_twice_or_joins_a_gene_to_itself():
    z = np.load(GRAPH)
    for key in [k for k in z.files if k.endswith("__a")]:
        label = key[:-3]
        a, b = z[key], z[f"{label}__b"]
        assert not (a == b).any(), f"{label} has a self-edge"
        assert len({(int(x), int(y)) for x, y in zip(a, b)}) == len(a), f"{label} has duplicates"
        assert (a < b).all(), f"{label} is not in canonical order"


@pytest.mark.skipif(not (os.path.exists(GRAPH) and os.path.exists(NODES)),
                    reason="Plasmodium cache not built")
def test_every_index_names_a_row_of_the_plasmodium_table():
    """The two graphs are separate files because an index means a row of ONE table."""
    z = np.load(GRAPH)
    n = len(pd.read_parquet(NODES, columns=["gene_id"]))
    for key in [k for k in z.files if k.endswith(("__a", "__b"))]:
        assert z[key].min() >= 0 and z[key].max() < n, key


@pytest.mark.skipif(not os.path.exists(GRAPH), reason="Plasmodium graph not built")
def test_the_domain_weight_records_how_many_domains_a_pair_shares():
    z = np.load(GRAPH)
    w = z["domain__w"]
    assert w.min() == 1.0 and w.max() > 1.0, "the weight is not carrying the shared-domain count"


def test_build_includes_the_correlation_layer_when_the_stages_are_there():
    """The category layers and the correlation layer are assembled by different code paths."""
    rng = np.random.default_rng(2)
    base = rng.normal(size=len(G.STAGE_COLUMNS))
    rows = [list(base), list(base + 0.001), list(rng.normal(size=len(G.STAGE_COLUMNS)))]
    frame = pd.concat([_series(3, rows)] * 25, ignore_index=True)
    said = []
    layers = G.build(frame, log=said.append)
    assert {"coexpression__a", "coexpression__b", "coexpression__w"} <= set(layers)
    assert any("coexpression" in line for line in said)


# --------------------------------------------------------------------------- measured contacts
XL = os.path.join(ROOT, "datasets", "reference", "plasmodb", G.CROSSLINK[0], G.CROSSLINK[1],
                  G.CROSSLINK[2])


def _xl_root(tmp_path, rows=None, sheet=None):
    rows = rows if rows is not None else [
        ("sp|Q1|PF3D7_010000", "sp|Q2|PF3D7_010001", 5),   # parasite-parasite: kept
        ("sp|P02549|SPTA1", "sp|Q2|PF3D7_010001", 12),      # human at one end: dropped
        ("sp|P02549|SPTA1", "sp|P11277|SPTB", 7),           # human both ends: dropped
        ("sp|Q1|PF3D7_010000", "sp|Q1|PF3D7_010000", 3),    # homomeric: dropped
        ("sp|Q9|PF3D7_099999", "sp|Q2|PF3D7_010001", 2),    # not in the table: dropped
    ]
    folder = tmp_path / "reference" / "plasmodb" / G.CROSSLINK[0] / G.CROSSLINK[1]
    folder.mkdir(parents=True)
    pd.DataFrame(rows, columns=["Protein1", "Protein2", "Num_Crosslinks"]).to_excel(
        folder / G.CROSSLINK[2], sheet_name=sheet or G.CROSSLINK_SHEET, index=False)
    return str(tmp_path)


def _xl_nodes():
    return pd.DataFrame({"gene_id": ["PF3D7_010000", "PF3D7_010001", "PF3D7_010002"]})


def test_a_host_protein_at_either_end_is_not_a_parasite_edge(tmp_path):
    """The experiment crosslinked parasite inside erythrocyte, so a third of the pairs are human."""
    a, b, w = G.crosslink_edges(_xl_nodes(), _xl_root(tmp_path), log=lambda *a: None)
    assert list(zip(a, b)) == [(0, 1)]
    assert list(w) == [5.0]


def test_a_homomeric_crosslink_is_not_an_edge(tmp_path):
    """It is evidence a protein self-associates, and drawing it puts a zero-length line in the graph."""
    a, b, _w = G.crosslink_edges(_xl_nodes(), _xl_root(tmp_path), log=lambda *a: None)
    assert not any(x == y for x, y in zip(a, b))


def test_a_protein_absent_from_the_table_is_dropped(tmp_path):
    a, _b, _w = G.crosslink_edges(_xl_nodes(), _xl_root(tmp_path), log=lambda *a: None)
    assert len(a) == 1


def test_missing_or_unreadable_crosslinks_yield_no_layer(tmp_path):
    assert len(G.crosslink_edges(_xl_nodes(), str(tmp_path), log=lambda *a: None)[0]) == 0
    assert len(G.crosslink_edges(pd.DataFrame(), _xl_root(tmp_path), log=lambda *a: None)[0]) == 0
    assert len(G.crosslink_edges(_xl_nodes(), _xl_root(tmp_path / "b", sheet="Other"),
                                 log=lambda *a: None)[0]) == 0


def test_a_crosslink_sheet_without_the_protein_columns_is_refused(tmp_path):
    folder = tmp_path / "reference" / "plasmodb" / G.CROSSLINK[0] / G.CROSSLINK[1]
    folder.mkdir(parents=True)
    pd.DataFrame({"wrong": [1]}).to_excel(folder / G.CROSSLINK[2],
                                          sheet_name=G.CROSSLINK_SHEET, index=False)
    assert len(G.crosslink_edges(_xl_nodes(), str(tmp_path), log=lambda *a: None)[0]) == 0


def test_a_pair_with_no_crosslink_count_still_becomes_an_edge(tmp_path):
    root = _xl_root(tmp_path, rows=[("sp|Q1|PF3D7_010000", "sp|Q2|PF3D7_010001", None)])
    a, _b, w = G.crosslink_edges(_xl_nodes(), root, log=lambda *a: None)
    assert len(a) == 1 and w[0] == 1.0


def test_build_includes_the_crosslink_layer_only_when_given_a_dataset_root(tmp_path):
    nodes = _xl_nodes().assign(orthogroup=["OG6_1", "OG6_1", "OG6_2"])
    without = G.build(nodes, log=lambda *a: None)
    assert "xlms__a" not in without
    with_root = G.build(nodes, log=lambda *a: None, dataset_root=_xl_root(tmp_path))
    assert "xlms__a" in with_root


@pytest.mark.skipif(not os.path.exists(XL), reason="crosslinks not fetched")
def test_the_crosslinks_recover_complexes_that_have_to_be_there():
    """PTEX is three subunits that translocate proteins together; a crosslink map that misses them
    is not measuring contacts. EXP2, PTEX150 and HSP101 must appear joined to one another."""
    nodes = pd.read_parquet(NODES)
    a, b, _w = G.crosslink_edges(nodes, os.path.join(ROOT, "datasets"), log=lambda *a: None)
    genes = nodes["gene_id"].to_numpy()
    pairs = {frozenset((genes[i], genes[j])) for i, j in zip(a, b)}
    ptex = {"PF3D7_1471100": "EXP2", "PF3D7_1436300": "PTEX150", "PF3D7_1116800": "HSP101"}
    joined = [p for p in pairs if len(p & set(ptex)) == 2]
    assert joined, f"no PTEX subunit pair among {len(pairs)} crosslinked pairs"
    assert 20 < len(pairs) < 200


def test_a_crosslink_table_with_no_usable_pair_yields_no_layer(tmp_path):
    """Every row human or homomeric: readable, and nothing to draw."""
    root = _xl_root(tmp_path, rows=[("sp|P02549|SPTA1", "sp|P11277|SPTB", 7),
                                    ("sp|Q1|PF3D7_010000", "sp|Q1|PF3D7_010000", 3)])
    assert len(G.crosslink_edges(_xl_nodes(), root, log=lambda *a: None)[0]) == 0
