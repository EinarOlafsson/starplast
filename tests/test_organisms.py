"""The organism registry must say exactly what the code already assumes, before anything reads it.

Step R0 of instruction 53: the registry is declared alongside the literals it will replace, and this
test holds the two to each other. When a literal is replaced by a registry lookup, its check here
keeps passing by construction; while it is still a literal, a drift between the two fails here
rather than as a wrong table in a window.
"""
from __future__ import annotations

import os
import sys

import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from starplast import organisms as O  # noqa: E402


def test_the_registry_reproduces_the_slot_tables():
    from starplast import slots
    assert O.table_map() == slots.SPECIES_TABLES
    assert O.prefix_map() == slots.SPECIES_PREFIXES
    for code, bridges in slots.SPECIES_BRIDGE_TABLES.items():
        assert O.get(code).host_bridges == bridges["host"]


def test_the_registry_reproduces_the_calibration_targets():
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    import calibrate_strategies as CS
    for code in CS.TARGETS:
        assert O.get(code).targets == CS.TARGETS[code]
        assert O.get(code).numbers == CS.NUMBERS[code]


def test_the_registry_reproduces_the_slot_generator_stages():
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    import generate_slot_table as G
    for code, stages in G.STAGES_BY_ORGANISM.items():
        assert O.get(code).contexts == frozenset(stages)


def test_the_registry_reproduces_the_localization_abbreviations():
    from starplast import localization
    assert {O.get(c).orthomcl for c in O.codes()} <= set(localization.SPECIES)


def test_the_registry_reproduces_the_window_species():
    pytest.importorskip("PyQt6")
    from starplast import app
    for name, where in app.SPECIES.items():
        space = O.by_species(name)
        assert (space.code, space.nodes, space.graph) == (where["code"], where["nodes"],
                                                           where["graph"])


@pytest.mark.parametrize("code", O.codes())
def test_every_shipped_table_is_the_space_it_claims(code):
    path = O.nodes_path(code)
    if not os.path.exists(path):
        pytest.skip(f"{code} table not built here")
    ids = pd.read_parquet(path, columns=["gene_id"])["gene_id"].astype(str)
    assert ids.map(O.get(code).matches).all(), ids[~ids.map(O.get(code).matches)].head().tolist()
    assert O.detect(ids) == code


def test_detection_goes_by_majority_not_by_one_stray_accession():
    tg = [f"TGME49_{i:06d}" for i in range(100)]
    assert O.detect(tg + ["PF3D7_0100100"]) == "Tg"
    assert O.detect(["PF3D7_0100100"] * 10 + tg[:3]) == "Pf"
    assert O.detect(["ENSG00000141510", "nonsense"]) is None


def test_partners_are_registered_and_reciprocal_between_the_parasites():
    for code in O.codes():
        partner = O.get(code).partner
        if partner:
            assert partner in O.SPACES
    assert O.get("Tg").partner == "Pf" and O.get("Pf").partner == "Tg"


def test_a_duplicate_or_malformed_space_is_refused():
    with pytest.raises(ValueError, match="twice"):
        O.register(O.get("Tg"))
    bad = O.Space(code="Xx", species="X", reference="x", kind="symbiont", gene_regex="X",
                  prefixes=("X",), database="", record_url="", nodes="x.parquet")
    with pytest.raises(ValueError, match="kind"):
        O.register(bad)


def test_record_links_name_the_gene():
    assert O.get("Tg").record("TGME49_208830").endswith("/gene/TGME49_208830")
    assert "plasmodb.org" in O.get("Pf").record("PF3D7_1133400")
