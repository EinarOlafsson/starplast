"""Frozen, synthetic query fixtures test identity and context, not biology.

All accessions here are fabricated format-valid identifiers. No fixture is a
claim about a real gene, mapping, host relation or ground-truth measurement.
"""
import copy
import json
from dataclasses import replace

import pandas as pd
import pytest

from starplast import identity, organisms as O
from starplast.query import AliasRecord, BiologicalContext, EntityRef, GeneResolver, Query, Resolution


G1 = EntityRef(O.TOXOPLASMA, "gene", "TGME49_100001")
G2 = EntityRef(O.TOXOPLASMA, "gene", "TGME49_100002")
P1 = EntityRef(O.FALCIPARUM, "gene", "PF3D7_0100001")
HOST = EntityRef(O.HUMAN, "protein", "SYNTHETIC_PROTEIN_1")
CONTEXT = BiologicalContext(stage="tachyzoite", host=O.HUMAN, tissue="synthetic tissue",
                            condition="synthetic condition α", strain="synthetic strain")

# One canonical fixture for every entity/question type, including host protein.
FIXTURES = {
    "gene": Query(O.TOXOPLASMA, "gene", (G1,), context=CONTEXT, target="compartment"),
    "protein": Query(O.HUMAN, "protein", (HOST,), context=CONTEXT, output="evidence"),
    "class": Query(O.TOXOPLASMA, "class", target="compartment", values=("rhoptry",), context=CONTEXT),
    "label": Query(O.TOXOPLASMA, "label", target="compartment", context=CONTEXT, output="scorecard"),
    "gene_set": Query(O.TOXOPLASMA, "gene_set", (G1, G2), context=CONTEXT),
    "trait": Query(O.TOXOPLASMA, "trait", target="synthetic_continuous_trait", context=CONTEXT),
    "pair": Query(O.TOXOPLASMA, "pair", (G1, HOST), context=CONTEXT, output="agreement"),
}


@pytest.mark.parametrize("name", FIXTURES)
def test_frozen_entity_fixture_round_trip(name):
    original = FIXTURES[name]
    transported = json.loads(original.to_json())
    assert transported["kind"] == name
    assert transported["organism"] == original.organism
    assert transported["context"] == {
        "stage": "tachyzoite", "host": O.HUMAN, "tissue": "synthetic tissue",
        "condition": "synthetic condition α", "strain": "synthetic strain",
    }
    assert transported["schema_version"] == 1
    assert Query.from_dict(transported) == original
    assert Query.from_json(original.to_json()).to_json() == original.to_json()


def test_label_class_trait_and_context_have_distinct_addresses():
    label = FIXTURES["label"]
    same_output_class = replace(FIXTURES["class"], output=label.output)
    trait = replace(label, kind="trait")
    other_context = replace(label, context=replace(CONTEXT, condition="other condition"))
    assert len({q.to_json() for q in (label, same_output_class, trait, other_context)}) == 4


def test_multilabel_class_values_are_exact_unordered_addresses():
    query = replace(FIXTURES["class"], values=("rhoptry", "nucleus", "rhoptry"))
    assert query.values == ("nucleus", "rhoptry")
    assert Query.from_json(query.to_json()) == query
    assert replace(query, values=("nucleus", "Rhoptry")) != query


def test_gene_sets_are_unordered_but_pairs_preserve_direction_and_endpoint_species():
    assert replace(FIXTURES["gene_set"], entities=(G2, G1, G2)) == FIXTURES["gene_set"]
    forward = FIXTURES["pair"]
    reverse = replace(forward, organism=O.HUMAN, entities=(HOST, G1))
    assert forward.to_json() != reverse.to_json()
    assert Query.from_json(reverse.to_json()).entities == (HOST, G1)


@pytest.mark.parametrize("changes", [
    {"organism": "unknown"}, {"kind": "unknown"}, {"output": "unknown"},
    {"schema_version": 2}, {"schema_version": True}, {"schema_version": 1.0},
    {"schema_version": "1"}, {"context": {}}, {"target": " untrimmed"},
    {"entities": []}, {"entities": ("not an entity",)}, {"entities": ()},
    {"entities": (G1, G2)}, {"entities": (HOST,)}, {"entities": (P1,)},
    {"values": ["rhoptry"]}, {"values": ("rhoptry",)},
])
def test_invalid_gene_queries_fail_before_lookup(changes):
    with pytest.raises(ValueError):
        replace(FIXTURES["gene"], **changes)


@pytest.mark.parametrize("query,changes", [
    ("class", {"target": ""}), ("class", {"values": ()}),
    ("class", {"values": ("",)}), ("class", {"values": (1,)}),
    ("label", {"entities": (G1,)}), ("label", {"values": ("rhoptry",)}),
    ("trait", {"target": ""}), ("protein", {"entities": (G1,)}),
    ("gene_set", {"entities": ()}), ("gene_set", {"entities": (G1, P1)}),
    ("gene_set", {"entities": (HOST,)}), ("pair", {"entities": (G1,)}),
    ("pair", {"entities": (HOST, G1)}),
])
def test_entity_specific_contracts_are_checked(query, changes):
    with pytest.raises(ValueError):
        replace(FIXTURES[query], **changes)


@pytest.mark.parametrize("changes", [
    {"organism": O.FALCIPARUM}, {"organism": O.HUMAN}, {"organism": O.MOUSE},
    {"identifier": "ALIAS1"}, {"identifier": ""},
    {"identifier": " TGME49_100001"}, {"identifier": 1}, {"kind": "pair"},
])
def test_entity_identity_is_explicit_and_typed(changes):
    with pytest.raises(ValueError):
        replace(G1, **changes)


def test_host_addresses_do_not_assert_host_gene_pack_availability():
    host_gene = EntityRef(O.MOUSE, "gene", "SYNTHETIC_HOST_GENE_1")
    query = Query(O.MOUSE, "gene", (host_gene,), context=BiologicalContext(tissue="synthetic"))
    assert Query.from_json(query.to_json()) == query
    assert O.MOUSE not in O.SPACES  # references are not gene-space admission


@pytest.mark.parametrize("changes", [{"host": O.TOXOPLASMA}, {"host": "unknown"},
                                        {"stage": 1}, {"condition": " trailing "}])
def test_invalid_contexts_are_rejected(changes):
    with pytest.raises(ValueError):
        replace(CONTEXT, **changes)


@pytest.mark.parametrize("path,value", [
    ((), None), (("extra",), "ignored?"), (("schema_version",), 9),
    (("context", "extra"), "ignored?"), (("context",), []),
    (("entities",), "not an array"), (("values",), "not an array"),
    (("entities", 0, "extra"), "ignored?"),
])
def test_json_never_silently_discards_or_coerces_fields(path, value):
    data = copy.deepcopy(FIXTURES["gene"].to_dict())
    if not path:
        data = value
    else:
        parent = data
        for part in path[:-1]:
            parent = parent[part]
        parent[path[-1]] = value
    with pytest.raises(ValueError):
        Query.from_json(json.dumps(data))


def test_missing_and_duplicate_json_fields_are_rejected():
    data = FIXTURES["gene"].to_dict()
    del data["context"]
    with pytest.raises(ValueError):
        Query.from_dict(data)
    with pytest.raises(ValueError, match="Repeated"):
        Query.from_json('{"kind":"gene","kind":"class"}')


def test_resolution_preserves_ambiguities_sources_and_unknowns():
    records = [AliasRecord("SAME1", G1, "symbol", "synthetic-source-a:row1"),
               AliasRecord("SAME1", G2, "symbol", "synthetic-source-b:row2"),
               AliasRecord("SHORT", G1, "alias", "synthetic-source-a:row3"),
               AliasRecord("SHORT", G1, "alias", "synthetic-source-b:row4")]
    resolver = GeneResolver(O.TOXOPLASMA, records)
    ambiguous = resolver.resolve(" same-1 ")
    assert ambiguous.status == "ambiguous" and ambiguous.entity is None
    assert {r.entity for r in ambiguous.choices} == {G1, G2}
    assert {r.source for r in ambiguous.choices} == {records[0].source, records[1].source}
    sole_gene = resolver.resolve("short")
    assert sole_gene.status == "resolved" and sole_gene.entity == G1
    assert len(sole_gene.choices) == 2  # sources are not extra candidate genes
    for text in ("", "absent", "we measured SHORT", "SHORT-extra"):
        result = resolver.resolve(text)
        assert result.status == "unresolved" and result.entity is None
        assert result.organism == O.TOXOPLASMA and result.query == text


def test_same_symbol_in_another_organism_never_changes_the_result():
    tg = GeneResolver(O.TOXOPLASMA, [AliasRecord("SAME1", G1, "symbol", "synthetic-a")])
    pf = GeneResolver(O.FALCIPARUM, [AliasRecord("SAME1", P1, "symbol", "synthetic-b")])
    assert tg.resolve("SAME1").entity == G1
    assert pf.resolve("SAME1").entity == P1
    with pytest.raises(ValueError, match="explicit organism"):
        GeneResolver(O.TOXOPLASMA, [AliasRecord("SAME1", P1, "symbol", "synthetic-b")])


@pytest.mark.parametrize("changes", [{"source": ""}, {"source": 1}, {"kind": ""},
                                     {"alias": "---"}, {"entity": HOST}])
def test_mapping_sources_and_gene_identity_are_required(changes):
    with pytest.raises(ValueError):
        replace(AliasRecord("ALIAS1", G1, "symbol", "synthetic-source"), **changes)


def test_existing_index_adapter_restores_collisions_without_changing_prose_matching(tmp_path):
    path = tmp_path / "synthetic_identity.tsv"
    pd.DataFrame([{"gene_id": gene.identifier, "gene_name": "SAME1", "previous_ids": ""}
                  for gene in (G1, G2)]).to_csv(path, sep="\t", index=False)
    index = identity.build_index([G1.identifier, G2.identifier], str(path), log=lambda *_: None)
    resolver = GeneResolver.from_index(index, O.TOXOPLASMA, "synthetic-index:v1")
    result = resolver.resolve("SAME1")
    assert result.status == "ambiguous"
    assert {r.entity for r in result.choices} == {G1, G2}
    assert all(r.source == "synthetic-index:v1" and r.kind == "ambiguous_alias" for r in result.choices)
    assert resolver.resolve(G1.identifier.lower()).entity == G1
    assert list(index.find("SAME1")) == []  # literature extraction stays conservative
    with pytest.raises(ValueError, match="explicit organism"):
        GeneResolver.from_index(index, O.FALCIPARUM, "synthetic-index:v1")
    with pytest.raises(ValueError, match="explicit organism"):
        GeneResolver.from_index(index, O.HUMAN, "synthetic-index:v1")
    with pytest.raises(ValueError):
        GeneResolver.from_index(index, O.TOXOPLASMA, "")


def test_non_string_alias_is_an_explicit_error():
    with pytest.raises(ValueError, match="string"):
        GeneResolver(O.TOXOPLASMA, []).resolve(1)


def test_a_manually_constructed_resolution_cannot_hide_foreign_choices():
    with pytest.raises(ValueError, match="explicit organism"):
        Resolution(O.TOXOPLASMA, "ALIAS1", (AliasRecord("ALIAS1", P1, "alias", "synthetic"),))
