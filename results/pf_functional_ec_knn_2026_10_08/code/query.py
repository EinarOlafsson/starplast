"""Organism-qualified addresses for evidence, inference and scorecard questions.

This is a transport contract, not a claim that a question has data or a valid
benchmark. Host protein references can be addressed before host gene packs exist.
Canonical gene membership and alias provenance come from an explicitly supplied
resolver; no organism is inferred from a symbol or silently substituted.

Labels name a target (``compartment``); classes additionally name its values
(``rhoptry``). Class values and gene sets are unordered, whereas pair endpoints
are ordered and independently organism-qualified. Context is part of the address.
Versioned JSON retains these distinctions without depending on UI or data loading.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field

from . import organisms
from .identity import GeneIndex, norm

SCHEMA_VERSION = 1
QUERY_KINDS = frozenset({"gene", "protein", "gene_set", "label", "class", "trait", "pair"})
OUTPUTS = frozenset({"evidence", "inferences", "scorecard", "agreement", "meta_inference"})


def _text(value, name, optional=False):
    if not isinstance(value, str) or value != value.strip() or (not optional and not value):
        raise ValueError(f"{name} must be a {'possibly empty ' if optional else 'nonempty '}trimmed string")


def _organism(code):
    _text(code, "organism")
    if code not in organisms.SPACES and code not in organisms.HOST_TABLES:
        raise ValueError(f"Unknown organism {code!r}")


def _object(value, fields):
    if not isinstance(value, dict) or set(value) != set(fields):
        raise ValueError(f"Expected an object with exactly these fields: {', '.join(fields)}")
    return value


def _array(value, name):
    if not isinstance(value, list):
        raise ValueError(f"{name} must be a JSON array")
    return value


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Repeated JSON field {key!r}")
        result[key] = value
    return result


@dataclass(frozen=True, order=True)
class EntityRef:
    """A canonical gene or protein accession qualified by organism and unit.

    Registered parasite genes obey their reference accession format. Host gene
    formats and membership await their gene packs; protein accessions remain
    opaque and must be validated against the caller's source table.
    """

    organism: str
    kind: str
    identifier: str

    def __post_init__(self):
        _organism(self.organism)
        _text(self.kind, "entity kind")
        _text(self.identifier, "identifier")
        if self.kind not in {"gene", "protein"}:
            raise ValueError("Entity kind must be gene or protein")
        if self.kind == "gene" and self.organism in organisms.SPACES:
            if not organisms.get(self.organism).matches(self.identifier):
                raise ValueError("Gene accession does not match the explicit organism reference")
        if self.kind == "gene" and any(space.matches(self.identifier) and code != self.organism
                                       for code, space in organisms.SPACES.items()):
            raise ValueError("Gene accession belongs to another explicit organism reference")


@dataclass(frozen=True)
class BiologicalContext:
    """Explicit context axes; empty means unspecified, never an inferred default.

    Free-text axes preserve study-specific contexts beyond the current stage
    registry. A specified host must be a known host organism code.
    """

    stage: str = ""
    host: str = ""
    tissue: str = ""
    condition: str = ""
    strain: str = ""

    def __post_init__(self):
        for name, value in asdict(self).items():
            _text(value, name, optional=True)
        if self.host:
            _organism(self.host)
            if self.host not in organisms.HOST_TABLES and organisms.get(self.host).kind != organisms.HOST:
                raise ValueError("Context host must be a host organism")


@dataclass(frozen=True)
class Query:
    """One typed question with a biological context and requested output.

    Gene/protein questions have one endpoint; gene sets have one or more genes
    from the primary organism. Pairs have two ordered endpoints, with the first
    in the primary organism, allowing parasite/host questions. Target addresses
    have no endpoints; class addresses add one or more exact, case-sensitive
    values. A gene question may optionally name the target it asks about.
    """

    organism: str
    kind: str
    entities: tuple[EntityRef, ...] = ()
    target: str = ""
    values: tuple[str, ...] = ()
    context: BiologicalContext = field(default_factory=BiologicalContext)
    output: str = "inferences"
    schema_version: int = SCHEMA_VERSION

    def __post_init__(self):
        _organism(self.organism)
        _text(self.kind, "query kind")
        _text(self.output, "output")
        _text(self.target, "target", optional=True)
        if self.kind not in QUERY_KINDS or self.output not in OUTPUTS:
            raise ValueError("Unsupported query kind or output")
        if type(self.schema_version) is not int or self.schema_version != SCHEMA_VERSION:
            raise ValueError("Unsupported query schema version")
        if not isinstance(self.context, BiologicalContext):
            raise ValueError("context must be a BiologicalContext")
        if not isinstance(self.entities, tuple) or not all(isinstance(e, EntityRef) for e in self.entities):
            raise ValueError("entities must be a tuple of EntityRef addresses")
        if not isinstance(self.values, tuple):
            raise ValueError("values must be a tuple")
        for value in self.values:
            _text(value, "class value")
        if self.kind in {"label", "class", "trait"}:
            if not self.target or self.entities:
                raise ValueError("Label, class and trait queries need a target and no endpoints")
        else:
            expected = 2 if self.kind == "pair" else 1
            if self.kind == "gene_set":
                if not self.entities or any(e.kind != "gene" for e in self.entities):
                    raise ValueError("A gene set needs at least one gene")
                object.__setattr__(self, "entities", tuple(sorted(set(self.entities))))
            elif len(self.entities) != expected:
                raise ValueError(f"{self.kind} needs {expected} endpoint(s)")
            if self.kind in {"gene", "protein"} and self.entities[0].kind != self.kind:
                raise ValueError("Query and endpoint kinds must agree")
            if self.entities[0].organism != self.organism:
                raise ValueError("First endpoint must match the primary organism")
            if self.kind != "pair" and any(e.organism != self.organism for e in self.entities):
                raise ValueError("All gene-set endpoints must match the primary organism")
        if self.kind == "class":
            if not self.values:
                raise ValueError("A class needs at least one value")
            object.__setattr__(self, "values", tuple(sorted(set(self.values))))
        elif self.values:
            raise ValueError("Only class queries have class values")

    def to_dict(self) -> dict:
        """Return a JSON-compatible object with all context and schema fields."""
        result = asdict(self)
        result["entities"] = [asdict(e) for e in self.entities]
        result["values"] = list(self.values)
        return result

    def to_json(self) -> str:
        """Return deterministic JSON for this address, not an inference cache key.

        Artifacts must additionally declare code, sources, fitting population and
        settings under the later artifact contract.
        """
        return json.dumps(self.to_dict(), sort_keys=True, ensure_ascii=False, separators=(",", ":"))

    @classmethod
    def from_dict(cls, value: dict) -> Query:
        """Validate a complete schema object, refusing unknown or missing fields."""
        value = dict(_object(value, cls.__dataclass_fields__))
        value["context"] = BiologicalContext(**_object(value["context"], BiologicalContext.__dataclass_fields__))
        value["entities"] = tuple(EntityRef(**_object(e, EntityRef.__dataclass_fields__))
                                  for e in _array(value["entities"], "entities"))
        value["values"] = tuple(_array(value["values"], "values"))
        return cls(**value)

    @classmethod
    def from_json(cls, value: str) -> Query:
        """Parse versioned JSON without silently accepting duplicate object keys."""
        return cls.from_dict(json.loads(value, object_pairs_hook=_unique_object))


@dataclass(frozen=True)
class AliasRecord:
    """One verified alias-to-gene mapping with its kind and source artifact ID.

    The source is supplied by the caller (ideally a version/hash/row address).
    An index adapter cites its index artifact rather than inventing source rows.
    This records mapping provenance, not biological inference confidence.
    """

    alias: str
    entity: EntityRef
    kind: str
    source: str

    def __post_init__(self):
        for name in ("alias", "kind", "source"):
            _text(getattr(self, name), name)
        if not norm(self.alias):
            raise ValueError("Alias must contain an identifier")
        if not isinstance(self.entity, EntityRef) or self.entity.kind != "gene":
            raise ValueError("Alias records must address a gene")


@dataclass(frozen=True)
class Resolution:
    """A resolved, ambiguous or unresolved lookup retaining every mapping source."""

    organism: str
    query: str
    choices: tuple[AliasRecord, ...]

    def __post_init__(self):
        _organism(self.organism)
        if not isinstance(self.query, str):
            raise ValueError("Resolution query must be a string")
        if not isinstance(self.choices, tuple) or any(
                not isinstance(r, AliasRecord) or r.entity.organism != self.organism for r in self.choices):
            raise ValueError("Resolution choices must belong to the explicit organism")

    @property
    def status(self) -> str:
        """Distinguish no match, one canonical gene and multiple candidate genes."""
        count = len({record.entity for record in self.choices})
        return "unresolved" if count == 0 else "resolved" if count == 1 else "ambiguous"

    @property
    def entity(self) -> EntityRef | None:
        """Return the sole canonical entity, or None when a choice is still needed."""
        return self.choices[0].entity if self.status == "resolved" else None


class GeneResolver:
    """Resolve exact user-entered aliases inside one explicit organism.

    Unlike running-prose extraction, verified short or ordinary-word aliases may
    be supplied explicitly. Normalization follows the existing identity module;
    there is no substring matching, suffix guessing or cross-organism fallback.
    """

    def __init__(self, organism: str, records):
        _organism(organism)
        self.organism = organism
        self._lookup = {}
        for record in records:
            if not isinstance(record, AliasRecord) or record.entity.organism != organism:
                raise ValueError("Every alias record must belong to the explicit organism")
            self._lookup.setdefault(norm(record.alias), set()).add(record)

    def resolve(self, alias: str) -> Resolution:
        """Return all candidate mappings with sources; never select an ambiguity."""
        if not isinstance(alias, str):
            raise ValueError("Alias must be a string")
        choices = sorted(self._lookup.get(norm(alias), ()),
                         key=lambda r: (r.entity, r.kind, r.source, r.alias))
        return Resolution(self.organism, alias, tuple(choices))

    @classmethod
    def from_index(cls, index: GeneIndex, organism: str, source: str) -> GeneResolver:
        """Adapt an existing gene index, including its previously hidden collisions.

        The original index stores no mapping kind for collided aliases. Those
        choices are honestly tagged ``ambiguous_alias`` and cite the supplied
        index artifact. All canonical rows are checked against the organism.
        """
        _text(source, "source")
        entities = {gid: EntityRef(organism, "gene", gid) for gid in index.canonical}
        records = [AliasRecord(alias, entities[gid], kind, source)
                   for alias, (gid, kind) in index.lookup.items()]
        records.extend(AliasRecord(alias, entities[gid], "ambiguous_alias", source)
                       for alias, genes in index.ambiguous.items() for gid in genes)
        return cls(organism, records)
