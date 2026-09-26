"""The organism registry: one declaration per species space, instead of literals in fifty files.

Starplast grew from one organism to two by copying: `slots.SPECIES_TABLES`, `app.SPECIES`,
`strategies.Context.shipped`, the calibration targets, the slot generator's stages and the
localization abbreviations each carry their own `"Tg"`/`"Pf"` table. Adding the next species -- a
host, a vector, another apicomplexan -- would have meant editing every one of them, and missing one
fails silently (several lookups default to Toxoplasma). This module is the one place a space is
declared; everything else is to read from it (instruction 53, step R0 onward).

A space is one species' own gene table, graph, slots and calibration. Spaces are never merged:
orthology and host-pathogen interfaces are bridges between them, declared per space.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field

PARASITE, HOST, VECTOR = "parasite", "host", "vector"


@dataclass(frozen=True)
class Space:
    """One species' space: what it is, how its genes are named, and where its tables live."""
    code: str                   # short, unique; slot keys are f"{code}_{name}"
    species: str                # binomial, as the menus show it
    reference: str              # strain or assembly the gene ids refer to
    kind: str                   # parasite, host or vector
    gene_regex: str             # full match for a canonical gene id in the node table
    prefixes: tuple             # accession prefixes that identify the table (first = canonical)
    database: str               # the genome database that curates it
    record_url: str             # a gene's record page; {id} is replaced by the gene id
    nodes: str                  # node table file in the data directory
    graph: str = ""             # edge layers
    host_bridges: str = ""      # host-pathogen pairs whose parasite end is in this space
    partner: str = ""           # the space `Context.other()` crosses to by orthology
    orthomcl: str = ""          # OrthoMCL / VEuPathDB four-letter abbreviation
    contexts: frozenset = frozenset()   # life stages (parasite) or tissues (host, vector)
    targets: tuple = ()         # categorical columns calibration holds out
    numbers: tuple = ()         # numeric columns calibration predicts
    hosts: tuple = ()           # codes of host spaces
    vectors: tuple = ()         # codes of vector spaces
    distribution: str = "wheel"  # shipped in the wheel, or a downloadable pack

    def record(self, gene_id: str) -> str:
        """The URL of one gene's record in its genome database."""
        return self.record_url.format(id=gene_id)

    def matches(self, gene_id: str) -> bool:
        """Whether `gene_id` is a canonical id of this space."""
        return re.fullmatch(self.gene_regex, str(gene_id)) is not None


SPACES: dict = {}


def register(space: Space) -> Space:
    """Add a space. A code used twice, or a kind outside parasite/host/vector, is refused."""
    if space.code in SPACES:
        raise ValueError(f"space {space.code!r} registered twice")
    if space.kind not in (PARASITE, HOST, VECTOR):
        raise ValueError(f"{space.code}: kind must be parasite, host or vector")
    re.compile(space.gene_regex)
    SPACES[space.code] = space
    return space


register(Space(
    code="Tg", species="Toxoplasma gondii", reference="ME49", kind=PARASITE,
    gene_regex=r"TGME49_\d{6}", prefixes=("TGME49_", "TGGT1_"), database="ToxoDB",
    record_url="https://toxodb.org/toxo/app/record/gene/{id}",
    nodes="nodes.parquet", graph="graph.npz", host_bridges="host_bridges.parquet",
    partner="Pf", orthomcl="tgon",
    contexts=frozenset({"tachyzoite", "bradyzoite", "sporozoite", "oocyst", "merozoite", "sexual",
                        "enteric", "gametocyte"}),
    targets=("compartment", "lopit_unified", "dtm_class", "screen_any_phenotype",
             "stage_enriched_derived"),
    numbers=("fit_invitro_hff", "fit_invivo_PE", "expr_tachy")))

register(Space(
    code="Pf", species="Plasmodium falciparum", reference="3D7", kind=PARASITE,
    gene_regex=r"PF3D7_(?:\d{7}|API\d{5}|MIT\d{5})", prefixes=("PF3D7_",), database="PlasmoDB",
    record_url="https://plasmodb.org/plasmo/app/record/gene/{id}",
    nodes="pf_nodes.parquet", graph="pf_graph.npz", host_bridges="pf_host_bridges.parquet",
    partner="Tg", orthomcl="pfal",
    contexts=frozenset({"ring", "trophozoite", "schizont", "gametocyte", "ookinete", "sporozoite",
                        "merozoite", "asexual blood stage", "sexual", "liver stage", "oocyst"}),
    targets=("pb_transferred_phenotype", "stage_enriched_derived", "is_exported"),
    numbers=("piggybac_mis", "expr_schizont", "mean_plddt")))


# --------------------------------------------------------------------------- lookups
def get(code: str) -> Space:
    """One space by code, or a KeyError naming the ones that exist."""
    if code not in SPACES:
        raise KeyError(f"no space {code!r}; there are {', '.join(SPACES)}")
    return SPACES[code]


def codes(kind: str | None = None, available: bool = False) -> list:
    """Space codes in registry order, optionally of one kind or only those whose table is present."""
    return [c for c, s in SPACES.items() if (kind is None or s.kind == kind)
            and (not available or os.path.exists(nodes_path(c)))]


def by_species(name: str) -> Space:
    """The space whose species is `name` (as the menus show it)."""
    for s in SPACES.values():
        if s.species == name:
            return s
    raise KeyError(f"no space for species {name!r}")


def nodes_path(code: str) -> str:
    """Where a space's node table lives (honouring STARPLAST_CACHE, as every table does)."""
    from . import paths
    return paths.cache_file(get(code).nodes)


def graph_path(code: str) -> str:
    """Where a space's edge layers live."""
    from . import paths
    return paths.cache_file(get(code).graph)


def table_map() -> dict:
    """code -> node table file: what `slots.SPECIES_TABLES` holds."""
    return {c: s.nodes for c, s in SPACES.items()}


def prefix_map() -> dict:
    """code -> accession prefixes: what `slots.SPECIES_PREFIXES` holds."""
    return {c: s.prefixes for c, s in SPACES.items()}


def detect(gene_ids, sample: int = 500) -> str | None:
    """The space whose canonical ids MOST of `gene_ids` fully match; None if none reaches half.

    A majority of full matches rather than the first prefix that appears anywhere: one stray
    accession from another species must not decide what a whole table is.
    """
    ids = [str(g) for g in list(gene_ids)[:sample]]
    if not ids:
        return None
    best, share = None, 0.0
    for code, space in SPACES.items():
        rx = re.compile(space.gene_regex)
        s = sum(bool(rx.fullmatch(g)) for g in ids) / len(ids)
        if s > share:
            best, share = code, s
    return best if share > 0.5 else None
