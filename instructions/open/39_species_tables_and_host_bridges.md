# 39 — One table per species, host tables beside them, and bridges rather than merges

**Status: open. Design settled 2026-08-15 with the user; implementation not started.**

## The shape

Three kinds of table and two kinds of bridge. Nothing is merged; everything that crosses a boundary
crosses it through a slot that can be inspected, held out, and disbelieved.

| kind | tables | keyed by | `unit` |
|---|---|---|---|
| parasite | `Tg`, `Pf`, `Pb`, `Pv`, `Pk`, `Py` | parasite gene | `gene` |
| host | `Hs` (human), `Mm` (mouse), `Ag` (*Anopheles gambiae*) | host gene | `host_gene` |
| bridge · host | one per (parasite species, host species) | gene pair | `pair` |
| bridge · orthology | one per parasite species pair | gene pair | `ortholog_group` |

`unit` already exists in the catalog and already takes `gene` and `pair`. `host_gene` and
`ortholog_group` are new values, not new machinery.

**The `Toxo_` -> `Tg_` rename is DONE (2026-08-16).** The prefix is now an organism code like every
other one, so `Tg_transcription_tachyzoite` sits beside `Pf_transcription_asexual_blood_stage` and
the host and species codes below can join them without one odd name out. It touched nine Python
files, the catalog, and the three generated tables; `Toxo` survives only inside the word
*Toxoplasma*, which is prose rather than an identifier.

## Why per-species, and why orthology is a bridge

Gene identifiers do not align (`PF3D7_*`, `PBANKA_*`, `TGME49_*`) and neither does the data:
piggyBac essentiality is *falciparum*, PlasmoGEM barcoded knockouts and most liver-stage work are
*berghei*, relapse biology is *vivax*, field variation is *falciparum* and *vivax*. Merged into one
table the NaN pattern would encode **which species was convenient to work on**, and since
missingness is already a feature here (`na_policy="indicator"`) the embedding would cluster on
laboratory history and it would look like biology.

There is already a precedent for the right answer in this codebase: `localization · transferred`
carries orthoLOPIT's cross-species transfer as its own slot with its own `target_family`, rather
than quietly filling the measured column. Cross-species transfer reuses that exactly:
`Pf_fitness · transferred from Pb` is a slot, marked as a transfer, `derived_from` the *berghei*
column it came from.

### The leakage rule, which is the whole reason this is a bridge

**`target_family` closure must span species.** Fitness is one family whether the number came from
*Pb* or *Pf*. Transfer *berghei* fitness onto *falciparum*, hold out *falciparum* fitness, and
"recover" it, and you have measured orthology, not biology — the same circularity this project has
already published once and then corrected. Concretely:

* every transferred column carries `derived_from` naming its source column **and** source species;
* `search.excluded_for` must exclude the whole family across every species, not per table;
* a finding whose evidence is a transferred column is marked `circular` unless the target species
  has its own measurement of that family;
* a test must exist that fails if a transfer survives its own family's hold-out.

## Host tables

Rows are host genes. These do not enter a parasite embedding as feature columns and must not: a
parasite map whose rows were partly host genes would break every claim of the form "cluster 5 is
71% IMC", which assumes rows are parasite genes. They earn their place by answering the question the
map keeps raising and cannot settle — *is the host partner this parasite protein binds even present
in the tissue this stage lives in?*

**One slot per host species.** Not a "warm-blooded host" fill: for *Plasmodium* the mosquito and the
human are sequential hosts, not alternatives, and averaging them answers neither.

### Host-context nodes (the `context` hierarchy already exists — hang these off it)

    Tg_tachyzoite      / human fibroblast, human monocyte, mouse BMDM, mouse brain
    Tg_bradyzoite      / mouse brain, mouse muscle, human neuron
    Pf_merozoite       / human erythrocyte
    Pf_ring/troph/schizont / human erythrocyte
    Pf_gametocyte      / human bone marrow, human erythrocyte
    Pf_ookinete        / Anopheles midgut
    Pf_oocyst          / Anopheles midgut basal lamina
    Pf_sporozoite      / Anopheles salivary gland, human dermis, human hepatocyte
    Pb_liver stage     / mouse hepatocyte
    Pv_hypnozoite      / human hepatocyte

### Host slot families, per context node

* `host proteome · <tissue>` — `unit=host_gene`, policy `average`
* `host transcriptome · <tissue>` — `unit=host_gene`, policy `average`
* `host surface / receptor repertoire · <tissue>` — `unit=host_gene`, policy `fill`; this is what a
  binding claim gets checked against
* `host response to infection · <tissue>` — `unit=host_gene`, policy `separate`; the disagreement
  between infected and uninfected is the signal, so it must not be averaged away

### Where the host data comes from

Named as resources rather than accessions, because accessions must be **resolved** by the candidate
search and never typed from memory — the rule the last handoff set, which held for all 193 existing
citations.

| host context | resource |
|---|---|
| human tissue transcriptome / proteome | Human Protein Atlas; GTEx; PRIDE for tissue proteomes |
| human erythrocyte proteome | PRIDE — red-cell proteome and membrane proteome studies |
| human hepatocyte | HPA liver; primary hepatocyte RNA-seq in GEO/ArrayExpress |
| human dermis / skin | HPA skin; GTEx skin |
| mouse tissue | Tabula Muris; mouse ENCODE; PRIDE mouse tissue proteomes |
| *Anopheles* midgut, salivary gland | VectorBase; PRIDE salivary-gland and midgut proteomes |
| host genome / annotation | Ensembl, UniProt (host gene identifiers and orthology) |

## Bridge slots (parasite ↔ host), `unit=pair`

One family, several assays, each its own slot because they are not comparable:

* `proximity labelling at the PV/PVM` — TurboID/BioID; the study the user cited belongs here
* `host IP-MS / co-IP` — parasite bait, host prey
* `crosslinking MS across the interface`
* `dual perturb-seq` — parasite gene perturbed, host transcriptome read (already exists as
  `host transcriptional effect per effector`; make it a bridge slot rather than a parasite-gene slot)
* `curated parasite–host PPI` — literature-curated, lower tier, `fill` behind the measured ones
* `receptor–ligand binding` — direct assays (e.g. erythrocyte invasion receptor work)

A pair row is `(parasite species, parasite gene, host species, host gene, evidence, source)`. The
species belong in the key: the same parasite protein binding a mouse and a human orthologue are two
observations, and which one was made is exactly what a reader needs.

## Orthology bridge, `unit=ortholog_group`

* Source: OrthoMCL / VEuPathDB ortholog groups; keep the group id, not just a best hit.
* One slot per ordered species pair: `Pf_* transferred from Pb`, `Tg_* transferred from Nc`, etc.
* Record one-to-many explicitly. A *falciparum* gene with three *berghei* orthologues is not a
  transfer, it is a question, and averaging across the three silently invents a number.

## Acquisition tiers

Define slots for every species — an empty slot is informative, that is the point of the table — but
scope **acquisition** as:

1. **Now:** *P. falciparum* (3D7) and *P. berghei* (ANKA). Between them they carry nearly all the
   genome-scale data, and they are the pair where transfer is most useful: *berghei* has the fitness
   screens, *falciparum* has the clinical relevance.
2. **Next:** *P. vivax* (relapse, field variation), *P. knowlesi* (in vitro tractable).
3. **Hosts:** human and mouse first; *Anopheles* with the vector species chosen deliberately —
   *gambiae* and *stephensi* are used by different laboratories and their proteomes are not
   interchangeable, so the slot must name which.

## Acceptance

* Column partition holds per table: no column claimed by two slots, none claimed by none — the check
  that already passes for Toxoplasma, run per species.
* A test that a transferred column does not survive its own family's hold-out.
* A test that a pair row without both species named is refused.
* Host tables never appear as feature columns in a parasite embedding: assert it.
* 100% coverage, no `pragma`. Coverage is currently 99% and should be back at 100 before this lands
  on top of it.
