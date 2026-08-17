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

## The data for the first host bridge is downloaded and verified (2026-08-16)

`interaction · with host proteins` is the last unfilled pair slot in the Toxoplasma arm, and it is
blocked on this instruction's bridge rather than on acquisition. The dataset is in place:

**PXD016383, MYR1 immunoprecipitations**, `proteinGroups.txt` at
`datasets/quarantine/2026_08_16_pride/Tg/host_interaction/`. 1,092 protein groups, LFQ intensities
for two MYR1 IPs (M1, M2) and two controls (R1, R2). 425 rows name a Toxoplasma accession and 674 a
human one, which is what a parasite–host bridge needs: one experiment, both ends.

**Verified on the parasite side, which is the side that can be checked today.** Taking
log2(mean LFQ in MYR1 IP over mean LFQ in control), MYR1 is **rank 1 of 325** at +34.4 — the bait
tops its own pulldown — with MYR3 at rank 47 and GRA44, GRA7, GRA9, GRA52 and GRA50 in the top ten.
That is MYR1's known neighbourhood at the vacuole membrane, so the IP worked and the enrichment is
the right way round.

### CORRECTION: the host side is not verifiable either, so this is not only structural

The paragraph above says the parasite side verifies, and it does. The host side does not, and
claiming this slot was "blocked on architecture, not data" was wrong. Reading it properly:

* **Contaminants hide inside protein groups.** The four most enriched host groups are keratins.
  MaxQuant prefixes a group with `CON__` only when the LEADING entry is a contaminant, and keratin
  arrives in the middle of a group headed by `sp|`. Filtering on the prefix keeps them; filtering on
  `CON__` appearing ANYWHERE in the group drops 26 groups and all of the keratin.
* **After filtering** — no contaminant entry anywhere, and at least two unique peptides in BOTH MYR1
  replicates — 219 host groups remain and 112 are enriched over the control.
* **Those 112 are led by tubulin, filamin C, HSP90, SERCA2, ribophorin and PDCD6**: abundant
  cytoskeletal and ER proteins, which is what an IP background looks like as much as what a vacuole
  translocon's neighbourhood looks like.

**There is no known-positive host partner of MYR1 to test the list against.** The parasite side could
be checked because MYR1 must top its own pulldown and MYR3 must be near it; the host side has no
equivalent, and enrichment over one control IP does not separate a specific partner from an abundant
protein that sticks.

So a host table built from this deposit would carry an unverifiable layer into the map. The slot
needs either a host interactome with its own controls and a checkable positive, or a bridge whose
evidence is something other than one IP.

### What is still needed structurally, once the data question is settled

1. **A human gene table.** The 674 host rows carry UniProt accessions and the deposit's
   `Fasta headers` column is empty, so mapping them to genes needs UniProt's ID mapping — a new
   external source, and the first identity layer in this project that is not ToxoDB's.
2. **A bridge that is a pair across two tables.** `slots.is_filled` now looks in the table matching a
   slot's unit (see `UNIT_TABLES`, added for the metabolite table), but an edge whose two ends live
   in different tables is not an edge in `graph.npz` and has no representation yet.
3. **Filling semantics for a bridge slot** — what `coverage` means when the denominator is pairs
   drawn from two tables of different size.

The metabolite table built on 2026-08-16 is the worked example for step 1's shape: a second table
keyed by its own identifier, reached through slots that declare their unit, graded on its own
denominator. Steps 2 and 3 are genuinely new.
