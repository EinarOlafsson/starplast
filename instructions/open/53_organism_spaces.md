# 53 · One space per organism, then links between pathogen, vector and host

**Status: open (designed 2026-09-26; WP0 waits for the data-audit-b and UI branches to merge).**

The user, 2026-09-26: "we should also make the host and vector datasets more comprehensive!
(plasmodium, cryptosporidium, human, mouse, feline, anophelus, rat, datasets and informationslots
need to be added) each in its own space first then we can start connecting pathogen to vector and
host."

## STATE

Measured on nightly @446e20e.

* **There is no organism registry.** About 150 `"Tg"`/`"Pf"` literals sit in 50 files, and at
  least five separate tables map species to files or ID prefixes. Three naming schemes coexist:
  `Tg/Pf`, `tgon/pfal/cpar` (`localization.py:34`), and `"Toxoplasma gondii"` (`app.py:277`).
* **"host" is two species in one table.** `host_proteins.parquet` has 36,579 rows keyed on UniProt
  accession. Its human columns (fibroblast_tpm, 19,087 rows) and mouse columns (brain_tpm 9,763,
  bmdm_tpm 15,437) never overlap. That breaks instruction 39's rule of one table per species.
* **Host slots belong to the parasites.** `slots.json` has 290 slots. Tg has 116 gene, 26 host_gene,
  10 pair and 3 metabolite slots; Pf has 99, 24, 9 and 3. So "Tg_host transcriptome · human
  fibroblast" is a Toxoplasma slot.
* **Two lookups silently default to Tg.**
  * `datasets.organism_of` (`datasets.py:2670`) returns "Tg" for every `host_*` dataset.
  * `strategies._guess_organism` treats anything that is not PF as Tg.
* **How Pf was added:** a separate builder (`plasmodium.build_all`) plus `pf_graph.py` copying Tg's
  graph construction, mirrored slots (`_pf_mirror`), a guarded rebuild (`build_plasmodium.py`), and
  ID-space tests (`test_pf_graph.py`, `test_plasmodium.py`). Every new space copies this template.
* **Public data only.** The VEuPathDB service API needs a key we do not have. The public
  `common/downloads/Current_Release` files need none; `archive.discover_veupathdb` already walks
  them.

Where each organism assumption lives (file:line):

| Concern | Where |
|---|---|
| Table/prefix maps | `slots.py:44-66`, `table_organism` :69 (first prefix match wins) |
| App | `app.py:277-310` SPECIES/load; Species menu :1498; `open_species` :2644; record link :4243; `EDGE_TYPES` :60, `COLOR_MODES` :83 |
| Strategies | `strategies.py` `shipped`, `_guess_organism`, `Context.shipped`, `family_of`, graph file, `other()` (the other one of Tg/Pf), synthetic Pf partner |
| Orthology | `strategy_catalog.py` `_source_default`, `_through_orthologs` (joins on `orthogroup`; Tg and Pf share 2,526 OG6 groups) |
| Leakage | `search.py` `SAME_QUANTITY`, `LAYER_SOURCES` (keyed on column name only); `leakage.py` |
| Slot generator | `generate_slot_table.py` `ORGANISMS`, `HOST_CONTEXTS_TG/PF`, `HOST_FAMILIES`, `HOST_COLUMNS`, `STAGES_BY_ORGANISM`, `all_slots`, `_rows` |
| Datasets / deposits | `Dataset` has no organism field; `organism_of` works from the key prefix; `Deposit.organism` is one of {Tg, Pf, host} |
| Host | `host.py` table names, `uniprot_index` (human/mouse), GTEx, FANTOM5, `TISSUE_REFERENCES` |
| Calibration | `calibrate_strategies.py` TARGETS/NUMBERS/`--organism`; `calibration.write` **overwrites the whole file** |

## WHY IT MATTERS

Adding a species today means editing about 30 files. Every new species multiplies the hard-coded
literals and the chances for silent Tg defaults. Host and vector biology can only become inference
targets once each species has its own table, slots, leakage families and calibration.

## WHAT TO DO

### WP0, registry and refactor (serial; everything else waits for it)

Add `starplast/organisms.py` with one `Space` per organism. A space declares:

* `code` (Tg, Pf, Pb, Cp, Hs, Mm, Rn, Fc, Ag, As), `species`, `reference`, and `kind`
  (parasite, host or vector);
* `gene_regex`, `alias_prefixes`, `id_source`, `identity_file` and `symbol_prefix`;
* the file paths for `nodes`, `graph` and `mentions`;
* `orthology` and `orthomcl_abbrev`;
* `hosts`, `vectors` and `partner` (Tg↔Pf, Pb→Pf);
* `contexts` (stages or tissues), `record_url` and `pubmed_query`;
* `edge_labels` and `colour_modes`;
* `calibration` tier, `distribution` (wheel or pack) and `builder`.

Helpers: `get`, `spaces(kind)`, `detect(ids)` (a majority full-match vote), `nodes_path`,
`graph_path` and `display_name`.

Refactor in steps. The suite stays green after each one.

* **R0.** The registry for Tg and Pf only, plus `tests/test_organisms.py`. The test asserts that the
  registry reproduces every current literal.
* **R1.** Replace each constant with a derived alias that keeps its old name.
  `Context.other()` → `organisms.get(code).partner`.
* **R2.** Tests read the organism set from the registry. Add a lint test that allowlists the
  remaining `"Tg"`/`"Pf"` literals, and shrink the allowlist over time.
* **R3.** Add `Dataset.organism`. Split `Deposit.organism` "host" into "Hs" and "Mm".
* **R4.** The slot generator reads the registry. Acceptance: `slots.json` is byte-identical.

### After WP0, in parallel worktrees (owned files do not overlap)

| WP | What | Owns | Accept |
|---|---|---|---|
| 1 | Space framework and on-demand packs | `spaces/__init__.py`, `spaces/veupath.py`, `packs.py`, `paths.py`, `scripts/build_space.py`, `check_wheel.py` | A pack round-trips its sha256; a tampered pack is refused; a lost column refuses the build |
| 2 | Slot generator per space | `generate_slot_table.py`, `scripts/slot_spaces/*.py` | Tg/Pf output byte-identical; each new space emits prefixed slots |
| 3 | Human (Hs), mouse (Mm) | `spaces/vertebrate.py`, `hs.py`, `mm.py`, `dataset_registry/hosts.py` | ENSG/ENSMUSG unique; ≥95% of reviewed UniProt resolved; GTEx matches today's values |
| 4 | Rat (Rn), cat (Fc) | `spaces/rn.py`, `fc.py` | Every column Fc takes through orthology carries `derived_from`; each empty slot has a verdict |
| 5 | Anopheles gambiae (Ag), stephensi (As) | `spaces/anopheles.py`, `dataset_registry/vectors.py` | AGAP IDs; As assembly recorded; blood-meal slots `separate` |
| 6 | Cryptosporidium (Cp), P. berghei (Pb), then P. vivax (Pv) | `spaces/crypto.py`, `pberghei.py`, `dataset_registry/apicomplexa.py` | OG6 present; Cp LOPIT loaded; Pf `pb_transferred_*` equal to the projection of the Pb space |
| 7 | Leakage per space | `search.py` families, `LAYER_DERIVED_FROM`, `leakage.py` | Attack tests: GTEx→HPA, BioGRID→STRING, DepMap→ORCS, cross-space transfer |
| 8 | Strategies across N spaces | `strategies.py` (`partner`, `bridge`, `applicable`), `calibration.py` (merge per space) | Publishing one space leaves the others byte-identical; strategies that do not apply are recorded as N/A |
| 9 | UI | `app.py` Space menu (Parasites / Hosts / Vectors, a download for any space not installed) | Every installed space opens offscreen |

**Phase 2** follows once the spaces exist:

* orthology bridges (OG6, Ensembl Compara, OrthoDB; 1:n relations kept explicit);
* host–pathogen PPI re-keyed to `(space_a, id_a, space_b, id_b)`, using HPIDB and the IntAct subset;
* vector–parasite links (stage↔tissue);
* `Context.bridge(code, kind)`;
* a "Linked spaces" dock that carries a selection across windows.

### Data per space

Resolve every identifier through an API; never type one. Licence codes: **O** open, **SA**
share-alike, **V** verify the terms first, **X** local build only.

* **Hs.**
  * GTEx v10 (O). HPA tissue, cell type and subcellular (SA). DepMap summaries (O).
  * UniProt/GO, with experimental codes only for targets (O).
  * STRING, BioGRID and IntAct (O). HuRI, BioPlex, OpenCell and CORUM (V). AlphaFold (O).
  * iPTMnet (V). **PhosphoSitePlus (X)**. gnomAD constraint (V). Bgee (O). gene2pubmed (O).
* **Mm.**
  * Tabula Muris Senis (O). FANTOM5, already present. ENCODE (O). ImmGen (V).
  * IMPC (O). MGI (V). BioGRID ORCS (O). STRING. AlphaFold.
  * Tg and Pb infection series (GEO).
* **Rn.** Rat BodyMap (O). RatGTEx (V). Bgee. RGD (V). STRING. AlphaFold. Tg-infected brain (GEO).
* **Fc.** Little measured data. Use Bgee or Expression Atlas where present, AlphaFold and STRING.
  Look for feline intestinal organoid series with Tg sexual stages. Most slots will be transferred
  and declared empty, each with a verdict.
* **Ag.**
  * Tissue atlases. Blood-meal time courses.
  * Plasmodium-infected midgut and salivary gland (GEO).
  * Midgut and salivary proteomes (PRIDE). Ag1000G (V). Insecticide resistance (V). AlphaFold.
* **As.** Lab transmission data (far less).
* **Cp.**
  * CryptoDB IOWA-ATCC, with `cgd` aliases. OG6. LOPIT, already in the tree.
  * Life-cycle scRNA and time courses (GEO). Oocyst and sporozoite proteomes (PRIDE).
  * There is no known genome-wide CRISPR screen; search for one before assuming.
* **Pb.** PlasmoGEM screens get a native home instead of existing only as transfers into Pf.
* **Pf (expand).** Malaria Cell Atlas (V). MalariaGEN summaries (V). PRIDE stage proteomes.

**Size.** The wheel is 52.6 MB (PyPI limit 100 MB).

* **In the wheel:** Tg, Pf, Cp and Pb (under ~65 MB).
* **As packs** in `user_cache_dir()/spaces/<code>/<version>/` with a sha256 manifest and per-column
  licences: Hs (~15-25 MB), Mm, Rn, Fc, Ag, As and Pv.

**Calibration.**

* **Tier A (full sweep):** Tg, Pf and Hs (Hs capped at ~8k rows).
* **Tier B (3 seeds, defaults + best):** Mm, Pb, Cp and Ag.
* **Tier C (self-test only):** Rn, Fc, As and Pv.

## HOW TO KNOW IT WORKED

* WP0: the full suite passes; `slots.json` is byte-identical; the registry reproduces every old
  literal; the lint test passes.
* Each space: the ID-space tests from `test_pf_graph.py` carried over; each slot filled, or empty
  with a verdict; its leakage attack tests pass; its self-tests and calibration published without
  touching any other space's numbers.

## TRAPS

* The host table stacks human and mouse. Splitting it changes `host_columns` and every deposit with
  organism "host".
* `calibration.write` overwrites the whole file. Merge per space before calibrating a third one.
* HPA consensus integrates GTEx; BioGRID ORCS and OGEE reuse DepMap; STRING's combined channel
  contains BioGRID and IntAct. These are one family each: holding one out must hold out the others.
* The Pf host contexts in `generate_slot_table.py` already name Anopheles midgut and salivary gland.
  Those become bridge slots into the Ag space, not Pf columns.
