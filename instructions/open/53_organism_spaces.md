# 53 · One space per organism, then links between pathogen, vector and host

**Status: open (designed 2026-09-26).**

Crash recovery: all implementation commits through `05650fd` survived and were pushed. The 45
focused recovery checks passed. Documentation CI passed; test CI had 4,138 passes and one stale
packaging assertion that disallowed subpackages. The assertion now checks that every Starplast
source package is included in the single distribution, including `starplast.spaces`.
All 50 publishing/release checks pass with this correction. Follow-up remote CI remains pending.

### Continued implementation

The user authorized sustained work on 2026-09-27, prioritizing token efficiency. See
`WORK_LOG.md` for completed milestones, pending dependencies and actual blockers. Current scope:
WP8's lossless per-space calibration publishing and WP1's versioned pack/build framework.
Acceptance: publishing a subset preserves untouched results/provenance; corrupt existing files
are not overwritten; a pack round-trips hashes, rejects tampering/unsafe paths, installs atomically,
and refuses to build over lost identifiers, columns or measured values.

Pf display fix implemented: per-space display recipes resolve 118 Pf features for all 5,720 genes.
`inference_spec` explicitly retains the recipe used by the shipped calibration. The executed
`notebooks/pf_display_layout_2026_09_27.ipynb` compares 17 target-specific/no-target strategy groupings
and both species' compatibility input matrices against pre-change fingerprints: all unchanged.
Only Pf coordinates and layout metadata changed; every edge array, both node tables, Tg layout,
slots and calibration were preserved. Focused checks: 156 passed. UI/docstring checks: 305 passed,
1 skipped offscreen; the real shared-OpenGL test then passed under Xvfb. UI/optimizer checks:
219 passed. Search/leakage/graph/strategy checks: 405 passed, 2 skipped. A final fingerprint check
confirmed all 17 groupings and both compatibility matrices unchanged. The clean 53.4 MB wheel
passes validation and package modules parse under Python 3.10 syntax. The earlier full-suite
GUI/GPU baseline failures remain documented below; no full-suite rerun or release was performed.

WP8 publication implemented: updates merge by organism/strategy, retain each result's original
sweep provenance and atomically replace the file. Corrupt existing files fail before writing.
Generated pages include all published spaces. Calibration/registry checks: 37 passed; shipped
calibration values were not changed. Remaining WP8 strategy applicability work is still open.

WP1 framework implemented: `packs.py` verifies complete ZIP manifests and SHA256 values, rejects
unsafe members and unreviewed/local-only license codes, downloads with bounded size and timeout,
and activates immutable versions atomically. `spaces.build_pack` and `scripts/build_space.py`
validate canonical unique IDs, graph order/indices and per-cell preservation against the previous
table. Registry paths resolve active packs without downloading or silently falling back. CLI and
usage are in `docs/space-packs.md`. 380 relevant checks passed and a clean 53.4 MB wheel includes
the new subpackage. Source-specific builders, acquisition and published download URLs remain open.

### Implementation pass, 2026-09-27

R3 implemented and validated with the baseline test failures recorded below. All 161 datasets now declare their organism;
the four host deposits declare Hs/Mm. The mixed protein cache is replaced by
`hs_host_proteins.parquet` (20,989 rows) and `mm_host_proteins.parquet` (15,590 rows). No identifier,
name, measurement or bridge endpoint was lost. The 54 identity-only rows were resolved from local
UniProt mapping files and human bridge sources, not guessed from names. The executed migration is
`notebooks/split_host_tables_2026_09_27.ipynb`; its manifest records checksums and coverage. The
script refuses unassigned/conflicting rows, unknown columns and overwriting later acquisitions.
Readers and builders select one host species at a time, including the slot viewer and generated
atlas. `slots.json` and the generated atlas remain byte-identical. Full host gene spaces, packs and
the Space menu remain subsequent work packages.

**Correction to the Pf follow-up:** the handoff's claim that strategies are unaffected is false.
`Context.blocks()` calls `embedding.default_spec()`. A proposed per-species recipe expands the
display from 20 source columns to 118 resolved features, but changes inference grouping from 50
to 62 blocks without a target, and changes it for every shipped calibration holdout. The experiment
is recorded in `notebooks/pf_layout_dependency_2026_09_27.ipynb`; production recipes/layouts were not
changed during that first audit. The subsequent display fix above explicitly separates display
and calibrated inference recipes. Do not reuse the old scores if the inference groupings change.

Validation (Python 3.12, pandas 3, offscreen Qt, local GPU environment):

* Focused host/deposit/dataset/organism/slot checks: 246 passed. The new migration and reader
  regressions also pass, including identity-only rows, zero/false measurements, ambiguous species,
  and refusing to overwrite later acquisitions. Leakage/search/source checks: 178 passed, 2 skipped.
* Full suite: **4,145 passed, 9 skipped, 7 failed** in 19 minutes. All seven failure cases reproduce
  on untouched `f02d075`: import-column counts, gated point sizes, persisted display settings,
  sprite caching, CPU UMAP fallback with cuML available, offscreen shortcut focus, and attention
  colors. These are existing GUI state/test-order and GPU-environment issues, not a green suite.
  A baseline run of `test_app_controls.py`, `test_display.py`, `test_build_graph.py`,
  `test_help_search.py`, and `test_app_smoke.py` with `--randomly-seed=42` reproduced six of them
  (11 failed, 356 passed). The seventh is reproduced by running
  `test_searching_an_exact_gene_id_selects_it` immediately before
  `test_a_gated_set_recedes_the_rest_of_the_map_without_recoloring_it` in `test_app_controls.py`
  with `-p no:randomly`: the selected gene retains its larger point size.
* `add_deposits.py --dry-run` found no loss. The migration checked every retained cell and parquet
  round trip; all 36,579 protein identifiers and all human bridge endpoints remain present.
* Regenerating the slot catalogue and atlas produced byte-identical files. Parasite node tables,
  graphs and calibration are unchanged. Both species' slot audits have zero orphan columns,
  double ownership or missing citations.
* A clean wheel builds and passes `check_wheel.py` (53.4 MB). Negative checks reject a wheel with
  the old mixed cache or a missing mouse cache. `release.py check`, Python 3.11 syntax parsing and
  `git diff --check` pass. No release/version bump is part of this milestone.

* **R0 done:** `starplast/organisms.py` declares Tg and Pf, and `tests/test_organisms.py` holds it to every literal it replaces.
* **R1 done:** these now read the registry:
  * `slots.SPECIES_TABLES`, `SPECIES_PREFIXES` and `SPECIES_BRIDGE_TABLES`;
  * `app.SPECIES` and `DEFAULT_SPECIES`, and the gene record link;
  * `strategies.Context.shipped`, its graph file, `other()` (the space's partner) and `_guess_organism` (registry detection first);
  * the calibration `TARGETS`/`NUMBERS`;
  * the slot generator's `STAGES_BY_ORGANISM` (regenerated `slots.json` byte-identical).
* **R2 started:** `test_no_file_gains_an_organism_literal` is a ratchet. It measured 156 bare "Tg"/"Pf" literals in 32 files, and any file that gains one fails. Lower the numbers as files move to the registry.
* **R3 implemented:** explicit dataset organisms and separate Hs/Mm protein tables (see above).

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

### Raw data acquired (2026-09-26)

5.9 GB is downloaded under `<STARPLAST_DATA>/spaces/<code>/<kind>/<source>/`.

* **What is there:** Hs, Mm, Rn, Fc, Ag, As, Cp and Pb, plus `shared/` (gene2pubmed, gene_orthologs,
  BioGRID, BioGRID PTMs and ORCS, IntAct).
* **How each deposit is recorded:** URLS.txt holds the resolving query per file, and SHA256SUMS
  verifies it (0 mismatches).
* **Per space:** a `MANIFEST.json` with release, licence code, sanity check and target slot family.
* **Also there:** `spaces/README.md`, and the scripts and an annotated notebook under
  `spaces/_acquisition/`.
* **Structure:** mean pLDDT per model from the AlphaFold DB search API for every space. There is no
  per-proteome summary file.
* **Hs highlights:** Ensembl 116 plus UniProt, GTEx v10, HPA 25.1 (HPA is CC BY 4.0 now, not SA),
  Bgee 15.2, DepMap 24Q4 (from figshare+; the portal is behind Cloudflare), COMPARTMENTS, OpenCell,
  STRING 12.0, HuRI, BioPlex, CORUM (non-commercial), PaxDb, iPTMnet, gnomAD v4.1, GOA, and 9
  infection series (Tg, Pf, Pv, Pb, Cp).
* **Mm:** Tabula Muris Senis bulk, IMPC DR24 viability and phenotypes, ImmGen, and 9 infection
  series.
* **Rn:** BodyMap, RatGTEx.
* **Fc:** cat small-intestine epithelium transcriptome and phosphoproteome under Tg (PMIDs 37491273,
  38003154).
* **Ag:** MozAtlas, blood-meal time course, PRIDE salivary gland/hemolymph/saliva, 5 infection
  series.
* **Cp:** Walzer 2024 single-cell atlas, oocyst proteomes. No genome-wide screen exists.
* **Pb:** PlasmoGEM Bushell 2017 and Stanway 2019 (the real tables), Russell 2023, Malaria Cell
  Atlas.

Traps it found:

* **VEuPathDB is closed to anonymous downloads.** Release 71 public downloads return 404 and the
  service returns 401, so the Cp/Pb/Ag/As references come from Ensembl Genomes r63.
* **Two PlasmoGEM files in the old dataset tree are fakes.** They are 1.8 kB HTML "preparing to
  download" pages, and `reference/plasmodb/pb_transfer/` does not exist. The real tables are under
  `spaces/Pb/essentiality/`.
* **The FANTOM5 path `host/fantom5/24670764/` does not exist.** The file is now in
  `spaces/Mm/expression/FANTOM5_E-MTAB-3579/`.
* **Several id spaces are split or retired:**
  * Cat: three id spaces (retired Felis_catus_9.0, and ENSFCTG).
  * *A. gambiae*: NCBI moved to AGAMI1_ ids with no AGAP cross-reference; bridge through UniProt.
  * *A. stephensi*: four id spaces.
  * Rat: STRING, COMPARTMENTS and PaxDb use old ENSRNOP ids.
  * Cp: no open IOWA-ATCC↔cgd map.
  * Pb: some tables use old 6-digit PBANKA ids.
* **gene2pubmed barely covers the parasites** (Cp: 29 PMIDs).
* **Junk to delete by hand:** `spaces/As/abundance/PXD001647/` holds 48 kB of unrelated Hydra and
  Daphnia mzTabs, recorded as `rejected_unusable`.

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

* **The Pf display previously used only 20 features.** The display follow-up now uses its own
  slots (145 source columns, 118 resolved features). `inference_spec` preserves the calibrated
  grouping separately: broad raw-column coverage never proved that grouping changes were harmless.
  The two executed notebooks above record both the dependency and the final invariance checks.

* R3 resolved the mixed human/mouse host cache. Host readers and deposits must continue to select
  one explicit organism; full host gene spaces remain pending.
* `calibration.write` now merges measured organism/strategy entries. Unswept entries keep their
  original provenance; changed strategy recipes still require a fresh calibration.
* HPA consensus integrates GTEx; BioGRID ORCS and OGEE reuse DepMap; STRING's combined channel
  contains BioGRID and IntAct. These are one family each: holding one out must hold out the others.
* The Pf host contexts in `generate_slot_table.py` already name Anopheles midgut and salivary gland.
  Those become bridge slots into the Ag space, not Pf columns.
