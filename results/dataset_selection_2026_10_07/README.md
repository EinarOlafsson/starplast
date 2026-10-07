# Dataset-selection audit — 7 October 2026

The citation-rate preference is implemented and persisted for future dataset
admissions. This audit resolves available bibliometrics and searches every slot;
**it does not establish that every installed dataset is the best in the literature**.
Origin-paper/deposit validation and eligible replacements remain instruction 65.04.

| Audit measure | Result |
|---|---:|
| Registered sources accounted for | 162 / 162 |
| Publication identities resolved from recorded PMID/DOI | 95 |
| Remaining source roles/identity gaps | 67 |
| Slot declarations searched | 297 / 297 |
| Unique scoped queries / primary search requests | 216 / 432 |
| Failed literature requests | 0 |
| Retained publication/slot/order records | 5,058 |
| Distinct retained publication identities | 2,032 |
| Slot views with truncated searches | 231 |
| Slot views with no query hits | 28 |
| Diagnostic quantity/context flags triaged | 6 |
| Admitted replacements / promoted datasets | 0 / 0 |

The 67 remaining entries are **not 67 failed downloads**: 13 are local/derived
sources without their own publication, 19 reference resources lack resolved
publication/version lineage, 25 other sources lack resolved publication identity,
six accession text searches need origin-paper review, and four accession searches
found no publication in this index. Known zero citations remain distinct from
unknown counts. All registry PMID lookups resolved; their existence verifies the
publication identity, not the registry's biological source association.

## Reproducible artifacts

- `source_publications.csv`: every registry source, recorded identifiers, exact
  query URL, publication date, citations/year and identity/association gaps.
- `slot_comparisons.csv`: all 297 slots and installed numeric-output diagnostics.
  Every biological admission remains pending and **no preference is emitted**.
  Diagnostic scores use stored-value fractions, which are not assay completeness.
- `slot_searches.csv`: both search orders, hit counts, errors and truncation flags
  for every slot; 594 slot/order views reuse 432 unique queries.
- `literature_candidates.json`: abstracts and dated citation factors for retained
  hits. Comprehensiveness and selection scores remain unknown until admission.
- `source_review_queue.csv`: all 162 sources, inventory scopes/statuses, numeric
  slot links and explicit source/coverage/alternative review requirements. Numeric
  links are not a census of relationship or metabolite evidence; the inventory
  retains those units separately.
- `diagnostic_reviews.csv`: the six quantity/context traps and retention decisions.
- `host_progress.csv`, `host_followup.json`: measured host slot coverage and the
  publication-version/assay follow-ups described below.
- `requests/`: exact primary Europe PMC JSON responses, with URL, UTC retrieval
  time and failure status. `manifest.json` hashes the used responses, inputs and
  four core exports. `review_manifest.json` separately hashes review inputs/outputs.
- `audit.ipynb`, `review.ipynb`, `followup/review.ipynb`: executed analyses and
  exact-identity host follow-up. `verification.ipynb` records artifact reconciliation.

Run from the repository root with the existing Starplast Python environment:

```sh
python scripts/audit_dataset_selection.py --out results/dataset_selection_2026_10_07 --as-of 2026-10-07
python scripts/review_dataset_selection.py --out results/dataset_selection_2026_10_07
```

These commands reuse the frozen responses. A fresh citation snapshot must use a
new directory and its actual observation date; current counts cannot reconstruct
a historical count. Weights, the 30-day age floor, provider, input hashes and
Python/pandas versions are frozen in the manifest. Analysis used Python 3.12.13
and pandas 3.0.5 in the existing environment under a 4 GB memory cap.

## Search and scoring limits

Queries declare the species, assay/question and available context. Host baseline
queries do not require parasite mention; infection questions do. Tissue synonyms
and atlas names permit broader discovery but require sample-level review. Vector
queries name Anopheles rather than substituting its parasite. Multiword assay
terms are quoted. The recorded sorts are the REST API's `sort_cited:y` and
`sort_date:y`; SOAP-style sort text is not appended as free-text search.

Each query retains the first ten most cited and first ten newest hits. These
first pages **do not guarantee retrieval of the highest citations/year paper**:
middle-aged candidates can be missed. The 231 truncated slot views need deeper
queries/pagination and primary deposit review under 65.04. Reviews, irrelevant
species mentioned in abstracts and small targeted studies are still proposals,
never admitted datasets. No query hits means a search gap, not proof that no data
exist. Queries were refined after a stopped initial run; only the final used
response hashes are part of this audit.

Europe PMC counts are provider-specific, as its [search help](https://europepmc.org/help)
explains. Query syntax/order comes from its
[REST reference](https://dev.europepmc.org/RestfulWebService). Citation preference
is separate from accuracy/capacity metrics against ground truth. No runtime
default, installed measurement, slot, strategy or calibration was changed.

## Host status and concrete findings

| Host evidence | Present | Remaining |
|---|---|---|
| Human protein table | 7 source datasets; 20,989 reference rows | Independent human gene space and inference pack |
| Mouse protein table | 3 source datasets; 15,590 reference rows | Independent mouse gene space and inference pack |
| Toxoplasma host slot views | 8 / 26 | 18 unfilled |
| Plasmodium host slot views | 4 / 24 | 20 unfilled |

Thus **12/50 host slot views (24%) are populated**. This is storage/slot coverage,
not overall host-project progress. The host gene builders and inference packs
(64.21, 64.22, 64.28, 64.29) remain pending; protein accession rows must not be
silently treated as admitted host gene entities.

The installed [2026 RBC study](https://www.nature.com/articles/s41597-026-06792-5)
reports 5,264 proteins and fraction-specific PSM counts. The
[2017 quantitative RBC study](https://pubmed.ncbi.nlm.nih.gov/28689405/) reports
2,650 proteins, absolute copy numbers and validation using labeled standards.
The frozen Europe PMC snapshot gives **4.53 versus 29.12 citations/year**,
respectively. The older quantitative source deserves admission review as
complementary evidence; its copy numbers cannot silently replace fraction PSMs.
Both comprehensiveness and assay resolution matter, and these protein totals
are not a common mapped assay denominator.

A [September 2026 rhoptry-discharge journal paper](https://pubmed.ncbi.nlm.nih.gov/42791346/)
matches the installed preprint's authors, question and headline host screen.
Journal/preprint supplementary-table equivalence remains unverified, so the
original preprint provenance and measurements are retained. Citation counts are
not summed across versions, and a journal date does not reset first publication.
The GTEx v10 columns also need release-to-publication lineage review: a resolved
2020 atlas PMID alone cannot verify newer release-specific cell-type fields.

Five diagnostic flags compared different outputs of the same source (e.g. site
counts against complete boolean flags); the sixth mixed in vivo brain and in
vitro tachyzoite transcription while the incumbent's bibliometrics were unknown.
All six were rejected as automatic replacements. Future automatic preference
also stops when any admitted alternative has missing comparison metadata.

Next: resolve source/deposit lineage, expand truncated searches, admit biologically
comparable candidates, measure common assay coverage, and validate each proposed
addition/replacement with mapping/QC, leakage and affected benchmarks. Then update
layouts/calibration/tutorials for any actual data or strategy changes.
