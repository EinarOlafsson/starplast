# 64.03 · Trace each measurement to its source

Status: DONE — ✅, 2026-10-07. Provenance contract, complete gap census and exact host pilot verified.

Parent: [64 · Information space and inference atlas](../open/64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.02 |
| Estimated remaining engineering time | 0 h |
| Existing code to inspect | `starplast/datasets.py`, `starplast/deposits.py`, `starplast/identity.py`, `starplast/results.py` |

## Deliverable

A source link for every inventory measurement: paper/deposit, file hash, source species, unit, transformation, mapping and redistribution status.

## Controlled scope

Implement the provenance contract and audit the existing tables; new source acquisition is separately scoped.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] Unresolved source records stay visibly unresolved; citation identifiers come from verified records.
- [x] Gene↔protein and strain/assembly mappings record cardinality, ambiguity and mapping loss.
- [x] A selected measurement traces through its transforms to the source file; inferred/transferred labels remain distinct from directly measured evidence.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

`starplast/provenance.py` implements schema 2: content-addressed files and transforms,
registry exploration scope separated from actual assay species, measurement units,
publication/file association status, mapping ambiguity/cardinality/loss, license,
redistribution and visible gaps. Cycles and duplicate output addresses are refused.

The executed census in `results/measurement_provenance_recovered_2026_10_07/`
accounts for all 162 sources and 628 source/output addresses. The selected raw-to-installed
pilot is mouse BMDM baseline TPM: 35,275 unique source genes; 15,531 map unambiguously,
53 are ambiguous, 19,691 unmapped; 15,437 target proteins with 58 many-to-one targets.
All 15,437 installed nonmissing values are reproduced exactly, including the original
`%.10g` TSV serialization. The initial rounding failure and preliminary snapshots
are preserved with original code copies and supersession notes.

Source recovery is a bounded follow-up under instruction 66. The initial source census
retains all 162 locations and stale links. GTEx v10 is located by an exact URL/name/species
receipt and SHA256. The current PMC Cloud Service recovered 12 named missing inputs
by matching PMCID/PMID, a unique published version, exact filename and metadata MD5;
other failures and version/name mismatches stay visible. No runtime data is changed.

Checks: 509 focused inventory/query/dataset/evidence/provenance checks passed; a separate
380-check provenance/recovery/docs/organism run passed with one optional pdoc module
skipped locally. The executed raw-to-installed census passes; source/output hashes,
Python 3.10 syntax and whitespace checks are recorded with the artifacts.

The commit **Trace measurement provenance and recover missing published inputs** is
published on `origin/nightly`. The parent tracker and handoff move on to 64.04.

Limitations: 627 addresses still lack an audited measurement unit. Most original assay
species, origin-publication associations, historical parameters, strain/assembly mapping
and licenses require source-specific review. An installed cache is explicitly not a raw
input. A resolved publication ID is not verification that it generated the file. Unknown
mapping losses and missing values are not negatives. Completing this contract/gap census
does not complete biological source admission, host gene spaces or inference scorecards.
