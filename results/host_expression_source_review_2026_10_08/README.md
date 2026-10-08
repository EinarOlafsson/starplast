# Host expression source/version review · 2026-10-08

Bounded **GT-HOST-TX-01**, under 65.04: existing GTEx and FANTOM5/Expression
Atlas E-MTAB-3579 host sources. No installed numerical change or biological
benchmark admission. `review.ipynb` verifies archived input receipts, primary
deposit/resource/publication metadata and original field schemas.

The [official GTEx catalogue](https://www.gtexportal.org/home/downloads/adult-gtex#qtl)
lists v11 using GENCODE 47 with no new samples/donors relative to v10. Its exact
listed median-summary filename was discovered in the public primary GCS object
listing and downloaded with provider MD5/size validation. Original filename and
SHA receipt are in `sources.json`; the file remains external. A static JavaScript
display gap did not stop the object review or the independent Atlas task.

| GTEx version comparison | Count |
|---|---:|
| v10 gene rows | 59,033 |
| v11 gene rows | 74,628 |
| Shared stable gene identifiers, preserving PAR_Y | 56,474 |
| v10-only / v11-only identifiers | 2,559 / 18,154 |
| Current native v10 / v11 protein projections | 19,088 / 19,099 |
| Changed fibroblast / hepatocyte values on shared genes | 26,865 / 26,158 |

Both files contain the selected cultured-fibroblast and hepatocyte fields.
Shared-gene rank correlations are 0.980913/0.980147. More annotated gene rows
are not measured protein comprehensiveness or proof of improved accuracy.
`selection_review.json` records the user policy, current Europe PMC family
citation snapshot and pending admission: neither version receives a preference
score/promotion. The citation's indexed date does not establish release-specific
publication lineage; v11 must not receive a new citation-age clock.

The provider-checksummed primary LCM README documents 172 sequenced pilot
samples from five tissue types, including targeted liver hepatocytes. This
supports an enriched reference, with fixed-tissue/LCM context distinct from
infection-study donors/cultures. Original README remains external; exact
association is in `lcm_readme_receipt.json`/`lcm_metadata.ipynb`. No per-sample
or instrument files were downloaded.

Atlas `tpmss.tsv` was downloaded from the freshly verified accession resource
index, retaining its original filename. Every parsed value/column/order matches
the historical local file (18,835 genes × 170 tissue/stage columns), despite
different generated-header hashes. Original design/deposit receipts remain
verified; the historical local name was not renamed to imply a primary filename.
CAGE is separate from RNA-seq TPM; four adult brain regions and juvenile biceps
femoris remain distinct contexts.

`comparison_v3/comparison.ipynb` records complete selected gene rows and native
source replay through explicit temporary symlinks. Replay is **not exact**:

| Installed field | Installed values | Missing replay | Extra replay on installed universe | Numerical differences on overlap |
|---|---:|---:|---:|---:|
| GTEx fibroblast TPM | 19,087 | 3 | 2 | 0 |
| GTEx hepatocyte TPM | 19,087 | 3 | 2 | 0 |
| Atlas brain TPM | 9,763 | 8 | 1,342 | 1,741 |
| Atlas muscle TPM | 8,635 | 7 | 1,276 | 0 |

Every installed protein and both values are retained in replay rows; unmatched
keys are explicit. No mismatch is hidden by deleting observations. Current
and September 25 archived mapping hashes are identical; those two mapping
copies alone do not explain the legacy differences.

`cutoff_diagnostic.ipynb` tests the already documented **0.5 cutoff**, without
threshold search. The recovered primary export has a zero query cutoff and
selected-region minima below 0.5. Applying 0.5 **per region before averaging**
eliminates all 1,741 overlapping brain-value differences: 9,755 brain and 8,628
muscle values match exactly. Eight/seven missing and 22/18 extra projected values
remain. This explains fidelity; transform/mapping revision requires a controlled
decision and affected checks. Missing/non-detected expression is not a biological
negative or zero. Original gene values remain intact.

Diagnostics retained: `comparison/` first refused duplicate keys produced by
splitting at the first dot (PAR_Y collapsed); `comparison_v2/` correctly refused
non-exact installed replay. Original code/partial outputs remain. v3 preserves
PAR_Y and original versioned identifiers, records all discrepancies and continues
the other source. `verification.ipynb` replays every reported count, gene/value
retention, primary table equivalence, provider hash and policy refusal.
**534 relevant checks passed** (host, source recovery, selection audit/policy,
provenance, ground truth, organisms and docstrings).

Next bounded dependency: **GT-HOST-TX-02**, with 66.04/64.21/64.22: review
installed mapping lineage and PAR_Y; preserve per-gene evidence and the explicit
legacy cutoff before corrected protein projection or v11 promotion. Host gene
spaces/inference packs remain open. Continue independent tasks if publication
or content remains unavailable.
