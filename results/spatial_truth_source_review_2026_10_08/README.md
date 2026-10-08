# Spatial truth-source review · 2026-10-08

Bounded partition **GT-SPATIAL-01**, under 65.04/64.10; two existing registered
papers, no new organism or target promotion. Primary articles:
[Barylyuk 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7670262/), PMID 33053376;
[Chisholm 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13369866/), PMID 42218142.
Original articles, tables and figures remain in the external dataset archive.

`preflight.ipynb` verifies existing bindings, primary PMC versions, identities,
published filenames/MD5s and schemas. Three additional companions were retrieved:
Toxoplasma Table S3; Plasmodium Supplementary Data 2 and supplementary figures.
`table_audit_v2/audit.ipynb` retains complete marker/prediction fields and every
source row, with exact installed-ID matches. It recovers four further primary
inputs: Toxoplasma Figures 1/2 and Table S9, Plasmodium Figure 3. All seven
companions have public primary metadata/MD5 and local SHA receipts. No restricted
access is bypassed. Unavailable content is a gap; independent work continues.

| Source field | Source rows | Marker rows | Classes | Exact-ID marker rows | Unmapped markers |
|---|---:|---:|---:|---:|---:|
| Toxoplasma final markers | 3,832 | 718 | 26 | 716 | 2 |
| Plasmodium S1/S2 markers | 3,000 | 542 | 20 | 542 | 0 |
| Plasmodium S1/S2/S3 markers | 3,000 | 484 | 20 | 484 | 0 |

The two unmapped Toxoplasma markers are retained explicitly; no aliases are
invented. All 27 unmapped Toxoplasma source rows remain present. Neither table
has duplicate or missing accessions. These are protein-table counts, not an
assertion that every grouped protein measurement unambiguously measures one gene.

The Toxoplasma paper says its 718 final markers combine 656 prior markers and
all 62 newly measured microscopy outcomes. Those 62 cannot independently test
the published final classifier trained on them. Both papers' marker selection
mixes previous knowledge, inferred function/location and profile consistency;
marker status alone is not independent experimental truth.

Plasmodium has 140/119 literal marker versus final-prediction disagreements in
the two fields. Coarse versus refined classes are unresolved, so these are not
scored as biological errors. **Final location (S1-S2) includes rhoptry assignments
from S1/S2/S3**; the source field name must not imply exclusively S1/S2 provenance.
The raw source strings and separate marker contexts are preserved.

`microscopy_review/review.ipynb` records a manual visual transcription of the
nine primary Figure 3 panels. Shorthand identifiers match unique full IDs in
the same paper's verified primer table; full prefixes are never guessed. All
nine map exactly to installed genes, appear in the spatial table and are absent
from **both** training-marker fields. The complete 15-target attempted cohort
is retained; the other six IFA outcomes remain unknown, never negative truth.
Panel compartment and co-stain labels remain exact, with no automatic translation
to finer classifier niches. Single tagged-line experiment, profiling-based target
selection, assay context, source dependence and compartment taxonomy remain gaps.

**Zero biological benchmarks admitted, no replacement or installed numerical
change, no new performance/calibration claim.** The nine measured reporter
outcomes are promising assay candidates; population-wide accuracy is not
identified by this selected validation cohort. Host review and remaining
categorical adapters continue independently.

`verification.ipynb` reproduces every saved marker/prediction row, source order,
class count and exact-ID mapping; verifies four extra media hashes and the
9-positive/6-unknown microscopy distinction. **426 relevant checks passed**
(source/PMC/table recovery, ground truth, provenance, organisms and docstrings).
The first table audit's attachment discovery failed because the publisher used
`media` rather than `supplementary-material`. That diagnostic and original code
are preserved in `table_audit/`; v2 continues other figures and discovers the
exact attachment across primary XML. All numerical marker rows match v1 exactly.
Current acquisition/audit code and reusable download/notebook helpers are in
`code/`; the preflight and phase manifests retain original input/output hashes.
