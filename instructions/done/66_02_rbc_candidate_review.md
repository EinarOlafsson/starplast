# 66.02 · Review the quantitative RBC candidate

Status: DONE — ✅, 2026-10-07 (local date). Controlled source review complete;
quantitative pack admission remains under 65.04.

Primary PRIDE project metadata verifies human samples, PMID 28689405 and the
2017 paper DOI. The 28 MB processed MaxQuant archive matches the deposit's SHA1
and byte count. Two DOI-matched ACS Figshare supplements match publisher MD5/size.
URLs, hashes, licenses, failures and executed retrieval notebooks are retained.

The published Table S3 has **2,653 protein groups and four donors**, separately
for whole erythrocytes and white ghosts. Eight copy-count fields have no missing,
negative or nonfinite values. Published zeros remain numerical observations;
they are not certified absent proteins or negative biological labels.
Whole-cell donor rank correlations are 0.842–0.911 and ghost correlations
0.884–0.945. The source-group/reference audit records 2,195 unambiguous groups,
330 ambiguous groups and 128 unmapped groups; 59 target proteins receive multiple
groups. With ambiguity and reverse collisions withheld, **2,061 single-protein
projections** remain. Original groups and donor values are preserved externally.

Decision: **retain installed fraction PSMs and intact-cell surface profiling**.
Whole-cell/white-ghost abundance does not replace selective surface measurements,
and copy counts are not interchangeable with fraction PSMs. The 2017 quantitative
source remains a separate candidate with declared quantities and group evidence.
Its October 7 Europe PMC snapshot has 268 citations, 29.115 citations/year and
recency factor 1.01587. Common assay comprehensiveness remains unknown, so the
admission-gated policy does not fabricate a preference score. CC BY-NC publisher
constraints remain explicit; original numeric source/derivative tables stay
external and are not bundled in the wheel. Rights-compatible pack admission,+projection validation and benchmark evidence remain 65.04.

Artifacts: `results/host_candidate_deposits_2026_10_07/`,
`results/host_candidate_downloads_2026_10_07/` and the reviewed
`results/rbc_candidate_review_2026_10_07_v2/`. The preliminary comparison missed
the actual four surface columns; it is preserved with original code and superseded
by the registry-driven comparison of all six incumbent measurement columns.

Checks: **548** source/host/deposit/review checks and **625** registry/selection/
generated-catalogue/review/docstring/organism checks passed. Exact gene-version
comparison, ambiguous group preservation and swapped-fraction refusal have tests.
Executed artifacts and hashes were checked. No installed data, slot, graph,
strategy output or calibration changed. The coherent host-review commit is
pushed to `origin/nightly`; the parent tracker and handoff retain remaining work.
