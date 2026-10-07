# Measurement trace census and exact raw-to-host reproduction

Executed October 7, 2026 using Python 3.12.13 and a 4 GB systemd scope.
`provenance.ipynb` contains the executed complete inventory trace census and mouse
BMDM baseline pilot. `manifest.json` pins all inputs, located file identities and
three JSON outputs. `traces.json` uses schema 2 and validates content identities.

All 162 registered sources reconcile to 628 source/output addresses. 364 addresses
have a located processed input, mapping reference or installed-cache identity;
627 lack an audited measurement unit and mapping audit. Registry exploration
scope and actual assay species are distinct fields. Unknown species, licenses,
mapping losses, origin-paper associations and transformation parameters are
explicit gaps. A cache/file hash does not establish an experimental association.

The mouse BMDM pilot reproduces every installed nonmissing TPM: **15,437 exact
matches, no lost installed IDs, maximum difference zero**. Its source has 35,275
unique Ensembl genes: 15,531 map, 53 are ambiguous and 19,691 unmapped. The mapped
genes reach 15,437 reviewed UniProt proteins, including 58 many-to-one targets.
Three M0 replicates are converted from FPKM to TPM before averaging, mapping and
selection of the most expressed source row for each target. The original pipeline
serialized the result to TSV with `%.10g`; reproducing that step explains the
initial maximum pre-serialization difference of 0.0000037736972444690764.

The failure diagnosis is retained in
`../measurement_provenance_diagnostics_2026_10_07/diagnostic.ipynb`. Preliminary
schema-1 and schema-2 censuses are retained in the dated sibling directories with
exact original code copies. This snapshot supersedes both and includes the 12
newly recovered PMC supplement bindings. It changes no installed measurements,
strategy scores, source defaults or calibration artifacts.

Validation: 509 focused provenance/recovery/query/inventory/dataset/evidence
checks passed. Another 380 checks for current PMC metadata, provenance,
recovery, documentation, docstrings and organism invariants passed; one optional
local pdoc module skipped. Python 3.10 syntax and frozen input/output checks
were verified. Biological admission and missing raw inputs remain instructions
65.04 and 66; this is the completed provenance contract and gap census, 64.03.
