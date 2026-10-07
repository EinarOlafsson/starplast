# Literature verifier audit — 2026-10-07

The literature graphs do not provide the missing broad independent coverage. No Toxoplasma
candidate passes both the existing dependence and calibration limits. One Plasmodium candidate
passes the exploratory screen and reaches two additional in-range claims: it agrees with one
and disagrees with the other. It remains a candidate, with no published claims promoted.

The audit evaluated abstract and full-text propagation separately, using positive
attention-corrected residuals and keeping raw counts as controls. It used the shipped generator
for each of the 12 labels and the same five orthogroup folds, seed and labelled genes as the
track record. Identity, truth, seed and fold mismatches abort the comparison. Banned layers,
missing corrections and nonpositive evidence are refused without substituting another graph.

There are 48 combinations: 42 measured and six unavailable (Plasmodium full text). The saved
tables contain 92,354 held-out predictions and 234,730 calls or abstentions on unknown genes.

| Corrected candidate | Shared-mistake upper 95% bound | Additional in-range claims reached | Decision |
|---|---:|---:|---|
| Tg compartment, abstracts | 2.721 | 15 | Fails independence (<1.6 required) |
| Tg compartment, full text | 2.044 | 205 | Fails independence |
| Tg LOPIT unified, abstracts / full text | 2.264 / 2.412 | 15 / 205 | Both fail independence |
| Tg cell-cycle phase, abstracts / full text | 2.639 / 1.901 | 95 / 950 | Both fail independence |
| Tg phenotype screens, full text | Five pass the numerical dependence gate | 0 | Combined calibration error 0.064–0.100 exceeds 0.05 |
| Pf localization, abstracts | 3.350 | 6 | Fails independence |
| Pf transferred phenotype, abstracts | 1.331 | 2 | Exploratory candidate; requires further validation |
| Pf export, abstracts | 1.000 | 0 | No additional unknown-gene coverage |

For the phenotype candidate the combined calibration error is 0.0243, with 600 confident
held-out calls at 87% precision. Those are pooled recipe results, not a measured precision for
the two newly reached genes. Measure performance on that stratum and possible reuse of source
papers before promotion. The two extra verdicts are in [newly_reached.csv](newly_reached.csv).
Outside-range genes remain separate throughout; reaching them does not calibrate their certainty.

Read [audit.ipynb](audit.ipynb) for the executed measurements and [review.ipynb](review.ipynb) for
the decision and exact additional gene calls. [summary.csv](summary.csv) contains every candidate,
control, refusal, interval, agreement count and coverage count. [manifest.json](manifest.json)
records source/data SHA-256 hashes, numerical settings and package versions; the hashes matched
before and after the final run. `heldout.parquet` and `unknown_calls.parquet` preserve the rows
behind the report. The standalone notebook can be re-executed from within this checkout.

Reproduce in a fresh directory, under a memory cap and with a RAM lease:

```bash
systemd-run --user --scope -p MemoryMax=8G env CUDA_VISIBLE_DEVICES= \
  OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 \
  ~/anaconda3/envs/starplast/bin/python scripts/audit_claim_verifiers.py \
  --out results/claim_verifiers_repeat
```

Validation: **415 passed**, covering the new audit (23 checks), claims, track record, organism
registry and docstrings. Tests cover corrected-vs-raw reach, evidence refusals, stale record
alignment, outside-range handling, minimum calibration/precision/count/reach requirements, CLI
completion and replay from a clean notebook namespace. The first complete measurement attempt
hit a final console-reporting `NameError`; it was fixed, covered and followed by a full successful
rerun. No runtime strategies, calibration tables or shipped claims changed.

Next: evaluate orthology or separate experimental evidence for broader coverage. Frozen
prospective testing remains open under instruction 63.
