# Frozen native kNN conformal label sets

Fixed first pilot on the registered Toxoplasma compartment prediction candidate,
seed 17, alpha 0.1 and per-class thresholds. The previously frozen training rank
state is reused without selection from outer-test outcomes. Train/calibration/test
populations are 2,126/569/560. The separate fitted base, calibration digest and
held-out set artifact preserve their roles and dependencies. Installed data,
runtime strategies and calibration are unchanged.

| Measure | Conformal pilot | All-training-classes control |
|---|---:|---:|
| Retained test genes | 560 | 560 |
| Native model classes | 26 | 26 |
| Empirical set coverage | 0.973214 | 1.000000 |
| Mean set size | 24.455357 | 26 |
| Singleton calls | 0 | 0 |
| Native set efficiency | 0.061786 | 0 |

This is agreement with stored classifier output, not independent biological
accuracy. Exchangeability is unresolved. Eleven rare classes use the overall
fallback; class-specific guarantees are not demonstrated. Zero singleton calls
means all 560 genes abstain from categorical calling; called-only accuracy is
unavailable. Raw base scores remain available and are not calibrated gene-level
correctness probabilities. Broad coverage alone would conceal the weak capacity
to narrow the answer.

`pilot.ipynb` contains actual execution. `verification.ipynb` compares every
native call, set, score and threshold on identical frozen partitions, verifies
all record-derived cards/baselines and freezes the broad-set control. Only the
native partition chooser is temporarily replaced for diagnostic parity and is
restored; no released inference code is changed. Rare classes remain eligible,
without the released deployment runner's prefilter. Typed set members retain
their meaning independently of the native display separator. Unbounded
thresholds use explicit status/null encoding.

Input/output hashes are in `manifest.json`; original code is under `code/`.
The control has its own immutable artifact manifest. 906 relevant checks passed,
one optional pdoc documentation module skipped locally. Item 64.10 is still open:
other missing adapters, the logistic variant, full outer coverage and biological
source admission are pending.

Held-out artifact: `ce629ecf4177dc01c6af04e789dfdd53c1aa206caaa45c5db4420f226f1db917`.
Base model: `806bb3df0c83834eb8960e99b95b8a5df4acc1b5341f559b7b8461d219ef0683`.
Calibration: `ea882dbc3a32e17b7d53430e1e4e4f4f8e90ca4bac7b2d0d7ca987457565af2b`.

This is the initial prototype, superseded by `../conformal_label_pilot_2026_10_08_v2/` for post-load payload identity enforcement. Its executed verification and exact original code are retained.
