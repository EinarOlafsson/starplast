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
The control has its own immutable artifact manifest. 908 relevant checks passed,
one optional pdoc documentation module skipped locally. Item 64.10 is still open:
other missing adapters, the logistic variant, full outer coverage and biological
source admission are pending.

Held-out artifact: `a6062f243a4e4fdf4702e9fc9162c438139da54684bf2873948ee1aa04bfa4ad`.
Base model: `de4e8e5d2bfdf5308e5a4ea9952aeb80dcc100120a5949926025d6c7e122f831`.
Calibration: `cc5b0da9c6a5e767376e64cb331aa893bb7bc2d07e52024a5eb018280b7cdc38`.

This final revision also rechecks canonical base-model payload identity after loading, refusing in-memory mutations of its class vocabulary/state. The earlier prototype and original code remain in the sibling directory; numerical outcomes are unchanged.
