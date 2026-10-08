# Frozen categorical ortholog transfer

The existing Plasmodium spatial-prediction column is projected by native
orthogroup aggregation onto the frozen Toxoplasma compartment cohort. Receiver
training genes alone fit the native source-to-label mapping. The projection
context contains no receiver target column. No installed data, runtime strategy
or calibration is changed; this is one bounded adapter under open item 64.10.

Both source and receiver targets are prediction-grade. Stage/context semantic
equivalence and biological independence remain unresolved. This pilot measures
agreement with stored predictions, not independent biological accuracy.

| Quantity | Result |
|---|---:|
| Eligible receiver genes | 3,827 |
| Train/tune/calibration/test genes | 2,126 / 572 / 569 / 560 |
| Receiver genes with mapped donor values | 967 |
| Mapped training/test genes | 550 / 126 |
| Retained test rows | 560 |
| Calls / abstentions | 123 / 437 |
| Correct / wrong calls | 65 / 58 |
| All-hidden / called-only agreement | 0.116071 / 0.528455 |
| Same-cohort training-majority baseline | 0.230357 |

The native method provides no class-score matrix, support or calibrated
confidence; AUROC/AUPRC remain unavailable. Missing donor mappings and unsupported
categories remain explicit abstentions. Three mapped test values have no learned
source category. Per-class precision, recall/F1, confusion and group bootstrap
are derived from all retained rows.

`pilot.ipynb` executes acquisition of the stored tables, native projection,
fitting and freezing. `verification.ipynb` checks every input/output hash, exact
native aggregation/calls, record-derived metrics and both matched baselines.
The initial artifact preflight refused required-field omissions before creating
an artifact; corrected scope and typed categorical output are recorded. Original
matching code is under `code/`; immutable JSON receipts are in `manifest.json`.

Artifact identity: `65034bf84b7dc87aefee2e21622294cb33f0789b4e9dc41fb81243ee9ec33f07`.
Mapping identity: `e1da9a0c47cecdd47087c34026ddc5e5add6d4f950570a81574c481b40654d80`.
464 relevant automated checks passed. Missing adapters, full outer coverage and
independent biological truth remain open. This item receives no completion tick.
