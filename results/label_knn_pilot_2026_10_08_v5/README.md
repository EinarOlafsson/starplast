# Frozen label kNN pilot — label_knn_pilot_2026_10_08_v5

Current first native feature-kNN pilot. Training average-rank percentiles exactly match Context.matrix on training rows. Held-out features use frozen training distributions, missing values become zero, and all-missing training columns are withheld. This is an explicit inductive protocol, not the original whole-context feature ranking.

The source target is a stored localization prediction, not independent biological ground truth. No biological benchmark is admitted; no installed data, runtime strategy or calibration is changed. Complete label-adapter coverage remains open under 64.10.

Split: 3,827 eligible genes; train/tune/calibration/test 2,126/572/569/560. All 560 test entities and 26 native score columns are retained. Calls: 432; abstentions: 128. Surrogate all-hidden accuracy 0.3500; called-only accuracy 0.4537. Training-majority baseline 0.2304.

Held-out artifact identity: `4924dfa57c5fc39a6bd3f92718e031d60d23163b34da869ef102a5c1ec03e350`. Inputs and JSON outputs have SHA-256 receipts in `manifest.json`. Original matching Python inputs are preserved under `code/`; a changed live source must not be substituted silently for those identities.

`pilot.ipynb` records successful execution where available; verification/diagnostic notebooks document the associated checks. Earlier run directories remain immutable diagnostics; their old code identities are retained.
