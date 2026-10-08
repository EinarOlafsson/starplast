# Frozen label kNN pilot — label_knn_pilot_2026_10_08_v2

Z-score diagnostic with tuple/list JSON normalization. Payload rows and scores use pandas default 10-digit JSON rounding; the later full-precision check supersedes its serialization convention.

The source target is a stored localization prediction, not independent biological ground truth. No biological benchmark is admitted; no installed data, runtime strategy or calibration is changed. Complete label-adapter coverage remains open under 64.10.

Split: 3,827 eligible genes; train/tune/calibration/test 2,126/572/569/560. All 560 test entities and 26 native score columns are retained. Calls: 431; abstentions: 129. Surrogate all-hidden accuracy 0.3143; called-only accuracy 0.4084. Training-majority baseline 0.2304.

Held-out artifact identity: `62cd18402d06e3b46890115d2badb19288d73e56968e953c6e9ad120505e14ae`. Inputs and JSON outputs have SHA-256 receipts in `manifest.json`. Original matching Python inputs are preserved under `code/`; a changed live source must not be substituted silently for those identities.

`pilot.ipynb` records successful execution where available; verification/diagnostic notebooks document the associated checks. Earlier run directories remain immutable diagnostics; their old code identities are retained.
