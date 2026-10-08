# Frozen label kNN pilot — label_knn_pilot_2026_10_08_v4

Fully replayed z-score diagnostic, with full floating payloads and one-ULP-sensitive metric hashes. Feature pipeline review found native rank scaling was needed; this run is not the canonical native feature-transform pilot.

The source target is a stored localization prediction, not independent biological ground truth. No biological benchmark is admitted; no installed data, runtime strategy or calibration is changed. Complete label-adapter coverage remains open under 64.10.

Split: 3,827 eligible genes; train/tune/calibration/test 2,126/572/569/560. All 560 test entities and 26 native score columns are retained. Calls: 431; abstentions: 129. Surrogate all-hidden accuracy 0.3143; called-only accuracy 0.4084. Training-majority baseline 0.2304.

Held-out artifact identity: `bdf9c9eb83cac7463318e13e4ac6c9fb0e92de9ca3a644150389ddc41ec3775a`. Inputs and JSON outputs have SHA-256 receipts in `manifest.json`. Original matching Python inputs are preserved under `code/`; a changed live source must not be substituted silently for those identities.

`pilot.ipynb` records successful execution where available; verification/diagnostic notebooks document the associated checks. Earlier run directories remain immutable diagnostics; their old code identities are retained.
