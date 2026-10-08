# Offline functional recovery scorecards

`packaging.ipynb` executes the checksum-preserving compact view of the exactly
replayed [v3 pilot](../functional_ec_knn_pilot_2026_10_08_v3/summary.json).
Its original typed artifact identity is
`b7de54fd3a835738309a9262d543e6d49921b8d7c26cf6dd29e84736d2987034`.

All source artifact payload bytes are verified. Every displayed row, label card,
complete-profile class card, seven major enzyme-class cards and matched training
control retains its exact original source payload hash. The original model state
stays upstream, referenced by its path/hash/size. The compact application bundle
has an external SHA-256 pin:
`807c4ea23239930722742348ddaadb7d51bcf961ff2736de53e4cc28ba62dc16`.
Packaging source and a manifest are retained here; the shipped copy is
`starplast/data/functional_results.json`.

Discoveries links its source EC label to **Functional tests**, with label, class,
profile, strategy and training-control views. All 182 held-out native outcomes
remain available, including 11 abstentions; gene clicks use the existing evidence
panel. Class outcome rows retain missed calls and false calls from other profiles.
Control cards do not display native predictions as if they were control outcomes.
Missing or changed bundles show an explicit unavailable reason.

This is frozen **reference-annotation recovery**, not independent biological
activity accuracy. Biological class precision/recall and calibrated confidence
remain unavailable. Missing or unresolved annotations never become verified
negative labels, and no new verified unknown-gene claims are created. The earlier
v2 packaging prototype and relative-path preflight error are preserved separately.
