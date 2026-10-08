# FN-EC-01: frozen functional annotation-profile recovery

`pilot.ipynb` actually fitted and replayed native kNN against complete EC major-class
profiles. Source review: `../functional_source_review_2026_10_08_v2/`. Typed held-out
artifact: `349317afd8f24b860d94ae1e11280a6130bf451700b02b9d2145fab85fc99162`.

The fixed seed 23 split protects recorded homology groups: 679 train, 175 tuning,
190 calibration and 182 test genes. Only training genes fit ranks, feature selection
and imputation. Tuning/calibration roles were not used. Fixed k=15/minimum share=0.3;
349 eligible numeric inputs exclude source EC, domain, homology and publication/full-text
attention counts. Graph input is empty. Overlapping enzyme-class memberships stay
in the complete profile, including two multi-class test genes. One unseen test
profile stays in the evaluation denominator.

| Measure | Result |
|---|---:|
| Exact profile recovery, all 182 test genes | 71/182 = 39.0% |
| Answered genes | 171/182 = 94.0% |
| Exact recovery among answers | 71/171 = 41.5% |
| Abstentions | 11 |
| Train-only majority recovery | 66/182 = 36.3% |

Native predictions, all class scores and floats, model state, complete-profile and
seven overlapping member-class scorecards, and matched controls replay **exactly**.
Source/code hashes and frozen roles are retained. Biological precision/recall and
calibrated confidence remain unavailable: this recovers recorded annotations,
not independently measured activity. Unannotated/unresolved genes are unknown,
not verified negatives. No new unknown-gene deployment claims are admitted.
The result is weak and does not support a broad accuracy claim.

Original nodes, claims, recipes and legacy track records are unchanged. The first
missing-context diagnostic and v2/v3 prototypes are retained separately. The v2 to
v3 revision hardens malformed-source handling and removes attention counts;
v4 preserves original descriptions/per-class provenance in the browser and updates
the declared code identity without changing v3 predictions or metrics. These are
correctness revisions, not held-out setting search. The application view is
`../functional_ec_ui_bundle_2026_10_08_v4/`.
