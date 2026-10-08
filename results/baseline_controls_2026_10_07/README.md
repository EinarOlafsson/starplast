# Shared baseline and control fixtures

Executed `controls.ipynb`, synthetic fixtures only. All six task recipes retain
their denominator and interpretation limits. Training-only majority/prevalence
and mean/median use the identical frozen evaluation cohort; predictions retain
split and universe identities. Uniform random ranking pins the candidate universe
and query-seed exclusions.

Planted positives recover completely; reversed and null rankings are refused by
the synthetic harness. Positive-only retrieval leaves precision/AUROC/AUPRC
unavailable rather than treating unlabelled entities as negatives. The network
control performs 20 swaps preserving degrees and detection-stratum contacts;
mixing remains unestablished. A complete graph permits no swaps and is explicitly
an unavailable control.

No runtime strategy change or independent biological accuracy/capacity claim.
These shared definitions and helpers need task-specific adapter execution against
the same admitted truth cohort and split as each strategy. Grouped nulls,
map/replication adapters and verified biological negatives remain follow-up work.
