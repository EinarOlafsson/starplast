# Superseded ground-truth diagnostic

Use `../ground_truth_registry_2026_10_07_v2/` for the reviewed census.
This initial snapshot used assay kind to classify unreviewed CRISPR outputs and
missed the liver-stage orthology transfer. No entry was admitted and no runtime
data was changed. Review replaced that shortcut with explicit source/output
declarations, including derived in-vivo composites and predicted AlphaFold
confidence. Original script/module bytes are preserved in `code/` with hashes
matching this snapshot's manifest. This diagnostic must not support biological
accuracy claims.
