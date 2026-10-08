# Shared scorecard task views — software acceptance

**539 checks passed** across all six task contracts, ten synthetic aggregate cards, prediction-set and numeric-interval metadata, explicit controls, positive-only truth, unavailable states, evidence quality and individual outcomes.

The executed `audit.ipynb` computes every recorded input/card/check in its shared namespace. `exports.json` retains the exact native cards; both independently constructed pure consumers agree on all float/null values, HTML and JSON. Metric routes resolve their task-specific definitions, including distinct nominal set and interval coverage. Code/card/record/input identities and runtime versions are retained in `summary.json` and byte receipts in `SHA256SUMS`.

These are authored software fixtures, with no source acquisition, estimator fitting, biological benchmark admission or calibrated confidence. No Qt application ran. Actual desktop host checks and scientific source/split/calibration/task-admission gates remain separate. **64.17 remains open.**

Run from the repository with the existing environment and a 400 MB external memory cap:

```bash
systemd-run --user --scope -p MemoryMax=400M /home/carruthers/anaconda3/envs/starplast/bin/python scripts/audit_scorecard_task_views.py
```

The script refuses an existing output directory. Use `--out` with a new destination to reproduce; notebook cell timings may differ, while canonical cards and displayed values must agree.
