# Initial ridge pilot: preserved metric diagnostic

Native training ranks, predictions, fitted state and matching baselines replay
exactly for 1,101 held-out genes. However, constant nonbinary-float baselines have
small nonzero computed standard deviation because of rounding. The preceding
metric guard consequently emitted arbitrary top/bottom-decile recall for a
constant baseline (correlations returned null with scipy warnings). These cards
are preserved with original code/input identity and are not canonical. No model,
source value or prediction changed. The corrected snapshot uses exact nonzero
range to detect variation; see `ridge_value_pilot_2026_10_08_v2/` after verification.
Do not interpret the baseline decile fields in this initial packet as capacity.
