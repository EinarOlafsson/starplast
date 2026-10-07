# 65.01 · Versioned dataset-selection preferences

Status: DONE — ✅, 2026-10-07.

Parent: [65 · Dataset selection audit](../open/65_dataset_selection_audit.md).

Delivered `SelectionPolicy`, `DatasetCandidate`, factor-exposing `evaluate`,
scope/provider/snapshot-consistent `rank_candidates` and conservative `preferred`.
Admission of organism/context/assay/quantity/quality precedes citation preference.
Annual citations use actual first-publication age and a 30-day floor. Recency
bonus is at most 10%, with five-year decay; comprehensiveness bonus is at most 25%.
Unknown metadata does not become zero. Actual zero-citation datasets remain
eligible, with coverage/recency tie breaks. No scores are biological accuracy.

Tests: **389 passed**, policy + all module docstrings + organism invariants under
a 4 GB cap, Python 3.12/offscreen Qt. Tests verify normalized rather than raw
citation preference, bounded bonuses, zero/unknown distinctions, admission gates,
scope/provider/date consistency and malformed candidates/policies.

Rule persistence: AGENTS.md and NEXT_SESSION.md; documentation:
`docs/dataset_selection.md`, linked from `docs/API.md`.

Commit: **Establish reproducible citation-rate preferences for dataset selection**,
published to nightly. Actual publication associations, assay completeness and
candidate eligibility remain audit inputs. This policy-only milestone does not
change source tables, existing default slot outputs or strategies; replacing an
incumbent requires 65.04 validation and the affected calibration checks.
