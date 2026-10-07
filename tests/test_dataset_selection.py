"""Selection preferences must not choose the wrong biological question or invent metadata."""
from dataclasses import replace

import pytest

from starplast.dataset_selection import DatasetCandidate, SelectionPolicy, evaluate, preferred, rank_candidates


AS_OF = "2026-10-07"


def _candidate(key="synthetic", **kwargs):
    values = dict(dataset_id=key, scope="synthetic_same_assay_context", publication_id="synthetic-paper",
                  first_publication_date="2025-10-07", citations=20, comprehensiveness=.5,
                  citation_provider="synthetic-provider", citation_snapshot=AS_OF,
                  admission="accepted", admission_reason="Synthetic known-good assay fixture")
    return DatasetCandidate(**dict(values, **kwargs))


def test_annual_citation_rate_beats_total_citation_count():
    old = _candidate("older", first_publication_date="2006-10-07", citations=200)
    new = _candidate("newer", citations=40)
    assert old.citations > new.citations
    assert preferred([old, new], AS_OF) == "newer"


def test_recency_and_comprehensiveness_are_preferences_with_declared_bounds():
    baseline = _candidate("baseline", comprehensiveness=0)
    comprehensive = replace(baseline, dataset_id="comprehensive", comprehensiveness=1)
    a, b = evaluate(baseline, AS_OF), evaluate(comprehensive, AS_OF)
    assert preferred([baseline, comprehensive], AS_OF) == "comprehensive"
    assert b["score"] / a["score"] == pytest.approx(1.25)
    assert 1 <= a["recency_factor"] <= 1.10
    # Equal citation-rate fixtures using exactly double the age in days.
    old = _candidate("old", first_publication_date="2024-10-07", citations=40)
    assert evaluate(old, AS_OF)["citations_per_year"] == pytest.approx(a["citations_per_year"])
    assert a["recency_factor"] > evaluate(old, AS_OF)["recency_factor"]


def test_slight_bonuses_do_not_swamp_a_large_citation_rate_advantage():
    popular = _candidate("popular", citations=100, comprehensiveness=0)
    broad = _candidate("broad", citations=20, comprehensiveness=1)
    assert preferred([broad, popular], AS_OF) == "popular"


def test_zero_citation_new_datasets_are_rankable_and_unknown_counts_are_not_zero():
    zero = _candidate("known_zero", citations=0)
    unknown = _candidate("unknown", citations=None, comprehensiveness=1)
    ranked = rank_candidates([unknown, zero], AS_OF)
    assert ranked[0]["dataset_id"] == "known_zero" and ranked[0]["score"] == 0
    assert ranked[1]["status"] == "metadata_incomplete" and ranked[1]["score"] is None
    assert preferred([unknown], AS_OF) is None
    assert preferred([zero, unknown], AS_OF) is None


def test_age_floor_limits_unstable_first_day_rates():
    same_day = _candidate(first_publication_date=AS_OF, citations=1)
    row = evaluate(same_day, AS_OF)
    assert row["status"] == "ranked" and row["age_years"] == 0
    assert row["citations_per_year"] < 13


def test_admission_always_precedes_popularity():
    valid = _candidate("valid")
    unsuitable = _candidate("wrong_assay", citations=999999, admission="rejected",
                             admission_reason="Synthetic wrong stage and quantity")
    pending = _candidate("unreviewed", citations=999999, admission="pending")
    assert preferred([unsuitable, pending, valid], AS_OF) == "valid"
    assert preferred([pending, unsuitable], AS_OF) is None


@pytest.mark.parametrize("changes", [
    {"first_publication_date": ""}, {"first_publication_date": "2027-01-01"},
    {"first_publication_date": "2025"}, {"citation_snapshot": ""},
    {"citation_snapshot": "2027-01-01"}, {"citation_snapshot": "2024-01-01"},
    {"citation_provider": ""}, {"publication_id": ""}, {"comprehensiveness": None},
])
def test_incomplete_or_impossible_metadata_is_visible(changes):
    result = evaluate(_candidate(**changes), AS_OF)
    assert result["status"] == "metadata_incomplete" and result["score"] is None and result["gap"]


def test_comparisons_refuse_mixed_scopes_providers_and_snapshot_dates():
    candidate = _candidate()
    for changes in ({"scope": "other_assay"}, {"citation_provider": "other-provider"},
                    {"citation_snapshot": "2026-10-06"}):
        with pytest.raises(ValueError):
            rank_candidates([candidate, replace(candidate, dataset_id="other", **changes)], AS_OF)
    with pytest.raises(ValueError, match="unique"):
        rank_candidates([candidate, candidate], AS_OF)


def test_zero_rate_ties_prefer_comprehensive_then_newer_data():
    old = _candidate("older", citations=0, first_publication_date="2020-01-01")
    new = _candidate("newer", citations=0)
    assert preferred([old, new], AS_OF) == "newer"
    assert preferred([new, replace(old, comprehensiveness=1)], AS_OF) == "older"
    assert preferred([], AS_OF) is None


@pytest.mark.parametrize("changes", [{"citations": -1}, {"citations": True}, {"citations": 1.5},
                                     {"comprehensiveness": 2}, {"comprehensiveness": float("nan")},
                                     {"comprehensiveness": True}, {"admission": "maybe"},
                                     {"admission_reason": ""}, {"scope": ""}])
def test_bad_candidate_contracts_fail(changes):
    with pytest.raises(ValueError):
        _candidate(**changes)


@pytest.mark.parametrize("changes", [{"version": 2}, {"version": True}, {"minimum_age_days": 0},
                                     {"recency_bonus": -1}, {"recency_decay_years": 0},
                                     {"comprehensiveness_bonus": float("inf")}])
def test_bad_policy_contracts_fail(changes):
    with pytest.raises(ValueError):
        SelectionPolicy(**changes)


def test_policy_can_be_frozen_and_explicitly_changed():
    policy = SelectionPolicy(recency_bonus=.05)
    assert SelectionPolicy(**policy.to_dict()) == policy
    assert evaluate(_candidate(), AS_OF, policy)["recency_factor"] < evaluate(_candidate(), AS_OF)["recency_factor"]
