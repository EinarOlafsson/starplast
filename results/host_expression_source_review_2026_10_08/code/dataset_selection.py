"""Rank comparable, admitted datasets using dated publication and coverage evidence.

Citation rate is the main preference; bounded bonuses reward recency and
comprehensiveness. These are selection preferences, not accuracy estimates.
Every candidate must already pass a biological admission check for the same
question. Missing bibliometrics and coverage cannot silently become zeros.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import date
import math


@dataclass(frozen=True)
class SelectionPolicy:
    """Versioned defaults for the user's citation-rate, recency and coverage rule.

    Citation rate uses actual days since first publication, with a 30-day floor
    to prevent unstable rates for papers published only hours ago. Recency adds
    at most 10%, decaying over five years; comprehensiveness adds at most 25%.
    Weights are preferences and can be changed explicitly, retaining the policy.
    """

    version: int = 1
    minimum_age_days: float = 30.0
    recency_bonus: float = 0.10
    recency_decay_years: float = 5.0
    comprehensiveness_bonus: float = 0.25

    def __post_init__(self):
        if type(self.version) is not int or self.version != 1:
            raise ValueError("Unsupported selection-policy version")
        for name in ("minimum_age_days", "recency_bonus", "recency_decay_years", "comprehensiveness_bonus"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError("Policy weights must be finite numbers")
        if self.minimum_age_days <= 0 or self.recency_decay_years <= 0:
            raise ValueError("Age floor and recency decay must be positive")
        if self.recency_bonus < 0 or self.comprehensiveness_bonus < 0:
            raise ValueError("Preference bonuses must be nonnegative")

    def to_dict(self) -> dict:
        """Return the weights to freeze alongside an audit or selection decision."""
        return asdict(self)


@dataclass(frozen=True)
class DatasetCandidate:
    """One source assessed for one scope, with explicit biological admission.

    Comprehensiveness is a verified fraction in [0, 1] for a common eligible
    population/assay scope. Storage coverage may be used only when explicitly
    labeled as its proxy; context breadth or replicate depth needs a separately
    defined common measure. Publication metadata includes provider and snapshot.
    """

    dataset_id: str
    scope: str
    publication_id: str = ""
    first_publication_date: str = ""
    citations: int | None = None
    comprehensiveness: float | None = None
    citation_provider: str = ""
    citation_snapshot: str = ""
    admission: str = "pending"
    admission_reason: str = ""

    def __post_init__(self):
        if not self.dataset_id or not self.scope:
            raise ValueError("Dataset and comparison scope are required")
        if self.admission not in {"accepted", "pending", "rejected"}:
            raise ValueError("Admission must be accepted, pending or rejected")
        if self.admission != "pending" and not self.admission_reason.strip():
            raise ValueError("Admission decisions need a reason")
        if self.citations is not None and (type(self.citations) is not int or self.citations < 0):
            raise ValueError("Citation count must be a nonnegative integer or unknown")
        if self.comprehensiveness is not None:
            if isinstance(self.comprehensiveness, bool) or not isinstance(self.comprehensiveness, (int, float)):
                raise ValueError("Comprehensiveness must be a numeric fraction")
            if not math.isfinite(self.comprehensiveness) or not 0 <= self.comprehensiveness <= 1:
                raise ValueError("Comprehensiveness must be between zero and one")


def evaluate(candidate: DatasetCandidate, as_of: str, policy: SelectionPolicy | None = None) -> dict:
    """Return the score and all factors, or an explicit admission/metadata gap.

    Score = citations/year × (1 + recency bonus × exp(-age/decay))
    × (1 + comprehensiveness bonus × comprehensiveness).
    Zero-citation ties use comprehensiveness, then newer first publication.
    Counts observed after the comparison date cannot be backdated.
    """
    policy = policy or SelectionPolicy()
    comparison_date = date.fromisoformat(as_of)
    result = {**asdict(candidate), "policy_version": policy.version, "as_of": as_of,
              "status": candidate.admission, "age_years": None, "citations_per_year": None,
              "recency_factor": None, "comprehensiveness_factor": None, "score": None, "gap": ""}
    if candidate.admission != "accepted":
        result["gap"] = candidate.admission_reason or "Biological suitability is not yet reviewed"
        return result
    required = ("publication_id", "first_publication_date", "citations", "comprehensiveness",
                "citation_provider", "citation_snapshot")
    missing = [name for name in required if getattr(candidate, name) is None or getattr(candidate, name) == ""]
    if missing:
        result.update(status="metadata_incomplete", gap=", ".join(missing))
        return result
    try:
        publication_date = date.fromisoformat(candidate.first_publication_date)
        snapshot_date = date.fromisoformat(candidate.citation_snapshot)
    except ValueError:
        result.update(status="metadata_incomplete", gap="Publication or citation snapshot date is invalid")
        return result
    if publication_date > comparison_date or snapshot_date > comparison_date or snapshot_date < publication_date:
        result.update(status="metadata_incomplete", gap="Publication/snapshot dates do not fit comparison date")
        return result
    age_days = (comparison_date - publication_date).days
    age_years = age_days / 365.2425
    rate = candidate.citations / (max(age_days, policy.minimum_age_days) / 365.2425)
    recent = 1 + policy.recency_bonus * math.exp(-age_years / policy.recency_decay_years)
    coverage = 1 + policy.comprehensiveness_bonus * candidate.comprehensiveness
    result.update(status="ranked", age_years=age_years, citations_per_year=rate,
                  recency_factor=recent, comprehensiveness_factor=coverage, score=rate * recent * coverage)
    return result


def rank_candidates(candidates, as_of: str, policy: SelectionPolicy | None = None) -> list[dict]:
    """Rank one scope, retaining pending/rejected/incomplete candidates afterward.

    Rankings use one citation provider and one snapshot date. Mixing providers
    or historical counts silently would create a preference for index coverage.
    Dataset IDs must be unique; multiple papers from one dataset need an explicit
    dataset/publication association rather than summed citation counts.
    """
    candidates = list(candidates)
    if len({c.scope for c in candidates}) > 1:
        raise ValueError("Only candidates for the same biological scope can be ranked together")
    if len({c.dataset_id for c in candidates}) != len(candidates):
        raise ValueError("Candidate dataset IDs must be unique")
    rows = [evaluate(c, as_of, policy) for c in candidates]
    ranked = [r for r in rows if r["status"] == "ranked"]
    if len({(r["citation_provider"], r["citation_snapshot"]) for r in ranked}) > 1:
        raise ValueError("Comparison needs a common citation provider and snapshot date")
    return sorted(rows, key=lambda r: (
        r["status"] != "ranked", -(r["score"] or 0),
        -(r["comprehensiveness"] or 0), -date.fromisoformat(r["first_publication_date"]).toordinal()
        if r["status"] == "ranked" else 0, r["dataset_id"]))


def preferred(candidates, as_of: str, policy: SelectionPolicy | None = None) -> str | None:
    """Return a winner only when every admitted alternative has comparable metadata.

    An admitted source with unknown citations or coverage could outrank a fully
    described source. Report that gap instead of treating it as a loser.
    Pending and rejected sources remain outside the admitted comparison.
    """
    rows = rank_candidates(candidates, as_of, policy)
    if any(r["admission"] == "accepted" and r["status"] != "ranked" for r in rows):
        return None
    return rows[0]["dataset_id"] if rows and rows[0]["status"] == "ranked" else None
