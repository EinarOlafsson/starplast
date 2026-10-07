# Dataset selection

When several datasets answer the same biological question, Starplast prefers
annual citation rate, with modest recency and comprehensiveness bonuses. This is
the user's selection rule of 2026-10-07 and is implemented by
`starplast.dataset_selection.SelectionPolicy`, version 1.

First admit only datasets with the right organism/host, cell type, stage, assay,
quantity, units and mapping, usable source data and passing quality checks.
Reviews, retracted work, incompatible conditions and unsuitable measurements do
not win through popularity. Measured, transferred, predicted and derived evidence
retain their roles. Complementary datasets and replicates remain accessible.

For a comparable admitted dataset:

```
age_years = days_since_first_publication / 365.2425
citations_per_year = citations / (max(days_since_first_publication, 30) / 365.2425)
recency_factor = 1 + 0.10 * exp(-age_years / 5)
comprehensiveness_factor = 1 + 0.25 * comprehensiveness
selection_score = citations_per_year * recency_factor * comprehensiveness_factor
```

The 30-day floor avoids unstable first-day extrapolations. It does not replace
the actual publication age reported in the audit. Comprehensiveness must be a
verified fraction in [0, 1] with a common eligible population/assay scope. Record
when installed storage coverage is only a proxy; imputed values and an unmeasured
universe cannot establish assay coverage. Context breadth, replicate depth and
assay resolution need explicit common definitions if used instead.

Rankings retain total citations, first publication date, citation rate, both
bonuses, score, provider and citation snapshot date. Compare one provider and one
snapshot date at a time. Missing citations, dates or comprehensiveness stay
unknown, so automatic selection returns no winner until those gaps are resolved.
Known zero-citation datasets remain eligible; zero-score ties prefer more
comprehensive data, then more recent publication. Dataset IDs break exact ties
deterministically. Citations from several papers describing one dataset are not
summed without an explicit publication association policy.

`DatasetCandidate` records the biological admission decision and its reason.
`evaluate()` exposes each factor or an explicit gap; `rank_candidates()` retains
pending/rejected/incomplete choices; `preferred()` returns a dataset ID only when
an admitted, fully described candidate exists. Weights can be changed explicitly
and frozen with the comparison. Citation preference is not a measured accuracy or
inference-capacity score; those remain tied to the ground-truth tests.

Europe PMC's documented citation counts depend on its indexed citation coverage;
record the provider rather than treating these as universal counts. Publication
and citation metadata can be resolved from its [REST service](https://europepmc.org/RestfulWebService).
The service and coverage interpretation are documented in its
[search help](https://europepmc.org/help).

Current selection audit and promotion gates are tracked in
[instruction 65](../instructions/open/65_dataset_selection_audit.md).

The [October 7 frozen audit](https://github.com/EinarOlafsson/starplast/blob/nightly/results/dataset_selection_2026_10_07/README.md)
accounts for all 162 sources and searches all 297 slots. It resolves publication
metadata for 95 sources and retains explicit gaps for the remaining source roles.
All literature hits and installed-output diagnostics require biological admission;
no dataset or runtime default was replaced. First pages by total citations and
publication date do not guarantee finding the highest annual citation rate;
truncated searches and source/deposit reviews remain open.
