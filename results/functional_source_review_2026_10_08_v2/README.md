# Functional source and EC nomenclature review

Executed `review.ipynb` pins public ENZYME metadata (02-Sep-2026), source receipts
and every installed EC annotation in both parasite organisms. Original annotations
remain unchanged. A unique transfer chain may resolve a nomenclature identifier;
split transfers, deleted, preliminary, missing and malformed entries cannot create
new gene assignments. A profile is eligible only if every assigned source EC resolves
uniquely to an active entry. All major-class memberships are retained together.

Complete profiles: 1,226 Toxoplasma genes, 1,050 Plasmodium curated-source genes,
and 1,055 Plasmodium orthology-derived genes. Unannotated or unresolved profiles
remain unknown. The stricter v2 review excludes two malformed orthology profiles;
the earlier review is preserved. Metadata defines nomenclature; gene-level
experimental/curated/predicted grades, assignment release and homology/source
independence remain unresolved. **No independent activity benchmark admitted.**

This is acquisition and source review, not a claim that the installed gene datasets
are the best literature choices or a validated dataset promotion. Raw inputs stay
in the external functional-nomenclature archive; retained receipts record URLs,
retrieval dates, bytes and hashes. Original code is copied under `code/` for replay.
The separate current domain-name audit is
`../functional_domain_metadata_2026_10_08_v4/`.
