# Journal/preprint rhoptry screen equivalence

Executed table comparison verifies all 20,010 gene symbols and all three score
fields exactly: combined score, beta and Wald FDR. The journal workbook member,
sheet name and SHA256 are recorded. The historical loader and `%.10g` serialization
reproduce all 18,739 installed protein rows in all three fields exactly.

The independent mapping census records 18,700 unambiguous symbols, 39 ambiguous
symbols and 1,271 unmapped symbols. The legacy symbol dictionary silently selects
the first reviewed accession for those 39; correction is separately controlled
under 66.04. The source note now states that gap rather than presenting the rows
as a fully unambiguous host gene space. Raw gene-level scores remain preserved.

PubMed verifies journal PMID/DOI and UpdateOf the registered preprint. The publisher
Dataset EV1 link establishes file association. The primary bioRxiv history confirms
first publication 2025-10-17 and v2 2025-10-21; journal publication is 2026-09-25.
Journal/preprint versions are one experiment. Citation counts are retained
separately, without summing overlaps or resetting source age to 2026. The journal
and preprint primary snapshots both report zero citations. Lineage retrieval
timestamps crossed UTC midnight to October 8, still October 7 locally.

Decision: upgrade the canonical citation to the peer-reviewed journal while
retaining exact installed values, v2 file/download provenance and first-publication
lineage. Dataset catalogue and generated source script are refreshed. No mapping
or runtime measurement changes accompany the citation correction. The initial
resource-fork failure and source archive receipt are retained separately.
548 host/source checks and 625 affected registry/selection/catalogue/review checks
passed. Exact version comparison and missingness differences have regression tests.
