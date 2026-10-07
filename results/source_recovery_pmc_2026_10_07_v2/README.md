# Missing supplements recovered from the current PMC Cloud Service

Source: **NIH NLM NCBI PubMed Central (PMC) Article Datasets on AWS**, accessed
October 7, 2026 from <https://registry.opendata.aws/ncbi-pmc>. Current service
documentation: <https://pmc.ncbi.nlm.nih.gov/tools/pmcaws/>. The retired legacy
archive is not used; the exact current README is archived in this snapshot.

The executed notebook discovers all versions of each recorded PMCID, checks
the declared PMID, manuscript/retraction state, exact named supplemental file
and published MD5. It never picks the highest version automatically. S3 URLs
are translated to HTTPS for the identical public bucket/object. No login or
restricted access is involved. Primary JSON and version listings are retained.

**22 missing named inputs attempted; 12 recovered; 57 total processed-source
bindings** including the prior verified locations. Recovered sources are the
two human erythrocyte datasets; Plasmodium liver/phenotype transfers, IP-MS,
CDPK1 sites, lactylome, acetylome, myristoylome, palmitome, complexes and crosslinks.
Checksummed files and per-source URL/SHA256 receipts are in the external archive
under `recovered_processed_sources/`, preserving stale legacy paths unchanged.

Nine attempted sources have no unique exact filename/published version match.
`sexual_stages` has a declared PMID/PMCID mismatch and is refused. The JSON retains
every source's failure, article version and article-level license. Matching an
article does not prove its legacy transformation or grant wheel redistribution.
No installed data, quantities or strategy results were changed.

The first cloud parser diagnostic is retained separately. Seven focused current
cloud matching tests passed as part of the 380-check run; boolean flags, S3 URL
conversion, mismatched papers, foreign object paths, retractions and absent
checksums have regression coverage.
