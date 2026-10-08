# Primary host candidate discovery

Executed discovery/follow-up notebooks retain the PRIDE project, file metadata
and PubMed journal-to-preprint linkage. The first PRIDE file endpoint returned
an empty successful response; the follow-up explicitly rejects that as evidence
of an empty deposit. The correct project-files endpoint returned 97 records:
96 instrument RAW files and one processed MaxQuant archive. Only the processed
archive was selected for download.

The PubMed XML verifies the journal rhoptry DOI and UpdateOf the recorded preprint.
This establishes article lineage, not data equivalence. The subsequent direct
table comparison is in `../rhoptry_journal_review_2026_10_07_v2/`.

The old `rhoptry_publication_linkage.json` correctly states equivalence was unknown
at discovery time; it is preserved rather than rewritten after later evidence.
