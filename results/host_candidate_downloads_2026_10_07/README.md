# Verified RBC candidate acquisition

The executed retrieval notebook downloaded the PRIDE processed archive and two
DOI-matched ACS Figshare supplements. Deposit SHA1/byte count and publisher MD5/
byte counts match; SHA256 identities and requested URLs are in `downloads.json`.
Primary Figshare article JSON supplies the association and license. Original
files remain in the external source archive. Instrument RAW files were not fetched.

The inspection notebooks establish actual sheet, fraction, donor, identifier and
quantity fields before parsing. The publisher table includes four donor columns
for each of whole-cell and white-ghost copy counts. Numeric source tables are not
committed or bundled in the wheel. The reviewed mapping and incumbent comparison
are in `../rbc_candidate_review_2026_10_07_v2/`; no acquisition automatically
promotes or replaces an installed source.
