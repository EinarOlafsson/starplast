# Functional domain nomenclature metadata

Canonical metadata audit for instruction 68.02. InterPro release 110.0 and Pfam release 38.2 metadata provide current names/types for installed identifiers. These are nomenclature, not newly inferred gene functions or independent biological accuracy.

All 34,983 installed organism/target/gene/identifier memberships and source descriptions are retained exactly. 8,043 unique source identifiers: 7,926 current entries, 112 missing-current InterPro identifiers, and 5 officially withdrawn Pfam identifiers. Unknown names/types remain null. Withdrawal forwarding identifiers are preserved as history and never applied to genes.

The companion InterPro names.dat matches entry.list exactly. All selected metadata records were independently replayed from literal source rows/blocks; declared complete release counts match parsed tables (54,922 InterPro; 30,134 current and 1,121 retired Pfam). Original node files retain their input hashes.

The official [InterPro data license](https://interpro-documentation.readthedocs.io/en/latest/license.html) covers downloadable InterPro and Pfam data under CC0-1.0. The archived data version says Pfam used UniProtKB 2025_03, while its current about page says 2026_01: this discrepancy remains explicit. Original gene-assignment model releases and experimental activity evidence are unresolved.

Public metadata sources and headers/dates/hashes are in source_receipts.json and status_receipts.json. Original downloads remain external; this repository retains the installed-ID subset and complete membership audit. The executed notebook and frozen code reproduce the result. Earlier acquisition and replay diagnostics are preserved separately; the v2 independent comparison correctly caught 102 trailing-whitespace differences, corrected by retaining literal source descriptor text.
