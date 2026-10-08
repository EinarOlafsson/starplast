# Browse datasets

Open **Tools → Datasets**. The browser shares the catalog window with **Slots**.
Filter by organism or host, entity unit, evidence family, measured context and
availability. Search also matches source names, columns, labels and recorded gaps.

Selecting a source opens its coverage and provenance card alongside its stored
values. **Previous** and **Next** reach every original row in pages of 200.
Hover over a shortened cell or heading to read its full text. Double-click a gene
identifier to open that gene in the corresponding organism's map. Host proteins,
metabolites and pair records retain their own identifiers and units.

Coverage counts nonmissing stored values in the declared table. Zero and `False`
are preserved; missing values retain unknown interpretation. Pair counts retain
an unknown assay denominator where none was recorded. Source cards show declared
derivations, publication links and supplied lineage; missing units, mapping,
source dates, versions and calibration are explicit. Biological accuracy requires
its own benchmark.

Session imports appear as separate sources with their preprocessing record. The
original source-file identity remains unavailable when the import record lacks it.
Importing values refreshes an open browser without changing installed cache files.

For Python callers, `starplast.dataset_space.DatasetSpace` accepts an explicitly
qualified inventory and local tables. `filter_rows`, `facets`, `source_card` and
`entities(row, start=0, limit=200)` expose the same scopes, immutable evidence cards
and original values used by the desktop browser.

Acceptance is recorded by `python -m scripts.audit_dataset_browser --out NEW_DIRECTORY`.
The executed notebook freezes local input/code hashes and checks source cards,
filters, exact stored values, pagination, session imports and gene navigation.
Use a fresh directory for each replay.

Gene search resolves exact recorded aliases within the selected organism. Multiple
alias or text matches show a candidate list; select an accession to open its gene.
The gene card’s **All evidence by biological question** link opens a condensed
view. **Show all evidence** includes every stored column and missing value.
Each row shows its question, context, quantity unit, declared origin and gaps;
source selection opens the original coverage/provenance card. Label and class
buttons reach existing recorded results when the current table matches the
installed context. Changed tables require their own validation.

Python callers can use `gene_evidence.build_resolver` for exact scoped aliases
and `gene_evidence.evidence_rows` for original per-column records. Gene entry
acceptance runs with `python -m scripts.audit_gene_evidence --out NEW_DIRECTORY`.
