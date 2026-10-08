# Human source-gene foundation, 2026-10-08

Controlled first partition of 64.21, not full completion or an activated gene
space. `build.ipynb` records the actual build and strict independent replay.
It reuses the verified GTEx v10 processed median-TPM file and UniProt 2026_03
reviewed human Ensembl cross-references. This source-scoped universe contains
58,988 canonical Ensembl genes; it is not a claim of full genome coverage.
Installed protein reference rows are never renamed or used as the gene universe.

All 59,033 original gene records retain source position, versioned identifier,
qualifier, source-provided name and both selected context measurements. The
45 qualified pseudoautosomal Y records retain their 90 measured cells separately
and cannot overwrite ordinary-gene values. All 118,066 source cells replay exactly
from their text; ordinary-gene nodes retain 117,976 of them. Remaining qualified
cells are retained evidence, not discarded data. Genuine median-TPM zero stays
observed; missing and unmapped stay distinct. Version/source positions and mapping
counts are not eligible biological features. Node names are source-provided,
without claiming HGNC validation.

Reviewed mapping census: 22,320 explicit candidate pairs; 39,630 source genes
unmapped, 19,208 with one reviewed protein, 150 with multiple proteins. 1,713
source genes link to a protein shared with another reference gene. Across the
full reference, 1,615 proteins link to multiple genes and 2,685 mapping genes
are outside the source universe. Full source cross-reference identifiers and
isoform values remain available. No ambiguous mapping is chosen by row order,
and no protein measurement is projected or averaged into gene rows.

Per-column metadata and typed measurement provenance retain species, source
hash/version/context, median-TPM units, evidence grade, candidate benchmark
status and licensing/redistribution gaps. GTEx cultured adult fibroblast medians
are not validated HFF measurements. Independent biological accuracy is not
admitted; zero expression is not a negative phenotype label.

Primary context/reference checks:
[GTEx download catalogue](https://gtexportal.org/home/downloads/adult-gtex/)
identifies the v10 GENCODE39 gene model and v11's annotation update with no new
samples/donors. The [GENCODE FAQ](https://www.gencodegenes.org/pages/faq.html)
explains qualified PAR_Y identifiers for releases 25–43. The GTEx catalogue
search index describes use with acknowledgement; direct static catalogue/license
requests return a JavaScript shell. Specific pack redistribution clearance remains
unresolved. Never substitute an article's license or infer a data license from
public availability. [UniProt terms](https://www.uniprot.org/help/license) are a
separate reference-rights source, not permission for GTEx measurements.

GTEx v11 remains a comparison candidate under 65.04. This reproducibility
foundation preserves the currently acquired v10 source; it does not select or
promote v10 as the best future gene space. Suitable sources still require the
publication-rate/recency/comprehensiveness review and fixed eligibility denominators.

618 relevant tests passed; the earlier 436-check run overlaps and is not additive.
Controls cover exact source-text conversion, zero/missingness, qualified-Y
non-overwrite, dimensions, canonical IDs, annotation collisions, source order,
many-to-many/reverse mappings, unreviewed exclusions, truncated gzip refusal,
packs, host/source reviews, provenance, registry/docstrings, ground truth/splits,
artifacts and layout/publishing compatibility. One existing module-execution
warning remains. Two real-build wrapper errors occurred before output writes:
missing SourceFile role and an unsupported role name. Corrected to the existing
`processed_input`/`mapping_reference` contract, then build/replay passed. Original
wrapper code and diagnostics are preserved in the adjacent preflight packet.

All inspected installed host and parasite input hashes remain unchanged. Serial
CPU work ran under a 2 GB systemd cap, no GPU. Full reference/admission review,
graph construction/order, independent opening and distributable pack gates remain
open; no synthetic graph or deployment accuracy is asserted.
