# Host symbol correction · first builder partition

66.04 remains open. `host.uniprot_index` now withholds a primary symbol naming
multiple reviewed proteins, with the same public schema, missing/truncated-input
handling and unchanged Ensembl mapping. Repeated records of the same accession
remain valid; reviewed row order cannot decide a protein's gene score.

`deposits.k562_rhoptry_evidence` retains every original worksheet row/symbol and
all three source score fields before protein projection. `k562_rhoptry_screen`
uses this reader and projects only uniquely mapped symbols. A large ambiguous
hit is preserved as gene evidence, never duplicated/assigned to one protein.

`audit.ipynb` compares the new parser with the archived pre-fix implementation
on identical human/mouse mapping inputs. Every excluded symbol is genuinely
multi-accession; all retained symbol mappings and every Ensembl mapping match
exactly. The full primary screen verifies **20,010 original gene records**:
18,700 uniquely mapped, 39 ambiguous, 1,271 unmapped. Every original symbol,
source row and all three scores are exact. The corrected **candidate** protein
projection retains all 18,700 unambiguous rows, with all three values exactly
matching installed data. All 39 original ambiguous scores and their alternative
accessions remain in `rhoptry_gene_evidence.parquet`; legacy first-assigned
accessions are recorded in `ambiguous_legacy_projection_review.parquet`.

**574 relevant checks passed; two existing skips.** New tests establish mapping
order independence and source-gene non-loss when protein projection withholds
ambiguous/unmapped genes. Host/species/deposit/provenance/ground-truth, leakage
and documentation checks also pass. Executed source replay verifies full counts,
values and immutable input identities; exact implementation code is archived.

**Installed protein tables/caches are unchanged.** The 39 legacy wrong protein
projections remain pending an explicit migration with withdrawal provenance,
preserved gene evidence and affected non-loss/layout/benchmark/calibration checks.
GT-HOST-TX-02 separately records PAR_Y/Atlas-cutoff and installed mapping-lineage
gaps. No biological accuracy, host-pack completion or release claim.
