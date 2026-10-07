# Measurement provenance

`starplast.provenance` provides a content-addressed trace for every inventory
source/output, including explicit gaps. It does not certify legacy source
associations merely because a file or publication identifier exists.

Each `MeasurementTrace` records the registered exploration scope, output organism,
storage unit and column, actual assay species with its own verification status,
measurement unit, context, evidence grade, publication association, input files,
processing steps, mapping audits, upstream sources and redistribution status.
An installed cache remains an installed cache, rather than a raw experiment.
Registry organism is not necessarily the original species of an orthology transfer.

`SourceFile` pins SHA-256 and byte count and can verify the current content.
`Transform` pins implementation code and explicit processing parameters.
`MappingAudit` counts unique source identifiers, successful mappings, ambiguous
sources withheld, unmapped sources, target proteins and many-to-one targets.
Version and assembly/namespace references are required for audited mappings;
unavailable strain or assembly audits stay visible in trace gaps. Missing values
and mapping failures never become experimental negatives.

`trace_sources()` follows declared upstream columns within their registry scope,
rejects cycles and retains unresolved or multiply declared parents. `write_traces()`
rejects duplicate output addresses and existing snapshot paths. `read_traces()`
checks each trace's content identity, including gaps and processing parameters.

## October 7 audit

The executed audit in `results/measurement_provenance_recovered_2026_10_07/` covers all
162 registered sources and 628 source/output addresses. Of these, 364 have a
located input or installed-cache identity. Most historical unit, species,
publication lineage, mapping and license audits remain unresolved: 627 addresses
have no verified measurement unit. This is a complete accounting of the gaps,
not completion of source-specific biological validation.

The selected raw-to-installed pilot is mouse BMDM baseline expression. It pins the
GEO source and reviewed UniProt mapping files, normalizes all source FPKM values
to TPM per sample, averages the three unstimulated replicates, resolves identifiers,
keeps the most expressed row for duplicate targets, serializes using the historical
ten-significant-digit TSV format and reproduces all 15,437 installed TPM values
exactly. Of 35,275 source entities, 15,531 map unambiguously, 53 are ambiguous and
19,691 remain unmapped; 58 target proteins have multiple unambiguous source genes.
Those losses are not hidden by the final table's coverage.

The preliminary schema-1 snapshot and rounding diagnostic are retained separately
with their original code copies. Schema 2 separates registry scope from actual
assay species. New source downloads and admission decisions are tracked under
instruction 66 and dataset replacement under 65.04.

To reproduce the census using an existing source-binding snapshot:

```bash
python scripts/build_measurement_provenance.py \
  --root "$STARPLAST_DATA" \
  --bindings results/source_recovery_pmc_2026_10_07_v2/bindings.json \
  --out results/NEW_PROVENANCE_SNAPSHOT
```

Use the dataset archive itself as `STARPLAST_DATA`, rather than its parent.
The builder checks input hashes before and after the audit and records an
executed annotated notebook. It changes no installed measurements or inferences.
