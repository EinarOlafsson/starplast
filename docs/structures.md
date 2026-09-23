# Local AlphaFold 3 models

Starplast can index a local AF3 collection, map the models to genes, and add
confidence and geometry features to the gene table. The coordinates stay on your
drive. The package contains a small feature table and a source manifest.

The September 2026 inventory found sequence-verified models for **1,210 of 8,140
Toxoplasma genes** in Einar Olafsson's local collection. This is the portion that
could be matched exactly to the bundled coding sequences in the available synced
folders. It is not a claim that only these genes have been folded. Unmatched
accessions, sequence-version differences and unsynced files remain unresolved.

## Build an index

```bash
pip install 'starplast[structures]'
starplast-index-structures /path/to/af_output /path/to/af3_structures \
  --uniprot-map /path/to/uniprot.csv \
  --features /path/to/af3_features.parquet
```

The optional UniProt CSV has `gene_nr` and `uniprot` columns. Gene IDs are resolved
through the identity table's current and previous accessions. A model then needs
an exact full-length sequence match or an exact, uniquely placed fragment within
that gene's translated CDS. Repeated sequences and ambiguous gene assignments
are withheld. An accession match by itself is insufficient.

The default index is `~/.cache/starplast/structure_catalog.parquet`.
`--index` chooses a different destination; set `STARPLAST_STRUCTURE_INDEX` to that
file when launching the app. Selecting a gene in the structure viewer checks the
index before trying AlphaFold DB. Full models take precedence over fragments.
Local paths are machine-specific; recreate the index when the drive moves.

## What the features mean

| Feature | Interpretation |
| --- | --- |
| Sequence coverage | Fraction of reference residues present in matched models |
| Mean and lower-quartile pLDDT | Confidence among modelled residues |
| Confident/very-confident sequence fraction | Fraction of the full reference with pLDDT ≥70/≥90 |
| Low-confidence modelled fraction | Fraction of modelled residues with pLDDT <50 |
| Confident radius of gyration | C-alpha spread of residues with pLDDT ≥70, in Å |
| Confident contacts per residue | C-alpha neighbours within 8 Å, excluding sequence separation ≤4 |
| pTM | AF3's whole-model confidence summary when available |

The radius and contact summary are only calculated from full-length models.
Low global confidence can make relative domain placement unreliable even when
individual residues have high pLDDT; geometry is a descriptive feature, not proof
of a stable conformation. pLDDT is prediction confidence, not a direct disorder
measurement or experimental validation.

Overlapping fragments contribute once per residue, using the maximum available
pLDDT. That choice is explicitly optimistic: compare it against individual-model
confidence when interpreting candidates. Missing residues remain missing. Seed
replicates and byte-identical mirrors do not count as independent evidence.
Parasite–host complexes, peptide jobs and decoys remain in the inventory but do
not enter these single-protein features.

The shipped summaries are unpublished computed data from Einar Olafsson's AF3
collection. They do not provide a public download of the underlying coordinates.
The source manifest records file hashes and sequence matches so a local result
can be traced back to a specific model.

See the [AF3 output documentation](https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md)
for confidence definitions and [Gemmi](https://gemmi.readthedocs.io/en/stable/mol.html)
for coordinate parsing.
