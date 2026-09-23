# Protein sequence representations

Starplast includes 320 ESM-2 features for **8,064 Toxoplasma genes**. These are
computed from the bundled coding sequences using the frozen 8-million-parameter
ESM-2 model. They give prediction methods access to patterns in protein sequence
without requiring a known function or localization label for each protein.

An ESM coordinate is a learned representation, not a named biological trait.
It should enter a prediction model as part of a feature block. Its value alone
does not mean that a gene has a particular function.

The model is `facebook/esm2_t6_8M_UR50D`, pinned to revision
`c731040fcd8d73dceaa04b0a8e6329b345b0f5df`. Features are in
`starplast/data/esm_features.parquet`; the accompanying manifest records the
sequence hash, sequence length, model revision, and pooling recipe for each gene.
No Starplast target labels were used to train or fine-tune this model.

Long proteins are split into windows of up to 1,000 residues, overlapping by 128.
Each residue's representations are averaged across the windows that contain it,
then pooled across the full protein. Every residue contributes a total weight of
one. This preserves coverage, though a window cannot capture interactions with
residues outside its context. These features are not a full-length structure model.

Invalid coding sequences, internal stop codons and conflicting translations are
withheld. The 76 genes without a valid bundled translation retain missing sequence
features; they do not receive a zero vector.

## Regenerate or use the features

The bundled feature table needs only pandas to read. Generating new representations
requires optional dependencies and downloads the pinned model:

```bash
pip install 'starplast[sequence]'
starplast-encode-proteins --output /path/to/esm_features.parquet --device cpu
```

Use `--device cuda` for a compatible NVIDIA installation. Encoding runs in evaluation
mode with gradients disabled. Per-sequence caches under `~/.cache/starplast/esm`
allow interrupted runs to resume. A different sequence, model revision or pooling
recipe produces a different cache key.

```python
from pathlib import Path
import pandas as pd
import starplast

path = Path(starplast.__file__).parent / 'data'
features = pd.read_parquet(path / 'esm_features.parquet')
```

The small model keeps regeneration practical. Larger models or structure-derived
representations are candidates for further benchmarks, not assumed improvements.
See the [ESM model card](https://huggingface.co/facebook/esm2_t6_8M_UR50D) and
[official ESM repository](https://github.com/facebookresearch/esm).
