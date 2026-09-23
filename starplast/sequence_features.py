"""Frozen ESM protein representations with residue-complete window pooling.

The optional encoder downloads a pinned public model once. No task labels enter
these features. Long proteins are tiled rather than truncated; overlapping token
representations receive reciprocal coverage weights before residue averaging.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

MODEL = "facebook/esm2_t6_8M_UR50D"
REVISION = "c731040fcd8d73dceaa04b0a8e6329b345b0f5df"


def windows(length, size=1000, overlap=128):
    """Return half-open windows and reciprocal overlap weights for each residue."""
    if length < 1 or not 0 <= overlap < size or size > 1022:
        raise ValueError("require length > 0 and 0 <= overlap < size <= 1022")
    intervals = []
    start = 0
    while start < length:
        end = min(start + size, length)
        intervals.append((start, end))
        if end == length:
            break
        start += size - overlap
    coverage = np.zeros(length, dtype=int)
    for start, end in intervals:
        coverage[start:end] += 1
    return [(start, end, 1.0 / coverage[start:end]) for start, end in intervals]


def cache_key(sequence, model=MODEL, revision=REVISION, size=1000, overlap=128):
    """Identity includes exact sequence, weights revision and pooling parameters."""
    return hashlib.sha256(json.dumps([sequence, model, revision, size, overlap,
                                      "overlap-corrected-residue-mean-v1"]).encode()).hexdigest()


def encode(sequences, cache_dir, device="cpu", model=MODEL, revision=REVISION,
           size=1000, overlap=128, batch_tokens=2048, log=print):
    """Encode a gene-to-protein mapping; resume safely from per-sequence caches.

    Returns numeric gene features and a provenance table. The 8M model emits 320
    coordinates; these are latent features, not named biological properties.
    The GPU is optional and evaluation mode disables dropout and gradients.
    """
    try:
        import torch
        from transformers import AutoModel, AutoTokenizer
    except ImportError as exc:
        raise ImportError("Install starplast[sequence] to generate protein embeddings") from exc
    if not revision or revision in {"main", "latest"}:
        raise ValueError("pin a model revision to make cached features reproducible")
    if batch_tokens < size + 2:
        raise ValueError("batch_tokens must fit a complete window plus special tokens")
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    tokenizer = encoder = None
    rows, provenance = [], []
    for i, (gene, sequence) in enumerate(sorted(sequences.items())):
        if not sequence or any(c not in "ACDEFGHIKLMNPQRSTVWYXBZUO" for c in sequence):
            provenance.append({"gene_id": gene, "status": "unsupported_sequence"})
            continue
        key = cache_key(sequence, model, revision, size, overlap)
        path = cache / (key + ".npy")
        if path.is_file():
            vector = np.load(path, allow_pickle=False)
        else:
            if encoder is None:
                tokenizer = AutoTokenizer.from_pretrained(model, revision=revision)
                encoder = AutoModel.from_pretrained(model, revision=revision, add_pooling_layer=False).to(device).eval()
            pieces = windows(len(sequence), size, overlap)
            total = np.zeros(encoder.config.hidden_size, dtype=np.float64)
            batch_size = max(1, batch_tokens // (size + 2))
            for start in range(0, len(pieces), batch_size):
                batch = pieces[start:start + batch_size]
                inputs = tokenizer([sequence[a:b] for a, b, _ in batch], return_tensors="pt", padding=True)
                inputs = {name: tensor.to(device) for name, tensor in inputs.items()}
                with torch.inference_mode():
                    hidden = encoder(**inputs).last_hidden_state.detach().float().cpu().numpy()
                for row, (a, b, weights) in zip(hidden, batch):
                    total += np.sum(row[1:b-a+1] * weights[:, None], axis=0)
            vector = (total / len(sequence)).astype(np.float32)
            tmp = path.with_suffix(".tmp.npy")
            np.save(tmp, vector, allow_pickle=False)
            tmp.replace(path)
        if vector.ndim != 1 or not np.isfinite(vector).all():
            raise ValueError(f"invalid sequence feature cache: {path}")
        rows.append({"gene_id": gene, **{f"esm_{j:03}": float(v) for j, v in enumerate(vector)}})
        provenance.append({"gene_id": gene, "status": "encoded", "sequence_length": len(sequence),
                           "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                           "cache_key": key, "model": model, "revision": revision,
                           "window": size, "overlap": overlap, "residue_coverage": 1.0})
        if (i + 1) % 100 == 0:
            log(f"Encoded/cached {i+1}/{len(sequences)} proteins", flush=True)
    return pd.DataFrame(rows), pd.DataFrame(provenance)


def main():
    """Generate ESM features from the bundled Toxoplasma coding sequences."""
    import argparse
    from .structure_catalog import protein_sequences, DATA
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cache", type=Path, default=Path.home()/".cache/starplast/esm")
    parser.add_argument("--device", default="cpu")
    args = parser.parse_args()
    genes = set(pd.read_parquet(DATA/"nodes.parquet", columns=["gene_id"]).gene_id)
    sequences = {g:s for g,s in protein_sequences().items() if g in genes}
    features, manifest = encode(sequences, args.cache, device=args.device)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    features.to_parquet(args.output, index=False)
    manifest.to_parquet(args.output.with_name(args.output.stem+"_manifest.parquet"), index=False)
    print(f"Saved {len(features)} protein representations to {args.output}")


if __name__ == "__main__":
    main()
