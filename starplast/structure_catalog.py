"""Index local AF3 predictions without copying coordinates into the package.

Models are matched to the bundled coding sequences, never by a numeric suffix
alone. Fragments retain residue ranges; complexes and controls stay out of the
single-protein feature table. Confidence describes a prediction, not an assay.
"""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path
import re

import numpy as np
import pandas as pd

DATA = Path(__file__).parent / "data"


def catalog_path() -> Path:
    """User-selected index, or the ordinary user cache location."""
    return Path(os.environ.get("STARPLAST_STRUCTURE_INDEX",
                str(Path.home() / ".cache/starplast/structure_catalog.parquet")))


def protein_sequences(cds_path=DATA / "toxodb_cds.tsv.gz") -> dict[str, str]:
    """Translate complete CDS records; omit internal stops and conflicting isoforms.

    A terminal stop is removed. Ambiguous amino acids remain X and participate in
    exact sequence matching. No missing sequence is inferred from gene names.
    """
    try:
        from Bio.Seq import Seq
    except ImportError as exc:
        raise ImportError("Install starplast[structures] to translate protein sequences") from exc
    table = pd.read_csv(cds_path, sep="\t", dtype=str)
    sequences = defaultdict(set)
    for gene, dna in table.iloc[:, :2].itertuples(index=False, name=None):
        if not isinstance(dna, str) or len(dna) % 3:
            continue
        seq = str(Seq(dna.upper()).translate()).rstrip("*")
        if seq and "*" not in seq:
            sequences[gene].add(seq)
    return {gene: next(iter(values)) for gene, values in sequences.items() if len(values) == 1}


def uniprot_aliases(mapping_path, node_ids, identity_path=DATA / "toxodb_identity.tsv"):
    """Resolve explicit UniProt/gene pairs through current and previous accessions."""
    from .identity import build_index, norm
    index = build_index(node_ids, str(identity_path), log=lambda *a: None)
    table = pd.read_csv(mapping_path, dtype=str)
    aliases = defaultdict(set)
    for number, accession in table[["gene_nr", "uniprot"]].itertuples(index=False, name=None):
        if pd.isna(number) or pd.isna(accession):
            continue
        gene = number if "_" in number else "TGME49_" + number.zfill(6)
        hit = index.lookup.get(norm(gene))
        if hit:
            aliases[accession].add(hit[0])
    return {acc: next(iter(ids)) for acc, ids in aliases.items() if len(ids) == 1}


def model_files(roots):
    """Yield root-ranked AF3 models and flat mirrors, excluding seed replicates."""
    seen = set()
    for root in map(Path, roots):
        if not root.is_dir():
            continue
        for entry in sorted(root.iterdir()):
            paths = sorted(entry.glob("*_model.cif")) if entry.is_dir() else [entry]
            for path in paths:
                if path.suffix.lower() not in {".cif", ".pdb"}:
                    continue
                resolved = path.resolve()
                if resolved not in seen:
                    seen.add(resolved)
                    yield path


def match_sequence(sequence, sequences, exact, alias=None):
    """Return unique gene, zero-based interval and match status; ambiguity abstains."""
    hits = exact.get(sequence, [])
    if len(hits) == 1:
        return hits[0], 0, len(sequence), "exact_full_sequence"
    if len(hits) > 1:
        return None, None, None, "ambiguous_sequence"
    if alias in sequences:
        full = sequences[alias]
        start = full.find(sequence)
        if start >= 0 and full.find(sequence, start + 1) < 0:
            return alias, start, start + len(sequence), "exact_fragment_sequence"
    return None, None, None, "sequence_unresolved"


def read_chain(path):
    """Read one peptide chain and residue pLDDT from an AF3 model using Gemmi.

    C-alpha B factors are AF3's residue-level confidence proxy. Coordinates are
    only used for descriptive geometry; this function does not infer function.
    """
    try:
        import gemmi
    except ImportError as exc:
        raise ImportError("Install starplast[structures] to read AF3 coordinate files") from exc
    structure = gemmi.read_structure(str(path))
    chains = []
    for chain in structure[0]:
        polymer = chain.get_polymer()
        if polymer and polymer.check_polymer_type() in (gemmi.PolymerType.PeptideL,
                                                       gemmi.PolymerType.PeptideD):
            chains.append(polymer)
    if len(chains) != 1:
        raise ValueError(f"expected one peptide chain, found {len(chains)}")
    residues = list(chains[0])
    sequence = gemmi.one_letter_code([r.name for r in residues]).upper()
    atoms = [next((a for a in r if a.name == "CA"), None) for r in residues]
    plddt = np.array([a.b_iso if a is not None else np.nan for a in atoms])
    coords = np.array([[a.pos.x, a.pos.y, a.pos.z] if a else [np.nan] * 3 for a in atoms])
    return sequence, plddt, coords


def geometry(coords, plddt):
    """Describe only confident C-alpha coordinates (pLDDT >= 70)."""
    from scipy.spatial import cKDTree
    keep = np.isfinite(coords).all(axis=1) & (plddt >= 70)
    points = coords[keep]
    if len(points) < 3:
        return {"confident_rg_angstrom": np.nan, "confident_contacts_per_residue": np.nan}
    rg = np.sqrt(np.mean(np.sum((points - points.mean(axis=0)) ** 2, axis=1)))
    pairs = cKDTree(points).query_pairs(8, output_type="ndarray")
    positions = np.flatnonzero(keep)
    contacts = np.sum(np.abs(positions[pairs[:, 0]] - positions[pairs[:, 1]]) > 4) if len(pairs) else 0
    return {"confident_rg_angstrom": float(rg),
            "confident_contacts_per_residue": float(2 * contacts / len(points))}


def inventory(roots, sequences, aliases=None, log=print):
    """Return model inventory and overlap-aware gene features from local AF3 files.

    Every coordinate file is fingerprinted. Identical mirrors are recorded but
    contribute once. Overlapping fragments use the highest residue confidence;
    whole-protein geometry is computed only from complete models. Missing residues
    stay missing rather than receiving a pLDDT of zero.
    """
    exact = defaultdict(list)
    for gene, seq in sequences.items():
        exact[seq].append(gene)
    aliases = aliases or {}
    rows, confidence, full_geometry, hashes = [], {}, {}, set()
    for i, path in enumerate(model_files(roots)):
        if i % 100 == 0:
            log(f"Indexed {i} coordinate files", flush=True)
        job = path.parent.name if path.name.endswith("_model.cif") else path.stem
        row = {"path": str(path.resolve()), "job": job, "gene_id": None,
               "model_type": "other_job", "mapping_status": "not_single_protein",
               "evidence_status": "computed_AF3_prediction"}
        # Complexes, decoys and peptide jobs are discoverable but cannot silently
        # stand in for a folded whole protein. Include only the named AF3 runs.
        if not re.match(r"(?i)^af3[_-]", job):
            rows.append(row)
            continue
        row["model_type"] = "single_protein"
        try:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            row["sha256"] = digest
            sequence, plddt, coords = read_chain(path)
            accession = re.sub(r"(?i)^af3[_-]", "", job).split("_")[0]
            gene, start, end, status = match_sequence(sequence, sequences, exact, aliases.get(accession))
            row.update(gene_id=gene, mapping_status=status, uniprot=accession,
                       residue_start=start, residue_end=end, model_residues=len(sequence),
                       sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                       mean_plddt=float(np.nanmean(plddt)))
            summary_path = path.with_name(path.name.replace("_model.cif", "_summary_confidences.json"))
            if summary_path != path and summary_path.exists():
                summary = json.loads(summary_path.read_text())
                for key in ("ptm", "iptm", "ranking_score", "fraction_disordered", "has_clash"):
                    row[key] = summary.get(key)
            row["duplicate_content"] = digest in hashes
            hashes.add(digest)
            if gene is not None:
                full_length = len(sequences[gene])
                row["model_type"] = "full_protein" if start == 0 and end == full_length else "fragment"
                row["sequence_coverage"] = (end - start) / full_length
            if gene is not None and not row["duplicate_content"]:
                values = confidence.setdefault(gene, np.full(full_length, np.nan))
                values[start:end] = np.fmax(values[start:end], plddt)
                if row["model_type"] == "full_protein":
                    previous = full_geometry.get(gene)
                    if previous is None or row["mean_plddt"] > previous[0]:
                        full_geometry[gene] = (row["mean_plddt"], geometry(coords, plddt), row.get("ptm"))
        except (OSError, ValueError, RuntimeError, IndexError) as exc:
            row["mapping_status"] = "parse_error"
            row["error"] = f"{type(exc).__name__}: {exc}"
        rows.append(row)
    features = []
    for gene, values in sorted(confidence.items()):
        measured = values[np.isfinite(values)]
        row = {"gene_id": gene, "af3_sequence_coverage": len(measured) / len(values),
               "af3_mean_plddt": float(measured.mean()),
               "af3_plddt_q25": float(np.quantile(measured, .25)),
               "af3_confident_sequence_fraction": float(np.sum(values >= 70) / len(values)),
               "af3_very_confident_sequence_fraction": float(np.sum(values >= 90) / len(values)),
               "af3_low_confidence_modelled_fraction": float(np.mean(measured < 50))}
        if gene in full_geometry:
            _, shape, ptm = full_geometry[gene]
            row.update({"af3_" + k: v for k, v in shape.items()})
            row["af3_ptm"] = ptm
        features.append(row)
    return pd.DataFrame(rows), pd.DataFrame(features)


def local_model(gene_id, index_path=None):
    """Best available indexed single-protein file, preferring complete models."""
    path = Path(index_path) if index_path else catalog_path()
    if not path.is_file():
        return None
    table = _cached_catalog(str(path), path.stat().st_mtime_ns)
    hits = table[(table.gene_id == gene_id) & table.model_type.isin(["full_protein", "fragment"])]
    if hits.empty:
        return None
    for candidate in hits.sort_values(["sequence_coverage", "mean_plddt"], ascending=False).path:
        if Path(candidate).is_file():
            return candidate
    return None


from functools import lru_cache


@lru_cache(maxsize=2)
def _cached_catalog(path, stamp):
    return pd.read_parquet(path)


def main():
    """Build a portable feature table and machine-local structure index."""
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("roots", nargs="+", type=Path)
    parser.add_argument("--uniprot-map", type=Path)
    parser.add_argument("--index", type=Path, default=catalog_path())
    parser.add_argument("--features", type=Path, required=True)
    args = parser.parse_args()
    sequences = protein_sequences()
    aliases = uniprot_aliases(args.uniprot_map, sequences) if args.uniprot_map else {}
    models, features = inventory(args.roots, sequences, aliases)
    args.index.parent.mkdir(parents=True, exist_ok=True)
    args.features.parent.mkdir(parents=True, exist_ok=True)
    models.to_parquet(args.index, index=False)
    features.to_parquet(args.features, index=False)
    summary = {"models": len(models), "genes_with_features": len(features),
               "mapping_status": models.mapping_status.value_counts().to_dict(),
               "model_types": models.model_type.value_counts().to_dict()}
    args.features.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def attach_features(nodes, feature_path=DATA / "af3_features.parquet"):
    """Attach sequence-verified AF3 features by gene ID, preserving row order/index."""
    path = Path(feature_path)
    if not path.is_file():
        return nodes
    features = pd.read_parquet(path).set_index("gene_id")
    if not features.index.is_unique:
        raise ValueError("AF3 feature table contains duplicate gene IDs")
    result = nodes.copy()
    for column in features:
        result[column] = result.gene_id.map(features[column])
    return result


if __name__ == "__main__":
    main()
