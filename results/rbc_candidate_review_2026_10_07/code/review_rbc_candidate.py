"""Audit the published RBC copy-number tables without conflating protein groups.

The frozen publisher metadata establishes the file identity and license. This
review keeps donors, fractions, original groups, zeroes and ambiguous mappings;
it never assigns a group's abundance to its first accession or substitutes copy
numbers for installed spectral counts or cell-surface accessibility.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import sys
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from starplast import host, organisms as O, paths  # noqa: E402
from starplast.dataset_selection import DatasetCandidate, evaluate  # noqa: E402
from starplast.provenance import audit_mapping  # noqa: E402


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parse_copy_numbers(path: Path) -> pd.DataFrame:
    """Read verified multirow headers and preserve every source protein group/donor."""
    raw = pd.read_excel(path, sheet_name="Table S3", header=None)
    if raw.shape[1] != 39 or raw.iloc[1, 0] != "Protein IDs":
        raise ValueError("Unexpected published Table S3 shape or identifier field")
    names = {}
    for start, context in ((31, "whole_rbc"), (35, "white_ghost")):
        expected = "Whole erythrocyte" if context == "whole_rbc" else "White ghosts"
        if raw.iloc[1, start] != expected:
            raise ValueError("Copy-number context differs from reviewed header")
        for offset in range(4):
            index = start + offset
            if str(raw.iloc[2, index]).strip().lower() != f"copy number donor {offset + 1}":
                raise ValueError("Unexpected copy-number donor header")
            names[index] = f"{context}_copies_donor{offset + 1}"
    data = raw.iloc[3:].copy()
    if data[0].isna().any() or data[0].duplicated().any():
        raise ValueError("Missing or duplicate source protein-group identifier")
    out = data[[0, 1, 2, *names]].rename(columns={0: "protein_group", 1: "protein_name", 2: "source_gene_names", **names})
    out.insert(0, "source_excel_row", out.index + 1)
    for column in names.values():
        out[column] = pd.to_numeric(out[column], errors="raise")
        if not np.isfinite(out[column]).all() or (out[column] < 0).any():
            raise ValueError("Copy counts must be finite and nonnegative; missing values stay unresolved")
    return out.reset_index(drop=True)


def group_candidates(identifier: str, reviewed: set) -> tuple[str, ...]:
    """Resolve all explicitly listed reviewed accessions; contaminants are excluded."""
    return tuple(sorted({token.strip().split("-")[0] for token in identifier.split(";")
                         if token.strip().split("-")[0] in reviewed and not token.strip().startswith(("CON__", "REV__"))}))


def review(root: Path, download_snapshot: Path, output: Path) -> dict:
    """Freeze quantity, mapping, donor QC, incumbent comparison and admission gaps."""
    if output.exists():
        raise ValueError("Use a new candidate review snapshot")
    output.mkdir(parents=True)
    records = json.loads((download_snapshot / "downloads.json").read_text())
    source = next(record for record in records if record["filename"].endswith(".xlsx"))
    file = Path(source["path"])
    if _sha(file) != source["sha256"]:
        raise ValueError("RBC file no longer matches its verified download")
    metadata_path = ROOT / "results/host_candidate_deposits_2026_10_07/rbc_pride.json"
    metadata = json.loads(metadata_path.read_text())
    if [item["accession"] for item in metadata["organisms"]] != ["NEWT:9606"]:
        raise ValueError("Candidate is not the declared human-cell dataset")
    frame = parse_copy_numbers(file)
    reference = root / host.UNIPROT_ROOT / host.UNIPROT_IDMAP["human"][1]
    with gzip.open(reference, "rt") as stream:
        next(stream)
        reviewed = {line.split("\t", 1)[0] for line in stream}
    candidates = [group_candidates(group, reviewed) for group in frame.protein_group]
    frame["reviewed_candidates"] = [";".join(values) for values in candidates]
    frame["mapping_status"] = ["unmapped" if not values else "ambiguous_group" if len(values) > 1 else "unique_reviewed_candidate"
                               for values in candidates]
    pairs = [(group, accession) for group, choices in zip(frame.protein_group, candidates) for accession in choices]
    mapping = audit_mapping(frame.protein_group, pairs, source_organism=O.HUMAN, target_organism=O.HUMAN,
                            source_namespace="published UniProt protein group", target_namespace="reviewed UniProt protein",
                            source_version="publisher SHA256:" + _sha(file),
                            target_version="reviewed reference SHA256:" + _sha(reference), reference=str(reference))
    reverse = Counter(values[0] for values in candidates if len(values) == 1)
    frame["single_protein_projection"] = [values[0] if len(values) == 1 and reverse[values[0]] == 1 else None for values in candidates]
    columns = [column for column in frame if "_copies_donor" in column]
    numeric = frame[columns]
    correlations = numeric.corr(method="spearman").round(6)
    zeroes = {name: int(numeric[name].eq(0).sum()) for name in columns}
    installed_path = Path(paths.cache_file(O.HOST_TABLES[O.HUMAN]))
    installed = pd.read_parquet(installed_path).set_index("host_id")
    scope = "human erythrocyte absolute abundance, copies per cell; four donors, whole cells and white ghosts separate"
    doi = metadata["references"][0]["doi"]
    query = {"query": "DOI:" + doi, "format": "json", "resultType": "core", "pageSize": 10}
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(query)
    retrieved = datetime.now(timezone.utc).isoformat()
    bibliometrics = {"query": query, "url": url, "retrieved_utc": retrieved, "status": "pending"}
    publication = None
    try:
        with urllib.request.urlopen(url, timeout=25) as response:
            payload = response.read(4 << 20)
        hits = json.loads(payload)
        exact = [item for item in hits["resultList"]["result"] if item.get("doi", "").lower() == doi.lower()]
        if len(exact) != 1:
            raise ValueError("No unique DOI-matched primary publication")
        publication = exact[0]
        (output / "publication.json").write_bytes(payload)
        bibliometrics.update(status="retrieved_unique_DOI", sha256=hashlib.sha256(payload).hexdigest())
    except Exception as error:
        bibliometrics.update(status="metadata_failed", error=f"{type(error).__name__}: {error}")
    candidate = DatasetCandidate(metadata["accession"], scope, publication_id="DOI:" + doi,
                                 first_publication_date=publication.get("firstPublicationDate", "") if publication else "",
                                 citations=publication.get("citedByCount") if publication else None,
                                 citation_provider="Europe PMC", citation_snapshot=retrieved[:10], comprehensiveness=None,
                                 admission="pending", admission_reason="Complementary absolute quantity; common assay denominator and redistribution scope unresolved")
    factors = evaluate(candidate, retrieved[:10])
    # Keep bibliometric factors independently visible even though biological admission gates ranking.
    if publication:
        days = (datetime.fromisoformat(retrieved).date() - datetime.fromisoformat(candidate.first_publication_date).date()).days
        bibliometrics.update(first_publication_date=candidate.first_publication_date, citations=candidate.citations,
                            citations_per_year=candidate.citations / (max(days, 30) / 365.2425))
    comparison = {}
    projected = set(frame.single_protein_projection.dropna())
    for column in ("rbc_membrane_psms", "rbc_cytoplasm_psms", "rbc_surface_log2fc"):
        if column in installed:
            members = set(installed.index[installed[column].notna()])
            comparison[column] = {"incumbent_nonmissing_proteins": len(members), "shared_unambiguous_proteins": len(members & projected),
                                  "same_quantity": False, "reason": "Absolute copy counts differ from fraction spectral counts and surface-accessibility enrichment"}
    from dataclasses import asdict
    summary = {"source": metadata["accession"], "source_species": "Homo sapiens", "source_groups": len(frame),
               "quantity": "published copies per cell", "contexts": ["whole erythrocyte", "white ghost"], "donors": 4,
               "mapping": asdict(mapping), "unique_noncolliding_protein_projection": len(projected),
               "projection_policy": "withhold multi-reviewed groups and reverse collisions; no first accession or summed groups",
               "missing_numeric_cells": int(numeric.isna().sum().sum()), "zero_cells_by_donor": zeroes,
               "zero_semantics": "Published numeric zeros retained; not certified absence or experimental negative labels",
               "donor_spearman": correlations.to_dict(), "incumbent_comparison": comparison,
               "citation_factors": factors, "bibliometrics": bibliometrics,
               "license": source["license"], "redistribution": "CC BY-NC constraints retained; source table not bundled in the wheel",
               "decision": "retain installed fraction and surface evidence; keep quantitative source as a separate candidate, not an interchangeable replacement",
               "remaining_admission": ["common scope/assay comprehensiveness denominator", "rights-compatible quantitative pack",
                                       "group-to-protein projection biological validation and source-specific held-out truth"]}
    # Original source tables and derivative numerical measurements remain in the external archive.
    export = file.parent / "group_preserving_review.parquet"
    if export.exists():
        raise ValueError("Existing candidate derivative needs an independent review")
    frame.to_parquet(export, index=False)
    summary["external_group_table"] = {"path": str(export), "sha256": _sha(export)}
    (output / "review.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    inputs = (Path(__file__), file, reference, metadata_path, download_snapshot / "downloads.json", installed_path)
    manifest = {"input_sha256": {str(path): _sha(path) for path in inputs},
                "outputs": {p.name: _sha(p) for p in output.iterdir()}}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return summary


def main():
    """Execute the RBC source review with annotations, mapping losses and explicit gaps."""
    from notebook_runner import ExecutedNotebook

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--downloads", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    nb = ExecutedNotebook("Quantitative RBC candidate: groups, donors, quantities and replacement decision")
    nb.md("The frozen scope is publisher Table S3 copy counts, human reviewed protein reference and installed "
          "RBC evidence. Protein groups and four donors in whole cells/white ghosts are kept distinct. "
          "Ambiguity and reverse mapping collisions are withheld from single-protein projection. "
          "Citation factors cannot promote an incomparable quantity. Original numeric tables remain external.")
    nb.code("from pathlib import Path", "import sys", f"sys.path.insert(0, {str(ROOT / 'scripts')!r})",
            "from review_rbc_candidate import review",
            f"review(Path({str(args.root)!r}), Path({str(args.downloads)!r}), Path({str(args.out)!r}))")
    nb.write(str(args.out / "review.ipynb"))


if __name__ == "__main__":
    main()
