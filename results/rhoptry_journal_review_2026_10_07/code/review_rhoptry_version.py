"""Retrieve the linked journal supplement and compare the installed preprint source.

The archive URL is taken from the public article's Dataset EV1 link. Source
measurements, mapping losses and installed value reproduction are tested before
any citation or data version is promoted. Runtime measurements are not edited.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
import zipfile

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from starplast import deposits, host, organisms as O, paths  # noqa: E402
from starplast.provenance import audit_mapping  # noqa: E402
from recover_source_files import download  # noqa: E402

PAGE = "https://link.springer.com/article/10.1038/s44318-026-00911-z"
URL = "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs44318-026-00911-z/MediaObjects/44318_2026_911_MOESM2_ESM.zip"
FIELDS = ("Gene", "Combined Rhoptry Score", "greenN|beta", "greenN|wald-fdr")


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def compare_tables(before: pd.DataFrame, after: pd.DataFrame) -> dict:
    """Compare exact gene-level assay values, retaining missingness and changed fields."""
    for table in (before, after):
        if not set(FIELDS).issubset(table.columns) or table.Gene.isna().any() or table.Gene.duplicated().any():
            raise ValueError("Version comparison needs a unique gene-level screen table")
    a, b = [table.set_index("Gene").sort_index() for table in (before, after)]
    report = {"before_rows": len(a), "after_rows": len(b), "gene_sets_equal": a.index.equals(b.index),
              "common_genes": len(a.index.intersection(b.index)), "fields": {}}
    for name in FIELDS[1:]:
        one = pd.to_numeric(a[name], errors="raise")
        two = pd.to_numeric(b[name], errors="raise").reindex(a.index)
        equal = one.eq(two) | (one.isna() & two.isna())
        delta = (one - two).dropna().abs()
        report["fields"][name] = {"identical_cells": int(equal.sum()), "different_cells": int((~equal).sum()),
                                  "before_missing": int(one.isna().sum()), "after_missing": int(b[name].isna().sum()),
                                  "max_abs_difference": float(delta.max()) if not delta.empty else None}
    report["assay_values_identical"] = report["gene_sets_equal"] and all(
        item["different_cells"] == 0 for item in report["fields"].values())
    return report


def review(root: Path, output: Path) -> dict:
    """Keep journal/preprint equivalence, actual symbol mapping losses and reproduction."""
    if output.exists():
        raise ValueError("Use a new journal-version review snapshot")
    output.mkdir(parents=True)
    target = root / "host/k562_rhoptry_screen/journal_2026" / Path(URL).name
    record = {"source_page": PAGE, "requested_url": URL,
              "retrieved_utc": datetime.now(timezone.utc).isoformat()}
    try:
        download(URL, target)
        record.update(status="retrieved_validated_archive", path=str(target), sha256=_sha(target), bytes=target.stat().st_size)
        with zipfile.ZipFile(target) as archive:
            record["archive_members"] = [{"name": item.filename, "bytes": item.file_size} for item in archive.infolist()]
            candidates = []
            for item in archive.infolist():
                if not item.filename.lower().endswith(".xlsx") or item.file_size > 64 << 20:
                    continue
                payload = archive.read(item)
                book = pd.ExcelFile(io.BytesIO(payload))
                for sheet in book.sheet_names:
                    table = pd.read_excel(book, sheet_name=sheet)
                    if set(FIELDS).issubset(table.columns):
                        candidates.append((item.filename, sheet, hashlib.sha256(payload).hexdigest(), table))
            if len(candidates) != 1:
                raise ValueError("No unique journal workbook/sheet contains the required screen fields")
        member, sheet, digest, after = candidates[0]
        before_path = root.joinpath(*deposits.K562)
        before = pd.read_excel(before_path, sheet_name="Rhoptry discharge MAgECK-MLE")
        comparison = compare_tables(before, after)
        comparison.update(journal_member=member, journal_sheet=sheet, journal_member_sha256=digest,
                          preprint_file_sha256=_sha(before_path))
        entries = root / host.UNIPROT_ROOT / host.UNIPROT_IDMAP["human"][1]
        pairs = []
        with gzip.open(entries, "rt") as stream:
            next(stream)
            for line in stream:
                fields = line.rstrip("\n").split("\t")
                if len(fields) > 1 and fields[1]:
                    pairs.append((fields[1], fields[0]))
        mapping = audit_mapping(before.Gene.astype(str).str.strip(), pairs,
                                source_organism=O.HUMAN, target_organism=O.HUMAN,
                                source_namespace="published gene symbol", target_namespace="reviewed UniProt protein",
                                source_version="preprint SHA256:" + _sha(before_path),
                                target_version="reviewed human SHA256:" + _sha(entries), reference=str(entries))
        reproduced = deposits.k562_rhoptry_screen(str(root)).set_index("host_id")
        serialized = pd.read_csv(io.StringIO(reproduced.to_csv(sep="\t", float_format="%.10g")), sep="\t", index_col=0)
        installed = pd.read_parquet(paths.cache_file(O.HOST_TABLES[O.HUMAN])).set_index("host_id")
        checks = {}
        for column in ("rhoptry_discharge_score", "rhoptry_discharge_beta", "rhoptry_discharge_fdr"):
            expected = installed[column].dropna()
            actual = serialized[column].reindex(expected.index)
            checks[column] = {"installed_nonmissing": len(expected), "exact_matches": int(actual.eq(expected).sum()),
                              "missing_installed_ids": int(actual.isna().sum()),
                              "max_abs_difference": float((actual - expected).abs().max())}
        if any(c["exact_matches"] != c["installed_nonmissing"] for c in checks.values()):
            raise ValueError("Historical source/serialization does not reproduce all installed screen values")
        record.update(comparison=comparison, mapping_audit=asdict(mapping), installed_reproduction=checks,
                      decision="retain_measurements; journal citation upgrade eligible after independent metadata identity checks"
                      if comparison["assay_values_identical"] else "journal measurements differ; no automatic source replacement")
        if mapping.ambiguous_entities:
            record["mapping_followup"] = "Legacy symbol mapping picks first reviewed accession; audit before promoting any remapped data"
        (target.parent / "URLS.txt").write_text(URL + "\n")
        (target.parent / "SHA256SUMS.txt").write_text(record["sha256"] + "  " + target.name + "\n")
    except Exception as error:
        record.update(status="review_or_retrieval_failed", error=f"{type(error).__name__}: {error}")
    (output / "review.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    manifest = {"input_sha256": {str(Path(__file__)): _sha(__file__),
                                str(ROOT / "starplast/deposits.py"): _sha(ROOT / "starplast/deposits.py"),
                                str(ROOT / "starplast/host.py"): _sha(ROOT / "starplast/host.py")},
                "outputs": {"review.json": _sha(output / "review.json")}}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return record


def main():
    """Run the controlled journal comparison in an executed annotated notebook."""
    from notebook_runner import ExecutedNotebook

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    nb = ExecutedNotebook("Host rhoptry screen: published journal versus installed preprint")
    nb.md("The source identity was verified through PubMed UpdateOf and the public publisher's Dataset EV1 link. "
          "Compare gene sets, combined score, beta and Wald FDR exactly. Audit all symbol mapping collisions "
          "and reproduce the installed values through the original serialization. No values are promoted by this review.")
    nb.code("from pathlib import Path", "import sys", f"sys.path.insert(0, {str(ROOT / 'scripts')!r})",
            "from review_rhoptry_version import review", f"review(Path({str(args.root)!r}), Path({str(args.out)!r}))")
    nb.write(str(args.out / "review.ipynb"))


if __name__ == "__main__":
    main()
