"""Build a complete, offline review queue from the frozen dataset-selection audit.

This records source gaps and triages diagnostic comparisons without admitting
abstract matches, merging publication versions or replacing measured datasets.
"""
from __future__ import annotations

import argparse
from datetime import date
import hashlib
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from starplast import organisms as O  # noqa: E402


def diagnostic_review(slot: str) -> str:
    """Explain the observed quantity/context traps in this dated diagnostic sweep."""
    notes = {
        O.TOXOPLASMA + "_transcription_tachyzoite":
            "In vivo brain and in vitro tachyzoite transcription are different contexts; incumbent bibliometrics are incomplete. No replacement supported.",
        O.TOXOPLASMA + "_localization_measured":
            "Same source: confidence flags and probability are different quantities. Fully stored confidence flags do not imply complete localization assays.",
        O.FALCIPARUM + "_acetylation":
            "Same source: site count and binary detection are different outputs; filled false flags do not extend the measured proteome.",
        O.FALCIPARUM + "_lactylation_asexual_blood_stage":
            "Same source: site count and binary detection are different outputs; filled false flags do not extend the measured proteome.",
        O.FALCIPARUM + "_fitness_transferred_from_pb":
            "Same source: transferred growth and transfer confidence are different quantities. Confidence is not a fitness replacement.",
        O.FALCIPARUM + "_field_variation_and_resistance_markers":
            "Same source diagnostic preference compares adjusted pN/pS with variant fraction. Cross-source comparisons also require matched population, quantity and normalization.",
    }
    return notes.get(slot, "Diagnostic flag requires source, context, quantity and assay-coverage review; no preference admitted.")


def build_review(output: Path) -> dict:
    """Reconcile all sources and slots into explicit decisions and follow-up gaps."""
    inputs = [output / name for name in ("source_publications.csv", "slot_comparisons.csv", "manifest.json",
                                        "followup/publications.json")]
    inventory_path = ROOT / "results/information_inventory_2026_10_07/inventory.json"
    inputs += [inventory_path, Path(__file__)]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    sources = pd.read_csv(inputs[0], keep_default_na=False)
    slots = pd.read_csv(inputs[1], keep_default_na=False)
    inventory = json.loads(inventory_path.read_text())
    ledger = []
    for source in sources.to_dict("records"):
        rows = [r for r in inventory if r["source_id"] == source["dataset_id"]]
        matching = [r["slot"] for r in slots.to_dict("records")
                    if any(c["source_id"] == source["dataset_id"] for c in json.loads(r["candidate_factors"]))]
        ledger.append({"dataset_id": source["dataset_id"], "organism": source["organism"],
                       "inventory_statuses": json.dumps(sorted({r["status"] for r in rows})),
                       "scopes": json.dumps(sorted({(r["organism"], r["unit"]) for r in rows})),
                       "numeric_slot_views": json.dumps(matching),
                       "publication_status": source["publication_status"],
                       "publication_id": source["publication_id"], "citations_per_year": source["citations_per_year"],
                       "origin_paper_verified": False, "measured_comprehensiveness": "unknown",
                       "literature_best_verified": False, "replacement_admission": "pending",
                       "next_check": "Check origin paper and deposit, matching question/quantity/context, mapped assay denominator, QC and eligible alternatives. Derived/reference resources need version and upstream lineage review."})
    pd.DataFrame(ledger).to_csv(output / "source_review_queue.csv", index=False)
    diagnostics = [{"slot": r["slot"], "diagnostic_source": r["diagnostic_source"],
                    "diagnostic_pattern": r["diagnostic_pattern"],
                    "decision": "no_replacement_supported", "reason": diagnostic_review(r["slot"])}
                   for r in slots.to_dict("records") if r["diagnostic_default_change"]]
    pd.DataFrame(diagnostics).to_csv(output / "diagnostic_reviews.csv", index=False)
    host = slots[slots.unit.eq("host_gene")]
    host_rows = [{"parasite": code, "filled_slot_views": int((arm.candidate_groups > 0).sum()),
                  "total_slot_views": len(arm), "host_gene_space_status": "pending", "inference_pack_status": "pending"}
                 for code, arm in host.groupby("organism")]
    pd.DataFrame(host_rows).to_csv(output / "host_progress.csv", index=False)
    records = json.loads(inputs[3].read_text())
    manifest = json.loads(inputs[2].read_text())
    bibliometrics = []
    for record in records:
        published = record.get("firstPublicationDate")
        age = (date.fromisoformat(manifest["as_of"]) - date.fromisoformat(published)).days if published else None
        count = record.get("citedByCount")
        rate = count / (max(age, manifest["policy"]["minimum_age_days"]) / 365.2425) if age is not None and count is not None else None
        bibliometrics.append({"pmid": record["pmid"], "doi": record.get("doi"), "title": record["title"],
                              "publication_date": published, "citations": count, "citations_per_year": rate,
                              "provider": "Europe PMC", "as_of": manifest["as_of"], "admission": "pending"})
    followup = {"publication_factors": bibliometrics, "decisions": [
        {"source": "host_k562_rhoptry_screen", "decision": "publication_version_review_pending",
         "reason": "A September 2026 journal paper matches the preprint authors, question and headline screen. Supplementary-table equivalence and first-publication/version citation lineage are not verified; retain the preprint measurements and provenance."},
        {"source": "host_erythrocyte_proteome", "decision": "retain_fraction_data; quantitative_complement_review_pending",
         "reason": "The 2017 study has higher annual citation rate and reports 2,650 proteins, absolute copy numbers and validation with labeled standards. The installed 2026 study reports 5,264 proteins and fraction-specific PSM counts. These quantities are not interchangeable; broad protein counts lack a common assay denominator. Admit and validate a separately labeled quantitative source before comparing defaults."},
        {"source": "host_gtex_transcriptome", "decision": "release_lineage_review_pending",
         "reason": "The registry associates v10 data with the 2020 GTEx atlas publication. Preserve release identity and check the originating-publication association; citation rate of the atlas alone does not verify the v10 cell-type columns."},
    ], "promotions": 0}
    (output / "host_followup.json").write_text(json.dumps(followup, indent=2) + "\n")
    if len(ledger) != 162 or len(slots) != 297 or len(host) != 50:
        raise ValueError("Dated review scope does not reconcile")
    if slots.change_proposed.any() or slots.rankable_groups.any():
        raise ValueError("Unreviewed diagnostics became admitted preferences")
    summary = {"source_review_rows": len(ledger), "slot_review_scope": len(slots),
               "diagnostics_triaged": len(diagnostics), "host_slot_views": len(host), "promotions": 0}
    artifact_names = ("source_review_queue.csv", "diagnostic_reviews.csv", "host_progress.csv", "host_followup.json")
    (output / "review_manifest.json").write_text(json.dumps({"input_sha256": hashes, "summary": summary,
        "outputs": {name: hashlib.sha256((output / name).read_bytes()).hexdigest() for name in artifact_names}}, indent=2) + "\n")
    return summary


def main():
    """Write the executed offline review notebook alongside its frozen audit."""
    from notebook_runner import ExecutedNotebook

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    notebook = ExecutedNotebook("Dataset-selection review queue and host progress")
    notebook.md("Review is separate from popularity scoring. Every source remains in the queue. "
                "Six diagnostic flags are triaged for context and quantity traps; no new source is admitted. "
                "Host publication follow-ups use the separately executed exact-identity lookup.")
    notebook.code("from pathlib import Path", "import sys", f"sys.path.insert(0, {str(ROOT / 'scripts')!r})",
                  "from review_dataset_selection import build_review",
                  f"build_review(Path({str(args.out.resolve())!r}))")
    notebook.write(str(args.out / "review.ipynb"))


if __name__ == "__main__":
    main()
