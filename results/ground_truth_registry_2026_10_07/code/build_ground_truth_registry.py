"""Freeze current truth candidates and gaps in an executed offline census notebook."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from starplast import datasets as D, organisms as O, paths, scorecard as C, strategies as S  # noqa: E402
from starplast.ground_truth import BenchmarkEntry, StrategyBenchmark, cohort_digest, eligibility_mask, write_registry  # noqa: E402
from starplast.provenance import SourceFile, read_traces  # noqa: E402

# These declarations classify the output, not the underlying instrument. The
# basis is the current registry and loader; paper association stays unresolved.
SOURCE_GRADES = {
    "lopit_tgon": ("prediction", "localization.lopit_labels: native MAP/MCMC assignments from measured profiles"),
    "pf_spatial_proteome": ("prediction", "registry: SVM classifier outputs from hyperLOPIT profiles"),
    "pf_berghei_transfer": ("orthology_transfer", "registry and deposit: berghei mutant phenotypes carried to falciparum"),
    "plasmodb_pf3d7_exportpred": ("prediction", "registry: ExportPred sequence model evaluated at multiple thresholds"),
    "stage_enriched": ("derived_quantity", "registry derived_from: standardized expression argmax"),
    "xue_singlecell": ("derived_quantity", "registry: inferred transcriptional phase/pseudotime from measured RNA profiles"),
    "splitcas9_imaging_screen": ("direct_experiment", "registry and imaging loader: screened gene perturbation phenotypes"),
    "host_k562_rhoptry_screen": ("direct_experiment", "journal-equivalence audit: measured K562 CRISPR screen; gene-symbol projection unresolved for 39 proteins"),
    "host_bmdm_baseline": ("direct_experiment", "provenance pilot: mouse baseline RNA abundance with exact raw-to-installed reproduction"),
}


def classify(source, column):
    """Classify explicit legacy output declarations without guessing from column names."""
    if source is None:
        return "unresolved", "No organism-qualified registry source for this target"
    if source.key == "lopit_tgon" and column == "compartment_best":
        return "orthology_transfer", "localization.lopit_labels: native classifier call with cross-species fallback; mixed row provenance"
    if source.derived_from:
        return "derived_quantity", "Registry declares derived_from=" + repr(source.derived_from)
    if source.key in SOURCE_GRADES:
        return SOURCE_GRADES[source.key]
    if source.kind in {"CRISPR_screen", "piggyBac_screen"} and "transferred" not in source.provides.lower():
        return "direct_experiment", "Registry declares a perturbation assay; actual source association and field semantics require review"
    if source.kind in {"curated", "curation"}:
        return "curation", "Registry declares curated annotations; review scope and primary evidence before admission"
    return "unresolved", "Source output grade not reviewed; assay kind alone cannot classify a processed target"


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build(output: Path) -> dict:
    """Freeze eligible observations and strategy references without changing runtime data."""
    if output.exists():
        raise ValueError("Use a new snapshot directory")
    provenance_path = ROOT / "results/measurement_provenance_recovered_2026_10_07/traces.json"
    traces = read_traces(provenance_path)
    by_address = {(t.output_organism, t.column): t for t in traces}
    inputs = [Path(__file__), provenance_path, *(ROOT / "starplast" / name for name in
              ("ground_truth.py", "datasets.py", "strategies.py", "track_record.py", "localization.py", "search.py",
               "scorecard.py", "organisms.py", "strategy_catalog.py", "strategy_graph.py", "strategy_learning.py"))]
    tables = {}
    for organism in (O.TOXOPLASMA, O.FALCIPARUM, O.HUMAN, O.MOUSE):
        path = Path(paths.cache_file(O.HOST_TABLES[organism])) if organism in O.HOST_TABLES else Path(O.nodes_path(organism))
        inputs.append(path)
        tables[organism] = pd.read_parquet(path)
    hashes = {str(path): _sha(path) for path in inputs}
    output.mkdir(parents=True)
    entries, target_summary = [], []
    for organism, frame in tables.items():
        host = organism in O.HOST_TABLES
        unit, id_column = ("protein", "host_id") if host else ("gene", "gene_id")
        ids = frame[id_column]
        if ids.isna().any() or ids.astype(str).duplicated().any():
            raise ValueError("Truth freeze requires a unique explicitly identified universe")
        ids = ids.astype(str).tolist()
        if host:
            wanted = {column for source in D.REGISTRY if source.organism == organism for column in source.columns if column in frame}
            kinds = {column: "numeric" if pd.api.types.is_numeric_dtype(frame[column]) and not
                     pd.api.types.is_bool_dtype(frame[column]) else "categorical" for column in wanted}
        else:
            ctx = S.Context(frame, graph={}, organism=organism)
            wanted = set(ctx.categorical_columns()) | set(ctx.numeric_targets()) | set(O.get(organism).targets) | set(O.get(organism).numbers)
            kinds = {column: "numeric" if pd.api.types.is_numeric_dtype(frame[column]) and not
                     pd.api.types.is_bool_dtype(frame[column]) else "categorical" for column in wanted if column in frame}
        freeze = frame[[id_column, *sorted(kinds)]].copy()
        freeze_path = output / (organism + "_truth.parquet")
        freeze.to_parquet(freeze_path, index=False)
        universe_path = output / (organism + "_universe.json")
        universe_path.write_text(json.dumps(ids) + "\n")
        truth_file = SourceFile.inspect(freeze_path, "installed_cache", association="frozen_installed_targets_not_original_assay")
        universe_file = SourceFile.inspect(universe_path, "mapping_reference", association="exact_installed_entity_order")
        for column, kind in sorted(kinds.items()):
            source = D.provenance(column, organism=organism)
            trace = by_address.get((organism, column))
            grade, basis = classify(source, column)
            values = pd.Series(frame[column].to_numpy(), index=ids)
            mask = eligibility_mask(values, kind)
            gaps = ["source_assayed_population_unresolved", "independent_validation_protocol_pending_64.05",
                    "origin_paper_and_processed_output_association_requires_review"]
            if trace:
                gaps.extend(trace.gaps)
            if grade == "unresolved":
                gaps.append("output_evidence_grade_unresolved")
            if grade in {"prediction", "orthology_transfer", "derived_quantity"}:
                gaps.append("requires_independent_experimental_truth_for_biological_accuracy")
            if host:
                gaps.append("protein_reference_is_not_an_admitted_host_gene_space")
            if source and source.key == "host_k562_rhoptry_screen":
                gaps.append("39_arbitrary_symbol_projections_pending_66.04; stored_cohort_not_admitted")
            contexts = tuple(trace.context) if trace and trace.context else ("legacy_context_unresolved",)
            quantity = trace.quantity_unit if trace else "unresolved"
            tasks = (C.T_VALUES,) if kind == "numeric" else (C.T_LABEL, C.T_RANK, C.T_SET, C.T_CLUSTER)
            for task in tasks:
                key = organism + ":" + column + ":" + task.replace(" ", "_") + ":v1"
                entries.append(BenchmarkEntry(
                    key, organism, column, task, unit, grade, basis, (source.key,) if source else (),
                    truth_file, universe_file, np.packbits(mask.to_numpy()).tobytes().hex(),
                    cohort_digest(ids, mask.tolist()), len(ids), int(mask.sum()), None,
                    contexts, quantity,
                    ("Missing/unknown/nonfinite values excluded; no imputation", "Ambiguous historical mappings not yet fully audited; candidate only",
                     "Target family and source siblings must be withheld; split protocol pending"),
                    "Missing is unknown. Zero/False is a stored result, not certified biological absence. "
                    "No negative class admitted until assay-specific review.", gaps=tuple(sorted(set(gaps)))))
            target_summary.append({"organism": organism, "target": column, "kind": kind, "source_id": source.key if source else "",
                                   "evidence_grade": grade, "stored_population": len(ids), "eligible_population": int(mask.sum()),
                                   "measured_population": None, "status": "candidate"})
    refs = []
    relationship_keys = {"link_prediction", "unwritten_links", "network_training"}
    for organism in tables:
        for strategy in S.catalog():
            candidates = tuple(e.benchmark_id for e in entries if e.organism == organism and e.task == strategy.task)
            gaps = ["strategy_specific_protocol_and_source_exclusions_pending", "no_independent_biological_benchmark_admitted"]
            if strategy.key in relationship_keys:
                candidates = ()
                gaps.append("pair_truth_and_assayed_negative_universe_unresolved; gene_labels_cannot_validate_relationships")
            if strategy.task == C.T_REPL:
                candidates = ()
                gaps.append("disjoint_discovery_validation_findings_not_registered; split_half_is_not_external_replication")
            if organism in O.HOST_TABLES:
                gaps.append("host_gene_context_not_yet_admitted; protein_candidates_only")
            refs.append(StrategyBenchmark(organism, strategy.key, strategy.task, candidates, tuple(gaps)))
    write_registry(output / "registry.json", entries, refs)
    pd.DataFrame(target_summary).to_csv(output / "targets.csv", index=False)
    pd.DataFrame([asdict(row) for row in refs]).to_json(output / "strategy_gaps.json", orient="records", indent=2)
    for path, digest in hashes.items():
        if _sha(path) != digest:
            raise ValueError("Source changed during truth census: " + path)
    summary = {"target_addresses": len(target_summary), "benchmark_candidates": len(entries), "admitted_benchmarks": 0,
               "strategies_per_species": len(S.catalog()), "strategy_species_addresses": len(refs),
               "grades": pd.DataFrame(target_summary).evidence_grade.value_counts().to_dict(),
               "measured_population": "unresolved; installed nonmissing observations are not assay denominators",
               "runtime_changes": "none"}
    outputs = {path.name: _sha(path) for path in output.iterdir() if path.is_file()}
    (output / "manifest.json").write_text(json.dumps({"schema_version": 1, "base_commit": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(), "input_sha256": hashes,
        "outputs": outputs, "summary": summary,
        "software": {"python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__}}, indent=2) + "\n")
    return summary


def main():
    """Execute and retain the frozen pilot census with its interpretation limits."""
    from notebook_runner import ExecutedNotebook

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    notebook = ExecutedNotebook("Ground-truth candidate registry and strategy validation gaps")
    notebook.md("Frozen pilot: installed parasite targets and host protein fields. No fitting, acquisition or numerical changes. "
                "Eligibility is observational availability, not independent truth admission. Synthetic controls cannot close biological gaps.")
    notebook.code("from pathlib import Path", "import sys", f"sys.path.insert(0, {str(ROOT / 'scripts')!r})",
                  "from build_ground_truth_registry import build", f"build(Path({str(args.out.resolve())!r}))")
    notebook.md("Every strategy/species task has candidate references or a specific gap. Predictor/transfer/derived truth remains "
                "separate from direct assays. All measured assay denominators are unresolved; zero/False does not certify a negative. "
                "Every entry is a candidate and cannot substantiate an independent biological accuracy claim. "
                "Source-specific review, split protocols and the host mapping correction remain required.")
    notebook.write(str(args.out / "census.ipynb"))


if __name__ == "__main__":
    main()
