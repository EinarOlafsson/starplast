"""Execute frozen biological candidate splits and deliberate leakage guard probes."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
from itertools import combinations
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from starplast import organisms as O, strategies as S, scorecard as C  # noqa: E402
from starplast.ground_truth import read_registry  # noqa: E402
from starplast.splits import cold_node_pairs, known_node_pairs, make_exclusions, make_split, read_split, write_split  # noqa: E402


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _refusal(call):
    try:
        call()
    except ValueError as error:
        return str(error)
    raise AssertionError("Deliberate leakage was not refused")


def build(output):
    """Freeze two real candidate cohorts; exercise guards without model fitting."""
    output = Path(output)
    if output.exists():
        raise ValueError("Use a new split-pilot directory")
    registry = ROOT / "results/ground_truth_registry_2026_10_07_v2/registry.json"
    entries, _refs = read_registry(registry)
    inputs = [registry, Path(__file__), *(ROOT / "starplast" / name for name in
              ("splits.py", "ground_truth.py", "strategies.py", "search.py", "datasets.py", "organisms.py", "embedding.py"))]
    frames = {code: pd.read_parquet(O.nodes_path(code)) for code in (O.TOXOPLASMA, O.FALCIPARUM)}
    inputs.extend(Path(O.nodes_path(code)) for code in frames)
    hashes = {str(path): _sha(path) for path in inputs}
    output.mkdir(parents=True)
    report = []
    for organism, target in ((O.TOXOPLASMA, "fit_invitro_hff"), (O.FALCIPARUM, "piggybac_mis")):
        entry = next(e for e in entries if e.organism == organism and e.target == target and e.task == C.T_VALUES)
        ids = json.loads(Path(entry.universe_file.path).read_text())
        mask = entry.mask(ids)
        frame = frames[organism].set_index("gene_id", drop=False)
        if frame.index.tolist() != ids:
            raise ValueError("Candidate and installed entity universes no longer match")
        ctx = S.Context(frame.reset_index(drop=True), graph={}, organism=organism)
        groups = ctx.groups()[mask.to_numpy()].tolist()
        eligible = mask.index[mask].tolist()
        split = make_split(eligible, groups, organism=organism, benchmark_id=entry.benchmark_id, seed=17,
                           feature_access="transductive")
        training = frame.loc[list(split.entities("train"))].reset_index(drop=True)
        exclusions = make_exclusions(training, target, benchmark_id=entry.benchmark_id, split=split)
        probes = {
            "scaling_test_access": _refusal(lambda: split.guard_fit("scaling", split.entities("test"))),
            "tuning_calibration_access": _refusal(lambda: split.guard_fit("setting_selection", split.entities("calibration"))),
            "calibration_test_access": _refusal(lambda: split.guard_fit("calibration", split.entities("test"))),
            "target_input": _refusal(lambda: exclusions.guard_inputs(columns=[target])),
            "test_labels_for_source_selection": _refusal(lambda: make_exclusions(frame.loc[list(split.entities("test"))], target,
                                            benchmark_id=entry.benchmark_id, split=split)),
        }
        for stage in ("imputation", "scaling", "representation", "feature_selection", "model_fit"):
            split.guard_fit(stage, split.entities("train"))
        split.guard_fit("setting_selection", split.entities("tune"))
        split.guard_fit("calibration", split.entities("calibration"))
        filename = organism + "_" + target + ".json"
        write_split(output / filename, split, exclusions)
        assert read_split(output / filename) == (split, exclusions)
        report.append({"organism": organism, "target": target, "benchmark_id": entry.benchmark_id,
                       "cohort_sha256": entry.cohort_sha256, "split_identity": split.identity,
                       "stored_population": entry.stored_population, "eligible_population": len(eligible),
                       "assayed_population": entry.measured_population, "candidate_status": entry.status,
                       "partitions": {role: len(split.entities(role)) for role in ("train", "tune", "calibration", "test")},
                       "protected_groups": len(set(groups)), "group_policy": "existing Context.groups; missing orthogroups are singleton entities, unknown homology remains a gap",
                       "excluded_columns": len(exclusions.columns), "excluded_layers": list(exclusions.layers),
                       "deliberate_leakage_refusals": probes, "accuracy": "not_measured; no model fitting"})
    # Pair fixtures certify endpoint semantics only, never biological truth or negatives.
    nodes = [f"fixture_node_{i}" for i in range(16)]
    pairs = list(combinations(nodes, 2))
    known = known_node_pairs(pairs, organism=O.TOXOPLASMA, benchmark_id="synthetic_pair_regime:v1", seed=17)
    node_split = make_split(nodes, nodes, organism=O.TOXOPLASMA, benchmark_id="synthetic_cold_nodes:v1", group_kind="entity", seed=17)
    cold = cold_node_pairs(pairs, node_split)
    (output / "pair_regime_fixtures.json").write_text(json.dumps({"grade": "synthetic_control",
        "biological_admission": "not_ground_truth", "known_node": asdict(known), "cold_node": asdict(cold)}, indent=2) + "\n")
    (output / "pilot_report.json").write_text(json.dumps(report, indent=2) + "\n")
    for path, digest in hashes.items():
        if _sha(path) != digest:
            raise ValueError("Pilot input changed: " + path)
    summary = {"candidate_cohorts": len(report), "fitted_models": 0, "deliberate_leakage_probes_refused": sum(len(row["deliberate_leakage_refusals"]) for row in report),
               "runtime_changes": "none", "real_pair_truth": "unresolved; synthetic regime fixture only"}
    (output / "manifest.json").write_text(json.dumps({"input_sha256": hashes, "summary": summary,
        "outputs": {path.name: _sha(path) for path in output.iterdir() if path.is_file()}}, indent=2) + "\n")
    return report


def main():
    """Keep executed split generation and its scientific limits as an annotated notebook."""
    from notebook_runner import ExecutedNotebook

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    notebook = ExecutedNotebook("Nested hold-outs, training-only source closure and leakage refusals")
    notebook.md("Frozen pilot: two existing parasite numeric candidates, protected homology groups, four disjoint roles. "
                "No fitting, source acquisition or runtime changes. Exclusion selection reads training labels only. "
                "Transductive joint graph access is declared; it does not permit fitting on held-out labels.")
    notebook.code("from pathlib import Path", "import sys", f"sys.path.insert(0, {str(ROOT / 'scripts')!r})",
                  "from freeze_split_pilots import build", f"build(Path({str(args.out.resolve())!r}))")
    notebook.md("Known-node and cold-node pair fixtures test endpoint guards only. Dataset/context/time groups have software checks; "
                "external biological truth, pair negatives and unknown homology still need source-specific review. "
                "No new biological accuracy or capacity estimate is produced; candidate assay populations remain unresolved. "
                "Future benchmark adapters must call these guards for actual inputs before fitting, tuning and calibration.")
    notebook.write(str(args.out / "pilots.ipynb"))


if __name__ == "__main__":
    main()
