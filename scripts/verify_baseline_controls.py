"""Execute annotated synthetic checks of shared baseline and null-control contracts."""
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

from starplast import organisms as O  # noqa: E402
from starplast.baselines import (RECIPES, assert_planted_recovery, label_baselines, matched_network_control,
                                positive_only_retrieval, random_ranking, value_baselines)  # noqa: E402
from starplast.splits import make_split, write_split  # noqa: E402


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build(output):
    """Keep matched synthetic cohorts, planted failures and unavailable nulls reviewable."""
    output = Path(output)
    if output.exists():
        raise ValueError("Use a new baseline-control directory")
    inputs = [Path(__file__), *(ROOT / "starplast" / name for name in
                              ("baselines.py", "splits.py", "ground_truth.py", "scorecard.py", "organisms.py"))]
    hashes = {str(path): _sha(path) for path in inputs}
    output.mkdir(parents=True)
    ids = [f"synthetic_gene_{i}" for i in range(40)]
    split = make_split(ids, ids, organism=O.TOXOPLASMA, benchmark_id="synthetic_baselines:v1", group_kind="entity", seed=4)
    train = split.entities("train")
    label_truth = pd.Series(["common"] * (len(train) - 2) + ["rare"] * 2, index=train)
    numeric_truth = pd.Series([float(i) for i in range(len(train))], index=train)
    labels = label_baselines(label_truth, split, seed=5)
    values = value_baselines(numeric_truth, split)
    assert labels.index.equals(values.index)
    assert tuple(labels.index) == split.entities("test")
    assert labels.attrs["evaluation_n"] == values.attrs["evaluation_n"]
    assert values["mean"].eq(numeric_truth.mean()).all()
    write_split(output / "synthetic_split.json", split)
    labels.to_csv(output / "label_baselines.csv")
    values.to_csv(output / "value_baselines.csv")
    label_truth.rename("label").to_csv(output / "training_labels.csv")
    numeric_truth.rename("value").to_csv(output / "training_values.csv")
    candidates = [f"synthetic_candidate_{i}" for i in range(100)]
    positives = candidates[:10]
    random = random_ranking(candidates, seed=7, query_seeds=["explicit_seed_outside_candidate_universe"])
    random.to_csv(output / "random_ranking.csv", index=False)
    planted = assert_planted_recovery(candidates, positives, depth=10, minimum_recall=1.)
    failures = {}
    for name, ranking in (("reversed_scores", candidates[::-1]), ("null_ranking", random.entity.tolist())):
        try:
            assert_planted_recovery(ranking, positives, depth=10, minimum_recall=.9)
        except AssertionError as error:
            failures[name] = str(error)
        else:
            raise AssertionError("Synthetic harness failed to reject " + name)
    positive_only = positive_only_retrieval(random.entity.tolist(), positives, depth=10)
    assert positive_only["precision"] is positive_only["auroc"] is positive_only["auprc"] is None
    edges = [(f"synthetic_node_{i}", f"synthetic_node_{(i + 1) % 20}") for i in range(20)]
    strata = {f"synthetic_node_{i}": "detection_high" if i % 2 else "detection_low" for i in range(20)}
    network = matched_network_control(edges, strata, swaps=20, seed=8)
    assert network["successful_swaps"] == 20 and network["degree_preserved"] and network["detection_contacts_preserved"]
    complete_nodes = list("abcde")
    unavailable = matched_network_control(list(combinations(complete_nodes, 2)), dict.fromkeys(complete_nodes, "detected"),
                                         swaps=5, max_attempts=30)
    assert unavailable["status"] == "unavailable_no_valid_swaps"
    summary = {"truth_grade": "synthetic_control", "biological_admission": "not_ground_truth", "task_recipes": len(RECIPES),
               "labels": labels.attrs, "values": values.attrs, "ranking": random.attrs, "planted": planted,
               "harness_failure_refusals": failures, "positive_only": positive_only,
               "network": network, "unavailable_network": unavailable,
               "adapter_status": "shared definitions and helpers; task-specific strategy execution remains pending"}
    (output / "recipes.json").write_text(json.dumps({task: asdict(recipe) for task, recipe in RECIPES.items()}, indent=2) + "\n")
    (output / "report.json").write_text(json.dumps(summary, indent=2) + "\n")
    for path, digest in hashes.items():
        if _sha(path) != digest:
            raise ValueError("Baseline input changed: " + path)
    (output / "manifest.json").write_text(json.dumps({"input_sha256": hashes, "truth_grade": "synthetic_control",
        "software": {"python": sys.version.split()[0], "pandas": pd.__version__},
        "outputs": {path.name: _sha(path) for path in output.iterdir() if path.is_file()}}, indent=2) + "\n")
    return {"task_recipes": len(RECIPES), "evaluation_population": len(labels), "matched_cohort": True,
            "planted_recall": planted["known_positive_recall"], "failed_harnesses_refused": len(failures),
            "network_swaps": network["successful_swaps"], "network_mixing": network["mixing"],
            "no_valid_swaps": unavailable["status"], "biological_accuracy": "not_measured"}


def main():
    """Execute software controls and retain their limits in a notebook."""
    from notebook_runner import ExecutedNotebook

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    notebook = ExecutedNotebook("Shared baseline contracts and synthetic positive/null controls")
    notebook.md("Frozen scope: synthetic fixtures only. Training labels/values fit baselines; evaluation entities match exactly. "
                "No biological benchmark admission, source acquisition or runtime strategy change.")
    notebook.code("from pathlib import Path", "import sys", f"sys.path.insert(0, {str(ROOT / 'scripts')!r})",
                  "from verify_baseline_controls import build", f"build(Path({str(args.out.resolve())!r}))")
    notebook.md("Planted recovery passes, while reversed and null rankings fail the harness. Known-positive recall is not full "
                "biological recall; precision/AUROC/AUPRC remain unidentified without verified negatives. Network swaps preserve "
                "degrees and per-node detection-stratum contacts, but mixing remains unestablished. A complete graph refuses a "
                "usable rewiring null. These controls certify software behavior only; subsequent adapters must execute recipes "
                "on the same admitted cohort and split as their strategy.")
    notebook.write(str(args.out / "controls.ipynb"))


if __name__ == "__main__":
    main()
