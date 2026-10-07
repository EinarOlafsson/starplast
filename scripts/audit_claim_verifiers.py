#!/usr/bin/env python3
"""Measure literature verifiers against the shipped claims without changing published recipes.

    python scripts/audit_claim_verifiers.py --out results/claim_verifiers_2026_10_07

An executed notebook records the analysis, input hashes and every refusal. Abstract and full-text
graphs are tested separately, with positive attention residuals; raw counts are retained only as a
control. Each run uses the track record's orthogroup folds and reports incremental coverage over
the verifiers already shipped. Passing this screen is a candidate for further validation, not a
published independent verification or a new biological measurement.
"""
from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from starplast import claims as C  # noqa: E402
from starplast import organisms as O  # noqa: E402
from starplast import strategies as S  # noqa: E402
from starplast import track_record as T  # noqa: E402


def literature_context(ctx, target: str, layer: str, corrected: bool = True):
    """An isolated graph using positive attention residuals, or an explicit refusal reason.

    A banned or absent layer is never replaced with a measurement graph. Missing residuals are
    refused for corrected runs instead of falling back to raw popularity. The source context,
    its graph arrays and its caches remain untouched.
    """
    if layer not in S.LITERATURE_LAYERS:
        raise ValueError(f"not a literature layer: {layer}")
    if layer in ctx.banned_layers(target):
        return None, "banned layer"
    if layer not in ctx.layers():
        return None, "layer unavailable"
    graph = dict(ctx.graph)
    weight = f"{layer}__{'r' if corrected else 'w'}"
    if weight not in graph:
        return None, "attention residuals unavailable" if corrected else "raw weights unavailable"
    weights = np.asarray(graph[weight], dtype=float)
    if not np.isfinite(weights).all():
        return None, "nonfinite edge weights"
    graph[f"{layer}__w"] = np.maximum(weights, 0.0)
    if not graph[f"{layer}__w"].any():
        return None, "no positive edges"
    return S.Context(ctx.nodes, graph=graph, organism=ctx.organism, seed=ctx.seed,
                     log=ctx.log, should_stop=ctx.should_stop), "measured"


def validate_reference(reference: pd.DataFrame, measured: pd.DataFrame, generator: str):
    """Refuse comparisons whose gene identities, truth, folds or seeds differ from the reference.

    Every labelled gene, including abstentions, must occur once in each record. Positional gene
    numbers alone are insufficient: a rebuilt or reordered table must not silently join two genes.
    """
    base = reference[(reference["mode"] == "together") &
                     (reference["strategy"].astype(str) == generator)]
    keys = ["organism", "target", "gene", "gene_id", "truth", "fold", "seed"]
    if base.empty or measured.empty or base["gene"].duplicated().any() or measured["gene"].duplicated().any():
        raise ValueError("comparison needs one held-out row per gene in both records")
    left = base.sort_values("gene")[keys].astype(str).reset_index(drop=True)
    right = measured.sort_values("gene")[keys].astype(str).reset_index(drop=True)
    if not left.equals(right):
        raise ValueError("candidate and reference differ in identities, labels, folds or seeds")


def coverage(claims: pd.DataFrame, calls: pd.DataFrame) -> dict:
    """Coverage of all unknown genes and of existing claims, with new reach counted separately.

    Outside-range claims are reported separately and never counted as newly testable in-range
    claims. Coverage means a verifier spoke; it is not evidence that its answer was correct.
    """
    reached = set(calls.loc[calls["prediction"].notna(), "gene_id"].astype(str))
    claim_ids = claims["gene_id"].astype(str)
    newly = (claims["status"].astype(str) == "untested") & claim_ids.isin(reached)
    tested = claims["status"].astype(str) == "tested"
    outside = claims["status"].astype(str) == "outside tested range"
    return {
        "unknown_genes": len(calls), "unknown_reached": len(reached),
        "existing_claims": len(claims), "claims_reached": int(claim_ids.isin(reached).sum()),
        "already_tested": int(tested.sum()), "new_in_range": int(newly.sum()),
        "outside_range_reached": int((outside & claim_ids.isin(reached)).sum()),
        "tested_with_candidate": int(tested.sum() + newly.sum()),
    }


def audit_candidate(ctx, target: str, layer: str, corrected: bool, reference: pd.DataFrame,
                    recipe: pd.Series, claims: pd.DataFrame) -> tuple:
    """One candidate's held-out dependence, agreement precision, calibration and additional reach.

    The generator and existing verifiers are fixed to the shipped recipe. Candidate selection is
    exploratory; raw controls and dependent candidates are never eligible for promotion. Even a
    passing candidate needs a fresh recipe build and broader validation before publication.
    """
    generator = str(recipe["generator"])
    name = f"literature_{layer}_{'residual' if corrected else 'raw'}"
    row = {"organism": ctx.organism, "target": target, "generator": generator,
           "verifier": name, "layer": layer, "corrected": corrected,
           "independent": False, "eligible": False}
    candidate, reason = literature_context(ctx, target, layer, corrected)
    row["status"] = reason
    if candidate is None:
        return row, pd.DataFrame(), pd.DataFrame()
    settings = {"layer": layer, "restart": 0.5}
    measured = T.evaluate(candidate, "layer_propagation", target, settings=settings)
    if measured.empty:
        row["status"] = "no held-out predictions"
        return row, measured, pd.DataFrame()
    validate_reference(reference, measured, generator)
    measured["strategy"] = name
    measured["setting_key"] = measured["setting_key"].astype(str) + (
        ", weights=positive_attention_residual" if corrected else ", weights=raw_control")
    combined = pd.concat([reference, measured], ignore_index=True)
    row.update(C.shared_mistakes(combined, generator, name))
    row["independent"] = bool(np.isfinite(row["ratio_high"]) and
                              row["ratio_high"] < C.INDEPENDENT_BELOW)
    rates = C.verification(combined, generator, name)
    for outcome in ("base", "agrees", "disagrees", "silent"):
        row.update({f"{outcome}_{key}": value for key, value in rates[outcome].items()})
    answered = measured[~measured["abstained"].astype(bool)]
    row["verifier_answered"] = len(answered)
    row["verifier_accuracy"] = float(answered["correct"].astype(float).mean())
    truth = ctx.truth(target)
    kept = truth.where(truth.map(truth.value_counts()) >= S.MIN_CLASS)
    query = np.flatnonzero(truth.isna().to_numpy())
    pred, support = T._predict(candidate, S.get("layer_propagation"), kept, query, settings, target)
    calls = pd.DataFrame({"organism": ctx.organism, "target": target, "verifier": name,
                          "gene_id": np.asarray(ctx.gene_ids)[query],
                          "prediction": pd.Series(pred).iloc[query].to_numpy(),
                          "support": pd.Series(support).iloc[query].to_numpy()})
    row.update(coverage(claims, calls))
    if row["independent"] and corrected:
        old = tuple(v.strip() for v in str(recipe["verifiers"]).split(",") if v.strip())
        model = C.verified_model(combined, generator, old + (name,))
        row["combined_calibration_error"] = model.calibration_error
        row["combined_calibrated"] = model.calibrated and model.certainty.calibrated
        if hasattr(model, "_oof"):
            confident = model._oof >= C.CONFIDENT
            row["combined_confident"] = int(confident.sum())
            row["combined_confident_precision"] = (float(model._y[confident].mean())
                                                   if confident.any() else np.nan)
            row["eligible"] = bool(row["combined_calibrated"] and row["new_in_range"] > 0 and
                                   row["combined_confident"] >= 30 and
                                   row["combined_confident_precision"] >= C.CONFIDENT)
    return row, measured, calls


def audit() -> tuple:
    """Audit both literature layers on every shipped claim recipe, keeping refusals in the report."""
    rows, records, calls = [], [], []
    for organism in O.codes(available=True):
        ctx = S.Context.shipped(organism)
        reference = T.shipped(organism)
        existing = C.shipped(organism)
        for _, recipe in C.recipes(organism).iterrows():
            target = str(recipe["target"])
            led = reference[reference["target"].astype(str) == target]
            claims = existing[existing["target"].astype(str) == target]
            for layer in S.LITERATURE_LAYERS:
                for corrected in (True, False):
                    row, record, unknown = audit_candidate(ctx, target, layer, corrected,
                                                           led, recipe, claims)
                    rows.append(row)
                    if not record.empty:
                        records.append(record)
                    if not unknown.empty:
                        calls.append(unknown)
                    print(f"{organism}/{target}/{row['verifier']}: {row['status']}; "
                          f"ratio high={row.get('ratio_high', np.nan):.3f}; "
                          f"new in range={row.get('new_in_range', 0)}; eligible={row['eligible']}",
                          flush=True)
    return (pd.DataFrame(rows), pd.concat(records, ignore_index=True) if records else pd.DataFrame(),
            pd.concat(calls, ignore_index=True) if calls else pd.DataFrame())


def fingerprints() -> dict:
    """SHA-256 identities of the code, node tables, graphs, record and claims read by the audit."""
    paths = list((ROOT / "starplast").glob("*.py")) + [Path(__file__)]
    paths += [Path(C._data(name)) for name in ("track_record.parquet", "claims.parquet",
                                              "claim_recipes.parquet")]
    for code in O.codes(available=True):
        paths.extend([Path(O.nodes_path(code)), Path(O.graph_path(code))])
    hashes = {}
    for path in paths:
        digest = hashlib.sha256()
        with path.open("rb") as source:
            for chunk in iter(lambda: source.read(1024 * 1024), b""):
                digest.update(chunk)
        hashes[str(path.relative_to(ROOT) if path.is_relative_to(ROOT) else path)] = digest.hexdigest()
    return hashes


def main(argv=None) -> int:
    """Write the measured tables, provenance and an executed notebook in a fresh output directory."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=False)
    before = fingerprints()
    nb = ExecutedNotebook("Wider independent verification: the literature graphs")
    nb.md("The published generator and its existing verifiers stay fixed. Each candidate uses the "
          "same held-out genes and orthogroup folds as the published track record. Positive "
          "attention-corrected residuals are the candidates; raw co-mentions are controls only.",
          "Independence is measured by the shared-mistake ratio's upper 95% bootstrap bound, "
          f"which must be below {C.INDEPENDENT_BELOW}. Passing is an exploratory result, not "
          "prospective validation: publications may themselves discuss the target experiments.",
          "Coverage counts genes the verifier answers about, including disagreements. Outside-range "
          "genes stay separate; this audit cannot calibrate their certainty. No claims or recipe "
          "files are updated.")
    nb.code("from pathlib import Path", "import sys, json, pandas as pd",
            "root = next(p for p in [Path.cwd(), *Path.cwd().parents]",
            "            if (p / 'pyproject.toml').exists() and (p / 'starplast').is_dir())",
            "sys.path.insert(0, str(root / 'scripts'))",
            "from audit_claim_verifiers import audit, fingerprints",
            f"out = Path({str(args.out.resolve())!r})", "before = fingerprints()")
    nb.code("summary, heldout, calls = audit()",
            "summary.to_csv(out / 'summary.csv', index=False)",
            "heldout.to_parquet(out / 'heldout.parquet', index=False)",
            "calls.to_parquet(out / 'unknown_calls.parquet', index=False)",
            "summary[['organism', 'target', 'verifier', 'status', 'independent', 'eligible']]")
    nb.md("## Corrected candidates and extra coverage",
          "The complete report is summary.csv. New in-range counts currently untested claims "
          "that a candidate can reach; it does not promote them to tested claims. Calibration "
          "is measured only when a corrected candidate meets the dependence threshold.")
    nb.code("summary[summary.corrected & summary.status.eq('measured')][",
            "    ['organism', 'target', 'verifier', 'ratio', 'ratio_low', 'ratio_high',",
            "     'agrees_n', 'agrees_rate', 'disagrees_rate', 'unknown_genes', 'unknown_reached',",
            "     'new_in_range', 'outside_range_reached', 'eligible']].round(3)")
    after = fingerprints()
    if before != after:
        raise RuntimeError("audit inputs changed while running; do not use these results")
    manifest = {"inputs": before, "independent_below": C.INDEPENDENT_BELOW,
                "folds": T.FOLDS, "fold_seed": 0, "restart": 0.5, "bootstrap": C.BOOT,
                "python": sys.version, "packages": {name: version(name) for name in
                                                     ("numpy", "pandas", "scipy", "scikit-learn")}}
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    nb.md("## Exact inputs", "Source and data hashes below are checked again after the analysis.")
    nb.code("assert fingerprints() == before, 'inputs changed during the audit'", "before")
    nb.write(str(args.out / "audit.ipynb"))
    print(nb.ns["summary"][["organism", "target", "verifier", "independent", "eligible"]]
          .to_string(index=False))
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
