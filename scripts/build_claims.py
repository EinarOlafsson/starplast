#!/usr/bin/env python3
"""Build the shipped claims: for every recorded label, what the best testable recipe claims about genes
that have no label, with each claim's certainty, independent verdicts and status.

Writes `starplast/data/claims.parquet` and `starplast/data/claim_recipes.parquet`. Every number in them
comes from the shipped track record (`track_record.parquet`), so build that first. A label whose best
recipe is not calibrated end to end produces no claims; that is recorded in the recipes table.

    python scripts/build_claims.py
"""
from __future__ import annotations

import os
import sys
import time

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from starplast import claims as C  # noqa: E402
from starplast import organisms as O  # noqa: E402
from starplast import strategies as S  # noqa: E402
from starplast import track_record as T  # noqa: E402

OUT = os.path.join(ROOT, "starplast", "data", "claims.parquet")
RECIPES = os.path.join(ROOT, "starplast", "data", "claim_recipes.parquet")


def main(log=print) -> int:
    claims, recipes = [], []
    for code in O.codes(available=True):
        ctx = S.Context.shipped(code)
        for target in T.recorded_targets(code):
            t0 = time.monotonic()
            table = C.candidates(code, target)
            rec = C.recipe(code, target)
            found = C.generate(ctx, rec)
            best = table[table["generator"] == rec.generator]
            row = best.iloc[0].to_dict() if len(best) else {}
            recipes.append({
                "organism": code, "target": target, "generator": rec.generator,
                "verifiers": ", ".join(rec.verifiers), "proven": rec.proven,
                "certainty_error": rec.model.certainty.calibration_error,
                "verified_error": rec.model.calibration_error,
                "heldout_confident": row.get("confident"),
                "heldout_confident_precision": row.get("confident_precision"),
                "claims": int(len(found))})
            if len(found):
                keep = ["organism", "target", "gene_id", "claim", "generator", "certainty",
                        "verified_certainty", "status", "confidence", "prior", "lift", "distance"]
                keep += [c for c in found.columns if c.endswith(" verdict") or c.startswith("by ")]
                claims.append(found[keep])
            log(f"{code}/{target}: {rec.generator} -> {', '.join(rec.verifiers) or 'no independent check'}"
                f"; {len(found):,} claims ({time.monotonic() - t0:.0f}s)")
    out = pd.concat(claims, ignore_index=True)
    verdicts = [c for c in out.columns if c.endswith(" verdict") or c.startswith("by ")]
    for column in ["organism", "target", "claim", "generator", "status"] + verdicts:
        out[column] = out[column].astype("category")
    out.to_parquet(OUT, index=False, compression="zstd")
    pd.DataFrame(recipes).to_parquet(RECIPES, index=False)
    log(f"wrote {OUT}: {len(out):,} claims, {os.path.getsize(OUT) / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
