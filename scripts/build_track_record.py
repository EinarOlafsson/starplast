#!/usr/bin/env python3
"""Build the shipped track record: every labelled gene held out, for every strategy that can say.

Writes `starplast/data/track_record.parquet`: one row per (strategy, held-out gene) for the default
label of each organism, plus the set hold-outs (one per class, and a random-size degradation curve).
Those are what the application reads to answer, for one gene, "could this software have told me what
I already know?".

Other labels are not shipped: 17 of them per organism multiply the build by seventeen and the file
with it. `starplast.track_record.evaluate` computes any of them on demand.

    python scripts/build_track_record.py                 # both organisms, the default label
    python scripts/build_track_record.py --organism Tg --target compartment --out /tmp/t.parquet
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from starplast import organisms as O  # noqa: E402
from starplast import strategies as S  # noqa: E402
from starplast import track_record as T  # noqa: E402

OUT = os.path.join(ROOT, "starplast", "data", "track_record.parquet")

#: Sizes for the degradation curve: one gene alone is the easiest case there is, five hundred hidden
#: together is the case a real screen presents.
CURVE_SIZES = (1, 5, 20, 100, 500)
CURVE_REPEATS = 4

#: A strategy whose five folds take longer than this gets complete fold coverage but no set
#: hold-outs: the sets are another ~45 passes, and for the slowest two that is most of an
#: hour for a curve that is flat wherever it has been measured.
SET_BUDGET_SECONDS = 30.0


def build(organism: str, target: str | None = None, log=print) -> pd.DataFrame:
    """Every supported strategy, held out fold by fold and by set, for one organism."""
    ctx = S.Context.shipped(organism)
    target = target or S.default_category(ctx)
    parts = []
    for key in T.supported():
        t0 = time.monotonic()
        folds = T.evaluate(ctx, key, target)
        if not len(folds):
            log(f"  {key:24s} cannot speak about {target}; skipped")
            continue
        fold_seconds = time.monotonic() - t0
        parts.append(folds)
        # Every strategy gets complete fold coverage; the set hold-outs are another ~45 runs each,
        # which for the slowest two would be most of an hour for a curve already known to be flat.
        # They are run for the strategies that can afford them, and the ledger says plainly which.
        if fold_seconds <= SET_BUDGET_SECONDS:
            parts.append(T.evaluate_sets(ctx, key, target, T.class_sets(ctx, target)))
            parts.append(T.evaluate_sets(ctx, key, target,
                                         T.random_sets(ctx, target, CURVE_SIZES, CURVE_REPEATS,
                                                       seed=1)))
        else:
            log(f"  {key:24s} folds only: {fold_seconds:.0f}s a pass is too slow for the sets")
        answered = int((~folds["abstained"].astype(bool)).sum())
        right = int(folds["correct"].fillna(False).sum())
        log(f"  {key:24s} {time.monotonic() - t0:6.1f}s  {right:>5,}/{answered:<6,} right "
            f"({right / answered:.3f})" if answered else f"  {key:24s} no answers")
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=list(T.COLUMNS))


def main(argv=None, log=print) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--organism", action="append", choices=O.codes())
    ap.add_argument("--target", default=None)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args(argv)
    frames = []
    for code in args.organism or O.codes(available=True):
        log(f"{code}:")
        frames.append(build(code, args.target, log=log))
    out = pd.concat(frames, ignore_index=True)
    # Labels repeat enormously; as categories the file is a fraction of the size.
    for column in ("organism", "strategy", "target", "setting_key", "mode", "set_name",
                   "truth", "prediction"):
        out[column] = out[column].astype("category")
    out.to_parquet(args.out, index=False, compression="zstd")
    log(f"wrote {args.out}: {len(out):,} rows, {os.path.getsize(args.out) / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
