#!/usr/bin/env python3
"""Build the shipped track record: every labelled gene held out, for every strategy that can say.

Writes `starplast/data/track_record.parquet`: one row per (strategy, held-out gene) for the default
label of each organism, plus the set hold-outs (one per class, and a random-size degradation curve).
Those are what the application reads to answer, for one gene, "could this software have told me what
I already know?".

Every label in `track_record.labels` is recorded: the organism's declared targets
(`organisms.Space.targets`) and, for Toxoplasma, the cell-cycle phase and the four specific screen
phenotypes. Columns that
describe where a label came from rather than biology (`compartment_source`, `ortholopit_donors`,
`screen_scorers_agree`, `chromosome`, ...) are deliberately not: a strategy "recovering" which
dataset assigned a compartment would be a number with no meaning. `track_record.alone` answers for
any other label on demand.

    python scripts/build_track_record.py --workers 4     # both organisms, every recorded label
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


def labels(code: str) -> list:
    """The labels recorded for one organism, the default first (`track_record.labels`)."""
    return T.labels(S.Context.shipped(code))


def _task(args):
    code, target = args
    t0 = time.monotonic()
    lines = []
    frame = build(code, target, log=lines.append)
    return code, target, frame, lines, time.monotonic() - t0


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
    ap.add_argument("--workers", type=int, default=1,
                    help="labels built in parallel; each worker holds one organism's table")
    args = ap.parse_args(argv)
    codes = args.organism or O.codes(available=True)
    tasks = [(code, t) for code in codes
             for t in ([args.target] if args.target else labels(code))]
    log(f"{len(tasks)} labels: " + ", ".join(f"{c}/{t}" for c, t in tasks))
    frames = []
    if args.workers > 1:
        from concurrent.futures import ProcessPoolExecutor
        import multiprocessing
        with ProcessPoolExecutor(args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
            for code, target, frame, lines, seconds in pool.map(_task, tasks):
                log(f"{code}/{target}: {seconds:.0f}s")
                for line in lines:
                    log(line)
                frames.append(frame)
    else:
        for code, target in tasks:
            log(f"{code}/{target}:")
            frames.append(build(code, target, log=log))
    out = pd.concat([f for f in frames if len(f)], ignore_index=True)
    # Labels repeat enormously; as categories the file is a fraction of the size.
    for column in ("organism", "strategy", "target", "setting_key", "mode", "set_name",
                   "truth", "prediction"):
        out[column] = out[column].astype("category")
    out.to_parquet(args.out, index=False, compression="zstd")
    log(f"wrote {args.out}: {len(out):,} rows, {os.path.getsize(args.out) / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
