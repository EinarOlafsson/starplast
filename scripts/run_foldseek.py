#!/usr/bin/env python3
"""All-against-all structural comparison of a proteome's AlphaFold models, as a pair table.

The Toxoplasma arm carries a `struct` layer -- Foldseek TM-align over 6,900 models at TM >= 0.7 --
and it is the layer that reaches what the others cannot: structural similarity needs no orthology, so
it finds relatives among the lineage-specific proteins where sequence search gives nothing. The
Plasmodium arm had no such layer because nobody had run it.

This is a BUILD step, not part of the node build: it wants a few gigabytes of models and several
minutes of CPU, and it produces a table that the loader then reads in a second. That split is
deliberate -- the expensive thing is archived as data, and everything downstream stays testable.

    python scripts/run_foldseek.py --models <dir of .cif.gz> --out datasets/.../pf_struct_pairs.tsv

Every pair the search reports is written, with its TM-score and its alignment; the THRESHOLD is
applied by the loader rather than here, so raising or lowering it later does not mean re-running an
hour of compute.
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

#: What AlphaFold names a model: `AF-<accession>-F1-model_v6.cif.gz`. The accession is what the
#: identity layer joins on, so it is pulled out here rather than left in the file name.
MODEL = re.compile(r"AF-([A-Z0-9]+)-F(\d+)-model")

#: Foldseek's own output fields. `alntmscore` is the TM-score of the alignment, which is what the
#: Toxoplasma layer thresholds on; `lddt` and `prob` come along because they cost nothing and a pair
#: table nobody can re-judge is a pair table that has to be recomputed to be re-judged.
FIELDS = "query,target,alntmscore,qtmscore,ttmscore,lddt,prob,evalue,alnlen"


def run(models: str, out: str, binary: str = "foldseek", threads: int = 0,
        exhaustive: bool = False, log=print) -> int:
    """Search every model against every other, writing one row per reported pair.

    `exhaustive` skips the prefilter and TM-aligns all N^2 pairs. It is off by default because it
    does not finish: 5,168 models is 26.7 million alignments, and the run was still inside the
    tmalign stage when a 90-minute budget expired, having written nothing -- foldseek holds the
    output until the stage completes, so a killed exhaustive run costs the whole search and leaves
    no partial result to resume from. The prefilter is the right trade at THIS threshold rather
    than in general: the layer keeps pairs at TM >= 0.7, which is close structural similarity, and
    the 3Di k-mer prefilter is built to retain exactly those. Exhaustive search earns its cost in
    the remote-homology range around TM 0.4, which this layer does not ship.
    """
    if not os.path.isdir(models):
        log(f"no model directory at {models}")
        return 1
    if shutil.which(binary) is None and not os.path.exists(binary):
        log(f"foldseek not found at {binary!r} -- install it from github.com/steineggerlab/foldseek")
        return 1
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="foldseek-") as tmp:
        raw = os.path.join(tmp, "hits.tsv")
        command = [binary, "easy-search", models, models, raw, os.path.join(tmp, "work"),
                   # TM-align mode. Foldseek's default 3Di+AA alignment is far faster and reports a
                   # different quantity; the Toxoplasma layer is TM-score, and two arms whose
                   # "structural similarity" means two different numbers cannot be compared.
                   "--alignment-type", "1",
                   "--format-output", FIELDS,
                   # Foldseek's own default e-value, and a cap on how many targets per query reach
                   # the aligner. Both were measured rather than chosen: at `-e 10` -- ten thousand
                   # times looser, which an earlier version used -- almost every target survives to
                   # be TM-aligned, and 200 queries against this database did not finish in ten
                   # minutes, putting the whole proteome past four hours. At these settings the
                   # same 200 queries take 160 seconds, which is 69 minutes for all 5,168.
                   #
                   # The cap costs nothing here, and that is checked rather than assumed: the
                   # busiest of those 200 queries returned 159 rows, so no query is truncated at
                   # 300. Raise it if a proteome with large structural families says otherwise --
                   # the tell is queries returning exactly `--max-seqs` rows.
                   "-e", "0.001",
                   "--max-seqs", "300"]
        if exhaustive:
            command += ["--exhaustive-search", "1"]
        if threads:
            command += ["--threads", str(threads)]
        log(" ".join(command))
        # Streamed, not captured. Foldseek holds its output until a stage completes, so a run that
        # is killed part way -- which is how the first two attempts ended -- reported an empty log
        # and left nothing to diagnose. Letting it write through means a timeout at least says
        # which stage it died in.
        result = subprocess.run(command, stdout=sys.stdout, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            log(f"foldseek exited {result.returncode}; its output is above")
            return result.returncode
        kept = 0
        with open(raw, encoding="utf8") as fh, open(out, "w", encoding="utf8") as dest:
            dest.write("accession_a\taccession_b\t" + "\t".join(FIELDS.split(",")[2:]) + "\n")
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < len(FIELDS.split(",")):
                    continue
                a, b = MODEL.search(parts[0]), MODEL.search(parts[1])
                if not a or not b or a.group(1) == b.group(1):
                    continue
                dest.write("\t".join([a.group(1), b.group(1)] + parts[2:]) + "\n")
                kept += 1
    log(f"wrote {kept:,} pairs to {out}")
    return 0


def main(argv=None, log=print) -> int:
    """Parse arguments and run the search."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", required=True, help="directory of AlphaFold model files")
    parser.add_argument("--out", required=True, help="TSV to write")
    parser.add_argument("--binary", default="foldseek", help="path to the foldseek binary")
    parser.add_argument("--threads", type=int, default=0, help="0 leaves the choice to foldseek")
    parser.add_argument("--exhaustive", action="store_true",
                        help="TM-align all N^2 pairs, skipping the prefilter (does not finish "
                             "for a whole proteome -- see run()'s docstring)")
    args = parser.parse_args(argv)
    return run(args.models, args.out, binary=args.binary, threads=args.threads,
               exhaustive=args.exhaustive, log=log)


if __name__ == "__main__":
    sys.exit(main())
