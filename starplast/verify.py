#!/usr/bin/env python3
"""Check the shipped numbers against the repository they came from.

Task 05 asked whether the expression columns should be rebuilt from GEO's primary deposits instead of
from the papers' supplementary tables. The answer taken here is neither, because each source fails in a
way the other detects:

    supplement-derived   exactly what the authors published, but not re-derivable by a reader and not
                         checkable against anything
    GEO-primary          re-downloadable and automatable end to end, but a reprocessed series can
                         change under you silently, and GEO's matrix is not always the one the paper
                         used

So the GEO path is built and used as a CHECK. Where the two agree, the reader can re-derive the number,
which is the property a paper needs. Where they disagree, the disagreement is the finding, and it is
recorded rather than resolved by quietly preferring one.

That matters more than usual here because this exact series is the project's canonical example of a
silent error. `GSE108740_FPKM.xlsx` reaches **16,520** and is genuine FPKM; the columns derived from it
in the node table top out at **3.3** because they were logged upstream and kept the name. Log the first
again and one gene stops dominating every distance; fail to log it and one gene dominates every
distance. Neither mistake raises anything.

Comparison is on RANKS. The two forms of the same series are on different scales by construction -- that
is the whole point -- so comparing values would report a difference that is only a transform. Spearman
asks the question that actually matters: do the two sources order the genes the same way?
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# Which shipped node column corresponds to which GEO sample column, and in which series. Stated
# explicitly rather than matched by name: `expr_tachy` and `Tachyzoites_T2 [FPKM]` do not look alike,
# and a name-matching heuristic would silently compare the wrong pair and report a low correlation as
# a data problem.
CORRESPONDENCE = {
    "GSE108740": {
        "file": "transcription/RNAseq/GSE108740/GSE108740_FPKM.xlsx",
        "id_column": "ToxoDB ID",
        "pairs": {
            # shipped column   ->   GEO sample columns averaged to form it
            "expr_tachy": ("Tachyzoites_T2 [FPKM]", "Tachyzoites_T4 [FPKM]"),
            "expr_cyst": ("Tissue_cysts_A [FPKM]", "Tissue_cysts_B [FPKM]"),
        },
    },
    "GSE206344": {
        # Two workbooks ship for this series and only one is the source. `LFCs` holds seven sheets of
        # pairwise contrasts; `Normalised_data` holds the per-sample matrix on sheet "1", whose column
        # names -- Unsporulated R1, Sporulating R1 -- are exactly what the shipped rna206344_* columns
        # are named after. Declaring the wrong one would have produced a number that looks like
        # verification and is not, which is why this correspondence went unrecorded until the sheet
        # names settled it.
        "file": "transcription/RNAseq/GSE206344/GSE206344_Normalised_data_ToxoDB_Release68.xlsx",
        "sheet": "1",
        "id_column": "ToxoDB ID Release 68",
        "pairs": {
            "expr_sporulated": ("Unsporulated R1", "Unsporulated R2", "Sporulating R1",
                                "Sporulating R2", "Sporulated R1", "Sporulated R2"),
            "rna206344_Unsporulated_R1": ("Unsporulated R1",),
            "rna206344_Sporulating_R1": ("Sporulating R1",),
            "rna206344_Sporulated_R2": ("Sporulated R2",),
        },
    },
}

AGREEMENT = 0.90        # Spearman rho above which the two sources are treated as the same ordering


def _spearman(a: pd.Series, b: pd.Series) -> float:
    ok = a.notna() & b.notna()
    if ok.sum() < 30:
        return float("nan")
    return float(a[ok].rank().corr(b[ok].rank()))


def compare_series(series: str, nodes: pd.DataFrame, resolve=None, log=print) -> pd.DataFrame:
    """Compare one GEO series' primary matrix against the shipped columns derived from it."""
    from . import paths

    spec = CORRESPONDENCE[series]
    p = paths.find(spec["file"])
    if not p:
        log(f"{series}: primary matrix not found; cannot verify")
        return pd.DataFrame()

    geo = pd.read_excel(p, sheet_name=spec.get("sheet", 0))
    ids = geo[spec["id_column"]].astype(str).str.strip()
    if resolve is not None:
        ids = ids.map(lambda a: resolve(a) or a)
    geo = geo.assign(gene_id=ids).drop_duplicates("gene_id").set_index("gene_id")

    shipped = nodes.set_index("gene_id")
    rows = []
    for col, geo_cols in spec["pairs"].items():
        if col not in shipped.columns:
            continue
        have = [c for c in geo_cols if c in geo.columns]
        if not have:
            continue
        # Mean across the replicate samples, then log -- the transform the shipped column already had
        # applied upstream. Comparing on ranks makes the log irrelevant to the answer, which is the
        # point: the check should not depend on getting the transform right.
        primary = np.log2(geo[have].mean(axis=1).clip(lower=0) + 1)
        common = shipped.index.intersection(primary.index)
        rho = _spearman(shipped.loc[common, col], primary.loc[common])
        rows.append({"series": series, "column": col,
                     "geo_columns": ", ".join(have),
                     "n_common": len(common),
                     "spearman": round(rho, 4),
                     "agrees": bool(rho >= AGREEMENT),
                     "shipped_max": round(float(shipped[col].max()), 2),
                     "primary_max_raw": round(float(geo[have].max().max()), 1)})
        log(f"{series} {col}: rho={rho:.4f} over {len(common):,} genes "
            f"(shipped max {rows[-1]['shipped_max']}, GEO raw max {rows[-1]['primary_max_raw']})")
    return pd.DataFrame(rows)


def verify_all(nodes: pd.DataFrame, resolve=None, log=print) -> pd.DataFrame:
    """Every series with a recorded correspondence. The table is the result, agreement or not."""
    out = [compare_series(s, nodes, resolve=resolve, log=log) for s in CORRESPONDENCE]
    out = [t for t in out if not t.empty]
    if not out:
        return pd.DataFrame()
    t = pd.concat(out, ignore_index=True)
    n_ok = int(t.agrees.sum())
    log(f"verification: {n_ok} of {len(t)} shipped columns reproduce their GEO source "
        f"(Spearman >= {AGREEMENT})")
    if n_ok < len(t):
        log("  DISAGREEMENT is the finding here, not a reason to swap sources. Columns: "
            + ", ".join(t[~t.agrees].column))
    return t
