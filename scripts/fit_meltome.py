#!/usr/bin/env python3
"""PXD056075 / MAP-X (PMID 41315737): per-protein melting temperature across 7 IDC timepoints x 3 reps.

Input : Nat Microbiol 2025 Supplementary Table 1 (41564_2025_2173_MOESM3_ESM.xlsx, cached as
        supptable1.parquet) -- Proteome Discoverer protein abundances in 10 TMT channels = 10
        temperatures 37..73 C (4 C steps), one row per protein x timepoint x replicate, organism Pf/Hs.
Method: TPP-style (Savitski 2014), matching the Starplast Toxoplasma column's conventions.
    1. Pf rows, protein FDR confidence High or Medium (q < 0.05), all 10 channels present.
    2. Fraction soluble = abundance_T / abundance_37.
    3. Per run: median fraction per temperature over all proteins -> fit sigmoid -> per-temperature
       correction factor fitted/median (the authors' X.scale.meltcurves does the same median-curve
       correction, then min-max scales; we keep the 37 C anchor instead of min-max scaling).
    4. Per protein: f(T) = (1-p)/(1+exp(-(a/T - b))) + p; Tm solves f = 0.5.
       Keep R2 >= 0.8, plateau p < 0.3 (curve crosses 0.5), 30 <= Tm <= 80 (as the Tg column).
Output: derived/pxd056075_tm_long.tsv (one row per protein x run) and derived/pxd056075_tm_gene.tsv.
"""
import os, warnings
import numpy as np, pandas as pd
from scipy.optimize import curve_fit

#: Where the deposit and its fitted intermediate live under the dataset root. `starplast.deposits`
#: reads the fitted table; this script is what regenerates it from the published supplement.
FOLDER = os.path.join("translation", "proteomics", "41315737")
HERE = os.path.dirname(os.path.abspath(__file__))
TEMPS = np.array([37, 41, 45, 49, 53, 57, 61, 65, 69, 73], float)
TC = [str(t) for t in TEMPS.astype(int)]
warnings.filterwarnings("ignore")


def sig(T, a, b, p):
    return (1 - p) / (1 + np.exp(-(a / T - b))) + p


def fit(y):
    best = None
    for a0, b0 in [(550, 10), (800, 15), (1500, 28), (300, 5)]:
        try:
            popt, _ = curve_fit(sig, TEMPS, y, p0=[a0, b0, 0.0], bounds=([0, -50, -0.5], [1e5, 1e3, 1.5]), maxfev=5000)
        except Exception:
            continue
        r = y - sig(TEMPS, *popt)
        r2 = 1 - (r ** 2).sum() / ((y - y.mean()) ** 2).sum()
        if best is None or r2 > best[1]:
            best = (popt, r2)
    if best is None:
        return np.nan, np.nan, np.nan
    (a, b, p), r2 = best
    if p >= 0.5:
        return np.nan, r2, p
    den = b - np.log(0.5 / (0.5 - p))
    tm = a / den if den > 0 else np.nan
    return tm, r2, p


def main():
    t = pd.read_parquet(os.path.join(HERE, "supptable1.parquet"))
    t = t[(t.organism == "Pf") & t["Protein FDR Confidence: Combined"].isin(["High", "Medium"])].copy()
    t["gene_id"] = t.Accession.str.replace(r"\.\d+-p\d+$", "", regex=True)
    t = t.dropna(subset=TC)
    t = t[(t[TC] > 0).all(axis=1)]
    # 7 duplicate rows within a run: average them (as the authors' X.clean.meltcurves does)
    t = t.groupby(["timepoint", "replicate", "gene_id"], as_index=False)[TC + ["# PSMs", "# Unique Peptides"]].mean()
    rel = t[TC].div(t["37"], axis=0)
    rows = []
    for (tp, rep), idx in t.groupby(["timepoint", "replicate"]).groups.items():
        R = rel.loc[idx]
        med = R.median().values
        pm, _ = curve_fit(sig, TEMPS, med, p0=[550, 10, 0], bounds=([0, -50, -0.5], [1e5, 1e3, 1.5]), maxfev=5000)
        corr = sig(TEMPS, *pm) / med
        corr[0] = 1.0
        N = R * corr
        for i, y in zip(idx, N.values):
            tm, r2, p = fit(y)
            rows.append((t.at[i, "gene_id"], tp, rep, tm, r2, p, t.at[i, "# PSMs"], t.at[i, "# Unique Peptides"]))
        print(tp, rep, "proteins", len(idx), "median-curve Tm", round(pm[0] / (pm[1] - np.log(0.5 / (0.5 - pm[2]))), 2))
    L = pd.DataFrame(rows, columns=["gene_id", "timepoint", "replicate", "tm", "r2", "plateau", "psms", "unique_peptides"])
    L["hpi"] = L.timepoint.str.replace("hpi", "").astype(int)
    L["pass"] = (L.r2 >= 0.8) & (L.plateau < 0.3) & L.tm.between(30, 80)
    os.makedirs(os.path.join(HERE, "derived"), exist_ok=True)
    L.to_csv(os.path.join(HERE, "derived", "pxd056075_tm_long.tsv"), sep="\t", index=False)
    print("curves", len(L), "pass", L["pass"].sum(), "(%.1f%%)" % (100 * L["pass"].mean()))
    P = L[L["pass"]]
    # per gene x timepoint mean, then gene summary
    tpm = P.groupby(["gene_id", "hpi"]).tm.mean().unstack()
    tpm.columns = ["tm_%dhpi" % c for c in tpm.columns]
    g = P.groupby("gene_id").agg(tm_median=("tm", "median"), tm_mean=("tm", "mean"), tm_sd=("tm", "std"),
                                 n_curves=("tm", "size"), n_timepoints=("hpi", "nunique"))
    g = g.join(tpm)
    g.to_csv(os.path.join(HERE, "derived", "pxd056075_tm_gene.tsv"), sep="\t")
    print("genes with >=1 passing curve:", len(g), " >=2 curves:", (g.n_curves >= 2).sum(), " >=3:", (g.n_curves >= 3).sum())
    print(g.describe().round(2).to_string())


if __name__ == "__main__":
    main()
