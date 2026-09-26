#!/usr/bin/env python3
"""Derive the 2026-09 deposit tables and write the notebook that shows how.

Runs every derivation in `starplast.deposits` against the raw deposits under the dataset root, with
the sanity checks that admitted each one, and writes both the keyed tables
(`starplast/data/deposit_<key>.tsv`) and `notebooks/derive_deposits_2026_09.ipynb`, executed, so the
notebook in the tree is the run that produced the shipped numbers.

Then `python scripts/add_deposits.py` merges the tables into the caches.

Run:  STARPLAST_DATA=/path/to/datasets python scripts/derive_deposits.py
"""
from __future__ import annotations

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from notebook_runner import ExecutedNotebook  # noqa: E402

OUT = os.path.join(ROOT, "notebooks", "derive_deposits_2026_09.ipynb")


def build(dataset_root: str) -> str:
    nb = ExecutedNotebook("Deriving the 2026-09 deposits", {"__name__": "__notebook__"})
    nb.md("Every column added in the September 2026 data audit comes from a public deposit through "
          "a computation -- a mean of replicate screens, a log, a moderated contrast. This notebook "
          "runs each one from the raw file and shows the check that admitted it. The code is "
          "`starplast/deposits.py`; this notebook calls it, so what is shown is what ships.",
          "Raw deposits live under the dataset root, in the taxonomy "
          "`<level>/<kind>/<PMID or accession>/`; every folder carries `URLS.txt` and "
          "`SHA256SUMS.txt` from the download.")
    nb.code("import warnings; warnings.filterwarnings('ignore')",
            "import os, numpy as np, pandas as pd",
            "from scipy import stats",
            f"ROOT = {ROOT!r}",
            f"DATA = {dataset_root!r}",
            "import sys; sys.path.insert(0, ROOT)",
            "from starplast import deposits as D",
            "nodes = pd.read_parquet(os.path.join(ROOT, 'starplast', 'data', 'nodes.parquet'))",
            "product = nodes.set_index('gene_id')['product']",
            "ribosomal = product.str.contains('ribosomal protein', case=False, na=False) & "
            "~product.str.contains('mitochondrial|apicoplast|kinase|methyltransferase|"
            "acetyltransferase', case=False, na=False)",
            "def me49(ids): return 'TGME49_' + pd.Series(ids, dtype=str).str.extract(r'_(\\d{6})')[0]",
            "print(len(nodes), 'genes;', int(ribosomal.sum()), 'cytosolic ribosomal proteins')")

    nb.md("## 1. In vivo fitness in six tissues (Giuliano et al. 2024, PMID 38977907)",
          "The four tissues already shipped were cited to the 2019 platform paper. If they are this "
          "study's composite scores, the supplement must reproduce them value for value. It does, "
          "and it carries heart and brain as well.")
    nb.code("invivo = D.giuliano_invivo(DATA)",
            "invivo['me49'] = me49(invivo.gene_id).values",
            "whole = invivo[~invivo.gene_id.str.contains(r'\\d[A-Z]$')].drop_duplicates('me49').set_index('me49')",
            "rows = []",
            "for c in ['fit_invivo_PE', 'fit_invivo_liver', 'fit_invivo_spleen', 'fit_invivo_lung']:",
            "    j = nodes.set_index('gene_id')[[c]].join(whole[[c]], rsuffix='_supp').dropna()",
            "    rows.append({'column': c, 'genes': len(j), 'max |difference|': (j[c] - j[c + '_supp']).abs().max(),",
            "                 'spearman': stats.spearmanr(j[c], j[c + '_supp']).correlation})",
            "pd.DataFrame(rows)")
    nb.code("pd.DataFrame({t: invivo[f'fit_invivo_{t}'].describe() for t in ['PE', 'lung', 'heart', 'brain']}).round(2)")

    nb.md("## 2. Serum restriction (Bitew et al. 2025, PMID 41407671)",
          "Two independent genome-wide screens, each in 10% and 1% serum. A screen of fibroblast "
          "fitness must find the ribosome essential; the difference between sera must NOT be "
          "fibroblast fitness again, or it is no new axis. GRA38 is the paper's gene.")
    nb.code("serum = D.serum_restriction(DATA)",
            "serum['me49'] = me49(serum.gene_id).values",
            "s = serum[~serum.gene_id.str.contains(r'\\d[A-Z]$')].drop_duplicates('me49').set_index('me49')",
            "s = s.join(nodes.set_index('gene_id')[['fit_invitro_hff']])",
            "rib = ribosomal.reindex(s.index).fillna(False).astype(bool)",
            "pd.DataFrame({c: {'ribosomal median': s.loc[rib, c].median(), 'others median': s.loc[~rib, c].median(),",
            "                  'rho with fit_invitro_hff': stats.spearmanr(s[c], s.fit_invitro_hff, nan_policy='omit').correlation}",
            "              for c in s.columns if c.startswith('fit_') and c != 'fit_invitro_hff'}).T.round(3)")
    nb.code("serum.set_index('gene_id').loc[['TGGT1_312420']].T  # GRA38")

    nb.md("## 3. Translation efficiency, GSE302107",
          "The authors' TE is a linear ratio of footprint RPKM to RNA RPKM; it is shipped as log2 to "
          "sit on the scale of the other TE columns. The test is reproducibility, agreement with the "
          "earlier studies, and the ribosome being the high-TE class.")
    nb.code("te = D.riboseq_302107(DATA).set_index('gene_id')",
            "other = nodes.set_index('gene_id').filter(regex='^te(99395_intracellular|245775_parent_tachy)')",
            "{'replicate rho': stats.spearmanr(te.te302107_tachy_r1, te.te302107_tachy_r2, nan_policy='omit').correlation,",
            " 'rho with GSE99395 mean': stats.spearmanr(te.mean(axis=1), other.filter(like='99395').mean(axis=1).reindex(te.index), nan_policy='omit').correlation,",
            " 'rho with GSE245775 mean': stats.spearmanr(te.mean(axis=1), other.filter(like='245775').mean(axis=1).reindex(te.index), nan_policy='omit').correlation,",
            " 'ribosomal median log2 TE': te.mean(axis=1)[ribosomal.reindex(te.index).fillna(False).astype(bool)].median(),",
            " 'others median log2 TE': te.mean(axis=1)[~ribosomal.reindex(te.index).fillna(False).astype(bool)].median()}")

    nb.md("## 4. mRNA decay after actinomycin D, GSE329845",
          "Raw counts, no spike-in: the result is decay relative to the median transcript. The model "
          "is intercept + treatment + replicate, moderated as in limma. Replicate agreement is "
          "checked on the per-replicate paired log ratios; stable ribosomal mRNAs are the biology "
          "check.")
    nb.code("decay, fit = D.mrna_decay_329845(DATA, return_fit=True)",
            "norm = fit['normalized']",
            "paired = pd.DataFrame({i: np.log2(norm[f'cWT_ActD_REP{i}'] + 1) - np.log2(norm[f'cWT_vehicle_REP{i}'] + 1) for i in (1, 2, 3)}).loc[decay.gene_id]",
            "print('size factors', {k: round(v, 3) for k, v in fit['size_factors'].items()})",
            "print('prior df', round(fit['df_prior'], 1))",
            "paired.corr(method='spearman').round(3)")
    nb.code("d = decay.assign(me49=me49(decay.gene_id).values).set_index('me49')['mrna_log2_remaining_4h_actinomycin']",
            "rib = ribosomal.reindex(d.index).fillna(False).astype(bool)",
            "old = nodes.set_index('gene_id')['mrna_remaining_5h_actinomycin'].reindex(d.index)",
            "{'genes': len(d), 'ribosomal median': d[rib].median(), 'others median': d[~rib].median(),",
            " 'Mann-Whitney p': stats.mannwhitneyu(d[rib], d[~rib]).pvalue,",
            " 'rho with the shipped 412-gene column': stats.spearmanr(d, old, nan_policy='omit').correlation,",
            " 'shared genes': int(old.notna().sum())}")

    nb.md("## 5. Host responses and a baseline macrophage",
          "Host tables are keyed by reviewed UniProt accession through Ensembl ids "
          "(`reference/uniprot`), never by symbol. Each contrast is infected against uninfected "
          "with the design's blocks; the checks are genes whose behaviour is not in doubt.")
    nb.code("hff = D.hff_tg_infection(DATA).drop_duplicates('host_name').set_index('host_name')",
            "print(int((hff.hff_tg_infection_padj < 0.05).sum()), 'genes at padj < 0.05 of', len(hff))",
            "hff.reindex(['CXCL8', 'IL6', 'CXCL10', 'ISG15', 'EGR1', 'GAPDH', 'ACTB']).round(3)")
    nb.code("bmdm = D.bmdm_baseline(DATA).drop_duplicates('host_name').set_index('host_name')",
            "bmdm['rank'] = bmdm.bmdm_tpm.rank(ascending=False).astype(int)",
            "bmdm.reindex(['Lyz2', 'Cd68', 'Csf1r', 'Emr1', 'Itgam', 'Alb', 'Apoa1']).round(1)")
    nb.code("hep = D.hepatocyte_pf_infection(DATA).drop_duplicates('host_name').set_index('host_name')",
            "print(int((hep.hepatocyte_pf_infection_padj < 0.05).sum()), 'genes at padj < 0.05 of', len(hep))",
            "hep.sort_values('hepatocyte_pf_infection_log2fc', ascending=False).head(10).round(3)")

    nb.md("## 6. Write the tables",
          "Each derivation is written to `starplast/data/deposit_<key>.tsv`. "
          "`scripts/add_deposits.py` merges them into the node and host tables, refusing any merge "
          "that would lose a value.")
    nb.code("written = D.derive_all(DATA, ROOT, log=print)",
            "{k: v.shape for k, v in written.items()}")
    return nb.write(OUT)


def main(argv=None) -> int:
    from starplast import paths
    root = os.environ.get("STARPLAST_DATA") or paths.dataset_root()
    print(f"dataset root: {root}")
    print(f"wrote {build(root)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
