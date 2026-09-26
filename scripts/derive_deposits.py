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

    nb.md("## 6. Carbon-source withdrawal (Uboldi et al., bioRxiv 2025)",
          "The dependence column is the authors' contrast, glutamine-only minus glucose-only: "
          "negative means the gene is needed when glucose is absent. The genes of glutamine "
          "catabolism must come first, and the ribosome must be needed in complete medium.")
    nb.code("carbon = D.glucose_limitation(DATA)",
            "carbon['rank'] = carbon.fit_glucose_dependence.rank()",
            "known = {'TGGT1_249390': 'GDH1', 'TGGT1_289650': 'PEPCK'}",
            "print(carbon.set_index('gene_id').loc[list(known), ['fit_glucose_dependence', 'fit_glucose_dependence_fdr', 'rank']].rename(index=known).round(3))",
            "c = carbon.assign(me49=me49(carbon.gene_id).values).drop_duplicates('me49').set_index('me49')",
            "rib = ribosomal.reindex(c.index).fillna(False).astype(bool)",
            "{'ribosomal median, complete medium': c.loc[rib, 'fit_complete_medium_2025'].median(),",
            " 'others median, complete medium': c.loc[~rib, 'fit_complete_medium_2025'].median(),",
            " 'genes at FDR 0.05': int((carbon.fit_glucose_dependence_fdr < 0.05).sum())}")

    nb.md("## 7. 5' UTR architecture (Peters et al., bioRxiv 2025)",
          "Sequence features, admitted because they predict translation the way they must: more "
          "upstream AUGs, less translation; a better start context, more.")
    nb.code("utr = D.utr5_architecture(DATA).set_index('gene_id')",
            "te_mean = nodes.set_index('gene_id').filter(regex=r'^te\\d').mean(axis=1)",
            "{c: stats.spearmanr(utr[c], te_mean.reindex(utr.index), nan_policy='omit').correlation",
            " for c in ['utr5_n_uaugs', 'utr5_n_uorfs', 'utr5_length', 'utr5_kozak_score']}")

    nb.md("## 8. Bradyzoite subtypes in the brain (Ulu et al. 2026)",
          "Five groups, all of them bradyzoites: the bradyzoite markers high everywhere, the "
          "tachyzoite antigen SAG1 low everywhere, and SRS22A marking Group B.")
    nb.code("bz = D.bradyzoite_subtypes(DATA).set_index('gene_id')",
            "pct = bz.rank(pct=True)",
            "pct.loc[['TGME49_259020', 'TGME49_291040', 'TGME49_268860', 'TGME49_233460']].rename(index={'TGME49_259020': 'BAG1', 'TGME49_291040': 'LDH2', 'TGME49_268860': 'ENO1', 'TGME49_233460': 'SAG1'}).round(2)")

    nb.md("## 9. Iron depletion (Hanna et al., mBio 2026)",
          "Iron-sulfur proteins need the iron that was withdrawn, so they should shift down -- and "
          "they do, modestly; protein and transcript should move together, though not in lockstep.")
    nb.code("iron = D.iron_depletion(DATA)",
            "iron['me49'] = me49(iron.gene_id).values",
            "# Proteins arrive on GT1 accessions and transcripts on ME49 ones; pair them by gene.",
            "i = iron.groupby('me49')[['iron_depletion_protein_log2fc', 'iron_depletion_rna_log2fc']].first()",
            "# The paper's own iron-sulfur flag, rather than a guess from product names.",
            "s1 = pd.read_excel(os.path.join(DATA, *D.IRON, 'mbio.03788-25-s0002.xlsx'), sheet_name=0, header=1)",
            "fes = set(me49(s1.loc[s1['FeS'].notna(), 'Protein_Accessions']).dropna())",
            "flag = i.index.isin(fes)",
            "{'iron-sulfur proteins': int(flag.sum()),",
            " 'iron-sulfur median log2fc': i.loc[flag, 'iron_depletion_protein_log2fc'].median(),",
            " 'others median log2fc': i.loc[~flag, 'iron_depletion_protein_log2fc'].median(),",
            " 'protein vs RNA rho': stats.spearmanr(i.iron_depletion_protein_log2fc, i.iron_depletion_rna_log2fc, nan_policy='omit').correlation}")

    nb.md("## 10. Organelle surfaces (Parker & Huet, bioRxiv 2026)",
          "A bait on the outside of the mitochondrion should find mitochondrial proteins, and one on "
          "the ER should find ER proteins -- judged against hyperLOPIT, which measured location "
          "independently. The apicoplast bait is the weak one, and the notebook shows it.")
    nb.code("surf = D.organelle_surface(DATA).set_index('gene_id')",
            "comp = nodes.set_index('gene_id')['compartment'].reindex(surf.index)",
            "rows = []",
            "for bait, words in (('mitochondrion', 'mitochondri'), ('er', 'ER'), ('apicoplast', 'apicoplast')):",
            "    seen = surf[f'surface_{bait}_stringent'].notna() & comp.notna()",
            "    hit = surf[f'surface_{bait}_stringent'] == 1",
            "    lab = comp.astype(str).str.contains(words)",
            "    table = [[int((seen & hit & lab).sum()), int((seen & hit & ~lab).sum())],",
            "             [int((seen & ~hit & lab).sum()), int((seen & ~hit & ~lab).sum())]]",
            "    odds, p = stats.fisher_exact(table)",
            "    rows.append({'bait': bait, 'stringent hits': int(hit.sum()), 'odds ratio': odds, 'p': p})",
            "pd.DataFrame(rows)")

    nb.md("## 11. Host genes rhoptry discharge needs (Valleau et al., bioRxiv 2025)",
          "The screen must recover its own pathway: SLC35A2 first, the N-glycan genes near the top, "
          "and the glycosylation genes the paper rules out nowhere near it.")
    nb.code("k562 = D.k562_rhoptry_screen(DATA).drop_duplicates('host_name').set_index('host_name')",
            "rank = k562.rhoptry_discharge_score.rank(ascending=False)",
            "rank.reindex(['SLC35A2', 'RFT1', 'DPAGT1', 'GFPT1', 'MGAT1', 'MGAT2', 'B4GALT1', 'FUT8', 'ST6GAL1']).astype('Int64')")

    nb.md("## 12. Plasmodium: melting temperature, fertility, and the third wave",
          "Melting points are fitted by `scripts/fit_meltome.py`, since none is published. The "
          "checks: chaperones labile and glycolysis stable; HAP2 needed by males only; and the "
          "counts the papers report, reproduced.")
    nb.code("pf = pd.read_parquet(os.path.join(ROOT, 'starplast', 'data', 'pf_nodes.parquet')).set_index('gene_id')",
            "tm = D.pf_melting_temperature(DATA).set_index('gene_id')['melting_temperature_tm']",
            "prod = pf['product'].reindex(tm.index).fillna('')",
            "{'HSP70/90 median Tm': tm[prod.str.contains('heat shock protein 70|heat shock protein 90', case=False)].median(),",
            " 'glycolytic median Tm': tm[prod.str.contains('glyceraldehyde-3-phosphate|enolase|pyruvate kinase', case=False)].median(),",
            " 'proteins': len(tm)}")
    nb.code("D.pb_fertility(DATA).set_index('gene_id').loc[['PF3D7_1014200']].rename(index={'PF3D7_1014200': 'HAP2'})")
    nb.code("{'gametocyte proteins newly made (paper: 705)': int(D.pf_gametocyte_proteome(DATA).gametocyte_newly_made.sum()),",
            " 'latency classifier genes (paper: 200)': int(D.pf_latency(DATA).latency_classifier_member.sum()),",
            " 'H3K27ac high-confidence (paper: 99)': int(D.pf_chromatin_proxiome(DATA).chromprox_h3k27ac_hit.sum()),",
            " 'H3K4me3 high-confidence (paper: 48)': int(D.pf_chromatin_proxiome(DATA).chromprox_h3k4me3_hit.sum()),",
            " 'febrile unique sites up / down (paper rows: 143 / 53)': tuple(int(x) for x in D.pf_febrile_phospho(DATA)[['febrile_phospho_n_sites_up', 'febrile_phospho_n_sites_down']].sum()),",
            " 'm6A transcripts with a site': int((D.pf_m6a(DATA).m6a_n_canonical_sites > 0).sum()),",
            " 'proteins engaged by at least one antimalarial': int((D.pf_target_engagement(DATA).engaged_n_compounds_hit > 0).sum())}")

    nb.md("## 13. Write the tables",
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
