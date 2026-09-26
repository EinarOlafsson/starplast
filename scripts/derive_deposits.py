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

    nb.md("## 13. The spatial proteome of the schizont (Chisholm et al. 2026, PMID 42218142)",
          "The first subcellular localization in the Plasmodium table. Three numbers have to come "
          "back exactly -- 3,000 proteins mapped, 1,646 classified, 24 niches -- and then the "
          "markers have to be where cell biology has put them for thirty years.")
    nb.code("lopit = D.pf_spatial_proteome(DATA).set_index('gene_id')",
            "pfid = pd.read_csv(os.path.join(ROOT, 'starplast', 'data', 'plasmodb_identity.tsv'), sep='\\t')",
            "symbol = pfid.dropna(subset=['gene_name']).drop_duplicates('gene_name').set_index('gene_name')['gene_id']",
            "{'proteins mapped (paper: 3000)': int(lopit.lopit_pf_svm_score.notna().sum()),",
            " 'classified (paper: 1646)': int(lopit.lopit_pf_location.notna().sum()),",
            " 'niches (paper: 24)': int(lopit.lopit_pf_location.nunique())}")
    nb.code("markers = ['RAP1', 'MAHRP1', 'ACP', 'EXP2', 'HSP101', 'GAPDH', 'SERA5']",
            "pd.DataFrame([{'marker': m, 'gene': symbol.get(m),",
            "               'location': lopit.lopit_pf_location.get(symbol.get(m)),",
            "               'svm score': lopit.lopit_pf_svm_score.get(symbol.get(m))}",
            "              for m in markers])")

    nb.md("## 14. Field variation and between-species dN/dS, from the same paper",
          "Two different questions, which is why they are a separate deposit. The field pN/pS must "
          "NOT be a copy of the laboratory-strain SNP ratio the table already carries -- if it were, "
          "the field-variation slot would be answered by data the map already had.")
    nb.code("var = D.pf_field_variation(DATA).set_index('gene_id')",
            "base = pf[['snp_nonsyn_syn_ratio', 'piggybac_mis', 'ortholog_number']]",
            "j = var.join(base)",
            "rows = []",
            "for c in ['field_pnps_adj', 'field_variant_fraction', 'dnds_laverania', 'dnds_plasmodium']:",
            "    for other in ['snp_nonsyn_syn_ratio', 'piggybac_mis']:",
            "        ok = j[[c, other]].dropna()",
            "        rows.append({'column': c, 'against': other, 'genes': len(ok),",
            "                     'rho': stats.spearmanr(ok[c], ok[other]).correlation})",
            "pd.DataFrame(rows).round(3)")

    nb.md("## 15. mRNA synthesis and decay rates (Painter et al. 2018, PMID 29985403)",
          "Fluxes in transcripts per minute, not half-lives. The paper's claim is that transcription "
          "runs in every stage of the cycle; the deposit's own peak-stage labels say how many genes "
          "peak in each window, and the timing is then checked against an expression series this "
          "study had nothing to do with.")
    nb.code("dyn = D.pf_mrna_dynamics(DATA).set_index('gene_id')",
            "raw = pd.read_excel(os.path.join(DATA, *D.PF_4TU), sheet_name='Transcription Rate')",
            "raw.columns = ['gene_id', 'stage', 'rate']",
            "print(raw.stage.value_counts().to_dict())",
            "{'transcription rate > 0': int((dyn.transcription_rate_4tu > 0).sum()),",
            " 'of': int(dyn.transcription_rate_4tu.notna().sum()),",
            " 'decay rate < 0': int((dyn.mrna_decay_rate_4tu < 0).sum()),",
            " 'of ': int(dyn.mrna_decay_rate_4tu.notna().sum())}")
    nb.code("expr = pf[['expr_ring', 'expr_early_trophozoite', 'expr_late_trophozoite', 'expr_schizont']]",
            "peak = expr.idxmax(axis=1).reindex(raw.gene_id.values)",
            "ring = raw[raw.stage.isin(['Early Ring', 'Mid Ring', 'Late Ring'])]",
            "share = {s: float((peak.reindex(g.gene_id.values) == 'expr_ring').mean())",
            "         for s, g in ring.groupby('stage')}",
            "mx = expr.max(axis=1).reindex(dyn.index)",
            "{'share peaking in the shipped ring column': {k: round(v, 3) for k, v in share.items()},",
            " 'rho(transcription rate, max stage expression)': round(float(stats.spearmanr(dyn.transcription_rate_4tu, mx, nan_policy='omit').correlation), 3),",
            " 'rho(|decay rate|, max stage expression)': round(float(stats.spearmanr(dyn.mrna_decay_rate_4tu.abs(), mx, nan_policy='omit').correlation), 3)}")

    nb.md("## 16. The blood-stage proteome and Hsp90 dependence (bioRxiv 2026)",
          "The paper lists 131 chaperone-dependent proteins. Its methods give the rule, so the rule "
          "is applied to the abundance table and the answer is compared with their list -- not "
          "copied from it. Two inhibitors with unrelated scaffolds, so what agrees between them is "
          "the chaperone's.")
    nb.code("hsp = D.pf_hsp90_chemoproteome(DATA).set_index('gene_id')",
            "s1 = pd.read_excel(os.path.join(DATA, *D.PF_HSP90, D.PF_HSP90_TABLE), sheet_name=0, header=1)",
            "s2 = pd.read_excel(os.path.join(DATA, *D.PF_HSP90, D.PF_HSP90_HITS), sheet_name='A. PfHsp90-dependent proteins', header=1)",
            "num = lambda c: pd.to_numeric(s1[c], errors='coerce')",
            "rule = ((num('Average lg2fc (GA)') < D.PF_HSP90_MIN_DROP) & (num('p-value (GA)') < D.PF_HSP90_MAX_P)",
            "        & (num('Average lg2fc (XL)') < D.PF_HSP90_MIN_DROP) & (num('p-value (XL)') < D.PF_HSP90_MAX_P))",
            "mine = set(s1.loc[rule, 'UniProt ID'].astype(str).str.strip())",
            "theirs = set(s2['UniProt ID'].astype(str).str.strip())",
            "{'the paper\\'s list': len(theirs), 'the rule applied here': len(mine),",
            " 'in both': len(mine & theirs), 'genes after mapping': int(hsp.hsp90_dependent.sum()),",
            " 'rho(GA, XL)': round(float(stats.spearmanr(hsp.hsp90_inhibition_ga_log2fc, hsp.hsp90_inhibition_xl_log2fc).correlation), 3)}")
    nb.code("prod = pf['product'].reindex(hsp.index).fillna('')",
            "dependent = hsp.hsp90_dependent == 1",
            "{'proteins quantified': len(hsp),",
            " 'median log2 abundance, dependent': round(float(hsp.loc[dependent, 'proteome_blood_log2'].median()), 2),",
            " 'median log2 abundance, the rest': round(float(hsp.loc[~dependent, 'proteome_blood_log2'].median()), 2),",
            " 'ribosomal proteins among the dependent': int(prod[dependent].str.contains('ribosomal protein', case=False).sum())}")

    nb.md("## 17. RNA dependence (Hollin et al. 2024, PMID 38355719)",
          "A protein whose sedimentation moves to the light fractions once the RNA is digested was "
          "being held there by RNA. The paper reports 898 of them; the q-values in the deposit have "
          "to give that number, and the classes that should be RNA-dependent have to be.")
    nb.code("rdeep = D.pf_rna_dependence(DATA).set_index('gene_id')",
            "prod = pf['product'].reindex(rdeep.index).fillna('')",
            "dep = rdeep.rna_dependent == 1",
            "families = {'ribosomal protein': 'ribosomal protein', 'RNA helicase': 'helicase',",
            "            'tRNA ligase': 'tRNA ligase', 'proteasome': 'proteasome subunit'}",
            "rows = [{'quantified (paper: 3671)': len(rdeep), 'RNA-dependent (paper: 898)': int(dep.sum())}]",
            "for name, word in families.items():",
            "    hit = prod.str.contains(word, case=False)",
            "    table = [[int((hit & dep).sum()), int((hit & ~dep).sum())],",
            "             [int((~hit & dep).sum()), int((~hit & ~dep).sum())]]",
            "    odds, p = stats.fisher_exact(table)",
            "    rows.append({'family': name, 'n': int(hit.sum()), 'dependent': int((hit & dep).sum()),",
            "                 'odds': round(odds, 2), 'p': p})",
            "pd.DataFrame(rows)")

    nb.md("## 18. The sexually committed proteome (Venugopal et al. 2026, PMID 41482054)",
          "MSRP1 must be up, and that proves nothing -- it is the protein the sort was done on. The "
          "check is the paper's finding that merozoite surface proteins separate the two "
          "populations.")
    nb.code("committed = D.pf_committed_proteome(DATA).set_index('gene_id')",
            "wanted = ['MSRP1', 'MSP1', 'MSP2', 'MSP3', 'AP2-G', 'GEXP5']",
            "pd.DataFrame([{'protein': w, 'gene': symbol.get(w),",
            "               'log2FC': committed.committed_vs_asexual_log2fc.get(symbol.get(w)),",
            "               'FDR': committed.committed_vs_asexual_fdr.get(symbol.get(w))}",
            "              for w in wanted]).round(4)")

    nb.md("## 19. The resistome (Luth et al. 2024, PMID 39607932)",
          "The design reproduces from the deposit -- 724 clones, 118 compounds -- and then the two "
          "kinds of column part company. Counting mutations puts culture-adaptation loci on top; the "
          "paper's own classification puts the antimalarial targets there.")
    nb.code("res = D.pf_resistome(DATA).set_index('gene_id')",
            "snv = pd.read_excel(os.path.join(DATA, *D.PF_RESISTOME, D.PF_RESISTOME_SNVS), header=1)",
            "print({'clones (paper: 724)': snv['Clone Name'].nunique(), 'compounds (paper: 118)': snv['Compound'].nunique()})",
            "product = pf['product']",
            "top = res.sort_values('resistance_selection_compounds', ascending=False).head(6)",
            "top.assign(product=product.reindex(top.index))[['resistance_selection_compounds', 'resistance_selection_clones', 'product']]")
    nb.code("called = res[res.resistance_target_compounds.notna()].sort_values('resistance_target_compounds', ascending=False)",
            "called.assign(product=product.reindex(called.index))[['resistance_target_compounds', 'resistance_selection_clones', 'product']].head(10)")
    nb.code("{'genes with a selected coding mutation': int(res.resistance_selection_clones.notna().sum()),",
            " 'genes the paper calls a target or resistance gene': int(res.resistance_target_compounds.notna().sum()),",
            " 'field dN/dS genes': int(res.pf6_field_dnds.notna().sum()),",
            " 'rho(field dN/dS, selected compounds)': round(float(stats.spearmanr(res.pf6_field_dnds, res.resistance_selection_compounds, nan_policy='omit').correlation), 3)}")

    nb.md("## 20. Write the tables",
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
