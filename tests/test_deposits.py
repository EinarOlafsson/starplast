"""The 2026-09 deposits: the derivations say what they mean, and the merge cannot lose anything.

Planted tables first, because a derivation that cannot be shown to get a known answer right is not
evidence about a real one. Then the shipped tables: every column the audit added is there, in the
units its registry entry claims, and the three facts the audit turned on stay true -- the in vivo
scores are Giuliano's, the serum differential is a new axis rather than fibroblast fitness again,
and the genome-wide decay column is not the old 412-gene one.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import deposits as D  # noqa: E402


# --------------------------------------------------------------------------- the statistics
def test_moderated_t_finds_a_planted_difference_and_not_a_planted_null():
    rng = np.random.default_rng(0)
    Y = rng.normal(size=(400, 6))
    Y[:50, 3:] += 3.0                                  # fifty genes really do differ
    X = np.column_stack([np.ones(6), [0, 0, 0, 1, 1, 1]]).astype(float)
    coef, _t, p, info = D.moderated_t(Y, X, [0, 1])
    assert coef[:50].mean() > 2.5 and abs(coef[50:].mean()) < 0.3
    q = D.bh(p)
    # Three replicates against three is not much power: a three-sigma difference is found in about
    # three quarters of the planted genes, and expecting more would be expecting the method to beat
    # its own design. What matters is that the false discoveries stay near the rate asked for.
    assert (q[:50] < 0.05).mean() > 0.6
    assert (q[50:] < 0.05).mean() < 0.05
    assert info["df_residual"] == 4 and info["df_prior"] > 0


def test_squeezing_the_variance_is_what_makes_three_replicates_usable():
    """A gene whose three values happen to agree must not be declared certain.

    The point of the moderated test: one gene given an artificially tiny spread should not outrank
    genes with a much larger difference, which an ordinary t-test would let it do.
    """
    from scipy import stats
    rng = np.random.default_rng(1)
    Y = rng.normal(size=(200, 6))
    Y[:20, 3:] += 2.0                                       # real, ordinary spread
    Y[20] = [0, 0, 0, 0.2, 0.2, 0.2]                        # tiny difference, tiny spread
    X = np.column_stack([np.ones(6), [0, 0, 0, 1, 1, 1]]).astype(float)
    _c, t_mod, _p, _i = D.moderated_t(Y, X, [0, 1])
    t_plain = stats.ttest_ind(Y[:, :3], Y[:, 3:], axis=1).statistic
    assert abs(t_plain[20]) > np.abs(t_plain[:20]).max()     # the plain test is fooled
    assert abs(t_mod[20]) < np.abs(t_mod[:20]).max()         # the moderated one is not


def test_bh_keeps_nan_and_is_monotone():
    q = D.bh([0.001, 0.02, np.nan, 0.5])
    assert np.isnan(q[2]) and q[0] <= q[1] <= q[3] and q[3] <= 1.0


# --------------------------------------------------------------------------- the merge
def _table(n=40):
    return pd.DataFrame({"gene_id": [f"TGME49_{200000 + 10 * i}" for i in range(n)]})


def test_a_whole_locus_beats_the_halves_gt1_splits(tmp_path, monkeypatch):
    """TGGT1_224540 and TGGT1_224540A/B all resolve to one ME49 gene; the whole one wins.

    Averaging all three moved one gene's in vivo score from -3.99 to +0.20 once. The A/B halves are
    different gene models, not replicates of the whole, so the accession that IS the gene is taken
    and the halves are used only where nothing else reaches that gene.
    """
    frame = pd.DataFrame({"gene_id": ["TGGT1_224540", "TGGT1_224540A", "TGGT1_224540B",
                                      "TGGT1_300000A", "TGGT1_300000B"],
                          "value": [-4.0, 0.3, 4.3, 1.0, 3.0]})
    frame.to_csv(tmp_path / "deposit_planted.tsv", sep="\t", index=False)
    monkeypatch.setattr(D, "DEPOSITS", (D.Deposit("planted", "Tg", lambda root: frame),))
    monkeypatch.setattr(D, "table_path",
                        lambda base, key: str(tmp_path / f"deposit_{key}.tsv"))
    ids = ["TGME49_224540", "TGME49_300000"]
    out = D.parasite_columns(str(tmp_path), "Tg", ids,
                             resolve=lambda a: "TGME49_" + a.split("_")[1].rstrip("AB"),
                             log=lambda *a: None)
    assert out["value"].iloc[0] == -4.0                 # the whole locus, not the mean of three
    assert out["value"].iloc[1] == 2.0                  # halves only: their mean is all there is


def test_a_column_that_would_lose_values_is_refused(tmp_path, monkeypatch):
    import scripts.add_deposits as A
    nodes = _table()
    nodes["value"] = 1.0
    path = tmp_path / "nodes.parquet"
    nodes.to_parquet(path, index=False)
    thin = pd.DataFrame({"gene_id": nodes["gene_id"][:5], "value": 2.0})
    monkeypatch.setattr(D, "DEPOSITS", (D.Deposit("planted", "Tg", lambda root: thin),))
    monkeypatch.setattr(D, "_read", lambda base, key: thin)
    said = []
    assert A.merge_parasite(str(path), "Tg", str(tmp_path), log=said.append) == 2
    assert any("REFUSED" in line for line in said)
    assert pd.read_parquet(path)["value"].notna().sum() == len(nodes)   # untouched


def test_an_empty_column_is_refused_rather_than_written(tmp_path, monkeypatch):
    """A column of nothing flips a slot to filled and adds a feature nobody measured."""
    import scripts.add_deposits as A
    nodes = _table()
    nodes.to_parquet(tmp_path / "nodes.parquet", index=False)
    empty = pd.DataFrame({"gene_id": ["TGME49_999999"], "value": [1.0]})
    monkeypatch.setattr(D, "DEPOSITS", (D.Deposit("planted", "Tg", lambda root: empty),))
    monkeypatch.setattr(D, "_read", lambda base, key: empty)
    assert A.merge_parasite(str(tmp_path / "nodes.parquet"), "Tg", str(tmp_path),
                            log=lambda *a: None) == 2


def test_a_missing_deposit_leaves_what_is_shipped_alone(tmp_path):
    """A machine without the raw files must not erase what a machine with them derived."""
    written = D.derive_all(str(tmp_path / "nothing-here"), str(tmp_path), log=lambda *a: None)
    assert written == {}


def test_every_registered_deposit_is_in_the_dataset_registry():
    from starplast import datasets
    keys = {d.key for d in datasets.REGISTRY}
    assert {d.key for d in D.DEPOSITS} <= keys


def test_every_deposit_column_is_claimed_by_its_registry_entry():
    from starplast import datasets
    for dep in D.DEPOSITS:
        frame = D._read(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), dep.key)
        if frame.empty:
            continue
        declared = set(datasets.get(dep.key).columns)
        produced = set(frame.columns) - {"gene_id", "host_id", "host_name"}
        assert produced <= declared, (dep.key, produced - declared)


# --------------------------------------------------------------------------- the shipped tables
@pytest.fixture(scope="module")
def tg():
    from starplast import paths
    return pd.read_parquet(paths.cache_file("nodes.parquet"))


@pytest.fixture(scope="module")
def pf():
    from starplast import paths
    return pd.read_parquet(paths.cache_file("pf_nodes.parquet"))


@pytest.fixture(scope="module")
def hosts():
    from starplast import paths
    return pd.read_parquet(paths.cache_file("host_proteins.parquet"))


def test_the_new_in_vivo_tissues_are_there_and_agree_with_the_four_that_were(tg):
    """Heart and brain come from the same sheet as the four already shipped, so they must behave
    like them -- correlated, since a gene that cannot survive anywhere cannot survive here."""
    for c in ("fit_invivo_heart", "fit_invivo_brain"):
        assert tg[c].notna().sum() > 7000
    both = tg[["fit_invivo_PE", "fit_invivo_heart"]].dropna()
    assert both.corr(method="spearman").iloc[0, 1] > 0.3


def test_the_serum_differential_is_a_new_axis_and_the_arms_are_not(tg):
    """Fitness in either serum IS fibroblast fitness again; the difference between them is not."""
    from scipy import stats

    def rho(a, b):
        d = tg[[a, b]].dropna()
        return abs(stats.spearmanr(d[a], d[b]).correlation)
    assert rho("fit_lipid_rich_p8", "fit_invitro_hff") > 0.7
    assert rho("fit_serum_differential_p8", "fit_invitro_hff") < 0.3


def test_ribosomal_proteins_are_needed_in_every_fitness_screen_that_was_added(tg):
    """The direction check. A screen whose sign is inverted looks like data and ranks genes
    backwards; ribosomal proteins are the class that cannot be dispensable."""
    product = tg["product"].fillna("")
    ribosomal = (product.str.contains("ribosomal protein", case=False)
                 & ~product.str.contains("mitochondrial|apicoplast|kinase|methyltransferase",
                                         case=False))
    for c in ("fit_lipid_rich_p8", "fit_lipid_limited_p8", "fit_complete_medium_2025",
              "fit_no_glucose", "fit_no_glutamine"):
        assert tg.loc[ribosomal, c].median() < tg.loc[~ribosomal, c].median() - 0.5, c


def test_the_genome_wide_decay_column_is_not_the_old_412_gene_one(tg):
    """Two measurements in different units that do not agree: kept apart, never averaged."""
    from scipy import stats
    new, old = "mrna_log2_remaining_4h_actinomycin", "mrna_remaining_5h_actinomycin"
    assert tg[new].notna().sum() > 5000 and tg[old].notna().sum() < 500
    d = tg[[new, old]].dropna()
    assert len(d) > 200 and abs(stats.spearmanr(d[new], d[old]).correlation) < 0.2


def test_stable_transcripts_and_high_efficiency_land_where_biology_says(tg):
    product = tg["product"].fillna("")
    ribosomal = (product.str.contains("ribosomal protein", case=False)
                 & ~product.str.contains("mitochondrial|apicoplast|kinase", case=False))
    assert tg.loc[ribosomal, "mrna_log2_remaining_4h_actinomycin"].median() > \
        tg.loc[~ribosomal, "mrna_log2_remaining_4h_actinomycin"].median()
    te = tg[["te302107_tachy_r1", "te302107_tachy_r2"]].mean(axis=1)
    assert te[ribosomal].median() > te[~ribosomal].median()


def test_the_two_new_translation_efficiency_replicates_agree_better_than_the_old_ones(tg):
    """The reason this study was taken: it is the most reproducible TE in the organism."""
    from scipy import stats
    d = tg[["te302107_tachy_r1", "te302107_tachy_r2"]].dropna()
    assert stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1]).correlation > 0.95


def test_upstream_start_codons_go_with_low_translation_efficiency(tg):
    """The 5' UTR columns are sequence features, and this is the relationship that makes them
    features of translation rather than noise."""
    from scipy import stats
    te = tg.filter(regex=r"^te\d").mean(axis=1)
    d = pd.DataFrame({"u": tg["utr5_n_uaugs"], "te": te}).dropna()
    assert stats.spearmanr(d["u"], d["te"]).correlation < -0.3


def test_the_bradyzoite_subtypes_are_all_bradyzoites_and_differ_from_one_another(tg):
    groups = [f"bzsub_{g}_expr" for g in "ABCDE"]
    percentile = tg[groups].rank(pct=True)
    for gene in ("TGME49_259020", "TGME49_291040"):            # BAG1, LDH2
        row = percentile[tg["gene_id"] == gene]
        assert (row.to_numpy() > 0.9).all(), gene
    sag1 = percentile[tg["gene_id"] == "TGME49_233460"]         # SAG1, a tachyzoite antigen
    assert (sag1.to_numpy() < 0.5).all()
    # Five columns, not one column copied five times: the subtypes differ, which is the finding.
    values = tg[groups].dropna()
    assert (values.nunique(axis=1) > 1).mean() > 0.5
    off_diagonal = values.corr(method="spearman").to_numpy()[np.triu_indices(5, 1)]
    assert off_diagonal.max() < 0.999


def test_a_bait_that_never_saw_a_protein_leaves_its_flag_missing(tg):
    """Absence of detection is not a non-hit: only a protein the bait tested can be called."""
    for bait in ("apicoplast", "mitochondrion", "er"):
        seen = tg[f"surface_{bait}_log2fc"].notna()
        flagged = tg[f"surface_{bait}_stringent"].notna()
        assert (flagged == seen).all(), bait


def test_iron_depletion_moves_the_protein_and_the_transcript_together(tg):
    from scipy import stats
    d = tg[["iron_depletion_protein_log2fc", "iron_depletion_rna_log2fc"]].dropna()
    assert len(d) > 2000
    assert stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1]).correlation > 0.2


def test_the_plasmodium_melting_points_behave_like_melting_points(pf):
    """Complexes melt together; chaperones are labile and glycolysis is stable."""
    tm = pf.set_index("gene_id")["melting_temperature_tm"].dropna()
    assert 2000 < len(tm) and 45 < tm.median() < 62
    product = pf.set_index("gene_id")["product"].reindex(tm.index).fillna("")
    hsp = product.str.contains("heat shock protein 70|heat shock protein 90", case=False)
    glyco = product.str.contains("glyceraldehyde-3-phosphate|enolase|pyruvate kinase", case=False)
    assert tm[hsp].median() < tm[glyco].median()


def test_fertility_is_sex_specific_and_the_known_genes_prove_it(pf):
    """HAP2 fails in males only. A column that lost the distinction would fail here."""
    d = pf.set_index("gene_id")
    hap2 = d.loc["PF3D7_1014200"]                              # HAP2/GCS1
    assert hap2["fertility_male"] < -2 and hap2["fertility_female"] > -1


def test_the_host_screen_recovers_the_pathway_its_paper_is_about(hosts):
    """SLC35A2 first of 20,010, and the glycosylation genes the paper excludes are not hits."""
    h = hosts.dropna(subset=["rhoptry_discharge_score"]).set_index("host_name")
    ranks = h["rhoptry_discharge_score"].rank(ascending=False)
    assert ranks.get("SLC35A2", 1e9) <= 5
    panel = [g for g in ("SLC35A2", "RFT1", "DPAGT1", "GFPT1", "MGAT1", "MGAT2") if g in ranks]
    assert np.median([ranks[g] for g in panel]) < 0.1 * len(ranks)
    for absent in ("B4GALT1", "FUT8", "ST6GAL1"):
        if absent in ranks:
            assert ranks[absent] > 0.1 * len(ranks), absent


def test_the_host_table_kept_every_protein_it_had_and_its_names(hosts):
    assert len(hosts) > 36000
    assert hosts["host_name"].notna().mean() > 0.9
    assert hosts["host_id"].is_unique


# --------------------------------------------------------------------------- leakage
def test_the_new_columns_are_grouped_with_what_they_restate(tg, pf):
    from starplast import search
    # Fitness in either serum, and the complete-medium arm of the carbon screen, ARE the fibroblast
    # screen measured again.
    banned = search.excluded_for(tg, "fit_invitro_hff")
    assert {"fit_lipid_rich_p8", "fit_lipid_limited_p8", "fit_complete_medium_2025"} <= banned
    # Its differential is a different quantity and is not swept up with them.
    assert "fit_serum_differential_p8" not in banned
    # One perturbation read two ways.
    assert "iron_depletion_rna_log2fc" in search.excluded_for(
        tg, "iron_depletion_protein_log2fc")
    # The arms of the carbon screen and the contrast between them share their replicates.
    assert {"fit_no_glutamine", "fit_glucose_dependence"} <= search.excluded_for(
        tg, "fit_no_glucose")
    # One screen, two sexes.
    assert "fertility_male" in search.excluded_for(pf, "fertility_female")
    # The spread of a fit goes with the value it describes.
    assert "melting_temperature_sd" in search.excluded_for(pf, "melting_temperature_tm")


def test_holding_out_a_new_target_never_leaves_the_target_itself_behind(tg, pf):
    from starplast import search
    for nodes, targets in ((tg, ("fit_glucose_dependence", "bzsub_A_expr", "utr5_n_uaugs",
                                 "surface_er_log2fc", "iron_depletion_protein_log2fc")),
                           (pf, ("melting_temperature_tm", "fertility_female"))):
        for target in targets:
            assert target in search.excluded_for(nodes, target), target


# --------------------------------------------------------------------------- the third wave
def test_the_plasmodium_chromatin_baits_reproduce_the_paper_and_separate(pf):
    """The two euchromatin counts are the paper's own, and the baits must tell the two chromatin
    states apart.

    Separation is judged on the proteins whose location is not in doubt, not on a threshold for how
    much the hit lists may overlap: 26 of about 85 proteins are called by both the HP1 and the
    H3K27ac bait, which is neither obviously wrong -- boundary proteins and abundant nuclear
    proteins are reached by any nuclear bait -- nor something this project can adjudicate. What it
    can check is that the four proteins whose chromatin state is settled land on the right side.
    """
    counts = {c: int(pf[c].sum()) for c in pf.columns if c.startswith("chromprox_")
              and c.endswith("_hit")}
    assert counts["chromprox_h3k27ac_hit"] == 99
    assert counts["chromprox_h3k4me3_hit"] == 48
    d = pf.set_index("gene_id")
    for gene in ("PF3D7_1220900", "PF3D7_1008000"):        # HP1 itself, and HDA2
        assert d.loc[gene, "chromprox_hp1_hit"] == 1, gene
        assert d.loc[gene, "chromprox_h3k27ac_hit"] == 0, gene
    for gene in ("PF3D7_1033700", "PF3D7_0823300"):        # BDP1 and GCN5, active chromatin
        assert d.loc[gene, "chromprox_h3k27ac_hit"] == 1, gene
        assert d.loc[gene, "chromprox_hp1_hit"] == 0, gene


def test_a_chromatin_bait_that_missed_a_protein_calls_nothing_about_it(pf):
    for bait in ("hp1", "h3k27ac", "h3k4me3", "centromere"):
        seen = pf[f"chromprox_{bait}_log2fc"].notna()
        assert (pf[f"chromprox_{bait}_hit"].notna() == seen).all(), bait


def test_target_engagement_ships_its_denominator(pf):
    """A protein hit twice out of five assays is not the evidence of twice out of 25."""
    both = pf[["engaged_n_compounds_tested", "engaged_n_compounds_hit"]].dropna()
    assert (both["engaged_n_compounds_hit"] <= both["engaged_n_compounds_tested"]).all()
    assert both["engaged_n_compounds_tested"].max() <= 25


def test_the_gametocyte_translatome_is_the_papers_705(pf):
    assert int(pf["gametocyte_newly_made"].sum()) == 705


def test_febrile_phosphorylation_favours_exported_proteins(pf):
    """The response is specific, not a global heat wobble: the proteins whose phosphorylation rises
    are the ones the parasite sends into the host cell."""
    rising = pf["febrile_phospho_n_sites_up"] > 0
    quantified = pf["febrile_phospho_max_abs_log2fc"].notna()
    exported = pf["is_exported"].fillna(False).astype(bool) if "is_exported" in pf else None
    if exported is None:
        pytest.skip("no export call in this table")
    assert exported[rising].mean() > 3 * exported[quantified].mean()


def test_m6a_stoichiometry_is_a_fraction_and_only_where_there_are_sites(pf):
    stoich = pf["m6a_canonical_stoichiometry"].dropna()
    assert ((stoich >= 0) & (stoich <= 1)).all()
    assert (pf.loc[stoich.index, "m6a_n_canonical_sites"] > 0).all()


def test_the_latency_classifier_is_two_hundred_genes(pf):
    assert int(pf["latency_classifier_member"].sum()) == 200


def test_every_berghei_transfer_is_held_out_together(pf):
    """The audit's finding: the blood-stage transfer predicted the transmission transfer at 0.58,
    nine robust standard deviations out, because both are the same berghei library."""
    from starplast import search
    for target in ("pb_transferred_phenotype", "fertility_female"):
        banned = search.excluded_for(pf, target)
        assert {"pb_transferred_phenotype", "pb_transferred_growth_rate", "fertility_female",
                "fertility_male"} <= banned, target


def test_withdrawing_a_nutrient_is_held_out_with_fibroblast_fitness(tg):
    """0.76 predictability: same library and cells with one nutrient removed is not an independent
    measurement of essentiality. The contrast between the arms is, and stays out."""
    from starplast import search
    banned = search.excluded_for(tg, "fit_invitro_hff")
    assert {"fit_no_glucose", "fit_no_glutamine"} <= banned
    assert "fit_glucose_dependence" not in banned
    assert "fit_serum_differential_p8" not in banned
