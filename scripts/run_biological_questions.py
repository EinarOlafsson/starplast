"""Run the tested subset of `instructions/open/59_biological_questions.md`.

The instruction file writes 100 biological questions this software can answer. This script runs the
subset that was actually executed, saves every output table, and writes one summary table plus an
executed notebook so the whole thing reproduces.

    python scripts/run_biological_questions.py --out results/questions_2026_09_30

Each question is one `QUESTION`: the organism, the strategy key, the settings, and whether its
self-test is run too. Nothing here changes a strategy, a calibration number or a shipped table --
it only calls `starplast.strategies.run` and `.test` the way a user would from the Strategies tab.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings
from dataclasses import dataclass, field

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

warnings.filterwarnings("ignore")

from starplast import organisms  # noqa: E402
from starplast import strategies as S  # noqa: E402

# The two parasite spaces, in registry order, so this script names no organism itself.
TG, PF = [s.code for s in organisms.SPACES.values() if s.kind == organisms.PARASITE][:2]


@dataclass
class Question:
    qid: str
    organism: str
    strategy: str
    settings: dict
    question: str
    anchor: str
    also_test: bool = False
    tables: tuple = ()
    top: int = 15
    genes_from: tuple = ()          # (column, value): build `genes` from a shipped label
    genes_limit: int = 0
    settings_genes: tuple = ()      # an explicit gene list, given as accessions


def _gene_set(ctx, column, value, limit=0):
    mask = (ctx.nodes[column].astype("object") == value).to_numpy()
    ids = [str(g) for g in ctx.gene_ids[mask]]
    return ids[:limit] if limit else ids


QUESTIONS: list[Question] = [
    # ---- Toxoplasma: what the spatial proteome supports ------------------------------------
    Question("Q01", TG, "recoverability_atlas", dict(target="compartment"),
             "Which hyperLOPIT compartments can the rest of the data rediscover on its own?",
             "PMID:33053376", also_test=True, tables=("categories",)),
    Question("Q02", TG, "feature_knn", dict(target="compartment"),
             "What compartment do the hyperLOPIT-unassigned proteins behave like?",
             "PMID:39082802", also_test=True, tables=("calls",)),
    Question("Q03", TG, "structural_homology", dict(target="compartment"),
             "Which unannotated proteins can be placed by fold alone, where sequence homology fails?",
             "PMID:40066064", also_test=True, tables=("calls",)),
    Question("Q04", TG, "physical_partners", dict(target="compartment"),
             "Which proteins can be placed by the proteins they physically touch?",
             "PMID:40874616", tables=("calls",)),
    Question("Q05", TG, "cluster_guilt", dict(target="compartment"),
             "Which map clusters are compartment-pure, and which unassigned genes sit in them?",
             "PMID:38270431", tables=("enrichment", "calls")),
    Question("Q06", TG, "conformal_calls", dict(target="compartment", alpha=0.1),
             "Which localization calls come with a stated 10% error rate?",
             "PMID:41396985", also_test=True, tables=("calls",)),
    Question("Q07", TG, "stacking", dict(target="compartment"),
             "Which kind of evidence does a learner trust most for localization?",
             "PMID:33053376", tables=("trust by evidence", "calls")),
    Question("Q08", TG, "random_forest", dict(target="compartment"),
             "Which measurements actually define a compartment label?",
             "PMID:33053376", tables=("what defines the label", "calls")),
    Question("Q09", TG, "block_ablation", dict(target="compartment"),
             "Which kind of evidence carries localization, and which is redundant?",
             "PMID:33053376", tables=("evidence",)),
    Question("Q10", TG, "graph_convolution", dict(target="compartment"),
             "Does smoothing measurements along the networks improve localization calls?",
             "PMID:38900844", tables=("calls",)),
    Question("Q11", TG, "supervised_classifier", dict(target="compartment"),
             "What does a plain classifier trained on known compartments call the rest?",
             "PMID:33053376", tables=("calls",)),
    Question("Q12", TG, "triangulation", dict(target="compartment", min_agree=2),
             "Which localization calls survive agreement between three independent methods?",
             "PMID:39082802", tables=("calls",)),
    Question("Q13", TG, "understudied_first", dict(target="compartment", min_agree=2),
             "Which never-published proteins get a confident compartment call?",
             "PMID:41582196", tables=("calls",)),
    Question("Q14", TG, "label_outliers", dict(target="compartment", top=200),
             "Which proteins behave unlike their own compartment -- dual-localized or mis-assigned?",
             "PMID:40874616", also_test=True, tables=("surprising genes",)),
    Question("Q15", TG, "layer_propagation", dict(target="compartment", layer="xlms"),
             "How far does compartment information diffuse through the crosslink network alone?",
             "PMID:40874616", tables=("calls",)),
    Question("Q16", TG, "layer_vote", dict(target="compartment"),
             "Do all nine measured networks together beat any one of them for localization?",
             "PMID:38900844", tables=("source weights", "calls")),
    Question("Q17", TG, "map_neighbours", dict(target="compartment"),
             "Can a gene be localized from its neighbours on a single UMAP?",
             "PMID:33053376", tables=("calls",)),

    # ---- Toxoplasma: gene lists a biologist brings -----------------------------------------
    Question("Q18", TG, "set_enrichment", dict(),
             "What do the hyperLOPIT dense-granule proteins have in common across every layer?",
             "PMID:36515555", tables=("profile",), genes_from=("compartment", "dense granules")),
    Question("Q19", TG, "positive_unlabeled", dict(bags=15, top=300),
             "Which uncharacterised proteins most resemble the known dense-granule proteins?",
             "PMID:36515555", tables=("candidates",),
             genes_from=("compartment", "dense granules")),
    Question("Q20", TG, "seed_expansion", dict(restart=0.3, top=300),
             "Which genes does the network walk attach to the rhoptry-1 proteins?",
             "PMID:41975467", tables=("candidates",), genes_from=("compartment", "rhoptries 1")),
    Question("Q21", TG, "geneset_hunt", dict(),
             "Is there any map on which the micronemal proteins fall out as one cluster?",
             "PMID:42470170", tables=("maps",), genes_from=("compartment", "micronemes")),
    Question("Q22", TG, "neighbour_space", dict(k=10, top=300),
             "What sits in GRA17's neighbourhood when every layer is put in one space?",
             "PMID:37498952", tables=("neighbourhoods",),
             settings_genes=("TGME49_222170", "TGME49_254470")),

    # ---- Toxoplasma: links the field has not written down -----------------------------------
    Question("Q23", TG, "link_prediction", dict(layer="xlms", top=300),
             "Which contacts did the crosslinking interactome most likely miss?",
             "PMID:40874616", also_test=True, tables=("predicted links",)),
    Question("Q24", TG, "unwritten_links", dict(min_layers=3, top=500),
             "Which gene pairs does the data link in three layers that no paper mentions together?",
             "PMID:28276701", tables=("unwritten pairs",)),
    Question("Q25", TG, "attention_correction", dict(target="compartment", layer="comention",
                                                     top=200),
             "Which co-mentioned pairs are linked more than the two genes' fame predicts?",
             "PMID:41582196", tables=("corrected ranking",)),
    Question("Q26", TG, "network_training", dict(layer="xlms", top=300),
             "Which crosslink edges would a model trained on all the networks add?",
             "PMID:40874616", tables=("gaps",)),

    # ---- Toxoplasma: fitness and conditions -------------------------------------------------
    Question("Q27", TG, "trait_regression", dict(target="fit_invitro_hff", top=200),
             "Which genes are far more or less fitness-conferring than everything else predicts?",
             "PMID:27594426", also_test=True, tables=("surprising genes",)),
    Question("Q28", TG, "trait_regression", dict(target="fit_invivo_brain", top=200),
             "Which genes defy the prediction of in-vivo brain fitness?",
             "PMID:41935076", tables=("surprising genes",)),
    Question("Q29", TG, "trait_regression", dict(target="fit_density_dependence", top=200),
             "Do genes that matter only at high parasite density look like anything else in the data?",
             "PMID:34023299", also_test=True, tables=("surprising genes",)),
    Question("Q30", TG, "masked_imputation", dict(target="fit_ifng"),
             "Where can IFN-gamma macrophage fitness be filled in honestly for unmeasured genes?",
             "PMID:33067458", also_test=True, tables=("column reliability", "filled")),
    Question("Q31", TG, "conformal_values", dict(target="fit_invitro_hff", alpha=0.1),
             "Which fitness predictions come with an interval that actually holds 90% of the time?",
             "PMID:27594426", also_test=True, tables=("surprises", "predicted intervals")),
    Question("Q32", TG, "condition_shift", dict(condition="fit_ifng", baseline="fit_naive_bmdm"),
             "Which genes matter more under interferon-gamma than in a naive macrophage, and why?",
             "PMID:33067458", tables=("shifted genes",)),

    # ---- Toxoplasma: families, strata and combinations --------------------------------------
    Question("Q33", TG, "paralog_divergence", dict(target="compartment", min_shared=10),
             "Which paralogue pairs changed compartment -- a family that divided its labour?",
             "PMID:42113865", tables=("paralog pairs",)),
    Question("Q34", TG, "stratum_focus", dict(target="compartment", stratum="lineage-specific"),
             "Does inference still work on the genes orthology cannot reach?",
             "PMID:35196325", also_test=True, tables=("calls",)),
    Question("Q35", TG, "conjunctions", dict(a="compartment", b="stage_enriched_derived"),
             "Are there kinds of gene defined by compartment and life stage at once?",
             "PMID:37081202", tables=("conjunctions",)),
    Question("Q36", TG, "blind_battery", dict(),
             "Tune a map with no labels at all: what does it turn out to encode?",
             "PMID:33053376", tables=("held-out features",)),
    Question("Q37", TG, "holdout_search", dict(target="compartment"),
             "Is there a map setting under which a hidden compartment falls out as clusters?",
             "PMID:33053376", tables=("categories",)),
    Question("Q38", TG, "multiplex_modules", dict(target="compartment"),
             "Which communities do several measured networks agree on?",
             "PMID:38900844", tables=("modules",)),
    Question("Q39", TG, "ortholog_transfer", dict(target="fit_invitro_hff", source="piggybac_mis"),
             "Does Plasmodium essentiality predict Toxoplasma fitness through orthogroups?",
             "PMID:39913589", also_test=True, tables=("transferred",)),

    # ---- Plasmodium falciparum --------------------------------------------------------------
    Question("Q40", PF, "recoverability_atlas", dict(target="lopit_pf_location"),
             "Does the schizont spatial proteome fall out of the rest of the Plasmodium data?",
             "PMID:42218142", also_test=True, tables=("categories",)),
    Question("Q41", PF, "structural_homology", dict(target="lopit_pf_location"),
             "Which Plasmodium hypotheticals can be localized by fold?",
             "PMID:42218142", tables=("calls",)),
    Question("Q42", PF, "conformal_calls", dict(target="lopit_pf_location", alpha=0.1),
             "Which Plasmodium localization calls carry a stated error rate?",
             "PMID:41396985", also_test=True, tables=("calls",)),
    Question("Q43", PF, "link_prediction", dict(layer="ip_ms", top=300),
             "Which protein contacts did the Plasmodium pulldowns miss?",
             "PMID:41482054", also_test=True, tables=("predicted links",)),
    Question("Q44", PF, "condition_shift", dict(condition="expr_gametocyte_v",
                                                baseline="expr_schizont"),
             "What changes between schizont and mature gametocyte beyond overall expression?",
             "PMID:41482054", tables=("shifted genes",)),
    Question("Q45", PF, "trait_regression", dict(target="piggybac_mis", top=200),
             "Which Plasmodium genes are more or less essential than the data predicts?",
             "PMID:39913589", also_test=True, tables=("surprising genes",)),
    Question("Q46", PF, "seed_expansion", dict(restart=0.3, top=300),
             "Which genes join the measured apicoplast proteins when the networks are walked?",
             "PMID:36563485", tables=("candidates",),
             genes_from=("lopit_pf_location", "apicoplast")),
    Question("Q47", PF, "set_enrichment", dict(),
             "What do the exported Plasmodium proteins share across every layer?",
             "PMID:41035680", tables=("profile",), genes_from=("is_exported", True)),
    Question("Q48", PF, "feature_knn", dict(target="lopit_pf_location"),
             "Does nearest-neighbour localization work in Plasmodium the way it does in Toxoplasma?",
             "PMID:42218142", also_test=True, tables=("calls",)),
    Question("Q49", PF, "label_outliers", dict(target="lopit_pf_location", top=200),
             "Which Plasmodium proteins behave unlike their measured compartment?",
             "PMID:42218142", tables=("surprising genes",)),
    Question("Q50", PF, "ortholog_transfer", dict(target="piggybac_mis", source="fit_invitro_hff"),
             "Does Toxoplasma fitness predict Plasmodium essentiality through orthogroups?",
             "PMID:39913589", also_test=True, tables=("transferred",)),
]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/questions_2026_09_30")
    ap.add_argument("--only", default="", help="comma-separated question ids")
    ap.add_argument("--notebook", default="notebooks/biological_questions_2026_09_30.ipynb")
    ap.add_argument("--no-notebook", action="store_true")
    args = ap.parse_args(argv)

    os.makedirs(args.out, exist_ok=True)
    wanted = {q.strip() for q in args.only.split(",") if q.strip()}
    todo = [q for q in QUESTIONS if not wanted or q.qid in wanted]

    contexts = {}
    rows = []
    for q in todo:
        ctx = contexts.setdefault(q.organism, S.Context.shipped(q.organism))
        settings = dict(q.settings)
        if q.genes_from:
            settings["genes"] = _gene_set(ctx, *q.genes_from, limit=q.genes_limit)
        if q.settings_genes:
            settings["genes"] = list(q.settings_genes)
        row = {"qid": q.qid, "organism": q.organism, "strategy": q.strategy,
               "question": q.question, "anchor": q.anchor,
               "grade": S.calibration(q.strategy, q.organism).get("grade"),
               "settings": json.dumps({k: (f"{len(v)} genes" if k == "genes" else v)
                                       for k, v in settings.items()})}
        t0 = time.monotonic()
        try:
            result = S.get(q.strategy).run(ctx, **settings)
            row["seconds"] = round(time.monotonic() - t0, 1)
            row["summary"] = result.summary
            for name, df in result.tables.items():
                safe = name.replace(" ", "_").replace("/", "_")
                df.to_csv(os.path.join(args.out, f"{q.qid}_{q.strategy}_{safe}.csv"), index=False)
                row[f"rows[{name}]"] = len(df)
            lead = next((n for n in q.tables if n in result.tables), None) or \
                next(iter(result.tables), None)
            if lead is not None:
                row["lead_table"] = lead
                row["head"] = result.tables[lead].head(q.top).to_dict("records")
            row["status"] = "ran"
        except Exception as exc:                                   # noqa: BLE001
            row["seconds"] = round(time.monotonic() - t0, 1)
            row["status"] = "error"
            row["summary"] = f"{type(exc).__name__}: {exc}"
        if q.also_test and row["status"] == "ran":
            t0 = time.monotonic()
            try:
                test = S.get(q.strategy).test(ctx, **settings)
                row["verdict"] = test.verdict
                row["metric"] = test.metric
                row["observed"] = test.observed
                row["chance"] = test.null_mean
                row["bar"] = test.null_high
                row["p"] = test.p_value
                row["skill"] = test.skill
                row["n_hidden"] = test.n_hidden
                row["test_summary"] = test.summary()
                row["scorecard"] = {k: (round(v, 4) if isinstance(v, float) else v)
                                    for k, v in (test.scorecard or {}).items()}
                row["test_seconds"] = round(time.monotonic() - t0, 1)
            except Exception as exc:                               # noqa: BLE001
                row["test_summary"] = f"{type(exc).__name__}: {exc}"
        rows.append(row)
        print(f"{q.qid} {q.organism} {q.strategy} {row['status']} {row['seconds']}s "
              f"{str(row.get('summary'))[:150]}", flush=True)
        if "verdict" in row:
            print("     test: {} metric={} observed={} chance={} bar={} skill={}".format(
                row.get("verdict"), row.get("metric"), row.get("observed"),
                row.get("chance"), row.get("bar"), row.get("skill")), flush=True)

    with open(os.path.join(args.out, "summary.json"), "w") as fh:
        json.dump(rows, fh, indent=1, default=str)

    import pandas as pd
    flat = pd.DataFrame([{k: v for k, v in r.items() if k != "head"} for r in rows])
    flat.to_csv(os.path.join(args.out, "summary.csv"), index=False)
    print("wrote", args.out)

    if not args.no_notebook:
        write_notebook(args.notebook, args.out, [q.qid for q in todo])
    return rows


def write_notebook(path, out, qids):
    from notebook_runner import ExecutedNotebook

    nb = ExecutedNotebook("Biological questions, 2026-09-30")
    nb.md("# 100 biological questions, and the subset that was run",
          "`instructions/open/59_biological_questions.md` writes 100 questions about *Toxoplasma "
          "gondii* and *Plasmodium falciparum* that this software can answer, each with its entry "
          "point, its strategy, its settings and a resolved literature anchor. This notebook runs "
          "the subset that was executed and shows the real output -- the top calls, the metric and "
          "its chance level -- so every claim in `docs/questions.md` can be checked.",
          "Nothing here changes a strategy, a calibration number or a shipped table.")
    nb.code("import os, sys, warnings",
            "warnings.filterwarnings('ignore')",
            "sys.path.insert(0, os.path.abspath('scripts'))",
            "from run_biological_questions import QUESTIONS",
            "len(QUESTIONS)")
    nb.md("## The questions that were run",
          "Each row is one question: the organism, the strategy, its settings and the paper that "
          "makes it a real question.")
    nb.code("import pandas as pd",
            "pd.set_option('display.width', 200)",
            "pd.DataFrame([{'qid': q.qid, 'organism': q.organism, 'strategy': q.strategy,"
            " 'question': q.question, 'anchor': q.anchor} for q in QUESTIONS])")
    nb.md("## The recorded outputs",
          "`summary.csv` holds one row per question with the strategy's own summary sentence, the "
          "self-test verdict where a test was run, and the number of rows in each output table. "
          "The per-question CSVs beside it hold the calls themselves.")
    nb.code(f"summary = pd.read_csv({out + '/summary.csv'!r})",
            "summary[['qid', 'organism', 'strategy', 'grade', 'seconds', 'status']]")
    nb.code("cols = [c for c in ('qid', 'verdict', 'observed', 'chance', 'skill') "
            "if c in summary.columns]",
            "summary.loc[summary.get('verdict').notna(), cols] if 'verdict' in summary "
            "else 'no self-tests in this run'")
    nb.md("## Reproducing one question",
          "Every question is one call. This is Q01 -- which hyperLOPIT compartments the rest of the "
          "data can rediscover -- run again here so the notebook is not just a table of numbers.")
    nb.code("from starplast import strategies as S",
            "ctx = S.Context.shipped('Tg')",
            "q01 = S.run('recoverability_atlas', ctx=ctx, target='compartment')",
            "print(q01.summary)",
            "q01.tables['categories'].sort_values('auroc', ascending=False).head(12)")
    nb.md("To re-run the whole subset: "
          "`python scripts/run_biological_questions.py --out " + out + "`.")
    nb.write(path)
    print("wrote", path)


if __name__ == "__main__":
    main()
